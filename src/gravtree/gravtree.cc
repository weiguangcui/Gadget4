/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file gravtree.cc
 *
 *  \brief driver routines for computing the (short-range) gravity
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <atomic>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../gravity/ewald.h"
#include "../gravtree/gravtree.h"
#include "../gravtree/gwalk.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../pm/pm.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/*!
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local particles, and particles are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */

#ifdef PMGRID

/*! \brief initializes the short range table
 *
 * The short range table contains the complementary error function
 * needed for the computation of the short range part of the gravity
 * force/potential in case of the TreePM algorithm.
 */
template <typename partset>
void gravtree<partset>::short_range_init(void)
{
  for(int i = 0; i <= NTAB; i++)
    {
      double u = (RCUT / 2.0) / NTAB * i;

      /*
       *  fac0 contains    g = g(u) = (erfc(u) - 1)/u
       *  fac1 contains    g'
       *  fac2 contains    g'' - g'/u
       *  fac3 contains    g''' - 3g''/u + 3g'/u^2
       *  fac4 contains    g'''' - 6*g'''/u + 15*g''/u^2 - 15*g'/u^3
       *  fac5 contains    g''''' - 10*g''''/u + 45*g'''/u^2 - 105*g''/u^3 +  105*g'/u^4
       */

      shortrange_factors[i].fac0 = (u > 0) ? (erfc(u) - 1.0) / u : -2.0 / sqrt(M_PI);
      shortrange_factors[i].fac1 = (u > 0) ? ((1.0 - erfc(u)) - 2.0 * u / sqrt(M_PI) * exp(-u * u)) / (u * u) : 0;
#if(MULTIPOLE_ORDER >= 2)
      shortrange_factors[i].fac2 =
          (u > 0) ? (-3.0 * (1.0 - erfc(u)) + (6 * u + 4 * pow(u, 3)) / sqrt(M_PI) * exp(-u * u)) / pow(u, 3) : 0;
#endif
#if(MULTIPOLE_ORDER >= 3)
      shortrange_factors[i].fac3 =
          (u > 0) ? (15.0 * (1.0 - erfc(u)) - (30 * u + 20 * pow(u, 3) + 8 * pow(u, 5)) / sqrt(M_PI) * exp(-u * u)) / pow(u, 4) : 0;
#endif
#if(MULTIPOLE_ORDER >= 4)
      shortrange_factors[i].fac4 =
          (u > 0) ? (-105.0 * (1.0 - erfc(u)) +
                     (210.0 * u + 140.0 * pow(u, 3) + 56.0 * pow(u, 5) + 16.0 * pow(u, 7)) / sqrt(M_PI) * exp(-u * u)) /
                        pow(u, 5)
                  : 0;
#endif
#if(MULTIPOLE_ORDER >= 5)
      shortrange_factors[i].fac5 = (u > 0) ? (945.0 * (1.0 - erfc(u)) - (1890.0 * u + 1260.0 * pow(u, 3) + 504.0 * pow(u, 5) +
                                                                         144.0 * pow(u, 7) + 32.0 * pow(u, 9)) /
                                                                            sqrt(M_PI) * exp(-u * u)) /
                                                 pow(u, 6)
                                           : 0;
#endif
    }
}

template <typename partset>
void gravtree<partset>::set_mesh_factors(void)
{
  for(int i = 0; i < 2; i++)
    {
      mf[i].rcut2       = Tp->Rcut[i] * Tp->Rcut[i];
      double dblrcut[3] = {Tp->Rcut[i], Tp->Rcut[i],
                           Tp->Rcut[i]}; /* for stretched boxes, the conversion may be different in each dimension */
      Tp->pos_to_signedintpos(dblrcut, (MySignedIntPosType *)mf[i].intrcut);

      if(Tp->Asmth[i] > 0)
        mf[i].asmthinv1 = 0.5 / Tp->Asmth[i];
      else
        mf[i].asmthinv1 = 0;

      mf[i].asmthinv2 = mf[i].asmthinv1 * mf[i].asmthinv1;
#if(MULTIPOLE_ORDER >= 2)
      mf[i].asmthinv3 = mf[i].asmthinv2 * mf[i].asmthinv1;
#endif
#if(MULTIPOLE_ORDER >= 3)
      mf[i].asmthinv4 = mf[i].asmthinv3 * mf[i].asmthinv1;
#endif
#if(MULTIPOLE_ORDER >= 4)
      mf[i].asmthinv5 = mf[i].asmthinv4 * mf[i].asmthinv1;
#endif
#if(MULTIPOLE_ORDER >= 5)
      mf[i].asmthinv6 = mf[i].asmthinv5 * mf[i].asmthinv1;
#endif

      mf[i].asmthfac = mf[i].asmthinv1 * (NTAB / (RCUT / 2.0));
    }
}

#endif

template <typename partset>
void gravtree<partset>::gravity_exchange_forces(void)
{
  int *send_count  = (int *)Mem.mymalloc_movable(&send_count, "send_count", sizeof(int) * D->NTask);
  int *send_offset = (int *)Mem.mymalloc_movable(&send_offset, "send_offset", sizeof(int) * D->NTask);
  int *recv_count  = (int *)Mem.mymalloc_movable(&recv_count, "recv_count", sizeof(int) * D->NTask);
  int *recv_offset = (int *)Mem.mymalloc_movable(&recv_offset, "tecv_offset", sizeof(int) * D->NTask);

  /* now communicate the forces in ResultsActiveImported */
  for(int j = 0; j < D->NTask; j++)
    recv_count[j] = 0;

  int n = 0, k = 0;

  for(int i = 0; i < D->NTask; i++)
    for(int j = 0; j < Recv_count[i]; j++, n++) /* Note that we access Tree.Recv_count here */
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Points[n].ActiveFlag)
#endif
          {
            ResultsActiveImported[k].index = Points[n].index;
            recv_count[i]++;
            k++;
          }
      }
  myMPI_Alltoall(recv_count, 1, MPI_INT, send_count, 1, MPI_INT, D->Communicator);

  recv_offset[0] = 0;
  send_offset[0] = 0;

  int Nexport = 0;

  for(int j = 0; j < D->NTask; j++)
    {
      Nexport += send_count[j];
      if(j > 0)
        {
          send_offset[j] = send_offset[j - 1] + send_count[j - 1];
          recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
        }
    }

  resultsactiveimported_data *tmp_results =
      (resultsactiveimported_data *)Mem.mymalloc("tmp_results", Nexport * sizeof(resultsactiveimported_data));

  /* exchange  data */
  for(int ngrp = 1; ngrp < (1 << D->PTask); ngrp++)
    {
      int recvTask = D->ThisTask ^ ngrp;

      if(recvTask < D->NTask)
        {
          if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
            {
              myMPI_Sendrecv(&ResultsActiveImported[recv_offset[recvTask]], recv_count[recvTask] * sizeof(resultsactiveimported_data),
                             MPI_BYTE, recvTask, TAG_FOF_A, &tmp_results[send_offset[recvTask]],
                             send_count[recvTask] * sizeof(resultsactiveimported_data), MPI_BYTE, recvTask, TAG_FOF_A, D->Communicator,
                             MPI_STATUS_IGNORE);
            }
        }
    }
  for(int i = 0; i < Nexport; i++)
    {
      int target = tmp_results[i].index;

      for(int k = 0; k < 3; k++)
        Tp->P[target].GravAccel[k] += tmp_results[i].GravAccel[k];
#ifdef EVALPOTENTIAL
      Tp->P[target].Potential += tmp_results[i].Potential;
#endif

      if(MeasureCostFlag)
        Tp->P[target].GravCost += tmp_results[i].GravCost;
    }
  Mem.myfree(tmp_results);

  Mem.myfree(recv_offset);
  Mem.myfree(recv_count);
  Mem.myfree(send_offset);
  Mem.myfree(send_count);
}

/* make sure that we instantiate the template */
#include "../data/simparticles.h"
template class gravtree<simparticles>;
