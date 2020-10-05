/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_processing.cc
 *
 *  \brief main routines for processing a halo with the Subfind algorithm
 */

#include "gadgetconfig.h"

#ifdef SUBFIND

#include <gsl/gsl_math.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fof/fof.h"
#include "../gravtree/gravtree.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/cxxsort.h"
#include "../sort/parallel_sort.h"
#include "../sort/peano.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

template <typename partset>
void fof<partset>::subfind_processing(domain<partset> *SubDomain, domain_options mode)
{
  double t0 = Logs.second();

  if(mode == COLL_SUBFIND)
    {
      /* make a sanity check: We should have exactly 1 group, stored on the root of the processor subset if we collectively do a group
       */
      if(SubThisTask == 0 && Ngroups != 1)
        Terminate("Ngroups=%d != 1  SubNTask=%d SubThisTask=%d", Ngroups, SubNTask, SubThisTask);

      if(SubThisTask != 0 && Ngroups != 0)
        Terminate("Ngroups=%d != 0  SubNTask=%d SubThisTask=%d", Ngroups, SubNTask, SubThisTask);

      subfind_collective_printf("SUBFIND: root-task=%d: Collectively doing halo %lld of length  %lld  on  %d  processors.\n", ThisTask,
                                Group[0].GroupNr, (long long)Group[0].Len, SubNTask);

      if(SubThisTask == 0)
        {
          GroupNr = Group[0].GroupNr;
          Ascale  = Group[0].Ascale;
        }

      for(int i = 0; i < Tp->NumPart; i++)
        if(Tp->PS[i].GroupNr.get() == GroupNr)
          Tp->PS[i].DomainFlag = 1;
        else
          Tp->PS[i].DomainFlag = 0;

      /* tell everybody in the set the group number */
      MPI_Bcast(&GroupNr, 1, MPI_LONG_LONG, 0, SubComm);

      /* tell everybody in the set the group's scale factor */
      MPI_Bcast(&Ascale, 1, MPI_DOUBLE, 0, SubComm);
    }
  else
    {
      for(int i = 0; i < Tp->NumPart; i++)
        Tp->PS[i].DomainFlag = 1;

      if(SubNTask != 1)
        Terminate("Strange: SubNTask=%d  Ngroups=%d  SubThisTask=%d (expect to be a single processor here)", SubNTask, Ngroups,
                  SubThisTask);
    }

  /* Create a domain decomposition for the sub-communicator and the particles in it.
   * For the serial algorithm, this will be trivial, for collectively treated groups, only particles in the group get a gravity weight.
   * The outline of the toplevel tree nodes resulting from the domain decomposition can be used to create a gravity tree.
   */

  if(SubDomain->NumNodes != 0)
    Terminate("SubDomain.NumNodes=%d\n", SubDomain->NumNodes);

  double ta = Logs.second();
  SubDomain->domain_decomposition(mode);
  double tb = Logs.second();

  mpi_printf("SUBFIND: subdomain decomposition took %g sec\n", Logs.timediff(ta, tb));

  if(mode == COLL_SUBFIND)
    SubDomain->particle_exchange_based_on_PS(SubComm);

  for(int i = 0; i < Tp->NumPart; i++)
    Tp->PS[i].SubRankInGr = INT_MAX; /* set a default that is larger than any reasonable group number */

  /* now let us sort according to GroupNr and Density. This step will temporarily break the association with SphP[] and other arrays!
   */
  submp = (submp_data *)Mem.mymalloc_movable(&submp, "submp", sizeof(submp_data) * Tp->NumPart);
  for(int i = 0; i < Tp->NumPart; i++)
    {
      submp[i].index   = i;
      submp[i].GroupNr = Tp->PS[i].GroupNr.get();
#ifndef SUBFIND_HBT
      submp[i].DM_Density = Tp->PS[i].u.s.u.DM_Density;
#endif
    }
  mycxxsort(submp, submp + Tp->NumPart, subfind_compare_submp_GroupNr_DM_Density);

  /* In this index list, we store the indices of the local particles in the group that we currently need to process (it is for sure
   * shorter than NumPart) */
  IndexList = (int *)Mem.mymalloc_movable(&IndexList, "IndexList", Tp->NumPart * sizeof(int));

  MPI_Comm SingleComm;
  int thistask;
  MPI_Comm_rank(SubDomain->Communicator, &thistask);
  MPI_Comm_split(SubDomain->Communicator, thistask, thistask, &SingleComm);  // create a communicator for single ranks

  /* prepare a domain decomposition for unbinding locally independent ones */
  domain<partset> SingleDomain{SingleComm, Tp};

  if(SingleDomain.NumNodes != 0)
    Terminate("SubDomain.NumNodes=%d\n", SingleDomain.NumNodes);

  double taa = Logs.second();
  SingleDomain.domain_decomposition(SERIAL_SUBFIND);
  double tbb = Logs.second();

  mpi_printf("SUBFIND: serial subfind subdomain decomposition took %g sec\n", Logs.timediff(taa, tbb));

  double ta0 = Logs.second();
  if(mode == COLL_SUBFIND)
    {
      /* determine the number of local group particles, and fill the list of the indices */
      NumPartGroup = 0;
      for(int i = 0; i < Tp->NumPart; i++)
        if(Tp->PS[i].GroupNr.get() == GroupNr)
          IndexList[NumPartGroup++] = i;

          /* call the processing of the group */
#ifdef SUBFIND_HBT
      subfind_hbt_single_group(SubDomain, &SingleDomain, mode, 0);
#else
      subfind_process_single_group(SubDomain, &SingleDomain, mode, 0);
#endif
    }
  else
    {
      int i = 0;
      for(int gr = 0; gr < Ngroups; gr++) /* process all local groups */
        {
          GroupNr = Group[gr].GroupNr;
          Ascale  = Group[gr].Ascale;

          /* determine the number of local group particles, and set up the list of the indices */
          NumPartGroup = 0;
          for(; i < Tp->NumPart; i++)
            if(Tp->PS[submp[i].index].GroupNr.get() == GroupNr)
              IndexList[NumPartGroup++] = submp[i].index;
            else
              break;

              /* do local group with Group[] index 'gr' */
#ifdef SUBFIND_HBT
          subfind_hbt_single_group(SubDomain, &SingleDomain, mode, gr);
#else
          subfind_process_single_group(SubDomain, &SingleDomain, mode, gr);
#endif
        }
    }
  double tb0 = Logs.second();
  mpi_printf("SUBFIND: subfind_hbt_single_group() processing for Ngroups=%d took %g sec\n", Ngroups, Logs.timediff(ta0, tb0));

  SingleDomain.domain_free();
  MPI_Comm_free(&SingleComm);

  Mem.myfree(IndexList);
  Mem.myfree(submp);

  SubDomain->domain_free();

  double t1 = Logs.second();

  subfind_collective_printf("SUBFIND: root-task=%d: Collective processing of halo %d took %g\n", ThisTask, Group[0].GroupNr,
                            Logs.timediff(t0, t1));

  if(!(mode == COLL_SUBFIND) && ThisTask == 0)
    mpi_printf("SUBFIND: root-task=%d: Serial processing of halo %d took %g\n", ThisTask, Group[0].GroupNr, Logs.timediff(t0, t1));
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
