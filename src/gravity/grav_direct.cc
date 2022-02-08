/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file grav_direct.cc
 *
 *  \brief calculates forces through direct summation
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../gravtree/gravtree.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

#ifdef ALLOW_DIRECT_SUMMATION

/*! \brief This function computes the gravitational forces for all active particles through direct summation.
 *
 */
template <>
void gravtree<simparticles>::gravity_direct(simparticles *Sp, domain<simparticles> *D, int timebin)
{
  TIMER_START(CPU_TREEDIRECT);

  D->mpi_printf("GRAVDIRECT: direct summation.  (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());

  double tstart = Logs.second();

  int *Send_count  = (int *)Mem.mymalloc_movable(&Send_count, "Send_count", sizeof(int) * D->NTask);
  int *Send_offset = (int *)Mem.mymalloc_movable(&Send_offset, "Send_offset", sizeof(int) * D->NTask);
  int *Recv_count  = (int *)Mem.mymalloc_movable(&Recv_count, "Recv_count", sizeof(int) * D->NTask);
  int *Recv_offset = (int *)Mem.mymalloc_movable(&Recv_offset, "Recv_offset", sizeof(int) * D->NTask);

  struct directdata
  {
    MyIntPosType IntPos[3];
    MyDouble Mass;
    unsigned char Type;
#if NSOFTCLASSES > 1
    unsigned char SofteningClass;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
    unsigned char InsideOutsideFlag : 1;
#endif
  };
  directdata *DirectDataIn, *DirectDataAll;

  struct accdata
  {
    MyFloat Acc[3];
#ifdef EVALPOTENTIAL
    MyFloat Potential;
#endif
  };
  accdata *DirectAccOut, *DirectAccIn;

  DirectDataIn = (directdata *)Mem.mymalloc("DirectDataIn", Sp->TimeBinsGravity.NActiveParticles * sizeof(directdata));

  int nforces = 0;

  for(int idx = 0; idx < Sp->TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = Sp->TimeBinsGravity.ActiveParticleList[idx];

      for(int k = 0; k < 3; k++)
        DirectDataIn[nforces].IntPos[k] = Sp->P[i].IntPos[k];

      DirectDataIn[nforces].Mass = Sp->P[i].getMass();

      DirectDataIn[nforces].Type = Sp->P[i].getType();
#if NSOFTCLASSES > 1
      DirectDataIn[nforces].SofteningClass = Sp->P[i].getSofteningClass();
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      DirectDataIn[nforces].InsideOutsideFlag = Sp->P[i].InsideOutsideFlag;
#endif
      nforces++;
    }

  MPI_Allgather(&nforces, 1, MPI_INT, Recv_count, 1, MPI_INT, D->Communicator);

  int nimport    = 0;
  Recv_offset[0] = 0;

  for(int j = 0; j < D->NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
    }

  DirectDataAll = (directdata *)Mem.mymalloc("DirectDataAll", nimport * sizeof(directdata));

  for(int j = 0; j < D->NTask; j++)
    {
      Send_count[j]  = Recv_count[j] * sizeof(directdata);
      Send_offset[j] = Recv_offset[j] * sizeof(directdata);
    }

  myMPI_Allgatherv(DirectDataIn, nforces * sizeof(directdata), MPI_BYTE, DirectDataAll, Send_count, Send_offset, MPI_BYTE,
                   D->Communicator);

  /* subdivide the work evenly */
  int first, count;
  subdivide_evenly(nimport, D->NTask, D->ThisTask, &first, &count);

  DirectAccOut = (accdata *)Mem.mymalloc("DirectDataOut", count * sizeof(accdata));

  /* now calculate the forces */
  for(int i = 0; i < count; i++)
    {
      int target     = i + first;
      int result_idx = i;

      vector<double> acc = 0.0;
#ifdef EVALPOTENTIAL
      double pot = 0.0;
#endif

#if NSOFTCLASSES > 1
      double h_i = All.ForceSoftening[DirectDataAll[target].SofteningClass];
#else
      double h_i = All.ForceSoftening[0];
#endif

      for(int j = 0; j < nimport; j++)
        {
#if NSOFTCLASSES > 1
          double h_j = All.ForceSoftening[DirectDataAll[j].SofteningClass];
#else
          double h_j = All.ForceSoftening[0];
#endif
          double hmax = (h_j > h_i) ? h_j : h_i;

          vector<double> dxyz;
          Sp->nearest_image_intpos_to_pos(DirectDataAll[j].IntPos, DirectDataAll[target].IntPos,
                                          dxyz.da); /* converts the integer distance to floating point */

          double r2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];

          double mass = DirectDataAll[j].Mass;

          /* now evaluate the force component */

          double r = sqrt(r2);

          double rinv = (r > 0) ? 1.0 / r : 0;

          gravtree<simparticles>::gfactors gfac;

#ifdef PMGRID
          mesh_factors *mfp = &mf[LOW_MESH];
#if defined(PLACEHIGHRESREGION)
          if((DoPM & TREE_ACTIVE_CUTTOFF_HIGHRES_PM))
            {
              if(DirectDataAll[j].InsideOutsideFlag == FLAG_INSIDE && DirectDataAll[target].InsideOutsideFlag == FLAG_INSIDE)
                mfp = &mf[HIGH_MESH];
            }
#endif
          if((DoPM & (TREE_ACTIVE_CUTTOFF_BASE_PM + TREE_ACTIVE_CUTTOFF_HIGHRES_PM)))
            {
              if(modify_gfactors_pm_monopole(gfac, r, rinv, mfp))
                return;  // if we are outside the cut-off radius, we have no interaction
            }
#endif
          get_gfactors_monopole(gfac, r, hmax, rinv);

#ifdef EVALPOTENTIAL
          pot -= mass * gfac.fac0;
#endif
          acc -= (mass * gfac.fac1 * rinv) * dxyz;

          if(DoEwald)
            {
              // EWALD treatment, only done for periodic boundaries in case PM is not active

              ewald_data ew;
              Ewald.ewald_gridlookup(DirectDataAll[j].IntPos, DirectDataAll[target].IntPos, ewald::POINTMASS, ew);

#ifdef EVALPOTENTIAL
              pot += mass * ew.D0phi;
#endif
              acc += mass * ew.D1phi;
            }
        }

      DirectAccOut[result_idx].Acc[0] = acc[0];
      DirectAccOut[result_idx].Acc[1] = acc[1];
      DirectAccOut[result_idx].Acc[2] = acc[2];
#ifdef EVALPOTENTIAL
      DirectAccOut[result_idx].Potential = pot;
#endif
    }

  /* now send the forces to the right places */

  DirectAccIn = (accdata *)Mem.mymalloc("DirectDataIn", nforces * sizeof(accdata));

  MPI_Request *requests = (MPI_Request *)Mem.mymalloc_movable(&requests, "requests", 2 * D->NTask * sizeof(MPI_Request));
  int n_requests        = 0;

  int recvTask = 0;
  int sendTask = 0;
  int send_first, send_count;
  subdivide_evenly(nimport, D->NTask, sendTask, &send_first, &send_count);

  while(recvTask < D->NTask && sendTask < D->NTask) /* go through both lists */
    {
      while(send_first + send_count < Recv_offset[recvTask])
        {
          if(sendTask >= D->NTask - 1)
            Terminate("sendTask >= NTask  recvTask=%d sendTask=%d", recvTask, sendTask);

          sendTask++;
          subdivide_evenly(nimport, D->NTask, sendTask, &send_first, &send_count);
        }

      while(Recv_offset[recvTask] + Recv_count[recvTask] < send_first)
        {
          if(recvTask >= D->NTask - 1)
            Terminate("recvTask >= NTask  recvTask=%d sendTask=%d", recvTask, sendTask);

          recvTask++;
        }

      int start = std::max<int>(Recv_offset[recvTask], send_first);
      int next  = std::min<int>(Recv_offset[recvTask] + Recv_count[recvTask], send_first + send_count);

      if(next - start >= 1)
        {
          if(D->ThisTask == sendTask)
            MPI_Isend(DirectAccOut + start - send_first, (next - start) * sizeof(accdata), MPI_BYTE, recvTask, TAG_PDATA_SPH,
                      D->Communicator, &requests[n_requests++]);

          if(D->ThisTask == recvTask)
            MPI_Irecv(DirectAccIn + start - Recv_offset[recvTask], (next - start) * sizeof(accdata), MPI_BYTE, sendTask, TAG_PDATA_SPH,
                      D->Communicator, &requests[n_requests++]);
        }

      if(next == Recv_offset[recvTask] + Recv_count[recvTask])
        recvTask++;
      else
        {
          sendTask++;
          if(sendTask >= D->NTask)
            break;

          subdivide_evenly(nimport, D->NTask, sendTask, &send_first, &send_count);
        }
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
  Mem.myfree(requests);

  nforces = 0;

  for(int idx = 0; idx < Sp->TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = Sp->TimeBinsGravity.ActiveParticleList[idx];

      for(int k = 0; k < 3; k++)
        Sp->P[i].GravAccel[k] = DirectAccIn[nforces].Acc[k];

#ifdef EVALPOTENTIAL
      Sp->P[i].Potential = DirectAccIn[nforces].Potential;
#endif
      nforces++;
    }

  Mem.myfree(DirectAccIn);
  Mem.myfree(DirectAccOut);
  Mem.myfree(DirectDataAll);
  Mem.myfree(DirectDataIn);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  D->mpi_printf("GRAVDIRECT: force is done.\n");

  All.TotNumDirectForces += Sp->TimeBinsGravity.GlobalNActiveParticles;

  double tend = Logs.second();

  double timedirect, sumt;
  timedirect = tend - tstart;

  MPI_Reduce(&timedirect, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, D->Communicator);

  if(D->ThisTask == 0)
    {
      fprintf(Logs.FdTimings, "Nf=%9lld  timebin=%d  active part/task: avg=%g   total-Nf=%lld\n",
              Sp->TimeBinsGravity.GlobalNActiveParticles, timebin, ((double)Sp->TimeBinsGravity.GlobalNActiveParticles) / D->NTask,
              All.TotNumDirectForces);
      fprintf(Logs.FdTimings, "  (direct) took=%g sec part/sec:  %g   ia/sec: %g\n", timedirect,
              Sp->TimeBinsGravity.GlobalNActiveParticles / (sumt + 1.0e-20),
              Sp->TimeBinsGravity.GlobalNActiveParticles / (sumt + 1.0e-20) * Sp->TimeBinsGravity.GlobalNActiveParticles);
      myflush(Logs.FdTimings);
    }

  TIMER_STOP(CPU_TREEDIRECT);
}

#endif
