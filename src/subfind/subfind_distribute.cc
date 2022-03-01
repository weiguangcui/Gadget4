/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_distribute.cc
 *
 *  \brief code to distribute particle data for Subfind processing
 */

#include "gadgetconfig.h"

#ifdef SUBFIND

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fof/fof.h"
#include "../gravtree/gravtree.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/cxxsort.h"
#include "../subfind/subfind.h"

template <typename partset>
void fof<partset>::subfind_distribute_groups(void)
{
  double t0 = Logs.second();

  int *Send_count  = (int *)Mem.mymalloc_movable(&Send_count, "Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc_movable(&Send_offset, "Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc_movable(&Recv_count, "Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc_movable(&Recv_offset, "Recv_offset", sizeof(int) * NTask);

  /* count how many we have of each task */
  for(int i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(int i = 0; i < Ngroups; i++)
    {
      int target = Group[i].TargetTask;

      if(target < 0 || target >= NTask)
        Terminate("target < 0 || target >= NTask");

      if(target != ThisTask)
        Send_count[target]++;
    }

  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

  Recv_offset[0] = Send_offset[0] = 0;
  int nexport = 0, nimport = 0;

  for(int i = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  group_properties *send_Group =
      (group_properties *)Mem.mymalloc_movable(&send_Group, "send_Group", nexport * sizeof(group_properties));

  for(int i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(int i = 0; i < Ngroups; i++)
    {
      int target = Group[i].TargetTask;

      if(target != ThisTask)
        {
          send_Group[Send_offset[target] + Send_count[target]] = Group[i];
          Send_count[target]++;

          Group[i] = Group[Ngroups - 1];
          Ngroups--;
          i--;
        }
    }

  if(Ngroups + nimport > MaxNgroups)
    {
      MaxNgroups = Ngroups + nimport;
      Group      = (group_properties *)Mem.myrealloc_movable(Group, sizeof(group_properties) * MaxNgroups);
    }

  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              myMPI_Sendrecv(&send_Group[Send_offset[recvTask]], Send_count[recvTask] * sizeof(group_properties), MPI_BYTE, recvTask,
                           TAG_DENS_A, &Group[Ngroups + Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(group_properties),
                           MPI_BYTE, recvTask, TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);
            }
        }
    }

  Ngroups += nimport;

  Mem.myfree_movable(send_Group);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  double t1 = Logs.second();

  mpi_printf("SUBFIND: subfind_distribute_groups() took %g sec\n", Logs.timediff(t0, t1));
}

/* This function redistributes the particles in P[] and PS[] according to what is stored in
 * PS[].TargetTask, and PS[].TargetIndex. NOTE: The associated SphP[] is not moved, i.e. the
 * association is broken until the particles are moved back into the original order.
 */
template <typename partset>
void fof<partset>::subfind_distribute_particles(MPI_Comm Communicator)
{
  int *Send_count  = (int *)Mem.mymalloc_movable(&Send_count, "Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc_movable(&Send_offset, "Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc_movable(&Recv_count, "Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc_movable(&Recv_offset, "Recv_offset", sizeof(int) * NTask);

  int CommThisTask, CommNTask;
  MPI_Comm_size(Communicator, &CommNTask);
  MPI_Comm_rank(Communicator, &CommThisTask);

  for(int n = 0; n < CommNTask; n++)
    Send_count[n] = 0;

  for(int n = 0; n < Tp->NumPart; n++)
    {
      int target = Tp->PS[n].TargetTask;

      if(target != CommThisTask)
        {
          if(target < 0 || target >= CommNTask)
            Terminate("n=%d targettask=%d", n, target);

          Send_count[target]++;
        }
    }

  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

  int nimport = 0, nexport = 0;
  Recv_offset[0] = 0, Send_offset[0] = 0;

  for(int j = 0; j < CommNTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* for resize */
  int load = (Tp->NumPart + (nimport - nexport)), max_load;
  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, Communicator);

  typename partset::pdata *partBuf =
      (typename partset::pdata *)Mem.mymalloc_movable(&partBuf, "partBuf", nexport * sizeof(typename partset::pdata));
  subfind_data *subBuf = (subfind_data *)Mem.mymalloc_movable(&subBuf, "subBuf", nexport * sizeof(subfind_data));

  for(int i = 0; i < CommNTask; i++)
    Send_count[i] = 0;

  for(int n = 0; n < Tp->NumPart; n++)
    {
      int target = Tp->PS[n].TargetTask;

      if(target != CommThisTask)
        {
          partBuf[Send_offset[target] + Send_count[target]] = Tp->P[n];
          subBuf[Send_offset[target] + Send_count[target]]  = Tp->PS[n];

          Tp->P[n]  = Tp->P[Tp->NumPart - 1];
          Tp->PS[n] = Tp->PS[Tp->NumPart - 1];

          Send_count[target]++;
          Tp->NumPart--;
          n--;
        }
    }

  /* do resize */
  if(max_load > (1.0 - ALLOC_TOLERANCE) * Tp->MaxPart || max_load < (1.0 - 3 * ALLOC_TOLERANCE) * Tp->MaxPart)
    Tp->reallocate_memory_maxpart(max_load / (1.0 - 2 * ALLOC_TOLERANCE));

  Tp->PS = (subfind_data *)Mem.myrealloc_movable(Tp->PS, load * sizeof(subfind_data));

  for(int i = 0; i < CommNTask; i++)
    Recv_offset[i] += Tp->NumPart;

#ifdef ISEND_IRECV_IN_DOMAIN

  MPI_Request *requests = (MPI_Request *)Mem.mymalloc("requests", 8 * CommNTask * sizeof(MPI_Request));
  int n_requests        = 0;

  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Recv_count[target] > 0)
            {
              MPI_Irecv(Tp->P + Recv_offset[target], Recv_count[target] * sizeof(typename partset::pdata), MPI_BYTE, target, TAG_PDATA,
                        Communicator, &requests[n_requests++]);
              MPI_Irecv(Tp->PS + Recv_offset[target], Recv_count[target] * sizeof(subfind_data), MPI_BYTE, target, TAG_KEY,
                        Communicator, &requests[n_requests++]);
            }
        }
    }

  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Send_count[target] > 0)
            {
              MPI_Issend(partBuf + Send_offset[target], Send_count[target] * sizeof(typename partset::pdata), MPI_BYTE, target,
                         TAG_PDATA, Communicator, &requests[n_requests++]);
              MPI_Issend(subBuf + Send_offset[target], Send_count[target] * sizeof(subfind_data), MPI_BYTE, target, TAG_KEY,
                         Communicator, &requests[n_requests++]);
            }
        }
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
  Mem.myfree(requests);

#else
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Send_count[target] > 0 || Recv_count[target] > 0)
            {
              myMPI_Sendrecv(partBuf + Send_offset[target], Send_count[target] * sizeof(typename partset::pdata), MPI_BYTE, target,
                             TAG_PDATA, Tp->P + Recv_offset[target], Recv_count[target] * sizeof(typename partset::pdata), MPI_BYTE,
                             target, TAG_PDATA, Communicator, MPI_STATUS_IGNORE);

              myMPI_Sendrecv(subBuf + Send_offset[target], Send_count[target] * sizeof(subfind_data), MPI_BYTE, target, TAG_KEY,
                             Tp->PS + Recv_offset[target], Recv_count[target] * sizeof(subfind_data), MPI_BYTE, target, TAG_KEY,
                             Communicator, MPI_STATUS_IGNORE);
            }
        }
    }
#endif

  Tp->NumPart += nimport;
  Mem.myfree_movable(subBuf);
  Mem.myfree_movable(partBuf);

  /* finally, let's also address the desired local order according to PS[].TargetIndex */
  FoFDomain->reorder_P_PS(0, Tp->NumPart);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
