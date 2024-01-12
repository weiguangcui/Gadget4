/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  rearrange.cc
 *
 *  \brief routines to rearrange some of the lightcone data on disk to allow easier processing later on
 */

#include "gadgetconfig.h"

#ifdef MERGERTREE

#include <gsl/gsl_rng.h>
#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../io/hdf5_util.h"
#include "../io/snap_io.h"
#include "../lightcone/lightcone.h"
#include "../lightcone/lightcone_particle_io.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mergertree/io_readtrees_mbound.h"
#include "../sort/cxxsort.h"
#include "../sort/parallel_sort.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"

#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
void sim::rearrange_snapshot(int argc, char **argv)
{
  if(argc < 5)
    Terminate("too few arguments: Snapshot rearrange requires  <firstnum>  <lastnum>");

  int conenr   = 0;  // not needed here
  int firstnum = atoi(argv[3]);
  int lastnum  = atoi(argv[4]);

  All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
  Driftfac.init_drift_table();

  rearrange_generic<simparticles>(Sp, conenr, firstnum, lastnum);
}

/* specialization */
template <>
void sim::rearrange_read<simparticles>(simparticles &Tp, int num, int conenr)
{
  /* read the lightcone  */
  snap_io Snap{&Tp, Communicator, All.SnapFormat}; /* get an I/O object */
  Snap.read_snapshot(num, MOST_BOUND_PARTICLE_SNAPHOT);
}

/* specialization */
template <>
void sim::rearrange_write<simparticles>(simparticles &Tp, int num, int conenr)
{
  /* write the snapshot  */
  snap_io Snap{&Tp, &MergerTree, Communicator, All.SnapFormat};     /* get an I/O object */
  Snap.write_snapshot(num, MOST_BOUND_PARTICLE_SNAPHOT_REORDERED);  // true?
}
#endif

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES) && defined(REARRANGE_OPTION)

void sim::rearrange_lightcone(int argc, char **argv)
{
  if(argc < 6)
    Terminate("too few arguments: Lightcone rearrange requires  <conenr>  <firstnum>  <lastnum>");

  int conenr   = atoi(argv[3]);
  int firstnum = atoi(argv[4]);
  int lastnum  = atoi(argv[5]);

  All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
  Driftfac.init_drift_table();

  if(LightCone.lightcone_init_times())
    endrun();

#ifdef LIGHTCONE_MASSMAPS
  LightCone.lightcone_init_massmaps();
  if(LightCone.lightcone_massmap_report_boundaries())
    endrun();
#endif

  if(LightCone.lightcone_init_boxlist())
    endrun();


  double linklength = 0;
  LightCone.lightcone_init_intposconverter(linklength);

  rearrange_generic<lcparticles>(Lp, conenr, firstnum, lastnum);
}

/* specialization */
template <>
void sim::rearrange_read<lcparticles>(lcparticles &Tp, int num, int conenr)
{
  /* read the lightcone  */
  lightcone_particle_io LcIO{&Tp, &LightCone, &MergerTree, Communicator, All.SnapFormat}; /* get an I/O object */
  LcIO.lightcone_read(num, conenr);
}

/* specialization */
template <>
void sim::rearrange_write<lcparticles>(lcparticles &Tp, int num, int conenr)
{
  /* write the lightcone  */
  lightcone_particle_io LcIO{&Tp, &LightCone, &MergerTree, Communicator, All.SnapFormat}; /* get an I/O object */
  LcIO.lightcone_save(num, conenr, true);
}

#endif

template <typename partset>
void sim::rearrange_generic(partset &Tp, int conenr, int firstnum, int lastnum)
{
  /* read the merger tree mostbound halo IDs  */
  {
    readtrees_mbound_io MtreeIO{&MergerTree, Communicator, All.SnapFormat}; /* get an I/O object */
    MtreeIO.read_trees_mostbound();                                         // optimize this here to only read most-bound ID
  }

  /* let's now sort the tree data according to the IDs */
  mycxxsort_parallel(MergerTree.HaloIDdata, MergerTree.HaloIDdata + MergerTree.Nhalos, MergerTree.compare_HaloIDdata_ID, Communicator);

  /* establish some lists with the minimum and maximum IDs stored on each processor */
  long long halominID = (MergerTree.Nhalos > 0) ? MergerTree.HaloIDdata[0].SubMostBoundID : 0;
  long long halomaxID = (MergerTree.Nhalos > 0) ? MergerTree.HaloIDdata[MergerTree.Nhalos - 1].SubMostBoundID : 0;

  long long *list_halominID = (long long *)Mem.mymalloc("list_halominID", NTask * sizeof(long long));
  long long *list_halomaxID = (long long *)Mem.mymalloc("list_halomaxID", NTask * sizeof(long long));
  int *list_Nhalos          = (int *)Mem.mymalloc("list_Nhalos", NTask * sizeof(int));

  MPI_Allgather(&halominID, 1, MPI_LONG_LONG, list_halominID, 1, MPI_LONG_LONG, Communicator);
  MPI_Allgather(&halomaxID, 1, MPI_LONG_LONG, list_halomaxID, 1, MPI_LONG_LONG, Communicator);
  MPI_Allgather(&MergerTree.Nhalos, 1, MPI_INT, list_Nhalos, 1, MPI_INT, Communicator);

  typedef mergertree::treehalo_ids_type treehalo_ids_type;

  for(int num = firstnum; num <= lastnum; num++)
    {
      rearrange_read(Tp, num, conenr);

      mpi_printf("REARRANGE: On Task=%d, %d particles\n", ThisTask, Tp.NumPart);

      /* let's now sort the lightcone_particle_data according to ID */
      mycxxsort_parallel(Tp.P, Tp.P + Tp.NumPart, Tp.compare_ID, Communicator);

      /* establish some lists with the minimum and maximum IDs stored on each processor for the lightcone particles */
      long long coneminID = (Tp.NumPart > 0) ? Tp.P[0].ID.get() : 0;
      long long conemaxID = (Tp.NumPart > 0) ? Tp.P[Tp.NumPart - 1].ID.get() : 0;

      long long *list_coneminID = (long long *)Mem.mymalloc("list_coneminID", NTask * sizeof(long long));
      long long *list_conemaxID = (long long *)Mem.mymalloc("list_conemaxID", NTask * sizeof(long long));
      int *list_NumPart         = (int *)Mem.mymalloc("list_NumPart", NTask * sizeof(int));

      MPI_Allgather(&coneminID, 1, MPI_LONG_LONG, list_coneminID, 1, MPI_LONG_LONG, Communicator);
      MPI_Allgather(&conemaxID, 1, MPI_LONG_LONG, list_conemaxID, 1, MPI_LONG_LONG, Communicator);
      MPI_Allgather(&Tp.NumPart, 1, MPI_INT, list_NumPart, 1, MPI_INT, Communicator);

      // let's assign matchIDs to the degree possible and assign TreeIDs in accordance to that

      // default is no matching Tree
      for(int n = 0; n < Tp.NumPart; n++)
        Tp.P[n].TreeID = -1;

      MPI_Request *requests = (MPI_Request *)Mem.mymalloc("requests", NTask * sizeof(MPI_Request));

      int nreq = 0;

      for(int task = 0; task < NTask; task++)
        {
          if(MergerTree.Nhalos > 0 && list_NumPart[task] > 0)
            if(!(halomaxID < list_coneminID[task] || halominID > list_conemaxID[task]))
              {
                MPI_Issend(MergerTree.HaloIDdata, MergerTree.Nhalos * sizeof(mergertree::treehalo_ids_type), MPI_BYTE, task, TAG_N,
                           Communicator, &requests[nreq++]);
              }
        }

      for(int task = 0; task < NTask; task++)
        {
          if(list_Nhalos[task] > 0 && Tp.NumPart > 0)
            if(!(list_halomaxID[task] < coneminID || list_halominID[task] > conemaxID))
              {
                treehalo_ids_type *halos = (treehalo_ids_type *)Mem.mymalloc("halos", list_Nhalos[task] * sizeof(treehalo_ids_type));

                MPI_Recv(halos, list_Nhalos[task] * sizeof(treehalo_ids_type), MPI_BYTE, task, TAG_N, Communicator, MPI_STATUS_IGNORE);

                int i_halo = 0;
                int i_cone = 0;

                while(i_halo < list_Nhalos[task] && i_cone < Tp.NumPart)
                  {
                    if(halos[i_halo].SubMostBoundID == Tp.P[i_cone].ID.get())
                      {
                        Tp.P[i_cone++].TreeID = halos[i_halo++].TreeID;
                      }
                    else if(halos[i_halo].SubMostBoundID < Tp.P[i_cone].ID.get())
                      i_halo++;
                    else
                      i_cone++;
                  }

                Mem.myfree(halos);
              }
        }

      if(nreq)
        MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);

      Mem.myfree(requests);
      Mem.myfree(list_NumPart);
      Mem.myfree(list_conemaxID);
      Mem.myfree(list_coneminID);

      /* now we sort the lightcone_particle_data according to TreeID */
      mycxxsort_parallel(Tp.P, Tp.P + Tp.NumPart, Tp.compare_TreeID_ID, Communicator);

      rearrange_fill_treetable<partset>(Tp);

      /* write the lightcone  */
      rearrange_write(Tp, num, conenr);

      /* free the storage again */
      Tp.free_memory();
    }

  Mem.myfree(list_Nhalos);
  Mem.myfree(list_halomaxID);
  Mem.myfree(list_halominID);
}

template <typename partset>
void sim::rearrange_fill_treetable(partset &Tp)
{
  /*
     we use HaloCount just like ParticleCount, and FirstHalo becomes "ParticleFirst */

  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);

  for(int i = 0; i < NTask; i++)
    Send_count[i] = 0;

  long long *TreeID_list = (long long *)Mem.mymalloc("TreeID_list", sizeof(long long) * Tp.NumPart);

  long long maxTreeID = (MergerTree.Ntrees > 0) ? MergerTree.TreeTable[MergerTree.Ntrees - 1].TreeID : -1;

  long long *maxTreeID_list = (long long *)Mem.mymalloc("maxTreeID_list", sizeof(long long) * NTask);
  MPI_Allgather(&maxTreeID, sizeof(long long), MPI_BYTE, maxTreeID_list, sizeof(long long), MPI_BYTE, Communicator);

  int target_task = 0;

  for(int i = 0; i < Tp.NumPart; i++)
    {
      TreeID_list[i] = Tp.P[i].TreeID;

      while(target_task < NTask - 1 && TreeID_list[i] > maxTreeID_list[target_task])
        target_task++;

      Send_count[target_task]++;
    }

  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

  Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      if(j > 0)
        Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
    }

  for(int i = 0; i < MergerTree.Ntrees; i++)
    MergerTree.TreeTable[i].HaloCount = 0;

  /* exchange data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              long long *treeid_tmp = (long long *)Mem.mymalloc("treeid_tmp", sizeof(long long) * Recv_count[recvTask]);

              myMPI_Sendrecv(&TreeID_list[Send_offset[recvTask]], Send_count[recvTask] * sizeof(long long), MPI_BYTE, recvTask,
                             TAG_DENS_A, treeid_tmp, Recv_count[recvTask] * sizeof(long long), MPI_BYTE, recvTask, TAG_DENS_A,
                             Communicator, MPI_STATUS_IGNORE);

              for(int i = 0; i < Recv_count[recvTask]; i++)
                {
                  if(treeid_tmp[i] != -1)
                    {
                      int off = treeid_tmp[i] - MergerTree.TreeTable[0].TreeID;

                      if(off < 0 || off >= MergerTree.Ntrees)
                        Terminate(
                            "strange: off=%d  MergerTree.Ntrees=%d  treeid_tmp[i]=%lld  MergerTree.TreeTable[0].TreeID=%lld  i=%d  "
                            "Recv_count[recvTask]=%d  ",
                            off, MergerTree.Ntrees, treeid_tmp[i], MergerTree.TreeTable[0].TreeID, i, Recv_count[recvTask]);

                      MergerTree.TreeTable[off].HaloCount++;
                    }
                }

              Mem.myfree(treeid_tmp);
            }
        }
    }

  Mem.myfree(maxTreeID_list);
  Mem.myfree(TreeID_list);

  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  /* now also fix the cumulative count */

  if(MergerTree.Ntrees > 0)
    MergerTree.TreeTable[0].FirstHalo = 0;

  for(int i = 1; i < MergerTree.Ntrees; i++)
    MergerTree.TreeTable[i].FirstHalo = MergerTree.TreeTable[i - 1].FirstHalo + MergerTree.TreeTable[i - 1].HaloCount;

  long long cumul = 0;
  for(int i = 0; i < MergerTree.Ntrees; i++)
    cumul += MergerTree.TreeTable[i].HaloCount;

  long long *cumul_list = (long long *)Mem.mymalloc("cumul_list", sizeof(long long) * NTask);

  MPI_Allgather(&cumul, sizeof(long long), MPI_BYTE, cumul_list, sizeof(long long), MPI_BYTE, Communicator);

  cumul = 0;
  for(int i = 0; i < ThisTask; i++)
    cumul += cumul_list[i];

  for(int i = 0; i < MergerTree.Ntrees; i++)
    MergerTree.TreeTable[i].FirstHalo += cumul;

  Mem.myfree(cumul_list);
}

#endif
