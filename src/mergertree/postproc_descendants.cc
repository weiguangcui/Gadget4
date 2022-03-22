/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file postproc_descendants.cc
 *
 *  \brief code to determine the descendant subhalo of all subhalos in one catalogue in another one
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
#include "../data/mymalloc.h"
#include "../fof/fof.h"
#include "../fof/fof_io.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mergertree/io_readsnap.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

// this is used in FOF postprocessing to assign the previous subhalo length to particles
void mergertree::get_previous_size_of_subhalo_for_each_particle(int num)
{
  if(num >= 0)
    {
      mpi_printf(
          "SUBFIND / MERGERTREE: We are loading the previous group catalogue to set the size of previous subhalos for each "
          "particle\n");

      domain<simparticles> Domain{Communicator, Sp};
      fof<simparticles> FoF{Communicator, Sp, &Domain};

      /* load previous snapshot and group/subhalo catalogues */
      fof_io<simparticles> FoF_IO{&FoF, Communicator, All.SnapFormat};
      FoF_IO.fof_subfind_load_groups(num);

      readsnap_io Snap_IO{this, this->Communicator, All.SnapFormat};
      Snap_IO.mergertree_read_snap_ids(num);

      /* make sure that group catalog and snapshot are in file order */
      mycxxsort_parallel(MtrP, MtrP + MtrP_NumPart, compare_MtrP_FileOffset, Communicator);
      mycxxsort_parallel(FoF.Group, FoF.Group + FoF.Ngroups, compare_Group_FileOffset, Communicator);
      mycxxsort_parallel(FoF.Subhalo, FoF.Subhalo + FoF.Nsubhalos, compare_Subhalo_FileOffset, Communicator);

      /* now assign (global) group and subhalo numbers to each particle belonging to the particular structure */
      mergertree_assign_group_numbers(&FoF);

      Mem.myfree_movable(FoF.Subhalo);
      Mem.myfree_movable(FoF.Group);

      /* get PrevSizeOfSubhalo from PrevMtrP by matching IDs */
      mergertree_match_ids_of_current_snap();
    }
  else
    {
      for(int i = 0; i < Sp->NumPart; i++)
        {
          Sp->P[i].PrevSubhaloNr.set(HALONR_MAX);
          Sp->P[i].PrevSizeOfSubhalo.set(0);
        }
    }
}

void mergertree::descendants_in_postprocessing(int num)
{
  if(num - 1 < 0)
    Terminate("cannot execute for the given snapnum");

  domain<simparticles> Domain{Communicator, Sp};
  fof<simparticles> FoF{Communicator, Sp, &Domain};

  /* load previous snapshot and group/subhalo catalogues */

  fof_io<simparticles> FoF_IO{&FoF, Communicator, All.SnapFormat};
  FoF_IO.fof_subfind_load_groups(num - 1);

  readsnap_io Snap_IO{this, this->Communicator, All.SnapFormat};
  Snap_IO.mergertree_read_snap_ids(num - 1);

  /* make sure that group catalog and snapshot are in file order */
  mycxxsort_parallel(MtrP, MtrP + MtrP_NumPart, compare_MtrP_FileOffset, Communicator);
  mycxxsort_parallel(FoF.Group, FoF.Group + FoF.Ngroups, compare_Group_FileOffset, Communicator);
  mycxxsort_parallel(FoF.Subhalo, FoF.Subhalo + FoF.Nsubhalos, compare_Subhalo_FileOffset, Communicator);

  if(FoF_IO.LegacyFormat)
    {
      mpi_printf("\nFOF/SUBFIND: Legacy format from Arepo detected, trying to adjust for this.\n");
      FoF.subfind_redetermine_groupnr();
    }

  /* now assign (global) group and subhalo numbers to each particle belonging to the particular structure */
  mergertree_assign_group_numbers(&FoF);

  Mem.myfree_movable(FoF.Subhalo);
  Mem.myfree_movable(FoF.Group);

  /* save previous data */
  PrevTotNsubhalos = FoF.TotNsubhalos;
  PrevNsubhalos    = FoF.Nsubhalos;

  PrevMtrP_NumPart = MtrP_NumPart;

  PrevMtrP =
      (mergertree_particle_data *)Mem.mymalloc_movable(&PrevMtrP, "PrevMtrP", PrevMtrP_NumPart * sizeof(mergertree_particle_data));
  memcpy(PrevMtrP, MtrP, PrevMtrP_NumPart * sizeof(mergertree_particle_data));
  Mem.myfree_movable(MtrP);

  /* load new snapshot and group/subhalo catalogues */
  FoF_IO.fof_subfind_load_groups(num);
  CurrTotNsubhalos = FoF.TotNsubhalos;
  CurrNsubhalos    = FoF.Nsubhalos;

  Snap_IO.mergertree_read_snap_ids(num);

  /* make sure that group catalog and snapshot are in file order */
  mycxxsort_parallel(MtrP, MtrP + MtrP_NumPart, compare_MtrP_FileOffset, Communicator);
  mycxxsort_parallel(FoF.Group, FoF.Group + FoF.Ngroups, compare_Group_FileOffset, Communicator);
  mycxxsort_parallel(FoF.Subhalo, FoF.Subhalo + FoF.Nsubhalos, compare_Subhalo_FileOffset, Communicator);

  if(FoF_IO.LegacyFormat)
    {
      mpi_printf("\nFOF/SUBFIND: Legacy format from Arepo detected, trying to adjust for this.\n");
      FoF.subfind_redetermine_groupnr();
    }

  /* now assign (global) group and subhalo numbers to each particle belonging to the particular structure */
  mergertree_assign_group_numbers(&FoF);

  Mem.myfree_movable(FoF.Subhalo);
  Mem.myfree_movable(FoF.Group);

  /* assign the determined new subhalonr */
  for(int i = 0; i < MtrP_NumPart; i++)
    {
      MtrP[i].SubhaloNr         = MtrP[i].PrevSubhaloNr;
      MtrP[i].SubhaloLen        = MtrP[i].PrevSubhaloLen;
      MtrP[i].GroupNr           = MtrP[i].PrevGroupNr;
      MtrP[i].RankInSubhalo     = MtrP[i].PrevRankInSubhalo;
      MtrP[i].PrevSubhaloNr     = HALONR_MAX;
      MtrP[i].PrevGroupNr       = HALONR_MAX;
      MtrP[i].PrevSubhaloLen    = 0;
      MtrP[i].PrevRankInSubhalo = 0;
    }

  /* get PrevSubhaloNr from PrevMtrP by matching IDs */
  mergertree_match_ids_of_previous_snap();

  mergertree_determine_descendants_postproc(num);
}

void mergertree::mergertree_match_ids_of_previous_snap(void)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  mycxxsort_parallel(MtrP, MtrP + MtrP_NumPart, compare_MtrP_ID, Communicator);
  mycxxsort_parallel(PrevMtrP, PrevMtrP + PrevMtrP_NumPart, compare_MtrP_ID, Communicator);

  MyIDType *list_min_id = (MyIDType *)Mem.mymalloc("list_min_id", NTask * sizeof(MyIDType));
  MyIDType *list_max_id = (MyIDType *)Mem.mymalloc("list_max_id", NTask * sizeof(MyIDType));
  int *list_numpart     = (int *)Mem.mymalloc("list_mumpart", NTask * sizeof(int));

  MPI_Allgather(&MtrP[0].ID, sizeof(MyIDType), MPI_BYTE, list_min_id, sizeof(MyIDType), MPI_BYTE, Communicator);
  MPI_Allgather(&MtrP[MtrP_NumPart > 0 ? MtrP_NumPart - 1 : 0].ID, sizeof(MyIDType), MPI_BYTE, list_max_id, sizeof(MyIDType), MPI_BYTE,
                Communicator);
  MPI_Allgather(&MtrP_NumPart, sizeof(int), MPI_BYTE, list_numpart, sizeof(int), MPI_BYTE, Communicator);

  int nexport = 0, nimport = 0;
  mergertree_particle_data *import_data = NULL, *export_data = NULL;

  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target = 0;

      for(int i = 0; i < PrevMtrP_NumPart; i++)
        {
          while(target < NTask - 1 && (list_numpart[target] == 0 || PrevMtrP[i].ID > list_max_id[target]))
            target++;

          if(list_numpart[target] != 0)
            if(PrevMtrP[i].ID >= list_min_id[target] && PrevMtrP[i].ID <= list_max_id[target])
              {
                if(mode == 0)
                  Send_count[target]++;
                else
                  {
                    int off          = Send_offset[target] + Send_count[target]++;
                    export_data[off] = PrevMtrP[i];
                  }
              }
        }

      if(mode == 0)
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
          Recv_offset[0] = Send_offset[0] = 0;
          for(int j = 0; j < NTask; j++)
            {
              nimport += Recv_count[j];
              nexport += Send_count[j];
              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          export_data = (mergertree_particle_data *)Mem.mymalloc("export_data", nexport * sizeof(mergertree_particle_data));
          import_data = (mergertree_particle_data *)Mem.mymalloc("import_data", nimport * sizeof(mergertree_particle_data));
        }
    }

  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(mergertree_particle_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &import_data[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(mergertree_particle_data), MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                       MPI_STATUS_IGNORE);
    }

  /* incoming data should already be sorted, so now do the match */

  for(int i = 0, j = 0; i < MtrP_NumPart && j < nimport;)
    {
      if(MtrP[i].ID < import_data[j].ID)
        i++;
      else if(MtrP[i].ID > import_data[j].ID)
        j++;
      else
        {
          MtrP[i].PrevSubhaloNr     = import_data[j].PrevSubhaloNr;
          MtrP[i].PrevSubhaloLen    = import_data[j].PrevSubhaloLen;
          MtrP[i].PrevGroupNr       = import_data[j].PrevGroupNr;
          MtrP[i].PrevRankInSubhalo = import_data[j].PrevRankInSubhalo;
          i++;
          j++;
        }
    }

  Mem.myfree(import_data);
  Mem.myfree(export_data);
  Mem.myfree(list_numpart);
  Mem.myfree(list_max_id);
  Mem.myfree(list_min_id);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

void mergertree::mergertree_assign_group_numbers(fof<simparticles> *FoF)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* Tell everybody, how many particles are stored by each processor */

  FoF->fof_assign_group_offset();
  FoF->subfind_assign_subhalo_offsettype();

  int ntype_loc[NTYPES];            /* local particle number of each type */
  long long ntype_tot[NTYPES];      /* total particle numbers of each type */
  long long ntype_previous[NTYPES]; /* cumulative number of particles of each type on previous processors */

  for(int i = 0; i < NTYPES; i++)
    ntype_loc[i] = 0;

  for(int i = 0; i < MtrP_NumPart; i++)
    ntype_loc[MtrP[i].Type]++;

  /* collect a table with the particle numbers of each type on each processors */
  int *ntype_all = (int *)Mem.mymalloc("ntype_all", NTYPES * NTask * sizeof(int));
  MPI_Allgather(ntype_loc, NTYPES, MPI_INT, ntype_all, NTYPES, MPI_INT, Communicator);

  for(int i = 0; i < NTYPES; i++)
    ntype_tot[i] = 0;

  for(int i = 0; i < NTask; i++)
    for(int j = 0; j < NTYPES; j++)
      ntype_tot[j] += ntype_all[i * NTYPES + j];

  for(int i = 0; i < NTYPES; i++)
    ntype_previous[i] = 0;

  for(int i = 0; i < ThisTask; i++)
    for(int j = 0; j < NTYPES; j++)
      ntype_previous[j] += ntype_all[i * NTYPES + j];

  /* tell everybody how many groups each processor holds */

  int *gcount = (int *)Mem.mymalloc("gcount", NTask * sizeof(int));
  MPI_Allgather(&FoF->Ngroups, 1, MPI_INT, gcount, 1, MPI_INT, Communicator);

  /* determine the number of the first group we hold */
  long long first_groupnr = 0;
  for(int i = 0; i < ThisTask; i++)
    first_groupnr += gcount[i];

  /* tell everybody how many subhalos each processor holds */

  int *scount = (int *)Mem.mymalloc("scount", NTask * sizeof(int));
  MPI_Allgather(&FoF->Nsubhalos, 1, MPI_INT, scount, 1, MPI_INT, Communicator);

  /* determine the number of the first subhalo we hold */
  long long first_subhalonr = 0;
  for(int i = 0; i < ThisTask; i++)
    first_subhalonr += scount[i];

  /* let's now figure out which groups are needed by different processors to assign the group number information */

  struct group_info
  {
    long long GroupNr;
    long long First;
    MyLenType Len;
  };
  group_info *export_group_data = NULL, *import_group_data = NULL;

  for(int type = 0; type < NTYPES; type++)
    {
      int nexport = 0, nimport = 0;

      for(int mode = 0; mode < 2; mode++)
        {
          for(int i = 0; i < NTask; i++)
            Send_count[i] = 0;

          int target                = 0;
          long long first_in_target = 0; /* this the first particle (of this type) on the current target processor */

          for(int i = 0; i < FoF->Ngroups && target < NTask;)
            {
              int flag = 0;

              /* check whether we have an overlap */
              if(FoF->Group[i].OffsetType[type] + FoF->Group[i].LenType[type] > first_in_target &&
                 FoF->Group[i].OffsetType[type] < (first_in_target + ntype_all[target * NTYPES + type]))
                {
                  flag = 1;

                  if(mode == 0)
                    Send_count[target]++;
                  else
                    {
                      int off                        = Send_offset[target] + Send_count[target]++;
                      export_group_data[off].GroupNr = first_groupnr + i;
                      export_group_data[off].First   = FoF->Group[i].OffsetType[type];
                      export_group_data[off].Len     = FoF->Group[i].LenType[type];
                    }
                }

              if(FoF->Group[i].OffsetType[type] + FoF->Group[i].LenType[type] > first_in_target + ntype_all[target * NTYPES + type])
                {
                  first_in_target += ntype_all[target * NTYPES + type];
                  target++;
                }
              else
                {
                  if(i < FoF->Ngroups && flag == 0 && FoF->Group[i].LenType[type] > 0)
                    {
                      Terminate(
                          "strange: type=%d mode=%d  i=%d  first_in_target=%lld target=%d   FoF->Group[i].LenType[type]=%lld  "
                          "Ngroups=%d",
                          type, mode, i, first_in_target, target, (long long)FoF->Group[i].LenType[type], FoF->Ngroups);
                    }

                  i++;
                }
            }

          if(mode == 0)
            {
              myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
              Recv_offset[0] = Send_offset[0] = 0;
              for(int j = 0; j < NTask; j++)
                {
                  nimport += Recv_count[j];
                  nexport += Send_count[j];
                  if(j > 0)
                    {
                      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                    }
                }

              export_group_data = (group_info *)Mem.mymalloc("export_group_data", nexport * sizeof(group_info));
              import_group_data = (group_info *)Mem.mymalloc("import_group_data", nimport * sizeof(group_info));
            }
        }

      for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
        {
          int recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              myMPI_Sendrecv(&export_group_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(group_info), MPI_BYTE, recvTask,
                           TAG_DENS_B, &import_group_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(group_info), MPI_BYTE,
                           recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
        }

      /* now let's go through the local particles and assign the group numbers */

      int p  = ntype_previous[type];
      int gr = 0;

      for(int i = 0; i < MtrP_NumPart; i++)
        {
          if(MtrP[i].Type == type)
            {
              MtrP[i].PrevGroupNr = HALONR_MAX; /* default is not in a group */

              while(gr < nimport && p > import_group_data[gr].First + import_group_data[gr].Len - 1)
                gr++;

              if(gr < nimport && p >= import_group_data[gr].First && p < import_group_data[gr].First + import_group_data[gr].Len)
                MtrP[i].PrevGroupNr = import_group_data[gr].GroupNr;

              p++;
            }
        }

      Mem.myfree(import_group_data);
      Mem.myfree(export_group_data);
    }

  /* let's now figure out which groups are needed by different processors to assign the group number information */

  struct subhalo_info
  {
    long long SubhaloNr;
    long long First;
    MyLenType Len;
    MyLenType SubhaloLen;
  };
  subhalo_info *export_subhalo_data = NULL, *import_subhalo_data = NULL;

  for(int type = 0; type < NTYPES; type++)
    {
      int nexport = 0, nimport = 0;

      for(int mode = 0; mode < 2; mode++)
        {
          for(int i = 0; i < NTask; i++)
            Send_count[i] = 0;

          int target                = 0;
          long long first_in_target = 0; /* this the first particle (of this type) on the current target processor */

          for(int i = 0; i < FoF->Nsubhalos && target < NTask;)
            {
              /* check whether we have an overlap */
              if(FoF->Subhalo[i].OffsetType[type] + FoF->Subhalo[i].LenType[type] > first_in_target &&
                 FoF->Subhalo[i].OffsetType[type] < (first_in_target + ntype_all[target * NTYPES + type]))
                {
                  if(mode == 0)
                    Send_count[target]++;
                  else
                    {
                      int off                             = Send_offset[target] + Send_count[target]++;
                      export_subhalo_data[off].SubhaloNr  = first_subhalonr + i;
                      export_subhalo_data[off].First      = FoF->Subhalo[i].OffsetType[type];
                      export_subhalo_data[off].Len        = FoF->Subhalo[i].LenType[type];
                      export_subhalo_data[off].SubhaloLen = FoF->Subhalo[i].Len;
                    }
                }

              if(FoF->Subhalo[i].OffsetType[type] + FoF->Subhalo[i].LenType[type] >
                 first_in_target + ntype_all[target * NTYPES + type])
                {
                  first_in_target += ntype_all[target * NTYPES + type];
                  target++;
                }
              else
                i++;
            }

          if(mode == 0)
            {
              myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
              Recv_offset[0] = Send_offset[0] = 0;
              for(int j = 0; j < NTask; j++)
                {
                  nimport += Recv_count[j];
                  nexport += Send_count[j];
                  if(j > 0)
                    {
                      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                    }
                }

              export_subhalo_data = (subhalo_info *)Mem.mymalloc("export_subhalo_data", nexport * sizeof(subhalo_info));
              import_subhalo_data = (subhalo_info *)Mem.mymalloc("import_subhalo_data", nimport * sizeof(subhalo_info));
            }
        }

      for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
        {
          int recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              myMPI_Sendrecv(&export_subhalo_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(subhalo_info), MPI_BYTE,
                           recvTask, TAG_DENS_B, &import_subhalo_data[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(subhalo_info), MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                           MPI_STATUS_IGNORE);
        }

      /* now let's go through the local particles and assign group numbers */

      int p  = ntype_previous[type];
      int gr = 0;

      for(int i = 0; i < MtrP_NumPart; i++)
        {
          if(MtrP[i].Type == type)
            {
              MtrP[i].PrevSubhaloNr  = HALONR_MAX; /* default is not in a group */
              MtrP[i].PrevSubhaloLen = 0;

              while(gr < nimport && p > import_subhalo_data[gr].First + import_subhalo_data[gr].Len - 1)
                gr++;

              if(gr < nimport && p >= import_subhalo_data[gr].First && p < import_subhalo_data[gr].First + import_subhalo_data[gr].Len)
                {
                  MtrP[i].PrevSubhaloNr  = import_subhalo_data[gr].SubhaloNr;
                  MtrP[i].PrevSubhaloLen = import_subhalo_data[gr].SubhaloLen;

                  int rank = p - import_subhalo_data[gr].First;  // Note: This is the rank within particles of the same type

                  if(rank > UCHAR_MAX)  // restrict this to 1 byte (which is the storage we have reserved for this)
                    rank = UCHAR_MAX;

                  MtrP[i].PrevRankInSubhalo = static_cast<unsigned short>(rank);
                }
              p++;
            }
        }

      Mem.myfree(import_subhalo_data);
      Mem.myfree(export_subhalo_data);
    }

  Mem.myfree(scount);
  Mem.myfree(gcount);
  Mem.myfree(ntype_all);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

void mergertree::mergertree_match_ids_of_current_snap(void)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  assign_particle_data *AsP = (assign_particle_data *)Mem.mymalloc("AsP", sizeof(assign_particle_data) * (Sp->NumPart + 1));

  for(int i = 0; i < Sp->NumPart; i++)
    {
      AsP[i].OriginTask  = ThisTask;
      AsP[i].OriginIndex = i;
      AsP[i].ID          = Sp->P[i].ID.get();
    }

  mycxxsort_parallel(AsP, AsP + Sp->NumPart, compare_AssignP_ID, Communicator);
  mycxxsort_parallel(MtrP, MtrP + MtrP_NumPart, compare_MtrP_ID, Communicator);

  MyIDType *list_min_id = (MyIDType *)Mem.mymalloc("list_min_id", NTask * sizeof(MyIDType));
  MyIDType *list_max_id = (MyIDType *)Mem.mymalloc("list_max_id", NTask * sizeof(MyIDType));
  int *list_numpart     = (int *)Mem.mymalloc("list_mumpart", NTask * sizeof(int));

  MPI_Allgather(&AsP[0].ID, sizeof(MyIDType), MPI_BYTE, list_min_id, sizeof(MyIDType), MPI_BYTE, Communicator);
  MPI_Allgather(&AsP[Sp->NumPart > 0 ? Sp->NumPart - 1 : 0].ID, sizeof(MyIDType), MPI_BYTE, list_max_id, sizeof(MyIDType), MPI_BYTE,
                Communicator);
  MPI_Allgather(&Sp->NumPart, sizeof(int), MPI_BYTE, list_numpart, sizeof(int), MPI_BYTE, Communicator);

  int nexport = 0, nimport = 0;
  mergertree_particle_data *import_data = NULL, *export_data = NULL;

  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target = 0;

      for(int i = 0; i < MtrP_NumPart; i++)
        {
          while(target < NTask - 1 && (list_numpart[target] == 0 || MtrP[i].ID > list_max_id[target]))
            target++;

          if(list_numpart[target] != 0)
            if(MtrP[i].ID >= list_min_id[target] && MtrP[i].ID <= list_max_id[target])
              {
                if(mode == 0)
                  Send_count[target]++;
                else
                  {
                    int off          = Send_offset[target] + Send_count[target]++;
                    export_data[off] = MtrP[i];
                  }
              }
        }

      if(mode == 0)
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
          Recv_offset[0] = Send_offset[0] = 0;
          for(int j = 0; j < NTask; j++)
            {
              nimport += Recv_count[j];
              nexport += Send_count[j];
              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          export_data = (mergertree_particle_data *)Mem.mymalloc("export_data", nexport * sizeof(mergertree_particle_data));
          import_data = (mergertree_particle_data *)Mem.mymalloc("import_data", nimport * sizeof(mergertree_particle_data));
        }
    }

  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(mergertree_particle_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &import_data[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(mergertree_particle_data), MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                       MPI_STATUS_IGNORE);
    }

  /* incoming data should already be sorted, so now do the match */

  for(int i = 0, j = 0; i < Sp->NumPart && j < nimport;)
    {
      if(AsP[i].ID < import_data[j].ID)
        i++;
      else if(AsP[i].ID > import_data[j].ID)
        j++;
      else
        {
          AsP[i].PrevSubhaloNr  = import_data[j].PrevSubhaloNr;
          AsP[i].PrevSubhaloLen = import_data[j].PrevSubhaloLen;
          i++;
          j++;
        }
    }

  mycxxsort_parallel(AsP, AsP + Sp->NumPart, compare_AssignP_Origin, Communicator);

  for(int i = 0; i < Sp->NumPart; i++)
    {
      Sp->P[i].PrevSizeOfSubhalo.set(AsP[i].PrevSubhaloLen);
      Sp->P[i].PrevSubhaloNr.set(AsP[i].PrevSubhaloNr);
    }

  Mem.myfree(import_data);
  Mem.myfree(export_data);
  Mem.myfree(list_numpart);
  Mem.myfree(list_max_id);
  Mem.myfree(list_min_id);

  Mem.myfree(AsP);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

#endif
