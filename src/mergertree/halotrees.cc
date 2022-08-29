/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  halotrees.cc
 *
 *  \brief constructs trees that link together just the subhalos related by descendant relationships
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
#include "../mergertree/io_descendant.h"
#include "../mergertree/io_halotrees.h"
#include "../mergertree/io_progenitors.h"
#include "../mergertree/io_treelinks.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

/* This is the main function for constructing the halo trees.
 */
void mergertree::halotrees_construct(int lastsnapnr)
{
  domain<simparticles> Domain{Communicator, Sp};
  fof<simparticles> FoF{Communicator, Sp, &Domain};

  LastSnapShotNr = lastsnapnr;

  Cats     = (halo_catalogue *)Mem.mymalloc_clear("Cats", sizeof(halo_catalogue) * (LastSnapShotNr + 1));
  CatTimes = (times_catalogue *)Mem.mymalloc("CatTimes", sizeof(times_catalogue) * (LastSnapShotNr + 1));

  /* Let's load all catalogs. */
  halotrees_load_catalogues(&FoF);
  mpi_printf("\nMERGERTREE: Catalogues loaded\n");

  /* As preliminary work, let everybody know the number of subhalos on each rank, and assign a global SubhaloNr to each. */
  halotrees_assign_global_subhalonr_and_groupnr();
  mpi_printf("\nMERGERTREE: SubhaloNr assigned\n");

  /* Assign the halos to disjoint trees - initially, each halo is in its own tree, which will then be linked up to form bigger trees.
   * To get a reasonably memory-balanced distribution, we assign random tasks to them initially. */
  halotrees_initial_treeassignment();
  mpi_printf("\nMERGERTREE: Initial tree assignment\n");

  /* We now proceed by linking the halos according to the descendant information */
  halotrees_link_trees();
  mpi_printf("\nMERGERTREE: Halo tree linked\n");

  /* because the previously selected first progenitor is not necessarily the most significant branch, we make a new choice
   * by sorting the progenitors by their significance, and then redetermine the first progenitor among the set of progenitors */
  halotrees_determine_mainprogenitor();
  mpi_printf("\nMERGERTREE: Determined the main progenitor\n");

  /* Now the trees are linked, i.e. all halos belonging to the same tree have a common tree-id. We proceed by assigning new indices to
   * all halos within the same treeid. */
  halotrees_assign_new_treeindices();
  mpi_printf("\nMERGERTREE: New tree indices assigned\n");

  /* next, we need to remap the global descendant and first/next progenitor pointers to another set of such pointers that use the local
   * tree ids */
  halotrees_remap_treepointers();
  mpi_printf("\nMERGERTREE: Remapping done\n");

  /* now we can go ahead and move the halos belonging to the same tree together, and set up a table with the trees */
  halotrees_collect_treehalos();
  mpi_printf("\nMERGERTREE: Treehalos collected\n\n");

  /* finally, output the trees to files */
  halotrees_io TreesIO{this, this->Communicator, All.SnapFormat};
  TreesIO.halotrees_save_trees();

  mpi_printf("\nMERGERTREE: Trees saved:  A total number of %lld subhalos have been assembled in %lld trees.", TotNhalos, TotNtrees);
  mpi_printf("\nMERGERTREE: Largest tree contains %d halos.\n\n", LargestHaloCount);

  /* now also output links from the subhalo catalogues to the trees */
  halotrees_save_subhalo_treelinks();
  mpi_printf("\nMERGERTREE: Tree links saved\n\n");
}

/* This function is in charge of loading all the group catalogues.
 */
void mergertree::halotrees_load_catalogues(fof<simparticles> *FoF)
{
  /* let's first get the catalogs */
  for(int num = 0; num <= LastSnapShotNr; num++)
    {
      /* load the group catalog */
      fof_io<simparticles> FoF_IO{FoF, Communicator, All.SnapFormat};
      FoF_IO.fof_subfind_load_groups(num);

      /* bring all groups and subhalos into the order in which they were stored in the files, because the
       * parallel load routine may have mangled up the order.
       */
      mycxxsort_parallel(FoF->Group, FoF->Group + FoF->Ngroups, compare_Group_FileOffset, Communicator);
      mycxxsort_parallel(FoF->Subhalo, FoF->Subhalo + FoF->Nsubhalos, compare_Subhalo_FileOffset, Communicator);

      if(FoF_IO.LegacyFormat)
        {
          mpi_printf("\nFOF/SUBFIND: Legacy format from Arepo detected, trying to adjust for this.\n");
          FoF->subfind_redetermine_groupnr();
          FoF->fof_assign_group_offset();
          FoF->subfind_assign_subhalo_offsettype();
        }

      /* save the groups info for later use */
      Cats[num].Ngroups    = FoF->Ngroups;
      Cats[num].TotNgroups = FoF->TotNgroups;

      Cats[num].Nsubhalos    = FoF->Nsubhalos;
      Cats[num].TotNsubhalos = FoF->TotNsubhalos;

      CatTimes[num].Time     = FoF->Time;
      CatTimes[num].Redshift = FoF->Redshift;

      // we copy the group catalogue data here. Simply setting the pointer will violate the rule to not
      // copy a (movable) memory pointer into another variable, and would later give corruption once the IO object is destroyed...
      Cats[num].Group = (fof<simparticles>::group_properties *)Mem.mymalloc_movable(
          &Cats[num].Group, "Cats[num].Group", FoF->Ngroups * sizeof(fof<simparticles>::group_properties));
      memcpy(Cats[num].Group, FoF->Group, FoF->Ngroups * sizeof(fof<simparticles>::group_properties));
      Mem.myfree_movable(FoF->Group);

      Cats[num].Subhalo = (fof<simparticles>::subhalo_properties *)Mem.mymalloc_movable(
          &Cats[num].Subhalo, "Subhalo", Cats[num].Nsubhalos * sizeof(fof<simparticles>::subhalo_properties));
      memcpy(Cats[num].Subhalo, FoF->Subhalo, Cats[num].Nsubhalos * sizeof(fof<simparticles>::subhalo_properties));
      Mem.myfree_movable(FoF->Subhalo);

      /* allocate some storage for extra subhalo info */
      Cats[num].SubExt =
          (subhalo_extension *)Mem.mymalloc_movable(&Cats[num].SubExt, "Cats[num].SubExt", FoF->Nsubhalos * sizeof(subhalo_extension));
    }

  /* now let's load the descendant tree information (except for the last snapshot) */
  for(int num = 0; num < LastSnapShotNr; num++)
    {
      /* fetch the descendants for the catalogs  */
      descendant_io DescIO{this, this->Communicator, All.SnapFormat};
      DescIO.mergertree_read_descendants(num);

      /* if needed, restore proper order as in files */
      mycxxsort_parallel(Descendants, Descendants + DescIO.Nsubhalos, compare_Desc_FileOffset, Communicator);

      if(DescIO.TotNsubhalos != Cats[num].TotNsubhalos)
        Terminate("inconsistency: DescIO.TotNsubhalos=%lld != Cats[num].TotNsubhalos=%lld", DescIO.TotNsubhalos,
                  Cats[num].TotNsubhalos);

      /* it's possible that the local number of descendants, Nsubhalos, does not match Cats[num].Nsubhalos,
       * implying that we need to reshuffle how things are distributed over processors.
       */
      halotrees_reshuffle((char **)&Descendants, sizeof(desc_list), DescIO.Nsubhalos, Cats[num].Nsubhalos);

      /* retain info for later use
       * we copy the data here. Simply setting the pointer will violate the rule to not
       * copy a (movable) memory pointer into another variable, and would later give corruption once the IO object is destroyed...
       */
      Cats[num].Descendants =
          (desc_list *)Mem.mymalloc_movable(&Cats[num].Descendants, "Cats[num].Descendants", Cats[num].Nsubhalos * sizeof(desc_list));
      memcpy(Cats[num].Descendants, Descendants, Cats[num].Nsubhalos * sizeof(desc_list));
      Mem.myfree_movable(Descendants);
    }

  /* load the progenitor tree information (except for the first snapshot) */
  for(int num = 1; num <= LastSnapShotNr; num++)
    {
      /* fetch from files */
      progenitors_io ProgIO{this, this->Communicator, All.SnapFormat};
      ProgIO.mergertree_read_progenitors(num);

      /* if needed, restore proper order as in files */
      mycxxsort_parallel(Progenitors, Progenitors + ProgIO.Nsubhalos, compare_Prog_FileOffset, Communicator);

      if(ProgIO.TotNsubhalos != Cats[num].TotNsubhalos)
        Terminate("inconsistency: ProgIO.TotNsubhalos=%lld != Cats[num].TotNsubhalos=%lld", ProgIO.TotNsubhalos,
                  Cats[num].TotNsubhalos);

      /* it's possible that the local number of descendants, Nsubhalos, does not match Cats[num].Nsubhalos,
       * implying that we need to reshuffle how things are distributed over processors.
       */
      halotrees_reshuffle((char **)&Progenitors, sizeof(prog_list), ProgIO.Nsubhalos, Cats[num].Nsubhalos);

      /* retain info for later use
       * we copy the data here. Simply setting the pointer will violate the rule to not
       * copy a (movable) memory pointer into another variable, and would later give corruption once the IO object is destroyed...
       */
      Cats[num].Progenitors =
          (prog_list *)Mem.mymalloc_movable(&Cats[num].Progenitors, "Cats[num].Progenitors", Cats[num].Nsubhalos * sizeof(prog_list));
      memcpy(Cats[num].Progenitors, Progenitors, Cats[num].Nsubhalos * sizeof(prog_list));
      Mem.myfree_movable(Progenitors);
    }
}

void mergertree::halotrees_save_subhalo_treelinks(void)
{
  for(int num = 0; num <= LastSnapShotNr; num++)
    {
      treelinks_io TreeLinkIO{this, this->Communicator, All.SnapFormat};

      TreeLinkIO.Nsubhalos    = Cats[num].Nsubhalos;
      TreeLinkIO.TotNsubhalos = Cats[num].TotNsubhalos;

      TreeLink = (treelink_data *)Mem.mymalloc("TreeLink", TreeLinkIO.Nsubhalos * sizeof(treelink_data));
      for(int i = 0; i < TreeLinkIO.Nsubhalos; i++)
        {
          TreeLink[i].TreeID    = Cats[num].Subhalo[i].TreeID;
          TreeLink[i].TreeIndex = Cats[num].Subhalo[i].TreeIndex;
        }

      /* save the tree-link info */
      TreeLinkIO.treelinks_save(num);

      Mem.myfree(TreeLink);
    }
}

/* Here we collect the group and subhalo numbers stored on each processor for each snapshot, for communication purposes.
 * We also assign for the subhalos of each catalogue a running global subhalo number which is then a unique identifier within this
 * snapshot's catalogue.
 */
void mergertree::halotrees_assign_global_subhalonr_and_groupnr(void)
{
  long long grprev = 0;

  for(int num = 0; num <= LastSnapShotNr; num++)
    {
      Cats[num].TabNgroups = (int *)Mem.mymalloc("Cats[num].TabNgroups", NTask * sizeof(int));
      MPI_Allgather(&Cats[num].Ngroups, 1, MPI_INT, Cats[num].TabNgroups, 1, MPI_INT, Communicator);

      Cats[num].TabNsubhalos = (int *)Mem.mymalloc("Cats[num].TabNsubhalos", NTask * sizeof(int));
      MPI_Allgather(&Cats[num].Nsubhalos, 1, MPI_INT, Cats[num].TabNsubhalos, 1, MPI_INT, Communicator);

      long long subprev = 0;

      for(int i = 0; i < ThisTask; i++)
        subprev += Cats[num].TabNsubhalos[i];

      for(int i = 0; i < Cats[num].Nsubhalos; i++)
        Cats[num].Subhalo[i].SubhaloNr = subprev + i;

      /* Note: SubhaloNr should be now the quantity to which Descendant/FirstProgenitor/NextProgenitor refer to.
       */

      /* also set the Group[].GroupNr field (was not stored in the field that were read in)
       * In contrast, Subhalo[].GroupNr was read in  */
      long long nbefore = 0;

      for(int i = 0; i < ThisTask; i++)
        nbefore += Cats[num].TabNgroups[i];

      for(int i = 0; i < Cats[num].Ngroups; i++)
        Cats[num].Group[i].GroupNr = nbefore + i;

      /* define a GroupNr field for a unique group number */
      for(int i = 0; i < Cats[num].Nsubhalos; i++)
        Cats[num].Subhalo[i].UniqueGroupNr = Cats[num].Subhalo[i].GroupNr + grprev;

      grprev += Cats[num].TotNgroups;

      /* assign some special properties to the subhalos in the tree which come from the FOF group catalogue (like M200, etc.) */

      for(int i = 0; i < Cats[num].Nsubhalos; i++)
        Cats[num].Subhalo[i].M_Crit200 = 0;

      int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
      int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
      int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
      int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

      struct exch_data
      {
        long long GroupNr;
        MyFloat M_Crit200;
        int loc_index;
      };

      exch_data *export_data = NULL, *import_data = NULL;
      int nimport = 0, nexport = 0;

      /* for communication bookkeeping reasons, we traverse the counting pattern twice */
      for(int mode = 0; mode < 2; mode++)
        {
          for(int i = 0; i < NTask; i++)
            Send_count[i] = 0;

          int target                = 0;
          long long ngroup_previous = 0;

          for(int i = 0; i < Cats[num].Nsubhalos; i++)
            {
              /* select only the main subhalos */
              if(Cats[num].Subhalo[i].SubRankInGr == 0)
                {
                  while(target < NTask - 1 && Cats[num].Subhalo[i].GroupNr >= (ngroup_previous + Cats[num].TabNgroups[target]))
                    {
                      ngroup_previous += Cats[num].TabNgroups[target];
                      target++;
                    }

                  if(mode == 0)
                    Send_count[target]++;
                  else
                    {
                      int off = Send_offset[target] + Send_count[target]++;

                      export_data[off].loc_index = i;
                      export_data[off].GroupNr   = Cats[num].Subhalo[i].GroupNr;
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

              export_data = (exch_data *)Mem.mymalloc("export_data", nexport * sizeof(exch_data));
              import_data = (exch_data *)Mem.mymalloc("import_data", nimport * sizeof(exch_data));
            }
        }

      /* send data to target processors  */
      for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          int recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(exch_data), MPI_BYTE, recvTask,
                             TAG_DENS_B, &import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(exch_data), MPI_BYTE,
                             recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
        }

      long long firstgrnr = 0;
      for(int i = 0; i < ThisTask; i++)
        firstgrnr += Cats[num].TabNgroups[i];

      /* now read out the information we want from the groups */

      for(int i = 0; i < nimport; i++)
        {
          int index = import_data[i].GroupNr - firstgrnr;

          if(Cats[num].Group[index].GroupNr != import_data[i].GroupNr)
            Terminate(
                "bummer: num=%d i=%d Cats[num].Ngroups=%d nimport=%d  index=%d Cats[num].Group[index].GroupNr=%lld "
                "import_data[i].GroupNr=%lld\n",
                num, i, Cats[num].Ngroups, nimport, index, Cats[num].Group[index].GroupNr, import_data[i].GroupNr);

          import_data[i].M_Crit200 = Cats[num].Group[index].M_Crit200;
        }

      /* send the results back */
      for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
        {
          int recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              myMPI_Sendrecv(&import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(exch_data), MPI_BYTE, recvTask,
                             TAG_DENS_B, &export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(exch_data), MPI_BYTE,
                             recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
        }

      /* now read it out and assign the data */
      for(int i = 0; i < nexport; i++)
        Cats[num].Subhalo[export_data[i].loc_index].M_Crit200 = export_data[i].M_Crit200;

      Mem.myfree(import_data);
      Mem.myfree(export_data);

      Mem.myfree(Recv_offset);
      Mem.myfree(Recv_count);
      Mem.myfree(Send_offset);
      Mem.myfree(Send_count);
    }
}

/* Give initially each subhalo its own unique TreeID, and a randomly placed processor.
 */
void mergertree::halotrees_initial_treeassignment(void)
{
  long long previous = 0;
  for(int num = LastSnapShotNr; num >= 0; num--)
    {
      long long prevtask = 0;
      for(int i = 0; i < ThisTask; i++)
        prevtask += Cats[num].TabNsubhalos[i];

      for(int i = 0; i < Cats[num].Nsubhalos; i++)
        {
          Cats[num].Subhalo[i].TreeTask = (int)(get_random_number() * NTask);
          Cats[num].Subhalo[i].TreeID   = previous + prevtask + i;
        }

      previous += Cats[num].TotNsubhalos;
    }
}

void mergertree::halotrees_select_interior_min_newtreeid(int mode, tlink *treehalos, long long totnsubs)
{
  for(int backforth = -1; backforth <= 1; backforth += 2)
    {
      long long i;

      if(backforth == -1)
        i = 0;
      else
        i = totnsubs - 2;

      while(i >= 0 && i <= totnsubs - 2)
        {
          if((mode == 0 && treehalos[i].UniqueGroupNr == treehalos[i + 1].UniqueGroupNr) ||
             (mode == 1 && treehalos[i].TreeID == treehalos[i + 1].TreeID))
            {
              if(treehalos[i].NewTreeID > treehalos[i + 1].NewTreeID)
                {
                  treehalos[i].NewTreeID   = treehalos[i + 1].NewTreeID;
                  treehalos[i].NewTreeTask = treehalos[i + 1].NewTreeTask;
                }
              else if(treehalos[i].NewTreeID < treehalos[i + 1].NewTreeID)
                {
                  treehalos[i + 1].NewTreeID   = treehalos[i].NewTreeID;
                  treehalos[i + 1].NewTreeTask = treehalos[i].NewTreeTask;
                }
            }
          if(backforth == -1)
            i++;
          else
            i--;
        }
    }
}

long long mergertree::halotrees_join_trees_via_fof_or_mostboundid_bridges(int mode)
{
  long long totnsubs = 0;

  for(int num = 0; num <= LastSnapShotNr; num++)
    totnsubs += Cats[num].Nsubhalos;

  tlink *treehalos = (tlink *)Mem.mymalloc("treehalos", (totnsubs + 1) * sizeof(tlink));

  long long count = 0;

  for(int num = 0; num <= LastSnapShotNr; num++)
    for(int i = 0; i < Cats[num].Nsubhalos; i++)
      {
        treehalos[count].TreeID   = Cats[num].Subhalo[i].TreeID;
        treehalos[count].TreeTask = Cats[num].Subhalo[i].TreeTask;

        treehalos[count].NewTreeID   = Cats[num].Subhalo[i].TreeID;
        treehalos[count].NewTreeTask = Cats[num].Subhalo[i].TreeTask;

        treehalos[count].OrigTask   = ThisTask;
        treehalos[count].OrderIndex = count;

        if(mode == 0)
          treehalos[count].UniqueGroupNr = Cats[num].Subhalo[i].UniqueGroupNr;
        else
          treehalos[count].UniqueGroupNr = Cats[num].Subhalo[i].SubMostBoundID;

        count++;
      }

  long long *count_list = (long long *)Mem.mymalloc("count_list", NTask * sizeof(long long));
  MPI_Allgather(&totnsubs, sizeof(long long), MPI_BYTE, count_list, sizeof(long long), MPI_BYTE, Communicator);

  int next_task = -1;
  for(int i = ThisTask + 1; i < NTask; i++)
    if(count_list[i] > 0)
      {
        next_task = i;
        break;
      }

  int prev_task = -1;
  for(int i = ThisTask - 1; i >= 0; i--)
    if(count_list[i] > 0)
      {
        prev_task = i;
        break;
      }

  tlink *elem_first = (tlink *)Mem.mymalloc("elem_first", NTask * sizeof(tlink));
  tlink *elem_last  = (tlink *)Mem.mymalloc("elem_last", NTask * sizeof(tlink));

  for(int mode = 0; mode < 2; mode++)
    {
      if(mode == 0)
        {
          /* bring halos in the same group together */
          mycxxsort_parallel(treehalos, treehalos + totnsubs, compare_tlink_GroupNr, Communicator);
          halotrees_select_interior_min_newtreeid(mode, treehalos, totnsubs);
        }
      else
        {
          /* bring halos in the same tree together */
          mycxxsort_parallel(treehalos, treehalos + totnsubs, compare_tlink_TreeID, Communicator);
          halotrees_select_interior_min_newtreeid(mode, treehalos, totnsubs);
        }

      long long totchanges = 0;
      int iter             = 0;
      do
        {
          int changes = 0;

          /* note: the 0th element is guaranteed to be allocated even on ranks with zero totnsubs */
          MPI_Allgather(&treehalos[0], sizeof(tlink), MPI_BYTE, elem_first, sizeof(tlink), MPI_BYTE, Communicator);
          MPI_Allgather(&treehalos[totnsubs > 0 ? totnsubs - 1 : 0], sizeof(tlink), MPI_BYTE, elem_last, sizeof(tlink), MPI_BYTE,
                        Communicator);

          if(prev_task >= 0 && totnsubs > 0 &&
             ((mode == 0 && elem_last[prev_task].UniqueGroupNr == treehalos[0].UniqueGroupNr) ||
              (mode == 1 && elem_last[prev_task].TreeID == treehalos[0].TreeID)) &&
             elem_last[prev_task].NewTreeID < treehalos[0].NewTreeID)
            {
              treehalos[0].NewTreeID   = elem_last[prev_task].NewTreeID;
              treehalos[0].NewTreeTask = elem_last[prev_task].NewTreeTask;
              changes++;
            }

          if(next_task >= 0 && totnsubs > 0 &&
             ((mode == 0 && elem_first[next_task].UniqueGroupNr == treehalos[totnsubs - 1].UniqueGroupNr) ||
              (mode == 1 && elem_first[next_task].TreeID == treehalos[totnsubs - 1].TreeID)) &&
             elem_first[next_task].NewTreeID < treehalos[totnsubs - 1].NewTreeID)
            {
              treehalos[totnsubs - 1].NewTreeID   = elem_first[next_task].NewTreeID;
              treehalos[totnsubs - 1].NewTreeTask = elem_first[next_task].NewTreeTask;
              changes++;
            }

          halotrees_select_interior_min_newtreeid(mode, treehalos, totnsubs);

          sumup_large_ints(1, &changes, &totchanges, Communicator);

          if(++iter > MAXITER)
            Terminate("too many iterations. mode=%d changes=%d", mode, changes);
        }
      while(totchanges > 0);
    }

  Mem.myfree(elem_last);
  Mem.myfree(elem_first);
  Mem.myfree(count_list);

  /* now bring halos into original order */
  mycxxsort_parallel(treehalos, treehalos + totnsubs, compare_tlink_OrigTask_OrderIndex, Communicator);

  /* transfer new TreeIDs back */
  count = 0;

  long long flips = 0;

  for(int num = 0; num <= LastSnapShotNr; num++)
    for(int i = 0; i < Cats[num].Nsubhalos; i++)
      {
        if(treehalos[count].NewTreeID != treehalos[count].TreeID)
          {
            Cats[num].Subhalo[i].TreeID   = treehalos[count].NewTreeID;
            Cats[num].Subhalo[i].TreeTask = treehalos[count].NewTreeTask;

            flips++;
          }

        count++;
      }

  Mem.myfree(treehalos);

  return flips;
}

/* This function is in charge of linking the halos according to the descendant information.
 */
void mergertree::halotrees_link_trees(void)
{
  /*  we will build up the trees by assigning any two subhalos which are in the same tree the same TreeID, and the same TreeTask */

  long long changesA = 0, changesB = 0;
  int iter = 0;

  do
    {
      changesA = 0;

      /* propagate the TreeIDs from level num to level num-1, via the descendant relationship from num-1 to num */
      if(iter & 1)
        for(int num = 1; num <= LastSnapShotNr; num++)
          changesA += halotrees_join_via_descendants(num);
      else
        for(int num = LastSnapShotNr; num > 0; num--)
          changesA += halotrees_join_via_descendants(num);

      changesB = 0;

      /* propagate the TreeIDs from level num to level num+1, via the progenitor relationship from num+1 to num */
      if(iter & 1)
        for(int num = LastSnapShotNr - 1; num >= 0; num--)
          changesB += halotrees_join_via_progenitors(num);
      else
        for(int num = 0; num < LastSnapShotNr; num++)
          changesB += halotrees_join_via_progenitors(num);

      MPI_Allreduce(MPI_IN_PLACE, &changesA, 1, MPI_LONG_LONG, MPI_SUM, Communicator);
      MPI_Allreduce(MPI_IN_PLACE, &changesB, 1, MPI_LONG_LONG, MPI_SUM, Communicator);

      mpi_printf("MERGERTREE: iteration %d of joining trees via descendants and progenitor,  %lld   %lld  links\n", iter, changesA,
                 changesB);

      if(++iter > MAXITER)
        Terminate("too many iterations");
    }
  while(changesA + changesB > 0);

  /* Now link all subhalos that are in a common FOF halo, if they have most bound IDs in common, or if they are linked by a progenitor
   * relation, so that they appear in the same tree */
  for(int mode = 0; mode < 2; mode++)
    {
      long long totlinks = 0;
      int iter           = 0;
      do
        {
          long long links = halotrees_join_trees_via_fof_or_mostboundid_bridges(mode);

          MPI_Allreduce(&links, &totlinks, 1, MPI_LONG_LONG, MPI_SUM, Communicator);

          if(mode == 0)
            mpi_printf("MERGERTREE: iteration %d of joining trees via common FOF bridges yielded %lld links\n", iter, totlinks);
          else
            mpi_printf("MERGERTREE: iteration %d of joining trees via common MostboundID bridges yielded %lld links\n", iter,
                       totlinks);

          if(++iter > MAXITER)
            Terminate("too many iterations");
        }
      while(totlinks > 0);
    }
}

/* this functions brings the halos of trees together according to the TreeID, so that they are ready for output to file.
 * we are also assigning new TreeIDs to get a contiguous numbering, and we set up a table that gives the number of
 * halos for each tree, and a cumulative offset into the output file where the tree starts
 */
void mergertree::halotrees_collect_treehalos(void)
{
  Nhalos = 0;

  for(int num = 0; num <= LastSnapShotNr; num++)
    Nhalos += Cats[num].Nsubhalos;

  /* get the total number of all subhalos in all trees */
  sumup_large_ints(1, &Nhalos, &TotNhalos, Communicator);

  /* we now release some memory that is not needed any more in order to reduce peak memory usage */
  for(int num = LastSnapShotNr; num >= 0; num--)
    {
      Mem.myfree(Cats[num].TabNsubhalos);
      Cats[num].TabNsubhalos = NULL;

      Mem.myfree(Cats[num].TabNgroups);
      Cats[num].TabNgroups = NULL;
    }

  for(int num = LastSnapShotNr; num >= 1; num--)
    {
      Mem.myfree(Cats[num].Progenitors);
      Cats[num].Progenitors = NULL;
    }

  for(int num = LastSnapShotNr - 1; num >= 0; num--)
    {
      Mem.myfree(Cats[num].Descendants);
      Cats[num].Descendants = NULL;
    }

  for(int num = LastSnapShotNr; num >= 0; num--)
    {
      Mem.myfree_movable(Cats[num].Group);
      Cats[num].Group = NULL;
    }

  Halos = (treehalo_type *)Mem.mymalloc_movable(&Halos, "Halos", (Nhalos + 1) * sizeof(treehalo_type));

  /* set up the halo data for the tree output */
  long long off = 0;
  for(int num = 0; num <= LastSnapShotNr; num++)
    for(int i = 0; i < Cats[num].Nsubhalos; i++)
      {
        Halos[off].TreeID    = Cats[num].Subhalo[i].TreeID;
        Halos[off].TreeIndex = Cats[num].Subhalo[i].TreeIndex;

        Halos[off].TreeMainProgenitor  = Cats[num].SubExt[i].TreeMainProgenitor;
        Halos[off].TreeFirstProgenitor = Cats[num].SubExt[i].TreeFirstProgenitor;
        Halos[off].TreeNextProgenitor  = Cats[num].SubExt[i].TreeNextProgenitor;
        Halos[off].TreeDescendant      = Cats[num].SubExt[i].TreeDescendant;

        Halos[off].TreeFirstDescendant = Cats[num].SubExt[i].TreeFirstDescendant;
        Halos[off].TreeNextDescendant  = Cats[num].SubExt[i].TreeNextDescendant;
        Halos[off].TreeProgenitor      = Cats[num].SubExt[i].TreeProgenitor;

        Halos[off].SnapNum       = num;
        Halos[off].SubhaloNr     = Cats[num].Subhalo[i].SubhaloNr;
        Halos[off].GroupNr       = Cats[num].Subhalo[i].GroupNr;
        Halos[off].UniqueGroupNr = Cats[num].Subhalo[i].UniqueGroupNr;

        Halos[off].SubProp = Cats[num].Subhalo[i];

        off++;
      }

  /* release some memory to reduce peak memory usage */
  for(int num = LastSnapShotNr; num >= 0; num--)
    {
      Mem.myfree_movable(Cats[num].SubExt);
      Cats[num].SubExt = NULL;
    }

  long long *count_list = (long long *)Mem.mymalloc("count_list", NTask * sizeof(long long));
  long long totnsubs    = Nhalos;
  MPI_Allgather(&totnsubs, sizeof(long long), MPI_BYTE, count_list, sizeof(long long), MPI_BYTE, Communicator);

  /* first set the fields TreeFirstHaloInFOFgroup and TreeNextHaloInFOFgroup */
  mycxxsort_parallel(Halos, Halos + Nhalos, compare_Halos_UniqueGroupNr_SubhaloNr, Communicator);

  long long previous_UniqueGroupNr = HALONR_MAX;
  long long new_treeindex          = HALONR_MAX;
  for(int i = 0; i < Nhalos; i++)
    {
      if(Halos[i].UniqueGroupNr != previous_UniqueGroupNr)
        {
          previous_UniqueGroupNr = Halos[i].UniqueGroupNr;
          new_treeindex          = Halos[i].TreeIndex;
        }

      Halos[i].TreeFirstHaloInFOFgroup = new_treeindex;

      if(i < Nhalos - 1 && Halos[i].UniqueGroupNr == Halos[i + 1].UniqueGroupNr)
        Halos[i].TreeNextHaloInFOFgroup = Halos[i + 1].TreeIndex;
      else
        Halos[i].TreeNextHaloInFOFgroup = -1;
    }

  treehalo_type *efirst = (treehalo_type *)Mem.mymalloc("efirst", NTask * sizeof(treehalo_type));
  treehalo_type *elast  = (treehalo_type *)Mem.mymalloc("elast", NTask * sizeof(treehalo_type));
  /* note: the 0th element is guaranteed to be allocated even on ranks with zero totnsubs */
  MPI_Allgather(&Halos[0], sizeof(treehalo_type), MPI_BYTE, efirst, sizeof(treehalo_type), MPI_BYTE, Communicator);
  MPI_Allgather(&Halos[Nhalos > 0 ? Nhalos - 1 : 0], sizeof(treehalo_type), MPI_BYTE, elast, sizeof(treehalo_type), MPI_BYTE,
                Communicator);

  if(Nhalos > 0)
    for(int task = ThisTask + 1; task < NTask; task++)
      if(count_list[task] > 0)
        {
          if(Halos[Nhalos - 1].UniqueGroupNr == efirst[task].UniqueGroupNr)
            Halos[Nhalos - 1].TreeNextHaloInFOFgroup = efirst[task].TreeIndex;

          break;
        }

  if(Nhalos > 0)
    {
      long long previous_UniqueGroupNr = HALONR_MAX;
      long long new_treeindex          = HALONR_MAX;

      for(int task = ThisTask - 1; task >= 0; task--)
        {
          if(count_list[task] > 0)
            {
              if(elast[task].UniqueGroupNr == Halos[0].UniqueGroupNr)
                {
                  previous_UniqueGroupNr = elast[task].UniqueGroupNr;
                  new_treeindex          = elast[task].TreeFirstHaloInFOFgroup;
                }
            }
        }

      for(int i = 0; i < Nhalos; i++)
        {
          if(Halos[i].UniqueGroupNr != previous_UniqueGroupNr)
            {
              previous_UniqueGroupNr = Halos[i].UniqueGroupNr;
              new_treeindex          = Halos[i].TreeIndex;
              break;
            }

          Halos[i].TreeFirstHaloInFOFgroup = new_treeindex;
        }
    }

  Mem.myfree(elast);
  Mem.myfree(efirst);

  /* parallel sort according to treeid and treeindex - this brings the tree halos together, in ascending order */
  mycxxsort_parallel(Halos, Halos + Nhalos, compare_Halos_TreeID_TreeIndex, Communicator);

  treehalo_type *elem_last = (treehalo_type *)Mem.mymalloc("elem_last", NTask * sizeof(treehalo_type));
  /* now count the trees and overwrite the TreeIDs with new continuous IDs starting from zero */
  MPI_Allgather(&Halos[Nhalos > 0 ? Nhalos - 1 : 0], sizeof(treehalo_type), MPI_BYTE, elem_last, sizeof(treehalo_type), MPI_BYTE,
                Communicator);

  long long treeid_previous = -1;
  for(int task = ThisTask - 1; task >= 0; task--)
    {
      if(count_list[task] > 0)
        {
          treeid_previous = elem_last[task].TreeID;
          break;
        }
    }

  Ntrees = 0;
  for(int i = 0; i < Nhalos; i++)
    {
      if(Halos[i].TreeID != treeid_previous)
        Ntrees++;

      treeid_previous = Halos[i].TreeID;
      Halos[i].TreeID = Ntrees - 1;
    }

  int *ntrees_list = (int *)Mem.mymalloc("ntrees_list", NTask * sizeof(int));
  MPI_Allgather(&Ntrees, 1, MPI_INT, ntrees_list, 1, MPI_INT, Communicator);

  TotNtrees               = 0;
  long long ntrees_before = 0;
  for(int i = 0; i < NTask; i++)
    {
      TotNtrees += ntrees_list[i];
      if(i < ThisTask)
        ntrees_before += ntrees_list[i];
    }

  for(int i = 0; i < Nhalos; i++)
    Halos[i].TreeID += ntrees_before;

  /* let's now allocate the table of trees, with an extra [-1] element */
  TreeTable = (halotrees_table *)Mem.mymalloc_movable(&TreeTable, "TreeTable", (Ntrees + 1) * sizeof(halotrees_table));
  memset(TreeTable, 0, (Ntrees + 1) * sizeof(halotrees_table));
  TreeTable += 1;

  /* update what we have stored for the last element */
  MPI_Allgather(&Halos[Nhalos > 0 ? Nhalos - 1 : 0], sizeof(treehalo_type), MPI_BYTE, elem_last, sizeof(treehalo_type), MPI_BYTE,
                Communicator);

  treeid_previous = -1;
  for(int task = ThisTask - 1; task >= 0; task--)
    {
      if(count_list[task] > 0)
        {
          treeid_previous = elem_last[task].TreeID;
          break;
        }
    }

  Ntrees = 0;
  for(int i = 0; i < Nhalos; i++)
    {
      if(Halos[i].TreeID != treeid_previous)
        Ntrees++;

      treeid_previous = Halos[i].TreeID;
      TreeTable[Ntrees - 1].HaloCount++;
      TreeTable[Ntrees - 1].TreeID = Halos[i].TreeID;
    }

  /* in the TreeTable[-1] element we have counted segments of trees that do not start on the local processor */

  halotrees_table *elem_first = (halotrees_table *)Mem.mymalloc("elem_first", NTask * sizeof(halotrees_table));

  MPI_Allgather(&TreeTable[-1], sizeof(halotrees_table), MPI_BYTE, elem_first, sizeof(halotrees_table), MPI_BYTE, Communicator);

  if(Ntrees > 0)
    for(int task = ThisTask + 1; task < NTask; task++)
      {
        if(TreeTable[Ntrees - 1].TreeID == elem_first[task].TreeID)
          TreeTable[Ntrees - 1].HaloCount += elem_first[task].HaloCount;
        else
          break;
      }

  Mem.myfree(elem_first);

  long long sumhalos = 0;
  LargestHaloCount   = 0;
  for(int i = 0; i < Ntrees; i++)
    {
      sumhalos += TreeTable[i].HaloCount;
      if(TreeTable[i].HaloCount > LargestHaloCount)
        LargestHaloCount = TreeTable[i].HaloCount;
    }
  MPI_Allreduce(MPI_IN_PLACE, &LargestHaloCount, 1, MPI_INT, MPI_MAX, Communicator);

  long long *list_sumhalos = (long long *)Mem.mymalloc("list_sumhalos", NTask * sizeof(long long));
  MPI_Allgather(&sumhalos, sizeof(long long), MPI_BYTE, list_sumhalos, sizeof(long long), MPI_BYTE, Communicator);

  long long sum_check = 0;
  for(int i = 0; i < NTask; i++)
    sum_check += list_sumhalos[i];

  if(sum_check != TotNhalos)
    Terminate("TotNTrees=%lld  sum_check=%lld  !=  TotNhalos=%lld", (long long)TotNtrees, sum_check, (long long)TotNhalos);

  if(Ntrees > 0)
    {
      TreeTable[0].FirstHalo = 0;
      for(int task = 0; task < ThisTask; task++)
        TreeTable[0].FirstHalo += list_sumhalos[task];
    }

  Mem.myfree(list_sumhalos);

  for(int i = 1; i < Ntrees; i++)
    TreeTable[i].FirstHalo = TreeTable[i - 1].FirstHalo + TreeTable[i - 1].HaloCount;

  Mem.myfree_movable(ntrees_list);
  Mem.myfree_movable(elem_last);
  Mem.myfree_movable(count_list);
}

/*--------------------------------------------------------------------------------------------------------------*/

/* This function redistributes a subdivided global array with local pieces stored in an address pointed to by ptr.
 * On the local processor, there are currently 'ncurrent' elements, but we want 'ntarget' elements. The size of
 * one element is given by 'len'. The buffer, whose address is stored in *ptr, must be resizable and will be
 * reallocated to the new size of ntarget*len bytes.
 */
void mergertree::halotrees_reshuffle(char **ptr, size_t len, int ncurrent, int ntarget)
{
  int *Send_count  = (int *)Mem.mymalloc_movable(&Send_count, "Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc_movable(&Send_offset, "Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc_movable(&Recv_count, "Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc_movable(&Recv_offset, "Recv_offset", sizeof(int) * NTask);

  /* copy current data to an auxiliary buffer */
  char *buf = (char *)Mem.mymalloc_movable(&buf, "buf", ncurrent * len);
  memcpy(buf, *ptr, ncurrent * len);

  /* now resize source buffer of data to be able to accommodate the data of the desired target size */
  *ptr = (char *)Mem.myrealloc_movable(*ptr, ntarget * len);

  /* collect the current and target layout of the array */
  int *tab_ncurrent = (int *)Mem.mymalloc("tab_ncurrent", NTask * sizeof(int));
  MPI_Allgather(&ncurrent, 1, MPI_INT, tab_ncurrent, 1, MPI_INT, Communicator);

  int *tab_ntarget = (int *)Mem.mymalloc("tab_ntarget", NTask * sizeof(int));
  MPI_Allgather(&ntarget, 1, MPI_INT, tab_ntarget, 1, MPI_INT, Communicator);

  /* now work out where our local data should go */
  int nimport = 0;

  for(int i = 0; i < NTask; i++)
    Send_count[i] = 0;

  int nbefore = 0;
  for(int i = 0; i < ThisTask; i++)
    nbefore += tab_ncurrent[i];

  int target = 0, ncum = 0;
  for(int i = 0; i < ncurrent; i++)
    {
      while(target < NTask - 1 && nbefore + i >= ncum + tab_ntarget[target])
        ncum += tab_ntarget[target++];

      Send_count[target]++;
    }

  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
  Recv_offset[0] = Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  if(nimport != ntarget)
    Terminate("nimport != ntarget");

  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&buf[Send_offset[recvTask] * len], Send_count[recvTask] * len, MPI_BYTE, recvTask, TAG_DENS_B,
                         *ptr + Recv_offset[recvTask] * len, Recv_count[recvTask] * len, MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                         MPI_STATUS_IGNORE);
    }

  Mem.myfree(tab_ntarget);
  Mem.myfree(tab_ncurrent);
  Mem.myfree(buf);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

/* This function sets up new pointer for navigating within a given tree based on TreeIndex.
 * This is done by translating the old pointers to the new ones within a tree.
 */
void mergertree::halotrees_remap_treepointers(void)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* initialize new pointers to default values */
  for(int num = 0; num <= LastSnapShotNr; num++)
    for(int i = 0; i < Cats[num].Nsubhalos; i++)
      {
        Cats[num].SubExt[i].TreeMainProgenitor  = -1;
        Cats[num].SubExt[i].TreeFirstProgenitor = -1;
        Cats[num].SubExt[i].TreeNextProgenitor  = -1;
        Cats[num].SubExt[i].TreeDescendant      = -1;
        Cats[num].SubExt[i].TreeFirstDescendant = -1;
        Cats[num].SubExt[i].TreeNextDescendant  = -1;
        Cats[num].SubExt[i].TreeProgenitor      = -1;
      }

  /* the three pointers we have to set point either one snapshot back (e.g. FirstProgenitor, delta=-1),
   * they are in the same snapshot (NextProgenitor, delta=0), or they point into the future (Descendant, delta=+1)
   */
  for(int set = 0; set < 3; set++)
    for(int delta = 1; delta >= -1; delta--)
      {
        if(set == 2 && delta != -1)
          continue;

        /* note that we take advantage of the fact that the Subgroups are ordered according to their global SubhaloNr */

        /* in snapshot number 'base' we want to remap the pointers */
        for(int base = 0; base <= LastSnapShotNr; base++)
          {
            /* this is the snapshot we are pointing to */
            int num = base + delta;

            if(set == 1 && delta != 1 && num == 0)
              continue;

            if((delta == -1 && num >= 0) || (delta >= 0 && base < LastSnapShotNr))
              {
                /* lets get the minimum and maximum subhalo numbers for the target snapshot number as a function of the task */
                long long *list_min_subhalonr = (long long *)Mem.mymalloc("list_min_subhalonr", NTask * sizeof(long long));
                long long *list_max_subhalonr = (long long *)Mem.mymalloc("list_max_subhalonr", NTask * sizeof(long long));

                long long min_subhalonr =
                    Cats[num].Nsubhalos > 0 ? Cats[num].Subhalo[0].SubhaloNr : HALONR_MAX;  // HALONR_MAX means empty here
                long long max_subhalonr = Cats[num].Nsubhalos > 0 ? Cats[num].Subhalo[Cats[num].Nsubhalos - 1].SubhaloNr
                                                                  : HALONR_MAX;  // HALONR_MAX means empty here

                MPI_Allgather(&min_subhalonr, sizeof(long long), MPI_BYTE, list_min_subhalonr, sizeof(long long), MPI_BYTE,
                              Communicator);
                MPI_Allgather(&max_subhalonr, sizeof(long long), MPI_BYTE, list_max_subhalonr, sizeof(long long), MPI_BYTE,
                              Communicator);

                /* prepare a list of the current pointer values on the local processor */
                data_list *list = (data_list *)Mem.mymalloc("list", Cats[base].Nsubhalos * sizeof(data_list));
                int count       = 0;
                for(int i = 0; i < Cats[base].Nsubhalos; i++)
                  {
                    if(set == 0)
                      switch(delta)
                        {
                          case -1:
                            list[count].targetsubhalonr = Cats[base].Progenitors[i].FirstProgSubhaloNr;
                            break;
                          case 0:
                            list[count].targetsubhalonr = Cats[base].Descendants[i].NextProgSubhaloNr;
                            break;
                          case +1:
                            list[count].targetsubhalonr = Cats[base].Descendants[i].DescSubhaloNr;
                            break;
                        }
                    else if(set == 1)
                      switch(delta)
                        {
                          case -1:
                            list[count].targetsubhalonr = Cats[base].Progenitors[i].ProgSubhaloNr;
                            break;
                          case 0:
                            list[count].targetsubhalonr = Cats[base].Progenitors[i].NextDescSubhaloNr;
                            break;
                          case +1:
                            list[count].targetsubhalonr = Cats[base].Descendants[i].FirstDescSubhaloNr;
                            break;
                        }
                    else
                      switch(delta)
                        {
                          case -1:
                            list[count].targetsubhalonr = Cats[base].Progenitors[i].MainProgSubhaloNr;
                            break;
                        }

                    list[count].originsubhalonr = Cats[base].Subhalo[i].SubhaloNr;

                    if(list[count].targetsubhalonr >= 0 &&
                       list[count].targetsubhalonr < (long long)HALONR_MAX)  // need to only process existing pointers
                      {
                        if(list[count].targetsubhalonr >= Cats[num].TotNsubhalos)
                          Terminate("set=%d delta=%d num=%d base=%d list[count].targetsubhalonr=%lld >= Cats[num].TotNsubhalos=%lld\n",
                                    set, delta, num, base, list[count].targetsubhalonr, Cats[num].TotNsubhalos);

                        list[count].origin   = i;
                        list[count].intreeid = Cats[base].Subhalo[i].TreeID;
                        count++;
                      }
                  }

                /* sort it by current pointer value (which is the subhalonr) */
                mycxxsort(list, list + count, compare_data_list_subhalonnr);

                /* we now need to send them to other target processors since they may contain the corresponding subhalos */

                int nexport = 0, nimport = 0;

                remap_data *import_data = NULL, *export_data = NULL;

                for(int mode = 0; mode < 2; mode++)  // go through this twice to simplify bookkeeping
                  {
                    for(int i = 0; i < NTask; i++)
                      Send_count[i] = 0;

                    int target = 0;

                    for(int i = 0; i < count; i++)
                      {
                        while(target < NTask - 1 &&
                              (list_min_subhalonr[target] == HALONR_MAX || list[i].targetsubhalonr > list_max_subhalonr[target]))
                          target++;

                        if(list_min_subhalonr[target] != HALONR_MAX && list[i].targetsubhalonr >= list_min_subhalonr[target] &&
                           list[i].targetsubhalonr <= list_max_subhalonr[target])
                          {
                            if(mode == 0)
                              Send_count[target]++;
                            else
                              {
                                int off = Send_offset[target] + Send_count[target]++;

                                export_data[off].loc_index       = i;
                                export_data[off].targetsubhalonr = list[i].targetsubhalonr;
                                export_data[off].originsubhalonr = list[i].originsubhalonr;
                                export_data[off].intreeid        = list[i].intreeid;
                              }
                          }
                        else
                          Terminate(
                              "this shouldn't happen:  set=%d delta=%d num=%d base=%d delta=%d   i=%d|count=%d   "
                              "list[i].targetsubhalonr=%lld   "
                              "target=%d  "
                              "list_min_subhalonr[target]=%lld  list_max_subhalonr[target]=%lld\n",
                              set, delta, num, base, delta, i, count, (long long)list[i].targetsubhalonr, target,
                              list_min_subhalonr[target], list_max_subhalonr[target]);
                      }

                    if(mode == 0)  // prepare offset tables
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

                        export_data = (remap_data *)Mem.mymalloc("export_data", nexport * sizeof(remap_data));
                        import_data = (remap_data *)Mem.mymalloc("import_data", nimport * sizeof(remap_data));
                      }
                  }

                for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
                  {
                    /* note: here we also have a transfer from each task to itself (for ngrp=0) */
                    int recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask)
                      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                        myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(remap_data), MPI_BYTE,
                                       recvTask, TAG_DENS_B, &import_data[Recv_offset[recvTask]],
                                       Recv_count[recvTask] * sizeof(remap_data), MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                                       MPI_STATUS_IGNORE);
                  }

                /* incoming data is not necessarily be sorted according to subhalorn, that's why we need to sort it now */
                for(int i = 0; i < nimport; i++)
                  import_data[i].orig_index = i;

                mycxxsort(import_data, import_data + nimport, compare_remap_data_subhalonr);

                /* now we check the incoming data, and prepare new values for the pointers */
                for(int i = 0, j = 0; i < Cats[num].Nsubhalos && j < nimport;)
                  {
                    if(Cats[num].Subhalo[i].SubhaloNr < import_data[j].targetsubhalonr)
                      i++;
                    else if(Cats[num].Subhalo[i].SubhaloNr > import_data[j].targetsubhalonr)
                      {
                        Terminate("Can't find targetsubhalonr=%lld Cats[num].Subhalo[i].SubhaloNr=%lld  i=%d|%d j=%d|%d",
                                  import_data[j].targetsubhalonr, Cats[num].Subhalo[i].SubhaloNr, i, Cats[num].Nsubhalos, j, nimport);
                      }
                    else  // we have a match
                      {
                        import_data[j].new_treeindexptr = Cats[num].Subhalo[i].TreeIndex;
                        import_data[j].treeid           = Cats[num].Subhalo[i].TreeID;

                        if(Cats[num].Subhalo[i].TreeID != import_data[j].intreeid)
                          Terminate(
                              "\nWe are not in the same tree, which shouldn't be the case:  set=%d delta=%d  num=%d\n"
                              "Cats[num].Subhalo[i].SubhaloNr=%lld \nCats[num].Subhalo[i].DescSubhaloNr=%lld\n"
                              "Cats[num].Subhalo[i].NextProgSubhalor=%lld\nCats[num].Subhalo[i].FirstProgSubhalor=%lld\n"
                              "Cats[num].Subhalo[i].NextDescSubhalor=%lld\nCats[num].Subhalo[i].FirstDescSubhalor=%lld\n"
                              "Cats[num].Subhalo[i].ProgSubhalorNr=%lld\n"
                              "originsubgalonr=%lld targetsubgalonr=%lld Cats[num].Subhalo[i].TreeID=%lld   "
                              "import_data[j].intreeid=%lld   i=%d|%d  j=%d|%d  num=%d",
                              set, delta, num, Cats[num].Subhalo[i].SubhaloNr, Cats[num].Descendants[i].DescSubhaloNr,
                              Cats[num].Descendants[i].NextProgSubhaloNr, Cats[num].Progenitors[i].FirstProgSubhaloNr,
                              Cats[num].Progenitors[i].NextDescSubhaloNr, Cats[num].Descendants[i].FirstDescSubhaloNr,
                              Cats[num].Progenitors[i].ProgSubhaloNr, import_data[j].originsubhalonr, import_data[j].targetsubhalonr,
                              Cats[num].Subhalo[i].TreeID, import_data[j].intreeid, i, Cats[num].Nsubhalos, j, nimport, num);

                        j++;
                      }
                  }

                mycxxsort(import_data, import_data + nimport, compare_remap_data_orig_index);

                /* send the results back */
                for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
                  {
                    /* note: here we also have a transfer from each task to itself (for ngrp=0) */
                    int recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask)
                      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                        myMPI_Sendrecv(&import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(remap_data), MPI_BYTE,
                                       recvTask, TAG_DENS_B, &export_data[Send_offset[recvTask]],
                                       Send_count[recvTask] * sizeof(remap_data), MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                                       MPI_STATUS_IGNORE);
                  }

                for(int i = 0; i < nexport; i++)
                  {
                    int idx = list[export_data[i].loc_index].origin;

                    if(Cats[base].Subhalo[idx].TreeID != export_data[i].treeid)
                      Terminate("something is wrong: delta=%d  i=%d|%d  idx=%d|%d     %lld %lld %lld", delta, i, nexport, idx,
                                Cats[base].Nsubhalos, Cats[base].Subhalo[idx].TreeID, export_data[i].treeid, export_data[i].intreeid);

                    if(set == 0)
                      switch(delta)
                        {
                          case -1:
                            Cats[base].SubExt[idx].TreeFirstProgenitor = export_data[i].new_treeindexptr;
                            break;

                          case 0:
                            Cats[base].SubExt[idx].TreeNextProgenitor = export_data[i].new_treeindexptr;
                            break;

                          case +1:
                            Cats[base].SubExt[idx].TreeDescendant = export_data[i].new_treeindexptr;
                            break;
                        }
                    else if(set == 1)
                      switch(delta)
                        {
                          case -1:
                            Cats[base].SubExt[idx].TreeProgenitor = export_data[i].new_treeindexptr;
                            break;

                          case 0:
                            Cats[base].SubExt[idx].TreeNextDescendant = export_data[i].new_treeindexptr;
                            break;

                          case +1:
                            Cats[base].SubExt[idx].TreeFirstDescendant = export_data[i].new_treeindexptr;
                            break;
                        }
                    else
                      switch(delta)
                        {
                          case -1:
                            Cats[base].SubExt[idx].TreeMainProgenitor = export_data[i].new_treeindexptr;
                            break;
                        }
                  }

                Mem.myfree(import_data);
                Mem.myfree(export_data);
                Mem.myfree(list);
                Mem.myfree(list_max_subhalonr);
                Mem.myfree(list_min_subhalonr);
              }
          }
      }

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

/*--------------------------------------------------------------------------------------------------------------*/

/* This function creates a new subhalo numbering within each tree, as identified by the set of subhalos with the same treeid.
 */
void mergertree::halotrees_assign_new_treeindices(void)
{
  long long totnsubs = 0;

  for(int num = 0; num <= LastSnapShotNr; num++)
    totnsubs += Cats[num].Nsubhalos;

  assign_data *a_data = (assign_data *)Mem.mymalloc("a_data", (totnsubs + 1) * sizeof(assign_data));

  long long count = 0;
  for(int num = 0; num <= LastSnapShotNr; num++)
    for(int i = 0; i < Cats[num].Nsubhalos; i++)
      {
        a_data[count].origin_task  = ThisTask;
        a_data[count].origin_num   = num;
        a_data[count].origin_index = i;
        a_data[count].treeid       = Cats[num].Subhalo[i].TreeID;
        count++;
      }

  mycxxsort_parallel(a_data, a_data + totnsubs, compare_assign_data_treeid_origin_num_origin_task_origin_index, Communicator);

  long long newid = 0;

  for(long long i = 0, treeindex = 0; i < totnsubs; i++)
    {
      if(i == 0 || a_data[i].treeid != a_data[i - 1].treeid)
        treeindex = 0;

      a_data[i].treeindex = treeindex++;

      if(i > 0 && a_data[i].treeid != a_data[i - 1].treeid)
        newid++;

      a_data[i].newtreeid = newid;
    }

  long long *count_list = (long long *)Mem.mymalloc("count_list", NTask * sizeof(long long));
  long long *newid_list = (long long *)Mem.mymalloc("newid_list", NTask * sizeof(long long));
  MPI_Allgather(&totnsubs, sizeof(long long), MPI_BYTE, count_list, sizeof(long long), MPI_BYTE, Communicator);
  MPI_Allgather(&newid, sizeof(long long), MPI_BYTE, newid_list, sizeof(long long), MPI_BYTE, Communicator);

  assign_data *elem_last  = (assign_data *)Mem.mymalloc("elem_last", NTask * sizeof(assign_data));
  assign_data *elem_first = (assign_data *)Mem.mymalloc("elem_first", NTask * sizeof(assign_data));

  /* note: the 0th element is guaranteed to be allocated even on ranks with zero totnsubs */
  MPI_Allgather(&a_data[totnsubs > 0 ? totnsubs - 1 : 0], sizeof(assign_data), MPI_BYTE, elem_last, sizeof(assign_data), MPI_BYTE,
                Communicator);
  MPI_Allgather(&a_data[0], sizeof(assign_data), MPI_BYTE, elem_first, sizeof(assign_data), MPI_BYTE, Communicator);

  if(count_list[ThisTask] > 0)
    {
      long long count_before = 0;
      for(int task = 0; task < ThisTask; task++)
        if(count_list[task] > 0)
          {
            if(a_data[0].treeid == elem_last[task].treeid)
              count_before += (elem_last[task].treeindex + 1);
          }
      if(count_before > 0)
        {
          for(long long i = 0; i < totnsubs; i++)
            {
              if(a_data[i].treeid != a_data[0].treeid)
                break;

              a_data[i].treeindex += count_before;
            }
        }

      long long newidoff = 0;
      for(int task = 0; task < ThisTask; task++)
        newidoff += newid_list[task];

      /* now also need to check whether there are treeid changes that line up with task boundaries. Be aware that that some tasks could
       * be empty */

      int taleft = 0;

      do
        {
          while(taleft < ThisTask && count_list[taleft] == 0)
            taleft++;

          if(taleft < ThisTask)
            {
              int taright = taleft + 1;

              while(count_list[taright] == 0)
                taright++;

              // taright may be at most equal to ThisTask here

              if(elem_last[taleft].treeid != elem_first[taright].treeid)
                newidoff++;

              taleft = taright;
            }
        }
      while(taleft < ThisTask);

      /* now assignn new TreeIDs that form consecutive numbers without gaps */
      for(long long i = 0; i < totnsubs; i++)
        a_data[i].newtreeid += newidoff;
    }

  Mem.myfree(elem_first);
  Mem.myfree(elem_last);
  Mem.myfree(newid_list);
  Mem.myfree(count_list);

  mycxxsort_parallel(a_data, a_data + totnsubs, compare_assign_data_origin_task_origin_num_origin_index, Communicator);

  count = 0;
  for(int num = 0; num <= LastSnapShotNr; num++)
    for(int i = 0; i < Cats[num].Nsubhalos; i++)
      {
        Cats[num].Subhalo[i].TreeIndex = a_data[count].treeindex;
        Cats[num].Subhalo[i].TreeID    = a_data[count].newtreeid;
        count++;
      }

  Mem.myfree(a_data);
}

/*--------------------------------------------------------------------------------------------------------------*/

/* This function puts subhalos from snapshot "num" that are linked by a descendant relation from snapshot "num-1" into the same tree
 * by unifying their tree-id.
 */
int mergertree::halotrees_join_via_descendants(int num)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* note: we take advantage of the fact that the Subgroups are ordered according to their global SubhaloNr for each output number */

  /* the following lists store the minimum and maximum subhalo number to be found on a given processor */
  long long *list_min_subhalonr = (long long *)Mem.mymalloc("list_min_subhalonr", NTask * sizeof(long long));
  long long *list_max_subhalonr = (long long *)Mem.mymalloc("list_max_subhalonr", NTask * sizeof(long long));

  /* this value flags that there are no subhalos on the corresponding processor */
  long long empty = HALONR_MAX;

  MPI_Allgather(Cats[num].Nsubhalos > 0 ? &Cats[num].Subhalo[0].SubhaloNr : &empty, sizeof(long long), MPI_BYTE, list_min_subhalonr,
                sizeof(long long), MPI_BYTE, Communicator);
  MPI_Allgather(Cats[num].Nsubhalos > 0 ? &Cats[num].Subhalo[Cats[num].Nsubhalos - 1].SubhaloNr : &empty, sizeof(long long), MPI_BYTE,
                list_max_subhalonr, sizeof(long long), MPI_BYTE, Communicator);

  int nexport = 0, nimport = 0;

  halotrees_data *import_data = NULL, *export_data = NULL;

  /* for efficiency reasons, we need to traverse the local descendant pointers in increasing order of their target subhalo number, so
   * let's create an auxiliary list for facilitating this.
   */
  descnr_data *sorted_list = (descnr_data *)Mem.mymalloc("sorted_list", Cats[num - 1].Nsubhalos * sizeof(descnr_data));

  for(int i = 0; i < Cats[num - 1].Nsubhalos; i++)
    {
      sorted_list[i].DescSubhaloNr = Cats[num - 1].Descendants[i].DescSubhaloNr;
      sorted_list[i].TreeID        = Cats[num - 1].Subhalo[i].TreeID;
      sorted_list[i].TreeTask      = Cats[num - 1].Subhalo[i].TreeTask;
      sorted_list[i].orig_index    = i;
    }

  mycxxsort(sorted_list, sorted_list + Cats[num - 1].Nsubhalos, compare_sorted_list_descsubhalonr);

  /* for communication bookkeeping reasons, we traverse the counting pattern twice */
  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target = 0;

      for(int i = 0; i < Cats[num - 1].Nsubhalos; i++)
        {
          while(target < NTask - 1 &&
                (list_min_subhalonr[target] == empty || sorted_list[i].DescSubhaloNr > list_max_subhalonr[target]))
            target++;

          if(list_min_subhalonr[target] != empty && sorted_list[i].DescSubhaloNr >= list_min_subhalonr[target] &&
             sorted_list[i].DescSubhaloNr <= list_max_subhalonr[target])
            {
              if(mode == 0)
                Send_count[target]++;
              else
                {
                  int off = Send_offset[target] + Send_count[target]++;

                  export_data[off].loc_index    = sorted_list[i].orig_index;
                  export_data[off].descendantnr = sorted_list[i].DescSubhaloNr;
                  export_data[off].treeid       = sorted_list[i].TreeID;
                  export_data[off].treetask     = sorted_list[i].TreeTask;
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

          export_data = (halotrees_data *)Mem.mymalloc("export_data", nexport * sizeof(halotrees_data));
          import_data = (halotrees_data *)Mem.mymalloc("import_data", nimport * sizeof(halotrees_data));
        }
    }

  /* send data to those target processors that hold the descendant subhalos in order to fetch their tree-ids */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(halotrees_data), MPI_BYTE, recvTask,
                         TAG_DENS_B, &import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(halotrees_data), MPI_BYTE,
                         recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
    }

  /* the collection of incoming data is not necessarily sorted according to descendantnr, so we need to sort it for efficient matching
   */
  for(int i = 0; i < nimport; i++)
    import_data[i].orig_order = i;

  mycxxsort(import_data, import_data + nimport, compare_halotrees_data_descendantnr);

  int changes = 0;

  /* now do the matching */
  for(int i = 0, j = 0; i < Cats[num].Nsubhalos && j < nimport;)
    {
      if(Cats[num].Subhalo[i].SubhaloNr < import_data[j].descendantnr)
        i++;
      else if(Cats[num].Subhalo[i].SubhaloNr > import_data[j].descendantnr)
        j++;
      else
        {
          if(import_data[j].treeid > Cats[num].Subhalo[i].TreeID)
            {
              import_data[j].treeid   = Cats[num].Subhalo[i].TreeID;
              import_data[j].treetask = Cats[num].Subhalo[i].TreeTask;
              changes++;
            }
          else if(import_data[j].treeid < Cats[num].Subhalo[i].TreeID)
            {
              Cats[num].Subhalo[i].TreeID   = import_data[j].treeid;
              Cats[num].Subhalo[i].TreeTask = import_data[j].treetask;
              changes++;
            }
          j++;
        }
    }

  /* reestablish original order */
  mycxxsort(import_data, import_data + nimport, compare_halotrees_data_orig_order);

  /* send the results back */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(halotrees_data), MPI_BYTE, recvTask,
                         TAG_DENS_B, &export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(halotrees_data), MPI_BYTE,
                         recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
    }

  /* now read it out and assign the new treeid/treetask value to the halos in the previous output (which are the progenitors) */
  for(int i = 0; i < nexport; i++)
    {
      if(export_data[i].treeid != HALONR_MAX)
        {
          Cats[num - 1].Subhalo[export_data[i].loc_index].TreeID   = export_data[i].treeid;
          Cats[num - 1].Subhalo[export_data[i].loc_index].TreeTask = export_data[i].treetask;
        }
    }

  Mem.myfree(import_data);
  Mem.myfree(export_data);
  Mem.myfree(sorted_list);
  Mem.myfree(list_max_subhalonr);
  Mem.myfree(list_min_subhalonr);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  return changes;
}

int mergertree::halotrees_join_via_progenitors(int num)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* note: we take advantage of the fact that the Subgroups are ordered according to their global SubhaloNr for each output number */

  /* the following lists store the minimum and maximum subhalo number to be found on a given processor */
  long long *list_min_subhalonr = (long long *)Mem.mymalloc("list_min_subhalonr", NTask * sizeof(long long));
  long long *list_max_subhalonr = (long long *)Mem.mymalloc("list_max_subhalonr", NTask * sizeof(long long));

  /* this value flags that there are no subhalos on the corresponding processor */
  long long empty = HALONR_MAX;

  MPI_Allgather(Cats[num].Nsubhalos > 0 ? &Cats[num].Subhalo[0].SubhaloNr : &empty, sizeof(long long), MPI_BYTE, list_min_subhalonr,
                sizeof(long long), MPI_BYTE, Communicator);
  MPI_Allgather(Cats[num].Nsubhalos > 0 ? &Cats[num].Subhalo[Cats[num].Nsubhalos - 1].SubhaloNr : &empty, sizeof(long long), MPI_BYTE,
                list_max_subhalonr, sizeof(long long), MPI_BYTE, Communicator);

  int nexport = 0, nimport = 0;

  halotrees_data *import_data = NULL, *export_data = NULL;

  /* for efficiency reasons, we need to traverse the local progenitor pointers in increasing order of their target subhalo number, so
   * let's create an auxiliary list for facilitating this.
   */
  prognr_data *sorted_list = (prognr_data *)Mem.mymalloc("sorted_list", Cats[num + 1].Nsubhalos * sizeof(prognr_data));

  for(int i = 0; i < Cats[num + 1].Nsubhalos; i++)
    {
      sorted_list[i].ProgSubhaloNr = Cats[num + 1].Progenitors[i].ProgSubhaloNr;
      sorted_list[i].TreeID        = Cats[num + 1].Subhalo[i].TreeID;
      sorted_list[i].TreeTask      = Cats[num + 1].Subhalo[i].TreeTask;
      sorted_list[i].orig_index    = i;
    }

  mycxxsort(sorted_list, sorted_list + Cats[num + 1].Nsubhalos, compare_sorted_list_progsubhalonr);

  /* for communication bookkeeping reasons, we traverse the counting pattern twice */
  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target = 0;

      for(int i = 0; i < Cats[num + 1].Nsubhalos; i++)
        {
          while(target < NTask - 1 &&
                (list_min_subhalonr[target] == empty || sorted_list[i].ProgSubhaloNr > list_max_subhalonr[target]))
            target++;

          if(list_min_subhalonr[target] != empty && sorted_list[i].ProgSubhaloNr >= list_min_subhalonr[target] &&
             sorted_list[i].ProgSubhaloNr <= list_max_subhalonr[target])
            {
              if(mode == 0)
                Send_count[target]++;
              else
                {
                  int off = Send_offset[target] + Send_count[target]++;

                  export_data[off].loc_index    = sorted_list[i].orig_index;
                  export_data[off].progenitornr = sorted_list[i].ProgSubhaloNr;
                  export_data[off].treeid       = sorted_list[i].TreeID;
                  export_data[off].treetask     = sorted_list[i].TreeTask;
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

          export_data = (halotrees_data *)Mem.mymalloc("export_data", nexport * sizeof(halotrees_data));
          import_data = (halotrees_data *)Mem.mymalloc("import_data", nimport * sizeof(halotrees_data));
        }
    }

  /* send data to those target processors that hold the descendant subhalos in order to fetch their tree-ids */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(halotrees_data), MPI_BYTE, recvTask,
                         TAG_DENS_B, &import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(halotrees_data), MPI_BYTE,
                         recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
    }

  /* the collection of incoming data is not necessarily sorted according to descendantnr, so we need to sort it for efficient
   * matching
   */
  for(int i = 0; i < nimport; i++)
    import_data[i].orig_order = i;

  mycxxsort(import_data, import_data + nimport, compare_halotrees_data_progenitornr);

  int changes = 0;

  /* now do the matching */
  for(int i = 0, j = 0; i < Cats[num].Nsubhalos && j < nimport;)
    {
      if(Cats[num].Subhalo[i].SubhaloNr < import_data[j].progenitornr)
        i++;
      else if(Cats[num].Subhalo[i].SubhaloNr > import_data[j].progenitornr)
        j++;
      else
        {
          if(import_data[j].treeid > Cats[num].Subhalo[i].TreeID)
            {
              import_data[j].treeid   = Cats[num].Subhalo[i].TreeID;
              import_data[j].treetask = Cats[num].Subhalo[i].TreeTask;
              changes++;
            }
          else if(import_data[j].treeid < Cats[num].Subhalo[i].TreeID)
            {
              Cats[num].Subhalo[i].TreeID   = import_data[j].treeid;
              Cats[num].Subhalo[i].TreeTask = import_data[j].treetask;
              changes++;
            }
          j++;
        }
    }

  /* reestablish original order */
  mycxxsort(import_data, import_data + nimport, compare_halotrees_data_orig_order);

  /* send the results back */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(halotrees_data), MPI_BYTE, recvTask,
                         TAG_DENS_B, &export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(halotrees_data), MPI_BYTE,
                         recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
    }

  /* now read it out and assign the new treeid/treetask value to the halos in the previous output (which are the progenitors) */
  for(int i = 0; i < nexport; i++)
    {
      if(export_data[i].treeid != HALONR_MAX)
        {
          Cats[num + 1].Subhalo[export_data[i].loc_index].TreeID   = export_data[i].treeid;
          Cats[num + 1].Subhalo[export_data[i].loc_index].TreeTask = export_data[i].treetask;
        }
    }

  Mem.myfree(import_data);
  Mem.myfree(export_data);
  Mem.myfree(sorted_list);
  Mem.myfree(list_max_subhalonr);
  Mem.myfree(list_min_subhalonr);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  return changes;
}

void mergertree::halotrees_determine_mainprogenitor(void)
{
  for(int num = 1; num <= LastSnapShotNr; num++)
    for(int i = 0; i < Cats[num].Nsubhalos; i++)
      Cats[num].Progenitors[i].MainProgSubhaloNr = HALONR_MAX;

  /* initialize the maximum branch length field */
  for(int num = 0; num <= LastSnapShotNr; num++)
    for(int i = 0; i < Cats[num].Nsubhalos; i++)
      Cats[num].SubExt[i].MaxLenProgBranch = Cats[num].Subhalo[i].Len;

  /* propagate the branch length of subhalos from snapshot "num-1" to those in snapshot "num" via the
   * descendant information */
  for(int num = 1; num <= LastSnapShotNr; num++)
    halotrees_propagate_max_branch_length_descendants(num);

  /* propagate the branch length of subhalos from snapshot "num-1" to those in snapshot "num" via the
   * progenitor information */
  /* we disable this as it could cause a spurious jump to a larger progenitor that shed a few particles to a newly formed group */
  // for(int num = 1; num <= LastSnapShotNr; num++)
  //   halotrees_propagate_max_branch_length_progenitors(num);

  mpi_printf("MERGERTREE: determination of main progenitor branch done\n");
}

/* This function determines the maximum branch size of subhalos from snapshot "num-1" to those snapshot "num" via the
 * descendant information
 */
void mergertree::halotrees_propagate_max_branch_length_descendants(int num)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* note: we take advantage of the fact that the Subgroups are ordered according to their global SubhaloNr for each output number */

  /* the following lists store the minimum and maximum subhalo number to be found on a given processor */
  long long *list_min_subhalonr = (long long *)Mem.mymalloc("list_min_subhalonr", NTask * sizeof(long long));
  long long *list_max_subhalonr = (long long *)Mem.mymalloc("list_max_subhalonr", NTask * sizeof(long long));

  /* this value flags that there are no subhalos on the corresponding processor */
  long long empty = HALONR_MAX;

  MPI_Allgather(Cats[num].Nsubhalos > 0 ? &Cats[num].Subhalo[0].SubhaloNr : &empty, sizeof(long long), MPI_BYTE, list_min_subhalonr,
                sizeof(long long), MPI_BYTE, Communicator);
  MPI_Allgather(Cats[num].Nsubhalos > 0 ? &Cats[num].Subhalo[Cats[num].Nsubhalos - 1].SubhaloNr : &empty, sizeof(long long), MPI_BYTE,
                list_max_subhalonr, sizeof(long long), MPI_BYTE, Communicator);

  int nexport = 0, nimport = 0;

  halotrees_propagate_data *import_data = NULL, *export_data = NULL;

  /* for efficiency reasons, we need to traverse the local descendant pointers in increasing order of their target subhalo number, so
   * let's create an auxiliary list for facilitating this.
   */
  halotrees_propagate_data *sorted_list =
      (halotrees_propagate_data *)Mem.mymalloc("sorted_list", Cats[num - 1].Nsubhalos * sizeof(halotrees_propagate_data));

  for(int i = 0; i < Cats[num - 1].Nsubhalos; i++)
    {
      sorted_list[i].DescSubhaloNr    = Cats[num - 1].Descendants[i].DescSubhaloNr;
      sorted_list[i].SubhaloNr        = Cats[num - 1].Descendants[i].PrevSubhaloNr;
      sorted_list[i].MaxLenProgBranch = Cats[num - 1].SubExt[i].MaxLenProgBranch;
    }

  mycxxsort(sorted_list, sorted_list + Cats[num - 1].Nsubhalos, compare_halotrees_propagate_data_DescSubhaloNr);

  /* for communication bookkeeping reasons, we traverse the counting pattern twice */
  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target = 0;

      for(int i = 0; i < Cats[num - 1].Nsubhalos; i++)
        {
          while(target < NTask - 1 &&
                (list_min_subhalonr[target] == empty || sorted_list[i].DescSubhaloNr > list_max_subhalonr[target]))
            target++;

          if(list_min_subhalonr[target] != empty && sorted_list[i].DescSubhaloNr >= list_min_subhalonr[target] &&
             sorted_list[i].DescSubhaloNr <= list_max_subhalonr[target])
            {
              if(mode == 0)
                Send_count[target]++;
              else
                {
                  int off = Send_offset[target] + Send_count[target]++;

                  export_data[off].DescSubhaloNr    = sorted_list[i].DescSubhaloNr;
                  export_data[off].SubhaloNr        = sorted_list[i].SubhaloNr;
                  export_data[off].MaxLenProgBranch = sorted_list[i].MaxLenProgBranch;
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

          export_data = (halotrees_propagate_data *)Mem.mymalloc("export_data", nexport * sizeof(halotrees_propagate_data));
          import_data = (halotrees_propagate_data *)Mem.mymalloc("import_data", nimport * sizeof(halotrees_propagate_data));
        }
    }

  /* send data to those target processors that hold the descendant subhalos in order to fetch their tree-ids */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(halotrees_propagate_data), MPI_BYTE,
                         recvTask, TAG_DENS_B, &import_data[Recv_offset[recvTask]],
                         Recv_count[recvTask] * sizeof(halotrees_propagate_data), MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                         MPI_STATUS_IGNORE);
    }

  /* the collection of incoming data is not necessarily sorted according to DescSubhaloNr, so we need to sort it for efficient
   * matching
   */
  mycxxsort(import_data, import_data + nimport, compare_halotrees_propagate_data_DescSubhaloNr);

  /* now do the matching */
  for(int i = 0, j = 0; i < Cats[num].Nsubhalos && j < nimport;)
    {
      if(Cats[num].Subhalo[i].SubhaloNr < import_data[j].DescSubhaloNr)
        i++;
      else if(Cats[num].Subhalo[i].SubhaloNr > import_data[j].DescSubhaloNr)
        j++;
      else
        {
          if(import_data[j].MaxLenProgBranch > Cats[num].SubExt[i].MaxLenProgBranch - Cats[num].Subhalo[i].Len)
            {
              Cats[num].SubExt[i].MaxLenProgBranch       = import_data[j].MaxLenProgBranch + Cats[num].Subhalo[i].Len;
              Cats[num].Progenitors[i].MainProgSubhaloNr = import_data[j].SubhaloNr;
            }

          j++;
        }
    }

  Mem.myfree(import_data);
  Mem.myfree(export_data);
  Mem.myfree(sorted_list);
  Mem.myfree(list_max_subhalonr);
  Mem.myfree(list_min_subhalonr);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

/* This function determines the maximum branch size of subhalos from snapshot "num-1" to those snapshot "num" via the
 * progenitor information
 */
void mergertree::halotrees_propagate_max_branch_length_progenitors(int num)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* note: we take advantage of the fact that the Subgroups are ordered according to their global SubhaloNr for each output number */

  /* the following lists store the minimum and maximum subhalo number to be found on a given processor */
  long long *list_min_subhalonr = (long long *)Mem.mymalloc("list_min_subhalonr", NTask * sizeof(long long));
  long long *list_max_subhalonr = (long long *)Mem.mymalloc("list_max_subhalonr", NTask * sizeof(long long));

  /* this value flags that there are no subhalos on the corresponding processor */
  long long empty = HALONR_MAX;

  MPI_Allgather(Cats[num - 1].Nsubhalos > 0 ? &Cats[num - 1].Subhalo[0].SubhaloNr : &empty, sizeof(long long), MPI_BYTE,
                list_min_subhalonr, sizeof(long long), MPI_BYTE, Communicator);
  MPI_Allgather(Cats[num - 1].Nsubhalos > 0 ? &Cats[num - 1].Subhalo[Cats[num - 1].Nsubhalos - 1].SubhaloNr : &empty,
                sizeof(long long), MPI_BYTE, list_max_subhalonr, sizeof(long long), MPI_BYTE, Communicator);

  int nexport = 0, nimport = 0;

  halotrees_propagate_data *import_data = NULL, *export_data = NULL;

  /* for efficiency reasons, we need to traverse the local descendant pointers in increasing order of their target subhalo number, so
   * let's create an auxiliary list for facilitating this.
   */
  halotrees_propagate_data *sorted_list =
      (halotrees_propagate_data *)Mem.mymalloc("sorted_list", Cats[num].Nsubhalos * sizeof(halotrees_propagate_data));

  for(int i = 0; i < Cats[num].Nsubhalos; i++)
    {
      sorted_list[i].ProgSubhaloNr = Cats[num].Progenitors[i].ProgSubhaloNr;
      sorted_list[i].index         = i;
    }

  mycxxsort(sorted_list, sorted_list + Cats[num].Nsubhalos, compare_halotrees_propagate_data_ProgSubhaloNr);

  /* for communication bookkeeping reasons, we traverse the counting pattern twice */
  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target = 0;

      for(int i = 0; i < Cats[num].Nsubhalos; i++)
        {
          while(target < NTask - 1 &&
                (list_min_subhalonr[target] == empty || sorted_list[i].ProgSubhaloNr > list_max_subhalonr[target]))
            target++;

          if(list_min_subhalonr[target] != empty && sorted_list[i].ProgSubhaloNr >= list_min_subhalonr[target] &&
             sorted_list[i].ProgSubhaloNr <= list_max_subhalonr[target])
            {
              if(mode == 0)
                Send_count[target]++;
              else
                {
                  int off = Send_offset[target] + Send_count[target]++;

                  export_data[off].ProgSubhaloNr = sorted_list[i].ProgSubhaloNr;
                  export_data[off].index         = sorted_list[i].index;
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

          export_data = (halotrees_propagate_data *)Mem.mymalloc("export_data", nexport * sizeof(halotrees_propagate_data));
          import_data = (halotrees_propagate_data *)Mem.mymalloc("import_data", nimport * sizeof(halotrees_propagate_data));
        }
    }

  /* send data to those target processors that hold the descendant subhalos in order to fetch their tree-ids */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(halotrees_propagate_data), MPI_BYTE,
                         recvTask, TAG_DENS_B, &import_data[Recv_offset[recvTask]],
                         Recv_count[recvTask] * sizeof(halotrees_propagate_data), MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                         MPI_STATUS_IGNORE);
    }

  for(int i = 0; i < nimport; i++)
    import_data[i].orig_order = i;

  /* the collection of incoming data is not necessarily sorted according to ProgSubhaloNr, so we need to sort it for efficient
   * matching
   */
  mycxxsort(import_data, import_data + nimport, compare_halotrees_propagate_data_ProgSubhaloNr);

  /* now do the matching */
  for(int i = 0, j = 0; i < Cats[num - 1].Nsubhalos && j < nimport;)
    {
      if(Cats[num - 1].Subhalo[i].SubhaloNr < import_data[j].ProgSubhaloNr)
        i++;
      else if(Cats[num - 1].Subhalo[i].SubhaloNr > import_data[j].ProgSubhaloNr)
        j++;
      else
        {
          import_data[j].MaxLenProgBranch = Cats[num - 1].SubExt[i].MaxLenProgBranch;
          import_data[j].SubhaloNr        = Cats[num - 1].Subhalo[i].SubhaloNr;

          if(Cats[num - 1].Subhalo[i].SubhaloNr != import_data[j].ProgSubhaloNr)
            Terminate("Cats[num - 1].Subhalo[i].SubhaloNr != import_data[j].ProgSubhaloNr)");

          j++;
        }
    }

  /* reestablish original order */
  mycxxsort(import_data, import_data + nimport, compare_halotrees_propagate_data_orig_order);

  /* send data back */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&import_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(halotrees_propagate_data), MPI_BYTE,
                         recvTask, TAG_DENS_B, &export_data[Send_offset[recvTask]],
                         Send_count[recvTask] * sizeof(halotrees_propagate_data), MPI_BYTE, recvTask, TAG_DENS_B, Communicator,
                         MPI_STATUS_IGNORE);
    }

  /* now read it out and assign the new treeid/treetask value to the halos in the previous output (which are the progenitors) */
  for(int i = 0; i < nexport; i++)
    {
      int q = export_data[i].index;

      if(export_data[i].MaxLenProgBranch > Cats[num].SubExt[q].MaxLenProgBranch - Cats[num].Subhalo[q].Len)
        {
          Cats[num].SubExt[q].MaxLenProgBranch       = export_data[i].MaxLenProgBranch + Cats[num].Subhalo[q].Len;
          Cats[num].Progenitors[q].MainProgSubhaloNr = export_data[i].SubhaloNr;
        }
    }

  Mem.myfree(import_data);
  Mem.myfree(export_data);
  Mem.myfree(sorted_list);
  Mem.myfree(list_max_subhalonr);
  Mem.myfree(list_min_subhalonr);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

#endif
