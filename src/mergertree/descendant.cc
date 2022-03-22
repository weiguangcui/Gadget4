/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  descendant.cc
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
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mergertree/io_descendant.h"
#include "../mergertree/io_progenitors.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

/* This function allocates and fills the "Descendants" array, which gives the number of the descendant subhalo in the newly created
 * subhalo catalogue for every subhalo from the previous catalogue.
 */
void mergertree::mergertree_determine_descendants_postproc(int num)
{
  Descendants = (desc_list *)Mem.mymalloc_movable(&Descendants, "Descendants", (PrevNsubhalos + 1) * sizeof(desc_list));
  Progenitors = (prog_list *)Mem.mymalloc_movable(&Progenitors, "Progenitors", (CurrNsubhalos + 1) * sizeof(prog_list));

  Sp->NumPart = MtrP_NumPart;

  /* allocate a work structure */
  desc = (desc_partdata *)Mem.mymalloc_movable(&desc, "desc", sizeof(desc_partdata) * Sp->NumPart);

  /* let's fill in some relevant data into the work-structure */
  for(int i = 0; i < Sp->NumPart; i++)
    {
      desc[i].CurrSubhaloNr.set(MtrP[i].SubhaloNr);
      desc[i].CurrRankInSubhalo.set(MtrP[i].RankInSubhalo);

      desc[i].PrevSubhaloNr.set(MtrP[i].PrevSubhaloNr);
      desc[i].PrevRankInSubhalo.set(MtrP[i].PrevRankInSubhalo);

      if(desc[i].PrevSubhaloNr.get() >= PrevTotNsubhalos && desc[i].PrevSubhaloNr.get() != HALONR_MAX)
        Terminate("strange: i=%d  desc[i].PrevSubhaloNr=%lld  PrevTotNsubhalos=%lld", i, (long long)desc[i].PrevSubhaloNr.get(),
                  (long long)PrevTotNsubhalos);
    }

  mergertree_determine_descendants(num);

  mpi_printf("DONE!\n");
}

/* This function allocates and fills the "Descendants" array, which gives the number of the descendant subhalo in the newly created
 * subhalo catalogue for every subhalo from the previous catalogue.
 */
void mergertree::mergertree_determine_descendants_on_the_fly(int num)
{
  if(num == 0)  // for the first output, we don't yet have anything to link
    return;

  Descendants = (desc_list *)Mem.mymalloc_movable(&Descendants, "Descendants", (PrevNsubhalos + 1) * sizeof(desc_list));
  Progenitors = (prog_list *)Mem.mymalloc_movable(&Progenitors, "Progenitors", (CurrNsubhalos + 1) * sizeof(prog_list));

  /* allocate a work structure */
  desc = (desc_partdata *)Mem.mymalloc_movable(&desc, "desc", sizeof(desc_partdata) * Sp->NumPart);

  /* Let's fill in some relevant data into a work-structure defined for every particle
   * For each particle, we know the new and the old subhalo number, as well as its rank in the previous subhalo
   */
  for(int i = 0; i < Sp->NumPart; i++)
    {
      desc[i].CurrSubhaloNr     = Sp->PS[i].SubhaloNr;
      desc[i].CurrRankInSubhalo = Sp->PS[i].RankInSubhalo;

      desc[i].PrevSubhaloNr     = Sp->P[i].PrevSubhaloNr;
      desc[i].PrevRankInSubhalo = Sp->P[i].PrevRankInSubhalo;

      if(desc[i].PrevSubhaloNr.get() >= PrevTotNsubhalos && desc[i].PrevSubhaloNr.get() != HALONR_MAX)
        Terminate("strange: i=%d  desc[i].PrevSubhaloNr=%lld  PrevTotNsubhalos=%lld", i, (long long)desc[i].PrevSubhaloNr.get(),
                  (long long)PrevTotNsubhalos);
    }

  mergertree_determine_descendants(num);

  Mem.myfree(desc);
  Mem.myfree(Progenitors);
  Mem.myfree(Descendants);
}

void mergertree::mergertree_determine_descendants(int num)
{
  /* determine matching pieces of subhalos and their mutual scores */
  int nmatch = mergertree_find_matching_segments_and_scores();

  /* select the main descendants */
  mergertree_select_maximum_score_descendants(nmatch);

  /* select the main progenitors */
  mergertree_select_maximum_score_progenitors(nmatch);

  /* let's determine the next-progenitor fields, which chain up those subhalos that have the same descendant */
  mergertree_chain_up_progenitors_with_same_descendant();

  /* set first progenitor field for chaining up those subhalos that have the same descendant  */
  mergertree_set_first_progenitor_with_same_descendant();

  /* let's determine the next-descendant fields, which chain up those subhalos that have the same progenitor */
  mergertree_chain_up_descendants_with_same_progenitor();

  /* set first descendant field for chaining up those subhalos that have the same progenitor */
  mergertree_set_first_descendant_with_same_progenitor();

  /**** write stuff to files *****/

  descendant_io Desc(this, this->Communicator, All.SnapFormat); /* get an I/O object */
  Desc.mergertree_save_descendants(num);

  progenitors_io Prog(this, this->Communicator, All.SnapFormat);
  Prog.mergertree_save_progenitors(num);

  /* calculate and output some statistics to characterize linking */

  int count_links = 0;
  for(int i = 0; i < PrevNsubhalos; i++)
    if(Descendants[i].DescSubhaloNr != HALONR_MAX)
      count_links++;

  long long tot_count_links;
  sumup_large_ints(1, &count_links, &tot_count_links, Communicator);

  mpi_printf("MERGERTREE: Was able to identify descendants for %lld out of %lld subhalos, i.e. for a fraction of %g\n",
             tot_count_links, (long long)PrevTotNsubhalos, tot_count_links / (PrevTotNsubhalos + SMALLNUM));

  int count_splits = 0;
  for(int i = 0; i < CurrNsubhalos; i++)
    if(Progenitors[i].NextDescSubhaloNr != HALONR_MAX)
      count_splits++;

  long long tot_count_splits;
  sumup_large_ints(1, &count_splits, &tot_count_splits, Communicator);

  mpi_printf("MERGERTREE: We have found secondary descendants for %lld halos out of %lld subhalos, i.e. for a fraction of %g\n",
             tot_count_splits, (long long)CurrTotNsubhalos, tot_count_splits / (CurrTotNsubhalos + SMALLNUM));
}

int mergertree::mergertree_find_matching_segments_and_scores(void)
{
  /*  let's eliminate unbound particles from the work list */
  int nmatch = Sp->NumPart;

  for(int i = 0; i < nmatch; i++)
    if(desc[i].CurrSubhaloNr.get() == HALONR_MAX || desc[i].PrevSubhaloNr.get() == HALONR_MAX)
      {
        desc[i] = desc[nmatch - 1];
        nmatch--;
        i--;
      }

  /* let's do the scoring */
  for(int i = 0; i < nmatch; i++)
    {
      if(desc[i].PrevRankInSubhalo.get() < NUM_MOST_BOUND_PARTICLES_USED_FOR_TRACKING)
        desc[i].DescScore = 1.0 / (1 + pow(desc[i].PrevRankInSubhalo.get(), 0.5));
      else
        desc[i].DescScore = 0;

      if(desc[i].CurrRankInSubhalo.get() < NUM_MOST_BOUND_PARTICLES_USED_FOR_TRACKING)
        desc[i].ProgScore = 1.0 / (1 + pow(desc[i].CurrRankInSubhalo.get(), 0.5));
      else
        desc[i].ProgScore = 0;
    }

  /* Now we sort the list such that the old subhalos are grouped together, and the new subhalo numbers are consecutive within them
   */
  mycxxsort_parallel(desc, desc + nmatch, mergertree_compare_PrevSubNr_NewSubNr, Communicator);

  /* eliminate duplicate matched pairs on local processor and sum up the scores
   */
  int start = 0;
  int count = 0;

  if(nmatch > 0)
    count = 1;

  for(int i = 1; i < nmatch; i++)
    if(desc[i].PrevSubhaloNr == desc[start].PrevSubhaloNr && desc[i].CurrSubhaloNr == desc[start].CurrSubhaloNr)
      {
        desc[start].DescScore += desc[i].DescScore;
        desc[start].ProgScore += desc[i].ProgScore;
      }
    else
      {
        desc[count] = desc[i];
        start       = count;
        count++;
      }

  nmatch = count;

  /* now we consolidate duplicate matched pairs on different processors. The list is still ordered, but there could be gaps
   * (i.e. processors with only one or potentially zero entries)
   */

  /* obtain last and first element of each processor, and the counts from each processor */
  desc_partdata *desc_first = (desc_partdata *)Mem.mymalloc("desc_first", NTask * sizeof(desc_partdata));
  desc_partdata *desc_last  = (desc_partdata *)Mem.mymalloc("desc_last", NTask * sizeof(desc_partdata));
  int *nmatch_list          = (int *)Mem.mymalloc("nmatch_list", NTask * sizeof(int));

  MPI_Allgather(&desc[0], sizeof(desc_partdata), MPI_BYTE, desc_first, sizeof(desc_partdata), MPI_BYTE,
                Communicator); /* note: the 0th element is guaranteed to be allocated */
  MPI_Allgather(&desc[nmatch > 0 ? nmatch - 1 : 0], sizeof(desc_partdata), MPI_BYTE, desc_last, sizeof(desc_partdata), MPI_BYTE,
                Communicator);
  MPI_Allgather(&nmatch, 1, MPI_INT, nmatch_list, 1, MPI_INT, Communicator);

  /* We go from the back through the tasks, and for the first element in each case, we move it ahead in the list as far as possible.
   * Eliminated items are marked through DescScore = -1, and are then filtered out later.
   */
  for(int i = NTask - 1; i > 0; i--)
    if(nmatch_list[i] > 0)
      {
        int target = -1;

        for(int j = i - 1; j >= 0; j--)
          {
            if(nmatch_list[j] > 0)
              {
                if(nmatch_list[j] > 1)
                  {
                    if(desc_last[j].PrevSubhaloNr == desc_first[i].PrevSubhaloNr &&
                       desc_last[j].CurrSubhaloNr == desc_first[i].CurrSubhaloNr)
                      target = j;
                  }
                else /* nmatch_list[j] == 1 */
                  {
                    if(desc_first[j].PrevSubhaloNr == desc_first[i].PrevSubhaloNr &&
                       desc_first[j].CurrSubhaloNr == desc_first[i].CurrSubhaloNr)
                      target = j;
                  }

                break;
              }
          }

        if(target >= 0)
          {
            if(nmatch_list[target] > 1)
              {
                desc_last[target].DescScore += desc_first[i].DescScore;
                desc_last[target].ProgScore += desc_first[i].ProgScore;

                if(ThisTask == target)
                  {
                    desc[nmatch - 1].DescScore += desc_first[i].DescScore;
                    desc[nmatch - 1].ProgScore += desc_first[i].ProgScore;
                  }
              }
            else
              {
                desc_first[target].DescScore += desc_first[i].DescScore;
                desc_first[target].ProgScore += desc_first[i].ProgScore;

                if(ThisTask == target)
                  {
                    desc[0].DescScore += desc_first[i].DescScore;
                    desc[0].ProgScore += desc_first[i].ProgScore;
                  }
              }

            desc_first[i].DescScore = -1;
            desc_first[i].ProgScore = -1;

            if(ThisTask == i)
              {
                desc[0].DescScore = -1;
                desc[0].ProgScore = -1;
              }
          }
      }

  Mem.myfree(nmatch_list);
  Mem.myfree(desc_last);
  Mem.myfree(desc_first);

  /* now eliminate the ones with negative score
   */
  if(nmatch > 0 && desc[0].DescScore < 0)
    {
      nmatch--;
      memmove(desc, desc + 1, nmatch * sizeof(desc_partdata));
    }

  return nmatch;
}

void mergertree::mergertree_chain_up_descendants_with_same_progenitor(void)
{
  /* sort by progenitor to bring equal ones next to each other */
  mycxxsort_parallel(Progenitors, Progenitors + CurrNsubhalos, mergertree_compare_ProgSubhaloNr, Communicator);

  prog_list *elem_first = (prog_list *)Mem.mymalloc("elem_first", NTask * sizeof(prog_list));
  prog_list *elem_last  = (prog_list *)Mem.mymalloc("elem_last", NTask * sizeof(prog_list));

  /* note: the 0th element is guaranteed to be allocated even on ranks with zero CurrNsubhalos */
  MPI_Allgather(&Progenitors[0], sizeof(prog_list), MPI_BYTE, elem_first, sizeof(prog_list), MPI_BYTE, Communicator);
  MPI_Allgather(&Progenitors[CurrNsubhalos > 0 ? CurrNsubhalos - 1 : 0], sizeof(prog_list), MPI_BYTE, elem_last, sizeof(prog_list),
                MPI_BYTE, Communicator);

  /* get list of the subhalo count on each processor, and the cumulative number stored before */
  int *tab_CurrNsubhalos = (int *)Mem.mymalloc("tab_CurrNsubhalos", sizeof(int) * NTask);
  MPI_Allgather(&CurrNsubhalos, 1, MPI_INT, tab_CurrNsubhalos, 1, MPI_INT, Communicator);

  int next_task = -1;
  for(int i = ThisTask + 1; i < NTask; i++)
    if(tab_CurrNsubhalos[i] > 0)
      {
        next_task = i;
        break;
      }

  int prev_task = -1;
  for(int i = ThisTask - 1; i >= 0; i--)
    if(tab_CurrNsubhalos[i] > 0)
      {
        prev_task = i;
        break;
      }

  for(int i = 0; i < CurrNsubhalos; i++)
    {
      if(i < CurrNsubhalos - 1)
        {
          if(Progenitors[i].ProgSubhaloNr == Progenitors[i + 1].ProgSubhaloNr && Progenitors[i].ProgSubhaloNr != HALONR_MAX)
            Progenitors[i].NextDescSubhaloNr = Progenitors[i + 1].SubhaloNr;
          else
            Progenitors[i].NextDescSubhaloNr = HALONR_MAX;
        }
      else
        {
          if(next_task >= 0 && Progenitors[i].ProgSubhaloNr == elem_first[next_task].ProgSubhaloNr &&
             Progenitors[i].ProgSubhaloNr != HALONR_MAX)
            Progenitors[i].NextDescSubhaloNr = elem_first[next_task].SubhaloNr;
          else
            Progenitors[i].NextDescSubhaloNr = HALONR_MAX;
        }

      if(i > 0)
        {
          if(Progenitors[i].ProgSubhaloNr != Progenitors[i - 1].ProgSubhaloNr && Progenitors[i].ProgSubhaloNr != HALONR_MAX)
            Progenitors[i].FirstDescFlag = 1; /* flags the first progenitors */
          else
            Progenitors[i].FirstDescFlag = 0;
        }
      else
        {
          if(Progenitors[i].ProgSubhaloNr != HALONR_MAX &&
             (ThisTask == 0 || (prev_task >= 0 && Progenitors[i].ProgSubhaloNr != elem_last[prev_task].ProgSubhaloNr)))
            Progenitors[i].FirstDescFlag = 1; /* flags the first progenitors */
          else
            Progenitors[i].FirstDescFlag = 0;
        }
    }

  Mem.myfree(tab_CurrNsubhalos);
  Mem.myfree(elem_last);
  Mem.myfree(elem_first);

  /* bring back into original order */
  mycxxsort_parallel(Progenitors, Progenitors + CurrNsubhalos, mergertree_compare_SubhaloNr, Communicator);
}

/* This function determines the next progenitor field, which chains up those subhalos that have the same descendant */
void mergertree::mergertree_chain_up_progenitors_with_same_descendant(void)
{
  /* sort by descendant to bring equal ones next to each other */
  mycxxsort_parallel(Descendants, Descendants + PrevNsubhalos, mergertree_compare_DescSubhaloNr, Communicator);

  desc_list *elem_first = (desc_list *)Mem.mymalloc("elem_first", NTask * sizeof(desc_list));
  desc_list *elem_last  = (desc_list *)Mem.mymalloc("elem_last", NTask * sizeof(desc_list));

  /* note: the 0th element is guaranteed to be allocated even on ranks with zero PrevNsubhalos */
  MPI_Allgather(&Descendants[0], sizeof(desc_list), MPI_BYTE, elem_first, sizeof(desc_list), MPI_BYTE, Communicator);
  MPI_Allgather(&Descendants[PrevNsubhalos > 0 ? PrevNsubhalos - 1 : 0], sizeof(desc_list), MPI_BYTE, elem_last, sizeof(desc_list),
                MPI_BYTE, Communicator);

  /* get list of the subhalo count on each processor, and the cumulative number stored before */
  int *tab_PrevNsubhalos = (int *)Mem.mymalloc("tab_PrevNsubhalos", sizeof(int) * NTask);
  MPI_Allgather(&PrevNsubhalos, 1, MPI_INT, tab_PrevNsubhalos, 1, MPI_INT, Communicator);

  int next_task = -1;
  for(int i = ThisTask + 1; i < NTask; i++)
    if(tab_PrevNsubhalos[i] > 0)
      {
        next_task = i;
        break;
      }

  int prev_task = -1;
  for(int i = ThisTask - 1; i >= 0; i--)
    if(tab_PrevNsubhalos[i] > 0)
      {
        prev_task = i;
        break;
      }

  for(int i = 0; i < PrevNsubhalos; i++)
    {
      if(i < PrevNsubhalos - 1)
        {
          if(Descendants[i].DescSubhaloNr == Descendants[i + 1].DescSubhaloNr && Descendants[i].DescSubhaloNr != HALONR_MAX)
            Descendants[i].NextProgSubhaloNr = Descendants[i + 1].PrevSubhaloNr;
          else
            Descendants[i].NextProgSubhaloNr = HALONR_MAX;
        }
      else
        {
          if(next_task >= 0 && Descendants[i].DescSubhaloNr == elem_first[next_task].DescSubhaloNr &&
             Descendants[i].DescSubhaloNr != HALONR_MAX)
            Descendants[i].NextProgSubhaloNr = elem_first[next_task].PrevSubhaloNr;
          else
            Descendants[i].NextProgSubhaloNr = HALONR_MAX;
        }

      if(i > 0)
        {
          if(Descendants[i].DescSubhaloNr != Descendants[i - 1].DescSubhaloNr && Descendants[i].DescSubhaloNr != HALONR_MAX)
            Descendants[i].FirstProgFlag = 1; /* flags the first progenitors */
          else
            Descendants[i].FirstProgFlag = 0;
        }
      else
        {
          if(Descendants[i].DescSubhaloNr != HALONR_MAX &&
             (ThisTask == 0 || (prev_task >= 0 && Descendants[i].DescSubhaloNr != elem_last[prev_task].DescSubhaloNr)))
            Descendants[i].FirstProgFlag = 1; /* flags the first progenitors */
          else
            Descendants[i].FirstProgFlag = 0;
        }
    }

  Mem.myfree(tab_PrevNsubhalos);
  Mem.myfree(elem_last);
  Mem.myfree(elem_first);

  /* bring back into original order */
  mycxxsort_parallel(Descendants, Descendants + PrevNsubhalos, mergertree_compare_PrevSubhaloNr, Communicator);
}

/******** set first progenitor field  ****/
void mergertree::mergertree_set_first_progenitor_with_same_descendant(void)
{
  /* sort by descendant to bring equal ones next to each other */
  mycxxsort_parallel(Descendants, Descendants + PrevNsubhalos, mergertree_compare_DescSubhaloNr, Communicator);

  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  int *tab_CurrNsubhalos = (int *)Mem.mymalloc("tab_CurrNsubhalos", sizeof(int) * NTask);
  MPI_Allgather(&CurrNsubhalos, 1, MPI_INT, tab_CurrNsubhalos, 1, MPI_INT, Communicator);

  long long cumul_currnsubhalos = 0;
  for(int i = 0; i < ThisTask; i++)
    cumul_currnsubhalos += tab_CurrNsubhalos[i];

  struct pair_data
  {
    long long subhalonr;
    long long firstprognr;
  };

  pair_data *send_data = NULL;
  pair_data *recv_data = NULL;
  int nexport = 0, nimport = 0;

  for(int mode = 0; mode < 2; mode++)  // go through this twice to simplify bookkeeping
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int task        = 0;
      long long first = 0;

      for(int i = 0; i < PrevNsubhalos; i++)
        {
          if(Descendants[i].FirstProgFlag && Descendants[i].DescSubhaloNr != HALONR_MAX)
            {
              while(task < NTask - 1 && Descendants[i].DescSubhaloNr >= first + tab_CurrNsubhalos[task])
                {
                  first += tab_CurrNsubhalos[task];
                  task++;
                }

              if(mode == 0)
                Send_count[task]++;
              else
                {
                  int off = Send_offset[task] + Send_count[task]++;

                  send_data[off].subhalonr   = Descendants[i].DescSubhaloNr;
                  send_data[off].firstprognr = Descendants[i].PrevSubhaloNr;
                }
            }
        }

      if(mode == 0)  // prepare offset tables
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

          Recv_offset[0] = 0;
          Send_offset[0] = 0;

          for(int j = 0; j < NTask; j++)
            {
              nexport += Send_count[j];
              nimport += Recv_count[j];

              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          send_data = (pair_data *)Mem.mymalloc("pair_data", nexport * sizeof(pair_data));
          recv_data = (pair_data *)Mem.mymalloc("pair_data", nimport * sizeof(pair_data));
        }
    }

  /* exchange data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            myMPI_Sendrecv(&send_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(pair_data), MPI_BYTE, recvTask, TAG_DENS_A,
                         &recv_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(pair_data), MPI_BYTE, recvTask, TAG_DENS_A,
                         Communicator, MPI_STATUS_IGNORE);
        }
    }

  for(int i = 0; i < CurrNsubhalos; i++)
    Progenitors[i].FirstProgSubhaloNr = HALONR_MAX;

  for(int i = 0; i < nimport; i++)
    {
      int off = recv_data[i].subhalonr - cumul_currnsubhalos;

      if(off < 0 || off >= CurrNsubhalos)
        Terminate("off = %d  CurrNsubhalos = %d", off, CurrNsubhalos);

      Progenitors[off].FirstProgSubhaloNr = recv_data[i].firstprognr;
    }

  Mem.myfree(recv_data);
  Mem.myfree(send_data);

  Mem.myfree(tab_CurrNsubhalos);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  /* bring back into original order */
  mycxxsort_parallel(Descendants, Descendants + PrevNsubhalos, mergertree_compare_PrevSubhaloNr, Communicator);
}

/********** pick the progenitor with the maximum score *****/
void mergertree::mergertree_select_maximum_score_progenitors(int nmatch)
{
  mycxxsort_parallel(desc, desc + nmatch, mergertree_compare_NewSubNr_PrevSubNr, Communicator);

  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  int *tab_CurrNsubhalos = (int *)Mem.mymalloc("tab_CurrNsubhalos", sizeof(int) * NTask);
  MPI_Allgather(&CurrNsubhalos, 1, MPI_INT, tab_CurrNsubhalos, 1, MPI_INT, Communicator);

  long long cumul_currnsubhalos = 0;
  for(int i = 0; i < ThisTask; i++)
    cumul_currnsubhalos += tab_CurrNsubhalos[i];

  desc_partdata *send_data = NULL;
  desc_partdata *recv_data = NULL;
  int nexport = 0, nimport = 0;
  for(int mode = 0; mode < 2; mode++)  // go through this twice to simplify bookkeeping
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int task                 = 0;
      unsigned long long first = 0;
      for(int i = 0; i < nmatch; i++)
        {
          while(task < NTask - 1 && desc[i].CurrSubhaloNr.get() >= first + tab_CurrNsubhalos[task])
            {
              first += tab_CurrNsubhalos[task];
              task++;
            }

          if(mode == 0)
            Send_count[task]++;
          else
            {
              int off = Send_offset[task] + Send_count[task]++;

              send_data[off] = desc[i];
            }
        }

      if(mode == 0)  // prepare offset tables
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

          Recv_offset[0] = 0;
          Send_offset[0] = 0;

          for(int j = 0; j < NTask; j++)
            {
              nexport += Send_count[j];
              nimport += Recv_count[j];

              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          send_data = (desc_partdata *)Mem.mymalloc("send_data", nexport * sizeof(desc_partdata));
          recv_data = (desc_partdata *)Mem.mymalloc("recv_data", nimport * sizeof(desc_partdata));
        }
    }

  /* exchange data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            myMPI_Sendrecv(&send_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(desc_partdata), MPI_BYTE, recvTask,
                         TAG_DENS_A, &recv_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(desc_partdata), MPI_BYTE,
                         recvTask, TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);
        }
    }

  for(int i = 0; i < CurrNsubhalos; i++)
    {
      Progenitors[i].SubhaloNr     = cumul_currnsubhalos + i;
      Progenitors[i].ProgSubhaloNr = HALONR_MAX;
      Progenitors[i].MaxScore      = 0;
    }

  for(int i = 0; i < nimport; i++)
    {
      int index = recv_data[i].CurrSubhaloNr.get() - cumul_currnsubhalos;

      if(index < 0 || index >= CurrNsubhalos)
        Terminate("index=%d  CurrNsubhalos=%d", index, CurrNsubhalos);

      if(recv_data[i].ProgScore > Progenitors[index].MaxScore)
        {
          Progenitors[index].MaxScore      = recv_data[i].ProgScore;
          Progenitors[index].ProgSubhaloNr = recv_data[i].PrevSubhaloNr.get();
        }
    }

  Mem.myfree(recv_data);
  Mem.myfree(send_data);

  Mem.myfree(tab_CurrNsubhalos);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

/********** determine the descendant with the maximum score *****/
void mergertree::mergertree_select_maximum_score_descendants(int nmatch)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* get lists of the subhalo count on each processors, and the cumulative number stored before */
  int *tab_PrevNsubhalos = (int *)Mem.mymalloc("tab_PrevNsubhalos", sizeof(int) * NTask);
  MPI_Allgather(&PrevNsubhalos, 1, MPI_INT, tab_PrevNsubhalos, 1, MPI_INT, Communicator);

  long long cumul_prevnsubhalos = 0;
  for(int i = 0; i < ThisTask; i++)
    cumul_prevnsubhalos += tab_PrevNsubhalos[i];

  desc_partdata *send_data = NULL;
  desc_partdata *recv_data = NULL;
  int nexport = 0, nimport = 0;
  for(int mode = 0; mode < 2; mode++)  // go through this twice to simplify bookkeeping
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int task                 = 0;
      unsigned long long first = 0;
      for(int i = 0; i < nmatch; i++)
        {
          if(PrevTotNsubhalos < 1)
            Terminate("PrevTotNsubhalos = %lld", PrevTotNsubhalos);

          while(task < NTask - 1 && desc[i].PrevSubhaloNr.get() >= first + tab_PrevNsubhalos[task])
            {
              first += tab_PrevNsubhalos[task];
              task++;
            }

          if(mode == 0)
            Send_count[task]++;
          else
            {
              int off = Send_offset[task] + Send_count[task]++;

              send_data[off] = desc[i];
            }
        }

      if(mode == 0)  // prepare offset tables
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

          Recv_offset[0] = 0;
          Send_offset[0] = 0;

          for(int j = 0; j < NTask; j++)
            {
              nexport += Send_count[j];
              nimport += Recv_count[j];

              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          send_data = (desc_partdata *)Mem.mymalloc("send_data", nexport * sizeof(desc_partdata));
          recv_data = (desc_partdata *)Mem.mymalloc("recv_data", nimport * sizeof(desc_partdata));
        }
    }

  /* exchange data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            myMPI_Sendrecv(&send_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(desc_partdata), MPI_BYTE, recvTask,
                         TAG_DENS_A, &recv_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(desc_partdata), MPI_BYTE,
                         recvTask, TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);
        }
    }

  for(int i = 0; i < PrevNsubhalos; i++)
    {
      Descendants[i].PrevSubhaloNr = cumul_prevnsubhalos + i;
      Descendants[i].DescSubhaloNr = HALONR_MAX;
      Descendants[i].MaxScore      = 0;
    }

  for(int i = 0; i < nimport; i++)
    {
      int index = recv_data[i].PrevSubhaloNr.get() - cumul_prevnsubhalos;

      if(index < 0 || index >= PrevNsubhalos)
        Terminate(
            "index=%d i=%d  nimport=%d  PrevNsubhalos=%d  recv_data[i].PrevSubhaloNr=%lld  PrevTotNsubhalos=%lld "
            "cumul_prevnsubhalos=%lld",
            index, i, nimport, PrevNsubhalos, (long long)recv_data[i].PrevSubhaloNr.get(), PrevTotNsubhalos, cumul_prevnsubhalos);

      if(recv_data[i].DescScore > Descendants[index].MaxScore)
        {
          Descendants[index].MaxScore      = recv_data[i].DescScore;
          Descendants[index].DescSubhaloNr = recv_data[i].CurrSubhaloNr.get();
        }
    }

  Mem.myfree(recv_data);
  Mem.myfree(send_data);

  Mem.myfree(tab_PrevNsubhalos);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

void mergertree::mergertree_set_first_descendant_with_same_progenitor(void)
{
  /* sort by progenitor to bring equal ones next to each other */
  mycxxsort_parallel(Progenitors, Progenitors + CurrNsubhalos, mergertree_compare_ProgSubhaloNr, Communicator);

  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  int *tab_PrevNsubhalos = (int *)Mem.mymalloc("tab_PrevNsubhalos", sizeof(int) * NTask);
  MPI_Allgather(&PrevNsubhalos, 1, MPI_INT, tab_PrevNsubhalos, 1, MPI_INT, Communicator);

  long long cumul_prevnsubhalos = 0;
  for(int i = 0; i < ThisTask; i++)
    cumul_prevnsubhalos += tab_PrevNsubhalos[i];

  struct pair_data
  {
    long long subhalonr;
    long long firstdescnr;
  };

  pair_data *send_data = NULL;
  pair_data *recv_data = NULL;
  int nexport = 0, nimport = 0;

  for(int mode = 0; mode < 2; mode++)  // go through this twice to simplify bookkeeping
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int task        = 0;
      long long first = 0;

      for(int i = 0; i < CurrNsubhalos; i++)
        {
          if(Progenitors[i].FirstDescFlag && Progenitors[i].ProgSubhaloNr != HALONR_MAX)
            {
              while(task < NTask - 1 && Progenitors[i].ProgSubhaloNr >= first + tab_PrevNsubhalos[task])
                {
                  first += tab_PrevNsubhalos[task];
                  task++;
                }

              if(mode == 0)
                Send_count[task]++;
              else
                {
                  int off = Send_offset[task] + Send_count[task]++;

                  send_data[off].subhalonr   = Progenitors[i].ProgSubhaloNr;
                  send_data[off].firstdescnr = Progenitors[i].SubhaloNr;
                }
            }
        }

      if(mode == 0)  // prepare offset tables
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

          Recv_offset[0] = 0;
          Send_offset[0] = 0;

          for(int j = 0; j < NTask; j++)
            {
              nexport += Send_count[j];
              nimport += Recv_count[j];

              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          send_data = (pair_data *)Mem.mymalloc("pair_data", nexport * sizeof(pair_data));
          recv_data = (pair_data *)Mem.mymalloc("pair_data", nimport * sizeof(pair_data));
        }
    }

  /* exchange data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            myMPI_Sendrecv(&send_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(pair_data), MPI_BYTE, recvTask, TAG_DENS_A,
                         &recv_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(pair_data), MPI_BYTE, recvTask, TAG_DENS_A,
                         Communicator, MPI_STATUS_IGNORE);
        }
    }

  for(int i = 0; i < PrevNsubhalos; i++)
    Descendants[i].FirstDescSubhaloNr = HALONR_MAX;

  for(int i = 0; i < nimport; i++)
    {
      int off = recv_data[i].subhalonr - cumul_prevnsubhalos;

      if(off < 0 || off >= PrevNsubhalos)
        Terminate("off = %d  PrevNsubhalos = %d", off, PrevNsubhalos);

      Descendants[off].FirstDescSubhaloNr = recv_data[i].firstdescnr;
    }

  Mem.myfree(recv_data);
  Mem.myfree(send_data);

  Mem.myfree(tab_PrevNsubhalos);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  /* bring back into original order */
  mycxxsort_parallel(Progenitors, Progenitors + CurrNsubhalos, mergertree_compare_SubhaloNr, Communicator);
}

#endif
