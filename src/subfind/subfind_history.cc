/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_history.cc
 *
 *  \brief this implements the SUBFIND_HBT algorithm for substructure finding
 */

#include "gadgetconfig.h"

#ifdef SUBFIND
#if defined(MERGERTREE) && defined(SUBFIND_HBT)

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
void fof<partset>::subfind_hbt_single_group(domain<partset> *SubDomain, domain<partset> *SingleDomain, domain_options mode, int gr)
{
  /* get total group length */
  long long totgrouplen;
  sumup_large_ints(1, &NumPartGroup, &totgrouplen, SubComm);

  /******************* determine subhalo candidates based on previous subhalo catalogue ***************/

  hbt_pcand_t *pcand = (hbt_pcand_t *)Mem.mymalloc_movable(&pcand, "pcand", (NumPartGroup + 1) * sizeof(hbt_pcand_t));

  for(int k = 0; k < NumPartGroup; k++)
    {
      /* provisionally assign a new subhalo number based on the previous group catalogue - this will be modified in the following  */
      Tp->PS[IndexList[k]].SubhaloNr = Tp->P[IndexList[k]].PrevSubhaloNr;

      pcand[k].SubhaloNr         = Tp->PS[IndexList[k]].SubhaloNr;
      pcand[k].PrevSizeOfSubhalo = Tp->P[IndexList[k]].PrevSizeOfSubhalo;
      pcand[k].index = -1;  // just put here to signal that this is invalid at this stage (note that we'll do a parallel sort)
    }

  /* sort according to subhalonr  */
  mycxxsort_parallel(pcand, pcand + NumPartGroup, subfind_hbt_compare_pcand_subhalonr, SubComm);

  int *NumPartGroup_list = (int *)Mem.mymalloc_movable(&NumPartGroup_list, "NumPartGroup_list", SubNTask * sizeof(int));
  MPI_Allgather(&NumPartGroup, sizeof(int), MPI_BYTE, NumPartGroup_list, sizeof(int), MPI_BYTE, SubComm);

  /* get the last element of each task */
  hbt_pcand_t *elem_last = (hbt_pcand_t *)Mem.mymalloc_movable(&elem_last, "elem_last", SubNTask * sizeof(hbt_pcand_t));

  /* note: the 0th element is guaranteed to be allocated even on ranks with zero NumPartGroup */
  MPI_Allgather(&pcand[NumPartGroup > 0 ? NumPartGroup - 1 : 0], sizeof(hbt_pcand_t), MPI_BYTE, elem_last, sizeof(hbt_pcand_t),
                MPI_BYTE, SubComm);

  /* if a new section begins on the current processor, we register it on this processor as a candidate */
  /* initialize the number of current candidates */
  bool element_before_present = false;
  hbt_pcand_t element_before{};

  for(int task = SubThisTask - 1; task >= 0; task--)
    {
      if(NumPartGroup_list[task] > 0)
        {
          element_before_present = true;
          element_before         = elem_last[task];
          break;
        }
    }

  int marked = 0;
  count_cand = 0;

  for(int i = 0; i < NumPartGroup; i++)
    {
      if(i == 0 && !element_before_present)
        count_cand++;
      else
        {
          MyHaloNrType prevnr;

          if(i == 0 && element_before_present)
            prevnr = element_before.SubhaloNr;
          else
            prevnr = pcand[i - 1].SubhaloNr;

          if(pcand[i].SubhaloNr != prevnr)
            count_cand++;
        }
    }

  /* allocate a list to store local subhalo candidates for this group */
  hbt_subcand_t *loc_candidates =
      (hbt_subcand_t *)Mem.mymalloc_movable(&loc_candidates, "loc_candidates", count_cand * sizeof(hbt_subcand_t));

  count_cand = 0;

  for(int i = 0; i < NumPartGroup; i++)
    {
      if(i == 0 && !element_before_present)
        {
          loc_candidates[count_cand].SubhaloNr = pcand[i].SubhaloNr;
          count_cand++;
        }
      else
        {
          MyHaloNrType prevnr;

          if(i == 0 && element_before_present)
            prevnr = element_before.SubhaloNr;
          else
            prevnr = pcand[i - 1].SubhaloNr;

          if(pcand[i].SubhaloNr != prevnr)
            {
              loc_candidates[count_cand].SubhaloNr = pcand[i].SubhaloNr;
              count_cand++;
            }
        }
    }

  /* establish total number of candidates */
  long long totcand;
  sumup_large_ints(1, &count_cand, &totcand, SubComm);

  int nsubhalos_old = Nsubhalos;

  if(Nsubhalos + totcand + 1 > MaxNsubhalos)
    {
      // warn("Nsubhalos=%d  + totcand=%lld >= MaxNsubhalos=%d", Nsubhalos, totcand, MaxNsubhalos);

      MaxNsubhalos = 1.25 * (Nsubhalos + totcand + 1);

      Subhalo = (subhalo_properties *)Mem.myrealloc_movable(Subhalo, MaxNsubhalos * sizeof(subhalo_properties));
    }

  /* assemble a list of the candidates on all tasks */
  hbt_subcand_t *all_candidates =
      (hbt_subcand_t *)Mem.mymalloc_movable(&all_candidates, "all_candidates", (totcand + 1) * sizeof(hbt_subcand_t));

  int *countlist = (int *)Mem.mymalloc_movable(&countlist, "countlist", SubNTask * sizeof(int));
  int *offset    = (int *)Mem.mymalloc_movable(&offset, "offset", SubNTask * sizeof(int));

  int count = count_cand * sizeof(hbt_subcand_t); /* length in bytes */
  MPI_Allgather(&count, 1, MPI_INT, countlist, 1, MPI_INT, SubComm);

  offset[0] = 0;
  for(int i = 1; i < SubNTask; i++)
    offset[i] = offset[i - 1] + countlist[i - 1];

  myMPI_Allgatherv(loc_candidates, count, MPI_BYTE, all_candidates, countlist, offset, MPI_BYTE, SubComm);

  /* sort the candidates by subhalonr */
  mycxxsort(all_candidates, all_candidates + totcand, subfind_hbt_compare_subcand_subhalonr);

  /* now determine the size of the candidates */
  long long *size_list            = (long long *)Mem.mymalloc_clear("size_list", totcand * sizeof(long long));
  long long *summed_prevsize_list = (long long *)Mem.mymalloc_clear("summed_prevsize_list", totcand * sizeof(long long));

  int j = 0;
  for(int k = 0; k < NumPartGroup; k++)
    {
      if(j >= totcand)
        Terminate("can't be: k=%d   j=%d  NumPartGroup=%d totcand=%lld\n", k, j, NumPartGroup, totcand);

      while(j < totcand && pcand[k].SubhaloNr > all_candidates[j].SubhaloNr)
        j++;

      if(pcand[k].SubhaloNr != all_candidates[j].SubhaloNr)
        Terminate("can't be:  k=%d NumPartGroup=%d   pcand[k].SubhaloNr=%lld    j=%d  all_candidates[j].SubhaloNr=%lld\n", k,
                  NumPartGroup, (long long)pcand[k].SubhaloNr.get(), j, (long long)all_candidates[j].SubhaloNr.get());

      size_list[j]++;
      summed_prevsize_list[j] += pcand[k].PrevSizeOfSubhalo.get();
    }

  MPI_Allreduce(MPI_IN_PLACE, size_list, totcand, MPI_LONG_LONG, MPI_SUM, SubComm);
  MPI_Allreduce(MPI_IN_PLACE, summed_prevsize_list, totcand, MPI_LONG_LONG, MPI_SUM, SubComm);

  for(int i = 0; i < totcand; i++)
    {
      all_candidates[i].len           = size_list[i];
      all_candidates[i].summedprevlen = summed_prevsize_list[i];
    }

  Mem.myfree(summed_prevsize_list);
  Mem.myfree(size_list);

  /* do a sanity test */
  long long lensum = 0;
  for(int i = 0; i < totcand; i++)
    lensum += all_candidates[i].len;

  if(lensum != totgrouplen)
    Terminate("lensum=%lld != totgrouplen=%lld\n", lensum, totgrouplen);

  /*******************************************/

  /* find the group of previously unbound ones, and if this candidate exists, eliminate it */
  for(int i = 0; i < totcand; i++)
    if(all_candidates[i].SubhaloNr.get() == HALONR_MAX)
      {
        all_candidates[i] = all_candidates[totcand - 1];
        totcand--;
        break;
      }

  /* let's now eliminate small groups and flag the corresponding particles as unbound */
  {
    /* sort the candidates according to previous subhalonr */
    mycxxsort(all_candidates, all_candidates + totcand, subfind_hbt_compare_subcand_subhalonr);

    /* reestablish a locally sorted pcand */
    for(int k = 0; k < NumPartGroup; k++)
      {
        pcand[k].SubhaloNr = Tp->PS[IndexList[k]].SubhaloNr;
        pcand[k].index     = IndexList[k];
      }
    mycxxsort(pcand, pcand + NumPartGroup, subfind_hbt_compare_pcand_subhalonr);

    long long sum = 0;

    int p = 0;

    for(int i = 0; i < totcand; i++)
      if(all_candidates[i].len < All.DesLinkNgb)
        {
          while(p < NumPartGroup && pcand[p].SubhaloNr < all_candidates[i].SubhaloNr)
            p++;

          while(p < NumPartGroup && pcand[p].SubhaloNr == all_candidates[i].SubhaloNr)
            {
              if(Tp->PS[pcand[p].index].SubhaloNr != all_candidates[i].SubhaloNr)
                Terminate(
                    "we have an issue! p=%d  NumPartGroup=%d  pcand[p].index=%d  pcand[p].SubhaloNr=%lld  "
                    "all_candidates[i].SubhaloNr=%lld  Tp->P[pcand[p].index].SubhaloNr=%lld",
                    p, NumPartGroup, pcand[p].index, (long long)pcand[p].SubhaloNr.get(), (long long)all_candidates[i].SubhaloNr.get(),
                    (long long)Tp->PS[pcand[p].index].SubhaloNr.get());

              pcand[p].SubhaloNr.set(HALONR_MAX);
              Tp->PS[pcand[p].index].SubhaloNr.set(HALONR_MAX);
              p++;
              sum++;
            }
        }

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_LONG_LONG, MPI_SUM, SubComm);

    long long sum2 = 0;

    for(int i = 0; i < totcand; i++)
      if(all_candidates[i].len < All.DesLinkNgb)
        {
          sum2 += all_candidates[i].len;
          all_candidates[i] = all_candidates[totcand - 1];
          totcand--;
          i--;
        }

    if(sum != sum2)
      Terminate("consistency check failed sum = %lld  sum2 = %lld\n", sum, sum2);

    /* sort according to subhalonr  */
    mycxxsort_parallel(pcand, pcand + NumPartGroup, subfind_hbt_compare_pcand_subhalonr, SubComm);
  }

  /* if a largest one exists, eliminate it, because this is the one we will treat as background halo */
  if(totcand > 0)
    {
      /* sort the candidates by size */
      mycxxsort(all_candidates, all_candidates + totcand, subfind_hbt_compare_subcand_len);

      /* sort the candidates by summed previous length, as this is arguably a more robust decision of which one should be the largest
       */
      mycxxsort(all_candidates, all_candidates + totcand, subfind_hbt_compare_subcand_summedprevlen);

      totcand--;
      for(int k = 0; k < NumPartGroup; k++)
        if(Tp->PS[IndexList[k]].SubhaloNr == all_candidates[totcand].SubhaloNr)
          Tp->PS[IndexList[k]].SubhaloNr.set(HALONR_MAX);
    }

  /* sort the candidates according to previous subhalonr */
  mycxxsort(all_candidates, all_candidates + totcand, subfind_hbt_compare_subcand_subhalonr);

  subfind_collective_printf("SUBFIND: root-task=%d: total number of subhalo coll_candidates=%lld\n", ThisTask, totcand);

  /*******************************************/
  /* Let's now see which candidates can be treated with serial CPUs, which is more efficient than doing them all collectively.
   * We identify them with those candidates that are sufficiently small, which should be most of them by number.
   */

  int n_small_cand = 0;
  int max_length   = 0;

  int task         = 0;
  int task_scatter = 0;
  for(int i = 0; i < totcand; i++)
    {
      if(all_candidates[i].len < 0.20 * Tp->TotNumPart / NTask)  // small enough
        {
          all_candidates[i].DoIt        = true;
          all_candidates[i].TargetTask  = task++;
          all_candidates[i].TargetIndex = n_small_cand;

          if(task >= SubNTask)
            task = 0;

          if(all_candidates[i].len > max_length)
            max_length = all_candidates[i].len;

          n_small_cand++;
        }
      else
        {
          all_candidates[i].DoIt        = false;
          all_candidates[i].TargetTask  = task_scatter++;
          all_candidates[i].TargetIndex = INT_MAX;
          if(task_scatter >= SubNTask)
            task_scatter = 0;
        }
    }

  subfind_collective_printf(
      "SUBFIND: root-task=%d: number of subhalo candidates small enough to be done with one cpu: %d. (Largest size %d)\n", ThisTask,
      n_small_cand, max_length);

  /* now get target information to particles */
  for(int k = 0; k < NumPartGroup; k++)
    {
      pcand[k].SubhaloNr = Tp->PS[IndexList[k]].SubhaloNr;
      pcand[k].index     = IndexList[k];
    }

  // note: local serial sort sufficient here
  mycxxsort(pcand, pcand + NumPartGroup, subfind_hbt_compare_pcand_subhalonr);

  if(SubNTask > 1)
    {
      /* we only need to redistribute the particles if we are processing groups collectively */
      /* Note: setting of TargetIndex make sure that the particles are grouped together on the target task */

      /* set default values for current particle distribution */
      for(int i = 0; i < Tp->NumPart; i++)
        {
          Tp->PS[i].u.s.origintask  = SubThisTask;
          Tp->PS[i].u.s.originindex = i;

          Tp->PS[i].TargetTask  = SubThisTask;
          Tp->PS[i].TargetIndex = INT_MAX;
        }

      int i = 0;
      for(int k = 0; k < totcand; k++)
        if(all_candidates[k].DoIt)
          {
            while(i < NumPartGroup && pcand[i].SubhaloNr < all_candidates[k].SubhaloNr)
              i++;

            while(i < NumPartGroup && pcand[i].SubhaloNr == all_candidates[k].SubhaloNr)
              {
                Tp->PS[pcand[i].index].TargetTask  = all_candidates[k].TargetTask;
                Tp->PS[pcand[i].index].TargetIndex = all_candidates[k].TargetIndex;
                i++;
              }
          }

      /* assemble the particles on individual processors (note: IndexList[] becomes temporarily meaningless)  */
      subfind_distribute_particles(SubComm);
    }

  /* now do the serial unbinding */
  /*----------------------------------------------------*/

  {
    unbind_list = (int *)Mem.mymalloc("unbind_list", Tp->NumPart * sizeof(int));

    int i = 0;  // particle index

    for(int k = 0; k < totcand; k++)
      if(all_candidates[k].DoIt)
        if(all_candidates[k].TargetTask == SubThisTask)
          {
            int len = 0;

            if(SubNTask > 1)
              {
                while(i < Tp->NumPart && Tp->PS[i].SubhaloNr < all_candidates[k].SubhaloNr)
                  i++;

                while(i < Tp->NumPart && Tp->PS[i].SubhaloNr == all_candidates[k].SubhaloNr && Tp->PS[i].GroupNr.get() == GroupNr)
                  {
                    unbind_list[len] = i;
                    len++;
                    i++;
                  }
              }
            else
              {
                while(i < NumPartGroup && Tp->PS[pcand[i].index].SubhaloNr < all_candidates[k].SubhaloNr)
                  i++;

                while(i < NumPartGroup && Tp->PS[pcand[i].index].SubhaloNr == all_candidates[k].SubhaloNr &&
                      Tp->PS[pcand[i].index].GroupNr.get() == GroupNr)
                  {
                    unbind_list[len] = pcand[i].index;
                    len++;
                    i++;
                  }
              }

            if(len != all_candidates[k].len)
              Terminate("this is unexpected: k=%d   len=%lld != all_candidates[k].len=%lld) \n", k, (long long)len,
                        (long long)all_candidates[k].len);

            /* default is that all particles end up unbound */
            for(int n = 0; n < len; n++)
              Tp->PS[unbind_list[n]].SubhaloNr.set(HALONR_MAX);

            /* call the serial unbind function */
            len = subfind_unbind(SingleDomain, SingleDomain->Communicator, unbind_list, len);

            if(len >= All.DesLinkNgb)
              {
                /* set as provisional group number the previous group number */
                for(int n = 0; n < len; n++)
                  {
                    Tp->PS[unbind_list[n]].SubhaloNr = all_candidates[k].SubhaloNr;
                    Tp->PS[unbind_list[n]].SizeOfSubhalo.set(len);
                  }

                if(Nsubhalos >= MaxNsubhalos)
                  Terminate("no storage: Nsubhalos=%d  MaxNsubhalos=%d  nsubhalos_old=%d totcand=%lld\n", Nsubhalos, MaxNsubhalos,
                            nsubhalos_old, totcand);

                /* ok, we found a substructure */
                marked += subfind_determine_sub_halo_properties(unbind_list, len, &Subhalo[Nsubhalos], SingleDomain->Communicator);

                Subhalo[Nsubhalos].GroupNr       = GroupNr;
                Subhalo[Nsubhalos].SubParentRank = 0;
                Subhalo[Nsubhalos].SubhaloNr     = all_candidates[k].SubhaloNr.get();

                Nsubhalos++;
              }
          }

    Mem.myfree(unbind_list);
  }

  if(SubNTask > 1)
    {
      /* bring them back to their original processor */
      for(int i = 0; i < Tp->NumPart; i++)
        {
          Tp->PS[i].TargetTask  = Tp->PS[i].u.s.origintask;
          Tp->PS[i].TargetIndex = Tp->PS[i].u.s.originindex;
        }

      subfind_distribute_particles(SubComm);
    }

  double t0 = Logs.second();

  /**************************************************/
  /**************************************************/
  /*******  now do remaining ones collectively  *****/

  /* first, add a fiducial candidate which will be our background halo, swallowing all unbound particles */

  all_candidates[totcand].DoIt = false; /* marks collective ones */
  all_candidates[totcand].SubhaloNr.set(HALONR_MAX);
  totcand++;

  for(int k = 0; k < totcand; k++)
    if(all_candidates[k].DoIt == false)
      {
        domain<partset> SubUnbindDomain{SubComm, Tp};

        int *unbind_list;

        if(mode == COLL_SUBFIND)
          {
            for(int i = 0; i < Tp->NumPart; i++)
              {
                Tp->PS[i].u.s.origintask  = SubThisTask;
                Tp->PS[i].u.s.originindex = i;
                Tp->PS[i].DomainFlag      = 0;
              }

            /* mark the one to be unbound in PS[] */
            for(int i = 0; i < NumPartGroup; i++)
              if(Tp->PS[IndexList[i]].SubhaloNr == all_candidates[k].SubhaloNr)
                Tp->PS[IndexList[i]].DomainFlag = 1;

            SubUnbindDomain.domain_decomposition(mode);
            subfind_distribute_particles(SubComm);

            LocalLen = 0;
            for(int i = 0; i < Tp->NumPart; i++)
              if(Tp->PS[i].DomainFlag)
                LocalLen++;

            unbind_list = (int *)Mem.mymalloc_movable(&unbind_list, "unbind_list", LocalLen * sizeof(int));

            /* refill unbind_list */
            LocalLen = 0;
            for(int i = 0; i < Tp->NumPart; i++)
              if(Tp->PS[i].DomainFlag)
                unbind_list[LocalLen++] = i;
          }

        else
          {
            unbind_list = (int *)Mem.mymalloc_movable(&unbind_list, "unbind_list", NumPartGroup * sizeof(int));

            LocalLen = 0;
            for(int i = 0; i < NumPartGroup; i++)
              if(Tp->PS[IndexList[i]].SubhaloNr == all_candidates[k].SubhaloNr)
                unbind_list[LocalLen++] = IndexList[i];
          }

        /* default is that all particles end up unbound */
        for(int n = 0; n < LocalLen; n++)
          Tp->PS[unbind_list[n]].SubhaloNr.set(HALONR_MAX);

        if(mode == COLL_SUBFIND)
          LocalLen = subfind_unbind(&SubUnbindDomain, SubComm, unbind_list, LocalLen);
        else
          LocalLen = subfind_unbind(SingleDomain, SingleDomain->Communicator, unbind_list, LocalLen);

        int FullLen;
        MPI_Allreduce(&LocalLen, &FullLen, 1, MPI_INT, MPI_SUM, SubComm);

        if(FullLen >= All.DesLinkNgb)
          {
            if(all_candidates[k].SubhaloNr.get() == HALONR_MAX)
              all_candidates[k].SubhaloNr.set(HALONR_MAX - 1);

            /* set as provisional group number the previous group number */
            for(int n = 0; n < LocalLen; n++)
              {
                Tp->PS[unbind_list[n]].SubhaloNr = all_candidates[k].SubhaloNr;
#if defined(MERGERTREE)
                Tp->PS[unbind_list[n]].SizeOfSubhalo.set(FullLen);
#endif
              }

            if(Nsubhalos >= MaxNsubhalos)
              Terminate("no storage: Nsubhalos=%d  MaxNsubhalos=%d  nsubhalos_old=%d totcand=%lld\n", Nsubhalos, MaxNsubhalos,
                        nsubhalos_old, totcand);

            marked += subfind_determine_sub_halo_properties(unbind_list, LocalLen, &Subhalo[Nsubhalos], SubComm);

            if(SubThisTask == 0)
              {
                Subhalo[Nsubhalos].GroupNr       = GroupNr;
                Subhalo[Nsubhalos].SubParentRank = 0;
                Subhalo[Nsubhalos].SubhaloNr     = all_candidates[k].SubhaloNr.get();

                Nsubhalos++;
              }
          }

        Mem.myfree(unbind_list);

        if(mode == COLL_SUBFIND)
          {
            SubUnbindDomain.domain_free();

            for(int i = 0; i < Tp->NumPart; i++)
              {
                Tp->PS[i].TargetTask  = Tp->PS[i].u.s.origintask;
                Tp->PS[i].TargetIndex = Tp->PS[i].u.s.originindex;
              }

            double t0 = Logs.second();
            subfind_distribute_particles(SubComm); /* bring them back to their original processor */
            double t1 = Logs.second();

            subfind_collective_printf("SUBFIND: root-task=%d: bringing the independent ones back took %g sec\n", ThisTask,
                                      Logs.timediff(t0, t1));

            /* Since we reestablished the original order, we can use IndexList[] again */
          }
      }

  double t1 = Logs.second();
  subfind_collective_printf("SUBFIND: root-task=%d: the collective unbinding of remaining halos took %g sec\n", ThisTask,
                            Logs.timediff(t0, t1));

  /* get the total substructure count */
  int countloc = Nsubhalos - nsubhalos_old;
  int countall;

  MPI_Allreduce(&countloc, &countall, 1, MPI_INT, MPI_SUM, SubComm);

  MPI_Allreduce(MPI_IN_PLACE, &marked, 1, MPI_INT, MPI_SUM, SubComm);

  /* Now let's save some properties of the substructures */
  if(SubThisTask == 0)
    {
      Group[gr].Nsubs = countall;

#if defined(SUBFIND_ORPHAN_TREATMENT)
      Group[gr].LenPrevMostBnd += marked;
#endif
    }

  /* also need to set the field SubRankInGr in all found Subhalos */

  hbt_subhalo_t *subhalo_list = (hbt_subhalo_t *)Mem.mymalloc("subhalo_list", countloc * sizeof(hbt_subhalo_t));

  for(int n = 0; n < countloc; n++)
    {
      subhalo_list[n].Len       = Subhalo[n + nsubhalos_old].Len;
      subhalo_list[n].ThisTask  = SubThisTask;
      subhalo_list[n].ThisIndex = n;
      subhalo_list[n].SubhaloNr = Subhalo[n + nsubhalos_old].SubhaloNr;
    }

  mycxxsort_parallel(subhalo_list, subhalo_list + countloc, subfind_hbt_compare_subhalolist_len, SubComm);

  int *countloc_list = (int *)Mem.mymalloc("countloc_list", SubNTask * sizeof(int));
  MPI_Allgather(&countloc, 1, MPI_INT, countloc_list, 1, MPI_INT, SubComm);
  int npreviousranks = 0;
  for(int i = 0; i < SubThisTask; i++)
    npreviousranks += countloc_list[i];

  for(int n = 0; n < countloc; n++)
    subhalo_list[n].SubRankInGr = n + npreviousranks;

  mycxxsort_parallel(subhalo_list, subhalo_list + countloc, subfind_hbt_compare_subhalolist_thistask_thisindex, SubComm);

  /* rank and index of main subhalo */
  int taskmain  = 0;
  int indexmain = 0;

  for(int n = 0; n < countloc; n++)
    {
      Subhalo[n + nsubhalos_old].SubRankInGr = subhalo_list[n].SubRankInGr;

      if(subhalo_list[n].SubRankInGr == 0)
        {
          /* here we have the main subhalo */
          taskmain  = SubThisTask;
          indexmain = n + nsubhalos_old;
        }
    }

  /* now we need to fill in SubRankInGr to T->PS[] so that the particles can be sorted in subhalo order
   * we use subhalo_list[] as translation table
   */

  for(int k = 0; k < NumPartGroup; k++)
    {
      pcand[k].SubhaloNr = Tp->PS[IndexList[k]].SubhaloNr;
      pcand[k].index     = IndexList[k];
    }

  /* sort locally  */
  mycxxsort(pcand, pcand + NumPartGroup, subfind_hbt_compare_pcand_subhalonr);

  int sizelocsubhalolist = countloc * sizeof(hbt_subhalo_t); /* length in bytes */
  MPI_Allgather(&sizelocsubhalolist, 1, MPI_INT, countlist, 1, MPI_INT, SubComm);

  offset[0] = 0;
  for(int i = 1; i < SubNTask; i++)
    offset[i] = offset[i - 1] + countlist[i - 1];

  hbt_subhalo_t *all_subhalo_list = (hbt_subhalo_t *)Mem.mymalloc("all_subhalo_list", countall * sizeof(hbt_subhalo_t));

  myMPI_Allgatherv(subhalo_list, sizelocsubhalolist, MPI_BYTE, all_subhalo_list, countlist, offset, MPI_BYTE, SubComm);

  /* sort locally  */
  mycxxsort(all_subhalo_list, all_subhalo_list + countall, subfind_hbt_compare_subhalolist_prevsubhalonr);

  int n = 0;

  for(int k = 0; k < NumPartGroup; k++)
    {
      if(pcand[k].SubhaloNr.get() == HALONR_MAX)
        Tp->PS[pcand[k].index].SubRankInGr = INT_MAX;
      else
        {
          while(n < countall && all_subhalo_list[n].SubhaloNr < (long long)pcand[k].SubhaloNr.get())
            n++;

          if(n >= countall)
            Terminate("unexpected: n=%d countall=%d", n, countall);

          if(all_subhalo_list[n].SubhaloNr != (long long)pcand[k].SubhaloNr.get())
            Terminate("also unexpected: k=%d NumPartGroup=%d  all_subhalo_list[n].SubhaloNr=%lld != pcand[k].SubhaloNr=%lld\n", k,
                      NumPartGroup, (long long)all_subhalo_list[n].SubhaloNr, (long long)pcand[k].SubhaloNr.get());

          Tp->PS[pcand[k].index].SubRankInGr = all_subhalo_list[n].SubRankInGr;
        }
    }

  if(countall > 0)
    {
      MPI_Allreduce(MPI_IN_PLACE, &taskmain, 1, MPI_INT, MPI_MAX, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, &indexmain, 1, MPI_INT, MPI_MAX, SubComm);

      subhalo_properties MainSubhalo;

      if(taskmain == SubThisTask)
        {
          MainSubhalo = Subhalo[indexmain];

          if(taskmain != 0)
            MPI_Send(&MainSubhalo, sizeof(subhalo_properties), MPI_BYTE, 0, TAG_N, SubComm);
        }

      if(SubThisTask == 0)
        {
          if(taskmain != 0)
            MPI_Recv(&MainSubhalo, sizeof(subhalo_properties), MPI_BYTE, taskmain, TAG_N, SubComm, MPI_STATUS_IGNORE);

          for(int j = 0; j < 3; j++)
            {
              Group[gr].Pos[j]    = MainSubhalo.Pos[j];
              Group[gr].IntPos[j] = MainSubhalo.IntPos[j];
            }
        }
    }

  Mem.myfree(all_subhalo_list);
  Mem.myfree(countloc_list);
  Mem.myfree(subhalo_list);
  Mem.myfree(offset);
  Mem.myfree(countlist);
  Mem.myfree(all_candidates);
  Mem.myfree(loc_candidates);
  Mem.myfree(elem_last);
  Mem.myfree(NumPartGroup_list);
  Mem.myfree(pcand);

  subfind_collective_printf("SUBFIND: root-task=%d: found %d bound substructures in FoF group of length %lld\n", ThisTask, countall,
                            totgrouplen);
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
#endif
