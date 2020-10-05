/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_excursionset.cc
 *
 *  \brief main routines for processing a halo with the classic Subfind algorithm
 */

#include "gadgetconfig.h"

#ifdef SUBFIND
#ifndef SUBFIND_HBT

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

#define HIGHBIT (1 << 30)

template <typename partset>
void fof<partset>::subfind_process_single_group(domain<partset> *SubDomain, domain<partset> *SingleDomain, domain_options mode, int gr)
{
  /* set up an inverse look-up of the position of local particles in the index list, for later use
   */
  for(int i = 0; i < NumPartGroup; i++)
    Tp->PS[IndexList[i]].InvIndex = i;

  /* allocated storage for an auxiliary array needed for sorting */
  sd = (sort_density_data *)Mem.mymalloc_movable(&sd, "sd", NumPartGroup * sizeof(sort_density_data));

  /* construct a tree for all particles in the halo */

  FoFGravTree.treeallocate(Tp->NumPart, Tp, SubDomain);
  FoFGravTree.treebuild(NumPartGroup, IndexList);

  /* determine the radius that encloses a certain number of link particles */
  subfind_find_linkngb(SubDomain, NumPartGroup, IndexList);

  /* now determine the indices of the nearest two denser neighbours within this link region */

  /* first, allocate some auxialiary arrays for storing this info */

  Tp->R2Loc = (typename partset::nearest_r2_data *)Mem.mymalloc_movable(&Tp->R2Loc, "R2Loc",
                                                                        Tp->NumPart * sizeof(typename partset::nearest_r2_data));

  /* find the nearest two denser particles around each point in the group */
  subfind_find_nearesttwo(SubDomain, NumPartGroup, IndexList);

  /* create an array that we can conveniently sort according to density in group subsets */
  for(int i = 0; i < NumPartGroup; i++)
    {
      if(IndexList[i] >= Tp->NumPart || IndexList[i] < 0)
        Terminate("Really?");

      sd[i].density    = Tp->PS[IndexList[i]].u.s.u.DM_Density;
      sd[i].ngbcount   = Tp->PS[IndexList[i]].nearest.count;
      sd[i].index      = {SubThisTask, i};
      sd[i].ngb_index1 = Tp->PS[IndexList[i]].nearest.index[0];
      sd[i].ngb_index2 = Tp->PS[IndexList[i]].nearest.index[1];
#ifdef MERGERTREE
      sd[i].PrevSizeOfSubhalo = Tp->P[IndexList[i]].PrevSizeOfSubhalo;
#else
      sd[i].PrevSizeOfSubhalo.set(0);
#endif
    }
  Mem.myfree(Tp->R2Loc);

  /* sort sd according to densities (in parallel if needed) */
  mycxxsort_parallel(sd, sd + NumPartGroup, subfind_compare_densities, SubComm);

  subfind_collective_printf("SUBFIND: root-task=%d: parallel sort of 'sd' done.\n", ThisTask);

  /* can now release the tree */
  FoFGravTree.treefree();

  /* allocate and initialize distributed link list for storing subhalo connectivity */
  SFHead    = (location *)Mem.mymalloc_movable(&SFHead, "SFHead", NumPartGroup * sizeof(location));
  SFNext    = (location *)Mem.mymalloc_movable(&SFNext, "SFNext", NumPartGroup * sizeof(location));
  SFTail    = (location *)Mem.mymalloc_movable(&SFTail, "SFTail", NumPartGroup * sizeof(location));
  SFLen     = (MyLenType *)Mem.mymalloc_movable(&SFLen, "SFLen", NumPartGroup * sizeof(MyLenType));
  SFPrevLen = (double *)Mem.mymalloc_movable(&SFPrevLen, "SFPrevLen", NumPartGroup * sizeof(double));

  for(int i = 0; i < NumPartGroup; i++)
    {
      SFHead[i]    = {-1, -1};
      SFNext[i]    = {-1, -1};
      SFTail[i]    = {-1, -1};
      SFLen[i]     = 0;
      SFPrevLen[i] = 0;
    }

  /* allocate a list to store subhalo candidates */
  max_coll_candidates = std::max<int>((NumPartGroup / 50), 200);
  coll_candidates =
      (coll_cand_dat *)Mem.mymalloc_movable(&coll_candidates, "coll_candidates", max_coll_candidates * sizeof(coll_cand_dat));

  /* initialize the number of current candidates */
  count_cand = 0;

  /* get total group length */
  long long totgrouplen;
  sumup_large_ints(1, &NumPartGroup, &totgrouplen, SubComm);

  /* determine subhalo candidates */
  subfind_col_find_coll_candidates(totgrouplen);

  /* establish total number of candidates */
  long long totcand;
  sumup_large_ints(1, &count_cand, &totcand, SubComm);

  subfind_collective_printf("SUBFIND: root-task=%d: total number of subhalo coll_candidates=%lld\n", ThisTask, totcand);

  for(int i = 0; i < NumPartGroup; i++)
    SFTail[i] = {-1, -1};

  /* default is to be not nested */
  for(int i = 0; i < count_cand; i++)
    coll_candidates[i].parent = 0;

  int count_leaves     = 0;
  long long nremaining = totcand;

  do
    {
      /* Let's see which coll_candidates can be unbound independent from each other.
       * We identify them with those candidates that have no embedded other candidate, which are most of them by number.
       */
      double t0                          = Logs.second();
      coll_cand_dat *tmp_coll_candidates = 0;
      if(SubThisTask == 0)
        tmp_coll_candidates = (coll_cand_dat *)Mem.mymalloc("tmp_coll_candidates", totcand * sizeof(coll_cand_dat));

      int *countlist = (int *)Mem.mymalloc("countlist", SubNTask * sizeof(int));
      int *offset    = (int *)Mem.mymalloc("offset", SubNTask * sizeof(int));

      int count = count_cand * sizeof(coll_cand_dat);
      MPI_Allgather(&count, 1, MPI_INT, countlist, 1, MPI_INT, SubComm);

      offset[0] = 0;
      for(int i = 1; i < SubNTask; i++)
        offset[i] = offset[i - 1] + countlist[i - 1];

      /* assemble a list of the candidates on subtask 0 */
      MPI_Gatherv(coll_candidates, countlist[SubThisTask], MPI_BYTE, tmp_coll_candidates, countlist, offset, MPI_BYTE, 0, SubComm);

      if(SubThisTask == 0)
        {
          for(int k = 0; k < totcand; k++)
            {
              tmp_coll_candidates[k].nsub  = k;
              tmp_coll_candidates[k].subnr = k;
            }

          mycxxsort(tmp_coll_candidates, tmp_coll_candidates + totcand, subfind_compare_coll_candidates_rank);
          for(int k = 0; k < totcand; k++)
            {
              if(tmp_coll_candidates[k].parent >= 0)
                {
                  tmp_coll_candidates[k].parent = 0;

                  for(int j = k + 1; j < totcand; j++)
                    {
                      if(tmp_coll_candidates[j].rank > tmp_coll_candidates[k].rank + tmp_coll_candidates[k].len)
                        break;

                      if(tmp_coll_candidates[j].parent < 0) /* ignore these */
                        continue;

                      if(tmp_coll_candidates[k].rank + tmp_coll_candidates[k].len >=
                         tmp_coll_candidates[j].rank + tmp_coll_candidates[j].len)
                        {
                          tmp_coll_candidates[k].parent++; /* we here count the number of subhalos that are enclosed */
                        }
                      else
                        {
                          Terminate("k=%d|%lld has rank=%d and len=%d.  j=%d has rank=%d and len=%d\n", k, totcand,
                                    (int)tmp_coll_candidates[k].rank, (int)tmp_coll_candidates[k].len, j,
                                    (int)tmp_coll_candidates[j].rank, (int)tmp_coll_candidates[j].len);
                        }
                    }
                }
            }

          mycxxsort(tmp_coll_candidates, tmp_coll_candidates + totcand, subfind_compare_coll_candidates_subnr);
        }

      /* send the stuff back */
      MPI_Scatterv(tmp_coll_candidates, countlist, offset, MPI_BYTE, coll_candidates, countlist[SubThisTask], MPI_BYTE, 0, SubComm);

      Mem.myfree(offset);
      Mem.myfree(countlist);

      if(SubThisTask == 0)
        Mem.myfree(tmp_coll_candidates);

      count_leaves   = 0;
      int max_length = 0;
      for(int i = 0; i < count_cand; i++)
        if(coll_candidates[i].parent == 0) /* if it's not a nested one, candidate is eligible for the independent/parallel unbinding */
          {
            /* if it seems large (heuristic criterion), let's rather do it collectively */
            if(coll_candidates[i].len > 0.20 * Tp->TotNumPart / NTask)
              coll_candidates[i].parent++; /* this trick will ensure that it is not considered in this round */
            else
              {
                if(coll_candidates[i].len > max_length)
                  max_length = coll_candidates[i].len;

                count_leaves++;
              }
          }

      /* get total count of these eligible subhalos, and their maximum length */
      MPI_Allreduce(MPI_IN_PLACE, &count_leaves, 1, MPI_INT, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, &max_length, 1, MPI_INT, MPI_MAX, SubComm);

      double t1 = Logs.second();

      subfind_collective_printf(
          "SUBFIND: root-task=%d: number of subhalo coll_candidates that can be done independently=%d. (Largest size %d, finding took "
          "%g sec)\n",
          ThisTask, count_leaves, max_length, Logs.timediff(t0, t1));

      /* if there are none left, we break and do the reset collectively */
      if(count_leaves <= 0)
        {
          subfind_collective_printf("SUBFIND: root-task=%d: too few, let's do the rest of %d collectively\n", ThisTask, nremaining);
          break;
        }

      /* seems large, let's rather do it collectively */
      if(max_length > 0.5 * Tp->TotNumPart / NTask)
        {
          subfind_collective_printf("SUBFIND: root-task=%d: too big coll_candidates, I do the rest collectively\n", ThisTask);
          break;
        }

      nremaining -= count_leaves;

      /* set default values */
      for(int i = 0; i < Tp->NumPart; i++)
        {
          Tp->PS[i].u.s.origintask  = SubThisTask;
          Tp->PS[i].u.s.originindex = i;

          Tp->PS[i].TargetTask = SubThisTask;
          Tp->PS[i].submark    = HIGHBIT;
        }

      for(int i = 0; i < NumPartGroup; i++)
        {
          if(SFTail[i].index >= 0) /* this means this particle is already bound to a substructure */
            Tp->PS[IndexList[i]].u.s.origintask |= HIGHBIT;
        }

      /* we now mark the particles that are in subhalo candidates that can be processed independently in parallel */
      int nsubs = 0;
      for(int master = 0; master < SubNTask; master++)
        {
          int ncand = count_cand;
          MPI_Bcast(&ncand, sizeof(ncand), MPI_BYTE, master, SubComm);

          for(int k = 0; k < ncand; k++)
            {
              MyLenType len;
              int parent;
              if(SubThisTask == master)
                {
                  len    = coll_candidates[k].len;
                  parent = coll_candidates[k].parent; /* this is here actually the daughter count */
                }

              MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, SubComm);
              MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, SubComm);
              MPI_Barrier(SubComm);

              if(parent == 0)
                {
                  if(SubThisTask != master)
                    subfind_poll_for_requests();
                  else
                    {
                      location p = coll_candidates[k].head;
                      for(MyLenType i = 0; i < coll_candidates[k].len; i++)
                        {
                          subfind_distlinklist_mark_particle(p, master, nsubs);

                          if(p.index < 0)
                            Terminate("Bummer\n");

                          p = subfind_distlinklist_get_next(p);
                        }

                      /* now tell the others to stop polling */
                      for(int i = 0; i < SubNTask; i++)
                        if(i != SubThisTask)
                          MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
                    }

                  MPI_Barrier(SubComm);
                }

              nsubs++;
            }
        }

      if(mode == COLL_SUBFIND)
        {
          /* this will make sure that the particles are grouped by submark on the target task */
          for(int i = 0; i < Tp->NumPart; i++)
            Tp->PS[i].TargetIndex = Tp->PS[i].submark;

          /* assemble the particles on individual processors (note: IndexList[] becomes temporarily meaningless)  */
          subfind_distribute_particles(SubComm);

          PPS = (PPS_data *)Mem.mymalloc("PPS", Tp->NumPart * sizeof(PPS_data));

          for(int i = 0; i < Tp->NumPart; i++)
            PPS[i].index = i;
        }
      else
        {
          PPS = (PPS_data *)Mem.mymalloc("PPS", Tp->NumPart * sizeof(PPS_data));

          for(int i = 0; i < Tp->NumPart; i++)
            {
              PPS[i].submark = Tp->PS[i].submark;
              PPS[i].index   = i;
            }

          mycxxsort(PPS, PPS + Tp->NumPart, subfind_compare_PPS);
        }

      MPI_Barrier(SubComm);
      double ta = Logs.second();

      subfind_unbind_independent_ones(SingleDomain, count_cand);

      MPI_Barrier(SubComm);
      double tb = Logs.second();

      Mem.myfree(PPS);

      subfind_collective_printf("SUBFIND: root-task=%d: unbinding of independent ones took %g sec\n", ThisTask, Logs.timediff(ta, tb));

      if(mode == COLL_SUBFIND)
        {
          for(int i = 0; i < Tp->NumPart; i++)
            {
              Tp->PS[i].u.s.origintask &= (HIGHBIT - 1); /* clear high bit if set */
              Tp->PS[i].TargetTask  = Tp->PS[i].u.s.origintask;
              Tp->PS[i].TargetIndex = Tp->PS[i].u.s.originindex;
            }

          t0 = Logs.second();
          subfind_distribute_particles(SubComm); /* bring them back to their original processor */
          t1 = Logs.second();

          subfind_collective_printf("SUBFIND: root-task=%d: bringing the independent ones back took %g sec\n", ThisTask,
                                    Logs.timediff(t0, t1));

          /* Since we reestablihed the original order, we can use IndexList[] again */
        }

      /* now mark the bound particles */
      for(int i = 0; i < NumPartGroup; i++)
        if(Tp->PS[IndexList[i]].submark >= 0 && Tp->PS[IndexList[i]].submark < nsubs)
          SFTail[i].index = Tp->PS[IndexList[i]].submark; /* we use this to flag bound parts of substructures */

      for(int i = 0; i < count_cand; i++)
        if(coll_candidates[i].parent == 0)
          coll_candidates[i].parent = -1;
    }
  while(count_leaves > 0);

  /*
   * Now we do the unbinding of the subhalo candidates that contain other subhalo candidates.
   * This will be done with several CPUs if needed.
   */

  double t0 = Logs.second();

  for(int master = 0, nr = 0; master < SubNTask; master++)
    {
      int ncand = count_cand;
      MPI_Bcast(&ncand, sizeof(ncand), MPI_BYTE, master, SubComm);

      for(int k = 0; k < ncand; k++)
        {
          MyLenType len;
          int parent, nsubs;
          if(SubThisTask == master)
            {
              len    = coll_candidates[k].len;
              nsubs  = coll_candidates[k].nsub;
              parent = coll_candidates[k].parent; /* this is here actually the daughter count */
            }

          MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, SubComm);
          MPI_Barrier(SubComm);

          if(parent >= 0)
            {
              MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, SubComm);
              MPI_Bcast(&nsubs, sizeof(nsubs), MPI_BYTE, master, SubComm);

              subfind_collective_printf("SUBFIND: root-task=%d: collective unbinding of nr=%d (%d) of length=%d\n", ThisTask, nr,
                                        nremaining, (int)len);

              nr++;

              LocalLen = 0;

              double tt0 = Logs.second();

              unbind_list = (int *)Mem.mymalloc_movable(&unbind_list, "unbind_list", NumPartGroup * sizeof(int));

              if(SubThisTask != master)
                subfind_poll_for_requests();
              else
                {
                  location p = coll_candidates[k].head;
                  for(int i = 0; i < coll_candidates[k].len; i++)
                    {
                      if(p.index < 0)
                        Terminate("Bummer i=%d \n", i);

                      subfind_distlinklist_add_particle(p);

                      p = subfind_distlinklist_get_next(p);
                    }

                  /* now tell the others to stop polling */
                  for(int i = 0; i < SubNTask; i++)
                    if(i != SubThisTask)
                      MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
                }

              if(LocalLen > NumPartGroup)
                Terminate("LocalLen=%d  > NumPartGroup=%d", LocalLen, NumPartGroup);

              /* rewrite list of group indices to particle indices */
              for(int i = 0; i < LocalLen; i++)
                {
                  unbind_list[i] = IndexList[unbind_list[i]];
                  if(unbind_list[i] < 0 || unbind_list[i] >= Tp->NumPart)
                    Terminate("bad!  unbind_list[i]=%d\n", unbind_list[i]);
                }

              /* mark the one to be unbound in PS[] */
              for(int i = 0; i < Tp->NumPart; i++)
                {
                  Tp->PS[i].u.s.origintask  = SubThisTask;
                  Tp->PS[i].u.s.originindex = i;
                  Tp->PS[i].DomainFlag      = 0;
                }

              for(int i = 0; i < LocalLen; i++)
                Tp->PS[unbind_list[i]].DomainFlag = 1;

              Mem.myfree(unbind_list);

              domain<partset> SubUnbindDomain{SubComm, Tp};

              if(SubUnbindDomain.NumNodes != 0)
                Terminate("SubDomain.NumNodes=%d\n", SubUnbindDomain.NumNodes);

              SubUnbindDomain.domain_decomposition(mode);

              if(mode == COLL_SUBFIND)
                subfind_distribute_particles(SubComm);

              /* refill unbind_list */
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

              LocalLen = subfind_unbind(&SubUnbindDomain, SubComm, unbind_list, LocalLen);

              for(int i = 0; i < Tp->NumPart; i++)
                {
                  Tp->PS[i].DomainFlag  = 0;
                  Tp->PS[i].TargetTask  = Tp->PS[i].u.s.origintask;
                  Tp->PS[i].TargetIndex = Tp->PS[i].u.s.originindex;
                }

              for(int i = 0; i < LocalLen; i++)
                Tp->PS[unbind_list[i]].DomainFlag = 1;

              Mem.myfree(unbind_list);

              SubUnbindDomain.domain_free();

              double ta = Logs.second();
              subfind_distribute_particles(SubComm); /* bring them back to their original processor */
              double tb = Logs.second();

              unbind_list = (int *)Mem.mymalloc_movable(&unbind_list, "unbind_list", NumPartGroup * sizeof(int));

              /* refill unbind_list */
              LocalLen = 0;
              for(int i = 0; i < Tp->NumPart; i++)
                if(Tp->PS[i].DomainFlag)
                  {
                    if(LocalLen >= NumPartGroup)
                      Terminate("LocalLen=%d >= NumPartGroup=%d", LocalLen, NumPartGroup);
                    unbind_list[LocalLen++] = i;
                  }

              subfind_collective_printf("SUBFIND: root-task=%d: bringing the collective ones back took %g sec\n", ThisTask,
                                        Logs.timediff(ta, tb));

              /* go from particle indices back to group indices */
              for(int i = 0; i < LocalLen; i++)
                {
                  unbind_list[i] = Tp->PS[unbind_list[i]].InvIndex;
                  if(unbind_list[i] < 0 || unbind_list[i] >= NumPartGroup)
                    Terminate("ups, bad!  unbind_list[i]=%d   NumPartGroup=%d\n", unbind_list[i], NumPartGroup);
                }

              double tt1 = Logs.second();

              int oldlen = len;

              MPI_Allreduce(&LocalLen, &len, 1, MPI_INT, MPI_SUM, SubComm);

              subfind_collective_printf(
                  "SUBFIND: root-task=%d: collective unbinding of nr=%d (%d) of length=%d, bound length=%d    took %g sec\n", ThisTask,
                  nr - 1, nremaining, oldlen, (int)len, Logs.timediff(tt0, tt1));

              if(len >= All.DesLinkNgb)
                {
                  /* ok, we found a substructure */
                  for(int i = 0; i < LocalLen; i++)
                    SFTail[unbind_list[i]].index = nsubs; /* we use this to flag the substructures */

                  if(SubThisTask == master)
                    coll_candidates[k].bound_length = len;
                }
              else
                {
                  /* bound particle count too low or zero */
                  if(SubThisTask == master)
                    coll_candidates[k].bound_length = 0;
                }

              Mem.myfree(unbind_list);
            }
        }
    }
  double t1 = Logs.second();

  subfind_collective_printf("SUBFIND: root-task=%d: the collective unbinding of remaining halos took %g sec\n", ThisTask,
                            Logs.timediff(t0, t1));

  /* get the total substructure count */
  int countall = 0;
  for(int k = 0; k < count_cand; k++)
    if(coll_candidates[k].bound_length >= All.DesLinkNgb)
      {
        if(coll_candidates[k].len < All.DesLinkNgb)
          Terminate("coll_candidates[k=%d].len=%lld bound=%lld\n", k, (long long)coll_candidates[k].len,
                    (long long)coll_candidates[k].bound_length);

        countall++;
      }

  MPI_Allreduce(MPI_IN_PLACE, &countall, 1, MPI_INT, MPI_SUM, SubComm);

  subfind_collective_printf("SUBFIND: root-task=%d: found %d bound substructures in FoF group of length %lld\n", ThisTask, countall,
                            totgrouplen);

  /* now determine the parent subhalo for each candidate */
  t0 = Logs.second();
  mycxxsort_parallel(coll_candidates, coll_candidates + count_cand, subfind_compare_coll_candidates_boundlength, SubComm);

  coll_cand_dat *tmp_coll_candidates = 0;

  if(SubThisTask == 0)
    tmp_coll_candidates = (coll_cand_dat *)Mem.mymalloc("tmp_coll_candidates", totcand * sizeof(coll_cand_dat));

  int *countlist = (int *)Mem.mymalloc("countlist", SubNTask * sizeof(int));
  int *offset    = (int *)Mem.mymalloc("offset", SubNTask * sizeof(int));

  int count_size = count_cand * sizeof(coll_cand_dat);
  MPI_Allgather(&count_size, 1, MPI_INT, countlist, 1, MPI_INT, SubComm);

  offset[0] = 0;
  for(int i = 1; i < SubNTask; i++)
    offset[i] = offset[i - 1] + countlist[i - 1];

  MPI_Gatherv(coll_candidates, countlist[SubThisTask], MPI_BYTE, tmp_coll_candidates, countlist, offset, MPI_BYTE, 0, SubComm);

  if(SubThisTask == 0)
    {
      for(int k = 0; k < totcand; k++)
        {
          tmp_coll_candidates[k].subnr  = k;
          tmp_coll_candidates[k].parent = 0;
        }

      mycxxsort(tmp_coll_candidates, tmp_coll_candidates + totcand, subfind_compare_coll_candidates_rank);

      for(int k = 0; k < totcand; k++)
        {
          for(int j = k + 1; j < totcand; j++)
            {
              if(tmp_coll_candidates[j].rank > tmp_coll_candidates[k].rank + tmp_coll_candidates[k].len)
                break;

              if(tmp_coll_candidates[k].rank + tmp_coll_candidates[k].len >= tmp_coll_candidates[j].rank + tmp_coll_candidates[j].len)
                {
                  if(tmp_coll_candidates[k].bound_length >= All.DesLinkNgb)
                    tmp_coll_candidates[j].parent = tmp_coll_candidates[k].subnr;
                }
              else
                {
                  Terminate("k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d bound=%d\n", k, countall,
                            (int)tmp_coll_candidates[k].rank, (int)tmp_coll_candidates[k].len,
                            (int)tmp_coll_candidates[k].bound_length, (int)tmp_coll_candidates[j].rank,
                            (int)tmp_coll_candidates[j].len, (int)tmp_coll_candidates[j].bound_length);
                }
            }
        }

      mycxxsort(tmp_coll_candidates, tmp_coll_candidates + totcand, subfind_compare_coll_candidates_subnr);
    }

  MPI_Scatterv(tmp_coll_candidates, countlist, offset, MPI_BYTE, coll_candidates, countlist[SubThisTask], MPI_BYTE, 0, SubComm);

  Mem.myfree(offset);
  Mem.myfree(countlist);

  if(SubThisTask == 0)
    Mem.myfree(tmp_coll_candidates);

  t1 = Logs.second();

  subfind_collective_printf("SUBFIND: root-task=%d: determination of parent subhalo took %g sec (presently allocated %g MB)\n",
                            ThisTask, Logs.timediff(t0, t1), Mem.getAllocatedBytesInMB());

  /* Now let's save some properties of the substructures */
  if(SubThisTask == 0)
    Group[gr].Nsubs = countall;

  t0 = Logs.second();

  unbind_list = (int *)Mem.mymalloc_movable(&unbind_list, "unbind_list", NumPartGroup * sizeof(int));

  for(int master = 0, subnr = 0; master < SubNTask; master++)
    {
      int ncand = count_cand;
      MPI_Bcast(&ncand, sizeof(ncand), MPI_BYTE, master, SubComm);

      for(int k = 0; k < ncand; k++)
        {
          MyLenType len;
          int parent, nsubs;
          if(SubThisTask == master)
            {
              len    = coll_candidates[k].bound_length;
              nsubs  = coll_candidates[k].nsub;
              parent = coll_candidates[k].parent;
            }

          MPI_Bcast(&len, sizeof(len), MPI_BYTE, master, SubComm);
          MPI_Barrier(SubComm);

          if(len > 0)
            {
              MPI_Bcast(&nsubs, sizeof(nsubs), MPI_BYTE, master, SubComm);
              MPI_Bcast(&parent, sizeof(parent), MPI_BYTE, master, SubComm);

              LocalLen = 0;

              if(SubThisTask != master)
                subfind_poll_for_requests();
              else
                {
                  location p = coll_candidates[k].head;
                  for(MyLenType i = 0; i < coll_candidates[k].len; i++)
                    {
                      subfind_distlinklist_add_bound_particles(p, nsubs);
                      p = subfind_distlinklist_get_next(p);
                    }

                  /* now tell the others to stop polling */
                  for(int i = 0; i < SubNTask; i++)
                    if(i != SubThisTask)
                      MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
                }

              int max_nsubhalos;
              MPI_Allreduce(&Nsubhalos, &max_nsubhalos, 1, MPI_INT, MPI_MAX, SubComm);

              if(max_nsubhalos >= MaxNsubhalos)
                {
                  if(ThisTask == 0)
                    warn("Nsubhalos=%d >= MaxNsubhalos=%d", max_nsubhalos, MaxNsubhalos);

                  MaxNsubhalos = 1.25 * max_nsubhalos;

                  Subhalo = (subhalo_properties *)Mem.myrealloc_movable(Subhalo, MaxNsubhalos * sizeof(subhalo_properties));
                }

              for(int i = 0; i < LocalLen; i++)
                {
                  unbind_list[i] = IndexList[unbind_list[i]]; /* move to particle index list */

                  if(unbind_list[i] < 0 || unbind_list[i] >= Tp->NumPart)
                    Terminate("strange");
                }

              int marked = subfind_determine_sub_halo_properties(unbind_list, LocalLen, &Subhalo[Nsubhalos], SubComm);

              for(int i = 0; i < LocalLen; i++)
                {
                  unbind_list[i] = Tp->PS[unbind_list[i]].InvIndex; /* move back to group index list */

                  if(unbind_list[i] < 0 || unbind_list[i] >= NumPartGroup)
                    Terminate("also strange");
                }

              MPI_Allreduce(MPI_IN_PLACE, &marked, 1, MPI_INT, MPI_SUM, SubComm);

              if(SubThisTask == 0)
                {
                  if(subnr == 0)
                    {
                      for(int j = 0; j < 3; j++)
                        {
                          Group[gr].Pos[j]    = Subhalo[Nsubhalos].Pos[j];
                          Group[gr].IntPos[j] = Subhalo[Nsubhalos].IntPos[j];
                        }
                    }

#if defined(SUBFIND_ORPHAN_TREATMENT)
                  Group[gr].LenPrevMostBnd += marked;
#endif

                  Subhalo[Nsubhalos].GroupNr       = GroupNr;
                  Subhalo[Nsubhalos].SubRankInGr   = subnr;
                  Subhalo[Nsubhalos].SubParentRank = parent;

                  Nsubhalos++;
                }

              /* Let's now assign the subhalo number within the group */
              for(int i = 0; i < LocalLen; i++)
                {
                  Tp->PS[IndexList[unbind_list[i]]].SubRankInGr = subnr;
#if defined(MERGERTREE)
                  Tp->PS[IndexList[unbind_list[i]]].SizeOfSubhalo.set(len);
#endif
                }

              subnr++;
            }
        }
    }

  subfind_collective_printf("SUBFIND: root-task=%d: determining substructure properties done\n", ThisTask);

  Mem.myfree(unbind_list);
  Mem.myfree(coll_candidates);
  Mem.myfree(SFPrevLen);
  Mem.myfree(SFLen);
  Mem.myfree(SFTail);
  Mem.myfree(SFNext);
  Mem.myfree(SFHead);
  Mem.myfree(sd);
}

/* This function finds the subhalo candidates (i.e. locally overdense structures bounded by a saddle point).
 * They can be nested inside each other, and will later be subjected to an unbinding procedure.
 */
template <typename partset>
void fof<partset>::subfind_col_find_coll_candidates(long long totgrouplen)
{
  subfind_collective_printf("SUBFIND: root-task=%d: building distributed linked list. (presently allocated %g MB)\n", ThisTask,
                            Mem.getAllocatedBytesInMB());

  double t0 = Logs.second();

  /* Now find the subhalo candidates by building up-link lists for them.
   * We go through the processors in the group in sequence, starting with the first one holding the densest particles according to the
   * sd[] array In each iteration, the processor we currently deal with is called 'master', the others are listening to incoming
   * requests for information.
   */
  for(int master = 0; master < SubNTask; master++)
    {
      double tt0 = Logs.second();
      if(SubThisTask != master)
        subfind_poll_for_requests(); /* if we are not the master task, we react to incoming requests for information */
      else
        {
          /* we go through the sd[] indices stored on the master task, which means we start with the densest particle */
          for(int k = 0; k < NumPartGroup; k++)
            {
              int ngbcount        = sd[k].ngbcount;
              location ngb_index1 = sd[k].ngb_index1;
              location ngb_index2 = sd[k].ngb_index2;

              switch(ngbcount)
                /* treat the different possible cases */
                {
                  case 0: /* this appears to be a lonely maximum -> new group */
                    subfind_distlinklist_set_all(sd[k].index, sd[k].index, sd[k].index, 1, {-1, -1}, sd[k].PrevSizeOfSubhalo);
                    break;

                  case 1: /* the particle is attached to exactly one group */
                    {
                      if(ngb_index1.task < 0 || ngb_index1.task >= SubNTask)
                        Terminate("ngb_index1.task=%d  SubNTask=%d", ngb_index1.task, SubNTask);

                      location head = subfind_distlinklist_get_head(ngb_index1);

                      if(head.index == -1)
                        Terminate("We have a problem!  head=%d for k=%d on task=%d\n", head.index, k, SubThisTask);

                      location tail;
                      int retcode =
                          subfind_distlinklist_get_tail_set_tail_increaselen(head, tail, sd[k].index, sd[k].PrevSizeOfSubhalo);

                      if(!(retcode & 1))
                        subfind_distlinklist_set_headandnext(sd[k].index, head, {-1, -1});
                      if(!(retcode & 2))
                        subfind_distlinklist_set_next(tail, sd[k].index);
                    }
                    break;

                  case 2: /* the particle merges two groups together */
                    {
                      location head, head_attach;

                      if(ngb_index1.task < 0 || ngb_index1.task >= SubNTask)
                        Terminate("ngb_index1.task=%d  SubNTask=%d", ngb_index1.task, SubNTask);

                      if(ngb_index2.task < 0 || ngb_index2.task >= SubNTask)
                        Terminate("ngb_index2.task=%d  SubNTask=%d", ngb_index2.task, SubNTask);

                      if(ngb_index1.task == ngb_index2.task)
                        {
                          subfind_distlinklist_get_two_heads(ngb_index1, ngb_index2, head, head_attach);
                        }
                      else
                        {
                          head        = subfind_distlinklist_get_head(ngb_index1);
                          head_attach = subfind_distlinklist_get_head(ngb_index2);
                        }

                      if(head.index == -1 || head_attach.index == -1)
                        Terminate("We have a problem!  head=%d/%d head_attach=%d/%d for k=%d on task=%d\n", head.task, head.index,
                                  head_attach.task, head_attach.index, k, SubThisTask);

                      if(head != head_attach)
                        {
                          location tail, tail_attach;
                          MyLenType len, len_attach;
                          double prevlen, prevlen_attach;

                          subfind_distlinklist_get_tailandlen(head, tail, len, prevlen);
                          subfind_distlinklist_get_tailandlen(head_attach, tail_attach, len_attach, prevlen_attach);

                          bool swap_len     = false;
                          bool swap_prevlen = false;

                          if(len_attach > len || (len_attach == len && head_attach < head))
                            swap_len = true;

                          if(prevlen > 0 && prevlen_attach > 0 && len >= All.DesLinkNgb && len_attach >= All.DesLinkNgb)
                            {
                              if(prevlen_attach > prevlen || (prevlen_attach == prevlen && swap_len == true))
                                swap_prevlen = true;
                            }
                          else
                            swap_prevlen = swap_len;

                          /* if other group is longer, swap */
                          if(swap_prevlen)
                            {
                              location tmp = head;
                              head         = head_attach;
                              head_attach  = tmp;

                              tmp         = tail;
                              tail        = tail_attach;
                              tail_attach = tmp;

                              MyLenType tmplen = len;
                              len              = len_attach;
                              len_attach       = tmplen;

                              double tmpprevlen = prevlen;
                              prevlen           = prevlen_attach;
                              prevlen_attach    = tmpprevlen;
                            }

                          /* only in case the attached group is long enough we bother to register it
                           as a subhalo candidate */

                          if(len_attach >= All.DesLinkNgb && len >= All.DesLinkNgb)
                            {
                              count_decisions++;

                              if(swap_prevlen != swap_len)
                                {
                                  printf(
                                      "SUBFIND: TASK=%d:  made a different main trunk decision due to previous length: prevlen=%g "
                                      "prevlen_attach=%g   len=%g len_attach=%g\n",
                                      ThisTask, prevlen / len, prevlen_attach / len_attach, (double)len, (double)len_attach);
                                  fflush(stdout);
                                  count_different_decisions++;
                                }

                              if(count_cand < max_coll_candidates)
                                {
                                  coll_candidates[count_cand].len  = len_attach;
                                  coll_candidates[count_cand].head = head_attach;
                                  count_cand++;
                                }
                              else
                                Terminate("Task %d: count=%d, max=%d, npartgroup=%d\n", SubThisTask, count_cand, max_coll_candidates,
                                          NumPartGroup);
                            }

                          /* now join the two groups */
                          subfind_distlinklist_set_tailandlen(head, tail_attach, len + len_attach, prevlen + prevlen_attach);
                          subfind_distlinklist_set_next(tail, head_attach);

                          location ss = head_attach;
                          do
                            {
                              ss = subfind_distlinklist_set_head_get_next(ss, head);
                            }
                          while(ss.index >= 0);
                        }

                      /* finally, attach the particle to 'head' */
                      location tail;
                      int retcode =
                          subfind_distlinklist_get_tail_set_tail_increaselen(head, tail, sd[k].index, sd[k].PrevSizeOfSubhalo);

                      if(!(retcode & 1))
                        subfind_distlinklist_set_headandnext(sd[k].index, head, {-1, -1});
                      if(!(retcode & 2))
                        subfind_distlinklist_set_next(tail, sd[k].index);
                    }
                    break;
                }
            }

          myflush(stdout);

          /* now tell the others to stop polling */
          for(int k = 0; k < SubNTask; k++)
            if(k != SubThisTask)
              MPI_Send(&k, 1, MPI_INT, k, TAG_POLLING_DONE, SubComm);
        }

      MPI_Barrier(SubComm);
      double tt1 = Logs.second();

      subfind_collective_printf("SUBFIND: root-task=%d: ma=%d/%d took %g sec\n", ThisTask, master, SubNTask, Logs.timediff(tt0, tt1));
    }
  double t1 = Logs.second();

  subfind_collective_printf("SUBFIND: root-task=%d: identification of primary coll_candidates took %g sec\n", ThisTask,
                            Logs.timediff(t0, t1));

  /* Add the full group as the final subhalo candidate.
   */
  location head = {-1, -1};
  location prev = {-1, -1};
  for(int master = 0; master < SubNTask; master++)
    {
      if(SubThisTask != master)
        subfind_poll_for_requests();
      else
        {
          for(int i = 0; i < NumPartGroup; i++)
            {
              location index = {SubThisTask, i};

              if(SFHead[i] == index)
                {
                  location tail;
                  MyLenType len;
                  double prevlen;
                  subfind_distlinklist_get_tailandlen(SFHead[i], tail, len, prevlen);
                  location next = subfind_distlinklist_get_next(tail);

                  if(next.index == -1)
                    {
                      if(prev.index < 0)
                        head = index;

                      if(prev.index >= 0)
                        subfind_distlinklist_set_next(prev, index);

                      prev = tail;
                    }
                }
            }

          /* now tell the others to stop polling */
          for(int k = 0; k < SubNTask; k++)
            if(k != SubThisTask)
              MPI_Send(&k, 1, MPI_INT, k, TAG_POLLING_DONE, SubComm);
        }

      MPI_Barrier(SubComm);
      MPI_Bcast(&head, sizeof(head), MPI_BYTE, master, SubComm);
      MPI_Bcast(&prev, sizeof(prev), MPI_BYTE, master, SubComm);
    }

  if(SubThisTask == SubNTask - 1)
    {
      if(count_cand < max_coll_candidates)
        {
          coll_candidates[count_cand].len  = totgrouplen;
          coll_candidates[count_cand].head = head;
          count_cand++;
        }
      else
        Terminate("count_cand=%d >= max_coll_candidates=%d", count_cand, max_coll_candidates);
    }

  subfind_collective_printf("SUBFIND: root-task=%d: adding background as candidate\n", ThisTask);

  /* go through the whole chain once to establish a rank order. For the rank we use SFLen[]
   */
  double ta = Logs.second();

  int master = head.task;

  if(master < 0 || master >= SubNTask)
    Terminate("master=%d  SubNTask=%d\n", master, SubNTask);

  if(SubThisTask != master)
    subfind_poll_for_requests();
  else
    {
      location p     = head;
      MyLenType rank = 0;

      while(p.index >= 0)
        {
          p = subfind_distlinklist_setrank_and_get_next(p, rank);
        }

      /* now tell the others to stop polling */
      for(int i = 0; i < SubNTask; i++)
        if(i != master)
          MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
    }

  MPI_Barrier(SubComm);

  /* for each candidate, we now pull out the rank of its head */
  for(int master = 0; master < SubNTask; master++)
    {
      if(SubThisTask != master)
        subfind_poll_for_requests();
      else
        {
          for(int k = 0; k < count_cand; k++)
            coll_candidates[k].rank = subfind_distlinklist_get_rank(coll_candidates[k].head);

          /* now tell the others to stop polling */
          for(int i = 0; i < SubNTask; i++)
            if(i != SubThisTask)
              MPI_Send(&i, 1, MPI_INT, i, TAG_POLLING_DONE, SubComm);
        }
    }
  MPI_Barrier(SubComm);

  double tb = Logs.second();

  subfind_collective_printf(
      "SUBFIND: root-task=%d: establishing of rank order took %g sec (grouplen=%lld) presently allocated %g MB\n", ThisTask,
      Logs.timediff(ta, tb), totgrouplen, Mem.getAllocatedBytesInMB());
}

template <typename partset>
void fof<partset>::subfind_unbind_independent_ones(domain<partset> *SingleDomain, int count_cand_l)
{
  unbind_list = (int *)Mem.mymalloc("unbind_list", Tp->NumPart * sizeof(int));

  mycxxsort(coll_candidates, coll_candidates + count_cand_l, subfind_compare_coll_candidates_nsubs);

  for(int k = 0, ii = 0; k < count_cand_l; k++)
    if(coll_candidates[k].parent == 0)
      {
        int i = PPS[ii].index;

        while(Tp->PS[i].submark < coll_candidates[k].nsub)
          {
            ii++;
            i = PPS[ii].index;

            if(i >= Tp->NumPart)
              Terminate("i >= NumPart");
          }

        if(Tp->PS[i].submark >= 0 && Tp->PS[i].submark < HIGHBIT)
          {
            int len   = 0;
            int nsubs = Tp->PS[i].submark;

            if(nsubs != coll_candidates[k].nsub)
              Terminate("TASK=%d i=%d k=%d nsubs=%d coll_candidates[k].nsub=%d\n", SubThisTask, i, k, nsubs, coll_candidates[k].nsub);

            while(i < Tp->NumPart)
              {
                if(Tp->PS[i].submark == nsubs)
                  {
                    Tp->PS[i].submark = HIGHBIT;
                    if((Tp->PS[i].u.s.origintask & HIGHBIT) == 0)
                      {
                        unbind_list[len] = i;
                        len++;
                      }
                    ii++;
                    i = PPS[ii].index;
                  }
                else
                  break;
              }

            /* call the serial unbind function */
            len = subfind_unbind(SingleDomain, SingleDomain->Communicator, unbind_list, len);

            if(len >= All.DesLinkNgb)
              {
                /* ok, we found a substructure */
                coll_candidates[k].bound_length = len;

                for(int j = 0; j < len; j++)
                  Tp->PS[unbind_list[j]].submark = nsubs; /* we use this to flag the substructures */
              }
            else
              coll_candidates[k].bound_length = 0;
          }
      }

  Mem.myfree(unbind_list);
}

struct loc_compound0
{
  int index;
  location loc;
  approxlen prevlen;
};

struct loc_compound1
{
  int index;
  location loc;
};

struct loc_compound2
{
  location loc;
  MyLenType len;
  double prevlen;
};

struct loc_compound3
{
  int index;
  MyLenType len;
  location tail;
  double prevlen;
};

struct loc_compound4
{
  int index;
  location head;
  location next;
};

struct loc_compound5
{
  int index;
  MyLenType len;
  location head;
  location tail;
  location next;
  approxlen prevlen;
};

struct loc_compound6
{
  location loc;
  MyLenType len;
};

template <typename partset>
void fof<partset>::subfind_poll_for_requests(void)
{
  int tag;
  do
    {
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, SubComm, &status);

      int source = status.MPI_SOURCE;
      tag        = status.MPI_TAG;

      /* MPI_Get_count(&status, MPI_BYTE, &count); */
      switch(tag)
        {
          case TAG_GET_TWOHEADS:
            {
              int ibuf[2];
              MPI_Recv(ibuf, 2, MPI_INT, source, TAG_GET_TWOHEADS, SubComm, MPI_STATUS_IGNORE);
              location buf[2];
              buf[0] = SFHead[ibuf[0]];
              buf[1] = SFHead[ibuf[1]];
              MPI_Send(buf, 2 * sizeof(location), MPI_BYTE, source, TAG_GET_TWOHEADS_DATA, SubComm);
            }
            break;

          case TAG_SET_NEWTAIL:
            {
              loc_compound0 data;
              MPI_Recv(&data, sizeof(data), MPI_BYTE, source, TAG_SET_NEWTAIL, SubComm, MPI_STATUS_IGNORE);

              int index        = data.index;
              location newtail = data.loc;
              location oldtail = SFTail[index]; /* return old tail */
              SFTail[index]    = newtail;
              SFLen[index]++;
              SFPrevLen[index] += data.prevlen.get();

              if(newtail.task == SubThisTask)
                {
                  SFHead[newtail.index] = {SubThisTask, index};
                  SFNext[newtail.index] = {-1, -1};
                }

              if(oldtail.task == SubThisTask)
                {
                  SFNext[oldtail.index] = newtail;
                }

              MPI_Send(&oldtail, sizeof(location), MPI_BYTE, source, TAG_GET_OLDTAIL, SubComm);
            }
            break;

          case TAG_SET_ALL:
            {
              loc_compound5 data;
              MPI_Recv(&data, sizeof(data), MPI_BYTE, source, TAG_SET_ALL, SubComm, MPI_STATUS_IGNORE);
              int index        = data.index;
              SFLen[index]     = data.len;
              SFHead[index]    = data.head;
              SFTail[index]    = data.tail;
              SFNext[index]    = data.next;
              SFPrevLen[index] = data.prevlen.get();
            }
            break;

          case TAG_GET_TAILANDLEN:
            {
              int index;
              MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
              loc_compound2 data = {SFTail[index], SFLen[index], SFPrevLen[index]};
              MPI_Send(&data, sizeof(data), MPI_BYTE, source, TAG_GET_TAILANDLEN_DATA, SubComm);
            }
            break;

          case TAG_SET_TAILANDLEN:
            {
              loc_compound3 data;
              MPI_Recv(&data, sizeof(data), MPI_BYTE, source, TAG_SET_TAILANDLEN, SubComm, MPI_STATUS_IGNORE);
              int index        = data.index;
              SFTail[index]    = data.tail;
              SFLen[index]     = data.len;
              SFPrevLen[index] = data.prevlen;
            }
            break;

          case TAG_SET_HEADANDNEXT:
            {
              loc_compound4 data;
              MPI_Recv(&data, sizeof(data), MPI_BYTE, source, TAG_SET_HEADANDNEXT, SubComm, MPI_STATUS_IGNORE);
              int index     = data.index;
              SFHead[index] = data.head;
              SFNext[index] = data.next;
            }
            break;

          case TAG_SET_NEXT:
            {
              loc_compound1 data;
              MPI_Recv(&data, sizeof(data), MPI_BYTE, source, TAG_SET_NEXT, SubComm, MPI_STATUS_IGNORE);
              int index     = data.index;
              SFNext[index] = data.loc;
            }
            break;

          case TAG_SETHEADGETNEXT:
            {
              loc_compound1 data;
              MPI_Recv(&data, sizeof(data), MPI_BYTE, source, TAG_SETHEADGETNEXT, SubComm, MPI_STATUS_IGNORE);
              int index     = data.index;
              location head = data.loc;
              location next;
              int task;
              do
                {
                  SFHead[index] = head;
                  next          = SFNext[index];
                  task          = next.task;
                  index         = next.index;
                }
              while(next.index >= 0 && task == SubThisTask);
              MPI_Send(&next, sizeof(location), MPI_BYTE, source, TAG_SETHEADGETNEXT_DATA, SubComm);
            }
            break;

          case TAG_GET_NEXT:
            {
              int index;
              MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
              MPI_Send(&SFNext[index], sizeof(location), MPI_BYTE, source, TAG_GET_NEXT_DATA, SubComm);
            }
            break;

          case TAG_GET_HEAD:
            {
              int index;
              MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
              MPI_Send(&SFHead[index], sizeof(location), MPI_BYTE, source, TAG_GET_HEAD_DATA, SubComm);
            }
            break;

          case TAG_ADD_PARTICLE:
            {
              int index;
              MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
              if(SFTail[index].index < 0) /* consider only particles not already in substructures */
                {
                  unbind_list[LocalLen] = index;
                  if(index >= NumPartGroup)
                    Terminate("What: index=%d NumPartGroup=%d\n", index, NumPartGroup);
                  LocalLen++;
                }
            }
            break;

          case TAG_MARK_PARTICLE:
            {
              int ibuf[3];
              MPI_Recv(ibuf, 3, MPI_INT, source, TAG_MARK_PARTICLE, SubComm, MPI_STATUS_IGNORE);
              int index   = ibuf[0];
              int target  = ibuf[1];
              int submark = ibuf[2];

              if(Tp->PS[IndexList[index]].submark != HIGHBIT)
                Terminate("TasK=%d i=%d P[i].submark=%d?\n", SubThisTask, IndexList[index], Tp->PS[IndexList[index]].submark);

              Tp->PS[IndexList[index]].TargetTask = target;
              Tp->PS[IndexList[index]].submark    = submark;
            }
            break;

          case TAG_ADDBOUND:
            {
              int ibuf[2];
              MPI_Recv(ibuf, 2, MPI_INT, source, TAG_ADDBOUND, SubComm, &status);
              int index = ibuf[0];
              int nsub  = ibuf[1];
              if(SFTail[index].index == nsub) /* consider only particles in this substructure */
                {
                  unbind_list[LocalLen] = index;
                  LocalLen++;
                }
            }
            break;

          case TAG_SETRANK:
            {
              loc_compound6 data;
              MPI_Recv(&data, sizeof(data), MPI_BYTE, source, TAG_SETRANK, SubComm, MPI_STATUS_IGNORE);
              int index      = data.loc.index;
              MyLenType rank = data.len;
              location next;
              do
                {
                  SFLen[index] = rank++;
                  next         = SFNext[index];
                  if(next.index < 0)
                    break;
                  index = next.index;
                }
              while(next.task == SubThisTask);
              data.loc = next;
              data.len = rank;
              MPI_Send(&data, sizeof(data), MPI_BYTE, source, TAG_SETRANK_OUT, SubComm);
            }
            break;

          case TAG_GET_RANK:
            {
              int index;
              MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
              MyLenType rank = SFLen[index];
              MPI_Send(&rank, sizeof(MyLenType), MPI_BYTE, source, TAG_GET_RANK_DATA, SubComm);
            }
            break;

          case TAG_POLLING_DONE:
            {
              int index;
              MPI_Recv(&index, 1, MPI_INT, source, tag, SubComm, &status);
            }
            break;

          default:
            Terminate("tag not present in the switch");
            break;
        }
    }
  while(tag != TAG_POLLING_DONE);
}

template <typename partset>
location fof<partset>::subfind_distlinklist_setrank_and_get_next(location loc, MyLenType &rank)
{
  int task = loc.task;
  int i    = loc.index;

  location next;

  if(SubThisTask == task)
    {
      SFLen[i] = rank++;
      next     = SFNext[i];
    }
  else
    {
      loc_compound6 data = {loc, rank};
      MPI_Send(&data, sizeof(data), MPI_BYTE, task, TAG_SETRANK, SubComm);
      MPI_Recv(&data, sizeof(data), MPI_BYTE, task, TAG_SETRANK_OUT, SubComm, MPI_STATUS_IGNORE);
      next = data.loc;
      rank = data.len;
    }
  return next;
}

template <typename partset>
location fof<partset>::subfind_distlinklist_set_head_get_next(location loc, location head)
{
  int task = loc.task;
  int i    = loc.index;

  location next;

  if(SubThisTask == task)
    {
      SFHead[i] = head;
      next      = SFNext[i];
    }
  else
    {
      loc_compound1 data = {i, head};
      MPI_Send(&data, sizeof(data), MPI_BYTE, task, TAG_SETHEADGETNEXT, SubComm);
      MPI_Recv(&next, sizeof(location), MPI_BYTE, task, TAG_SETHEADGETNEXT_DATA, SubComm, MPI_STATUS_IGNORE);
    }

  return next;
}

template <typename partset>
void fof<partset>::subfind_distlinklist_set_next(location loc, location next)
{
  int task = loc.task;
  int i    = loc.index;

  if(SubThisTask == task)
    {
      SFNext[i] = next;
    }
  else
    {
      loc_compound1 data = {i, next};
      MPI_Send(&data, sizeof(data), MPI_BYTE, task, TAG_SET_NEXT, SubComm);
    }
}

template <typename partset>
void fof<partset>::subfind_distlinklist_add_particle(location loc)
{
  int task = loc.task;
  int i    = loc.index;

  if(SubThisTask == task)
    {
      if(SFTail[i].index < 0) /* consider only particles not already in substructures */
        {
          if(i >= NumPartGroup)
            Terminate("What: index=%d NumPartGroup=%d\n", i, NumPartGroup);

          unbind_list[LocalLen] = i;
          LocalLen++;
        }
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_ADD_PARTICLE, SubComm);
    }
}

template <typename partset>
void fof<partset>::subfind_distlinklist_mark_particle(location loc, int target, int submark)
{
  int task = loc.task;
  int i    = loc.index;

  if(SubThisTask == task)
    {
      if(Tp->PS[IndexList[i]].submark != HIGHBIT)
        Terminate("Tas=%d i=%d PS[i].submark=%d?\n", SubThisTask, i, Tp->PS[IndexList[i]].submark);

      Tp->PS[IndexList[i]].TargetTask = target;
      Tp->PS[IndexList[i]].submark    = submark;
    }
  else
    {
      int ibuf[3] = {i, target, submark};
      MPI_Send(ibuf, 3, MPI_INT, task, TAG_MARK_PARTICLE, SubComm);
    }
}

template <typename partset>
void fof<partset>::subfind_distlinklist_add_bound_particles(location loc, int nsub)
{
  int task = loc.task;
  int i    = loc.index;

  if(SubThisTask == task)
    {
      if(SFTail[i].index == nsub) /* consider only particles not already in substructures */
        {
          unbind_list[LocalLen] = i;
          LocalLen++;
        }
    }
  else
    {
      int ibuf[2] = {i, nsub};
      MPI_Send(ibuf, 2, MPI_INT, task, TAG_ADDBOUND, SubComm);
    }
}

template <typename partset>
location fof<partset>::subfind_distlinklist_get_next(location loc)
{
  int task = loc.task;
  int i    = loc.index;

  location next;

  if(SubThisTask == task)
    {
      next = SFNext[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_NEXT, SubComm);
      MPI_Recv(&next, sizeof(location), MPI_BYTE, task, TAG_GET_NEXT_DATA, SubComm, MPI_STATUS_IGNORE);
    }

  return next;
}

template <typename partset>
MyLenType fof<partset>::subfind_distlinklist_get_rank(location loc)
{
  int task = loc.task;
  int i    = loc.index;

  MyLenType rank;

  if(SubThisTask == task)
    {
      rank = SFLen[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_RANK, SubComm);
      MPI_Recv(&rank, sizeof(MyLenType), MPI_BYTE, task, TAG_GET_RANK_DATA, SubComm, MPI_STATUS_IGNORE);
    }

  return rank;
}

template <typename partset>
location fof<partset>::subfind_distlinklist_get_head(location loc)
{
  int task = loc.task;
  int i    = loc.index;

  location head;

  if(SubThisTask == task)
    {
      head = SFHead[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_HEAD, SubComm);
      MPI_Recv(&head, sizeof(location), MPI_BYTE, task, TAG_GET_HEAD_DATA, SubComm, MPI_STATUS_IGNORE);
    }

  return head;
}

template <typename partset>
void fof<partset>::subfind_distlinklist_get_two_heads(location ngb_index1, location ngb_index2, location &head, location &head_attach)
{
  if(ngb_index1.task != ngb_index2.task)
    Terminate("ngb_index1.task != ngb_index2.task");

  int task = ngb_index1.task;
  int i1   = ngb_index1.index;
  int i2   = ngb_index2.index;

  if(SubThisTask == task)
    {
      head        = SFHead[i1];
      head_attach = SFHead[i2];
    }
  else
    {
      int ibuf[2] = {i1, i2};
      MPI_Send(ibuf, 2, MPI_INT, task, TAG_GET_TWOHEADS, SubComm);
      location buf[2];
      MPI_Recv(buf, 2 * sizeof(location), MPI_BYTE, task, TAG_GET_TWOHEADS_DATA, SubComm, MPI_STATUS_IGNORE);
      head        = buf[0];
      head_attach = buf[1];
    }
}

template <typename partset>
void fof<partset>::subfind_distlinklist_set_headandnext(location loc, location head, location next)
{
  int task = loc.task;
  int i    = loc.index;

  if(SubThisTask == task)
    {
      SFHead[i] = head;
      SFNext[i] = next;
    }
  else
    {
      loc_compound4 data = {i, head, next};
      MPI_Send(&data, sizeof(data), MPI_BYTE, task, TAG_SET_HEADANDNEXT, SubComm);
    }
}

template <typename partset>
int fof<partset>::subfind_distlinklist_get_tail_set_tail_increaselen(location loc, location &tail, location newtail, approxlen prevlen)
{
  int task = loc.task;
  int i    = loc.index;

  int retcode = 0;

  if(SubThisTask == task)
    {
      location oldtail = SFTail[i];
      SFTail[i]        = newtail;
      SFLen[i]++;
      SFPrevLen[i] += prevlen.get();
      tail = oldtail;

      if(newtail.task == SubThisTask)
        {
          SFHead[newtail.index] = loc;
          SFNext[newtail.index] = {-1, -1};
          retcode |= 1;
        }

      if(oldtail.task == SubThisTask)
        {
          SFNext[oldtail.index] = newtail;
          retcode |= 2;
        }
    }
  else
    {
      loc_compound0 data = {i, newtail, prevlen};
      MPI_Send(&data, sizeof(data), MPI_BYTE, task, TAG_SET_NEWTAIL, SubComm);
      location oldtail;
      MPI_Recv(&oldtail, sizeof(location), MPI_BYTE, task, TAG_GET_OLDTAIL, SubComm, MPI_STATUS_IGNORE);
      tail = oldtail;

      if(newtail.task == task)
        retcode |= 1;
      if(oldtail.task == task)
        retcode |= 2;
    }

  return retcode;
}

template <typename partset>
void fof<partset>::subfind_distlinklist_set_tailandlen(location loc, location tail, MyLenType len, double prevlen)
{
  int task = loc.task;
  int i    = loc.index;

  if(SubThisTask == task)
    {
      SFTail[i]    = tail;
      SFLen[i]     = len;
      SFPrevLen[i] = prevlen;
    }
  else
    {
      loc_compound3 data = {i, len, tail, prevlen};
      MPI_Send(&data, sizeof(data), MPI_BYTE, task, TAG_SET_TAILANDLEN, SubComm);
    }
}

template <typename partset>
void fof<partset>::subfind_distlinklist_get_tailandlen(location loc, location &tail, MyLenType &len, double &prevlen)
{
  int task = loc.task;
  int i    = loc.index;

  if(SubThisTask == task)
    {
      tail    = SFTail[i];
      len     = SFLen[i];
      prevlen = SFPrevLen[i];
    }
  else
    {
      MPI_Send(&i, 1, MPI_INT, task, TAG_GET_TAILANDLEN, SubComm);

      loc_compound2 data;
      MPI_Recv(&data, sizeof(data), MPI_BYTE, task, TAG_GET_TAILANDLEN_DATA, SubComm, MPI_STATUS_IGNORE);
      tail    = data.loc;
      len     = data.len;
      prevlen = data.prevlen;
    }
}

template <typename partset>
void fof<partset>::subfind_distlinklist_set_all(location loc, location head, location tail, MyLenType len, location next,
                                                approxlen prevlen)
{
  int task = loc.task;
  int i    = loc.index;

  if(SubThisTask == task)
    {
      SFHead[i]    = head;
      SFTail[i]    = tail;
      SFNext[i]    = next;
      SFLen[i]     = len;
      SFPrevLen[i] = prevlen.get();
    }
  else
    {
      loc_compound5 data = {i, len, head, tail, next, prevlen};
      MPI_Send(&data, sizeof(data), MPI_BYTE, task, TAG_SET_ALL, SubComm);
    }
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
