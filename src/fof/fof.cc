/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file fof.cc
 *
 *  \brief main routines of parallel FoF group finder
 */

#include "gadgetconfig.h"

#ifdef FOF

#include <mpi.h>

#include <algorithm>
#include <climits>
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
#include "../fof/fof_io.h"
#include "../gravtree/gravtree.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../sort/cxxsort.h"
#include "../sort/parallel_sort.h"
#include "../sort/peano.h"
#include "../subfind/subfind.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

using namespace std;

/*! Computation of a FOF group catalogue.
 *
 * \param num
 * If called with -1 as argument, only FOF is carried out and no group catalogs are saved to disk, for
 *  >= 0, the code will store the group/subhalo catalogs, and bring the particles into output order.
 * In this case, the calling routine (which is normally write_snapshot()) will need to free PS[] and bring
 * the particles back into the original order.
 */
template <typename partset>
void fof<partset>::fof_fof(int num, const char *grpcat_basename, const char *grpcat_dirbasename, double inner_distance)
{
  TIMER_START(CPU_FOF);

  mpi_printf("FOF: Begin to compute FoF group catalogue...  (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());

  double ta = Logs.second();

  /* determine linking length */
  Tp->LinkL = fof_get_comoving_linking_length();
  mpi_printf("FOF: Comoving linking length: %g\n", Tp->LinkL);

  /* allocate link lists for book-keeping local particle sets */
#if defined(LIGHTCONE_PARTICLES_GROUPS)
  Tp->DistanceOrigin = (double *)Mem.mymalloc("DistanceOrigin", Tp->NumPart * sizeof(double));
#endif

  Tp->MinID     = (MyIDStorage *)Mem.mymalloc("MinID", Tp->NumPart * sizeof(MyIDStorage));  // smallest particle ID within FOF group
  Tp->MinIDTask = (int *)Mem.mymalloc("MinIDTask", Tp->NumPart * sizeof(int));              // processor on which this ID is stored
  Tp->Head = (int *)Mem.mymalloc("Head", Tp->NumPart * sizeof(int));  // first particle in chaining list if local FOF group segment
  Tp->Next = (int *)Mem.mymalloc("Next", Tp->NumPart * sizeof(int));  // next particle in chaining list
  Tp->Tail = (int *)Mem.mymalloc("Tail", Tp->NumPart * sizeof(int));  // points to last particle in chaining list
  Tp->Len  = (int *)Mem.mymalloc("Len", Tp->NumPart * sizeof(int));  // length of local FOF group segment (note: 32 bit enough even for
                                                                     // huge groups because they are split across processors)

  int *numpart_list = (int *)Mem.mymalloc("numpart_list", NTask * sizeof(int));

  MPI_Allgather(&Tp->NumPart, 1, MPI_INT, numpart_list, 1, MPI_INT, Communicator);

  long long NumPartTot = 0;
  for(int i = 0; i < NTask; i++)
    NumPartTot += numpart_list[i];

  mpi_printf("FOF: NumPartTot=%lld\n", NumPartTot);

  if(NumPartTot > (long long)ID_MAX)
    Terminate("The chosen ID data type is not sufficiently big to store unique IDs for NumPartTot=%lld particles\n", NumPartTot);

  MyIDType id = 0;
  for(int i = 0; i < ThisTask; i++)
    id += numpart_list[i];

  Mem.myfree(numpart_list);

  /* initialize link-lists, each particle is in a group of its own initially */
  for(int i = 0; i < Tp->NumPart; i++)
    {
      Tp->Head[i] = Tp->Tail[i] = i;
      Tp->Len[i]                = 1;
      Tp->Next[i]               = -1;
      Tp->MinID[i].set(id++);  // we use new IDs here instead of P[].ID to make sure that also for lightcone group finding MinID is
                               // unique for all box replicas
      Tp->MinIDTask[i] = ThisTask;

#if defined(LIGHTCONE_PARTICLES_GROUPS)
      Tp->DistanceOrigin[i] = fof_distance_to_origin(i);
      Tp->P[i].setFlagSaveDistance();
#endif
    }

    /* make a list with the particles that are of the primary link type(s) */
#ifdef LEAN
  int npart       = Tp->NumPart;
  int *targetlist = NULL;
#else
  int npart       = 0;
  int *targetlist = (int *)Mem.mymalloc("targetlist", Tp->NumPart * sizeof(int));
  for(int i = 0; i < Tp->NumPart; i++)
    if(is_type_primary_link_type(Tp->P[i].getType()))
      targetlist[npart++] = i;
#endif

  /* build neighbour tree with primary link particles only */
  FoFNgbTree.treeallocate(Tp->NumPart, Tp, FoFDomain);
  FoFNgbTree.treebuild(npart, targetlist);

  /* call routine to find groups made up of the primary particle types */
  double cputime = fof_find_groups();
  mpi_printf("FOF: primary group finding took = %g sec\n", cputime);

  /* call routine to attach secondary particles/cells to primary groups */
  cputime = fof_find_nearest_dmparticle();
  mpi_printf("FOF: attaching gas and star particles to nearest dm particles took = %g sec\n", cputime);

  /* free some arrays that are not needed any more */
  FoFNgbTree.treefree();
#ifndef LEAN
  Mem.myfree(targetlist);
#endif
  Mem.myfree(Tp->Len);
  Mem.myfree(Tp->Tail);
  Mem.myfree(Tp->Next);

#if defined(LIGHTCONE_PARTICLES_GROUPS)
  double save_distance = 0;
  if(inner_distance > 0)
    save_distance = inner_distance + Tp->LinkL;
#endif

  double t0 = Logs.second();

  /* transfer the still required link list information to a somewhat smaller structure, "FOF_PList"
   * particles that are in the same group have the same MinID. The MinIDTask variable informs about the processor, on which
   * the particle which has ID=MinID is stored (and this particle is also member of the group).
   */
  FOF_PList = (fof_particle_list *)Mem.mymalloc_movable(&FOF_PList, "FOF_PList", Tp->NumPart * sizeof(fof_particle_list));
  for(int i = 0; i < Tp->NumPart; i++)
    {
      FOF_PList[i].MinID     = Tp->MinID[Tp->Head[i]];
      FOF_PList[i].MinIDTask = Tp->MinIDTask[Tp->Head[i]];
      FOF_PList[i].Pindex    = i;

#if defined(LIGHTCONE_PARTICLES_GROUPS)
      FOF_PList[i].DistanceOrigin = Tp->DistanceOrigin[Tp->Head[i]];

      if(Tp->DistanceOrigin[Tp->Head[i]] < save_distance)
        Tp->P[i].clearFlagSaveDistance();
#endif
    }

  /* free the rest of the original link lists */
  Mem.myfree_movable(Tp->Head);
  Mem.myfree_movable(Tp->MinIDTask);
  Mem.myfree_movable(Tp->MinID);
#if defined(LIGHTCONE_PARTICLES_GROUPS)
  Mem.myfree_movable(Tp->DistanceOrigin);
#endif

  /* Allocate a list of group pieces in FOF_Glist, one entry for each local group segment.
   * If a group is split across processor boundaries, it will appear on each of the processors with
   * a segment.
   */
  FOF_GList = (fof_group_list *)Mem.mymalloc_movable(&FOF_GList, "FOF_GList", sizeof(fof_group_list) * Tp->NumPart);

  fof_compile_catalogue(inner_distance);

  /* now determine total group number, largest group size, and output some log messages */
  double t1 = Logs.second();
  mpi_printf("FOF: compiling local group data and catalogue took = %g sec\n", Logs.timediff(t0, t1));

  sumup_large_ints(1, &Ngroups, &TotNgroups, Communicator);
  sumup_longs(1, &Nids, &TotNids, Communicator);

  /* determine the largest group size */
  long long largestgroup = 0;
  if(TotNgroups > 0)
    {
      long long largestloc = 0;

      for(int i = 0; i < NgroupsExt; i++)
        if(FOF_GList[i].Count > largestloc)
          largestloc = FOF_GList[i].Count;
      MPI_Allreduce(&largestloc, &largestgroup, 1, MPI_LONG_LONG, MPI_MAX, Communicator);
    }

  mpi_printf("FOF: Total number of FOF groups with at least %d particles: %lld\n", FOF_GROUP_MIN_LEN, TotNgroups);
  mpi_printf("FOF: Largest FOF group has %lld particles.\n", largestgroup);
  mpi_printf("FOF: Total number of particles in FOF groups: %lld\n", TotNids);

  t0 = Logs.second();

  /* now allocate some storage for the group catalogue, and begin to fill it in with partial
   * group properties
   */
  MaxNgroups = NgroupsExt;
  Group      = (group_properties *)Mem.mymalloc_movable(&Group, "Group", sizeof(group_properties) * MaxNgroups);

  mpi_printf("FOF: group properties are now allocated.. (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());

  /* Sort FOF_GList according to MinID. We are going to match this with FOF_PList (still ordered according to MinID) to get
   * access to the particles making up each group piece.
   */
  mycxxsort(FOF_GList, FOF_GList + NgroupsExt, fof_compare_FOF_GList_MinID);

  /* compute partial group properties for each local group segment */
  long long count_nids = 0;
  for(int i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID.get() < FOF_GList[i].MinID.get())
        {
          start++;
          if(start > Tp->NumPart)
            Terminate("start > Tp->NumPart");
        }

      if(FOF_PList[start].MinID.get() != FOF_GList[i].MinID.get())
        Terminate("ID mismatch");

      int lenloc = 0;
      for(lenloc = 0; start + lenloc < Tp->NumPart;)
        if(FOF_PList[start + lenloc].MinID.get() == FOF_GList[i].MinID.get())
          lenloc++;
        else
          break;

      Group[i].MinID     = FOF_GList[i].MinID.get();
      Group[i].MinIDTask = FOF_GList[i].MinIDTask;
      Group[i].Len       = FOF_GList[i].Count;

      /* calculate velocity dispersion etc for a local group segment */
      fof_compute_group_properties(i, start, lenloc);

      start += lenloc;
      count_nids += lenloc;
    }

  Mem.myfree_movable(FOF_GList);
  FOF_GList = NULL;

  /* do a sanity check */
  long long totNids;
  sumup_longs(1, &count_nids, &totNids, Communicator);
  if(totNids != TotNids)
    Terminate("Task=%d Nids=%lld count_nids=%lld totNids=%lld TotNids=%lld\n", ThisTask, Nids, count_nids, totNids, TotNids);

  /* add in group properties for each external group segment */
  fof_add_in_properties_of_group_segments();

  t1 = Logs.second();
  mpi_printf("FOF: computation of group properties took = %g sec\n", Logs.timediff(t0, t1));

  /* now we assign group numbers */
  fof_assign_group_numbers();

  Mem.myfree_movable(FOF_PList);
  FOF_PList = NULL;

  /* finalize the computation of the properties of the groups, and prune the list to just contain the main (local) segment */
  fof_finish_group_properties();

  /* Sort the groups in parallel according to group-number, which will be our output order. */
  mycxxsort_parallel(Group, Group + Ngroups, fof_compare_Group_GroupNr, Communicator);

  fof_assign_group_offset();

  double tb = Logs.second();

  mpi_printf("FOF: Finished computing FoF groups.  Complete work took %g sec  (presently allocated=%g MB)\n", Logs.timediff(ta, tb),
             Mem.getAllocatedBytesInMB());

#ifdef SUBFIND
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);

      char catname[MAXLEN_PATH_EXTRA];
      snprintf(catname, MAXLEN_PATH_EXTRA, "%s_subhalo_tab", grpcat_basename);

      subfind_find_subhalos(num, catname, grpcat_dirbasename);

      TIMER_START(CPU_FOF);
    }
#else
  Nsubhalos    = 0;
  TotNsubhalos = 0;
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);
      TIMER_START(CPU_SNAPSHOT);

      fof_io<partset> FoF_IO{this, this->Communicator, All.SnapFormat};

      char catname[MAXLEN_PATH_EXTRA];
      snprintf(catname, MAXLEN_PATH_EXTRA, "%s_tab", grpcat_basename);

      FoF_IO.fof_subfind_save_groups(num, catname, grpcat_dirbasename);

      TIMER_STOP(CPU_SNAPSHOT);
      TIMER_START(CPU_FOF);
    }
#endif

  Mem.myfree_movable(Group);

  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);
      TIMER_START(CPU_SNAPSHOT);

      /* now distribute the particles into output order */
      t0 = Logs.second();
      /* distribute particles such that FOF groups (and SUBFIND halos) will appear in consecutive way in snapshot files */
      fof_prepare_output_order();
      t1 = Logs.second();
      mpi_printf("FOF: preparing output order of particles took %g sec\n", Logs.timediff(t0, t1));

      TIMER_STOP(CPU_SNAPSHOT);
      TIMER_START(CPU_FOF);
    }

  TIMER_STOP(CPU_FOF);
}

template <typename partset>
void fof<partset>::fof_prepare_output_order(void)
{
  int ntype[NTYPES];
  for(int i = 0; i < NTYPES; i++)
    ntype[i] = 0;

  for(int i = 0; i < Tp->NumPart; i++)
    {
#ifndef LEAN
      Tp->PS[i].Type = Tp->P[i].getType();
#endif
      ntype[Tp->P[i].getType()]++;
    }

#if defined(MERGERTREE) && defined(SUBFIND)

  /* we determine a continuously increasing subhalo number (not starting again at zero in every Group), and also assign the ranking
   * within a subhalo
   */
  mycxxsort_parallel(Tp->PS, Tp->PS + Tp->NumPart, fof_compare_subfind_data_GroupNr_SubRankInNr_BndEgy, Communicator);

  int *num_list = (int *)Mem.mymalloc("num_list", NTask * sizeof(int));
  MPI_Allgather(&Tp->NumPart, 1, MPI_INT, num_list, 1, MPI_INT, Communicator);

  int prev_non_empty_task = ThisTask - 1;

  while(prev_non_empty_task >= 0)
    if(num_list[prev_non_empty_task] == 0)
      prev_non_empty_task--;
    else
      break;

  /* obtain last element of each processor */
  subfind_data *aux_last_element = (subfind_data *)Mem.mymalloc("aux_last_element", NTask * sizeof(subfind_data));
  MPI_Allgather(&Tp->PS[Tp->NumPart > 0 ? Tp->NumPart - 1 : 0], sizeof(subfind_data), MPI_BYTE, aux_last_element, sizeof(subfind_data),
                MPI_BYTE, Communicator);

  int first_index = INT_MAX;
  int first_task  = NTask;
  /* find out the very first particle in a subhalo (which is not necessarily in the first FOF group) */
  for(int i = 0; i < Tp->NumPart; i++)
    if(Tp->PS[i].GroupNr.get() < HALONR_MAX && Tp->PS[i].SubRankInGr < INT_MAX) /* particle is in a subhalo */
      {
        first_index = i;
        break;
      }

  int *first_list = (int *)Mem.mymalloc("first_list", NTask * sizeof(int));
  MPI_Allgather(&first_index, 1, MPI_INT, first_list, 1, MPI_INT, Communicator);
  for(int n = 0; n < NTask; n++)
    {
      if(first_list[n] != INT_MAX)
        {
          first_index = first_list[n];
          first_task  = n;
          break;
        }
    }
  Mem.myfree(first_list);

  long long subnr = 0;
  long long rank  = 0;

  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(Tp->PS[i].GroupNr.get() < HALONR_MAX && Tp->PS[i].SubRankInGr < INT_MAX) /* particle is in a subhalo */
        {
          if(i == 0)
            {
              if(prev_non_empty_task >= 0)
                {
                  if(Tp->PS[i].GroupNr.get() != aux_last_element[prev_non_empty_task].GroupNr.get() ||
                     Tp->PS[i].SubRankInGr !=
                         aux_last_element[prev_non_empty_task].SubRankInGr) /* we are the first particle of a new subhalo */
                    {
                      if(ThisTask != first_task || i != first_index)  // to prevent that we start a new subhalo for the very first one
                        {
                          subnr++;
                          rank = 0;
                        }
                    }
                }
            }
          else if(Tp->PS[i].GroupNr.get() != Tp->PS[i - 1].GroupNr.get() ||
                  Tp->PS[i].SubRankInGr != Tp->PS[i - 1].SubRankInGr) /* we are the first particle of a new subhalo */
            {
              if(ThisTask != first_task || i != first_index)  // to prevent that we start a new subhalo for the very first one
                {
                  subnr++;
                  rank = 0;
                }
            }

          Tp->PS[i].SubhaloNr.set(subnr);
          Tp->PS[i].RankInSubhalo.set(rank++);

          if(subnr < 0 || subnr >= TotNsubhalos)
            Terminate("i=%d  NumPart=%d  subnr=%lld  PS[i].SubhaloNr.get()=%lld >= TotNsubhalos=%lld", i, Tp->NumPart, subnr,
                      (long long)Tp->PS[i].SubhaloNr.get(), TotNsubhalos);
        }
      else
        {
          Tp->PS[i].SubhaloNr.set(HALONR_MAX);
          Tp->PS[i].RankInSubhalo.set(INT_MAX);
        }
    }

  long long *subnr_list = (long long *)Mem.mymalloc("subnr_list", NTask * sizeof(long long));
  MPI_Allgather(&subnr, 1, MPI_LONG_LONG, subnr_list, 1, MPI_LONG_LONG, Communicator);

  long long subnr_prev = 0;
  for(int i = 0; i < ThisTask; i++)
    subnr_prev += subnr_list[i];

  for(int i = 0; i < Tp->NumPart; i++)
    if(Tp->PS[i].GroupNr.get() < HALONR_MAX && Tp->PS[i].SubRankInGr < INT_MAX) /* particle is in a subhalo */
      Tp->PS[i].SubhaloNr.set(Tp->PS[i].SubhaloNr.get() + subnr_prev);

  /* obtain previous element of each processor */
  long long *rank_list = (long long *)Mem.mymalloc("rank_list", NTask * sizeof(long long));
  MPI_Allgather(&rank, 1, MPI_LONG_LONG, rank_list, 1, MPI_LONG_LONG, Communicator);

  if(prev_non_empty_task >= 0)
    {
      long long rank = rank_list[prev_non_empty_task];

      for(int i = 0; i < Tp->NumPart; i++)
        {
          if(Tp->PS[i].GroupNr.get() < HALONR_MAX && Tp->PS[i].SubRankInGr < INT_MAX)
            {
              if(i == 0)
                {
                  if(prev_non_empty_task >= 0)
                    if(Tp->PS[i].GroupNr.get() != aux_last_element[prev_non_empty_task].GroupNr.get() ||
                       Tp->PS[i].SubRankInGr !=
                           aux_last_element[prev_non_empty_task].SubRankInGr) /* we are the first particle of a new subhalo */
                      break;
                }
              else if(Tp->PS[i].GroupNr.get() != Tp->PS[i - 1].GroupNr.get() ||
                      Tp->PS[i].SubRankInGr != Tp->PS[i - 1].SubRankInGr) /* we are the first particle of a new subhalo */
                break;

              Tp->PS[i].RankInSubhalo.set(rank++);
            }
        }
    }

  Mem.myfree(rank_list);
  Mem.myfree(subnr_list);
  Mem.myfree(aux_last_element);
  Mem.myfree(num_list);

  /* now bring back into starting order */
  mycxxsort_parallel(Tp->PS, Tp->PS + Tp->NumPart, fof_compare_subfind_data_OriginTask_OriginIndex, Communicator);
#endif

  /* note: the following will destroy the value of the Potential, which is not needed any more at this point */
  for(int i = 0; i < Tp->NumPart; i++)
    {
#if defined(RANDOMIZE_DOMAINCENTER) && defined(PERIODIC)
      Tp->PS[i].u.Key =
          peano_hilbert_key(Tp->P[i].IntPos[0] - Tp->CurrentShiftVector[0], Tp->P[i].IntPos[1] - Tp->CurrentShiftVector[1],
                            Tp->P[i].IntPos[2] - Tp->CurrentShiftVector[2], BITS_FOR_POSITIONS);
#else
      Tp->PS[i].u.Key = peano_hilbert_key(Tp->P[i].IntPos[0], Tp->P[i].IntPos[1], Tp->P[i].IntPos[2], BITS_FOR_POSITIONS);
#endif
#ifdef SUBFIND
      /* make sure that for particles not in group we have no binding energy and no subrankgr, so that they are ordered only by the key
       */
      if(Tp->PS[i].GroupNr.get() == HALONR_MAX)
        {
          Tp->PS[i].SubRankInGr        = 0;
          Tp->PS[i].v.DM_BindingEnergy = 0;
        }
#endif
    }

#ifndef LEAN
  mycxxsort(Tp->PS, Tp->PS + Tp->NumPart, fof_compare_subfind_data_Type);
#endif

  for(int i = 0, off = 0; i < NTYPES; i++)
    {
      mycxxsort_parallel(Tp->PS + off, Tp->PS + off + ntype[i], fof_compare_subfind_data_GroupNr_SubNr_Egy_Key, Communicator);

      off += ntype[i];
    }

  for(int i = 0; i < Tp->NumPart; i++)
    {
      Tp->PS[i].TargetTask  = ThisTask;
      Tp->PS[i].TargetIndex = i;
    }

  /* now bring back into starting order */
  mycxxsort_parallel(Tp->PS, Tp->PS + Tp->NumPart, fof_compare_subfind_data_OriginTask_OriginIndex, Communicator);

  /* finally, reorder both P[] and PS[] */
  FoFDomain->particle_exchange_based_on_PS(Communicator);
}

/* calculate linkling length based on mean particle separation */
template <typename partset>
double fof<partset>::fof_get_comoving_linking_length(void)
{
  int ndm = 0;
  long long ndmtot;
  double mass = 0, masstot;

  for(int i = 0; i < Tp->NumPart; i++)
    if(is_type_primary_link_type(Tp->P[i].getType()))
      {
        ndm++;
        mass += Tp->P[i].getMass();
      }
  sumup_large_ints(1, &ndm, &ndmtot, Communicator);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, Communicator);

  double rhodm;
  if(Tp->TotNumGas > 0)
    rhodm = (All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  else
    rhodm = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  return FOF_LINKLENGTH * pow(masstot / ndmtot / rhodm, 1.0 / 3);
}

template <typename partset>
void fof<partset>::fof_compile_catalogue(double inner_distance)
{
  /* sort according to MinID, this brings particles belonging to the same group together */
  mycxxsort(FOF_PList, FOF_PList + Tp->NumPart, fof_compare_FOF_PList_MinID);

  /* we now use the auxiliary FOF_GList structure to determine the group lengths.
   * LocCount will count the length of the group piece that is on the principal processor
   * of that group, while other group pieces on other processors are counted through ExtCount
   */
  for(int i = 0; i < Tp->NumPart; i++)
    {
      FOF_GList[i].MinID     = FOF_PList[i].MinID;
      FOF_GList[i].MinIDTask = FOF_PList[i].MinIDTask;
      FOF_GList[i].Count     = 1;
#if defined(LIGHTCONE_PARTICLES_GROUPS)
      FOF_GList[i].DistanceOrigin = FOF_PList[i].DistanceOrigin;
#endif
    }

  /* Now we are going to eliminate duplicates in FOF_GList with respect to MinID, i.e. each group
   * piece on the local processor will be squashed to one entry. Some groups may be present as
   * pieces on several different processors.
   */
  if(Tp->NumPart)
    NgroupsExt = 1;
  else
    NgroupsExt = 0;

  for(int i = 1, start = 0; i < Tp->NumPart; i++)
    {
      if(FOF_GList[i].MinID.get() == FOF_GList[start].MinID.get())
        {
          if(FOF_GList[i].MinIDTask != FOF_GList[start].MinIDTask)
            Terminate("FOF_GList[i].MinIDTask=%d != FOF_GList[start].MinIDTask=%d", FOF_GList[i].MinIDTask,
                      FOF_GList[start].MinIDTask);

          FOF_GList[start].Count += FOF_GList[i].Count;

#if defined(LIGHTCONE_PARTICLES_GROUPS)
          if(FOF_GList[start].DistanceOrigin != FOF_GList[i].DistanceOrigin)
            Terminate(
                "start=%d i=%d  Tp->NumPart=%d  FOF_GList[start].DistanceOrigin=%g != FOF_GList[i].DistanceOrigin=%g  MinID=%lld  "
                "MinIDTask=%d\n",
                start, i, Tp->NumPart, FOF_GList[start].DistanceOrigin, FOF_GList[i].DistanceOrigin,
                (long long)FOF_GList[i].MinID.get(), FOF_GList[i].MinIDTask);
#endif
        }
      else
        {
          start            = NgroupsExt;
          FOF_GList[start] = FOF_GList[i];
          NgroupsExt++;
        }
    }

  /* we resize FOF_GList, which has shrunk */
  FOF_GList = (fof_group_list *)Mem.myrealloc_movable(FOF_GList, sizeof(fof_group_list) * NgroupsExt);

  /* sort the group pieces according to task */
  mycxxsort(FOF_GList, FOF_GList + NgroupsExt, fof_compare_FOF_GList_MinIDTask);

  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* count how many group pieces we have for each task */
  for(int i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(int i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  /* inform everybody about how much they have to receive */
  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

  /* count how many we get and prepare offset tables */
  int nimport    = 0;
  Recv_offset[0] = 0, Send_offset[0] = 0;
  for(int j = 0; j < NTask; j++)
    {
      if(j == ThisTask) /* we will not exchange the ones that are local */
        Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* allocate some temporary storage for foreign group pieces */
  fof_group_list *get_FOF_GList = (fof_group_list *)Mem.mymalloc("get_FOF_GList", nimport * sizeof(fof_group_list));

  /* get them */
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              myMPI_Sendrecv(&FOF_GList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(fof_group_list), MPI_BYTE, recvTask,
                             TAG_DENS_A, &get_FOF_GList[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(fof_group_list),
                             MPI_BYTE, recvTask, TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);
            }
        }
    }

  /* for the incoming pieces, we re-purpose MinIDTask and set it in ascending order in order to use this later to reestablish the order
   * we had */
  for(int i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = i;

  /* sort both the local group pieces and the incoming group pieces so that we can match them efficiently */
  mycxxsort(FOF_GList, FOF_GList + NgroupsExt, fof_compare_FOF_GList_MinID);
  mycxxsort(get_FOF_GList, get_FOF_GList + nimport, fof_compare_FOF_GList_MinID);

  /* Merge the imported group pieces with the local group pieces.
   * For all local group pieces in FOF_GList, the total number of particles on other processors will
   * be established in ExtCount.
   */
  int start = 0;
  for(int i = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID.get() < get_FOF_GList[i].MinID.get())
        {
          start++;
          if(start >= NgroupsExt)
            Terminate("start=%d >= NgroupsExt=%d", start, NgroupsExt);
        }

      if(FOF_GList[start].MinIDTask != ThisTask)
        Terminate("FOF_GList[start].MinIDTask=%d != ThisTask=%d", FOF_GList[start].MinIDTask, ThisTask);

      if(FOF_GList[start].MinID.get() != get_FOF_GList[i].MinID.get())
        Terminate(
            "FOF_GList[start].MinID != get_FOF_GList[i].MinID start=%d i=%d FOF_GList[start].MinID=%llu get_FOF_GList[i].MinID=%llu\n",
            start, i, (long long)FOF_GList[start].MinID.get(), (long long)get_FOF_GList[i].MinID.get());

      FOF_GList[start].Count += get_FOF_GList[i].Count;
    }

  /* copy the size information back into the received list, to inform the group pieces on originating processors */
  start = 0;
  for(int i = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID.get() < get_FOF_GList[i].MinID.get())
        {
          start++;
          if(start >= NgroupsExt)
            Terminate("start >= NgroupsExt");
        }

      get_FOF_GList[i].Count = FOF_GList[start].Count;
    }

  /* Sort the imported list according to MinIDTask. This reestablishes the order we had before the previous.
   * We also sort the local group pieces according to MinIDTask, so that we can fill in the exported info
   * at the right place again.
   */
  mycxxsort(get_FOF_GList, get_FOF_GList + nimport, fof_compare_FOF_GList_MinIDTask);
  mycxxsort(FOF_GList, FOF_GList + NgroupsExt, fof_compare_FOF_GList_MinIDTask);

  /* fix the value of MinIDTask again that we had temporarily overwritten */
  for(int i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = ThisTask;

  /* bring the data back to the originating processors */
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              myMPI_Sendrecv(&get_FOF_GList[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(fof_group_list), MPI_BYTE, recvTask,
                             TAG_DENS_A, &FOF_GList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(fof_group_list), MPI_BYTE,
                             recvTask, TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);
            }
        }
    }

  /* free our temporary list */
  Mem.myfree(get_FOF_GList);

  /* Now we determine how many groups we have above the group size limit.
   * Groups that are too small, or that are not guaranteed to be in the current lightcone segment
   * are eliminated. Local groups are counted in length to get the total size of particles in stored groups.
   */
#if defined(LIGHTCONE_PARTICLES_GROUPS)
  double save_distance = 0;
  if(inner_distance > 0)
    save_distance = inner_distance + Tp->LinkL;
#endif

  Ngroups = 0, Nids = 0;
  for(int i = 0; i < NgroupsExt; i++)
    {
      if(FOF_GList[i].Count < FOF_GROUP_MIN_LEN)
        {
          FOF_GList[i] = FOF_GList[NgroupsExt - 1];
          NgroupsExt--;
          i--;
        }
#if defined(LIGHTCONE_PARTICLES_GROUPS)
      else if(FOF_GList[i].DistanceOrigin < save_distance)
        {
          FOF_GList[i] = FOF_GList[NgroupsExt - 1];
          NgroupsExt--;
          i--;
        }
#endif
      else
        {
          if(FOF_GList[i].MinIDTask == ThisTask)
            {
              Ngroups++;
              Nids += FOF_GList[i].Count;
            }
        }
    }

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  /* we resize FOF_GList again, which has shrunk */
  FOF_GList = (fof_group_list *)Mem.myrealloc_movable(FOF_GList, sizeof(fof_group_list) * NgroupsExt);
}

template <typename partset>
void fof<partset>::fof_assign_group_numbers(void)
{
  mpi_printf("FOF: start assigning group numbers\n");

  double t0 = Logs.second();

  for(int i = 0; i < NgroupsExt; i++)
    Group[i].OriginTask = ThisTask;

  // carry out a parallel sort that brings the group segments into order of the total group length, with the primary segment coming
  // first
  mycxxsort_parallel(Group, Group + NgroupsExt, fof_compare_Group_Len_MinID_DiffOriginTaskMinIDTask, Communicator);

  /* assign group numbers to all local groups */
  int ngr = 0;
  for(int i = 0; i < NgroupsExt; i++)
    {
      if(Group[i].OriginTask == Group[i].MinIDTask)
        ngr++;

      Group[i].GroupNr = ngr - 1;
    }

  /* now need to count how many groups there are on earlier CPUs, so that we can
   * increase the assigned group numbers accordingly
   */
  int *ngr_list = (int *)Mem.mymalloc("ngr_list", sizeof(int) * NTask);
  MPI_Allgather(&ngr, 1, MPI_INT, ngr_list, 1, MPI_INT, Communicator);

  long long ngr_sum = 0;
  for(int j = 0; j < ThisTask; j++)
    ngr_sum += ngr_list[j];

  Mem.myfree(ngr_list);

  /* increase the group numbers with the cumulative count on earlier processors */
  for(int i = 0; i < NgroupsExt; i++)
    Group[i].GroupNr += ngr_sum;

  /* check that have consistent group numbers */
  sumup_large_ints(1, &ngr, &ngr_sum, Communicator);
  if(ngr_sum != TotNgroups)
    Terminate("inconsistency  ngr_sum=%lld\n", ngr_sum);

  /* bring the group list back into the original order */
  mycxxsort_parallel(Group, Group + NgroupsExt, fof_compare_Group_OriginTask_MinID, Communicator);

  /* Let's now mark all particles that are not in any group by assigning them a fiducial maximum (unused)
   * group number
   */
  for(int i = 0; i < Tp->NumPart; i++)
    Tp->PS[i].GroupNr.set(HALONR_MAX);

  long long Nids_old = Nids;

  /* we now aim to assign the group number also to the particles that make up groups,
   * in the auxiliary PS array
   */
  Nids = 0;
  for(int i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID.get() < Group[i].MinID)
        {
          start++;
          if(start > Tp->NumPart)
            Terminate("start > Tp->NumPart");
        }

      if(FOF_PList[start].MinID.get() != Group[i].MinID)
        Terminate("FOF_PList[start=%d].MinID=%lld != Group[i=%d].MinID=%lld", start, (long long)FOF_PList[start].MinID.get(), i,
                  (long long)Group[i].MinID);

      int lenloc;
      for(lenloc = 0; start + lenloc < Tp->NumPart;)
        if(FOF_PList[start + lenloc].MinID.get() == Group[i].MinID)
          {
            Tp->PS[FOF_PList[start + lenloc].Pindex].GroupNr.set(Group[i].GroupNr);
            Nids++;
            lenloc++;
          }
        else
          break;

      start += lenloc;
    }

  long long totNids;
  sumup_longs(1, &Nids, &totNids, Communicator);

  if(totNids != TotNids)
    Terminate("Task=%d Nids=%lld  Nids_old=%lld  totNids=%lld TotNids=%lld\n", ThisTask, Nids, Nids_old, totNids, TotNids);

  double t1 = Logs.second();

  mpi_printf("FOF: Assigning of group numbers took = %g sec\n", Logs.timediff(t0, t1));
}

/* This function computes the Group[].OffsetType variable, which gives for each type how many
 * particles are before it of the same time in earlier groups in the array.
 */
template <typename partset>
void fof<partset>::fof_assign_group_offset(void)
{
  /* Tell everybody, how many particles the groups stored by each processor contain */

  int gtype_loc[NTYPES]; /* particles of each type associated with locally stored groups */
  long long gtype_previous[NTYPES];

  for(int i = 0; i < NTYPES; i++)
    gtype_loc[i] = 0;

  for(int i = 0; i < Ngroups; i++)
    for(int j = 0; j < NTYPES; j++)
      gtype_loc[j] += Group[i].LenType[j];

  int *gtype_all = (int *)Mem.mymalloc("gtype_all", NTYPES * NTask * sizeof(int));
  MPI_Allgather(gtype_loc, NTYPES, MPI_INT, gtype_all, NTYPES, MPI_INT, Communicator);

  for(int i = 0; i < NTYPES; i++)
    gtype_previous[i] = 0;

  for(int i = 0; i < ThisTask; i++)
    for(int j = 0; j < NTYPES; j++)
      gtype_previous[j] += gtype_all[i * NTYPES + j];

  if(Ngroups > 0)
    for(int j = 0; j < NTYPES; j++)
      Group[0].OffsetType[j] = gtype_previous[j];

  for(int i = 1; i < Ngroups; i++)
    for(int j = 0; j < NTYPES; j++)
      Group[i].OffsetType[j] = Group[i - 1].OffsetType[j] + Group[i - 1].LenType[j];

  Mem.myfree(gtype_all);
}

template <typename partset>
void fof<partset>::fof_compute_group_properties(int gr, int start, int len)
{
  int start_index = FOF_PList[start].Pindex;

  Group[gr].Mass   = 0;
  Group[gr].Ascale = 0;
#ifdef STARFORMATION
  Group[gr].Sfr = 0;
#endif
#if defined(SUBFIND_ORPHAN_TREATMENT)
  Group[gr].LenPrevMostBnd = 0;
#endif

  for(int k = 0; k < 3; k++)
    {
      Group[gr].CM[k]          = 0;
      Group[gr].Vel[k]         = 0;
      Group[gr].FirstIntPos[k] = Tp->P[start_index].IntPos[k];
    }

  for(int k = 0; k < NTYPES; k++)
    {
      Group[gr].LenType[k]  = 0;
      Group[gr].MassType[k] = 0;
    }

  for(int k = 0; k < len; k++)
    {
      int index = FOF_PList[start + k].Pindex;

      Group[gr].Mass += Tp->P[index].getMass();
      int type = Tp->P[index].getType();

      Group[gr].Ascale += Tp->P[index].getMass() * Tp->P[index].getAscale();

      Group[gr].LenType[type]++;
      Group[gr].MassType[type] += Tp->P[index].getMass();

#if defined(SUBFIND_ORPHAN_TREATMENT)
      if(Tp->P[index].ID.is_previously_most_bound())
        Group[gr].LenPrevMostBnd++;
#endif

#ifdef STARFORMATION
      if(Tp->P[index].getType() == 0)
        Group[gr].Sfr += Tp->SphP[index].Sfr;
#endif

      double xyz[3];
      Tp->nearest_image_intpos_to_pos(Tp->P[index].IntPos, Tp->P[start_index].IntPos,
                                      xyz); /* converts the integer distance to floating point */

      for(int j = 0; j < 3; j++)
        {
          Group[gr].CM[j] += Tp->P[index].getMass() * xyz[j];
          Group[gr].Vel[j] += Tp->P[index].getMass() * Tp->P[index].Vel[j];
        }
    }
}

template <typename partset>
void fof<partset>::fof_add_in_properties_of_group_segments(void)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * NTask);

  /* sort the groups according to task */
  mycxxsort(Group, Group + NgroupsExt, fof_compare_Group_MinIDTask);

  /* count how many we have of each task */
  for(int i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(int i = 0; i < NgroupsExt; i++)
    Send_count[Group[i].MinIDTask]++;

  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

  int nimport    = 0;
  Recv_offset[0] = 0, Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      if(j == ThisTask) /* we will not exchange the ones that are local */
        Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  group_properties *get_Group = (group_properties *)Mem.mymalloc("get_Group", sizeof(group_properties) * nimport);

  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group data */
              myMPI_Sendrecv(&Group[Send_offset[recvTask]], Send_count[recvTask] * sizeof(group_properties), MPI_BYTE, recvTask,
                             TAG_DENS_A, &get_Group[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(group_properties), MPI_BYTE,
                             recvTask, TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);
            }
        }
    }

  /* sort the groups again according to MinID */
  mycxxsort(Group, Group + NgroupsExt, fof_compare_Group_MinID);
  mycxxsort(get_Group, get_Group + nimport, fof_compare_Group_MinID);

  int start = 0;
  /* now add in the partial imported group data to the main ones */
  for(int i = 0; i < nimport; i++)
    {
      while(Group[start].MinID < get_Group[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            Terminate("start >= NgroupsExt");
        }

      Group[start].Mass += get_Group[i].Mass;
      Group[start].Ascale += get_Group[i].Ascale;

      for(int j = 0; j < NTYPES; j++)
        {
          Group[start].LenType[j] += get_Group[i].LenType[j];
          Group[start].MassType[j] += get_Group[i].MassType[j];
        }
#if defined(SUBFIND_ORPHAN_TREATMENT)
      Group[start].LenPrevMostBnd += get_Group[i].LenPrevMostBnd;
#endif
#ifdef STARFORMATION
      Group[start].Sfr += get_Group[i].Sfr;
#endif

      double tmpxyz[3] = {get_Group[i].CM[0] / get_Group[i].Mass, get_Group[i].CM[1] / get_Group[i].Mass,
                          get_Group[i].CM[2] / get_Group[i].Mass};

      MyIntPosType delta[3];
      Tp->pos_to_signedintpos(tmpxyz, (MySignedIntPosType *)delta);

      delta[0] += get_Group[i].FirstIntPos[0];
      delta[1] += get_Group[i].FirstIntPos[1];
      delta[2] += get_Group[i].FirstIntPos[2];

      Tp->constrain_intpos(delta); /* will only do something if we have a stretched box */

      double xyz[3];
      Tp->nearest_image_intpos_to_pos(delta, Group[start].FirstIntPos, xyz); /* converts the integer distance to floating point */

      for(int j = 0; j < 3; j++)
        {
          Group[start].CM[j] += get_Group[i].Mass * xyz[j];
          Group[start].Vel[j] += get_Group[i].Vel[j];
        }
    }

  Mem.myfree(get_Group);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

template <typename partset>
void fof<partset>::fof_finish_group_properties(void)
{
  for(int i = 0; i < NgroupsExt; i++)
    {
      if(Group[i].MinIDTask == ThisTask)
        {
          double cm[3];
          for(int j = 0; j < 3; j++)
            {
              Group[i].Vel[j] /= Group[i].Mass;
              cm[j] = Group[i].CM[j] / Group[i].Mass;
            }

          Group[i].Ascale /= Group[i].Mass;

          MyIntPosType delta[3];
          Tp->pos_to_signedintpos(cm, (MySignedIntPosType *)delta);

          delta[0] += Group[i].FirstIntPos[0];
          delta[1] += Group[i].FirstIntPos[1];
          delta[2] += Group[i].FirstIntPos[2];

          Tp->constrain_intpos(delta); /* will only do something if we have a stretched box */

          fof_get_halo_position(delta, cm);

          Group[i].CM[0] = cm[0];
          Group[i].CM[1] = cm[1];
          Group[i].CM[2] = cm[2];

          /* define group position as CM . This will be overwritten in case Subfind is used with
           * the position of the potential minimum
           */
          for(int j = 0; j < 3; j++)
            Group[i].Pos[j] = Group[i].CM[j];

          for(int j = 0; j < 3; j++)
            Group[i].IntPos[j] = delta[j];
        }
    }

  int ngr = NgroupsExt;

  /* eliminate the non-local groups */
  for(int i = 0; i < ngr; i++)
    {
      if(Group[i].MinIDTask != ThisTask)
        {
          Group[i] = Group[ngr - 1];
          i--;
          ngr--;
        }
    }

  if(ngr != Ngroups)
    Terminate("ngr != Ngroups");

  mycxxsort(Group, Group + Ngroups, fof_compare_Group_MinID);
}

#if defined(LIGHTCONE_PARTICLES_GROUPS)

template <>
double fof<simparticles>::fof_distance_to_origin(int i)
{
  return 0;  // dummy information for ordinary timeslices
}

template <>
double fof<lcparticles>::fof_distance_to_origin(int i)
{
  return Tp->signedintpos_to_distanceorigin((MySignedIntPosType *)Tp->P[i].IntPos);
}

#endif

template <>
void fof<simparticles>::fof_get_halo_position(MyIntPosType *intpos, double *pos)
{
  Tp->intpos_to_pos(intpos, pos);
}

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
template <>
void fof<lcparticles>::fof_get_halo_position(MyIntPosType *intpos, double *pos)
{
  MyIntPosType origin[3] = {0, 0, 0};

  Tp->nearest_image_intpos_to_pos(intpos, origin, pos);
}
#endif

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif /* of FOF */
