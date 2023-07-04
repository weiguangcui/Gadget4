/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  parallel_sort.h
 *
 *  \brief parallel sort routine that leaves the number of elements per processor invariant
 */

#ifndef PARALLEL_SORT_H
#define PARALLEL_SORT_H

#include "gadgetconfig.h"

#include "../data/mymalloc.h"
#include "../sort/cxxsort.h"

//#define CHECK_LOCAL_RANK

template <typename It, typename Comp>
class IdxComp__
{
 private:
  It begin;
  Comp comp;

 public:
  IdxComp__(It begin_, Comp comp_) : begin(begin_), comp(comp_) {}
  bool operator()(std::size_t a, std::size_t b) const { return comp(*(begin + a), *(begin + b)); }
};

/*! Performs an indirect sort on the supplied iterator range and returns in
    \a idx a \a vector containing the indices of the smallest, second smallest,
    third smallest, etc. element, according to \a comp. */
template <typename It, typename T2, typename Comp>
inline void buildIndex(It begin, It end, T2 *idx, Comp comp)
{
  using namespace std;
  T2 num = end - begin;
  for(T2 i = 0; i < num; ++i)
    idx[i] = i;
  mycxxsort(idx, idx + num, IdxComp__<It, Comp>(begin, comp));
}

template <typename T, typename Comp>
void get_local_rank(const T &element, std::size_t tie_breaking_rank, const T *base, size_t nmemb, size_t noffs_thistask,
                    long long left, long long right, size_t *loc, Comp comp)
{
  if(right < left)
    Terminate("right < left");

  if(left == 0 && right == (int)nmemb + 1)
    {
      if(comp(base[nmemb - 1], element))
        {
          *loc = nmemb;
          return;
        }
      else if(comp(element, base[0]))
        {
          *loc = 0;
          return;
        }
    }

  if(right == left) /* looks like we already converged to the proper rank */
    {
      *loc = left;
    }
  else
    {
      if(comp(base[right - 1], element)) /* the last element is smaller, hence all elements are on the left */
        *loc = (right - 1) + 1;
      else if(comp(element, base[left])) /* the first element is already larger, hence no element is on the left */
        *loc = left;
      else
        {
          while(right > left)
            {
              long long mid = ((right - 1) + left) / 2;

              int cmp = comp(base[mid], element) ? -1 : (comp(element, base[mid]) ? +1 : 0);
              if(cmp == 0)
                {
                  if(mid + noffs_thistask < tie_breaking_rank)
                    cmp = -1;
                  else if(mid + noffs_thistask > tie_breaking_rank)
                    cmp = +1;
                }

              if(cmp == 0) /* element has exactly been found */
                {
                  *loc = mid;
                  break;
                }

              if((right - 1) == left) /* elements is not on this CPU */
                {
                  if(cmp < 0)
                    *loc = mid + 1;
                  else
                    *loc = mid;
                  break;
                }

              if(cmp < 0)
                {
                  left = mid + 1;
                }
              else
                {
                  if((right - 1) == left + 1)
                    {
                      if(mid != left)
                        Terminate("Can't be: -->left=%lld  right=%lld\n", left, right);

                      *loc = left;
                      break;
                    }

                  right = mid;
                }
            }
        }
    }
}

#ifdef CHECK_LOCAL_RANK
template <typename T, typename Comp>
inline void check_local_rank(const T &element,                /* element of which we want the rank */
                             size_t tie_breaking_rank,        /* the initial global rank of this element (needed for breaking ties) */
                             const T *base,                   /* base address of local data */
                             size_t nmemb,                    /* number and size of local data */
                             size_t noffs_thistask,           /* cumulative length of data on lower tasks */
                             long long left, long long right, /* range of elements on local task that may hold the element */
                             size_t loc, Comp comp)           /* user-specified  comparison function */
{
  long long count = 0;

  for(size_t i = 0; i < nmemb; i++)
    {
      int cmp = comp(base[i], element) ? -1 : (comp(element, base[i]) ? +1 : 0);

      if(cmp == 0)
        {
          if(noffs_thistask + i < tie_breaking_rank)
            cmp = -1;
        }

      if(cmp < 0)
        count++;
    }

  if(count != (long long)loc)
    Terminate("Inconsistency: loc=%lld count=%lld  left=%lld right=%lld  nmemb=%lld\n", (long long)loc, count, left, right,
              (long long)nmemb);
}
#endif

template <typename T, typename Comp>
inline double mycxxsort_parallel(T *begin, T *end, Comp comp, MPI_Comm comm)
{
  const int MAX_ITER_PARALLEL_SORT = 500;
  int ranks_not_found, Local_ThisTask, Local_NTask, Color, new_max_loc;
  size_t tie_breaking_rank, new_tie_breaking_rank, rank;
  MPI_Comm MPI_CommLocal;

  double ta    = Logs.second();
  size_t nmemb = end - begin;
  size_t size  = sizeof(T);
  /* do a serial sort of the local data up front */
  mycxxsort(begin, end, comp);

  /* we create a communicator that contains just those tasks with nmemb > 0. This makes
   *  it easier to deal with CPUs that do not hold any data.
   */
  if(nmemb)
    Color = 1;
  else
    Color = 0;

  int thistask;
  MPI_Comm_rank(comm, &thistask);

  MPI_Comm_split(comm, Color, thistask, &MPI_CommLocal);
  MPI_Comm_rank(MPI_CommLocal, &Local_ThisTask);
  MPI_Comm_size(MPI_CommLocal, &Local_NTask);

  if(Local_NTask > 1 && Color == 1)
    {
      size_t *nlist = (size_t *)Mem.mymalloc("nlist", Local_NTask * sizeof(size_t));
      size_t *noffs = (size_t *)Mem.mymalloc("noffs", Local_NTask * sizeof(size_t));

      MPI_Allgather(&nmemb, sizeof(size_t), MPI_BYTE, nlist, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

      noffs[0] = 0;
      for(int i = 1; i < Local_NTask; i++)
        noffs[i] = noffs[i - 1] + nlist[i - 1];

      T *element_guess                  = (T *)Mem.mymalloc("element_guess", Local_NTask * size);
      size_t *element_tie_breaking_rank = (size_t *)Mem.mymalloc("element_tie_breaking_rank", Local_NTask * sizeof(size_t));
      size_t *desired_glob_rank         = (size_t *)Mem.mymalloc("desired_glob_rank", Local_NTask * sizeof(size_t));
      size_t *current_glob_rank         = (size_t *)Mem.mymalloc("current_glob_rank", Local_NTask * sizeof(size_t));
      size_t *current_loc_rank          = (size_t *)Mem.mymalloc("current_loc_rank", Local_NTask * sizeof(size_t));
      long long *range_left             = (long long *)Mem.mymalloc("range_left", Local_NTask * sizeof(long long));
      long long *range_right            = (long long *)Mem.mymalloc("range_right", Local_NTask * sizeof(long long));
      int *max_loc                      = (int *)Mem.mymalloc("max_loc", Local_NTask * sizeof(int));

      size_t *list           = (size_t *)Mem.mymalloc("list", Local_NTask * sizeof(size_t));
      size_t *range_len_list = (size_t *)Mem.mymalloc("range_len_list", Local_NTask * sizeof(long long));
      T median_element;
      T *median_element_list                = (T *)Mem.mymalloc("median_element_list", Local_NTask * size);
      size_t *tie_breaking_rank_list        = (size_t *)Mem.mymalloc("tie_breaking_rank_list", Local_NTask * sizeof(size_t));
      int *index_list                       = (int *)Mem.mymalloc("index_list", Local_NTask * sizeof(int));
      int *max_loc_list                     = (int *)Mem.mymalloc("max_loc_list", Local_NTask * sizeof(int));
      size_t *source_range_len_list         = (size_t *)Mem.mymalloc("source_range_len_list", Local_NTask * sizeof(long long));
      size_t *source_tie_breaking_rank_list = (size_t *)Mem.mymalloc("source_tie_breaking_rank_list", Local_NTask * sizeof(long long));
      T *source_median_element_list         = (T *)Mem.mymalloc("source_median_element_list", Local_NTask * size);
      T new_element_guess;

      for(int i = 0; i < Local_NTask - 1; i++)
        {
          desired_glob_rank[i] = noffs[i + 1];
          current_glob_rank[i] = 0;
          range_left[i]        = 0;     /* first element that it can be */
          range_right[i]       = nmemb; /* first element that it can not be */
        }

      /* now we determine the first split element guess, which is the same for all divisions in the first iteration */

      /* find the median of each processor, and then take the median among those values.
       * This should work reasonably well even for extremely skewed distributions
       */
      long long range_len = range_right[0] - range_left[0];

      if(range_len >= 1)
        {
          long long mid     = (range_left[0] + range_right[0]) / 2;
          median_element    = begin[mid];
          tie_breaking_rank = mid + noffs[Local_ThisTask];
        }

      MPI_Gather(&range_len, sizeof(long long), MPI_BYTE, range_len_list, sizeof(long long), MPI_BYTE, 0, MPI_CommLocal);
      MPI_Gather(&median_element, size, MPI_BYTE, median_element_list, size, MPI_BYTE, 0, MPI_CommLocal);
      MPI_Gather(&tie_breaking_rank, sizeof(size_t), MPI_BYTE, tie_breaking_rank_list, sizeof(size_t), MPI_BYTE, 0, MPI_CommLocal);

      if(Local_ThisTask == 0)
        {
          for(int j = 0; j < Local_NTask; j++)
            max_loc_list[j] = j;

          /* eliminate the elements that are undefined because the corresponding CPU has zero range left */
          int nleft = Local_NTask;

          for(int j = 0; j < nleft; j++)
            {
              if(range_len_list[j] < 1)
                {
                  range_len_list[j] = range_len_list[nleft - 1];
                  if(range_len_list[nleft - 1] >= 1 && j != (nleft - 1))
                    {
                      median_element_list[j]    = median_element_list[nleft - 1];
                      tie_breaking_rank_list[j] = tie_breaking_rank_list[nleft - 1];
                      max_loc_list[j]           = max_loc_list[nleft - 1];
                    }

                  nleft--;
                  j--;
                }
            }

          /* do a serial sort of the remaining elements (indirectly, so that we have the order of tie breaking list as well) */
          buildIndex(median_element_list, median_element_list + nleft, index_list, comp);

          /* now select the median of the medians */
          int mid                      = nleft / 2;
          element_guess[0]             = median_element_list[index_list[mid]];
          element_tie_breaking_rank[0] = tie_breaking_rank_list[index_list[mid]];
          max_loc[0]                   = max_loc_list[index_list[mid]];
        }

      MPI_Bcast(element_guess, size, MPI_BYTE, 0, MPI_CommLocal);
      MPI_Bcast(&element_tie_breaking_rank[0], sizeof(size_t), MPI_BYTE, 0, MPI_CommLocal);
      MPI_Bcast(&max_loc[0], 1, MPI_INT, 0, MPI_CommLocal);

      for(int i = 1; i < Local_NTask - 1; i++)
        {
          element_guess[i]             = element_guess[0];
          element_tie_breaking_rank[i] = element_tie_breaking_rank[0];
          max_loc[i]                   = max_loc[0];
        }

      int iter = 0;

      do
        {
          for(int i = 0; i < Local_NTask - 1; i++)
            {
              if(current_glob_rank[i] != desired_glob_rank[i])
                {
                  get_local_rank(element_guess[i], element_tie_breaking_rank[i], begin, nmemb, noffs[Local_ThisTask], range_left[i],
                                 range_right[i], &current_loc_rank[i], comp);

#ifdef CHECK_LOCAL_RANK
                  check_local_rank(element_guess[i], element_tie_breaking_rank[i], begin, nmemb, noffs[Local_ThisTask], range_left[i],
                                   range_right[i], current_loc_rank[i], comp);
#endif
                }
            }

          /* now compute the global ranks by summing the local ranks */
          /* Note: the last element in current_loc_rank is not defined. It will be summed by the last processor, and stored in the last
           * element of current_glob_rank */
          myMPI_Alltoall(current_loc_rank, sizeof(size_t), MPI_BYTE, list, sizeof(size_t), MPI_BYTE, MPI_CommLocal);
          rank = 0;
          for(int j = 0; j < Local_NTask; j++)
            rank += list[j];
          MPI_Allgather(&rank, sizeof(size_t), MPI_BYTE, current_glob_rank, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

          ranks_not_found = 0;
          for(int i = 0; i < Local_NTask - 1; i++)
            {
              if(current_glob_rank[i] != desired_glob_rank[i]) /* here we're not yet done */
                {
                  ranks_not_found++;

                  if(current_glob_rank[i] < desired_glob_rank[i])
                    {
                      range_left[i] = current_loc_rank[i];

                      if(Local_ThisTask == max_loc[i])
                        range_left[i]++;
                    }

                  if(current_glob_rank[i] > desired_glob_rank[i])
                    range_right[i] = current_loc_rank[i];
                }
            }

          /* now we need to determine new element guesses */
          for(int i = 0; i < Local_NTask - 1; i++)
            {
              if(current_glob_rank[i] != desired_glob_rank[i]) /* here we're not yet done */
                {
                  /* find the median of each processor, and then take the median among those values.
                   * This should work reasonably well even for extremely skewed distributions
                   */
                  source_range_len_list[i] = range_right[i] - range_left[i];

                  if(source_range_len_list[i] >= 1)
                    {
                      long long middle                 = (range_left[i] + range_right[i]) / 2;
                      source_median_element_list[i]    = begin[middle];
                      source_tie_breaking_rank_list[i] = middle + noffs[Local_ThisTask];
                    }
                }
            }

          myMPI_Alltoall(source_range_len_list, sizeof(long long), MPI_BYTE, range_len_list, sizeof(long long), MPI_BYTE,
                         MPI_CommLocal);
          myMPI_Alltoall(source_median_element_list, size, MPI_BYTE, median_element_list, size, MPI_BYTE, MPI_CommLocal);
          myMPI_Alltoall(source_tie_breaking_rank_list, sizeof(size_t), MPI_BYTE, tie_breaking_rank_list, sizeof(size_t), MPI_BYTE,
                         MPI_CommLocal);

          if(Local_ThisTask < Local_NTask - 1)
            {
              if(current_glob_rank[Local_ThisTask] !=
                 desired_glob_rank[Local_ThisTask]) /* in this case we're not yet done for this split point */
                {
                  for(int j = 0; j < Local_NTask; j++)
                    max_loc_list[j] = j;

                  /* eliminate the elements that are undefined because the corresponding CPU has zero range left */
                  int nleft = Local_NTask;

                  for(int j = 0; j < nleft; j++)
                    {
                      if(range_len_list[j] < 1)
                        {
                          range_len_list[j] = range_len_list[nleft - 1];
                          if(range_len_list[nleft - 1] >= 1 && j != (nleft - 1))
                            {
                              median_element_list[j]    = median_element_list[nleft - 1];
                              tie_breaking_rank_list[j] = tie_breaking_rank_list[nleft - 1];
                              max_loc_list[j]           = max_loc_list[nleft - 1];
                            }

                          nleft--;
                          j--;
                        }
                    }

                  if((iter & 1))
                    {
                      size_t max_range = 0, maxj = 0;

                      for(int j = 0; j < nleft; j++)
                        if(range_len_list[j] > max_range)
                          {
                            max_range = range_len_list[j];
                            maxj      = j;
                          }

                      /* now select the median element from the task which has the largest range */
                      new_element_guess     = median_element_list[maxj];
                      new_tie_breaking_rank = tie_breaking_rank_list[maxj];
                      new_max_loc           = max_loc_list[maxj];
                    }
                  else
                    {
                      /* do a serial sort of the remaining elements (indirectly, so that we have the order of tie breaking list as
                       * well) */
                      buildIndex(median_element_list, median_element_list + nleft, index_list, comp);

                      /* now select the median of the medians */
                      int mid               = nleft / 2;
                      new_element_guess     = median_element_list[index_list[mid]];
                      new_tie_breaking_rank = tie_breaking_rank_list[index_list[mid]];
                      new_max_loc           = max_loc_list[index_list[mid]];
                    }
                }
              else
                {
                  /* in order to preserve existing guesses */
                  new_element_guess     = element_guess[Local_ThisTask];
                  new_tie_breaking_rank = element_tie_breaking_rank[Local_ThisTask];
                  new_max_loc           = max_loc[Local_ThisTask];
                }
            }

          MPI_Allgather(&new_element_guess, size, MPI_BYTE, element_guess, size, MPI_BYTE, MPI_CommLocal);
          MPI_Allgather(&new_tie_breaking_rank, sizeof(size_t), MPI_BYTE, element_tie_breaking_rank, sizeof(size_t), MPI_BYTE,
                        MPI_CommLocal);
          MPI_Allgather(&new_max_loc, 1, MPI_INT, max_loc, 1, MPI_INT, MPI_CommLocal);

          iter++;

          if(iter > (MAX_ITER_PARALLEL_SORT - 100) && Local_ThisTask == 0)
            {
              printf("PSORT: iter=%d: ranks_not_found=%d  Local_NTask=%d\n", iter, ranks_not_found, Local_NTask);
              myflush(stdout);
              if(iter > MAX_ITER_PARALLEL_SORT)
                Terminate("can't find the split points. That's odd");
            }
        }
      while(ranks_not_found);

      Mem.myfree(source_median_element_list);
      Mem.myfree(source_tie_breaking_rank_list);
      Mem.myfree(source_range_len_list);
      Mem.myfree(max_loc_list);
      Mem.myfree(index_list);
      Mem.myfree(tie_breaking_rank_list);
      Mem.myfree(median_element_list);

      /* At this point we have found all the elements corresponding to the desired split points */
      /* we can now go ahead and determine how many elements of the local CPU have to go to each other CPU */

      size_t *send_count  = (size_t *)Mem.mymalloc("send_count", Local_NTask * sizeof(size_t));
      size_t *recv_count  = (size_t *)Mem.mymalloc("recv_count", Local_NTask * sizeof(size_t));
      size_t *send_offset = (size_t *)Mem.mymalloc("send_offset", Local_NTask * sizeof(size_t));
      size_t *recv_offset = (size_t *)Mem.mymalloc("recv_offset", Local_NTask * sizeof(size_t));

      for(int i = 0; i < Local_NTask; i++)
        send_count[i] = 0;

      int target = 0;

      for(size_t i = 0; i < nmemb; i++)
        {
          while(target < Local_NTask - 1)
            {
              int cmp = comp(begin[i], element_guess[target]) ? -1 : (comp(element_guess[target], begin[i]) ? +1 : 0);
              if(cmp == 0)
                {
                  if(i + noffs[Local_ThisTask] < element_tie_breaking_rank[target])
                    cmp = -1;
                  else if(i + noffs[Local_ThisTask] > element_tie_breaking_rank[target])
                    cmp = +1;
                }
              if(cmp >= 0)
                target++;
              else
                break;
            }
          send_count[target]++;
        }

      myMPI_Alltoall(send_count, sizeof(size_t), MPI_BYTE, recv_count, sizeof(size_t), MPI_BYTE, MPI_CommLocal);

      size_t nimport = 0;

      recv_offset[0] = 0;
      send_offset[0] = 0;
      for(int j = 0; j < Local_NTask; j++)
        {
          nimport += recv_count[j];

          if(j > 0)
            {
              send_offset[j] = send_offset[j - 1] + send_count[j - 1];
              recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
            }
        }

      if(nimport != nmemb)
        Terminate("nimport=%lld != nmemb=%lld", (long long)nimport, (long long)nmemb);

      for(int j = 0; j < Local_NTask; j++)
        {
          send_count[j] *= size;
          recv_count[j] *= size;

          send_offset[j] *= size;
          recv_offset[j] *= size;
        }

      T *basetmp = (T *)Mem.mymalloc("basetmp", nmemb * size);

      /* exchange the data */
      myMPI_Alltoallv(begin, send_count, send_offset, basetmp, recv_count, recv_offset, sizeof(char), 1, MPI_CommLocal);

      memcpy(static_cast<void *>(begin), static_cast<void *>(basetmp), nmemb * size);
      Mem.myfree(basetmp);

      mycxxsort(begin, begin + nmemb, comp);

      Mem.myfree(recv_offset);
      Mem.myfree(send_offset);
      Mem.myfree(recv_count);
      Mem.myfree(send_count);

      Mem.myfree(range_len_list);
      Mem.myfree(list);
      Mem.myfree(max_loc);
      Mem.myfree(range_right);
      Mem.myfree(range_left);
      Mem.myfree(current_loc_rank);
      Mem.myfree(current_glob_rank);
      Mem.myfree(desired_glob_rank);
      Mem.myfree(element_tie_breaking_rank);
      Mem.myfree(element_guess);
      Mem.myfree(noffs);
      Mem.myfree(nlist);
    }

  MPI_Comm_free(&MPI_CommLocal);

  double tb = Logs.second();
  return Logs.timediff(ta, tb);
}

#endif
