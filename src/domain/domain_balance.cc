/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file domain_balance.cc
 *
 *  \brief contains routines to improve the domain balance by combining several patches per MPI rank
 */

#include "gadgetconfig.h"

#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/** This function uses the cumulative cost function (which weights work-load and memory-load equally) to subdivide
 *  the list of top-level leaf nodes into pieces that are (approximately) equal in size.
 */

#ifdef DOMAIN_SPECIAL_CHECK

template <typename partset>
void domain<partset>::domain_special_check(int mode, int ndomains)
{
  double *cost_data = (double *)Mem.mymalloc_clear("cost_data", sizeof(double) * 2 * ndomains * (NumTimeBinsToBeBalanced + 1));

  double *load         = cost_data;
  double *loadsph      = cost_data + ndomains;
  double *binGravCost  = cost_data + 2 * ndomains;
  double *binHydroCost = cost_data + 2 * ndomains + ndomains * NumTimeBinsToBeBalanced;

  for(int i = 0; i < Tp->NumPart; i++)
    {
      int no = n_to_no(i);

      int n;

      if(mode == 0)
        n = no;
      else
        n = TaskOfLeaf[no];

      if(n < 0 || n >= ndomains)
        Terminate("strange");

#ifdef SELFGRAVITY
      for(int k = 0; k < NumTimeBinsToBeBalanced; k++)
        {
          int bin = ListOfTimeBinsToBeBalanced[k];

          if(bin >= Tp->P[i].getTimeBinGrav())
            binGravCost[k * ndomains + n] += GravCostNormFactors[k] * Tp->P[i].getGravCost();
        }
#endif

      load[n] += NormFactorLoad;

      if(Tp->P[i].getType() == 0)
        {
          for(int k = 0; k < NumTimeBinsToBeBalanced; k++)
            {
              int bin = ListOfTimeBinsToBeBalanced[k];

              if(bin >= Tp->P[i].getTimeBinHydro())
                binHydroCost[k * ndomains + n] += HydroCostNormFactors[k];
            }

          loadsph[n] += NormFactorLoadSph;
        }
    }

  MPI_Allreduce(MPI_IN_PLACE, cost_data, 2 * ndomains * (NumTimeBinsToBeBalanced + 1), MPI_DOUBLE, MPI_SUM, Communicator);

  if(All.NumCurrentTiStep == 0 || All.NumCurrentTiStep == 2 || All.NumCurrentTiStep == 4)
    {
      if(ThisTask == 0)
        {
          char buf[MAXLEN_PATH_EXTRA];
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s/domain_data_%d_step%d.txt", All.OutputDir, mode, All.NumCurrentTiStep);
          FILE *fd = fopen(buf, "w");
          fprintf(fd, "%d %d\n", ndomains, NumTimeBinsToBeBalanced);
          for(int n = 0; n < ndomains; n++)
            {
              fprintf(fd, "%g  ", load[n]);
              for(int k = 0; k < NumTimeBinsToBeBalanced; k++)
                fprintf(fd, "%g  ", binGravCost[k * ndomains + n]);
              fprintf(fd, "\n");
            }
          fclose(fd);
        }
    }

  Mem.myfree(cost_data);
}

#endif

/** This function assigns the domain pieces to individual MPI tasks with the goal to balance the work-load
 *  on different timebins. The algorithm used works as follows:
 *
 *  The domains are assigned to the CPUs in sequence of decreasing "effective load", which is a simple combined measure of
 *  relative total gravity, hydro and memory load. For each assignment, a number of possible target CPUs are evaluated, and
 *  the assignment leading to the lowest total runtime is adopted.
 *  The set of target CPUs that is tested in each step is the one that
 *  consists of the CPUs that currently have the lowest load in the set of primary tasks that are examined.
 */
template <typename partset>
void domain<partset>::domain_combine_multipledomains(void)
{
  double t0 = Logs.second();

  /* we first determine the detailed cost of all the domain pieces (which are the top leaves in the tree), so that we can combine them
   * later on efficiently for different choices of nextra
   */

  double *cost_data = (double *)Mem.mymalloc_clear("cost_data", sizeof(double) * 2 * NTopleaves * (NumTimeBinsToBeBalanced + 1));

  double *load         = cost_data;
  double *loadsph      = cost_data + NTopleaves;
  double *binGravCost  = cost_data + 2 * NTopleaves;
  double *binHydroCost = cost_data + 2 * NTopleaves + NTopleaves * NumTimeBinsToBeBalanced;

  for(int i = 0; i < Tp->NumPart; i++)
    {
      int no = n_to_no(i);  // get the leave node

#ifdef SELFGRAVITY
      for(int k = 0; k < NumTimeBinsToBeBalanced; k++)
        {
          int bin = ListOfTimeBinsToBeBalanced[k];

          if(bin >= Tp->P[i].getTimeBinGrav())
            binGravCost[k * NTopleaves + no] += GravCostNormFactors[k] * Tp->P[i].getGravCost();
        }
#endif

      load[no] += NormFactorLoad;

      if(Tp->P[i].getType() == 0)
        {
          for(int k = 0; k < NumTimeBinsToBeBalanced; k++)
            {
              int bin = ListOfTimeBinsToBeBalanced[k];

              if(bin >= Tp->P[i].getTimeBinHydro())
                binHydroCost[k * NTopleaves + no] += HydroCostNormFactors[k];
            }

          loadsph[no] += NormFactorLoadSph;
        }
    }

  allreduce_sum<double>(cost_data, 2 * NTopleaves * (NumTimeBinsToBeBalanced + 1), Communicator);
  /*
  MPI_Allreduce(MPI_IN_PLACE, cost_data, 2 * NTopleaves * (NumTimeBinsToBeBalanced + 1), MPI_DOUBLE, MPI_SUM, Communicator);
*/

#ifdef DOMAIN_SPECIAL_CHECK
  if(All.NumCurrentTiStep == 0 || All.NumCurrentTiStep == 2 || All.NumCurrentTiStep == 4)
    {
      if(ThisTask == 0)
        {
          char buf[MAXLEN_PATH_EXTRA];
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s/domain_data_0_step%d.txt", All.OutputDir, All.NumCurrentTiStep);
          FILE *fd = fopen(buf, "w");
          fprintf(fd, "%d %d\n", NTopleaves, NumTimeBinsToBeBalanced);
          for(int n = 0; n < NTopleaves; n++)
            {
              fprintf(fd, "%g  ", load[n]);
              for(int k = 0; k < NumTimeBinsToBeBalanced; k++)
                fprintf(fd, "%g  ", binGravCost[k * NTopleaves + n]);
              fprintf(fd, "\n");
            }
          fclose(fd);
        }
    }
#endif

  /* let's now find the optimum combination */

  /* first, enumerate the possibilities that we are going to try */

  int cnt_try                   = 0;
  balance_try_data *balance_try = NULL;

  for(int rep = 0; rep < 2; rep++)  // we repeat this twice, first just counting, then allocating and filling the list
    {
      cnt_try = 0;

      double fac      = 0;
      int nextra      = 0;
      int nextra_prev = -1;

      while(nextra <= NTask)
        {
          if(nextra != nextra_prev)
            {
              double base_balance = ((double)(MultipleDomains * NTask)) / (MultipleDomains * NTask - nextra);

              double excess_balance = 0.01;

              while(base_balance + excess_balance < 2.5)
                {
                  if(rep == 1)
                    {
                      balance_try[cnt_try].nextra      = nextra;
                      balance_try[cnt_try].try_balance = base_balance + excess_balance;
                    }

                  cnt_try++;

                  excess_balance *= 1.25;
                }

              nextra_prev = nextra;
            }

          if(fac == 0)
            fac = 0.01;
          else
            fac *= 1.25;

          nextra = fac * NTask;
        }

      if(rep == 0)
        balance_try = (balance_try_data *)Mem.mymalloc("balance_try", cnt_try * sizeof(balance_try_data));
    }

  if(NumNodes == 0)
    determine_compute_nodes();

  domain_printf("DOMAIN: we are going to try at most %d different settings for combining the domains on tasks=%d, nnodes=%d\n",
                cnt_try, NTask, NumNodes);

  /* sort the combinations such that we first try those yielding a lower imbalance */
  mycxxsort(balance_try, balance_try + cnt_try, domain_compare_trybalance);

  int start_try = 0;
  int completed = 0;
  int niter     = 0;

  while(completed == 0 && start_try < cnt_try)
    {
      double glob_max_cost = 0;

      if(start_try + ThisTask < cnt_try)
        {
          int nextra         = balance_try[start_try + ThisTask].nextra;
          double try_balance = balance_try[start_try + ThisTask].try_balance;

          int ndomains = MultipleDomains * NTask - nextra;

          domain_segments_data *domainAssign =
              (domain_segments_data *)Mem.mymalloc_clear("domainAssign", ndomains * sizeof(domain_segments_data));

          /* consolidate the finely split domains pieces into larger chunks, exactly ndomains of them */
          {
            double max_cost       = 0;
            double costhalfnode   = (0.5 * TotalCost) / NTopleaves;
            double costavg        = TotalCost / ndomains;
            double cost_before    = 0;
            double costavg_before = 0;
            int start             = 0;

            double total_cost = 0, total_load = 0;

            for(int n = 0; n < ndomains; n++)
              {
                int end = start;

                double cost = domain_leaf_cost[end].Cost;

                while((cost + cost_before + (end + 1 < NTopleaves ? domain_leaf_cost[end + 1].Cost : 0) <
                       costavg + costavg_before + costhalfnode) ||
                      (n == ndomains - 1 && end < NTopleaves - 1))
                  {
                    if((NTopleaves - end) > (ndomains - n))
                      end++;
                    else
                      break;

                    cost += domain_leaf_cost[end].Cost;
                  }

                domainAssign[n].start = start;
                domainAssign[n].end   = end;

                /* let's also determine the grav-cost and hydro-cost separately for each timebin of all the domain-pieces */
                for(int no = domainAssign[n].start; no <= domainAssign[n].end; no++)
                  {
                    domainAssign[n].load += load[no];
                    domainAssign[n].loadsph += loadsph[no];

                    total_load += load[no] + loadsph[no];
                    total_cost += load[no] + loadsph[no];

                    for(int i = 0; i < NumTimeBinsToBeBalanced; i++)
                      {
                        domainAssign[n].bin_GravCost[i] += binGravCost[i * NTopleaves + no];
                        domainAssign[n].bin_HydroCost[i] += binHydroCost[i * NTopleaves + no];

                        total_cost += binGravCost[i * NTopleaves + no] + binHydroCost[i * NTopleaves + no];
                      }
                  }

                cost_before += cost;
                costavg_before += costavg;

                start = end + 1;

                if(max_cost < cost)
                  max_cost = cost;
              }

            domain_printf("DOMAIN: total_cost=%g  total_load=%g\n", total_cost, total_load);
          }

          /* now start to map the domain pieces onto different tasks
           */

          struct tasklist_data
          {
            double bin_GravCost[TIMEBINS];
            double bin_HydroCost[TIMEBINS];
            double load;
            double loadsph;
          };

          tasklist_data *tasklist = (tasklist_data *)Mem.mymalloc_clear("tasklist", NTask * sizeof(tasklist_data));

          int n_cost_items = 0;
          cost_queue_data *cost_queues[2 * TIMEBINS + 2];
          int first_unusued_in_cost_queue[2 * TIMEBINS + 2];

          for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
            {
              if(GravCostPerListedTimeBin[n] > 0.0)
                {
                  cost_queues[n_cost_items] = (cost_queue_data *)Mem.mymalloc("cost_queues", ndomains * sizeof(cost_queue_data));
                  for(int i = 0; i < ndomains; i++)
                    {
                      cost_queues[n_cost_items][i].value = domainAssign[i].bin_GravCost[n];
                      cost_queues[n_cost_items][i].index = i;
                    }
#ifdef SIMPLE_DOMAIN_AGGREGATION
                  domain_determinate_aggregated_value(cost_queues[n_cost_items], ndomains);
#endif
                  mycxxsort(cost_queues[n_cost_items], cost_queues[n_cost_items] + ndomains, domain_sort_cost_queue_data);
                  first_unusued_in_cost_queue[n_cost_items] = 0;

                  n_cost_items++;
                }

              if(HydroCostNormFactors[n] > 0.0)
                {
                  cost_queues[n_cost_items] = (cost_queue_data *)Mem.mymalloc("cost_queues", ndomains * sizeof(cost_queue_data));
                  for(int i = 0; i < ndomains; i++)
                    {
                      cost_queues[n_cost_items][i].value = domainAssign[i].bin_HydroCost[n];
                      cost_queues[n_cost_items][i].index = i;
                    }
#ifdef SIMPLE_DOMAIN_AGGREGATION
                  domain_determinate_aggregated_value(cost_queues[n_cost_items], ndomains);
#endif
                  mycxxsort(cost_queues[n_cost_items], cost_queues[n_cost_items] + ndomains, domain_sort_cost_queue_data);
                  first_unusued_in_cost_queue[n_cost_items] = 0;

                  n_cost_items++;
                }
            }

          if(NormFactorLoad > 0.0)
            {
              cost_queues[n_cost_items] = (cost_queue_data *)Mem.mymalloc("cost_queues", ndomains * sizeof(cost_queue_data));
              for(int i = 0; i < ndomains; i++)
                {
                  cost_queues[n_cost_items][i].value = domainAssign[i].load;
                  cost_queues[n_cost_items][i].index = i;
                }
#ifdef SIMPLE_DOMAIN_AGGREGATION
              domain_determinate_aggregated_value(cost_queues[n_cost_items], ndomains);
#endif
              mycxxsort(cost_queues[n_cost_items], cost_queues[n_cost_items] + ndomains, domain_sort_cost_queue_data);
              first_unusued_in_cost_queue[n_cost_items] = 0;

              n_cost_items++;
            }

          if(NormFactorLoadSph > 0.0)
            {
              cost_queues[n_cost_items] = (cost_queue_data *)Mem.mymalloc("cost_queues", ndomains * sizeof(cost_queue_data));
              for(int i = 0; i < ndomains; i++)
                {
                  cost_queues[n_cost_items][i].value = domainAssign[i].loadsph;
                  cost_queues[n_cost_items][i].index = i;
                }
#ifdef SIMPLE_DOMAIN_AGGREGATION
              domain_determinate_aggregated_value(cost_queues[n_cost_items], ndomains);
#endif
              mycxxsort(cost_queues[n_cost_items], cost_queues[n_cost_items] + ndomains, domain_sort_cost_queue_data);
              first_unusued_in_cost_queue[n_cost_items] = 0;

              n_cost_items++;
            }

          int nextqueue     = 0;
          int ndomains_left = ndomains;
          int target        = 0;

          for(target = 0; target < NTask && ndomains_left > 0; target++)
            {
              int count    = 0;  // number of pieces added to this target task
              int failures = 0;

              while(ndomains_left > 0)
                {
                  int k = first_unusued_in_cost_queue[nextqueue];
                  while(domainAssign[cost_queues[nextqueue][k].index].used)
                    {
                      if(k == ndomains - 1)
                        Terminate("target=%d   nextqueue=%d  ndomains_left=%d  k == ndomains - 1", target, nextqueue, ndomains_left);

                      k++;
                      first_unusued_in_cost_queue[nextqueue]++;
                    }

                  /* this is our next candidate for adding */
                  int n = cost_queues[nextqueue][k].index;

                  nextqueue = (nextqueue + 1) % n_cost_items;

                  // Let's see what imbalance we would get
                  double max_cost = 0;

                  if(max_cost < tasklist[target].load + domainAssign[n].load)
                    max_cost = tasklist[target].load + domainAssign[n].load;

                  if(max_cost < tasklist[target].loadsph + domainAssign[n].loadsph)
                    max_cost = tasklist[target].loadsph + domainAssign[n].loadsph;

                  for(int bin = 0; bin < NumTimeBinsToBeBalanced; bin++)
                    {
                      if(max_cost < tasklist[target].bin_GravCost[bin] + domainAssign[n].bin_GravCost[bin])
                        max_cost = tasklist[target].bin_GravCost[bin] + domainAssign[n].bin_GravCost[bin];

                      if(max_cost < tasklist[target].bin_HydroCost[bin] + domainAssign[n].bin_HydroCost[bin])
                        max_cost = tasklist[target].bin_HydroCost[bin] + domainAssign[n].bin_HydroCost[bin];
                    }

                  if(count > 0)
                    {
                      if(max_cost * NTask > try_balance)
                        {
                          failures++;

                          if(failures > n_cost_items)
                            break;
                          else
                            continue;
                        }
                    }

                  if(max_cost > glob_max_cost)
                    glob_max_cost = max_cost;

                  domainAssign[n].task = target;
                  domainAssign[n].used = 1;

                  tasklist[target].load += domainAssign[n].load;
                  tasklist[target].loadsph += domainAssign[n].loadsph;
                  for(int bin = 0; bin < NumTimeBinsToBeBalanced; bin++)
                    {
                      tasklist[target].bin_GravCost[bin] += domainAssign[n].bin_GravCost[bin];
                      tasklist[target].bin_HydroCost[bin] += domainAssign[n].bin_HydroCost[bin];
                    }

                  ndomains_left--;
                  count++;

                  // if the following condition holds, no reason any further to try to stuff more than one piece together
                  if(ndomains_left < NTask - target)
                    break;
                }

              // do an extra skip, so that we do not typically start with the same queue
              nextqueue = (nextqueue + 1) % n_cost_items;
            }

          if(ndomains_left == 0)
            {
              domain_printf("DOMAIN: combining multiple-domains succeeded, target=%d  NTask=%d\n", target, NTask);
              completed = 1;

#ifdef DOMAIN_SPECIAL_CHECK
              if(All.NumCurrentTiStep == 0 || All.NumCurrentTiStep == 2 || All.NumCurrentTiStep == 4)
                {
                  char buf[MAXLEN_PATH_EXTRA];
                  snprintf(buf, "%s/domain_data_1_step%d_task%d.txt", All.OutputDir, All.NumCurrentTiStep, ThisTask);
                  FILE *fd = fopen(buf, "w");
                  fprintf(fd, "%d %d\n", ndomains, NumTimeBinsToBeBalanced);
                  for(int n = 0; n < ndomains; n++)
                    {
                      fprintf(fd, "%g  ", domainAssign[n].load);
                      for(int k = 0; k < NumTimeBinsToBeBalanced; k++)
                        fprintf(fd, "%g  ", domainAssign[n].bin_GravCost[k]);
                      fprintf(fd, "\n");
                    }
                  fclose(fd);
                }
              if(All.NumCurrentTiStep == 0 || All.NumCurrentTiStep == 2 || All.NumCurrentTiStep == 4)
                {
                  char buf[MAXLEN_PATH_EXTRA];
                  snprintf(buf, MAXLEN_PATH_EXTRA, "%s/domain_data_2_step%d_task%d.txt", All.OutputDir, All.NumCurrentTiStep,
                           ThisTask);
                  FILE *fd = fopen(buf, "w");
                  fprintf(fd, "%d %d\n", NTask, NumTimeBinsToBeBalanced);
                  for(int n = 0; n < NTask; n++)
                    {
                      fprintf(fd, "%g  ", tasklist[n].load);
                      for(int k = 0; k < NumTimeBinsToBeBalanced; k++)
                        fprintf(fd, "%g  ", tasklist[n].bin_GravCost[k]);
                      fprintf(fd, "\n");
                    }
                  fprintf(fd, "%g\n", glob_max_cost * NTask);
                  fclose(fd);
                }
#endif

              /* store the mapping of the topleaves to tasks */
              for(int n = 0; n < ndomains; n++)
                {
                  for(int i = domainAssign[n].start; i <= domainAssign[n].end; i++)
                    TaskOfLeaf[i] = domainAssign[n].task;
                }
            }
          else
            {
              glob_max_cost = MAX_DOUBLE_NUMBER;
            }

          for(int num = n_cost_items - 1; num >= 0; num--)
            Mem.myfree(cost_queues[num]);

          Mem.myfree(tasklist);
          Mem.myfree(domainAssign);
        }
      else
        glob_max_cost = MAX_DOUBLE_NUMBER; /* this processor was not doing anything */

      struct
      {
        double cost;
        int rank;
      } global = {glob_max_cost, ThisTask};

      MPI_Allreduce(MPI_IN_PLACE, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);

      MPI_Allreduce(MPI_IN_PLACE, &completed, 1, MPI_INT, MPI_MAX, Communicator);

      niter++;

      if(completed)
        {
          domain_printf(
              "DOMAIN: best solution found after %d iterations by task=%d for nextra=%d, reaching maximum imbalance of %g|%g\n", niter,
              global.rank, balance_try[start_try + global.rank].nextra, global.cost * NTask,
              balance_try[start_try + global.rank].try_balance);

          MPI_Bcast(TaskOfLeaf, NTopleaves, MPI_INT, global.rank, Communicator);
          break;
        }

      start_try += NTask;
    }

  if(completed == 0)
    Terminate("domain balancing failed");

  Mem.myfree(balance_try);

  Mem.myfree(cost_data);

  double t1 = Logs.second();
  domain_printf("DOMAIN: combining multiple-domains took %g sec\n", Logs.timediff(t0, t1));
}

#ifdef SIMPLE_DOMAIN_AGGREGATION
template <typename partset>
void domain<partset>::domain_determinate_aggregated_value(cost_queue_data *data, int ndomains)
{
  if(NumNodes < 1)
    Terminate("NumNodes=%d\n", NumNodes);

  int nbase      = ndomains / NumNodes;
  int additional = ndomains % NumNodes;
  int start      = 0;
  int end        = 0;

  for(int i = 0; i < NumNodes; i++)
    {
      start = end;
      end   = start + nbase;

      if(additional > 0)
        {
          end++;
          additional--;
        }

      double aggregated_value = 0;
      for(int n = start; n < end; n++)
        aggregated_value += data[n].value;

      for(int n = start; n < end; n++)
        data[n].aggregated_value = aggregated_value;
    }

  if(end != ndomains)
    Terminate("end != ndomains");
}
#endif

/** This function prepares the measurement of the total cost on each domain.
 *  In particular, we determine how the timebins are mapped to the explicit measurements
 *  of the gravity cost stored in the P.GravCost[] array (which in general will only be available for a subset
 *  of all timebins). For the unmatched timebins, a closest bin is selected that is the most similar in terms
 *  of particle number on the bin. Finally, the routine also determines how often each timebin is executed in
 *  one cycle associated with the highest occupied timebin.
 */
template <>
void domain<simparticles>::domain_init_sum_cost(void)
{
  long long tot_count[TIMEBINS], tot_count_sph[TIMEBINS];

  sumup_large_ints(TIMEBINS, Tp->TimeBinsGravity.TimeBinCount, tot_count, Communicator);
  sumup_large_ints(TIMEBINS, Tp->TimeBinsHydro.TimeBinCount, tot_count_sph, Communicator);

  for(int i = 0; i < TIMEBINS; i++)
    {
      domain_to_be_balanced[i] = 0;
      domain_grav_weight[i]    = 1;
      domain_hydro_weight[i]   = 1;
    }

  domain_to_be_balanced[All.HighestActiveTimeBin] = 1;
  domain_grav_weight[All.HighestActiveTimeBin]    = 1;
  domain_hydro_weight[All.HighestActiveTimeBin]   = 1;

  ListOfTimeBinsToBeBalanced[0] = All.HighestActiveTimeBin;
  NumTimeBinsToBeBalanced       = 1;

  if(Mode == COLL_SUBFIND)
    return;

#ifdef HIERARCHICAL_GRAVITY

  for(int j = All.HighestActiveTimeBin - 1; j >= All.LowestOccupiedTimeBin; j--)
    {
      if(tot_count[j] > 0 || tot_count_sph[j] > 0)
        {
          ListOfTimeBinsToBeBalanced[NumTimeBinsToBeBalanced++] = j;

          domain_to_be_balanced[j] = 1;
        }

      domain_grav_weight[j] += 2;
    }

  for(int i = All.SmallestTimeBinWithDomainDecomposition - 1, weight = 1; i >= All.LowestOccupiedTimeBin; i--, weight *= 2)
    {
      if(tot_count[i] > 0)
        {
          domain_grav_weight[i] = weight;

          for(int j = i - 1; j >= All.LowestOccupiedTimeBin; j--)
            domain_grav_weight[j] += 2 * weight;
        }

      if(tot_count_sph[i] > 0)
        domain_hydro_weight[i] = weight;
    }

#else

  for(int i = All.SmallestTimeBinWithDomainDecomposition - 1, weight = 1; i >= All.LowestOccupiedTimeBin; i--, weight *= 2)
    {
      if(tot_count[i] > 0 || tot_count_sph[i] > 0)
        {
          ListOfTimeBinsToBeBalanced[NumTimeBinsToBeBalanced++] = i;
          domain_to_be_balanced[i]                              = 1;
        }

      if(tot_count[i] > 0)
        domain_grav_weight[i] = weight;

      if(tot_count_sph[i] > 0)
        domain_hydro_weight[i] = weight;
    }

#endif
}

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES)

template <>
void domain<lcparticles>::domain_init_sum_cost(void)
{
  for(int i = 0; i < TIMEBINS; i++)
    {
      domain_to_be_balanced[i] = 0;
      domain_grav_weight[i]    = 1;
      domain_hydro_weight[i]   = 1;
    }

  domain_to_be_balanced[0] = 1;
}

#endif

#include "../data/simparticles.h"
template class domain<simparticles>;

#ifdef LIGHTCONE_PARTICLES
#include "../data/lcparticles.h"
template class domain<lcparticles>;
#endif
