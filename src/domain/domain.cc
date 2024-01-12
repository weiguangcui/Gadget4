/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file domain.cc
 *
 *  \brief code for work-load balanced domain decomposition
 *
 *  This class contains the code for the domain decomposition of the
 *  simulation volume.  The domains are constructed from disjoint subsets
 *  of leaves of a fiducial top-level tree that covers the full
 *  simulation volume. Domain boundaries hence run along tree-node
 *  divisions of a fiducial global BH oct-tree. As a result of this method, the
 *  gravity forces are in principle strictly independent of the way the domains
 *  are cut. The domain decomposition can be carried out for an arbitrary
 *  number of CPUs. Individual domain pieces are not cubical, but spatially
 *  coherent since the leaves are traversed in a Peano-Hilbert order and
 *  individual domains form segments along this order.  This also ensures
 *  that each domain has a small surface to volume ratio, which reduces
 *  communication.
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
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../sort/peano.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/*! This is the main routine for the domain decomposition.  It acts as a
 *  driver routine that allocates various temporary buffers, maps the
 *  particles back onto the periodic box if needed, and then does the
 *  domain decomposition, followed by a final Peano-Hilbert order of all particles
 *  as a tuning measure.
 */
template <typename partset>
void domain<partset>::domain_decomposition(domain_options mode)
{
  Mode = mode;

  TIMER_START(CPU_DOMAIN);

  double t0 = Logs.second();

  domain_printf("DOMAIN: Begin domain decomposition (sync-point %d).\n", All.NumCurrentTiStep);

  /* map the particles by a shift vector back if desired */
  do_box_wrapping();

  /* determine which time bins need to be balanced */
  domain_init_sum_cost();

  /* find total cost factors, including a determination of MultipleDomains */
  domain_find_total_cost();

  /* allocate some fields that need to stay */
  domain_allocate();

  /* allocate some arrays we are going to use */
  domain_key = (peanokey *)Mem.mymalloc_movable(&domain_key, "domain_key", (sizeof(peanokey) * Tp->MaxPart));
  domain_leaf_cost =
      (domain_cost_data *)Mem.mymalloc_movable(&domain_leaf_cost, "domain_leaf_cost", (MaxTopNodes * sizeof(domain_cost_data)));

  topNodes = (local_topnode_data *)Mem.mymalloc_movable(&topNodes, "topNodes", (MaxTopNodes * sizeof(local_topnode_data)));

  /* determine top-level tree */
  domain_determineTopTree();

  /* combine on each MPI task several of the domains (namely the number MultipleDomains) */
  domain_combine_multipledomains();

  Mem.myfree(topNodes);
  Mem.myfree(domain_leaf_cost);

  /* move stars out of gas block if present */
  domain_rearrange_particle_sequence();

  /* finally, carry out the actual particle exchange */
  if(Mode == STANDARD)
    domain_exchange();
  else if(Mode == COLL_SUBFIND)
    domain_coll_subfind_prepare_exchange();

  double t1 = Logs.second();

  domain_printf("DOMAIN: domain decomposition done. (took in total %g sec)\n", Logs.timediff(t0, t1));

  if(Mode == STANDARD)
    {
      TIMER_STOPSTART(CPU_DOMAIN, CPU_PEANO);

      peano_hilbert_order(domain_key);

      TIMER_STOPSTART(CPU_PEANO, CPU_DOMAIN);
    }

  Mem.myfree(domain_key);

  Mem.myfree(ListOfTopleaves);

  TaskOfLeaf = (int *)Mem.myrealloc_movable(TaskOfLeaf, NTopleaves * sizeof(int));
  TopNodes   = (topnode_data *)Mem.myrealloc_movable(TopNodes, NTopnodes * sizeof(topnode_data));

  ListOfTopleaves = (int *)Mem.mymalloc_movable(&ListOfTopleaves, "ListOfTopleaves", (NTopleaves * sizeof(int)));

  memset(NumTopleafOfTask, 0, NTask * sizeof(int));

  for(int i = 0; i < NTopleaves; i++)
    NumTopleafOfTask[TaskOfLeaf[i]]++;

  FirstTopleafOfTask[0] = 0;
  for(int i = 1; i < NTask; i++)
    FirstTopleafOfTask[i] = FirstTopleafOfTask[i - 1] + NumTopleafOfTask[i - 1];

  memset(NumTopleafOfTask, 0, NTask * sizeof(int));

  for(int i = 0; i < NTopleaves; i++)
    {
      int task             = TaskOfLeaf[i];
      int off              = FirstTopleafOfTask[task] + NumTopleafOfTask[task]++;
      ListOfTopleaves[off] = i;
    }

  if(Mode == STANDARD)
    {
      /* the following will reconstruct the timebins and report the balance
       * in case we are dealing with simparticles, for lcparticles nothing will be done
       */
      domain_report_balance();
    }

#ifdef DOMAIN_SPECIAL_CHECK
  if(ThisTask == 0 && All.NumCurrentTiStep == 4)
    Terminate("stop");
#endif

  TIMER_STOP(CPU_DOMAIN);
}

/*! This function allocates all the stuff that will be required for the tree-construction/walk later on */
template <typename partset>
void domain<partset>::domain_allocate(int maxtopnodes)
{
  MaxTopNodes = maxtopnodes;

  if(FirstTopleafOfTask)
    Terminate("domain storage already allocated");

  FirstTopleafOfTask = (int *)Mem.mymalloc_movable(&FirstTopleafOfTask, "FirstTopleafOfTask", NTask * sizeof(int));
  NumTopleafOfTask   = (int *)Mem.mymalloc_movable(&NumTopleafOfTask, "NumTopleafOfTask", NTask * sizeof(int));
  TopNodes           = (topnode_data *)Mem.mymalloc_movable(&TopNodes, "TopNodes", (MaxTopNodes * sizeof(topnode_data)));
  TaskOfLeaf         = (int *)Mem.mymalloc_movable(&TaskOfLeaf, "TaskOfLeaf", (MaxTopNodes * sizeof(int)));
  ListOfTopleaves    = (int *)Mem.mymalloc_movable(&ListOfTopleaves, "DomainListOfLocalTopleaves", (MaxTopNodes * sizeof(int)));
}

/*! This function allocates all the stuff that will be required for the tree-construction/walk later on */
template <typename partset>
void domain<partset>::domain_allocate(void)
{
  int maxtopnodes = All.TopNodeAllocFactor * std::max<int>(All.TopNodeFactor * MultipleDomains * NTask, BASENUMBER);
  domain_allocate(maxtopnodes);
}

template <typename partset>
void domain<partset>::domain_free(void)
{
  if(!FirstTopleafOfTask)
    Terminate("domain storage not allocated");

  Mem.myfree_movable(ListOfTopleaves);
  Mem.myfree_movable(TaskOfLeaf);
  Mem.myfree_movable(TopNodes);
  Mem.myfree_movable(NumTopleafOfTask);
  Mem.myfree_movable(FirstTopleafOfTask);

  ListOfTopleaves    = NULL;
  TaskOfLeaf         = NULL;
  TopNodes           = NULL;
  NumTopleafOfTask   = NULL;
  FirstTopleafOfTask = NULL;
}

template <typename partset>
void domain<partset>::domain_printf(char *buf)
{
  if(All.RestartFlag == RST_BEGIN || All.RestartFlag == RST_RESUME || All.RestartFlag == RST_STARTFROMSNAP)
    {
      fprintf(Logs.FdDomain, "%s", buf);
    }
}

template <>
void domain<simparticles>::domain_find_total_cost(void)
{
  /* for each timebin that should be balanced, collect the gravity cost of
   * the particles active in that timebin
   */

  for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
    {
      GravCostPerListedTimeBin[n]    = 0.0;
      MaxGravCostPerListedTimeBin[n] = 0.0;
      /* do the same for the hydrodynamics
       */
      HydroCostPerListedTimeBin[n] = 0.0;
    }

  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(Tp->P[i].GravCost == 0)
        Tp->P[i].GravCost = 1;

      for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
        {
          int bin = ListOfTimeBinsToBeBalanced[n];

          if(bin >= Tp->P[i].TimeBinGrav)
            {
              GravCostPerListedTimeBin[n] += Tp->P[i].GravCost;
              if(MaxGravCostPerListedTimeBin[n] < Tp->P[i].GravCost)
                MaxGravCostPerListedTimeBin[n] = Tp->P[i].GravCost;
            }

          if(Tp->P[i].getType() == 0)
            {
              if(bin >= Tp->P[i].getTimeBinHydro())
                HydroCostPerListedTimeBin[n] += 1.0;
            }
        }
    }

  long long sum[2] = {Tp->NumPart, Tp->NumGas};

  MPI_Allreduce(MPI_IN_PLACE, sum, 2, MPI_LONG_LONG, MPI_SUM, Communicator);

  NormFactorLoad    = 1.0 / sum[0];
  NormFactorLoadSph = sum[1] > 0.0 ? 1.0 / sum[1] : 0.0;

  MultipleDomains = 0;
  TotalCost       = 0.0;

  if(NormFactorLoad > 0.0)
    {
      MultipleDomains += 1;
      TotalCost += 1.0;
    }

  if(NormFactorLoadSph > 0.0)
    {
      MultipleDomains += 1;
      TotalCost += 1.0;
    }

  MPI_Allreduce(MPI_IN_PLACE, GravCostPerListedTimeBin, NumTimeBinsToBeBalanced, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(MPI_IN_PLACE, HydroCostPerListedTimeBin, NumTimeBinsToBeBalanced, MPI_DOUBLE, MPI_SUM, Communicator);

  MPI_Allreduce(MPI_IN_PLACE, MaxGravCostPerListedTimeBin, NumTimeBinsToBeBalanced, MPI_DOUBLE, MPI_MAX, Communicator);

  double limit = 1.0 / (All.TopNodeFactor * MultipleDomains * NTask);

  for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
    {
      if(GravCostPerListedTimeBin[n] > 0.0)
        MultipleDomains += 1;

      if(HydroCostPerListedTimeBin[n] > 0.0)
        MultipleDomains += 1;

      GravCostNormFactors[n]  = GravCostPerListedTimeBin[n] > 0.0 ? 1.0 / GravCostPerListedTimeBin[n] : 0.0;
      HydroCostNormFactors[n] = HydroCostPerListedTimeBin[n] > 0.0 ? 1.0 / HydroCostPerListedTimeBin[n] : 0.0;

      if(MaxGravCostPerListedTimeBin[n] * GravCostNormFactors[n] > limit)
        GravCostNormFactors[n] = limit / MaxGravCostPerListedTimeBin[n];

      if(HydroCostNormFactors[n] > limit)
        HydroCostNormFactors[n] = limit;

      TotalCost += GravCostPerListedTimeBin[n] * GravCostNormFactors[n];
      TotalCost += HydroCostPerListedTimeBin[n] * HydroCostNormFactors[n];
    }
}

#ifdef LIGHTCONE_PARTICLES
template <>
void domain<lcparticles>::domain_find_total_cost(void)
{
  long long sum[2] = {Tp->NumPart, Tp->NumGas};

  MPI_Allreduce(MPI_IN_PLACE, sum, 2, MPI_LONG_LONG, MPI_SUM, Communicator);

  NormFactorLoad    = 1.0 / sum[0];
  NormFactorLoadSph = sum[1] > 0.0 ? 1.0 / sum[1] : 0.0;

  MultipleDomains = 0;
  TotalCost       = 0.0;

  if(NormFactorLoad > 0.0)
    {
      MultipleDomains += 1;
      TotalCost += 1.0;
    }

  if(NormFactorLoadSph > 0.0)
    {
      MultipleDomains += 1;
      TotalCost += 1.0;
    }

  NumTimeBinsToBeBalanced = 0;
}
#endif

template <>
void domain<simparticles>::domain_rearrange_particle_sequence(void)
{
  if(Mode != STANDARD)
    return;

  for(int i = 0; i < Tp->NumGas; i++)
    if(Tp->P[i].getType() != 0) /*If not a gas particle, swap to the end of the list */
      {
        particle_data psave = Tp->P[i];
        peanokey key        = domain_key[i];

        Tp->P[i]      = Tp->P[Tp->NumGas - 1];
        Tp->SphP[i]   = Tp->SphP[Tp->NumGas - 1];
        domain_key[i] = domain_key[Tp->NumGas - 1];

        Tp->P[Tp->NumGas - 1]      = psave;
        domain_key[Tp->NumGas - 1] = key;

        Tp->NumGas--;
        i--;
      }
  /*Now we have rearranged the particles,
   *we don't need to do it again unless there are more stars*/
}

template <>
void domain<simparticles>::domain_report_balance(void)
{
  Tp->reconstruct_timebins();

  TIMER_STOPSTART(CPU_DOMAIN, CPU_LOGS);

  if(Mode != STANDARD)
    return;

  /* now execute code to report balance */

  /* get total particle counts */
  long long loc_count[2 * TIMEBINS], glob_count[2 * TIMEBINS];

  for(int i = 0; i < TIMEBINS; i++)
    {
      loc_count[i]            = Tp->TimeBinsGravity.TimeBinCount[i];
      loc_count[TIMEBINS + i] = Tp->TimeBinsHydro.TimeBinCount[i];
    }

  MPI_Reduce(loc_count, glob_count, 2 * TIMEBINS, MPI_LONG_LONG_INT, MPI_SUM, 0, Communicator);

  double loc_max_data[2 * TIMEBINS + 3], glob_max_data[2 * TIMEBINS + 3];
  loc_max_data[2 * TIMEBINS + 0] = Tp->NumPart;
  loc_max_data[2 * TIMEBINS + 1] = Tp->NumGas;
  loc_max_data[2 * TIMEBINS + 2] = Tp->NumPart - Tp->NumGas;

  double glob_sum_data[2 * TIMEBINS];

  double *loc_HydroCost  = &loc_max_data[0];
  double *loc_GravCost   = &loc_max_data[TIMEBINS];
  double *max_HydroCost  = &glob_max_data[0];
  double *max_GravCost   = &glob_max_data[TIMEBINS];
  double *glob_HydroCost = &glob_sum_data[0];
  double *glob_GravCost  = &glob_sum_data[TIMEBINS];

  for(int i = 0; i < TIMEBINS; i++)
    {
      loc_GravCost[i]  = 0;
      loc_HydroCost[i] = 0;
    }

#ifdef SELFGRAVITY
  for(int i = 0; i < Tp->NumPart; i++)
    {
      for(int bin = Tp->P[i].TimeBinGrav; bin <= All.HighestOccupiedTimeBin; bin++)
        {
          loc_GravCost[bin] += MIN_FLOAT_NUMBER + domain_grav_weight[bin] * Tp->P[i].GravCost;
        }
    }
#endif

  for(int i = 0; i < Tp->NumPart; i++)
    if(Tp->P[i].getType() == 0)
      loc_HydroCost[Tp->P[i].getTimeBinHydro()] += 1.0;

  /* now determine the cumulative cost for the hydrodynamics */
  for(int i = 1; i <= All.HighestOccupiedTimeBin; i++)
    loc_HydroCost[i] += loc_HydroCost[i - 1];

  MPI_Reduce(loc_max_data, glob_sum_data, 2 * TIMEBINS, MPI_DOUBLE, MPI_SUM, 0, Communicator);
  MPI_Reduce(loc_max_data, glob_max_data, 2 * TIMEBINS + 3, MPI_DOUBLE, MPI_MAX, 0, Communicator);

  if(ThisTask == 0)
    {
      double max_tot = glob_max_data[2 * TIMEBINS + 0];
      double max_sph = glob_max_data[2 * TIMEBINS + 1];
      double max_dm  = glob_max_data[2 * TIMEBINS + 2];

      long long *tot_count     = &glob_count[0];
      long long *tot_count_sph = &glob_count[TIMEBINS];

      long long tot_cumulative[TIMEBINS];
      tot_cumulative[0] = tot_count[0];

      for(int i = 1; i < TIMEBINS; i++)
        tot_cumulative[i] = tot_count[i] + tot_cumulative[i - 1];

      double tot_gravcost = 0, max_gravcost = 0, tot_hydrocost = 0, max_hydrocost = 0;

      for(int i = 0; i < TIMEBINS; i++)
        {
          tot_gravcost += domain_to_be_balanced[i] * glob_GravCost[i] / NTask;
          max_gravcost += domain_to_be_balanced[i] * max_GravCost[i];

          tot_hydrocost += domain_to_be_balanced[i] * glob_HydroCost[i] / NTask;
          max_hydrocost += domain_to_be_balanced[i] * max_HydroCost[i];
        }

      double bal_grav_bin[TIMEBINS], bal_grav_bin_rel[TIMEBINS];
      double bal_hydro_bin[TIMEBINS], bal_hydro_bin_rel[TIMEBINS];

      for(int i = 0; i < TIMEBINS; i++)
        {
          if(tot_count[i] > 0)
            {
              bal_grav_bin[i] = max_GravCost[i] / (glob_GravCost[i] / NTask + SMALLNUM);
              bal_grav_bin_rel[i] =
                  (tot_gravcost + domain_to_be_balanced[i] * (max_GravCost[i] - glob_GravCost[i] / NTask)) / (tot_gravcost + SMALLNUM);
            }
          else
            {
              bal_grav_bin[i]     = 0.0;
              bal_grav_bin_rel[i] = 0.0;
            }

          if(tot_count_sph[i] > 0)
            {
              bal_hydro_bin[i]     = max_HydroCost[i] / (glob_HydroCost[i] / NTask + SMALLNUM);
              bal_hydro_bin_rel[i] = (tot_hydrocost + domain_to_be_balanced[i] * (max_HydroCost[i] - glob_HydroCost[i] / NTask)) /
                                     (tot_hydrocost + SMALLNUM);
            }
          else
            {
              bal_hydro_bin[i]     = 0.0;
              bal_hydro_bin_rel[i] = 0.0;
            }
        }

      char buf[MAXLEN_PATH];
      snprintf(buf, MAXLEN_PATH, "\nDOMAIN BALANCE, Sync-Point %d, Time: %g\n", All.NumCurrentTiStep, All.Time);
      domain_printf(buf);
      snprintf(buf, MAXLEN_PATH, "Timebins:       Gravity       Hydro  cumulative      grav-balance       hydro-balance\n");
      domain_printf(buf);

      long long tot = 0, tot_sph = 0;

      for(int i = TIMEBINS - 1; i >= 0; i--)
        {
#if(defined(SELFGRAVITY) || defined(EXTERNALGRAVITY))
          if(tot_count_sph[i] > 0 || tot_count[i] > 0)
#else
          if(tot_count[i] > 0)
            tot += tot_count[i];

          if(tot_count_sph[i] > 0)
#endif
            {
              char buf[MAXLEN_PATH];
              snprintf(buf, MAXLEN_PATH, "%c%cbin=%2d     %10llu  %10llu  %10llu    %6.3f |%6.3f  %c   %6.3f |%6.3f\n",
                       i == All.HighestActiveTimeBin ? '>' : ' ', i >= All.SmallestTimeBinWithDomainDecomposition ? '|' : ' ', i,
                       tot_count[i], tot_count_sph[i], tot_cumulative[i], bal_grav_bin[i], bal_grav_bin_rel[i],
                       domain_to_be_balanced[i] > 0 ? '*' : ' ', bal_hydro_bin[i], bal_hydro_bin_rel[i]);
              domain_printf(buf);

              tot += tot_count[i];
              tot_sph += tot_count_sph[i];
            }
        }

      snprintf(buf, MAXLEN_PATH, "-------------------------------------------------------------------------------------\n");
      domain_printf(buf);
      snprintf(buf, MAXLEN_PATH, "BALANCE,  LOAD:  %6.3f      %6.3f      %6.3f  WORK:     %6.3f              %6.3f\n",
               max_dm / (tot - tot_sph + SMALLNUM) * NTask, max_sph / (tot_sph + SMALLNUM) * NTask, max_tot / (tot + SMALLNUM) * NTask,
               max_gravcost / (tot_gravcost + SMALLNUM), max_hydrocost / (tot_hydrocost + SMALLNUM));
      domain_printf(buf);
      snprintf(buf, MAXLEN_PATH, "-------------------------------------------------------------------------------------\n");
      domain_printf(buf);
      snprintf(buf, MAXLEN_PATH, "\n");
      domain_printf(buf);
      myflush(Logs.FdDomain);
    }

  TIMER_STOPSTART(CPU_LOGS, CPU_DOMAIN);
}

#ifdef LIGHTCONE_PARTICLES

template <>
void domain<lcparticles>::domain_rearrange_particle_sequence(void)
{
}

template <>
void domain<lcparticles>::domain_report_balance(void)
{
}

#endif

#include "../data/simparticles.h"
template class domain<simparticles>;

#ifdef LIGHTCONE_PARTICLES
#include "../data/lcparticles.h"
template class domain<lcparticles>;
#endif
