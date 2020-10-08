/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file domain_toplevel.cc
 *
 *  \brief construction of top-level for subdividing the volume into domains
 */

#include "gadgetconfig.h"

#include <mpi.h>
#include <algorithm>
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
#include "../sort/peano.h"
#include "../system/system.h"

template <typename partset>
void domain<partset>::domain_do_local_refine(int n, int *list)
{
  double *worklist     = (double *)Mem.mymalloc("worklist", 8 * n * sizeof(double));
  long long *countlist = (long long *)Mem.mymalloc("countlist", 8 * n * sizeof(long long));

  /* create the new nodes */
  for(int k = 0; k < n; k++)
    {
      int i                = list[k];
      TopNodes[i].Daughter = NTopnodes;
      NTopnodes += 8;
      NTopleaves += 7;

      for(int j = 0; j < 8; j++)
        {
          int sub = TopNodes[i].Daughter + j;

          TopNodes[sub].Daughter = -1;
          topNodes[sub].Level    = topNodes[i].Level + 1;
          topNodes[sub].StartKey = topNodes[i].StartKey + get_peanokey_offset(j, (3 * (BITS_FOR_POSITIONS - topNodes[sub].Level)));
          topNodes[sub].PIndex   = topNodes[i].PIndex;
          topNodes[sub].Count    = 0;
          topNodes[sub].Cost     = 0;
        }

      int sub = TopNodes[i].Daughter;

      for(int p = topNodes[i].PIndex, j = 0; p < topNodes[i].PIndex + topNodes[i].Count; p++)
        {
          if(j < 7)
            while(mp[p].key >= topNodes[sub + 1].StartKey)
              {
                j++;
                sub++;
                topNodes[sub].PIndex = p;
                if(j >= 7)
                  break;
              }

          topNodes[sub].Cost += mp[p].cost;
          topNodes[sub].Count++;
        }

      for(int j = 0; j < 8; j++)
        {
          int sub              = TopNodes[i].Daughter + j;
          worklist[k * 8 + j]  = topNodes[sub].Cost;
          countlist[k * 8 + j] = topNodes[sub].Count;
        }
    }

  allreduce_sum<double>(worklist, 8 * n, Communicator);
  allreduce_sum<long long>(countlist, 8 * n, Communicator);

  /* store the results in the corresponding top nodes */
  for(int k = 0; k < n; k++)
    {
      int i = list[k];

      for(int j = 0; j < 8; j++)
        {
          int sub                = TopNodes[i].Daughter + j;
          topNodes[sub].Cost     = worklist[k * 8 + j];
          topNodes[sub].CountTot = countlist[k * 8 + j];
        }
    }

  Mem.myfree(countlist);
  Mem.myfree(worklist);
}

/*! This function walks the global top tree in order to establish the
 *  number of leaves it has, and for assigning the leaf numbers along the
 *  Peano-Hilbert Curve. These leaves are later combined to domain pieces,
 *  which are distributed to different processors.
 */
template <typename partset>
void domain<partset>::domain_walktoptree(int no)
{
  if(TopNodes[no].Daughter == -1)
    {
      TopNodes[no].Leaf = NTopleaves;

      domain_leaf_cost[TopNodes[no].Leaf].Cost = topNodes[no].Cost;

      NTopleaves++;
    }
  else
    {
      for(int i = 0; i < 8; i++)
        domain_walktoptree(TopNodes[no].Daughter + i);
    }
}

template <>
double domain<simparticles>::domain_get_cost_summed_over_timebins(int i)
{
  double cost = 0;

  for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
    {
      int bin = ListOfTimeBinsToBeBalanced[n];

      if(bin >= Tp->P[i].TimeBinGrav)
        cost += GravCostNormFactors[n] * Tp->P[i].GravCost;

      if(Tp->P[i].getType() == 0)
        {
#ifdef SUBFIND
          if(Mode == COLL_SUBFIND)
            {
              if(Tp->PS[i].DomainFlag)
                cost += HydroCostNormFactors[n];
            }
          else
#endif
            {
              if(bin >= Tp->P[i].getTimeBinHydro())
                cost += HydroCostNormFactors[n];
            }
        }
    }

  return cost;
}

#ifdef LIGHTCONE_PARTICLES
template <>
double domain<lcparticles>::domain_get_cost_summed_over_timebins(int i)
{
  return 0;
}
#endif

/*! This function constructs the global top-level tree node that is used
 *  for the domain decomposition. This is done by considering the string of
 *  Peano-Hilbert keys for all particles, which is recursively chopped off
 *  in pieces of eight segments until each segment holds at most a certain
 *  number of particles.
 */
template <typename partset>
void domain<partset>::domain_determineTopTree(void)
{
  double t0           = Logs.second();
  int message_printed = 0;

  mp = (domain_peano_hilbert_data *)Mem.mymalloc_movable(&mp, "mp", sizeof(domain_peano_hilbert_data) * Tp->NumPart);

  int count  = 0;
  double sum = 0.0;

  for(int i = 0; i < Tp->NumPart; i++)
    {
      mp[i].key = domain_key[i] = peano_hilbert_key(Tp->P[i].IntPos[0], Tp->P[i].IntPos[1], Tp->P[i].IntPos[2], BITS_FOR_POSITIONS);
      mp[i].cost                = 0;
      count++;

      mp[i].cost += domain_get_cost_summed_over_timebins(i);

      mp[i].cost += NormFactorLoad;

      if(Tp->P[i].getType() == 0)
        mp[i].cost += NormFactorLoadSph;

      sum += mp[i].cost;
    }

  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  domain_printf("DOMAIN: Sum=%g  TotalCost=%g  NumTimeBinsToBeBalanced=%d  MultipleDomains=%d\n", sum, TotalCost,
                NumTimeBinsToBeBalanced, MultipleDomains);

  mycxxsort(mp, mp + Tp->NumPart, domain_compare_key);

  NTopnodes            = 1;
  NTopleaves           = 1;
  TopNodes[0].Daughter = -1;
  topNodes[0].Level    = 0;
  topNodes[0].StartKey = {0, 0, 0};
  topNodes[0].PIndex   = 0;
  topNodes[0].Count    = count; /* this is the local number */
  topNodes[0].CountTot = Tp->TotNumPart;
  topNodes[0].Cost     = sum;

  /* in list[], we store the node indices hat should be refined */
  int *list = (int *)Mem.mymalloc_movable(&list, "list", MaxTopNodes * sizeof(int));

  double limit = 1.0 / (All.TopNodeFactor * NTask);

  int iter = 0;

  do
    {
      count = 0;

      for(int n = 0; n < NTopnodes; n++)
        if(TopNodes[n].Daughter == -1)  // consider only leaf nodes
          {
            if(topNodes[n].CountTot > 1)  // only split nodes with at list 1 particles
              if(topNodes[n].Cost > limit)
                {
                  while(NTopnodes + 8 * (count + 1) > MaxTopNodes)
                    {
                      domain_printf("DOMAIN: Increasing TopNodeAllocFactor=%g  ", All.TopNodeAllocFactor);
                      All.TopNodeAllocFactor *= 1.3;
                      domain_printf("new value=%g\n", All.TopNodeAllocFactor);
                      if(All.TopNodeAllocFactor > 1000)
                        Terminate("something seems to be going seriously wrong here. Stopping.\n");

                      MaxTopNodes = All.TopNodeAllocFactor * std::max<int>(All.TopNodeFactor * MultipleDomains * NTask, BASENUMBER);

                      topNodes   = (local_topnode_data *)Mem.myrealloc_movable(topNodes, (MaxTopNodes * sizeof(local_topnode_data)));
                      TopNodes   = (topnode_data *)Mem.myrealloc_movable(TopNodes, (MaxTopNodes * sizeof(topnode_data)));
                      TaskOfLeaf = (int *)Mem.myrealloc_movable(TaskOfLeaf, (MaxTopNodes * sizeof(int)));
                      domain_leaf_cost =
                          (domain_cost_data *)Mem.myrealloc_movable(domain_leaf_cost, (MaxTopNodes * sizeof(domain_cost_data)));
                      list = (int *)Mem.myrealloc_movable(list, MaxTopNodes * sizeof(int));
                    }

                  if(topNodes[n].Level >= BITS_FOR_POSITIONS - 2)
                    {
                      if(message_printed == 0)
                        {
                          domain_printf(
                              "DOMAIN: Note: we would like to refine top-tree beyond the level allowed by the selected positional "
                              "accuracy.\n");
                          message_printed = 1;
                        }
                    }
                  else
                    {
                      list[count] = n;
                      count++;
                    }
                }
          }

      if(count > 0)
        {
          domain_do_local_refine(count, list);
          iter++;
        }
    }
  while(count > 0);

  Mem.myfree(list);
  Mem.myfree(mp);

  /* count the number of top leaves */
  NTopleaves = 0;
  domain_walktoptree(0);

  double t1 = Logs.second();

  domain_printf("DOMAIN: NTopleaves=%d, determination of top-level tree involved %d iterations and took %g sec\n", NTopleaves, iter,
                Logs.timediff(t0, t1));
}

#include "../data/simparticles.h"
template class domain<simparticles>;

#ifdef LIGHTCONE_PARTICLES
#include "../data/lcparticles.h"
template class domain<lcparticles>;
#endif
