/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  hydra.cc
 *
 *  \brief computation of SPH forces and rate of entropy generation
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../sort/cxxsort.h"
#include "../sph/kernel.h"
#include "../sph/sph.h"
#include "../system/system.h"

/*! This file contains the "second SPH loop", where the SPH forces are
 *  computed, and where the rate of change of entropy due to the shock heating
 *  (via artificial viscosity) is computed.
 */

inline int sph::sph_hydro_evaluate_particle_node_opening_criterion(pinfo &pdat, ngbnode *nop)
{
  if(nop->level <= LEVEL_ALWAYS_OPEN)  // always open the root node (note: full node length does not fit in the integer type)
    return NODE_OPEN;

  if(nop->Ti_Current != All.Ti_Current)
    nop->drift_node(All.Ti_Current, Tp);

  MyNgbTreeFloat dist = std::max<MyNgbTreeFloat>(nop->MaxHsml, pdat.hsml);

  MyIntPosType search_min[3], search_range[3];

  MyIntPosType inthsml = dist * Tp->FacCoordToInt;

  for(int i = 0; i < 3; i++)
    {
      search_min[i]   = pdat.searchcenter[i] - inthsml;
      search_range[i] = inthsml + inthsml;
    }

  MyIntPosType left[3], right[3];

  left[0]  = Tp->nearest_image_intpos_to_intpos_X(nop->center_offset_min[0] + nop->center[0], search_min[0]);
  right[0] = Tp->nearest_image_intpos_to_intpos_X(nop->center_offset_max[0] + nop->center[0], search_min[0]);

  /* check whether we can stop walking along this branch */
  if(left[0] > search_range[0] && right[0] > left[0])
    return NODE_DISCARD;

  left[1]  = Tp->nearest_image_intpos_to_intpos_Y(nop->center_offset_min[1] + nop->center[1], search_min[1]);
  right[1] = Tp->nearest_image_intpos_to_intpos_Y(nop->center_offset_max[1] + nop->center[1], search_min[1]);

  /* check whether we can stop walking along this branch */
  if(left[1] > search_range[1] && right[1] > left[1])
    return NODE_DISCARD;

  left[2]  = Tp->nearest_image_intpos_to_intpos_Z(nop->center_offset_min[2] + nop->center[2], search_min[2]);
  right[2] = Tp->nearest_image_intpos_to_intpos_Z(nop->center_offset_max[2] + nop->center[2], search_min[2]);

  /* check whether we can stop walking along this branch */
  if(left[2] > search_range[2] && right[2] > left[2])
    return NODE_DISCARD;

  return NODE_OPEN;
}

inline void sph::sph_hydro_check_particle_particle_interaction(pinfo &pdat, int p, int p_type, unsigned char shmrank)
{
#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  if(skip_actual_force_computation)
    return;
#endif

  if(p_type == NODE_TYPE_LOCAL_PARTICLE) /* local particle */
    {
      particle_data *P        = get_Pp(p, shmrank);
      sph_particle_data *SphP = get_SphPp(p, shmrank);

      if(P->getType() > 0)
        return;

      if(P->get_Ti_Current() != All.Ti_Current)
        Tp->drift_particle(P, SphP, All.Ti_Current);  // this function avoids race conditions

      MyNgbTreeFloat dist   = std::max<MyNgbTreeFloat>(SphP->Hsml, pdat.hsml);
      MyNgbTreeFloat distsq = dist * dist;

      double posdiff[3];
      Tp->nearest_image_intpos_to_pos(P->IntPos, pdat.searchcenter, posdiff); /* converts the integer distance to floating point */

      double rad2 = posdiff[0] * posdiff[0] + posdiff[1] * posdiff[1] + posdiff[2] * posdiff[2];
      if(rad2 > distsq || rad2 == 0)
        return;

      if(pdat.numngb >= MAX_NGBS)
        Terminate("pdat.numngb >= MAX_NGBS");

      int n = pdat.numngb++;

      Ngbhydrodat[n].SphCore = SphP;
      Ngbhydrodat[n].IntPos  = P->IntPos;
      Ngbhydrodat[n].Mass    = P->getMass();
#ifndef LEAN
      Ngbhydrodat[n].TimeBinHydro = P->TimeBinHydro;
#endif
    }
  else if(p_type == NODE_TYPE_FETCHED_PARTICLE)
    {
      foreign_sphpoint_data *foreignpoint = get_foreignpointsp(p - EndOfForeignNodes, shmrank);

      MyNgbTreeFloat dist   = std::max<MyNgbTreeFloat>(foreignpoint->SphCore.Hsml, pdat.hsml);
      MyNgbTreeFloat distsq = dist * dist;

      /* converts the integer distance to floating point */
      double posdiff[3];
      Tp->nearest_image_intpos_to_pos(foreignpoint->IntPos, pdat.searchcenter, posdiff);

      double rad2 = posdiff[0] * posdiff[0] + posdiff[1] * posdiff[1] + posdiff[2] * posdiff[2];
      if(rad2 > distsq || rad2 == 0)
        return;

      if(pdat.numngb >= MAX_NGBS)
        Terminate("pdat.numngb >= MAX_NGBS");

      int n = pdat.numngb++;

      Ngbhydrodat[n].SphCore      = &foreignpoint->SphCore;
      Ngbhydrodat[n].IntPos       = foreignpoint->IntPos;
      Ngbhydrodat[n].Mass         = foreignpoint->Mass;
      Ngbhydrodat[n].TimeBinHydro = foreignpoint->TimeBinHydro;
    }
  else
    Terminate("unexpected");
}

inline void sph::sph_hydro_open_node(pinfo &pdat, ngbnode *nop, int mintopleafnode, int committed)
{
  /* open node */
  int p                 = nop->nextnode;
  unsigned char shmrank = nop->nextnode_shmrank;

  while(p != nop->sibling || (shmrank != nop->sibling_shmrank && nop->sibling >= MaxPart + D->NTopnodes))
    {
      if(p < 0)
        Terminate(
            "p=%d < 0  nop->sibling=%d nop->nextnode=%d shmrank=%d nop->sibling_shmrank=%d nop->foreigntask=%d  "
            "first_nontoplevelnode=%d",
            p, nop->sibling, nop->nextnode, shmrank, nop->sibling_shmrank, nop->OriginTask, MaxPart + D->NTopnodes);

      int next;
      unsigned char next_shmrank;
      char type;

      if(p < MaxPart) /* a local particle */
        {
          /* note: here shmrank cannot change */
          next         = get_nextnodep(shmrank)[p];
          next_shmrank = shmrank;
          type         = NODE_TYPE_LOCAL_PARTICLE;
        }
      else if(p < MaxPart + MaxNodes) /* an internal node  */
        {
          ngbnode *nop = get_nodep(p, shmrank);
          next         = nop->sibling;
          next_shmrank = nop->sibling_shmrank;
          type         = NODE_TYPE_LOCAL_NODE;
        }
      else if(p >= ImportedNodeOffset && p < EndOfTreePoints) /* an imported Treepoint particle  */
        {
          Terminate("not expected for SPH");
        }
      else if(p >= EndOfTreePoints && p < EndOfForeignNodes) /* an imported tree node */
        {
          ngbnode *nop = get_nodep(p, shmrank);
          next         = nop->sibling;
          next_shmrank = nop->sibling_shmrank;
          type         = NODE_TYPE_FETCHED_NODE;
        }
      else if(p >= EndOfForeignNodes) /* an imported particle below an imported tree node */
        {
          foreign_sphpoint_data *foreignpoint = get_foreignpointsp(p - EndOfForeignNodes, shmrank);

          next         = foreignpoint->Nextnode;
          next_shmrank = foreignpoint->Nextnode_shmrank;
          type         = NODE_TYPE_FETCHED_PARTICLE;
        }
      else
        {
          /* a pseudo point */
          Terminate(
              "should not happen: p=%d MaxPart=%d MaxNodes=%d  ImportedNodeOffset=%d  EndOfTreePoints=%d  EndOfForeignNodes=%d "
              "shmrank=%d",
              p, MaxPart, MaxNodes, ImportedNodeOffset, EndOfTreePoints, EndOfForeignNodes, shmrank);
        }

      sph_hydro_interact(pdat, p, type, shmrank, mintopleafnode, committed);

      p       = next;
      shmrank = next_shmrank;
    }
}

inline void sph::sph_hydro_interact(pinfo &pdat, int no, char no_type, unsigned char shmrank, int mintopleafnode, int committed)
{
  if(no_type <= NODE_TYPE_FETCHED_PARTICLE)  // we are interacting with a particle
    {
      sph_hydro_check_particle_particle_interaction(pdat, no, no_type, shmrank);
    }
  else  // we are interacting with a node
    {
      ngbnode *nop = get_nodep(no, shmrank);

      if(nop->not_empty == 0)
        return;

      if(no < MaxPart + MaxNodes)                // we have a top-levelnode
        if(nop->nextnode >= MaxPart + MaxNodes)  // if the next node is not a top-level, we have a leaf node
          mintopleafnode = no;

      int openflag = sph_hydro_evaluate_particle_node_opening_criterion(pdat, nop);

      if(openflag == NODE_OPEN) /* we need to open it */
        {
          if(nop->cannot_be_opened_locally.load(std::memory_order_acquire))
            {
              // are we in the same shared memory node?
              if(Shmem.GetNodeIDForSimulCommRank[nop->OriginTask] == Shmem.GetNodeIDForSimulCommRank[D->ThisTask])
                {
                  Terminate("this should not happen any more");
                }
              else
                {
                  tree_add_to_fetch_stack(nop, no, shmrank);  // will only add unique copies

                  tree_add_to_work_stack(pdat.target, no, shmrank, mintopleafnode);
                }
            }
          else
            {
              int min_buffer_space =
                  std::min<int>(MaxOnWorkStack - (NumOnWorkStack + NewOnWorkStack), MaxOnFetchStack - NumOnFetchStack);

              if(min_buffer_space >= committed + 8 * TREE_NUM_BEFORE_NODESPLIT)
                sph_hydro_open_node(pdat, nop, mintopleafnode, committed + 8 * TREE_NUM_BEFORE_NODESPLIT);
              else
                tree_add_to_work_stack(pdat.target, no, shmrank, mintopleafnode);
            }
        }
    }
}

void sph::hydro_forces_determine(int ntarget, int *targetlist)
{
  TIMER_STORE;
  TIMER_START(CPU_HYDRO);

  D->mpi_printf("SPH-HYDRO: Begin hydro-force calculation.  (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());
  D->mpi_printf("SPH-HYDRO: global Nhydro=%llu (task zero: NumGas=%d, Nhydro=%d)\n", Tp->TimeBinsHydro.GlobalNActiveParticles,
                Tp->NumGas, ntarget);

  double ta = Logs.second();

  // let's grab at most half the still available memory for imported points and nodes
  int nspace = (0.33 * Mem.FreeBytes) / (sizeof(ngbnode) + 8 * sizeof(foreign_sphpoint_data));

  MaxForeignNodes  = nspace;
  MaxForeignPoints = 8 * nspace;
  NumForeignNodes  = 0;
  NumForeignPoints = 0;

  sum_NumForeignNodes  = 0;
  sum_NumForeignPoints = 0;

  /* the following two arrays hold imported tree nodes and imported points to augment the local tree */
  Foreign_Nodes  = (ngbnode *)Mem.mymalloc_movable(&Foreign_Nodes, "Foreign_Nodes", MaxForeignNodes * sizeof(ngbnode));
  Foreign_Points = (foreign_sphpoint_data *)Mem.mymalloc_movable(&Foreign_Points, "Foreign_Points",
                                                                 MaxForeignPoints * sizeof(foreign_sphpoint_data));

  tree_initialize_leaf_node_access_info();

  max_ncycles = 0;

  prepare_shared_memory_access();

  if(All.ComovingIntegrationOn)
    {
      fac_mu       = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;
      fac_vsic_fix = All.cf_hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);
    }
  else
    {
      fac_mu       = 1.0;
      fac_vsic_fix = 1.0;
    }

  Ngbhydrodat = (ngbdata_hydro *)Mem.mymalloc("Ngbhydrodat", MAX_NGBS * sizeof(ngbdata_hydro));

  NumOnWorkStack         = 0;
  AllocWorkStackBaseLow  = std::max<int>(1.5 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
  AllocWorkStackBaseHigh = AllocWorkStackBaseLow + TREE_EXPECTED_CYCLES * TREE_MIN_WORKSTACK_SIZE;
  MaxOnWorkStack         = AllocWorkStackBaseLow;

  WorkStack = (workstack_data *)Mem.mymalloc("WorkStack", AllocWorkStackBaseHigh * sizeof(workstack_data));

  for(int i = 0; i < ntarget; i++)
    {
      int target = targetlist[i];

      clear_hydro_result(&Tp->SphP[target]);

      WorkStack[NumOnWorkStack].Target         = target;
      WorkStack[NumOnWorkStack].Node           = MaxPart;
      WorkStack[NumOnWorkStack].ShmRank        = Shmem.Island_ThisTask;
      WorkStack[NumOnWorkStack].MinTopLeafNode = MaxPart + D->NTopnodes;
      NumOnWorkStack++;
    }

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  workstack_data *WorkStackBak = (workstack_data *)Mem.mymalloc("WorkStackBak", NumOnWorkStack * sizeof(workstack_data));
  int NumOnWorkStackBak        = NumOnWorkStack;
  memcpy(WorkStackBak, WorkStack, NumOnWorkStack * sizeof(workstack_data));
#endif

  // set a default size of the fetch stack equal to half the work stack (this may still be somewhat too large)
  MaxOnFetchStack = std::max<int>(0.1 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
  StackToFetch    = (fetch_data *)Mem.mymalloc_movable(&StackToFetch, "StackToFetch", MaxOnFetchStack * sizeof(fetch_data));

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  for(int rep = 0; rep < 2; rep++)
    {
      if(rep == 0)
        {
          skip_actual_force_computation = true;
        }
      else
        {
          skip_actual_force_computation = false;
          NumOnWorkStack                = NumOnWorkStackBak;
          memcpy(WorkStack, WorkStackBak, NumOnWorkStack * sizeof(workstack_data));
        }
#endif

      while(NumOnWorkStack > 0)  // repeat until we are out of work
        {
          NewOnWorkStack  = 0;  // gives the new entries
          NumOnFetchStack = 0;
          MaxOnWorkStack  = std::min<int>(AllocWorkStackBaseLow + max_ncycles * TREE_MIN_WORKSTACK_SIZE, AllocWorkStackBaseHigh);

          TIMER_START(CPU_HYDROWALK);

          int item = 0;

          while(item < NumOnWorkStack)
            {
              int committed = 8 * TREE_NUM_BEFORE_NODESPLIT;
              int min_buffer_space =
                  std::min<int>(MaxOnWorkStack - (NumOnWorkStack + NewOnWorkStack), MaxOnFetchStack - NumOnFetchStack);
              if(min_buffer_space >= committed)
                {
                  int target     = WorkStack[item].Target;
                  int no         = WorkStack[item].Node;
                  int shmrank    = WorkStack[item].ShmRank;
                  int mintopleaf = WorkStack[item].MinTopLeafNode;
                  item++;

                  pinfo pdat;
                  get_pinfo(target, pdat);

                  if(no == MaxPart)
                    {
                      // we have a pristine particle that's processed for the first time
                      sph_hydro_interact(pdat, no, NODE_TYPE_LOCAL_NODE, shmrank, mintopleaf, committed);
                    }
                  else
                    {
                      // we have a node that we previously could not open
                      ngbnode *nop = get_nodep(no, shmrank);

                      if(nop->cannot_be_opened_locally)
                        {
                          Terminate("item=%d:  no=%d  now we should be able to open it!", item, no);
                        }
                      else
                        sph_hydro_open_node(pdat, nop, mintopleaf, committed);
                    }

                  hydro_evaluate_kernel(pdat);
                }
              else
                break;
            }

          if(item == 0 && NumOnWorkStack > 0)
            Terminate("Can't even process a single particle");

          TIMER_STOP(CPU_HYDROWALK);

          TIMER_START(CPU_HYDROFETCH);

          tree_fetch_foreign_nodes(FETCH_SPH_HYDRO);

          TIMER_STOP(CPU_HYDROFETCH);

          /* now reorder the workstack such that we are first going to do residual pristine particles, and then
           * imported nodes that hang below the first leaf nodes */
          NumOnWorkStack = NumOnWorkStack - item + NewOnWorkStack;
          memmove(WorkStack, WorkStack + item, NumOnWorkStack * sizeof(workstack_data));

          /* now let's sort such that we can go deep on top-level node branches, allowing us to clear them out eventually */
          mycxxsort(WorkStack, WorkStack + NumOnWorkStack, compare_workstack);

          max_ncycles++;
        }

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
    }
#endif

  Mem.myfree(StackToFetch);
#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  Mem.myfree(WorkStackBak);
#endif
  Mem.myfree(WorkStack);
  Mem.myfree(Ngbhydrodat);

  /* now factor in a prefactor for the computed rates */
  for(int i = 0; i < ntarget; i++)
    {
      int target = targetlist[i];

      double fac = GAMMA_MINUS1 / (All.cf_atime2_hubble_a * pow(Tp->SphP[target].Density, GAMMA_MINUS1));

      Tp->SphP[target].DtEntropy *= fac;
    }

  /* Now the tree-based hydrodynamical force computation is finished,
   * output some performance metrics
   */

  TIMER_START(CPU_HYDROIMBALANCE);

  MPI_Allreduce(MPI_IN_PLACE, &max_ncycles, 1, MPI_INT, MPI_MAX, D->Communicator);

  TIMER_STOP(CPU_HYDROIMBALANCE);

  cleanup_shared_memory_access();

  /* free temporary buffers */
  Mem.myfree(Foreign_Points);
  Mem.myfree(Foreign_Nodes);

  double tb = Logs.second();

  TIMER_STOPSTART(CPU_HYDRO, CPU_LOGS);

  D->mpi_printf("SPH-HYDRO: hydro-force computation done. took %8.3f\n", Logs.timediff(ta, tb));

  struct detailed_timings
  {
    double tree, wait, fetch, all;
    double numnodes;
    double NumForeignNodes, NumForeignPoints;
    double fillfacFgnNodes, fillfacFgnPoints;
  };
  detailed_timings timer, tisum, timax;

  timer.tree             = TIMER_DIFF(CPU_HYDROWALK);
  timer.wait             = TIMER_DIFF(CPU_HYDROIMBALANCE);
  timer.fetch            = TIMER_DIFF(CPU_HYDROFETCH);
  timer.all              = timer.tree + timer.wait + timer.fetch + TIMER_DIFF(CPU_HYDRO);
  timer.numnodes         = NumNodes;
  timer.NumForeignNodes  = NumForeignNodes;
  timer.NumForeignPoints = NumForeignPoints;
  timer.fillfacFgnNodes  = NumForeignNodes / ((double)MaxForeignNodes);
  timer.fillfacFgnPoints = NumForeignPoints / ((double)MaxForeignPoints);

  MPI_Reduce((double *)&timer, (double *)&tisum, (int)(sizeof(detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_SUM, 0,
             D->Communicator);
  MPI_Reduce((double *)&timer, (double *)&timax, (int)(sizeof(detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_MAX, 0,
             D->Communicator);

  All.TotNumHydro += Tp->TimeBinsHydro.GlobalNActiveParticles;

  if(D->ThisTask == 0)
    {
      fprintf(Logs.FdHydro, "Nf=%9lld  highest active timebin=%d  total-Nf=%lld\n", Tp->TimeBinsHydro.GlobalNActiveParticles,
              All.HighestActiveTimeBin, All.TotNumHydro);
      fprintf(Logs.FdHydro, "   work-load balance: %g   part/sec: raw=%g, effective=%g\n",
              timax.tree / ((tisum.tree + 1e-20) / D->NTask), Tp->TimeBinsGravity.GlobalNActiveParticles / (tisum.tree + 1.0e-20),
              Tp->TimeBinsGravity.GlobalNActiveParticles / ((timax.tree + 1.0e-20) * D->NTask));
      fprintf(Logs.FdHydro,
              "   maximum number of nodes: %g, filled: %g  NumForeignNodes: max=%g avg=%g fill=%g NumForeignPoints: max=%g avg=%g "
              "fill=%g  cycles=%d\n",
              timax.numnodes, timax.numnodes / MaxNodes, timax.NumForeignNodes, tisum.NumForeignNodes / D->NTask,
              timax.fillfacFgnNodes, timax.NumForeignPoints, tisum.NumForeignPoints / D->NTask, timax.fillfacFgnPoints, max_ncycles);
      fprintf(Logs.FdHydro, "   avg times: <all>=%g  <tree>=%g  <wait>=%g  <fetch>=%g  sec\n", tisum.all / D->NTask,
              tisum.tree / D->NTask, tisum.wait / D->NTask, tisum.fetch / D->NTask);
      myflush(Logs.FdHydro);
    }

  TIMER_STOP(CPU_LOGS);
}

#ifdef EXPLICIT_VECTORIZATION
void sph::hydro_evaluate_kernel(pinfo &pdat)
{
#ifndef LEAN
  particle_data *P_i        = &Tp->P[pdat.target];
  sph_particle_data *SphP_i = &Tp->SphP[pdat.target];

  /* the particles needs to be active */
  if(P_i->getTimeBinHydro() > All.HighestSynchronizedTimeBin)
    Terminate("bummer");

  double shinv, shinv3, shinv4;
  kernel_hinv(SphP_i->Hsml, &shinv, &shinv3, &shinv4);

  Vec4d hinv(shinv);
  Vec4d hinv3(shinv3);
  Vec4d hinv4(shinv4);

  Vec4d dwnorm(NORM * shinv3);
  Vec4d dwknorm(NORM * shinv4);

  Vec4d rho_i(SphP_i->Density);

#ifdef PRESSURE_ENTROPY_SPH
  Vec4d p_over_rho2_i((double)SphP_i->Pressure / ((double)SphP_i->PressureSphDensity * (double)SphP_i->PressureSphDensity));
#else
  Vec4d p_over_rho2_i((double)SphP_i->Pressure / ((double)SphP_i->Density * (double)SphP_i->Density));
#endif

  Vec4d sound_i(SphP_i->Csnd);
  Vec4d h_i(SphP_i->Hsml);

  Vec4d v_i[3];
  for(int i = 0; i < NUMDIMS; i++)
    {
      v_i[i] = SphP_i->VelPred[i];
    }
  Vec4d DhsmlDensityFactor_i(SphP_i->DhsmlDensityFactor);
#ifdef PRESSURE_ENTROPY_SPH
  Vec4d DhsmlDerivedDensityFactor_i(SphP_i->DhsmlDerivedDensityFactor);
  Vec4d EntropyToInvGammaPred_i(SphP_i->EntropyToInvGammaPred);
#endif

#if !defined(NO_SHEAR_VISCOSITY_LIMITER) && !defined(TIMEDEP_ART_VISC)
  Vec4d f_i(fabs(SphP_i->DivVel) / (fabs(SphP_i->DivVel) + SphP_i->CurlVel + 0.0001 * SphP_i->Csnd / SphP_i->Hsml / fac_mu));
#endif

#ifdef TIMEDEP_ART_VISC
  Vec4d alpha_i(SphP_i->Alpha);
#endif
  /* Now start the actual SPH computation for this particle */

  double dacc[3]     = {0};
  double dentr       = 0;
  Vec4d MaxSignalVel = sound_i;

  const int vector_length = 4;
  const int array_length  = (pdat.numngb + vector_length - 1) & (-vector_length);

  for(int n = pdat.numngb; n < array_length; n++) /* fill up neighbour array so that sensible data is accessed */
    Ngbhydrodat[n] = Ngbhydrodat[0];

  for(int n = 0; n < array_length; n += vector_length)
    {
      sph_particle_data_hydrocore *ngb0 = Ngbhydrodat[n + 0].SphCore;
      sph_particle_data_hydrocore *ngb1 = Ngbhydrodat[n + 1].SphCore;
      sph_particle_data_hydrocore *ngb2 = Ngbhydrodat[n + 2].SphCore;
      sph_particle_data_hydrocore *ngb3 = Ngbhydrodat[n + 3].SphCore;

      ngbdata_hydro *P0_j = &Ngbhydrodat[n + 0];
      ngbdata_hydro *P1_j = &Ngbhydrodat[n + 1];
      ngbdata_hydro *P2_j = &Ngbhydrodat[n + 2];
      ngbdata_hydro *P3_j = &Ngbhydrodat[n + 3];

      /* converts the integer distance to floating point */
      Vec4d dpos[NUMDIMS];
      double posdiff[array_length][3];
      for(int i = 0; i < 4; i++)
        {
          Tp->nearest_image_intpos_to_pos(P_i->IntPos, Ngbhydrodat[n + i].IntPos, &(posdiff[i][0]));
        }

      for(int i = 0; i < NUMDIMS; i++)
        {
          dpos[i] = Vec4d(posdiff[0][i], posdiff[1][i], posdiff[2][i], posdiff[3][i]);
        }

      Vec4d r2(0);

      for(int i = 0; i < NUMDIMS; i++)
        {
          r2 += dpos[i] * dpos[i];
        }

      Vec4d r = sqrt(r2);

      Vec4d v_j[NUMDIMS];
      for(int i = 0; i < NUMDIMS; i++)
        {
          v_j[i] = Vec4d(ngb0->VelPred[i], ngb1->VelPred[i], ngb2->VelPred[i], ngb3->VelPred[i]);
        }

      Vec4d pressure(ngb0->Pressure, ngb1->Pressure, ngb2->Pressure, ngb3->Pressure);
      Vec4d rho_j(ngb0->Density, ngb1->Density, ngb2->Density, ngb3->Density);
#ifdef PRESSURE_ENTROPY_SPH
      Vec4d rho_press_j(ngb0->PressureSphDensity, ngb1->PressureSphDensity, ngb2->PressureSphDensity, ngb3->PressureSphDensity);
      Vec4d p_over_rho2_j = pressure / (rho_press_j * rho_press_j);
#else
      Vec4d p_over_rho2_j = pressure / (rho_j * rho_j);
#endif

      Vec4d wk_i, dwk_i;
      Vec4d u = r * hinv;
      kernel_main_vector(u, dwnorm, dwknorm, &wk_i, &dwk_i);
      Vec4db decision = (r < h_i);
      Vec4d fac       = select(decision, 1., 0.);
      wk_i *= fac;
      dwk_i *= fac;

      Vec4d h_j(ngb0->Hsml, ngb1->Hsml, ngb2->Hsml, ngb3->Hsml);
      Vec4d hinv_j = 1 / h_j;
#ifdef THREEDIMS
      Vec4d hinv3_j = hinv_j * hinv_j * hinv_j;
#endif

#ifdef TWODIMS
      Vec4d hinv3_j = hinv_j * hinv_j;
#endif

#ifdef ONEDIMS
      Vec4d hinv3_j = hinv_j;
#endif
      Vec4d hinv4_j = hinv3_j * hinv_j;

      Vec4d wk_j, dwk_j;
      u = r * hinv_j;
      kernel_main_vector(u, NORM * hinv3_j, NORM * hinv4_j, &wk_j, &dwk_j);
      decision = (r < h_j);
      fac      = select(decision, 1., 0.);
      wk_j *= fac;
      dwk_j *= fac;

      Vec4d sound_j(ngb0->Csnd, ngb1->Csnd, ngb2->Csnd, ngb3->Csnd);
      Vec4d vsig = sound_i + sound_j;
      if(n + vector_length > pdat.numngb)
        {
          wk_i.cutoff(vector_length - (array_length - pdat.numngb));
          dwk_i.cutoff(vector_length - (array_length - pdat.numngb));
          wk_j.cutoff(vector_length - (array_length - pdat.numngb));
          dwk_j.cutoff(vector_length - (array_length - pdat.numngb));
          vsig.cutoff(vector_length - (array_length - pdat.numngb));
        }

      Vec4d dwk_ij = 0.5 * (dwk_i + dwk_j);

      MaxSignalVel = max(MaxSignalVel, vsig);

      Vec4d visc(0);

      Vec4d dv[NUMDIMS];
      for(int i = 0; i < NUMDIMS; i++)
        {
          dv[i] = v_i[i] - v_j[i];
        }

      Vec4d vdotr2(0);
      for(int i = 0; i < NUMDIMS; i++)
        {
          vdotr2 += dv[i] * dpos[i];
        }

      if(All.ComovingIntegrationOn)
        vdotr2 += All.cf_atime2_hubble_a * r2;

      decision            = (vdotr2 < 0);
      Vec4d viscosity_fac = select(decision, 1, 0);

      /* ... artificial viscosity */

      Vec4d mu_ij = fac_mu * vdotr2 / r;

      vsig -= 3 * mu_ij;

#if defined(NO_SHEAR_VISCOSITY_LIMITER) || defined(TIMEDEP_ART_VISC)
      Vec4d f_i(1);
      Vec4d f_j(1);
#else
      Vec4d DivVel_j(ngb0->DivVel, ngb1->DivVel, ngb2->DivVel, ngb3->DivVel);
      Vec4d CurlVel_j(ngb0->CurlVel, ngb1->CurlVel, ngb2->CurlVel, ngb3->CurlVel);
      Vec4d f_j = abs(DivVel_j) / (abs(DivVel_j) + CurlVel_j + 0.0001 * sound_j / fac_mu * hinv_j);
#endif

#ifdef TIMEDEP_ART_VISC
      Vec4d alpha_j(ngb0->Alpha, ngb1->Alpha, ngb2->Alpha, ngb3->Alpha);
      Vec4d BulkVisc_ij = 0.5 * (alpha_i + alpha_j);

#else
      Vec4d BulkVisc_ij(All.ArtBulkViscConst);
#endif
      Vec4d rho_ij_inv = 2.0 / (rho_i + rho_j);
      visc             = 0.25 * BulkVisc_ij * vsig * (-mu_ij) * rho_ij_inv * (f_i + f_j);
      Vec4d mass_j(P0_j->Mass, P1_j->Mass, P2_j->Mass, P3_j->Mass);
#ifdef VISCOSITY_LIMITER_FOR_LARGE_TIMESTEPS
      Vec4i timeBin_i(P_i->TimeBinHydro);
      Vec4i timeBin_j(P0_j->TimeBinHydro, P1_j->TimeBinHydro, P2_j->TimeBinHydro, P3_j->TimeBinHydro);

      Vec4i timebin = max(timeBin_i, timeBin_j);
      Vec4i integer_time(((integertime)1) << timebin[0], ((integertime)1) << timebin[1], ((integertime)1) << timebin[2],
                         ((integertime)1) << timebin[3]);

      Vec4ib decision_i    = (timebin != 0);
      Vec4i factor_timebin = select(decision_i, Vec4i(1), Vec4i(0));
      Vec4d dt             = to_double(2 * integer_time * factor_timebin) * All.Timebase_interval;

      decision = (dt > 0 && dwk_ij < 0);

      Vec4d visc_alternavtive = 0.5 * fac_vsic_fix * vdotr2 / ((P_i->getMass() + mass_j) * dwk_ij * r * dt);

      Vec4d visc2 = select(decision, visc_alternavtive, visc);
      visc        = min(visc, visc2);
#endif

      Vec4d hfc_visc = mass_j * visc * dwk_ij / r * viscosity_fac;

#ifndef PRESSURE_ENTROPY_SPH
      /* Formulation derived from the Lagrangian */
      dwk_i *= DhsmlDensityFactor_i;
      Vec4d DhsmlDensityFactor_j(ngb0->DhsmlDensityFactor, ngb1->DhsmlDensityFactor, ngb2->DhsmlDensityFactor,
                                 ngb3->DhsmlDensityFactor);
      dwk_j *= DhsmlDensityFactor_j;

      Vec4d hfc = mass_j * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r + hfc_visc;
#else
      Vec4d EntropyToInvGammaPred_j(ngb0->EntropyToInvGammaPred, ngb1->EntropyToInvGammaPred, ngb2->EntropyToInvGammaPred,
                                    ngb3->EntropyToInvGammaPred);
      Vec4d DhsmlDerivedDensityFactor_j(ngb0->DhsmlDerivedDensityFactor, ngb1->DhsmlDerivedDensityFactor,
                                        ngb2->DhsmlDerivedDensityFactor, ngb3->DhsmlDerivedDensityFactor);
      /* leading order term */
      Vec4d hfc = mass_j *
                  (p_over_rho2_i * dwk_i * EntropyToInvGammaPred_j / EntropyToInvGammaPred_i +
                   p_over_rho2_j * dwk_j * EntropyToInvGammaPred_i / EntropyToInvGammaPred_j) /
                  r;

      /* grad-h term */
      hfc += mass_j *
             (p_over_rho2_i * dwk_i * SphP_i->DhsmlDerivedDensityFactor + p_over_rho2_j * dwk_j * DhsmlDerivedDensityFactor_j) / r;

      /* add viscous term */
      hfc += hfc_visc;
#endif

      for(int i = 0; i < NUMDIMS; i++)
        {
          dacc[i] += horizontal_add(-hfc * dpos[i]);
        }
      dentr += horizontal_add(0.5 * (hfc_visc)*vdotr2);
    }

  SphP_i->HydroAccel[0] += dacc[0];
  SphP_i->HydroAccel[1] += dacc[1];
  SphP_i->HydroAccel[2] += dacc[2];
  SphP_i->DtEntropy += dentr;

  for(int i = 0; i < 4; i++)
    {
      if(SphP_i->MaxSignalVel < MaxSignalVel[i])
        SphP_i->MaxSignalVel = MaxSignalVel[i];
    }
#endif
}

#else

/*! This function is the 'core' of the SPH force computation. A target
 *  particle is specified which may either be local, or reside in the
 *  communication buffer.
 */
void sph::hydro_evaluate_kernel(pinfo &pdat)
{
#ifndef LEAN
  particle_data *P_i        = &Tp->P[pdat.target];
  sph_particle_data *SphP_i = &Tp->SphP[pdat.target];

  /* the particles needs to be active */
  if(P_i->getTimeBinHydro() > All.HighestSynchronizedTimeBin)
    Terminate("bummer");

#ifdef PRESSURE_ENTROPY_SPH
  double p_over_rho2_i = (double)SphP_i->Pressure / ((double)SphP_i->PressureSphDensity * (double)SphP_i->PressureSphDensity);
#else
  double p_over_rho2_i = (double)SphP_i->Pressure / ((double)SphP_i->Density * (double)SphP_i->Density);
#endif

  kernel_hydra kernel;

  kernel.sound_i = SphP_i->Csnd;
  kernel.h_i     = SphP_i->Hsml;

  /* Now start the actual SPH computation for this particle */

  double daccx        = 0;
  double daccy        = 0;
  double daccz        = 0;
  double dentr        = 0;
  double MaxSignalVel = kernel.sound_i;

  for(int n = 0; n < pdat.numngb; n++)
    {
      sph_particle_data_hydrocore *SphP_j = Ngbhydrodat[n].SphCore;
      ngbdata_hydro *P_j                  = &Ngbhydrodat[n];

      /* converts the integer distance to floating point */
      double posdiff[3];
      Tp->nearest_image_intpos_to_pos(P_i->IntPos, P_j->IntPos, posdiff);

      kernel.dx = posdiff[0];
      kernel.dy = posdiff[1];
      kernel.dz = posdiff[2];

      double r2  = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;
      kernel.h_j = SphP_j->Hsml;

      if(r2 < kernel.h_i * kernel.h_i || r2 < kernel.h_j * kernel.h_j)
        {
          kernel.r = sqrt(r2);
          if(kernel.r > 0)
            {
#ifdef PRESSURE_ENTROPY_SPH
              double p_over_rho2_j =
                  (double)SphP_j->Pressure / ((double)SphP_j->PressureSphDensity * (double)SphP_j->PressureSphDensity);
#else
              double p_over_rho2_j = (double)SphP_j->Pressure / ((double)SphP_j->Density * (double)SphP_j->Density);
#endif

              kernel.sound_j = SphP_j->Csnd;

              kernel.dvx        = SphP_i->VelPred[0] - SphP_j->VelPred[0];
              kernel.dvy        = SphP_i->VelPred[1] - SphP_j->VelPred[1];
              kernel.dvz        = SphP_i->VelPred[2] - SphP_j->VelPred[2];
              kernel.vdotr2     = kernel.dx * kernel.dvx + kernel.dy * kernel.dvy + kernel.dz * kernel.dvz;
              kernel.rho_ij_inv = 2.0 / (SphP_i->Density + SphP_j->Density);

              if(All.ComovingIntegrationOn)
                kernel.vdotr2 += All.cf_atime2_hubble_a * r2;

              double hinv, hinv3, hinv4;
              if(kernel.r < kernel.h_i)
                {
                  kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
                  double u = kernel.r * hinv;
                  kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, COMPUTE_DWK);
                }
              else
                {
                  kernel.dwk_i = 0;
                  kernel.wk_i  = 0;
                }

              if(kernel.r < kernel.h_j)
                {
                  kernel_hinv(kernel.h_j, &hinv, &hinv3, &hinv4);
                  double u = kernel.r * hinv;
                  kernel_main(u, hinv3, hinv4, &kernel.wk_j, &kernel.dwk_j, COMPUTE_DWK);
                }
              else
                {
                  kernel.dwk_j = 0;
                  kernel.wk_j  = 0;
                }

              kernel.dwk_ij = 0.5 * (kernel.dwk_i + kernel.dwk_j);

              kernel.vsig = kernel.sound_i + kernel.sound_j;

              if(kernel.vsig > MaxSignalVel)
                MaxSignalVel = kernel.vsig;

              double visc = 0;

              if(kernel.vdotr2 < 0) /* ... artificial viscosity */
                {
                  double mu_ij = fac_mu * kernel.vdotr2 / kernel.r;

                  kernel.vsig -= 3 * mu_ij;

#if defined(NO_SHEAR_VISCOSITY_LIMITER) || defined(TIMEDEP_ART_VISC)
                  double f_i         = 1.;
                  double f_j         = 1.;
#else
                  double f_i =
                      fabs(SphP_i->DivVel) / (fabs(SphP_i->DivVel) + SphP_i->CurlVel + 0.0001 * SphP_i->Csnd / SphP_i->Hsml / fac_mu);

                  double f_j =
                      fabs(SphP_j->DivVel) / (fabs(SphP_j->DivVel) + SphP_j->CurlVel + 0.0001 * kernel.sound_j / fac_mu / kernel.h_j);
#endif

#ifdef TIMEDEP_ART_VISC
                  double BulkVisc_ij = 0.5 * (SphP_i->Alpha + SphP_j->Alpha);

#else
                  double BulkVisc_ij = All.ArtBulkViscConst;
#endif

                  visc        = 0.25 * BulkVisc_ij * kernel.vsig * (-mu_ij) * kernel.rho_ij_inv * (f_i + f_j);
#ifdef VISCOSITY_LIMITER_FOR_LARGE_TIMESTEPS
                  int timebin = std::max<int>(P_i->TimeBinHydro, P_j->TimeBinHydro);

                  double dt = 2 * (timebin ? (((integertime)1) << timebin) : 0) * All.Timebase_interval;

                  if(dt > 0 && kernel.dwk_ij < 0)
                    {
                      visc = std::min<double>(
                          visc, 0.5 * fac_vsic_fix * kernel.vdotr2 / ((P_i->getMass() + P_j->Mass) * kernel.dwk_ij * kernel.r * dt));
                    }
#endif
                }

              double hfc_visc = P_j->Mass * visc * kernel.dwk_ij / kernel.r;

#ifndef PRESSURE_ENTROPY_SPH
              /* Formulation derived from the Lagrangian */
              kernel.dwk_i *= SphP_i->DhsmlDensityFactor;
              kernel.dwk_j *= SphP_j->DhsmlDensityFactor;

              double hfc = P_j->Mass * (p_over_rho2_i * kernel.dwk_i + p_over_rho2_j * kernel.dwk_j) / kernel.r + hfc_visc;
#else
              /* leading order term */
              double hfc = P_j->Mass *
                           (p_over_rho2_i * kernel.dwk_i * SphP_j->EntropyToInvGammaPred / SphP_i->EntropyToInvGammaPred +
                            p_over_rho2_j * kernel.dwk_j * SphP_i->EntropyToInvGammaPred / SphP_j->EntropyToInvGammaPred) /
                           kernel.r;

              /* grad-h term */
              hfc += P_j->Mass *
                     (p_over_rho2_i * kernel.dwk_i * SphP_i->DhsmlDerivedDensityFactor +
                      p_over_rho2_j * kernel.dwk_j * SphP_j->DhsmlDerivedDensityFactor) /
                     kernel.r;

              /* add viscous term */
              hfc += hfc_visc;
#endif

              daccx += (-hfc * kernel.dx);
              daccy += (-hfc * kernel.dy);
              daccz += (-hfc * kernel.dz);
              dentr += (0.5 * (hfc_visc)*kernel.vdotr2);
            }
        }
    }

  SphP_i->HydroAccel[0] += daccx;
  SphP_i->HydroAccel[1] += daccy;
  SphP_i->HydroAccel[2] += daccz;
  SphP_i->DtEntropy += dentr;

  if(SphP_i->MaxSignalVel < MaxSignalVel)
    SphP_i->MaxSignalVel = MaxSignalVel;
#endif
}
#endif

/* this routine clears the fields in the SphP particle structure that are additively computed by the SPH density loop
 * by summing over neighbours
 */
inline void sph::clear_hydro_result(sph_particle_data *SphP)
{
  for(int k = 0; k < 3; k++)
    SphP->HydroAccel[k] = 0;

  SphP->DtEntropy    = 0;
  SphP->MaxSignalVel = 0;
}
