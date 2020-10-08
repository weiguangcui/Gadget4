/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file timestep_treebased.cc
 *
 *  \brief routines to find the timestep by checking for the arrival of the first waves from anywhere
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../data/simparticles.h"
#include "../domain/domain.h"
#include "../gravtree/gravtree.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

inline int sph::sph_treetimestep_evaluate_particle_node_opening_criterion(pinfo &pdat, ngbnode *nop)
{
  if(nop->level <= LEVEL_ALWAYS_OPEN)  // always open the root node (note: full node length does not fit in the integer type)
    return NODE_OPEN;

  if(nop->Ti_Current != All.Ti_Current)
    nop->drift_node(All.Ti_Current, Tp);

  sph_particle_data *SphP = &Tp->SphP[pdat.target];

  double vsig = SphP->Csnd + nop->MaxCsnd;

  MyIntPosType new_range_min[3], new_range_max[3];
  MySignedIntPosType left[3], right[3];
  double vleft, vright;

  // ----------- x-checks

  vleft = (-vsig + (nop->vmin[0] - SphP->VelPred[0]));
  if(fabs(vleft) * SphP->CurrentMaxTiStep > MaxBoxDist)
    SphP->CurrentMaxTiStep = MaxBoxDist / fabs(vleft);

  vright = (vsig + (nop->vmax[0] - SphP->VelPred[0]));
  if(fabs(vright) * SphP->CurrentMaxTiStep > MaxBoxDist)
    SphP->CurrentMaxTiStep = MaxBoxDist / fabs(vright);

  new_range_min[0] = nop->center_offset_min[0] + nop->center[0] +
                     (MyIntPosType)Tp->pos_to_signedintpos(SphP->CurrentMaxTiStep * vleft + 2 * pdat.hsml);
  new_range_max[0] = nop->center_offset_max[0] + nop->center[0] +
                     (MyIntPosType)Tp->pos_to_signedintpos(SphP->CurrentMaxTiStep * vright - 2 * pdat.hsml);

  left[0]  = (MySignedIntPosType)Tp->nearest_image_intpos_to_intpos_X(new_range_min[0], pdat.searchcenter[0]);
  right[0] = (MySignedIntPosType)Tp->nearest_image_intpos_to_intpos_X(new_range_max[0], pdat.searchcenter[0]);

  /* check whether we can stop walking along this branch */
  if(right[0] < 0 || left[0] > 0)
    return NODE_DISCARD;

  // ----------- y-checks

  vleft = (-vsig + (nop->vmin[1] - SphP->VelPred[1]));
  if(fabs(vleft) * SphP->CurrentMaxTiStep > MaxBoxDist)
    SphP->CurrentMaxTiStep = MaxBoxDist / fabs(vleft);

  vright = (vsig + (nop->vmax[1] - SphP->VelPred[1]));
  if(fabs(vright) * SphP->CurrentMaxTiStep > MaxBoxDist)
    SphP->CurrentMaxTiStep = MaxBoxDist / fabs(vright);

  new_range_min[1] = nop->center_offset_min[1] + nop->center[1] +
                     (MyIntPosType)Tp->pos_to_signedintpos(SphP->CurrentMaxTiStep * vleft + 2 * pdat.hsml);
  new_range_max[1] = nop->center_offset_max[1] + nop->center[1] +
                     (MyIntPosType)Tp->pos_to_signedintpos(SphP->CurrentMaxTiStep * vright - 2 * pdat.hsml);

  left[1]  = (MySignedIntPosType)Tp->nearest_image_intpos_to_intpos_Y(new_range_min[1], pdat.searchcenter[1]);
  right[1] = (MySignedIntPosType)Tp->nearest_image_intpos_to_intpos_Y(new_range_max[1], pdat.searchcenter[1]);

  /* check whether we can stop walking along this branch */
  if(right[1] < 0 || left[1] > 0)
    return NODE_DISCARD;

  // ----------- z-checks

  vleft = (-vsig + (nop->vmin[2] - SphP->VelPred[2]));
  if(fabs(vleft) * SphP->CurrentMaxTiStep > MaxBoxDist)
    SphP->CurrentMaxTiStep = MaxBoxDist / fabs(vleft);

  vright = (vsig + (nop->vmax[2] - SphP->VelPred[2]));
  if(fabs(vright) * SphP->CurrentMaxTiStep > MaxBoxDist)
    SphP->CurrentMaxTiStep = MaxBoxDist / fabs(vright);

  new_range_min[2] = nop->center_offset_min[2] + nop->center[2] +
                     (MyIntPosType)Tp->pos_to_signedintpos(SphP->CurrentMaxTiStep * vleft + 2 * pdat.hsml);
  new_range_max[2] = nop->center_offset_max[2] + nop->center[2] +
                     (MyIntPosType)Tp->pos_to_signedintpos(SphP->CurrentMaxTiStep * vright - 2 * pdat.hsml);

  left[2]  = (MySignedIntPosType)Tp->nearest_image_intpos_to_intpos_Z(new_range_min[2], pdat.searchcenter[2]);
  right[2] = (MySignedIntPosType)Tp->nearest_image_intpos_to_intpos_Z(new_range_max[2], pdat.searchcenter[2]);

  /* check whether we can stop walking along this branch */
  if(right[2] < 0 || left[2] > 0)
    return NODE_DISCARD;

  return NODE_OPEN;
}

inline void sph::sph_treetimestep_check_particle_particle_interaction(pinfo &pdat, int p, int p_type, unsigned char shmrank)
{
#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  if(skip_actual_force_computation)
    return;
#endif

  sph_particle_data *SphP_i = &Tp->SphP[pdat.target];

  if(p_type == NODE_TYPE_LOCAL_PARTICLE) /* local particle */
    {
      particle_data *P_j        = get_Pp(p, shmrank);
      sph_particle_data *SphP_j = get_SphPp(p, shmrank);

      if(P_j->getType() > 0)
        return;

      if(P_j->get_Ti_Current() != All.Ti_Current)
        Tp->drift_particle(P_j, SphP_j, All.Ti_Current);  // this function avoids race conditions

      double dxyz[3];
      Tp->nearest_image_intpos_to_pos(P_j->IntPos, pdat.searchcenter, dxyz); /* converts the integer distance to floating point */

      double dist2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];

      if(dist2 > 0)
        {
          double dist = sqrt(dist2);
          double vsig = SphP_i->Csnd + SphP_j->Csnd -
                        ((SphP_j->VelPred[0] - SphP_i->VelPred[0]) * dxyz[0] + (SphP_j->VelPred[1] - SphP_i->VelPred[1]) * dxyz[1] +
                         (SphP_j->VelPred[2] - SphP_i->VelPred[2]) * dxyz[2]) /
                            dist;

          if(vsig > 0)
            {
              dist += 2 * SphP_i->Hsml; /* take one smoothing length as minimum distance in order to protect against unreasonably
                              small steps if two particles are very close */

              double dt = dist / vsig;
              if(SphP_i->CurrentMaxTiStep > dt)
                SphP_i->CurrentMaxTiStep = dt;
            }
        }
    }
  else if(p_type == NODE_TYPE_FETCHED_PARTICLE)
    {
      foreign_sphpoint_data *foreignpoint = get_foreignpointsp(p - EndOfForeignNodes, shmrank);

      sph_particle_data_hydrocore *SphP_j = &foreignpoint->SphCore;

      /* converts the integer distance to floating point */
      double dxyz[3];
      Tp->nearest_image_intpos_to_pos(foreignpoint->IntPos, pdat.searchcenter, dxyz);

      double dist2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];

      if(dist2 > 0)
        {
          double dist = sqrt(dist2);
          double vsig = SphP_i->Csnd + SphP_j->Csnd -
                        ((SphP_j->VelPred[0] - SphP_i->VelPred[0]) * dxyz[0] + (SphP_j->VelPred[1] - SphP_i->VelPred[1]) * dxyz[1] +
                         (SphP_j->VelPred[2] - SphP_i->VelPred[2]) * dxyz[2]) /
                            dist;

          if(vsig > 0)
            {
              dist += 2 * SphP_i->Hsml; /* take one smoothing length as minimum distance in order to protect against unreasonably
                              small steps if two particles are very close */

              double dt = dist / vsig;
              if(SphP_i->CurrentMaxTiStep > dt)
                SphP_i->CurrentMaxTiStep = dt;
            }
        }
    }
  else
    Terminate("unexpected");
}

inline void sph::sph_treetimestep_open_node(pinfo &pdat, ngbnode *nop, int mintopleafnode, int committed)
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

      sph_treetimestep_interact(pdat, p, type, shmrank, mintopleafnode, committed);

      p       = next;
      shmrank = next_shmrank;
    }
}

inline void sph::sph_treetimestep_interact(pinfo &pdat, int no, char no_type, unsigned char shmrank, int mintopleafnode, int committed)
{
  if(no_type <= NODE_TYPE_FETCHED_PARTICLE)  // we are interacting with a particle
    {
      sph_treetimestep_check_particle_particle_interaction(pdat, no, no_type, shmrank);
    }
  else  // we are interacting with a node
    {
      ngbnode *nop = get_nodep(no, shmrank);

      if(nop->not_empty == 0)
        return;

      if(no < MaxPart + MaxNodes)                // we have a top-levelnode
        if(nop->nextnode >= MaxPart + MaxNodes)  // if the next node is not a top-level, we have a leaf node
          mintopleafnode = no;

      int openflag = sph_treetimestep_evaluate_particle_node_opening_criterion(pdat, nop);

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
                sph_treetimestep_open_node(pdat, nop, mintopleafnode, committed + 8 * TREE_NUM_BEFORE_NODESPLIT);
              else
                tree_add_to_work_stack(pdat.target, no, shmrank, mintopleafnode);
            }
        }
    }
}

void sph::tree_based_timesteps(void)
{
  if(Tp->TotNumGas > 0)
    {
      TIMER_START(CPU_TREE_TIMESTEPS);

      D->mpi_printf("TIMESTEP-TREEWALK: Begin\n");

      MaxBoxDist = 0.25 * Tp->RegionLen / (1 << LEVEL_ALWAYS_OPEN);

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

      NumOnWorkStack         = 0;
      AllocWorkStackBaseLow  = std::max<int>(1.5 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
      AllocWorkStackBaseHigh = AllocWorkStackBaseLow + TREE_EXPECTED_CYCLES * TREE_MIN_WORKSTACK_SIZE;
      MaxOnWorkStack         = AllocWorkStackBaseLow;

      WorkStack = (workstack_data *)Mem.mymalloc("WorkStack", AllocWorkStackBaseHigh * sizeof(workstack_data));

      for(int i = 0; i < Tp->TimeBinsHydro.NActiveParticles; i++)
        {
          int target = Tp->TimeBinsHydro.ActiveParticleList[i];

          if(Tp->SphP[target].Csnd < MIN_FLOAT_NUMBER)
            Tp->SphP[target].Csnd = MIN_FLOAT_NUMBER;

          Tp->SphP[target].CurrentMaxTiStep = 2.0 * Tp->SphP[target].Hsml / (Tp->SphP[target].MaxSignalVel + MIN_FLOAT_NUMBER);

          /* note: for cosmological integration, CurrentMaxTiStep stores  1/a^2 times the maximum allowed physical timestep */
          if(Tp->SphP[target].CurrentMaxTiStep >= All.MaxSizeTimestep / All.cf_atime2_hubble_a / All.CourantFac)
            Tp->SphP[target].CurrentMaxTiStep = All.MaxSizeTimestep / All.cf_atime2_hubble_a / All.CourantFac;

          WorkStack[NumOnWorkStack].Target         = target;
          WorkStack[NumOnWorkStack].Node           = MaxPart;
          WorkStack[NumOnWorkStack].ShmRank        = Shmem.Island_ThisTask;
          WorkStack[NumOnWorkStack].MinTopLeafNode = MaxPart + D->NTopnodes;
          NumOnWorkStack++;
        }

      // set a default size of the fetch stack equal to half the work stack (this may still be somewhat too large)
      MaxOnFetchStack = std::max<int>(0.1 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
      StackToFetch    = (fetch_data *)Mem.mymalloc_movable(&StackToFetch, "StackToFetch", MaxOnFetchStack * sizeof(fetch_data));

      while(NumOnWorkStack > 0)  // repeat until we are out of work
        {
          NewOnWorkStack  = 0;  // gives the new entries
          NumOnFetchStack = 0;
          MaxOnWorkStack  = std::min<int>(AllocWorkStackBaseLow + max_ncycles * TREE_MIN_WORKSTACK_SIZE, AllocWorkStackBaseHigh);

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
                      sph_treetimestep_interact(pdat, no, NODE_TYPE_LOCAL_NODE, shmrank, mintopleaf, committed);
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
                        sph_treetimestep_open_node(pdat, nop, mintopleaf, committed);
                    }
                }
              else
                break;
            }

          if(item == 0 && NumOnWorkStack > 0)
            Terminate("Can't even process a single particle");

          tree_fetch_foreign_nodes(FETCH_SPH_TREETIMESTEP);

          /* now reorder the workstack such that we are first going to do residual pristine particles, and then
           * imported nodes that hang below the first leaf nodes */
          NumOnWorkStack = NumOnWorkStack - item + NewOnWorkStack;
          memmove(WorkStack, WorkStack + item, NumOnWorkStack * sizeof(workstack_data));

          /* now let's sort such that we can go deep on top-level node branches, allowing us to clear them out eventually */
          mycxxsort(WorkStack, WorkStack + NumOnWorkStack, compare_workstack);

          max_ncycles++;
        }

      Mem.myfree(StackToFetch);
      Mem.myfree(WorkStack);

      /* now multiply the determined values with the specified CourantFactor */
      for(int i = 0; i < Tp->TimeBinsHydro.NActiveParticles; i++)
        {
          int target = Tp->TimeBinsHydro.ActiveParticleList[i];

          Tp->SphP[target].CurrentMaxTiStep *= All.CourantFac;
        }

      MPI_Allreduce(MPI_IN_PLACE, &max_ncycles, 1, MPI_INT, MPI_MAX, D->Communicator);

      cleanup_shared_memory_access();

      /* free temporary buffers */
      Mem.myfree(Foreign_Points);
      Mem.myfree(Foreign_Nodes);

      double tb = Logs.second();

      TIMER_STOPSTART(CPU_TREE_TIMESTEPS, CPU_LOGS);

      D->mpi_printf("TIMESTEP-TREEWALK: took %g sec, max_ncycles = %d, part/sec = %g\n", Logs.timediff(ta, tb), max_ncycles,
                    Tp->TimeBinsHydro.GlobalNActiveParticles / (D->NTask * Logs.timediff(ta, tb) + MIN_FLOAT_NUMBER));

      TIMER_STOP(CPU_LOGS);
    }
}
