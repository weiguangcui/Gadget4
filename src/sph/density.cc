/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  density.cc
 *
 *  \brief SPH density computation and smoothing length determination
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
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../sort/cxxsort.h"
#include "../sph/kernel.h"
#include "../sph/sph.h"
#include "../system/system.h"

/*! This file contains the function for the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  iteratively corrects the smoothing length to the desired value.
 */

/* This function checks whether there is a spatial overlap between the (rectangular) enclosing box
 * of the particles contained in a node, and the search region.
 */
inline int sph::sph_density_evaluate_particle_node_opening_criterion(pinfo &pdat, ngbnode *nop)
{
  if(nop->level <= LEVEL_ALWAYS_OPEN)  // always open the root node (note: full node length does not fit in the integer type)
    return NODE_OPEN;

  if(nop->Ti_Current != All.Ti_Current)
    nop->drift_node(All.Ti_Current, Tp);

  MyIntPosType left[3], right[3];

  left[0]  = Tp->nearest_image_intpos_to_intpos_X(nop->center_offset_min[0] + nop->center[0], pdat.search_min[0]);
  right[0] = Tp->nearest_image_intpos_to_intpos_X(nop->center_offset_max[0] + nop->center[0], pdat.search_min[0]);

  /* check whether we can stop walking along this branch */
  if(left[0] > pdat.search_range[0] && right[0] > left[0])
    return NODE_DISCARD;

  left[1]  = Tp->nearest_image_intpos_to_intpos_Y(nop->center_offset_min[1] + nop->center[1], pdat.search_min[1]);
  right[1] = Tp->nearest_image_intpos_to_intpos_Y(nop->center_offset_max[1] + nop->center[1], pdat.search_min[1]);

  /* check whether we can stop walking along this branch */
  if(left[1] > pdat.search_range[1] && right[1] > left[1])
    return NODE_DISCARD;

  left[2]  = Tp->nearest_image_intpos_to_intpos_Z(nop->center_offset_min[2] + nop->center[2], pdat.search_min[2]);
  right[2] = Tp->nearest_image_intpos_to_intpos_Z(nop->center_offset_max[2] + nop->center[2], pdat.search_min[2]);

  /* check whether we can stop walking along this branch */
  if(left[2] > pdat.search_range[2] && right[2] > left[2])
    return NODE_DISCARD;

  return NODE_OPEN;
}

/* Check whether the potential neighbor referenced by p/p_type/shmrank is inside the smoothing length, and if yes
 * add it to the interaction list built up for this particle.
 */
inline void sph::sph_density_check_particle_particle_interaction(pinfo &pdat, int p, int p_type, unsigned char shmrank)
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

      double posdiff[3];
      Tp->nearest_image_intpos_to_pos(P->IntPos, pdat.searchcenter, posdiff); /* converts the integer distance to floating point */

      if(posdiff[0] * posdiff[0] + posdiff[1] * posdiff[1] + posdiff[2] * posdiff[2] > pdat.hsml2)
        return;

      if(pdat.numngb >= MAX_NGBS)
        Terminate("pdat.numngb >= MAX_NGBS");

      int n = pdat.numngb++;

      Ngbdensdat[n].IntPos  = P->IntPos;
      Ngbdensdat[n].VelPred = SphP->VelPred;
      Ngbdensdat[n].Mass    = P->getMass();
#ifdef PRESSURE_ENTROPY_SPH
      Ngbdensdat[n].EntropyToInvGammaPred = SphP->EntropyToInvGammaPred;
#endif
#ifdef TIMEDEP_ART_VISC
      Ngbdensdat[n].Csnd = SphP->Csnd;
#endif
    }
  else if(p_type == NODE_TYPE_FETCHED_PARTICLE)
    {
      foreign_sphpoint_data *foreignpoint = get_foreignpointsp(p - EndOfForeignNodes, shmrank);

      /* converts the integer distance to floating point */
      double posdiff[3];
      Tp->nearest_image_intpos_to_pos(foreignpoint->IntPos, pdat.searchcenter, posdiff);

      if(posdiff[0] * posdiff[0] + posdiff[1] * posdiff[1] + posdiff[2] * posdiff[2] > pdat.hsml2)
        return;

      if(pdat.numngb >= MAX_NGBS)
        Terminate("pdat.numngb >= MAX_NGBS");

      int n = pdat.numngb++;

      Ngbdensdat[n].IntPos  = foreignpoint->IntPos;
      Ngbdensdat[n].VelPred = foreignpoint->SphCore.VelPred;
      Ngbdensdat[n].Mass    = foreignpoint->Mass;
#ifdef PRESSURE_ENTROPY_SPH
      Ngbdensdat[n].EntropyToInvGammaPred = foreignpoint->SphCore.EntropyToInvGammaPred;
#endif
#ifdef TIMEDEP_ART_VISC
      Ngbdensdat[n].Csnd = foreignpoint->SphCore.Csnd;
#endif
    }
  else
    Terminate("unexpected");
}

/* Continues to walk the tree for the particle referenced in pdat by opening a node.
 */
inline void sph::sph_density_open_node(pinfo &pdat, ngbnode *nop, int mintopleafnode, int committed)
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
              "shmrank=%d  nop->nextnode=%d  nop->cannot_be_opened_locally=%d  nop->not_empty=%d  nop-TopNodes=%d",
              p, MaxPart, MaxNodes, ImportedNodeOffset, EndOfTreePoints, EndOfForeignNodes, shmrank, nop->nextnode,
              (int)nop->cannot_be_opened_locally, (int)nop->not_empty, (int)(nop - TopNodes));
        }

      sph_density_interact(pdat, p, type, shmrank, mintopleafnode, committed);

      p       = next;
      shmrank = next_shmrank;
    }
}

/* Take care of SPH density interaction between the particle referenced in pdat, and the node
 * referenced through no/shmrank. The node can either be a normal node or an imported node from another
 * shared memory machine, and the node can either be on the present MPI rank, or from another MPI rank in the
 * local shared memory machine.
 */
inline void sph::sph_density_interact(pinfo &pdat, int no, char no_type, unsigned char shmrank, int mintopleafnode, int committed)
{
  if(no_type <= NODE_TYPE_FETCHED_PARTICLE)  // we are interacting with a particle
    {
      sph_density_check_particle_particle_interaction(pdat, no, no_type, shmrank);
    }
  else  // we are interacting with a node
    {
      ngbnode *nop = get_nodep(no, shmrank);

      if(nop->not_empty == 0)
        return;

      if(no < MaxPart + MaxNodes)                // we have a top-levelnode
        if(nop->nextnode >= MaxPart + MaxNodes)  // if the next node is not a top-level, we have a leaf node
          mintopleafnode = no;

      int openflag = sph_density_evaluate_particle_node_opening_criterion(pdat, nop);

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
                sph_density_open_node(pdat, nop, mintopleafnode, committed + 8 * TREE_NUM_BEFORE_NODESPLIT);
              else
                tree_add_to_work_stack(pdat.target, no, shmrank, mintopleafnode);
            }
        }
    }
}

/* Internal driver routine to compute densities for the particles with indices listed in the targetlist
 * array.
 */
void sph::densities_determine(int ntarget, int *targetlist)
{
  Ngbdensdat = (ngbdata_density *)Mem.mymalloc("Ngbdensdat", MAX_NGBS * sizeof(ngbdata_density));

  NumOnWorkStack         = 0;
  AllocWorkStackBaseLow  = std::max<int>(1.5 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
  AllocWorkStackBaseHigh = AllocWorkStackBaseLow + TREE_EXPECTED_CYCLES * TREE_MIN_WORKSTACK_SIZE;
  MaxOnWorkStack         = AllocWorkStackBaseLow;

  WorkStack = (workstack_data *)Mem.mymalloc("WorkStack", AllocWorkStackBaseHigh * sizeof(workstack_data));

  for(int i = 0; i < ntarget; i++)
    {
      int target = targetlist[i];

      clear_density_result(&Tp->SphP[target]);

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

          TIMER_START(CPU_DENSWALK);

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
                      sph_density_interact(pdat, no, NODE_TYPE_LOCAL_NODE, shmrank, mintopleaf, committed);
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
                        sph_density_open_node(pdat, nop, mintopleaf, committed);
                    }

                  density_evaluate_kernel(pdat);
                }
              else
                break;
            }

          if(item == 0 && NumOnWorkStack > 0)
            Terminate("Can't even process a single particle");

          TIMER_STOP(CPU_DENSWALK);

          TIMER_START(CPU_DENSFETCH);

          tree_fetch_foreign_nodes(FETCH_SPH_DENSITY);

          TIMER_STOP(CPU_DENSFETCH);

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
  Mem.myfree(Ngbdensdat);
}

/* This first makes sure that all active SPH particles are drifted to the current time,
 * and then calls the SPH density computation for them.
 */
void sph::compute_densities(void)
{
  if(Tp->TimeBinsHydro.GlobalNActiveParticles > 0)
    {
      /* now drift the active hydro particles if not done already */
      for(int i = 0; i < Tp->TimeBinsHydro.NActiveParticles; i++)
        {
          particle_data *P        = &Tp->P[Tp->TimeBinsHydro.ActiveParticleList[i]];
          sph_particle_data *SphP = &Tp->SphP[Tp->TimeBinsHydro.ActiveParticleList[i]];

          Tp->drift_particle(P, SphP, All.Ti_Current);
#ifdef TIMEDEP_ART_VISC
          SphP->DivVelOld = SphP->DivVel;
#endif
        }

      /* compute density (and updates pressure) */
      density(Tp->TimeBinsHydro.ActiveParticleList, Tp->TimeBinsHydro.NActiveParticles);
    }
}

/* Compute the SPH densities for all particles listed in the index-array. The routine finds a suitable
 * smoothing length by a bisection algorithm.
 */
void sph::density(int *list, int ntarget)
{
  TIMER_STORE;
  TIMER_START(CPU_DENSITY);

  D->mpi_printf("SPH-DENSITY: Begin density calculation. (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());
  D->mpi_printf("SPH-DENSITY: Ndensities=%llu (task zero has: NumGas=%d, Ndensities=%d)\n", Tp->TimeBinsHydro.GlobalNActiveParticles,
                Tp->NumGas, ntarget);

  double ta = Logs.second();

  /* Create list of targets. We do this here to simplify the treatment later on */
  int *targetList = (int *)Mem.mymalloc("TargetList", Tp->NumGas * sizeof(int));
  MyFloat *Left   = (MyFloat *)Mem.mymalloc("Left", Tp->NumGas * sizeof(MyFloat));
  MyFloat *Right  = (MyFloat *)Mem.mymalloc("Right", Tp->NumGas * sizeof(MyFloat));

  int ndensities = ntarget;

  for(int i = 0; i < ntarget; i++)
    {
      int target    = list[i];
      targetList[i] = target;
      Left[target] = Right[target] = 0.0;
    }

  int iter = 0;

  // let's grab at most half the still available memory for imported points and nodes
  int nspace = (0.5 * Mem.FreeBytes) / (sizeof(ngbnode) + 8 * sizeof(foreign_sphpoint_data));

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

  do
    {
      double t0 = Logs.second();

      /* now do the primary work with this call */

      densities_determine(ndensities, targetList);

      /* do final operations on results */
      int npleft = 0;

      for(int i = 0; i < ndensities; i++)
        {
          int target = targetList[i];

          if(target >= 0)
            {
              if(Tp->P[target].getType() != 0)
                Terminate("P[target].getType() != 0");

              sph_particle_data *SphP = Tp->SphP;
              if(SphP[target].Density > 0)
                {
#ifdef WENDLAND_BIAS_CORRECTION
                  SphP[target].Density -= get_density_bias(SphP[target].Hsml, Tp->P[target].getMass(), All.DesNumNgb);
#endif
                  SphP[target].DhsmlDensityFactor *= SphP[target].Hsml / (NUMDIMS * SphP[target].Density);
                  if(SphP[target].DhsmlDensityFactor >
                     -0.9) /* note: this would be -1 if only a single particle at zero lag is found */
                    SphP[target].DhsmlDensityFactor = 1 / (1 + SphP[target].DhsmlDensityFactor);
                  else
                    SphP[target].DhsmlDensityFactor = 1;

#ifndef IMPROVED_VELOCITY_GRADIENTS
                  SphP[target].CurlVel = sqrt(SphP[target].Rot[0] * SphP[target].Rot[0] + SphP[target].Rot[1] * SphP[target].Rot[1] +
                                              SphP[target].Rot[2] * SphP[target].Rot[2]) /
                                         SphP[target].Density;

                  SphP[target].DivVel /= SphP[target].Density;
#else
                  SphP[target].set_velocity_gradients();
#endif
                  SphP[target].DtHsml    = (1.0 / NUMDIMS) * SphP[target].DivVel * SphP[target].Hsml;
                  SphP[target].DtDensity = -SphP[target].DivVel * SphP[target].Density;

#ifndef PRESSURE_ENTROPY_SPH
                  SphP[target].set_thermodynamic_variables();
#endif
                }

#ifdef PRESSURE_ENTROPY_SPH
              if(SphP[target].EntropyToInvGammaPred > 0 && SphP[target].PressureSphDensity > 0)
                {
                  SphP[target].DhsmlDerivedDensityFactor *=
                      SphP[target].Hsml / (NUMDIMS * SphP[target].Density * SphP[target].EntropyToInvGammaPred);
                  SphP[target].DhsmlDerivedDensityFactor *= -SphP[target].DhsmlDensityFactor;
                  SphP[target].PressureSphDensity /= SphP[target].EntropyToInvGammaPred;
#ifdef WENDLAND_BIAS_CORRECTION /* Dehnen & Aly 2012, eq (18), (19) */
                  SphP[target].PressureSphDensity -= get_density_bias(SphP[target].Hsml, Tp->P[target].getMass(), All.DesNumNgb);
#endif
                  SphP[target].DtPressureSphDensity = -SphP[target].DivVel * SphP[target].PressureSphDensity;
                  SphP[target].set_thermodynamic_variables();
                }
              else
                {
                  SphP[target].DhsmlDerivedDensityFactor = 0;
                  SphP[target].EntropyToInvGammaPred     = 0;
                  SphP[target].PressureSphDensity        = 0;
                }

#endif
#ifdef ADAPTIVE_HYDRO_SOFTENING
              Tp->P[target].setSofteningClass(Tp->get_softeningtype_for_hydro_particle(target));
#endif
              /* now check whether we had enough neighbours */
              double desnumngb    = All.DesNumNgb;
              double desnumngbdev = All.MaxNumNgbDeviation;

              double hfac = 1;
              for(int i = 0; i < NUMDIMS; i++)
                {
                  hfac *= SphP[target].Hsml;
                }

              SphP[target].NumNgb = NORM_COEFF * hfac * SphP[target].Density / Tp->P[target].getMass();
              if(SphP[target].NumNgb < (desnumngb - desnumngbdev) || (SphP[target].NumNgb > (desnumngb + desnumngbdev)))
                {
                  if(Left[target] > 0 && Right[target] > 0)
                    if((Right[target] - Left[target]) < 1.0e-3 * Left[target])
                      {
                        /* this one should be ok */
                        continue;
                      }

                  /* need to redo this particle */
                  targetList[npleft++] = target;

                  if(SphP[target].NumNgb < (desnumngb - desnumngbdev))
                    Left[target] = std::max<double>(SphP[target].Hsml, Left[target]);
                  else
                    {
                      if(Right[target] != 0)
                        {
                          if(SphP[target].Hsml < Right[target])
                            Right[target] = SphP[target].Hsml;
                        }
                      else
                        Right[target] = SphP[target].Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      double pos[3];
                      Tp->intpos_to_pos(Tp->P[target].IntPos, pos); /* converts the integer coordinates to floating point */

                      printf("target=%d Hsml=%g  task=%d ID=%llu Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", target,
                             SphP[target].Hsml, D->ThisTask, (unsigned long long)Tp->P[target].ID.get(), Left[target], Right[target],
                             SphP[target].NumNgb, Right[target] - Left[target], pos[0], pos[1], pos[2]);
                      myflush(stdout);
                    }

                  if(Right[target] > 0 && Left[target] > 0)
                    SphP[target].Hsml = pow(0.5 * (pow(Left[target], 3) + pow(Right[target], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[target] == 0 && Left[target] == 0)
                        Terminate("Right[i] == 0 && Left[i] == 0 SphP[i].Hsml=%g\n", SphP[target].Hsml);

                      if(Right[target] == 0 && Left[target] > 0)
                        {
                          if(Tp->P[target].getType() == 0 && fabs(SphP[target].NumNgb - desnumngb) < 0.5 * desnumngb)
                            {
                              double fac = 1 - (SphP[target].NumNgb - desnumngb) / (NUMDIMS * SphP[target].NumNgb) *
                                                   SphP[target].DhsmlDensityFactor;

                              if(fac < 1.26)
                                SphP[target].Hsml *= fac;
                              else
                                SphP[target].Hsml *= 1.26;
                            }
                          else
                            SphP[target].Hsml *= 1.26;
                        }

                      if(Right[target] > 0 && Left[target] == 0)
                        {
                          if(Tp->P[target].getType() == 0 && fabs(SphP[target].NumNgb - desnumngb) < 0.5 * desnumngb && iter < 4)
                            {
                              double fac = 1 - (SphP[target].NumNgb - desnumngb) / (NUMDIMS * SphP[target].NumNgb) *
                                                   SphP[target].DhsmlDensityFactor;

                              if(fac > 1 / 1.26)
                                SphP[target].Hsml *= fac;
                              else
                                SphP[target].Hsml /= 1.26;
                            }
                          else
                            SphP[target].Hsml /= 1.26;
                        }
                    }
                }
            }
        }

      ndensities = npleft;

      double t1 = Logs.second();

      if(npleft > 0)
        {
          iter++;

          D->mpi_printf("SPH-DENSITY: ngb iteration %4d: took %8.3f  , need to repeat for %012lld local particles.\n", iter,
                        Logs.timediff(t0, t1), npleft);

          if(iter > MAXITER)
            Terminate("failed to converge in neighbour iteration in density()\n");
        }
      else
        D->mpi_printf("SPH-DENSITY: ngb iteration %4d: took %8.3f\n", ++iter, Logs.timediff(t0, t1));
    }
  while(ndensities > 0);

#ifdef TIMEDEP_ART_VISC
  for(int i = 0; i < ntarget; i++)
    {
      int target = list[i];
      double dt =
          (Tp->P[target].getTimeBinHydro() ? (((integertime)1) << Tp->P[target].getTimeBinHydro()) : 0) * All.Timebase_interval;
      double dtime = All.cf_atime * dt / All.cf_atime_hubble_a;
      Tp->SphP[target].set_viscosity_coefficient(dtime);
    }
#endif

  TIMER_START(CPU_DENSIMBALANCE);

  MPI_Allreduce(MPI_IN_PLACE, &max_ncycles, 1, MPI_INT, MPI_MAX, D->Communicator);

  TIMER_STOP(CPU_DENSIMBALANCE);

  cleanup_shared_memory_access();

  /* free temporary buffers */

  Mem.myfree(Foreign_Points);
  Mem.myfree(Foreign_Nodes);

  Mem.myfree(Right);
  Mem.myfree(Left);
  Mem.myfree(targetList);

  double tb = Logs.second();

  TIMER_STOPSTART(CPU_DENSITY, CPU_LOGS);

  D->mpi_printf("SPH-DENSITY: density computation done. took %8.3f\n", Logs.timediff(ta, tb));

  struct detailed_timings
  {
    double tree, wait, fetch, all;
    double numnodes;
    double NumForeignNodes, NumForeignPoints;
    double fillfacFgnNodes, fillfacFgnPoints;
  };
  detailed_timings timer, tisum, timax;

  timer.tree             = TIMER_DIFF(CPU_DENSWALK);
  timer.wait             = TIMER_DIFF(CPU_DENSIMBALANCE);
  timer.fetch            = TIMER_DIFF(CPU_DENSFETCH);
  timer.all              = timer.tree + timer.wait + timer.fetch + TIMER_DIFF(CPU_DENSITY);
  timer.numnodes         = NumNodes;
  timer.NumForeignNodes  = NumForeignNodes;
  timer.NumForeignPoints = NumForeignPoints;
  timer.fillfacFgnNodes  = NumForeignNodes / ((double)MaxForeignNodes);
  timer.fillfacFgnPoints = NumForeignPoints / ((double)MaxForeignPoints);

  MPI_Reduce((double *)&timer, (double *)&tisum, (int)(sizeof(detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_SUM, 0,
             D->Communicator);
  MPI_Reduce((double *)&timer, (double *)&timax, (int)(sizeof(detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_MAX, 0,
             D->Communicator);

  All.TotNumDensity += Tp->TimeBinsHydro.GlobalNActiveParticles;

  if(D->ThisTask == 0)
    {
      fprintf(Logs.FdDensity, "Nf=%9lld  highest active timebin=%d  total-Nf=%lld\n", Tp->TimeBinsHydro.GlobalNActiveParticles,
              All.HighestActiveTimeBin, All.TotNumDensity);
      fprintf(Logs.FdDensity, "   work-load balance: %g   part/sec: raw=%g, effective=%g\n",
              timax.tree / ((tisum.tree + 1e-20) / D->NTask), Tp->TimeBinsGravity.GlobalNActiveParticles / (tisum.tree + 1.0e-20),
              Tp->TimeBinsGravity.GlobalNActiveParticles / ((timax.tree + 1.0e-20) * D->NTask));
      fprintf(Logs.FdDensity,
              "   maximum number of nodes: %g, filled: %g  NumForeignNodes: max=%g avg=%g fill=%g NumForeignPoints: max=%g avg=%g "
              "fill=%g  cycles=%d\n",
              timax.numnodes, timax.numnodes / MaxNodes, timax.NumForeignNodes, tisum.NumForeignNodes / D->NTask,
              timax.fillfacFgnNodes, timax.NumForeignPoints, tisum.NumForeignPoints / D->NTask, timax.fillfacFgnPoints, max_ncycles);
      fprintf(Logs.FdDensity, "   avg times: <all>=%g  <tree>=%g  <wait>=%g  <fetch>=%g  sec\n", tisum.all / D->NTask,
              tisum.tree / D->NTask, tisum.wait / D->NTask, tisum.fetch / D->NTask);
      myflush(Logs.FdDensity);
    }

  TIMER_STOP(CPU_LOGS);
}

#ifdef EXPLICIT_VECTORIZATION

/* Main SPH compute kernel, in a version that is explicitely vectorized. Several neighbours are
 * processed through one vector. The calculation should be semantically equivalent to the standard
 * looped version without explicit vector instructions.
 */
void sph::density_evaluate_kernel(pinfo &pdat)
{
  particle_data *targetP        = &Tp->P[pdat.target];
  sph_particle_data *targetSphP = &Tp->SphP[pdat.target];

  double shinv, shinv3, shinv4;
  kernel_hinv(targetSphP->Hsml, &shinv, &shinv3, &shinv4);

  Vec4d hinv(shinv);
  Vec4d hinv3(shinv3);
  Vec4d hinv4(shinv4);

  Vec4d dwnorm(NORM * shinv3);
  Vec4d dwknorm(NORM * shinv4);

  Vec4d v_i[NUMDIMS];
  for(int i = 0; i < NUMDIMS; i++)
    {
      v_i[i] = targetSphP->VelPred[i];
    }

  Vec4d cs_i(targetSphP->Csnd);
  const int vector_length = 4;
  const int array_length  = (pdat.numngb + vector_length - 1) & (-vector_length);

  for(int n = pdat.numngb; n < array_length; n++) /* fill up neighbour array so that sensible data is accessed */
    Ngbdensdat[n] = Ngbdensdat[0];

  Vec4d decayVel_loc(0);

  for(int n = 0; n < array_length; n += vector_length)
    {
      struct ngbdata_density *ngb0 = &Ngbdensdat[n + 0];
      struct ngbdata_density *ngb1 = &Ngbdensdat[n + 1];
      struct ngbdata_density *ngb2 = &Ngbdensdat[n + 2];
      struct ngbdata_density *ngb3 = &Ngbdensdat[n + 3];

      Vec4d dpos[NUMDIMS];
#if defined(LONG_X_BITS) || defined(LONG_Y_BITS) || defined(LONG_Z_BITS)
      double posdiff[array_length][3];
      for(int i = 0; i < 4; i++)
        {
          Tp->nearest_image_intpos_to_pos(targetP->IntPos, Ngbdensdat[n + i].IntPos, &(posdiff[i][0]));
        }

      for(int i = 0; i < NUMDIMS; i++)
        {
          dpos[i] = Vec4d(posdiff[0][i], posdiff[1][i], posdiff[2][i], posdiff[3][i]);
        }
#else
      for(int i = 0; i < NUMDIMS; i++)
        {
          dpos[i] = Tp->nearest_image_intpos_to_doublepos_vectorial(targetP->IntPos[i], ngb0->IntPos[i], ngb1->IntPos[i],
                                                                    ngb2->IntPos[i], ngb3->IntPos[i]);
        }
#endif

      Vec4d v_j[NUMDIMS];
      for(int i = 0; i < NUMDIMS; i++)
        {
          v_j[i] = Vec4d(ngb0->VelPred[i], ngb1->VelPred[i], ngb2->VelPred[i], ngb3->VelPred[i]);
        }

      Vec4d mass_j(ngb0->Mass, ngb1->Mass, ngb2->Mass, ngb3->Mass);
      Vec4d r2(0);

      for(int i = 0; i < NUMDIMS; i++)
        {
          r2 += dpos[i] * dpos[i];
        }

      Vec4d r = sqrt(r2);

      Vec4d u = r * hinv;

      /* now calculate the kernel */
      Vec4d wk, dwk;
      kernel_main_vector(u, dwnorm, dwknorm, &wk, &dwk);

      if(n + vector_length > pdat.numngb) /* we have excess elements */
        {
          mass_j.cutoff(vector_length - (array_length - pdat.numngb));
          wk.cutoff(vector_length - (array_length - pdat.numngb));
        }

      Vec4d mj_wk = mass_j * wk;

      targetSphP->Density += horizontal_add(mj_wk);

#ifdef PRESSURE_ENTROPY_SPH
      Vec4d entr_j(ngb0->EntropyToInvGammaPred, ngb1->EntropyToInvGammaPred, ngb2->EntropyToInvGammaPred, ngb3->EntropyToInvGammaPred);

      targetSphP->PressureSphDensity += horizontal_add(mj_wk * entr_j);

      targetSphP->DhsmlDerivedDensityFactor += horizontal_add(-mass_j * entr_j * (NUMDIMS * hinv * wk + u * dwk));
#endif

      targetSphP->DhsmlDensityFactor += horizontal_add(-mass_j * (NUMDIMS * hinv * wk + u * dwk));

      Vec4db decision_r_gt_0 = (r > 0);

      r = select(decision_r_gt_0, r, 1.0); /* note, for r=0, we have dwk=0 */

      Vec4d mj_dwk_r = mass_j * dwk / r;

      Vec4d dv[NUMDIMS];
      for(int i = 0; i < NUMDIMS; i++)
        {
          dv[i] = v_i[i] - v_j[i];
        }

#ifndef IMPROVED_VELOCITY_GRADIENTS
      Vec4d dpos_times_dv(0);
      for(int i = 0; i < NUMDIMS; i++)
        dpos_times_dv += dpos[i] * dv[i];

      targetSphP->DivVel += horizontal_add(-mj_dwk_r * dpos_times_dv);

#ifdef TWODIMS
      targetSphP->Rot[2] += horizontal_add(mj_dwk_r * (dpos[1] * dv[0] - dpos[0] * dv[1]));
#endif
#ifdef THREEDIMS
      targetSphP->Rot[0] += horizontal_add(mj_dwk_r * (dpos[2] * dv[1] - dpos[1] * dv[2]));
      targetSphP->Rot[1] += horizontal_add(mj_dwk_r * (dpos[0] * dv[2] - dpos[2] * dv[0]));
      targetSphP->Rot[2] += horizontal_add(mj_dwk_r * (dpos[1] * dv[0] - dpos[0] * dv[1]));
#endif
#else
      for(int i = 0; i < NUMDIMS; i++)
        {
          for(int j = 0; j < NUMDIMS; j++)
            {
              targetSphP->dvel[i][j] -= horizontal_add(mj_dwk_r * dv[i] * dpos[j]);
            }
        }
      targetSphP->dpos.dx_dx -= horizontal_add(mj_dwk_r * dpos[0] * dpos[0]);
      targetSphP->dpos.dx_dy -= horizontal_add(mj_dwk_r * dpos[0] * dpos[1]);
      targetSphP->dpos.dx_dz -= horizontal_add(mj_dwk_r * dpos[0] * dpos[2]);
      targetSphP->dpos.dy_dy -= horizontal_add(mj_dwk_r * dpos[1] * dpos[1]);
      targetSphP->dpos.dy_dz -= horizontal_add(mj_dwk_r * dpos[1] * dpos[2]);
      targetSphP->dpos.dz_dz -= horizontal_add(mj_dwk_r * dpos[2] * dpos[2]);
#endif
#ifdef TIMEDEP_ART_VISC
      Vec4d vdotr2 = 0;
      for(int i = 0; i < NUMDIMS; i++)
        {
          vdotr2 += dpos[i] * dv[i];
        }

      if(All.ComovingIntegrationOn)
        vdotr2 += All.cf_atime2_hubble_a * r2;

      Vec4d mu_ij = vdotr2 / (All.cf_afac3 * All.cf_atime * r);

      Vec4d cs_j(ngb0->Csnd, ngb1->Csnd, ngb2->Csnd, ngb3->Csnd);
      Vec4d cs_sum = cs_i + cs_j;

      Vec4d decay_vel_2 = cs_sum - mu_ij;

      Vec4db decision = (vdotr2 > 0);

      Vec4d decay_vel = select(decision, cs_sum, decay_vel_2);

      Vec4d fac_decay_vel = select(decision_r_gt_0, 1, 0);

      decay_vel = decay_vel * fac_decay_vel;

      decayVel_loc = max(decayVel_loc, decay_vel);

#endif
    }

#ifdef TIMEDEP_ART_VISC
  for(int i = 0; i < vector_length; i++)
    {
      if(decayVel_loc[i] > targetSphP->decayVel)
        targetSphP->decayVel = decayVel_loc[i];
    }
#endif
}

#else

/* Main SPH compute kernel. This function carries out the SPH computations for the neighbouring particle
 * data stored in the Ngbdensdat[] array. The results are added to the particle referenced through the pdat
 * structure.
 */
void sph::density_evaluate_kernel(pinfo &pdat)
{
  particle_data *targetP = &Tp->P[pdat.target];
  sph_particle_data *targetSphP = &Tp->SphP[pdat.target];

  kernel_density kernel;

  kernel_hinv(targetSphP->Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

  for(int n = 0; n < pdat.numngb; n++)
    {
      struct ngbdata_density *ngb = &Ngbdensdat[n];

      /*  note: in periodic case, closest image will be found through integer wrap around  */

      double posdiff[3];
      Tp->nearest_image_intpos_to_pos(targetP->IntPos, ngb->IntPos, posdiff); /* converts the integer distance to floating point */

      for(int i = 0; i < NUMDIMS; i++)
        {
          kernel.dpos[i] = posdiff[i];
        }

      double r2 = 0;

      for(int i = 0; i < NUMDIMS; i++)
        {
          r2 += kernel.dpos[i] * kernel.dpos[i];
        }

      kernel.r = sqrt(r2);

      double u = kernel.r * kernel.hinv;

      kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, COMPUTE_WK_AND_DWK);

      double mass_j = ngb->Mass;
      kernel.mj_wk = (mass_j * kernel.wk);

      targetSphP->Density += kernel.mj_wk;

#ifdef PRESSURE_ENTROPY_SPH
      targetSphP->PressureSphDensity += kernel.mj_wk * ngb->EntropyToInvGammaPred;
      targetSphP->DhsmlDerivedDensityFactor +=
          (-mass_j * ngb->EntropyToInvGammaPred * (NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk));

#endif

      targetSphP->DhsmlDensityFactor += (-mass_j * (NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk));

      if(kernel.r > 0)
        {
          kernel.mj_dwk_r = mass_j * kernel.dwk / kernel.r;

          for(int i = 0; i < NUMDIMS; i++)
            {
              kernel.dv[i] = targetSphP->VelPred[i] - ngb->VelPred[i];
            }

#ifndef IMPROVED_VELOCITY_GRADIENTS
          double dpos_times_dv = 0;
          for(int i = 0; i < NUMDIMS; i++)
            dpos_times_dv += kernel.dpos[i] * kernel.dv[i];

          targetSphP->DivVel += (-kernel.mj_dwk_r * (dpos_times_dv));
#ifdef TWODIMS
          targetSphP->Rot[2] += (kernel.mj_dwk_r * (kernel.dpos[1] * kernel.dv[0] - kernel.dpos[0] * kernel.dv[1]));
#endif
#ifdef THREEDIMS
          targetSphP->Rot[0] += (kernel.mj_dwk_r * (kernel.dpos[2] * kernel.dv[1] - kernel.dpos[1] * kernel.dv[2]));
          targetSphP->Rot[1] += (kernel.mj_dwk_r * (kernel.dpos[0] * kernel.dv[2] - kernel.dpos[2] * kernel.dv[0]));
          targetSphP->Rot[2] += (kernel.mj_dwk_r * (kernel.dpos[1] * kernel.dv[0] - kernel.dpos[0] * kernel.dv[1]));
#endif
#else
          for(int i = 0; i < NUMDIMS; i++)
            {
              for(int j = 0; j < NUMDIMS; j++)
                {
                  targetSphP->dvel[i][j] -= kernel.mj_dwk_r * kernel.dv[i] * kernel.dpos[j];
                }
            }

          targetSphP->dpos.dx_dx -= kernel.mj_dwk_r * kernel.dpos[0] * kernel.dpos[0];
          targetSphP->dpos.dx_dy -= kernel.mj_dwk_r * kernel.dpos[0] * kernel.dpos[1];
          targetSphP->dpos.dx_dz -= kernel.mj_dwk_r * kernel.dpos[0] * kernel.dpos[2];
          targetSphP->dpos.dy_dy -= kernel.mj_dwk_r * kernel.dpos[1] * kernel.dpos[1];
          targetSphP->dpos.dy_dz -= kernel.mj_dwk_r * kernel.dpos[1] * kernel.dpos[2];
          targetSphP->dpos.dz_dz -= kernel.mj_dwk_r * kernel.dpos[2] * kernel.dpos[2];

#endif
#ifdef TIMEDEP_ART_VISC
          double vdotr2 = 0;
          for(int i = 0; i < NUMDIMS; i++)
            {
              vdotr2 += kernel.dpos[i] * kernel.dv[i];
            }

          if(All.ComovingIntegrationOn)
            vdotr2 += All.cf_atime2_hubble_a * r2;

          double mu_ij = vdotr2 / (All.cf_afac3 * All.cf_atime * kernel.r);
          double decay_vel;

          if(vdotr2 < 0)
            decay_vel = targetSphP->Csnd + ngb->Csnd - mu_ij;
          else
            decay_vel = targetSphP->Csnd + ngb->Csnd;

          if(decay_vel > targetSphP->decayVel)
            targetSphP->decayVel = decay_vel;
#endif
        }
    }
}
#endif

/* this routine clears the fields in the SphP particle structure that are additively computed by the SPH density loop
 * by summing over neighbours
 */
inline void sph::clear_density_result(sph_particle_data *SphP)
{
  SphP->Density            = 0;
  SphP->DhsmlDensityFactor = 0;
  SphP->DivVel             = 0;

  for(int k = 0; k < 3; k++)
    SphP->Rot[k] = 0;

#ifdef PRESSURE_ENTROPY_SPH
  SphP->PressureSphDensity        = 0;
  SphP->DhsmlDerivedDensityFactor = 0;
#endif
#ifdef IMPROVED_VELOCITY_GRADIENTS
  SphP->dpos = {0};
  for(int i = 0; i < NUMDIMS; i++)
    {
      for(int j = 0; j < NUMDIMS; j++)
        {
          SphP->dvel[i][j] = 0;
        }
    }
#endif
#ifdef TIMEDEP_ART_VISC
  SphP->decayVel = 0;
#endif
}
