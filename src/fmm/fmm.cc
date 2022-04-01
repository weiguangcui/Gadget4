/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file fmm.cc
 *
 *  \brief main routines for gravitational force computation with the fast multipole method (FMM)
 */

#include "gadgetconfig.h"

#ifdef FMM

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fmm/fmm.h"
#include "../gravity/ewald.h"
#include "../gravtree/gravtree.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../pm/pm.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

void fmm::fmm_force_passdown(int no, unsigned char no_shmrank, taylor_data taylor_current)
{
  if(no >= MaxPart && no < MaxPart + MaxNodes) /* an internal node  */
    {
      /* first, let's add the external field expansion to the local one accumulated for this node */
#ifdef EVALPOTENTIAL
      taylor_current.coeff.phi += TaylorCoeff[no].coeff.phi;
#endif
      taylor_current.coeff.dphi += TaylorCoeff[no].coeff.dphi;
#if(MULTIPOLE_ORDER >= 2)
      taylor_current.coeff.d2phi += TaylorCoeff[no].coeff.d2phi;
#endif
#if(MULTIPOLE_ORDER >= 3)
      taylor_current.coeff.d3phi += TaylorCoeff[no].coeff.d3phi;
#endif
#if(MULTIPOLE_ORDER >= 4)
      taylor_current.coeff.d4phi += TaylorCoeff[no].coeff.d4phi;
#endif
#if(MULTIPOLE_ORDER >= 5)
      taylor_current.coeff.d5phi += TaylorCoeff[no].coeff.d5phi;
#endif

      taylor_current.coeff.interactions += TaylorCoeff[no].coeff.interactions;
    }
  else
    Terminate("this is not an internal node, which should not happen\n");

  gravnode *node_no = get_nodep(no, no_shmrank);

  int p = node_no->nextnode; /* open cell */

  unsigned char shmrank = node_no->nextnode_shmrank;

  while(p != node_no->sibling || (shmrank != node_no->sibling_shmrank && node_no->sibling >= MaxPart + D->NTopnodes))
    {
      if(p < MaxPart || (p >= ImportedNodeOffset && p < EndOfTreePoints)) /* we have found a single particle */
        {
          if(shmrank != Shmem.Island_ThisTask)
            Terminate("odd");

          int m, mp;
          MyIntPosType *intpos;

          if(p >= ImportedNodeOffset) /* an imported Treepoint particle  */
            {
              m      = p - ImportedNodeOffset;
              intpos = Points[m].IntPos;
              mp     = -1;
              p      = get_nextnodep(shmrank)[p - MaxNodes];
            }
          else
            {
              m      = p;
              intpos = Tp->P[p].IntPos;
              mp     = p;
              p      = get_nextnodep(shmrank)[p];
            }

          /* apply expansion to particle */

          vector<MyReal> dxyz;
          Tp->nearest_image_intpos_to_pos(intpos, node_no->s.da, dxyz.da);

#ifdef EVALPOTENTIAL
          MyReal pot = taylor_current.coeff.phi + taylor_current.coeff.dphi * dxyz;
#endif
          vector<MyReal> dphi = taylor_current.coeff.dphi;

#if(MULTIPOLE_ORDER >= 2)
          vector<MyReal> d2phidxyz = taylor_current.coeff.d2phi * dxyz;
          dphi += d2phidxyz;
#ifdef EVALPOTENTIAL
          pot += 0.5f * (d2phidxyz * dxyz);
#endif
#endif
#if(MULTIPOLE_ORDER >= 3)
          vector<MyReal> d3phidxyz2 = contract_twice(taylor_current.coeff.d3phi, dxyz);
          dphi += 0.5f * d3phidxyz2;
#ifdef EVALPOTENTIAL
          pot += static_cast<MyReal>(1.0 / 6) * (dxyz * d3phidxyz2);
#endif
#endif
#if(MULTIPOLE_ORDER >= 4)
          vector<MyReal> d4phidxyz3 = contract_thrice(taylor_current.coeff.d4phi, dxyz);
          dphi += static_cast<MyReal>(1.0 / 6) * d4phidxyz3;
#ifdef EVALPOTENTIAL
          pot += static_cast<MyReal>(1.0 / 24) * (dxyz * d4phidxyz3);
#endif
#endif
#if(MULTIPOLE_ORDER >= 5)
          vector<MyReal> d5phidxyz4 = contract_fourtimes(taylor_current.coeff.d5phi, dxyz);
          dphi += static_cast<MyReal>(1.0 / 24) * d5phidxyz4;
#ifdef EVALPOTENTIAL
          pot += static_cast<MyReal>(1.0 / 120) * (dxyz * d5phidxyz4);
#endif
#endif

          if(mp >= 0)
            {
#ifndef HIERARCHICAL_GRAVITY
              if(Tp->TimeBinSynchronized[Tp->P[mp].TimeBinGrav])
#endif
                {
                  Tp->P[mp].GravAccel -= dphi;
#ifdef EVALPOTENTIAL
                  Tp->P[mp].Potential += pot;
#endif

                  if(MeasureCostFlag)
                    Tp->P[mp].GravCost += taylor_current.coeff.interactions;

                  interactioncountEffective += taylor_current.coeff.interactions;
                }
            }
          else
            {
#ifndef HIERARCHICAL_GRAVITY
              if(Points[m].ActiveFlag)
#endif
                {
                  int idx = ResultIndexList[m];
                  ResultsActiveImported[idx].GravAccel -= dphi;
#ifdef EVALPOTENTIAL
                  ResultsActiveImported[idx].Potential += pot;
#endif
                  if(MeasureCostFlag)
                    ResultsActiveImported[idx].GravCost += taylor_current.coeff.interactions;

                  interactioncountEffective += taylor_current.coeff.interactions;
                }
            }
        }
      else if(p < MaxPart + MaxNodes) /* an internal node  */
        {
          gravnode *node_p = get_nodep(p, shmrank);

          if(fmm_depends_on_local_mass(p, shmrank))
            {
              taylor_data taylor_sub = taylor_current;

              vector<MyReal> dxyz;
              Tp->nearest_image_intpos_to_pos(node_p->s.da, node_no->s.da, dxyz.da);

              /* now shift the expansion center */

#ifdef EVALPOTENTIAL
              taylor_sub.coeff.phi += taylor_current.coeff.dphi * dxyz;
#endif

#if(MULTIPOLE_ORDER >= 2)
              vector<MyReal> delta_dphi = taylor_current.coeff.d2phi * dxyz;
              taylor_sub.coeff.dphi += delta_dphi;
#ifdef EVALPOTENTIAL
              taylor_sub.coeff.phi += 0.5f * (delta_dphi * dxyz);
#endif
#endif
#if(MULTIPOLE_ORDER >= 3)
              symtensor2<MyReal> delta_d2phi = taylor_current.coeff.d3phi * dxyz;

              taylor_sub.coeff.d2phi += delta_d2phi;

              delta_dphi = delta_d2phi * dxyz;

              taylor_sub.coeff.dphi += 0.5f * delta_dphi;
#ifdef EVALPOTENTIAL
              taylor_sub.coeff.phi += static_cast<MyReal>(1.0 / 6) * (delta_dphi * dxyz);
#endif
#endif
#if(MULTIPOLE_ORDER >= 4)
              symtensor3<MyReal> delta_d3phi = taylor_current.coeff.d4phi * dxyz;
              taylor_sub.coeff.d3phi += delta_d3phi;

              delta_d2phi = delta_d3phi * dxyz;
              taylor_sub.coeff.d2phi += 0.5f * delta_d2phi;

              delta_dphi = delta_d2phi * dxyz;
              taylor_sub.coeff.dphi += static_cast<MyReal>(1.0 / 6) * delta_dphi;
#ifdef EVALPOTENTIAL
              taylor_sub.coeff.phi += static_cast<MyReal>(1.0 / 24) * (delta_dphi * dxyz);
#endif
#endif
#if(MULTIPOLE_ORDER >= 5)
              symtensor4<MyReal> delta_d4phi = taylor_current.coeff.d5phi * dxyz;
              taylor_sub.coeff.d4phi += delta_d4phi;

              delta_d3phi = delta_d4phi * dxyz;
              taylor_sub.coeff.d3phi += 0.5f * delta_d3phi;

              delta_d2phi = delta_d3phi * dxyz;
              taylor_sub.coeff.d2phi += static_cast<MyReal>(1.0 / 6) * delta_d2phi;

              delta_dphi = delta_d2phi * dxyz;
              taylor_sub.coeff.dphi += static_cast<MyReal>(1.0 / 24) * delta_dphi;
#ifdef EVALPOTENTIAL
              taylor_sub.coeff.phi += static_cast<MyReal>(1.0 / 120) * (delta_dphi * dxyz);
#endif
#endif
              fmm_force_passdown(p, shmrank, taylor_sub);
            }

          p       = node_p->sibling;
          shmrank = node_p->sibling_shmrank;
        }
      else if(p >= EndOfTreePoints && p < EndOfForeignNodes) /* an imported tree node */
        {
          Terminate("A");
          gravnode *nop = get_nodep(p, shmrank);
          p             = nop->sibling;
          shmrank       = nop->sibling_shmrank;
        }
      else if(p >= EndOfForeignNodes) /* an imported particle below an imported tree node */
        {
          Terminate("B");
          foreign_gravpoint_data *foreignpoint = get_foreignpointsp(p - EndOfForeignNodes, shmrank);
          p                                    = foreignpoint->Nextnode;
          shmrank                              = foreignpoint->Nextnode_shmrank;
        }
      else
        {
          /* a pseudo point */
          Terminate(
              "should not happen: p=%d MaxPart=%d MaxNodes=%d  ImportedNodeOffset=%d  EndOfTreePoints=%d  EndOfForeignNodes=%d "
              "shmrank=%d",
              p, MaxPart, MaxNodes, ImportedNodeOffset, EndOfTreePoints, EndOfForeignNodes, shmrank);
        }
    }
}

inline void fmm::fmm_open_both(gravnode *noptr_sink, gravnode *noptr_source, int mintopleafnode, int committed)
{
  int self_flag = 0;
  if(noptr_sink == noptr_source)
    self_flag = 1;

  /* open node */
  int p_sink                 = noptr_sink->nextnode;
  unsigned char shmrank_sink = noptr_sink->nextnode_shmrank;

  while(p_sink != noptr_sink->sibling ||
        (shmrank_sink != noptr_sink->sibling_shmrank && noptr_sink->sibling >= MaxPart + D->NTopnodes))
    {
      int next_sink;
      unsigned char next_shmrank_sink;
      char type_sink;

      if(p_sink < MaxPart) /* a local particle */
        {
          /* note: here shmrank cannot change */
          next_sink         = get_nextnodep(shmrank_sink)[p_sink];
          next_shmrank_sink = shmrank_sink;
          type_sink         = NODE_TYPE_LOCAL_PARTICLE;
        }
      else if(p_sink < MaxPart + MaxNodes) /* an internal node  */
        {
          gravnode *nop     = get_nodep(p_sink, shmrank_sink);
          next_sink         = nop->sibling;
          next_shmrank_sink = nop->sibling_shmrank;
          type_sink         = NODE_TYPE_LOCAL_NODE;
        }
      else if(p_sink >= ImportedNodeOffset && p_sink < EndOfTreePoints) /* an imported Treepoint particle */
        {
          /* note: here shmrank cannot change */
          next_sink         = get_nextnodep(shmrank_sink)[p_sink - MaxNodes];
          next_shmrank_sink = shmrank_sink;
          type_sink         = NODE_TYPE_TREEPOINT_PARTICLE;
        }
      else if(p_sink >= EndOfTreePoints && p_sink < EndOfForeignNodes) /* an imported tree node */
        {
          gravnode *nop     = get_nodep(p_sink, shmrank_sink);
          next_sink         = nop->sibling;
          next_shmrank_sink = nop->sibling_shmrank;
          type_sink         = NODE_TYPE_FETCHED_NODE;
        }
      else if(p_sink >= EndOfForeignNodes)
        {
          foreign_gravpoint_data *foreignpoint = get_foreignpointsp(p_sink - EndOfForeignNodes, shmrank_sink);
          next_sink                            = foreignpoint->Nextnode;
          next_shmrank_sink                    = foreignpoint->Nextnode_shmrank;
          type_sink                            = NODE_TYPE_FETCHED_PARTICLE;
        }
      else
        {
          /* a pseudo point */
          next_sink = 0;
          type_sink = 0;
          Terminate("pseudo particle - should not happen");
        }

      int p_source                 = noptr_source->nextnode; /* open cell */
      unsigned char shmrank_source = noptr_source->nextnode_shmrank;

      while(p_source != noptr_source->sibling ||
            (shmrank_source != noptr_source->sibling_shmrank && noptr_source->sibling >= MaxPart + D->NTopnodes))
        {
          int next_source;
          unsigned char next_shmrank_source;
          char type_source;

          if(p_source < MaxPart) /* a local particle */
            {
              /* note: here shmrank cannot change */
              next_source         = get_nextnodep(shmrank_source)[p_source];
              next_shmrank_source = shmrank_source;
              type_source         = NODE_TYPE_LOCAL_PARTICLE;
            }
          else if(p_source < MaxPart + MaxNodes) /* an internal node  */
            {
              gravnode *nop       = get_nodep(p_source, shmrank_source);
              next_source         = nop->sibling;
              next_shmrank_source = nop->sibling_shmrank;
              type_source         = NODE_TYPE_LOCAL_NODE;
            }
          else if(p_source >= ImportedNodeOffset && p_source < EndOfTreePoints) /* an imported Treepoint particle */
            {
              /* note: here shmrank cannot change */
              next_source         = get_nextnodep(shmrank_source)[p_source - MaxNodes];
              next_shmrank_source = shmrank_source;
              type_source         = NODE_TYPE_TREEPOINT_PARTICLE;
            }
          else if(p_source >= EndOfTreePoints && p_source < EndOfForeignNodes) /* an imported tree node */
            {
              gravnode *nop       = get_nodep(p_source, shmrank_source);
              next_source         = nop->sibling;
              next_shmrank_source = nop->sibling_shmrank;
              type_source         = NODE_TYPE_FETCHED_NODE;
            }
          else if(p_source >= EndOfForeignNodes)
            {
              foreign_gravpoint_data *foreignpoint = get_foreignpointsp(p_source - EndOfForeignNodes, shmrank_source);
              next_source                          = foreignpoint->Nextnode;
              next_shmrank_source                  = foreignpoint->Nextnode_shmrank;
              type_source                          = NODE_TYPE_FETCHED_PARTICLE;
            }
          else
            {
              /* a pseudo point */
              next_source = 0;
              type_source = 0;
              Terminate("pseudo particle - should not happen");
            }

          if(self_flag == 0 || p_source >= p_sink)
            {
              if(type_sink >= NODE_TYPE_LOCAL_NODE && type_source <= NODE_TYPE_FETCHED_PARTICLE)
                {
                  /* in this case we have node-particle interaction, which we swap into a particle-node interaction */
                  fmm_force_interact(p_source, p_sink, type_source, type_sink, shmrank_source, shmrank_sink, mintopleafnode,
                                     committed);
                }
              else
                {
                  fmm_force_interact(p_sink, p_source, type_sink, type_source, shmrank_sink, shmrank_source, mintopleafnode,
                                     committed);
                }
            }

          p_source       = next_source;
          shmrank_source = next_shmrank_source;
        }

      p_sink       = next_sink;
      shmrank_sink = next_shmrank_sink;
    }
}

inline void fmm::fmm_open_node(int no_particle, gravnode *nop, char type_particle, unsigned char shmrank_particle, int mintopleafnode,
                               int committed)
{
  int p                 = nop->nextnode;
  unsigned char shmrank = nop->nextnode_shmrank;

  while(p != nop->sibling || (shmrank != nop->sibling_shmrank && nop->sibling >= MaxPart + D->NTopnodes))
    {
      if(p < 0)
        Terminate("p=%d < 0", p);

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
          gravnode *nop = get_nodep(p, shmrank);
          next          = nop->sibling;
          next_shmrank  = nop->sibling_shmrank;
          type          = NODE_TYPE_LOCAL_NODE;
        }
      else if(p >= ImportedNodeOffset && p < EndOfTreePoints) /* an imported Treepoint particle  */
        {
          /* note: here shmrank cannot change */
          next         = get_nextnodep(shmrank)[p - MaxNodes];
          next_shmrank = shmrank;
          type         = NODE_TYPE_TREEPOINT_PARTICLE;
        }
      else if(p >= EndOfTreePoints && p < EndOfForeignNodes) /* an imported tree node */
        {
          gravnode *nop = get_nodep(p, shmrank);
          next          = nop->sibling;
          next_shmrank  = nop->sibling_shmrank;
          type          = NODE_TYPE_FETCHED_NODE;
        }
      else if(p >= EndOfForeignNodes) /* an imported particle below an imported tree node */
        {
          foreign_gravpoint_data *foreignpoint = get_foreignpointsp(p - EndOfForeignNodes, shmrank);

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

      fmm_force_interact(no_particle, p, type_particle, type, shmrank_particle, shmrank, mintopleafnode, committed);

      p       = next;
      shmrank = next_shmrank;
    }
}

inline void fmm::fmm_particle_particle_interaction(int no_sink, int no_source, int type_sink, int type_source,
                                                   unsigned char shmrank_sink, unsigned char shmrank_source)
{
#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  if(skip_actual_force_computation)
    return;
#endif

  MyIntPosType *intpos_n, *intpos_m;
  MyReal mass_n, mass_m;
  MyReal h_n, h_m;
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  int test_n, test_m;
#endif

  /* in one of the following three cases we have a single particle on the sink side */
  if(type_sink == NODE_TYPE_LOCAL_PARTICLE)
    {
      particle_data *P = get_Pp(no_sink, shmrank_sink);

      intpos_n = P->IntPos;
      mass_n   = P->getMass();
      h_n      = All.ForceSoftening[P->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_n = P->InsideOutsideFlag;
#endif
    }
  else if(type_sink == NODE_TYPE_TREEPOINT_PARTICLE)
    {
      gravpoint_data *Pointp = get_pointsp(no_sink - ImportedNodeOffset, shmrank_sink);

      intpos_n = Pointp->IntPos;
      mass_n   = Pointp->Mass;
      h_n      = All.ForceSoftening[Pointp->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_n = Pointp->InsideOutsideFlag;
#endif
    }
  else /* a point that was fetched */
    {
      foreign_gravpoint_data *foreignpoint = get_foreignpointsp(no_sink - EndOfForeignNodes, shmrank_sink);

      intpos_n = foreignpoint->IntPos;
      mass_n   = foreignpoint->Mass;
      h_n      = All.ForceSoftening[foreignpoint->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_n = foreignpoint->InsideOutsideFlag;
#endif
    }

  /* in one of the following three cases we have a single particle on the source side */
  if(type_source == NODE_TYPE_LOCAL_PARTICLE)
    {
      particle_data *P = get_Pp(no_source, shmrank_source);

      intpos_m = P->IntPos;
      mass_m   = P->getMass();
      h_m      = All.ForceSoftening[P->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_m = P->InsideOutsideFlag;
#endif
    }
  else if(type_source == NODE_TYPE_TREEPOINT_PARTICLE)
    {
      gravpoint_data *Pointp = get_pointsp(no_source - ImportedNodeOffset, shmrank_source);

      intpos_m = Pointp->IntPos;
      mass_m   = Pointp->Mass;
      h_m      = All.ForceSoftening[Pointp->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_m = Pointp->InsideOutsideFlag;
#endif
    }
  else
    {
      foreign_gravpoint_data *foreignpoint = get_foreignpointsp(no_source - EndOfForeignNodes, shmrank_source);

      intpos_m = foreignpoint->IntPos;
      mass_m   = foreignpoint->Mass;
      h_m      = All.ForceSoftening[foreignpoint->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_m = foreignpoint->InsideOutsideFlag;
#endif
    }

  MyReal h_max = (h_m > h_n) ? h_m : h_n;

  vector<MyReal> dxyz;
  Tp->nearest_image_intpos_to_pos(intpos_m, intpos_n, dxyz.da); /* converts the integer distance to floating point */

  MyReal r2   = dxyz.r2();
  MyReal r    = sqrt(r2);
  MyReal rinv = (r > 0) ? 1 / r : 0;

  gfactors gfac;

#ifdef PMGRID
  if(DoPM)
    {
      mesh_factors *mfp = &mf[LOW_MESH];

#ifdef PLACEHIGHRESREGION
      if((DoPM & TREE_ACTIVE_CUTTOFF_HIGHRES_PM))
        {
          if(test_m == FLAG_INSIDE && test_n == FLAG_INSIDE)
            mfp = &mf[HIGH_MESH];
        }
#endif
      if(modify_gfactors_pm_monopole(gfac, r, rinv, mfp)) /* if we are outside the cut-off radius, we have no interaction */
        return;
    }
#endif

  get_gfactors_monopole(gfac, r, h_max, rinv);

#ifdef EVALPOTENTIAL
  MyReal D0 = gfac.fac0;
#endif
  vector<MyReal> D1 = gfac.fac1 * rinv * dxyz;

  if(DoEwald)
    {
      ewald_data ew;
      Ewald.ewald_gridlookup(intpos_m, intpos_n, ewald::POINTMASS, ew);

      D1 -= ew.D1phi;

#ifdef EVALPOTENTIAL
      D0 -= ew.D0phi;
#endif
    }

  if(shmrank_sink == Shmem.Island_ThisTask)
    {
      if(type_sink == NODE_TYPE_LOCAL_PARTICLE)
        {
#ifndef HIERARCHICAL_GRAVITY
          if(Tp->TimeBinSynchronized[Tp->P[no_sink].TimeBinGrav])
#endif
            {
              Tp->P[no_sink].GravAccel -= mass_m * D1;
#ifdef EVALPOTENTIAL
              Tp->P[no_sink].Potential -= mass_m * D0;
#endif
              if(MeasureCostFlag)
                Tp->P[no_sink].GravCost++;

              interactioncountPP += 1;
            }
        }
      else if(type_sink == NODE_TYPE_TREEPOINT_PARTICLE)
        {
#ifndef HIERARCHICAL_GRAVITY
          if(Points[no_sink - ImportedNodeOffset].ActiveFlag)
#endif
            {
              int idx = ResultIndexList[no_sink - ImportedNodeOffset];
              ResultsActiveImported[idx].GravAccel -= mass_m * D1;
#ifdef EVALPOTENTIAL
              ResultsActiveImported[idx].Potential -= mass_m * D0;
#endif
              if(MeasureCostFlag)
                ResultsActiveImported[idx].GravCost++;

              interactioncountPP += 1;
            }
        }
    }

  if(shmrank_source == Shmem.Island_ThisTask)
    {
      if(type_source == NODE_TYPE_LOCAL_PARTICLE)
        {
#ifndef HIERARCHICAL_GRAVITY
          if(Tp->TimeBinSynchronized[Tp->P[no_source].TimeBinGrav])
#endif
            {
              Tp->P[no_source].GravAccel += mass_n * D1;
#ifdef EVALPOTENTIAL
              Tp->P[no_source].Potential -= mass_n * D0;
#endif

              if(MeasureCostFlag)
                Tp->P[no_source].GravCost++;

              interactioncountPP += 1;
            }
        }
      else if(type_source == NODE_TYPE_TREEPOINT_PARTICLE)
        {
#ifndef HIERARCHICAL_GRAVITY
          if(Points[no_source - ImportedNodeOffset].ActiveFlag)
#endif
            {
              int idx = ResultIndexList[no_source - ImportedNodeOffset];
              ResultsActiveImported[idx].GravAccel += mass_n * D1;
#ifdef EVALPOTENTIAL
              ResultsActiveImported[idx].Potential -= mass_n * D0;
#endif
              if(MeasureCostFlag)
                ResultsActiveImported[idx].GravCost++;

              interactioncountPP += 1;
            }
        }
    }
}

inline void fmm::fmm_particle_node_interaction(int no_sink, int no_source, int type_sink, int type_source, unsigned char shmrank_sink,
                                               unsigned char shmrank_source, gravnode *noptr_source, vector<MyReal> &dxyz, MyReal &r2)
{
#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  if(skip_actual_force_computation)
    return;
#endif

  /* 'sink' is a particle
   * 'source' node is a node with multipole moments.
   * 'dxyz' is the distance vector, pointing from sink to source, i.e. dxyz = pos(source) - pos(sink)
   */

  MyReal mass_i, h_i;
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  int test_point;
#endif

  MyIntPosType *intpos_i;

  if(type_sink == NODE_TYPE_LOCAL_PARTICLE)
    {
      particle_data *P = get_Pp(no_sink, shmrank_sink);

      intpos_i = P->IntPos;
      mass_i   = P->getMass();
      h_i      = All.ForceSoftening[P->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_point = P->InsideOutsideFlag;
#endif
    }
  else if(type_sink == NODE_TYPE_TREEPOINT_PARTICLE)
    {
      gravpoint_data *Pointp = get_pointsp(no_sink - ImportedNodeOffset, shmrank_sink);

      intpos_i = Pointp->IntPos;
      mass_i   = Pointp->Mass;
      h_i      = All.ForceSoftening[Pointp->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_point = Pointp->InsideOutsideFlag;
#endif
    }
  else /* a point that was fetched */
    {
      foreign_gravpoint_data *foreignpoint = get_foreignpointsp(no_sink - EndOfForeignNodes, shmrank_sink);

      intpos_i = foreignpoint->IntPos;
      mass_i   = foreignpoint->Mass;
      h_i      = All.ForceSoftening[foreignpoint->getSofteningClass()];
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_point = foreignpoint->InsideOutsideFlag;
#endif
    }

  MyReal h_j   = All.ForceSoftening[noptr_source->getSofteningClass()];
  MyReal h_max = (h_j > h_i) ? h_j : h_i;

  /* do cell-particle interaction, node can be used */
  MyReal r = sqrt(r2);

  MyReal rinv = (r > 0) ? 1 / r : 0;

  gfactors gfac;

#ifdef PMGRID
  if(DoPM)
    {
      mesh_factors *mfp = &mf[LOW_MESH];

#ifdef PLACEHIGHRESREGION
      if((DoPM & TREE_ACTIVE_CUTTOFF_HIGHRES_PM))
        {
          int test_node = noptr_source->overlap_flag;
          if(test_node == FLAG_INSIDE && test_point == FLAG_INSIDE)
            mfp = &mf[HIGH_MESH];
        }
#endif

      if(modify_gfactors_pm_multipole(gfac, r, rinv, mfp)) /* if we are outside the cut-off radius, we have no interaction */
        return;
    }
#endif

  get_gfactors_multipole(gfac, r, h_max, rinv);

#ifdef EVALPOTENTIAL
  MyReal g0 = gfac.fac0;
  MyReal D0 = g0;
#endif

  MyReal g1         = gfac.fac1 * rinv;
  vector<MyReal> D1 = g1 * dxyz;

#if(MULTIPOLE_ORDER >= 2)
  MyReal g2               = gfac.fac2 * gfac.rinv2;
  symtensor2<MyReal> aux2 = dxyz % dxyz;  // construct outer product of the two vectors
  symtensor2<MyReal> D2   = g2 * aux2;
  D2[qXX] += g1;
  D2[qYY] += g1;
  D2[qZZ] += g1;
#endif
#if(MULTIPOLE_ORDER >= 3)
  MyReal g3 = gfac.fac3 * gfac.rinv3;
  symtensor3<MyReal> D3;
  symtensor3<MyReal> aux3;
  setup_D3(INIT, D3, dxyz, aux2, aux3, g2, g3);
#endif

#if(MULTIPOLE_ORDER >= 4)
  MyReal g4 = gfac.fac4 * gfac.rinv2 * gfac.rinv2;
  symtensor4<MyReal> D4;
  symtensor4<MyReal> aux4;
  setup_D4(INIT, D4, dxyz, aux2, aux3, aux4, g2, g3, g4);
#endif

#if(MULTIPOLE_ORDER >= 5)
  MyReal g5 = gfac.fac5 * gfac.rinv3 * gfac.rinv2;
  symtensor5<MyReal> D5;
  symtensor5<MyReal> aux5;
  setup_D5(INIT, D5, dxyz, aux3, aux4, aux5, g3, g4, g5);
#endif

  if(DoEwald)
    {
      ewald_data ew;
      Ewald.ewald_gridlookup(noptr_source->s.da, intpos_i, ewald::MULTIPOLES, ew);

#ifdef EVALPOTENTIAL
      D0 -= ew.D0phi;
#endif
      D1 -= ew.D1phi;
#if(MULTIPOLE_ORDER >= 2)
      D2 -= ew.D2phi;
#endif
#if(MULTIPOLE_ORDER >= 3)
      D3 -= ew.D3phi;
#endif
#if(MULTIPOLE_ORDER >= 4)
      D4 -= ew.D4phi;
#endif
#if(MULTIPOLE_ORDER >= 5)
      D5 -= ew.D5phi;
#endif
    }

  /* finally store the force on the particle */
  if(shmrank_sink == Shmem.Island_ThisTask)
    if(type_sink == NODE_TYPE_LOCAL_PARTICLE || type_sink == NODE_TYPE_TREEPOINT_PARTICLE)
      {
        MyReal mass_j = noptr_source->mass;

#if(MULTIPOLE_ORDER >= 3) || ((MULTIPOLE_ORDER >= 2) && defined(EXTRAPOTTERM))
        symtensor2<MyDouble> &Q2_j = noptr_source->Q2Tensor;
#endif
#if(MULTIPOLE_ORDER >= 4) || ((MULTIPOLE_ORDER >= 3) && defined(EXTRAPOTTERM))
        symtensor3<MyDouble> &Q3_j = noptr_source->Q3Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5) || ((MULTIPOLE_ORDER >= 4) && defined(EXTRAPOTTERM))
        symtensor4<MyDouble> &Q4_j = noptr_source->Q4Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5) && defined(EXTRAPOTTERM)
        symtensor5<MyDouble> &Q5_j = noptr_source->Q5Tensor;
#endif

#ifdef EVALPOTENTIAL
        MyReal pot = -mass_j * D0;
#if(MULTIPOLE_ORDER >= 3) || ((MULTIPOLE_ORDER >= 2) && defined(EXTRAPOTTERM))
        pot -= 0.5f * (D2 * Q2_j);
#endif
#if(MULTIPOLE_ORDER >= 4) || ((MULTIPOLE_ORDER >= 3) && defined(EXTRAPOTTERM))
        pot -= static_cast<MyReal>(1.0 / 6) * (D3 * Q3_j);
#endif
#if(MULTIPOLE_ORDER >= 5) || ((MULTIPOLE_ORDER >= 4) && defined(EXTRAPOTTERM))
        pot -= static_cast<MyReal>(1.0 / 24) * (D4 * Q4_j);
#endif
#if(MULTIPOLE_ORDER >= 5) && defined(EXTRAPOTTERM)

        pot -= static_cast<MyReal>(1.0 / 120) * (D5 * Q5_j);
#endif
#endif

        vector<MyReal> dphi = mass_j * D1;

#if(MULTIPOLE_ORDER >= 3)
        dphi += static_cast<MyReal>(0.5) * (D3 * Q2_j);
#endif
#if(MULTIPOLE_ORDER >= 4)
        dphi += static_cast<MyReal>(1.0 / 6) * (D4 * Q3_j);
#endif
#if(MULTIPOLE_ORDER >= 5)
        dphi += static_cast<MyReal>(1.0 / 24) * (D5 * Q4_j);
#endif

        if(type_sink == NODE_TYPE_LOCAL_PARTICLE)
          {
#ifndef HIERARCHICAL_GRAVITY
            if(Tp->TimeBinSynchronized[Tp->P[no_sink].TimeBinGrav])
#endif
              {
                Tp->P[no_sink].GravAccel -= dphi;
#ifdef EVALPOTENTIAL
                Tp->P[no_sink].Potential += pot;
#endif
                if(MeasureCostFlag)
                  Tp->P[no_sink].GravCost++;

                interactioncountPN += 1;
                interactioncountEffective += 1;
              }
          }
        else
          {
#ifndef HIERARCHICAL_GRAVITY
            if(Points[no_sink - ImportedNodeOffset].ActiveFlag)
#endif
              {
                int idx = ResultIndexList[no_sink - ImportedNodeOffset];
                ResultsActiveImported[idx].GravAccel -= dphi;
#ifdef EVALPOTENTIAL
                ResultsActiveImported[idx].Potential += pot;
#endif
                if(MeasureCostFlag)
                  ResultsActiveImported[idx].GravCost++;

                interactioncountPN += 1;
                interactioncountEffective += 1;
              }
          }
      }

  if(fmm_depends_on_local_mass(no_source, shmrank_source))
    if(type_source == NODE_TYPE_LOCAL_NODE)
      {
        /* mediating field expansion of point particle on source node */
#ifdef EVALPOTENTIAL
        TaylorCoeff[no_source].coeff.phi += (-mass_i) * D0;
#endif
        TaylorCoeff[no_source].coeff.dphi += (-mass_i) * D1;
#if(MULTIPOLE_ORDER >= 2)
        TaylorCoeff[no_source].coeff.d2phi += (-mass_i) * D2;
#endif
#if(MULTIPOLE_ORDER >= 3)
        TaylorCoeff[no_source].coeff.d3phi += (-mass_i) * D3;
#endif
#if(MULTIPOLE_ORDER >= 4)
        TaylorCoeff[no_source].coeff.d4phi += (-mass_i) * D4;
#endif
#if(MULTIPOLE_ORDER >= 5)
        TaylorCoeff[no_source].coeff.d5phi += (-mass_i) * D5;
#endif
        TaylorCoeff[no_source].coeff.interactions += 1;

        interactioncountPN += 1;
      }
}

inline void fmm::fmm_node_node_interaction(int no_sink, int no_source, int type_sink, int type_source, unsigned char shmrank_sink,
                                           unsigned char shmrank_source, gravnode *noptr_sink, gravnode *noptr_source,
                                           vector<MyReal> &dxyz, MyReal &r2)
{
#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  if(skip_actual_force_computation)
    return;
#endif

  /* 'sink' is a node with multipole moments
   * 'source' node is a node with multipole moments
   * 'dxyz' is the distance vector, pointing from sink to source, i.e. dxyz = pos(source) - pos(sink)
   */

  /* now do the node-node interaction */
  MyReal r = sqrt(r2);

#if(MULTIPOLE_ORDER >= 3) || ((MULTIPOLE_ORDER >= 2) && defined(EXTRAPOTTERM))
  symtensor2<MyDouble> &Q2_m = noptr_source->Q2Tensor;
  symtensor2<MyDouble> &Q2_n = noptr_sink->Q2Tensor;
#endif
#if(MULTIPOLE_ORDER >= 4) || ((MULTIPOLE_ORDER >= 3) && defined(EXTRAPOTTERM))
  symtensor3<MyDouble> &Q3_m = noptr_source->Q3Tensor;
  symtensor3<MyDouble> &Q3_n = noptr_sink->Q3Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5) || ((MULTIPOLE_ORDER >= 4) && defined(EXTRAPOTTERM))
  symtensor4<MyDouble> &Q4_m = noptr_source->Q4Tensor;
  symtensor4<MyDouble> &Q4_n = noptr_sink->Q4Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5) && defined(EXTRAPOTTERM) && defined(EVALPOTENTIAL)
  symtensor5<MyDouble> &Q5_m = noptr_source->Q5Tensor;
  symtensor5<MyDouble> &Q5_n = noptr_sink->Q5Tensor;
#endif

  MyReal mass_m = noptr_source->mass;
  MyReal mass_n = noptr_sink->mass;

  MyReal rinv = (r > 0) ? 1 / r : 0;

  MyReal h_n   = All.ForceSoftening[noptr_sink->getSofteningClass()];
  MyReal h_m   = All.ForceSoftening[noptr_source->getSofteningClass()];
  MyReal h_max = (h_m > h_n) ? h_m : h_n;

  gfactors gfac;

#ifdef PMGRID
  if(DoPM)
    {
      mesh_factors *mfp = &mf[LOW_MESH];

#ifdef PLACEHIGHRESREGION
      if((DoPM & TREE_ACTIVE_CUTTOFF_HIGHRES_PM))
        {
          if(noptr_source->overlap_flag == FLAG_INSIDE && noptr_sink->overlap_flag == FLAG_INSIDE)
            mfp = &mf[HIGH_MESH];
        }
#endif

      if(modify_gfactors_pm_multipole(gfac, r, rinv, mfp)) /* if we are outside the cut-off radius, we have no interaction */
        return;
    }
#endif

  get_gfactors_multipole(gfac, r, h_max, rinv);

#ifdef EVALPOTENTIAL
  MyReal g0 = gfac.fac0;
  MyReal D0 = g0;
#endif

  MyReal g1         = gfac.fac1 * rinv;
  vector<MyReal> D1 = g1 * dxyz;

#if(MULTIPOLE_ORDER >= 2)
  MyReal g2               = gfac.fac2 * gfac.rinv2;
  symtensor2<MyReal> aux2 = dxyz % dxyz;  // construct outer product of the two vectors
  symtensor2<MyReal> D2   = g2 * aux2;
  D2[qXX] += g1;
  D2[qYY] += g1;
  D2[qZZ] += g1;
#endif

#if(MULTIPOLE_ORDER >= 3)
  MyReal g3 = gfac.fac3 * gfac.rinv3;
  symtensor3<MyReal> D3;
  symtensor3<MyReal> aux3;
  setup_D3(INIT, D3, dxyz, aux2, aux3, g2, g3);
#endif

#if(MULTIPOLE_ORDER >= 4)
  MyReal g4 = gfac.fac4 * gfac.rinv2 * gfac.rinv2;
  symtensor4<MyReal> D4;
  symtensor4<MyReal> aux4;
  setup_D4(INIT, D4, dxyz, aux2, aux3, aux4, g2, g3, g4);
#endif

#if(MULTIPOLE_ORDER >= 5)
  MyReal g5 = gfac.fac5 * gfac.rinv3 * gfac.rinv2;
  symtensor5<MyReal> D5;
  symtensor5<MyReal> aux5;
  setup_D5(INIT, D5, dxyz, aux3, aux4, aux5, g3, g4, g5);
#endif

  if(DoEwald)
    {
      ewald_data ew;
      Ewald.ewald_gridlookup(noptr_source->s.da, noptr_sink->s.da, ewald::MULTIPOLES, ew);

#ifdef EVALPOTENTIAL
      D0 -= ew.D0phi;
#endif
      D1 -= ew.D1phi;
#if(MULTIPOLE_ORDER >= 2)
      D2 -= ew.D2phi;
#endif
#if(MULTIPOLE_ORDER >= 3)
      D3 -= ew.D3phi;
#endif
#if(MULTIPOLE_ORDER >= 4)
      D4 -= ew.D4phi;
#endif
#if(MULTIPOLE_ORDER >= 5)
      D5 -= ew.D5phi;
#endif
    }

  if(fmm_depends_on_local_mass(no_sink, shmrank_sink))
    if(type_sink == NODE_TYPE_LOCAL_NODE)
      {
#ifdef EVALPOTENTIAL
        TaylorCoeff[no_sink].coeff.phi += (-mass_m) * D0;
#if(MULTIPOLE_ORDER >= 3) || ((MULTIPOLE_ORDER >= 2) && defined(EXTRAPOTTERM))
        TaylorCoeff[no_sink].coeff.phi += static_cast<MyReal>(-0.5) * (D2 * Q2_m);
#endif
#if(MULTIPOLE_ORDER >= 4) || ((MULTIPOLE_ORDER >= 3) && defined(EXTRAPOTTERM))
        TaylorCoeff[no_sink].coeff.phi += static_cast<MyReal>(-1.0 / 6) * (D3 * Q3_m);
#endif
#if(MULTIPOLE_ORDER >= 5) || ((MULTIPOLE_ORDER >= 4) && defined(EXTRAPOTTERM))
        TaylorCoeff[no_sink].coeff.phi += static_cast<MyReal>(-1.0 / 24) * (D4 * Q4_m);
#endif
#if(MULTIPOLE_ORDER >= 5) && defined(EXTRAPOTTERM)
        TaylorCoeff[no_sink].coeff.phi += static_cast<MyReal>(-1.0 / 120) * (D5 * Q5_m);
#endif
#endif

        TaylorCoeff[no_sink].coeff.dphi += mass_m * D1;
#if(MULTIPOLE_ORDER >= 2)
        TaylorCoeff[no_sink].coeff.d2phi += (-mass_m) * D2;
#endif
#if(MULTIPOLE_ORDER >= 3)
        TaylorCoeff[no_sink].coeff.dphi += static_cast<MyReal>(0.5) * (D3 * Q2_m);
        TaylorCoeff[no_sink].coeff.d3phi += mass_m * D3;
#endif
#if(MULTIPOLE_ORDER >= 4)
        TaylorCoeff[no_sink].coeff.dphi += static_cast<MyReal>(1.0 / 6) * (D4 * Q3_m);
        TaylorCoeff[no_sink].coeff.d2phi += static_cast<MyReal>(-0.5) * (D4 * Q2_m);
        TaylorCoeff[no_sink].coeff.d4phi += (-mass_m) * D4;
#endif
#if(MULTIPOLE_ORDER >= 5)
        TaylorCoeff[no_sink].coeff.dphi += static_cast<MyReal>(1.0 / 24) * (D5 * Q4_m);
        TaylorCoeff[no_sink].coeff.d2phi += static_cast<MyReal>(-1.0 / 6) * (D5 * Q3_m);
        TaylorCoeff[no_sink].coeff.d3phi += static_cast<MyReal>(0.5) * (D5 * Q2_m);
        TaylorCoeff[no_sink].coeff.d5phi += mass_m * D5;
#endif

        TaylorCoeff[no_sink].coeff.interactions += 1;
        interactioncountNN += 1;
      }

  if(fmm_depends_on_local_mass(no_source, shmrank_source))
    if(type_source == NODE_TYPE_LOCAL_NODE)
      {
#ifdef EVALPOTENTIAL
        TaylorCoeff[no_source].coeff.phi += (-mass_n) * D0;
#if(MULTIPOLE_ORDER >= 3) || ((MULTIPOLE_ORDER >= 2) && defined(EXTRAPOTTERM))
        TaylorCoeff[no_source].coeff.phi += static_cast<MyReal>(-0.5) * (D2 * Q2_n);
#endif
#if(MULTIPOLE_ORDER >= 4) || ((MULTIPOLE_ORDER >= 3) && defined(EXTRAPOTTERM))
        TaylorCoeff[no_source].coeff.phi += static_cast<MyReal>(1.0 / 6) * (D3 * Q3_n);
#endif
#if(MULTIPOLE_ORDER >= 5) || ((MULTIPOLE_ORDER >= 4) && defined(EXTRAPOTTERM))
        TaylorCoeff[no_source].coeff.phi += static_cast<MyReal>(-1.0 / 24) * (D4 * Q4_n);
#endif
#if(MULTIPOLE_ORDER >= 5) && defined(EXTRAPOTTERM)
        TaylorCoeff[no_source].coeff.phi += static_cast<MyReal>(1.0 / 120) * (D5 * Q5_n);
#endif
#endif

        TaylorCoeff[no_source].coeff.dphi += (-mass_n) * D1;
#if(MULTIPOLE_ORDER >= 2)
        TaylorCoeff[no_source].coeff.d2phi += (-mass_n) * D2;
#endif
#if(MULTIPOLE_ORDER >= 3)
        TaylorCoeff[no_source].coeff.dphi += static_cast<MyReal>(-0.5) * (D3 * Q2_n);
        TaylorCoeff[no_source].coeff.d3phi += (-mass_n) * D3;
#endif
#if(MULTIPOLE_ORDER >= 4)
        TaylorCoeff[no_source].coeff.dphi += static_cast<MyReal>(1.0 / 6) * (D4 * Q3_n);
        TaylorCoeff[no_source].coeff.d2phi += static_cast<MyReal>(-0.5) * (D4 * Q2_n);
        TaylorCoeff[no_source].coeff.d4phi += (-mass_n) * D4;
#endif
#if(MULTIPOLE_ORDER >= 5)
        TaylorCoeff[no_source].coeff.dphi += static_cast<MyReal>(-1.0 / 24) * (D5 * Q4_n);
        TaylorCoeff[no_source].coeff.d2phi += static_cast<MyReal>(1.0 / 6) * (D5 * Q3_n);
        TaylorCoeff[no_source].coeff.d3phi += static_cast<MyReal>(-0.5) * (D5 * Q2_n);
        TaylorCoeff[no_source].coeff.d5phi += (-mass_n) * D5;
#endif

        TaylorCoeff[no_source].coeff.interactions += 1;
        interactioncountNN += 1;
      }
}

inline int fmm::fmm_evaluate_node_node_opening_criterion(gravnode *noptr_sink, gravnode *noptr_source, vector<MyReal> &dxyz,
                                                         MyReal &r2)
{
  if(noptr_source->level != noptr_sink->level)
    Terminate("This shouldn't happen: noptr_level=%d   noptr_sink->level=%d ", noptr_source->level, noptr_sink->level);

  if(noptr_source->level <=
     1)  // always open the root node, and the next level (note: full node length does not fit in the integer type)
    return NODE_OPEN;

  /* Note: we will always have noptr_sink->len == noptr_source->len in our algorithm! */

  MyIntPosType halflen = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - noptr_sink->level);
  MyIntPosType intlen  = halflen << 1;

#ifndef TREE_NO_SAFETY_BOX
  // We always open adjacent nodes to protect against worst-case force errors
  MyIntPosType twolens = intlen + (intlen - 1);
  MyIntPosType dist[3];
  Tp->nearest_image_intpos_to_absolute_intdist(noptr_source->center.da, noptr_sink->center.da, dist);

  if(dist[0] < twolens && dist[1] < twolens && dist[2] < twolens)
    return NODE_OPEN;
#endif

  /* converts the integer distance of the centers of mass to floating point */
  Tp->nearest_image_intpos_to_pos(noptr_source->s.da, noptr_sink->s.da, dxyz.da);

  r2 = dxyz.r2();

#ifdef PMGRID
  mesh_factors *mfp = &mf[LOW_MESH];

#ifdef PLACEHIGHRESREGION
  if((DoPM & TREE_ACTIVE_CUTTOFF_HIGHRES_PM))
    {
      int test_source = noptr_source->overlap_flag;
      int test_sink   = noptr_sink->overlap_flag;

      if((test_source == FLAG_BOUNDARYOVERLAP && test_sink != FLAG_OUTSIDE) ||
         (test_sink == FLAG_BOUNDARYOVERLAP && test_source != FLAG_OUTSIDE))
        {
          Terminate("this shouldn't happen any more");
          return NODE_OPEN;
        }
      else
        {
          if(test_source == FLAG_INSIDE && test_sink == FLAG_INSIDE)
            mfp = &mf[HIGH_MESH];
        }
    }
#endif

  if(DoPM && r2 > mfp->rcut2 && noptr_sink->level > 1)
    {
      /* check whether we can ignore any interactions along this branch */
      MyIntPosType dist_x = noptr_source->center[0] - noptr_sink->center[0];
      dist_x              = (((MySignedIntPosType)dist_x) >= 0) ? dist_x : -dist_x;
      if(dist_x > mfp->intrcut[0] + intlen)
        return NODE_DISCARD;

      MyIntPosType dist_y = noptr_source->center[1] - noptr_sink->center[1];
      dist_y              = (((MySignedIntPosType)dist_y) >= 0) ? dist_y : -dist_y;
      if(dist_y > mfp->intrcut[1] + intlen)
        return NODE_DISCARD;

      MyIntPosType dist_z = noptr_source->center[2] - noptr_sink->center[2];
      dist_z              = (((MySignedIntPosType)dist_z) >= 0) ? dist_z : -dist_z;
      if(dist_z > mfp->intrcut[2] + intlen)
        return NODE_DISCARD;
    }
#endif

  /* evaluate generalized opening criterion */

  MyReal len  = intlen * Tp->FacIntToCoord;
  MyReal len2 = len * len;

  if(All.RelOpeningCriterionInUse == 0) /* check Barnes-Hut opening criterion */
    {
      if(4 * len2 > r2 * errTolTheta2)
        return NODE_OPEN;
    }
  else
    {
      if(4 * len2 > r2 * errTolThetaMax2)
        return NODE_OPEN;

      MyReal mmax = (noptr_source->mass < noptr_sink->mass) ? noptr_sink->mass : noptr_source->mass;
      MyReal amin =
          errTolForceAcc * ((noptr_sink->MinOldAcc < noptr_source->MinOldAcc) ? noptr_sink->MinOldAcc : noptr_source->MinOldAcc);
#if(MULTIPOLE_ORDER == 1)
      if(mmax > r2 * amin)
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 2)
      if(square(mmax * len) > r2 * square(r2 * amin))
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 3)
      if(mmax * len2 > r2 * r2 * amin)
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 4)
      if(square(mmax * len * len2) > r2 * square(r2 * r2 * amin))
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 5)
      if(mmax * len2 * len2 > r2 * r2 * r2 * amin)
        return NODE_OPEN;
#endif
    }

#if NSOFTCLASSES > 1
  MyReal h_m   = All.ForceSoftening[noptr_sink->getSofteningClass()];
  MyReal h_n   = All.ForceSoftening[noptr_source->getSofteningClass()];
  MyReal h_max = (h_m > h_n) ? h_m : h_n;

  if(r2 < h_max * h_max)
    {
      if(All.ForceSoftening[noptr_source->minsofttype] < All.ForceSoftening[noptr_source->maxsofttype] &&
         h_max > All.ForceSoftening[noptr_sink->minsofttype])
        return NODE_OPEN;
      else if(All.ForceSoftening[noptr_sink->minsofttype] < All.ForceSoftening[noptr_sink->maxsofttype] &&
              h_max > All.ForceSoftening[noptr_source->minsofttype])
        return NODE_OPEN;
    }
#endif

  return NODE_USE;
}

inline int fmm::fmm_evaluate_particle_node_opening_criterion(int no_sink, char type_sink, unsigned char shmrank_sink,
                                                             gravnode *nop_source, vector<MyReal> &dxyz, MyReal &r2)
{
  MyIntPosType *intpos_n;
  MyReal mass_n;
  MyReal aold;
#if NSOFTCLASSES > 1
  MyReal h_n;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  int test_point;
#endif

  if(type_sink == NODE_TYPE_LOCAL_PARTICLE)
    {
      particle_data *P = get_Pp(no_sink, shmrank_sink);

      intpos_n = P->IntPos;
      mass_n   = P->getMass();
      aold     = P->OldAcc;
#if NSOFTCLASSES > 1
      h_n = All.ForceSoftening[P->getSofteningClass()];
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_point = P->InsideOutsideFlag;
#endif
    }
  else if(type_sink == NODE_TYPE_TREEPOINT_PARTICLE)
    {
      gravpoint_data *Pointp = get_pointsp(no_sink - ImportedNodeOffset, shmrank_sink);

      intpos_n = Pointp->IntPos;
      mass_n   = Pointp->Mass;
      aold     = Pointp->OldAcc;
#if NSOFTCLASSES > 1
      h_n = All.ForceSoftening[Pointp->getSofteningClass()];
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_point = Pointp->InsideOutsideFlag;
#endif
    }
  else /* a point that was fetched */
    {
      foreign_gravpoint_data *foreignpoint = get_foreignpointsp(no_sink - EndOfForeignNodes, shmrank_sink);

      intpos_n = foreignpoint->IntPos;
      mass_n   = foreignpoint->Mass;
      aold     = foreignpoint->OldAcc;
#if NSOFTCLASSES > 1
      h_n = All.ForceSoftening[foreignpoint->getSofteningClass()];
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      test_point = foreignpoint->InsideOutsideFlag;
#endif
    }

  if(nop_source->level <= LEVEL_ALWAYS_OPEN)  // always open the root node (note: full node length does not fit in the integer type)
    return NODE_OPEN;

  MyIntPosType halflen = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - nop_source->level);
  MyIntPosType intlen  = halflen << 1;

#ifndef TREE_NO_SAFETY_BOX
  MyIntPosType dist[3];
  Tp->nearest_image_intpos_to_absolute_intdist(nop_source->center.da, intpos_n, dist);
  // if we are close to the node, and therefore open it to protect against worst-case force errors
  if(dist[0] < intlen && dist[1] < intlen && dist[2] < intlen)
    return NODE_OPEN;
#endif

  /* check a variant of the classic opening criterion */
  Tp->nearest_image_intpos_to_pos(nop_source->s.da, intpos_n, dxyz.da); /* converts the integer distance to floating point */

  r2 = dxyz.r2();

#ifdef PMGRID
  mesh_factors *mfp = &mf[LOW_MESH];

#ifdef PLACEHIGHRESREGION
  if((DoPM & TREE_ACTIVE_CUTTOFF_HIGHRES_PM))
    {
      int test_node = nop_source->overlap_flag;

      if(test_point == FLAG_INSIDE && test_node == FLAG_BOUNDARYOVERLAP)
        {
          Terminate("this shouldn't happen any more");
          return NODE_OPEN;
        }
      else
        {
          if(test_node == FLAG_INSIDE && test_point == FLAG_INSIDE)
            mfp = &mf[HIGH_MESH];
        }
    }
#endif

  if(DoPM && r2 > mfp->rcut2)
    {
      /* check whether we can ignore any interactions along this branch */
      MyIntPosType dist_x = nop_source->center[0] - intpos_n[0];
      dist_x              = (((MySignedIntPosType)dist_x) >= 0) ? dist_x : -dist_x;
      if(dist_x > mfp->intrcut[0] + halflen)
        return NODE_DISCARD;

      MyIntPosType dist_y = nop_source->center[1] - intpos_n[1];
      dist_y              = (((MySignedIntPosType)dist_y) >= 0) ? dist_y : -dist_y;
      if(dist_y > mfp->intrcut[1] + halflen)
        return NODE_DISCARD;

      MyIntPosType dist_z = nop_source->center[2] - intpos_n[2];
      dist_z              = (((MySignedIntPosType)dist_z) >= 0) ? dist_z : -dist_z;
      if(dist_z > mfp->intrcut[2] + halflen)
        return NODE_DISCARD;
    }
#endif

  MyReal len  = intlen * Tp->FacIntToCoord;
  MyReal len2 = len * len;

  if(All.RelOpeningCriterionInUse == 0) /* check Barnes-Hut opening criterion */
    {
      if(4 * len2 > r2 * errTolTheta2)
        return NODE_OPEN;
    }
  else /* check relative opening criterion */
    {
      if(4 * len2 > r2 * errTolThetaMax2)
        return NODE_OPEN;

      MyReal mmax = (nop_source->mass < mass_n) ? mass_n : nop_source->mass;
      MyReal amin = errTolForceAcc * ((aold < nop_source->MinOldAcc) ? aold : nop_source->MinOldAcc);

#if(MULTIPOLE_ORDER == 1)
      if(mmax > r2 * amin)
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 2)
      if(square(mmax * len) > r2 * square(r2 * amin))
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 3)
      if(mmax * len2 > r2 * r2 * amin)
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 4)
      if(square(mmax * len * len2) > r2 * square(r2 * r2 * amin))
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 5)
      if(mmax * len2 * len2 > r2 * r2 * r2 * amin)
        return NODE_OPEN;
#endif
    }

#if NSOFTCLASSES > 1
  MyReal h_m = All.ForceSoftening[nop_source->getSofteningClass()];

  if(h_m > h_n)
    {
      if(r2 < h_m * h_m)
        if(All.ForceSoftening[nop_source->minsofttype] < All.ForceSoftening[nop_source->maxsofttype])
          {
            return NODE_OPEN;
          }
    }
#endif

  return NODE_USE;
}

/* function to account for interaction of two nodes in the tree */
void fmm::fmm_force_interact(int no_sink, int no_source, char type_sink, char type_source, unsigned char shmrank_sink,
                             unsigned char shmrank_source, int mintopleafnode, int committed)
{
  if(type_sink <= NODE_TYPE_FETCHED_PARTICLE && type_source <= NODE_TYPE_FETCHED_PARTICLE) /* particle-particle interaction */
    {
      /* nothing to be done, or if we do not deal with at least one local particle */
      if(type_sink == NODE_TYPE_FETCHED_PARTICLE && type_source == NODE_TYPE_FETCHED_PARTICLE)
        return;

      if(no_sink != no_source || shmrank_source != shmrank_sink)  // exclude self-interaction
        fmm_particle_particle_interaction(no_sink, no_source, type_sink, type_source, shmrank_sink, shmrank_source);
    }
  else if(!(type_sink > NODE_TYPE_FETCHED_PARTICLE && type_source > NODE_TYPE_FETCHED_PARTICLE)) /* cell-particle interaction */
    {
      /* we have arranged it such that the particle is always on the sink side, the node on the source side */

      gravnode *noptr_source = get_nodep(no_source, shmrank_source);

      /* noting to be done if we do not deal with any local mass */
      if(fmm_depends_on_local_mass(no_source, shmrank_source) == false &&
         (type_sink == NODE_TYPE_FETCHED_PARTICLE ||
          (type_sink < NODE_TYPE_FETCHED_PARTICLE && shmrank_sink != Shmem.Island_ThisTask)))
        return;

      if(noptr_source->not_empty == 0)
        return;

      if(no_source < MaxPart + MaxNodes)                  // we have a top-levelnode
        if(noptr_source->nextnode >= MaxPart + MaxNodes)  // if the next node is not a top-level, we have a leaf node
          mintopleafnode = no_source;

      MyReal r2;
      vector<MyReal> dxyz;

      int openflag = fmm_evaluate_particle_node_opening_criterion(no_sink, type_sink, shmrank_sink, noptr_source, dxyz, r2);

      if(openflag == NODE_USE)
        {
          fmm_particle_node_interaction(no_sink, no_source, type_sink, type_source, shmrank_sink, shmrank_source, noptr_source, dxyz,
                                        r2);
        }
      else if(openflag == NODE_OPEN) /* open cell in a cell-particle interaction */
        {
          if(noptr_source->cannot_be_opened_locally)
            {
              // are we in the same shared memory node?
              if(Shmem.GetNodeIDForSimulCommRank[noptr_source->OriginTask] == Shmem.GetNodeIDForSimulCommRank[D->ThisTask])
                {
                  Terminate("this should not happen any more");
                }
              else
                {
                  tree_add_to_fetch_stack(noptr_source, no_source, shmrank_source);

                  fmm_add_to_work_stack(no_sink, no_source, shmrank_sink, shmrank_source, mintopleafnode);
                }
            }
          else
            {
              int min_buffer_space =
                  std::min<int>(MaxOnWorkStack - (NumOnWorkStack + NewOnWorkStack), MaxOnFetchStack - NumOnFetchStack);

              if(min_buffer_space >= committed + 8 * TREE_NUM_BEFORE_NODESPLIT)
                fmm_open_node(no_sink, noptr_source, type_sink, shmrank_sink, mintopleafnode,
                              committed + 8 * TREE_NUM_BEFORE_NODESPLIT);
              else
                fmm_add_to_work_stack(no_sink, no_source, shmrank_sink, shmrank_source, mintopleafnode);
            }
        }
    }
  else /* cell - cell interaction */
    {
      gravnode *noptr_sink   = get_nodep(no_sink, shmrank_sink);
      gravnode *noptr_source = get_nodep(no_source, shmrank_source);

      /* at least one of the cells needs to depend on local particles */
      if(fmm_depends_on_local_mass(no_sink, shmrank_sink) || fmm_depends_on_local_mass(no_source, shmrank_source))
        {
          /* both cells need to be non-empty */
          if(noptr_sink->not_empty != 0 && noptr_source->not_empty != 0)
            {
              if(noptr_sink == noptr_source) /* self-interaction */
                {
                  if(no_source < MaxPart + MaxNodes)                  // we have a top-levelnode
                    if(noptr_source->nextnode >= MaxPart + MaxNodes)  // if the next node is not a top-level, we have a leaf node
                      mintopleafnode = no_source;

                  if(noptr_sink->cannot_be_opened_locally)
                    {
                      Terminate("should not happen because we have a self-interaction of a supposedly local node");
                    }
                  else
                    {
                      int min_buffer_space =
                          std::min<int>(MaxOnWorkStack - (NumOnWorkStack + NewOnWorkStack), MaxOnFetchStack - NumOnFetchStack);

                      if(min_buffer_space >= committed + 8 * 8 * TREE_NUM_BEFORE_NODESPLIT * TREE_NUM_BEFORE_NODESPLIT)
                        fmm_open_both(noptr_sink, noptr_sink, mintopleafnode,
                                      committed + 8 * 8 * TREE_NUM_BEFORE_NODESPLIT * TREE_NUM_BEFORE_NODESPLIT);
                      else
                        fmm_add_to_work_stack(no_source, no_sink, shmrank_source, shmrank_sink, mintopleafnode);
                    }
                }
              else
                {
                  MyReal r2;
                  vector<MyReal> dxyz;

                  int openflag = fmm_evaluate_node_node_opening_criterion(noptr_sink, noptr_source, dxyz, r2);

                  if(openflag == NODE_USE)
                    {
                      /* evaluate the interaction */
                      fmm_node_node_interaction(no_sink, no_source, type_sink, type_source, shmrank_sink, shmrank_source, noptr_sink,
                                                noptr_source, dxyz, r2);
                    }
                  else if(openflag == NODE_OPEN)
                    {
                      /* open both */

                      if(no_source < MaxPart + MaxNodes)                  // we have a top-levelnode
                        if(noptr_source->nextnode >= MaxPart + MaxNodes)  // if the next node is not a top-level, we have a leaf node
                          mintopleafnode = std::min<int>(mintopleafnode, no_source);

                      if(no_sink < MaxPart + MaxNodes)                  // we have a top-levelnode
                        if(noptr_sink->nextnode >= MaxPart + MaxNodes)  // if the next node is not a top-level, we have a leaf node
                          mintopleafnode = std::min<int>(mintopleafnode, no_sink);

                      if(noptr_source->cannot_be_opened_locally || noptr_sink->cannot_be_opened_locally)
                        {
                          if(noptr_source->cannot_be_opened_locally && noptr_sink->cannot_be_opened_locally)
                            Terminate("this should not happen, because then both nodes would be foreign");

                          if(noptr_source->cannot_be_opened_locally)
                            tree_add_to_fetch_stack(noptr_source, no_source, shmrank_source);

                          if(noptr_sink->cannot_be_opened_locally)
                            tree_add_to_fetch_stack(noptr_sink, no_sink, shmrank_sink);

                          fmm_add_to_work_stack(no_source, no_sink, shmrank_source, shmrank_sink, mintopleafnode);
                        }
                      else
                        {
                          int min_buffer_space =
                              std::min<int>(MaxOnWorkStack - (NumOnWorkStack + NewOnWorkStack), MaxOnFetchStack - NumOnFetchStack);

                          if(min_buffer_space >= committed + 8 * 8 * TREE_NUM_BEFORE_NODESPLIT * TREE_NUM_BEFORE_NODESPLIT)
                            fmm_open_both(noptr_sink, noptr_source, mintopleafnode,
                                          committed + 8 * 8 * TREE_NUM_BEFORE_NODESPLIT * TREE_NUM_BEFORE_NODESPLIT);
                          else
                            fmm_add_to_work_stack(no_source, no_sink, shmrank_source, shmrank_sink, mintopleafnode);
                        }
                    }
                }
            }
        }
    }
}

void fmm::fmm_determine_nodes_with_local_mass(int no, int sib)
{
  gravnode *nop = get_nodep(no, Shmem.Island_ThisTask);

  int p = nop->nextnode;

  /* if the next node is not a top-level node, we have reached a leaf node, and we need to do nothing */
  if(p < MaxPart || p >= FirstNonTopLevelNode)
    return;

  unsigned char depends_on_local_mass = 0;

  while(p != nop->sibling)
    {
      if(p >= 0)
        {
          if(p >= MaxPart && p < MaxPart + MaxNodes) /* we have an internal node */
            {
              int nextsib = get_nodep(p, Shmem.Island_ThisTask)->sibling;

              fmm_determine_nodes_with_local_mass(p, nextsib);
            }

          if(p < MaxPart) /* a particle */
            {
              Terminate("stop");
            }
          else if(p < MaxPart + MaxNodes) /* an internal node  */
            {
              depends_on_local_mass |= Topnode_depends_on_local_mass[p];

              p = get_nodep(p, Shmem.Island_ThisTask)->sibling;
            }
          else if(p < MaxPart + MaxNodes + D->NTopleaves) /* a pseudo particle */
            {
              /* we are processing a local leaf-node which does not have any particles.
               * can continue to the next element, which should end the work.
               */
              Terminate("stop");
            }
          else
            {
              Terminate("stop");
            }
        }
    }

  Topnode_depends_on_local_mass[no] = depends_on_local_mass;
}

void fmm::gravity_fmm(int timebin)
{
  interactioncountPP        = 0;
  interactioncountPN        = 0;
  interactioncountNN        = 0;
  interactioncountEffective = 0;

  TIMER_STORE;
  TIMER_START(CPU_TREE);

  D->mpi_printf("FMM: Begin tree force. timebin=%d  (presently allocated=%g MB)\n", timebin, All.ErrTolTheta,
                Mem.getAllocatedBytesInMB());

#ifdef PMGRID
  set_mesh_factors();
#endif

  Topnode_depends_on_local_mass = (char *)Mem.mymalloc_clear("Topnode_depends_on_local_mass", D->NTopnodes * sizeof(char));
  Topnode_depends_on_local_mass -= MaxPart;

  for(int n = 0; n < D->NTopleaves; n++)
    {
      if(D->TaskOfLeaf[n] == D->ThisTask)
        {
          int no = NodeIndex[n];

          if(TopNodes[no].not_empty)
            Topnode_depends_on_local_mass[no] = 1;
        }
    }

  fmm_determine_nodes_with_local_mass(MaxPart, -1);

  TIMER_START(CPU_TREESTACK);

  NumOnWorkStack         = 0;
  AllocWorkStackBaseLow  = std::max<int>(1.5 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
  AllocWorkStackBaseHigh = AllocWorkStackBaseLow + TREE_EXPECTED_CYCLES * TREE_MIN_WORKSTACK_SIZE;
  MaxOnWorkStack         = std::max<int>(AllocWorkStackBaseLow, 2 * 8 * 8 * TREE_NUM_BEFORE_NODESPLIT * TREE_NUM_BEFORE_NODESPLIT);

  FMM_WorkStack   = (fmm_workstack_data *)Mem.mymalloc("FMM_WorkStack", AllocWorkStackBaseHigh * sizeof(fmm_workstack_data));
  ResultIndexList = (int *)Mem.mymalloc("ResultIndexList", NumPartImported * sizeof(int));

  for(int i = 0; i < Tp->TimeBinsGravity.NActiveParticles; i++)
    {
      int target = Tp->TimeBinsGravity.ActiveParticleList[i];

      /* let's do a safety check here to protect against accidental use of zero softening lengths */
      int softtype = Tp->P[target].getSofteningClass();
      if(All.ForceSoftening[softtype] == 0)
        Terminate("Particle with ID=%lld of type=%d and softening type=%d was assigned zero softening\n",
                  (long long)Tp->P[target].ID.get(), Tp->P[target].getType(), softtype);
    }

  int ncount = 0;

  for(int i = 0; i < NumPartImported; i++)
    {
#ifndef HIERARCHICAL_GRAVITY
      if(Points[i].ActiveFlag)
#endif
        {
          ResultIndexList[i] = ncount++;
        }
    }

  NumOnWorkStack = 0;
  NewOnWorkStack = 0;

  /* for starting, request the self-interaction between the root node */
  fmm_add_to_work_stack(MaxPart, MaxPart, Shmem.Island_ThisTask, Shmem.Island_ThisTask, MaxPart + D->NTopnodes);

  NumOnWorkStack = NewOnWorkStack;

  ResultsActiveImported =
      (resultsactiveimported_data *)Mem.mymalloc_clear("ResultsActiveImported", ncount * sizeof(resultsactiveimported_data));

  TaylorCoeff = (taylor_data *)Mem.mymalloc_clear("TaylorCoeff", NumNodes * sizeof(taylor_data));
  TaylorCoeff -= MaxPart;

  /******************************************/
  /* now execute the tree walk calculations */
  /******************************************/

  errTolForceAcc  = All.ErrTolForceAcc;
  errTolThetaMax2 = All.ErrTolThetaMax * All.ErrTolThetaMax;
  errTolTheta2    = All.ErrTolTheta * All.ErrTolTheta;

  sum_NumForeignNodes  = 0;
  sum_NumForeignPoints = 0;

  // set a default size of the fetch stack equal to half the work stack (this may still be somewhat too large)
  MaxOnFetchStack = std::max<int>(0.1 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
  MaxOnFetchStack = std::max<int>(MaxOnFetchStack, 2 * 8 * 8 * TREE_NUM_BEFORE_NODESPLIT * TREE_NUM_BEFORE_NODESPLIT);
  StackToFetch    = (fetch_data *)Mem.mymalloc_movable(&StackToFetch, "StackToFetch", MaxOnFetchStack * sizeof(fetch_data));

  // let's grab at most half the still available memory for imported points and nodes
  int nspace = (0.33 * Mem.FreeBytes) / (sizeof(gravnode) + 8 * sizeof(foreign_gravpoint_data));

  MaxForeignNodes  = nspace;
  MaxForeignPoints = 8 * nspace;
  NumForeignNodes  = 0;
  NumForeignPoints = 0;

  /* the following two arrays hold imported tree nodes and imported points to augment the local tree */
  Foreign_Nodes  = (gravnode *)Mem.mymalloc_movable(&Foreign_Nodes, "Foreign_Nodes", MaxForeignNodes * sizeof(gravnode));
  Foreign_Points = (foreign_gravpoint_data *)Mem.mymalloc_movable(&Foreign_Points, "Foreign_Points",
                                                                  MaxForeignPoints * sizeof(foreign_gravpoint_data));

  tree_initialize_leaf_node_access_info();

  TIMER_STOP(CPU_TREESTACK);

  double t0       = Logs.second();
  int max_ncycles = 0;

  prepare_shared_memory_access();

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

          NumOnWorkStack = 0;
          NewOnWorkStack = 0;

          /* for starting, request the self-interaction between the root node */
          fmm_add_to_work_stack(MaxPart, MaxPart, Shmem.Island_ThisTask, Shmem.Island_ThisTask, MaxPart + D->NTopnodes);

          NumOnWorkStack = NewOnWorkStack;
        }
#endif

      while(NumOnWorkStack > 0)  // repeat until we are out of work
        {
          NewOnWorkStack  = 0;  // gives the new entries
          NumOnFetchStack = 0;
          MaxOnWorkStack  = std::min<int>(AllocWorkStackBaseLow + max_ncycles * TREE_MIN_WORKSTACK_SIZE, AllocWorkStackBaseHigh);
          MaxOnWorkStack  = std::max<int>(MaxOnWorkStack, 2 * 8 * 8 * TREE_NUM_BEFORE_NODESPLIT * TREE_NUM_BEFORE_NODESPLIT);

          TIMER_START(CPU_TREEWALK);

          int item = 0;

          while(item < NumOnWorkStack)
            {
              int min_buffer_space =
                  std::min<int>(MaxOnWorkStack - (NumOnWorkStack + NewOnWorkStack), MaxOnFetchStack - NumOnFetchStack);

              int committed = 8 * 8 * TREE_NUM_BEFORE_NODESPLIT * TREE_NUM_BEFORE_NODESPLIT;

              if(min_buffer_space >= committed)
                {
                  int no1                = FMM_WorkStack[item].Node1;
                  int no2                = FMM_WorkStack[item].Node2;
                  unsigned char shmrank1 = FMM_WorkStack[item].ShmRank1;
                  unsigned char shmrank2 = FMM_WorkStack[item].ShmRank2;
                  int mintopleaf         = FMM_WorkStack[item].MinTopLeafNode;
                  item++;

                  char type1 = 0, type2 = 0;

                  if(no1 < MaxPart) /* a local particle */
                    type1 = NODE_TYPE_LOCAL_PARTICLE;
                  else if(no1 < MaxPart + MaxNodes) /* an internal node  */
                    type1 = NODE_TYPE_LOCAL_NODE;
                  else if(no1 >= ImportedNodeOffset && no1 < EndOfTreePoints) /* an imported Treepoint particle  */
                    type1 = NODE_TYPE_TREEPOINT_PARTICLE;
                  else if(no1 >= EndOfTreePoints && no1 < EndOfForeignNodes) /* an imported LET node */
                    type1 = NODE_TYPE_FETCHED_NODE;
                  else if(no1 >= EndOfForeignNodes) /* an imported LED particle */
                    type1 = NODE_TYPE_FETCHED_PARTICLE;

                  if(no2 < MaxPart) /* a local particle */
                    type2 = NODE_TYPE_LOCAL_PARTICLE;
                  else if(no2 < MaxPart + MaxNodes) /* an internal node  */
                    type2 = NODE_TYPE_LOCAL_NODE;
                  else if(no2 >= ImportedNodeOffset && no2 < EndOfTreePoints) /* an imported Treepoint particle  */
                    type2 = NODE_TYPE_TREEPOINT_PARTICLE;
                  else if(no2 >= EndOfTreePoints && no2 < EndOfForeignNodes) /* an imported LET node */
                    type2 = NODE_TYPE_FETCHED_NODE;
                  else if(no2 >= EndOfForeignNodes) /* an imported LED particle */
                    type2 = NODE_TYPE_FETCHED_PARTICLE;

                  if(no1 == MaxPart && no2 == MaxPart)
                    {
                      // we have the interaction between the two root nodes
                      fmm_force_interact(no1, no2, type1, type2, shmrank1, shmrank2, mintopleaf, committed);
                    }
                  else
                    {
                      if(type1 > NODE_TYPE_FETCHED_PARTICLE && type2 > NODE_TYPE_FETCHED_PARTICLE)
                        {
                          /* node-node interaction */

                          gravnode *nop1 = get_nodep(no1, shmrank1);
                          gravnode *nop2 = get_nodep(no2, shmrank2);

                          if(nop1->cannot_be_opened_locally || nop2->cannot_be_opened_locally)
                            Terminate("how can this be");

                          fmm_open_both(nop1, nop2, mintopleaf, committed);
                        }
                      else
                        {
                          /* particle-node interaction, particle should on the sink side */

                          // we have a node that we previously could not open
                          gravnode *nop2 = get_nodep(no2, shmrank2);

                          if(nop2->cannot_be_opened_locally)
                            Terminate("now we should be able to open it!");

                          fmm_open_node(no1, nop2, type1, shmrank1, mintopleaf, committed);
                        }
                    }
                }
              else
                break;
            }

          if(item == 0 && NumOnWorkStack > 0)
            Terminate("Can't even process a single particle");

          TIMER_STOP(CPU_TREEWALK);

          TIMER_START(CPU_TREEFETCH);

          tree_fetch_foreign_nodes(FETCH_GRAVTREE);

          TIMER_STOP(CPU_TREEFETCH);

          TIMER_START(CPU_TREESTACK);

          /* now reorder the workstack such that we are first going to do residual pristine particles, and then
           * imported nodes that hang below the first leaf nodes */
          NumOnWorkStack = NumOnWorkStack - item + NewOnWorkStack;
          memmove(FMM_WorkStack, FMM_WorkStack + item, NumOnWorkStack * sizeof(fmm_workstack_data));

          /* now let's sort such that we can go deep on top-level node branches, allowing us to clear them out eventually */
          mycxxsort(FMM_WorkStack, FMM_WorkStack + NumOnWorkStack, compare_fmm_workstack);

          TIMER_STOP(CPU_TREESTACK);

          max_ncycles++;
        }

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
    }
#endif

  TIMER_START(CPU_TREEIMBALANCE);

  MPI_Allreduce(MPI_IN_PLACE, &max_ncycles, 1, MPI_INT, MPI_MAX, D->Communicator);

  TIMER_STOP(CPU_TREEIMBALANCE);

  cleanup_shared_memory_access();

  /* free temporary buffers */

  Mem.myfree(Foreign_Points);
  Mem.myfree(Foreign_Nodes);
  Mem.myfree(StackToFetch);

  TIMER_START(CPU_TREEWALK);

  taylor_data taylor_current{};  // note: the curly braces initialize this to zero in this case

  /* propagate node expansions to particles */
  if(fmm_depends_on_local_mass(MaxPart, Shmem.Island_ThisTask))
    fmm_force_passdown(MaxPart, Shmem.Island_ThisTask, taylor_current);

  TIMER_STOP(CPU_TREEWALK);

  Mem.myfree(TaylorCoeff + MaxPart);

  double t1 = Logs.second();

  D->mpi_printf("FMM: Forces calculated, with %d cycles took %g sec\n", max_ncycles, Logs.timediff(t0, t1));

  /* now communicate the forces in ResultsActiveImported */
  gravity_exchange_forces();

  Mem.myfree(ResultsActiveImported);
  Mem.myfree(ResultIndexList);
  Mem.myfree(FMM_WorkStack);

  TIMER_STOP(CPU_TREE);

  D->mpi_printf("FMM: tree-force is done.\n");

  /*  gather some diagnostic information */

  TIMER_START(CPU_LOGS);

  struct detailed_timings
  {
    double all, tree, wait, fetch, stack, lastpm;
    double costtotal, numnodes;
    double interactioncountEffective;
    double interactioncountPP, interactioncountPN, interactioncountNN;
    double NumForeignNodes, NumForeignPoints;
    double fillfacFgnNodes, fillfacFgnPoints;
    double sumcost;
  };
  detailed_timings timer, tisum, timax;

  memset(&timer, 0, sizeof(detailed_timings));

  if(MeasureCostFlag)
    {
      double sum = 0;
      for(int i = 0; i < Tp->NumPart; i++)
        if(Tp->TimeBinSynchronized[Tp->P[i].TimeBinGrav])
          sum += Tp->P[i].GravCost;

      timer.sumcost = sum;
    }

  timer.tree                      = TIMER_DIFF(CPU_TREEWALK);
  timer.wait                      = TIMER_DIFF(CPU_TREEIMBALANCE);
  timer.fetch                     = TIMER_DIFF(CPU_TREEFETCH);
  timer.stack                     = TIMER_DIFF(CPU_TREESTACK);
  timer.all                       = timer.tree + timer.wait + timer.fetch + timer.stack + TIMER_DIFF(CPU_TREE);
  tisum.lastpm                    = All.CPUForLastPMExecution;
  timer.costtotal                 = interactioncountPP + interactioncountEffective;
  timer.interactioncountPP        = interactioncountPP;
  timer.interactioncountPN        = interactioncountPN;
  timer.interactioncountNN        = interactioncountNN;
  timer.interactioncountEffective = interactioncountEffective;
  timer.NumForeignNodes           = NumForeignNodes;
  timer.NumForeignPoints          = NumForeignPoints;
  timer.fillfacFgnNodes           = NumForeignNodes / ((double)MaxForeignNodes);
  timer.fillfacFgnPoints          = NumForeignPoints / ((double)MaxForeignPoints);
  timer.numnodes                  = NumNodes;

  MPI_Reduce((double *)&timer, (double *)&tisum, (int)(sizeof(detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_SUM, 0,
             D->Communicator);
  MPI_Reduce((double *)&timer, (double *)&timax, (int)(sizeof(detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_MAX, 0,
             D->Communicator);

  All.TotNumOfForces += Tp->TimeBinsGravity.GlobalNActiveParticles;

  if(D->ThisTask == 0)
    {
      fprintf(Logs.FdTimings, "Nf=%9lld FMM  timebin=%d  total-Nf=%lld\n", Tp->TimeBinsGravity.GlobalNActiveParticles, timebin,
              All.TotNumOfForces);
      fprintf(Logs.FdTimings, "   work-load balance: %g   part/sec: raw=%g, effective=%g     ia/part: avg=%g   (%g|%g|%g)\n",
              timax.tree / ((tisum.tree + 1e-20) / D->NTask), Tp->TimeBinsGravity.GlobalNActiveParticles / (tisum.tree + 1.0e-20),
              Tp->TimeBinsGravity.GlobalNActiveParticles / ((timax.tree + 1.0e-20) * D->NTask),
              tisum.costtotal / (Tp->TimeBinsGravity.GlobalNActiveParticles + 1.0e-20),
              tisum.interactioncountPP / (Tp->TimeBinsGravity.GlobalNActiveParticles + 1.0e-20),
              tisum.interactioncountPN / (Tp->TimeBinsGravity.GlobalNActiveParticles + 1.0e-20),
              tisum.interactioncountNN / (Tp->TimeBinsGravity.GlobalNActiveParticles + 1.0e-20));
      fprintf(Logs.FdTimings,
              "   maximum number of nodes: %g, filled: %g  NumForeignNodes: max=%g avg=%g fill=%g NumForeignPoints: max=%g avg=%g "
              "fill=%g  cycles=%d\n",
              timax.numnodes, timax.numnodes / MaxNodes, timax.NumForeignNodes, tisum.NumForeignNodes / D->NTask,
              timax.fillfacFgnNodes, timax.NumForeignPoints, tisum.NumForeignPoints / D->NTask, timax.fillfacFgnPoints, max_ncycles);
      fprintf(Logs.FdTimings,
              "   avg times: <all>=%g  <tree>=%g  <wait>=%g  <fetch>=%g  <stack>=%g  "
              "(lastpm=%g) sec\n",
              tisum.all / D->NTask, tisum.tree / D->NTask, tisum.wait / D->NTask, tisum.fetch / D->NTask, tisum.stack / D->NTask,
              tisum.lastpm / D->NTask);
      fprintf(Logs.FdTimings, "   total interaction cost: %g  (imbalance=%g)  total cost measure: %g %g\n", tisum.costtotal,
              timax.costtotal / (tisum.costtotal / D->NTask), tisum.sumcost,
              tisum.interactioncountPP + tisum.interactioncountEffective);
      myflush(Logs.FdTimings);
    }

  Mem.myfree(Topnode_depends_on_local_mass + MaxPart);

  TIMER_STOP(CPU_LOGS);
}

#endif
