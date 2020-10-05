/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  gwalk.cc
 *
 *  \brief implements the routines for walking the gravity tree and accumulating forces
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../gravity/ewald.h"
#include "../gravtree/gravtree.h"
#include "../gravtree/gwalk.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../sort/cxxsort.h"
#include "../sort/peano.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/*! This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local particles, and particles are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */

inline void gwalk::evaluate_particle_particle_interaction(const pinfo &pdat, const int no, const char jtype, int shmrank)
{
#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  if(skip_actual_force_computation)
    return;
#endif

  /* first, let's get the relevant information of our partner particle, which can be in three different lists */

  MyIntPosType *intpos;
  MyReal mass;
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  int overlap_check;
#endif

  MyReal hmax = pdat.h_i;

  if(jtype == NODE_TYPE_LOCAL_PARTICLE)
    {
      particle_data *P = get_Pp(no, shmrank);

      intpos = P->IntPos;
      mass   = P->getMass();

#if NSOFTCLASSES > 1
      MyReal h_j = All.ForceSoftening[P->getSofteningClass()];
      if(h_j > hmax)
        hmax = h_j;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      overlap_check = P->InsideOutsideFlag;
#endif
    }
  else if(jtype == NODE_TYPE_TREEPOINT_PARTICLE)
    {
      int n = no - ImportedNodeOffset;

      gravpoint_data *Pointp = get_pointsp(n, shmrank);

      intpos = Pointp->IntPos;
      mass   = Pointp->Mass;

#if NSOFTCLASSES > 1
      MyReal h_j = All.ForceSoftening[Pointp->getSofteningClass()];
      if(h_j > hmax)
        hmax = h_j;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      overlap_check = Pointp->InsideOutsideFlag;
#endif
    }
  else /* a point that was fetched */
    {
      int n = no - EndOfForeignNodes;

      foreign_gravpoint_data *foreignpoint = get_foreignpointsp(n, shmrank);

      intpos = foreignpoint->IntPos;
      mass   = foreignpoint->Mass;
#if NSOFTCLASSES > 1
      MyReal h_j = All.ForceSoftening[foreignpoint->getSofteningClass()];
      if(h_j > hmax)
        hmax = h_j;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
      overlap_check = foreignpoint->InsideOutsideFlag;
#endif
    }

#ifdef PMGRID
  mesh_factors *mfp = &mf[LOW_MESH];
#if defined(PLACEHIGHRESREGION)
  if((DoPM & TREE_ACTIVE_CUTTOFF_HIGHRES_PM))
    {
      if(overlap_check == FLAG_INSIDE && pdat.InsideOutsideFlag == FLAG_INSIDE)
        mfp = &mf[HIGH_MESH];
    }
#endif
#endif

  vector<MyReal> dxyz;
  Tp->nearest_image_intpos_to_pos(intpos, pdat.intpos, dxyz.da);

  MyReal r2   = dxyz.r2();
  MyReal r    = sqrt(r2);
  MyReal rinv = (r > 0) ? 1 / r : 0;

  gravtree<simparticles>::gfactors gfac;

#ifdef PMGRID
  if((DoPM & (TREE_ACTIVE_CUTTOFF_BASE_PM + TREE_ACTIVE_CUTTOFF_HIGHRES_PM)))
    {
      if(modify_gfactors_pm_monopole(gfac, r, rinv, mfp))
        return;  // if we are outside the cut-off radius, we have no interaction
    }
#endif

  get_gfactors_monopole(gfac, r, hmax, rinv);

#ifdef EVALPOTENTIAL
  *pdat.pot -= mass * gfac.fac0;
#endif
  *pdat.acc -= (mass * gfac.fac1 * rinv) * dxyz;

  if(DoEwald)
    {
      // EWALD treatment, only done for periodic boundaries in case PM is not active

      ewald_data ew;
      Ewald.ewald_gridlookup(intpos, pdat.intpos, ewald::POINTMASS, ew);

#ifdef EVALPOTENTIAL
      *pdat.pot += mass * ew.D0phi;
#endif
      *pdat.acc += mass * ew.D1phi;
    }

  if(MeasureCostFlag)
    *pdat.GravCost += 1;

  interactioncountPP += 1;
}

inline int gwalk::evaluate_particle_node_opening_criterion_and_interaction(const pinfo &pdat, gravnode *nop)
{
  if(nop->level == 0)  // always open the root node (note: full node length does not fit in the integer type)
    return NODE_OPEN;

  MyIntPosType halflen = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - nop->level);
  MyIntPosType intlen  = halflen << 1;

#if defined(PMGRID) || !defined(TREE_NO_SAFETY_BOX)
  MyIntPosType dist[3];
  Tp->nearest_image_intpos_to_absolute_intdist(nop->center.da, pdat.intpos, dist);
#endif

#ifndef TREE_NO_SAFETY_BOX
  // if we are close to the node, and therefore open it to protect against worst-case force errors
  if(dist[0] < intlen && dist[1] < intlen && dist[2] < intlen)
    return NODE_OPEN;
#endif

#ifdef PMGRID
  mesh_factors *mfp = &mf[LOW_MESH];

  if((DoPM & (TREE_ACTIVE_CUTTOFF_BASE_PM + TREE_ACTIVE_CUTTOFF_HIGHRES_PM)))
    {
#ifdef PLACEHIGHRESREGION
      if((DoPM & TREE_ACTIVE_CUTTOFF_HIGHRES_PM))
        {
          int overlap_check = nop->overlap_flag;

          if(pdat.InsideOutsideFlag == FLAG_INSIDE && overlap_check == FLAG_BOUNDARYOVERLAP)
            Terminate(
                "this case should not happen:  node center=(%g|%g|%g)  len=%g   particle=(%g|%g|%g)  "
                "rel-dist=(%g|%g|%g)\n",
                nop->center[0] * Tp->FacIntToCoord, nop->center[1] * Tp->FacIntToCoord, nop->center[2] * Tp->FacIntToCoord,
                intlen * Tp->FacIntToCoord, pdat.intpos[0] * Tp->FacIntToCoord, pdat.intpos[1] * Tp->FacIntToCoord,
                pdat.intpos[2] * Tp->FacIntToCoord, pdat.intpos[0] * Tp->FacIntToCoord - nop->center[0] * Tp->FacIntToCoord,
                pdat.intpos[1] * Tp->FacIntToCoord - nop->center[1] * Tp->FacIntToCoord,
                pdat.intpos[2] * Tp->FacIntToCoord - nop->center[2] * Tp->FacIntToCoord);

          if(overlap_check == FLAG_INSIDE && pdat.InsideOutsideFlag == FLAG_INSIDE)
            mfp = &mf[HIGH_MESH];
        }
#endif

      /* check whether we can stop walking along this branch */
      if(dist[0] > mfp->intrcut[0] + halflen)
        return NODE_DISCARD;

      /* check whether we can stop walking along this branch */
      if(dist[1] > mfp->intrcut[1] + halflen)
        return NODE_DISCARD;

      /* check whether we can stop walking along this branch */
      if(dist[2] > mfp->intrcut[2] + halflen)
        return NODE_DISCARD;
    }
#endif

  /* converts the integer distance to floating point */
  vector<MyReal> dxyz;
  Tp->nearest_image_intpos_to_pos(nop->s.da, pdat.intpos, dxyz.da);

  MyReal r2 = dxyz.r2();

  MyReal mass = nop->mass;

  MyReal len  = intlen * Tp->FacIntToCoord;
  MyReal len2 = len * len;

  if(All.RelOpeningCriterionInUse == 0) /* check Barnes-Hut opening criterion */
    {
      if(len2 > r2 * theta2)
        return NODE_OPEN;
    }
  else /* check relative opening criterion */
    {
#if(MULTIPOLE_ORDER <= 2)
      if(mass * len2 > r2 * r2 * errTolForceAcc * pdat.aold)
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 3)
      if(square(mass * len * len2) > r2 * square(r2 * r2 * errTolForceAcc * pdat.aold))
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 4)
      if(mass * len2 * len2 > r2 * r2 * r2 * errTolForceAcc * pdat.aold)
        return NODE_OPEN;
#elif(MULTIPOLE_ORDER == 5)
      if(square(mass * len2 * len2 * len) > r2 * square(r2 * r2 * r2 * errTolForceAcc * pdat.aold))
        return NODE_OPEN;
#endif
      // carry out an additional test to protect against pathological force errors for very large opening angles
      if(len2 > r2 * thetamax2)
        return NODE_OPEN;
    }

  MyReal hmax = pdat.h_i;

#if NSOFTCLASSES > 1
  MyReal h_j = All.ForceSoftening[nop->maxsofttype];

  if(h_j > hmax)
    {
      if(r2 < h_j * h_j)
        {
          if(All.ForceSoftening[nop->minsofttype] < All.ForceSoftening[nop->maxsofttype])
            return NODE_OPEN;
        }
      hmax = h_j;
    }
#endif

    /**************************/

    // now evaluate the multipole moment interaction

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  if(skip_actual_force_computation)
    return NODE_USE;
#endif

  MyReal r    = sqrt(r2);
  MyReal rinv = (r > 0) ? 1 / r : 0;

  gravtree<simparticles>::gfactors gfac;

#ifdef PMGRID
  if((DoPM & (TREE_ACTIVE_CUTTOFF_BASE_PM + TREE_ACTIVE_CUTTOFF_HIGHRES_PM)))
    {
      if(modify_gfactors_pm_multipole(gfac, r, rinv, mfp))
        return NODE_DISCARD;  // if we are outside the cut-off radius, we have no interaction
    }
#endif

  get_gfactors_multipole(gfac, r, hmax, rinv);

#ifdef EVALPOTENTIAL
  MyReal g0 = gfac.fac0;
  *pdat.pot -= mass * g0;  //                       monopole potential
#endif

  MyReal g1 = gfac.fac1 * rinv;
  *pdat.acc -= (mass * g1) * dxyz;  //              monopole force

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
  MyReal g2             = gfac.fac2 * gfac.rinv2;
  vector<MyReal> Q2dxyz = nop->Q2Tensor * dxyz;
  MyReal Q2dxyz2        = Q2dxyz * dxyz;
  MyReal Q2trace        = nop->Q2Tensor.trace();
#if(MULTIPOLE_ORDER >= 3)
  MyReal g3 = gfac.fac3 * gfac.rinv3;
  *pdat.acc -= static_cast<MyReal>(0.5) * (g2 * Q2trace + g3 * Q2dxyz2) * dxyz + g2 * Q2dxyz;  //  quadrupole force
#endif
#ifdef EVALPOTENTIAL
  *pdat.pot -= static_cast<MyReal>(0.5) * (g1 * Q2trace + g2 * Q2dxyz2);  //  quadrupole potential
#endif
#endif

#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
  symtensor3<MyDouble> &Q3 = nop->Q3Tensor;

  symtensor2<MyDouble> Q3dxyz = Q3 * dxyz;
  vector<MyDouble> Q3dxyz2    = Q3dxyz * dxyz;
  MyReal Q3dxyz3              = Q3dxyz2 * dxyz;
  MyReal Q3dxyzTrace          = Q3dxyz.trace();

  vector<MyDouble> Q3vec;
  Q3vec[0] = Q3[dXXX] + Q3[dXYY] + Q3[dXZZ];
  Q3vec[1] = Q3[dYXX] + Q3[dYYY] + Q3[dYZZ];
  Q3vec[2] = Q3[dZXX] + Q3[dZYY] + Q3[dZZZ];

#if(MULTIPOLE_ORDER >= 4)
  MyReal g4 = gfac.fac4 * gfac.rinv2 * gfac.rinv2;
  *pdat.acc -=
      static_cast<MyReal>(0.5) *
      (g2 * Q3vec + g3 * Q3dxyz2 + (static_cast<MyReal>(1.0 / 3) * g4 * Q3dxyz3 + g3 * Q3dxyzTrace) * dxyz);  //      octupole force
#endif
#ifdef EVALPOTENTIAL
  *pdat.pot -= static_cast<MyReal>(1.0 / 6) * (3 * g2 * Q3dxyzTrace + g3 * Q3dxyz3);  //       octupole potential
#endif
#endif

#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
  // now compute the hexadecupole force
  symtensor4<MyDouble> &Q4 = nop->Q4Tensor;

  symtensor3<MyDouble> Q4dxyz  = Q4 * dxyz;
  symtensor2<MyDouble> Q4dxyz2 = Q4dxyz * dxyz;
  vector<MyDouble> Q4dxyz3     = Q4dxyz2 * dxyz;
  MyReal Q4dxyz4               = Q4dxyz3 * dxyz;
  MyReal Q4dxyz2trace          = Q4dxyz2.trace();

  symtensor2<MyDouble> QT;
  QT[qXX]        = Q4[sXXXX] + Q4[sXXYY] + Q4[sXXZZ];
  QT[qYY]        = Q4[sYYXX] + Q4[sYYYY] + Q4[sYYZZ];
  QT[qZZ]        = Q4[sZZXX] + Q4[sZZYY] + Q4[sZZZZ];
  QT[qXY]        = Q4[sXYXX] + Q4[sXYYY] + Q4[sXYZZ];
  QT[qXZ]        = Q4[sXZXX] + Q4[sXZYY] + Q4[sXZZZ];
  QT[qYZ]        = Q4[sYZXX] + Q4[sYZYY] + Q4[sYZZZ];
  MyReal QTtrace = QT.trace();

#if(MULTIPOLE_ORDER >= 5)
  vector<MyDouble> QTdxyz = QT * dxyz;
  MyReal g5               = gfac.fac5 * gfac.rinv2 * gfac.rinv3;
  *pdat.acc -=
      static_cast<MyReal>(1.0 / 24) * (g3 * (3 * QTtrace * dxyz + 12 * QTdxyz) + g4 * (6 * Q4dxyz2trace * dxyz + 4 * Q4dxyz3) +
                                       g5 * Q4dxyz4 * dxyz);  //  hexadecupole force
#endif
#ifdef EVALPOTENTIAL
  *pdat.pot -= static_cast<MyReal>(1.0 / 24) * (g2 * 3 * QTtrace + g3 * 6 * Q4dxyz2trace + g4 * Q4dxyz4);  //  hexadecupole potential
#endif
#endif

#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM) && defined(EVALPOTENTIAL))
  symtensor5<MyDouble> &Q5 = nop->Q5Tensor;

  symtensor4<MyDouble> Q5dxyz  = Q5 * dxyz;
  symtensor3<MyDouble> Q5dxyz2 = Q5dxyz * dxyz;
  symtensor2<MyDouble> Q5dxyz3 = Q5dxyz2 * dxyz;
  vector<MyDouble> Q5dxyz4     = Q5dxyz3 * dxyz;
  MyReal Q5dxyz5               = Q5dxyz4 * dxyz;

  MyReal Q5dxyzTtrace = Q5dxyz[sXXXX] + Q5dxyz[sYYYY] + Q5dxyz[sZZZZ] + 2 * (Q5dxyz[sXXYY] + Q5dxyz[sXXZZ] + Q5dxyz[sYYZZ]);

  // now compute the triakontadipole  potential term
  *pdat.pot -= static_cast<MyReal>(1.0 / 120) * (g3 * 15 * Q5dxyzTtrace + g4 * 10 * Q5dxyz3.trace() + g5 * Q5dxyz5);
#endif

  if(DoEwald)
    {
      // EWALD treatment, only done for periodic boundaries in case PM is not active

      ewald_data ew;
      Ewald.ewald_gridlookup(nop->s.da, pdat.intpos, ewald::MULTIPOLES, ew);

#ifdef EVALPOTENTIAL
      *pdat.pot += mass * ew.D0phi;
#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
      *pdat.pot += static_cast<MyReal>(0.5) * (nop->Q2Tensor * ew.D2phi);
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
      *pdat.pot += static_cast<MyReal>(1.0 / 6) * (nop->Q3Tensor * ew.D3phi);
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
      *pdat.pot += static_cast<MyReal>(1.0 / 24) * (nop->Q4Tensor * ew.D4phi);
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM) && defined(EVALPOTENTIAL))
      *pdat.pot += static_cast<MyReal>(1.0 / 120) * (nop->Q5Tensor * ew.D5phi);
#endif
#endif
      *pdat.acc += mass * ew.D1phi;
#if(MULTIPOLE_ORDER >= 3)
      *pdat.acc += static_cast<MyReal>(0.5) * (ew.D3phi * nop->Q2Tensor);
#endif
#if(MULTIPOLE_ORDER >= 4)
      *pdat.acc += static_cast<MyReal>(1.0 / 6) * (ew.D4phi * nop->Q3Tensor);
#endif
#if(MULTIPOLE_ORDER >= 5)
      *pdat.acc += static_cast<MyReal>(1.0 / 24) * (ew.D5phi * nop->Q4Tensor);
#endif
    }

  interactioncountPN += 1;

  if(MeasureCostFlag)
    *pdat.GravCost += 1;

  return NODE_USE;
}

inline void gwalk::gwalk_open_node(const pinfo &pdat, int i, char ptype, gravnode *nop, int mintopleafnode, int committed)
{
  /* open node */
  int p                 = nop->nextnode;
  unsigned char shmrank = nop->nextnode_shmrank;

  while(p != nop->sibling || (shmrank != nop->sibling_shmrank && nop->sibling >= MaxPart + D->NTopnodes))
    {
      if(p < 0)
        Terminate(
            "p=%d < 0  nop->sibling=%d nop->nextnode=%d shmrank=%d nop->sibling_shmrank=%d nop->foreigntask=%d  mass=%g  "
            "first_nontoplevelnode=%d",
            p, nop->sibling, nop->nextnode, shmrank, nop->sibling_shmrank, nop->OriginTask, nop->mass, MaxPart + D->NTopnodes);

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

      gravity_force_interact(pdat, i, p, ptype, type, shmrank, mintopleafnode, committed);

      p       = next;
      shmrank = next_shmrank;
    }
}

void gwalk::gravity_force_interact(const pinfo &pdat, int i, int no, char ptype, char no_type, unsigned char shmrank,
                                   int mintopleafnode, int committed)
{
  if(no_type <= NODE_TYPE_FETCHED_PARTICLE)  // we are interacting with a particle
    {
      evaluate_particle_particle_interaction(pdat, no, no_type, shmrank);
    }
  else  // we are interacting with a node
    {
      gravnode *nop = get_nodep(no, shmrank);

      if(nop->not_empty == 0)
        return;

      if(no < MaxPart + MaxNodes)                // we have a top-level node
        if(nop->nextnode >= MaxPart + MaxNodes)  // if the next node is not a top-level, we have a leaf node
          {
            mintopleafnode = no;

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
            // if the leaf node is on this shared memory, we have all the data, so we for sure don't need to import anything on this
            // branch
            if(skip_actual_force_computation)
              if(Shmem.GetNodeIDForSimulCommRank[nop->OriginTask] == Shmem.GetNodeIDForSimulCommRank[D->ThisTask])
                return;
#endif
          }

      int openflag = evaluate_particle_node_opening_criterion_and_interaction(pdat, nop);

      if(openflag == NODE_OPEN) /* cell can't be used, need to open it */
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

                  tree_add_to_work_stack(i, no, shmrank, mintopleafnode);
                }
            }
          else
            {
              int min_buffer_space =
                  std::min<int>(MaxOnWorkStack - (NumOnWorkStack + NewOnWorkStack), MaxOnFetchStack - NumOnFetchStack);

              if(min_buffer_space >= committed + 8 * TREE_NUM_BEFORE_NODESPLIT)
                gwalk_open_node(pdat, i, ptype, nop, mintopleafnode, committed + 8 * TREE_NUM_BEFORE_NODESPLIT);
              else
                tree_add_to_work_stack(i, no, shmrank, mintopleafnode);
            }
        }
    }
}

/*

*/

/*! \brief This function computes the gravitational forces for all active particles.
 *
 * The tree walk is done in two phases: First the local part of the force tree is processed (gravity_primary_loop() ).
 * Whenever an external node is encountered during the walk, this node is saved on a list.
 * This node list along with data about the particles is then exchanged among tasks.
 * In the second phase (gravity_secondary_loop() ) each task now continues the tree walk for
 * the imported particles. Finally the resulting partial forces are send back to the original task
 * and are summed up there to complete the tree force calculation.
 *
 * Particles are only exported to other processors when really needed, thereby allowing a
 * good use of the communication buffer. Every particle is sent at most once to a given processor
 * together with the complete list of relevant tree nodes to be checked on the other task.
 *
 * Particles which drifted into the domain of another task are sent to this task for the force computation.
 * Afterwards the resulting force is sent back to the originating task.
 *
 * In order to improve the work load balancing during a domain decomposition, the work done by each
 * node/particle is measured. The work is measured for the interaction partners (i.e. the nodes or particles)
 * and not for the particles itself that require a force computation. This way, work done for imported
 * particles is accounted for at the task where the work actually incurred. The cost measurement is
 * only done for the "GRAVCOSTLEVELS" highest occupied time bins. The variable #MeasureCostFlag will state whether a
 * measurement is done at the present time step.
 *
 * The tree imbalance can be further reduced using chunking. The particles requiring a force computation
 * are split into chunks of size #Nchunksize. A set of every #Nchunk -th chunk is processed first.
 * Then the process is repeated, processing the next set of chunks. This way the amount of exported particles
 * is more balanced, as communication heavy regions are mixed with less communication intensive regions.
 *
 */

void gwalk::gravity_tree(int timebin)
{
  interactioncountPP = 0;
  interactioncountPN = 0;

  TIMER_STORE;
  TIMER_START(CPU_TREE);

  D->mpi_printf("GRAVTREE: Begin tree force. timebin=%d (presently allocated=%g MB)\n", timebin, Mem.getAllocatedBytesInMB());

#ifdef PMGRID
  set_mesh_factors();
#endif

  TIMER_START(CPU_TREESTACK);

  // Create list of targets (the work queue). There are initially two possible sources of points, local ones, and imported ones.

  NumOnWorkStack         = 0;
  AllocWorkStackBaseLow  = std::max<int>(1.5 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
  AllocWorkStackBaseHigh = AllocWorkStackBaseLow + TREE_EXPECTED_CYCLES * TREE_MIN_WORKSTACK_SIZE;
  MaxOnWorkStack         = AllocWorkStackBaseLow;

  WorkStack       = (workstack_data *)Mem.mymalloc("WorkStack", AllocWorkStackBaseHigh * sizeof(workstack_data));
  ResultIndexList = (int *)Mem.mymalloc("ResultIndexList", NumPartImported * sizeof(int));

  for(int i = 0; i < Tp->TimeBinsGravity.NActiveParticles; i++)
    {
      int target = Tp->TimeBinsGravity.ActiveParticleList[i];

      // if we have exported particles, we need to explicitly check whether this particle is among them
      if(NumPartExported > 0)
        {
          MyIntPosType xxb       = Tp->P[target].IntPos[0];
          MyIntPosType yyb       = Tp->P[target].IntPos[1];
          MyIntPosType zzb       = Tp->P[target].IntPos[2];
          MyIntPosType mask      = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - 1));
          unsigned char shiftx   = (BITS_FOR_POSITIONS - 3);
          unsigned char shifty   = (BITS_FOR_POSITIONS - 2);
          unsigned char shiftz   = (BITS_FOR_POSITIONS - 1);
          unsigned char level    = 0;
          unsigned char rotation = 0;

          int no = 0;
          while(D->TopNodes[no].Daughter >= 0) /* walk down top tree to find correct leaf */
            {
              unsigned char pix     = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                                   ((unsigned char)((zzb & mask) >> (shiftz--))));
              unsigned char subnode = peano_incremental_key(pix, &rotation);
              mask >>= 1;
              level++;
              no = D->TopNodes[no].Daughter + subnode;
            }

          no       = D->TopNodes[no].Leaf;
          int task = D->TaskOfLeaf[no];

          if(task == D->ThisTask)
            {
              WorkStack[NumOnWorkStack].Target         = target;
              WorkStack[NumOnWorkStack].Node           = MaxPart;
              WorkStack[NumOnWorkStack].ShmRank        = Shmem.Island_ThisTask;
              WorkStack[NumOnWorkStack].MinTopLeafNode = MaxPart + D->NTopnodes;
              NumOnWorkStack++;
            }
        }
      else
        {
          WorkStack[NumOnWorkStack].Target         = target;
          WorkStack[NumOnWorkStack].Node           = MaxPart;
          WorkStack[NumOnWorkStack].ShmRank        = Shmem.Island_ThisTask;
          WorkStack[NumOnWorkStack].MinTopLeafNode = MaxPart + D->NTopnodes;
          NumOnWorkStack++;
        }

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

          WorkStack[NumOnWorkStack].Target         = i + ImportedNodeOffset;
          WorkStack[NumOnWorkStack].Node           = MaxPart;
          WorkStack[NumOnWorkStack].ShmRank        = Shmem.Island_ThisTask;
          WorkStack[NumOnWorkStack].MinTopLeafNode = MaxPart + D->NTopnodes;
          NumOnWorkStack++;
        }
    }

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  workstack_data *WorkStackBak = (workstack_data *)Mem.mymalloc("WorkStackBak", NumOnWorkStack * sizeof(workstack_data));
  int NumOnWorkStackBak        = NumOnWorkStack;
  memcpy(WorkStackBak, WorkStack, NumOnWorkStack * sizeof(workstack_data));
#endif

  ResultsActiveImported =
      (resultsactiveimported_data *)Mem.mymalloc_clear("ResultsActiveImported", ncount * sizeof(resultsactiveimported_data));

  /******************************************/
  /* now execute the tree walk calculations */
  /******************************************/

  theta2         = All.ErrTolTheta * All.ErrTolTheta;
  thetamax2      = All.ErrTolThetaMax * All.ErrTolThetaMax;
  errTolForceAcc = All.ErrTolForceAcc;

  sum_NumForeignNodes  = 0;
  sum_NumForeignPoints = 0;

  // set a default size of the fetch stack equal to half the work stack (this may still be somewhat too large)
  MaxOnFetchStack = std::max<int>(0.1 * (Tp->NumPart + NumPartImported), TREE_MIN_WORKSTACK_SIZE);
  StackToFetch    = (fetch_data *)Mem.mymalloc_movable(&StackToFetch, "StackToFetch", MaxOnFetchStack * sizeof(fetch_data));

  // let's grab at most half the still available memory for imported points and nodes
  int nspace = (0.5 * Mem.FreeBytes) / (sizeof(gravnode) + 8 * sizeof(foreign_gravpoint_data));

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
          NumOnWorkStack                = NumOnWorkStackBak;
          memcpy(WorkStack, WorkStackBak, NumOnWorkStack * sizeof(workstack_data));
        }
#endif

      while(NumOnWorkStack > 0)  // repeat until we are out of work
        {
          NewOnWorkStack  = 0;  // gives the new entries
          NumOnFetchStack = 0;
          MaxOnWorkStack  = std::min<int>(AllocWorkStackBaseLow + max_ncycles * TREE_MIN_WORKSTACK_SIZE, AllocWorkStackBaseHigh);

          TIMER_START(CPU_TREEWALK);

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
                  int ptype = get_pinfo(target, pdat);

                  if(no == MaxPart)
                    {
                      // we have a pristine particle that's processed for the first time
                      gravity_force_interact(pdat, target, no, ptype, NODE_TYPE_LOCAL_NODE, shmrank, mintopleaf, committed);
                    }
                  else
                    {
                      // we have a node that we previously could not open
                      gravnode *nop = get_nodep(no, shmrank);

                      if(nop->cannot_be_opened_locally)
                        {
                          Terminate("item=%d:  no=%d  now we should be able to open it!", item, no);
                        }
                      else
                        gwalk_open_node(pdat, target, ptype, nop, mintopleaf, committed);
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
          memmove(WorkStack, WorkStack + item, NumOnWorkStack * sizeof(workstack_data));

          /* now let's sort such that we can go deep on top-level node branches, allowing us to clear them out eventually */
          mycxxsort(WorkStack, WorkStack + NumOnWorkStack, compare_workstack);

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

  double t1 = Logs.second();

  D->mpi_printf("GRAVTREE: tree-forces are calculated, with %d cycles took %g sec\n", max_ncycles, Logs.timediff(t0, t1));

  /* now communicate the forces in ResultsActiveImported */
  gravity_exchange_forces();

  Mem.myfree(ResultsActiveImported);
#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  Mem.myfree(WorkStackBak);
#endif
  Mem.myfree(ResultIndexList);
  Mem.myfree(WorkStack);

  TIMER_STOP(CPU_TREE);

  D->mpi_printf("GRAVTREE: tree-force is done.\n");

  /*  gather some diagnostic information */

  TIMER_START(CPU_LOGS);

  struct detailed_timings
  {
    double tree, wait, fetch, stack, all, lastpm;
    double costtotal, numnodes;
    double interactioncountPP, interactioncountPN;
    double NumForeignNodes, NumForeignPoints;
    double fillfacFgnNodes, fillfacFgnPoints;
  };
  detailed_timings timer, tisum, timax;

  timer.tree               = TIMER_DIFF(CPU_TREEWALK);
  timer.wait               = TIMER_DIFF(CPU_TREEIMBALANCE);
  timer.fetch              = TIMER_DIFF(CPU_TREEFETCH);
  timer.stack              = TIMER_DIFF(CPU_TREESTACK);
  timer.all                = timer.tree + timer.wait + timer.fetch + timer.stack + TIMER_DIFF(CPU_TREE);
  timer.lastpm             = All.CPUForLastPMExecution;
  timer.costtotal          = interactioncountPP + interactioncountPN;
  timer.numnodes           = NumNodes;
  timer.interactioncountPP = interactioncountPP;
  timer.interactioncountPN = interactioncountPN;
  timer.NumForeignNodes    = NumForeignNodes;
  timer.NumForeignPoints   = NumForeignPoints;
  timer.fillfacFgnNodes    = NumForeignNodes / ((double)MaxForeignNodes);
  timer.fillfacFgnPoints   = NumForeignPoints / ((double)MaxForeignPoints);

  MPI_Reduce((double *)&timer, (double *)&tisum, (int)(sizeof(detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_SUM, 0,
             D->Communicator);
  MPI_Reduce((double *)&timer, (double *)&timax, (int)(sizeof(detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_MAX, 0,
             D->Communicator);

  All.TotNumOfForces += Tp->TimeBinsGravity.GlobalNActiveParticles;

  if(D->ThisTask == 0)
    {
      fprintf(Logs.FdTimings, "Nf=%9lld  timebin=%d  total-Nf=%lld\n", Tp->TimeBinsGravity.GlobalNActiveParticles, timebin,
              All.TotNumOfForces);
      fprintf(Logs.FdTimings, "   work-load balance: %g   part/sec: raw=%g, effective=%g     ia/part: avg=%g   (%g|%g)\n",
              timax.tree / ((tisum.tree + 1e-20) / D->NTask), Tp->TimeBinsGravity.GlobalNActiveParticles / (tisum.tree + 1.0e-20),
              Tp->TimeBinsGravity.GlobalNActiveParticles / ((timax.tree + 1.0e-20) * D->NTask),
              tisum.costtotal / (Tp->TimeBinsGravity.GlobalNActiveParticles + 1.0e-20),
              tisum.interactioncountPP / (Tp->TimeBinsGravity.GlobalNActiveParticles + 1.0e-20),
              tisum.interactioncountPN / (Tp->TimeBinsGravity.GlobalNActiveParticles + 1.0e-20));
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
      fprintf(Logs.FdTimings, "   total interaction cost: %g  (imbalance=%g)\n", tisum.costtotal,
              timax.costtotal / (tisum.costtotal / D->NTask));
      myflush(Logs.FdTimings);
    }

  TIMER_STOP(CPU_LOGS);
}

/* make sure that we instantiate the template */
#include "../data/simparticles.h"
template class gravtree<simparticles>;
