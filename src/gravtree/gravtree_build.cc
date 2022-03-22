/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file gravtree_build.cc
 *
 *  \brief routines for building the gravitational tree
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
#include "../domain/domain.h"
#include "../gravtree/gravtree.h"
#include "../io/io.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../pm/pm.h"
#include "../sort/cxxsort.h"
#include "../sort/peano.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/*!
 *  This file contains the construction of the tree used for calculating the gravitational force.
 *  The type tree implemented is a geometrical oct-tree, starting from a cube encompassing
 *  all particles. This cube is automatically found in the domain decomposition, which also
 *  splits up the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. In this version of the code, the tree
 *  construction may be repeated every timestep without a renewed domain decomposition.
 *  If particles are on the "wrong" processor because a new domain decomposition has not been
 *  carried out, they are sent as temporary points to the right insertion processor according
 *  to the layout of the top-level nodes. In addition, the mapping of the top-level nodes to
 *  processors may be readjusted in order to improve work-load balance for the current time step.
 *
 */

template <typename partset>
void gravtree<partset>::report_log_message(void)
{
  TIMER_START(CPU_LOGS);

  int max_imported;
  long long tot_imported;
  sumup_large_ints(1, &NumPartImported, &tot_imported, D->Communicator);
  MPI_Reduce(&NumPartImported, &max_imported, 1, MPI_INT, MPI_MAX, 0, D->Communicator);

  MyReal numnodes = NumNodes, tot_numnodes;
  MPI_Reduce(&numnodes, &tot_numnodes, 1, MPI_DOUBLE, MPI_SUM, 0, D->Communicator);

  if(Ninsert == Tp->NumPart)
    D->mpi_printf(
        "GRAVTREE: Tree construction done. took %g sec  <numnodes>=%g  NTopnodes=%d NTopleaves=%d tree-build-scalability=%g\n",
        Buildtime, tot_numnodes / D->NTask, D->NTopnodes, D->NTopleaves,
        ((double)((tot_numnodes - D->NTask * ((double)D->NTopnodes)) + D->NTopnodes)) / (tot_numnodes > 0 ? tot_numnodes : 1));

  TIMER_STOP(CPU_LOGS);
}

template <>
void gravtree<simparticles>::fill_in_export_points(gravpoint_data *exp_point, int i, int no)
{
  /* this point has to go to another task */
  for(int j = 0; j < 3; j++)
    exp_point->IntPos[j] = Tp->P[i].IntPos[j];

  exp_point->Mass   = Tp->P[i].getMass();
  exp_point->OldAcc = Tp->P[i].OldAcc;
  exp_point->index  = i;
  exp_point->Type   = Tp->P[i].getType();
#if NSOFTCLASSES > 1
  exp_point->SofteningClass = Tp->P[i].getSofteningClass();
#endif
  exp_point->no = no;
#ifndef HIERARCHICAL_GRAVITY
  if(Tp->TimeBinSynchronized[Tp->P[i].TimeBinGrav])
    exp_point->ActiveFlag = 1;
  else
    exp_point->ActiveFlag = 0;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  exp_point->InsideOutsideFlag = Tp->P[i].InsideOutsideFlag;
#endif
}

#if defined(LIGHTCONE) && (defined(LIGHTCONE_PARTICLES_GROUPS) || defined(LIGHTCONE_IMAGE_COMP_HSML_VELDISP))

template <>
void gravtree<lcparticles>::fill_in_export_points(gravpoint_data *exp_point, int i, int no)
{
  /* this point has to go to another task */
  for(int j = 0; j < 3; j++)
    exp_point->IntPos[j] = Tp->P[i].IntPos[j];

  exp_point->Mass   = Tp->P[i].getMass();
  exp_point->OldAcc = 0;
  exp_point->index  = i;
  exp_point->Type   = Tp->P[i].getType();
#if NSOFTCLASSES > 1
  exp_point->SofteningClass = Tp->P[i].getSofteningClass();
#endif
  exp_point->no = no;
}

#endif

/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
template <typename partset>
void gravtree<partset>::exchange_topleafdata(void)
{
  struct leafnode_data
  {
    vector<MyIntPosType> s;
    MyDouble mass;
#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
    symtensor2<MyDouble> Q2Tensor;
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
    symtensor3<MyDouble> Q3Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
    symtensor4<MyDouble> Q4Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
    symtensor5<MyDouble> Q5Tensor;
#endif
#ifdef FMM
    float MinOldAcc;
#endif
    unsigned char not_empty;
#if NSOFTCLASSES > 1
    unsigned char maxsofttype;
    unsigned char minsofttype;
#endif
    unsigned char level;
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
    unsigned char overlap_flag : 2;
#endif
  };
  leafnode_data *glob_leaf_node_data, *loc_leaf_node_data;

  glob_leaf_node_data = (leafnode_data *)Mem.mymalloc("glob_leaf_node_data", D->NTopleaves * sizeof(leafnode_data));

  /* share the pseudo-particle data accross CPUs */
  int *recvcounts = (int *)Mem.mymalloc("recvcounts", sizeof(int) * D->NTask);
  int *recvoffset = (int *)Mem.mymalloc("recvoffset", sizeof(int) * D->NTask);
  int *bytecounts = (int *)Mem.mymalloc("bytecounts", sizeof(int) * D->NTask);
  int *byteoffset = (int *)Mem.mymalloc("byteoffset", sizeof(int) * D->NTask);

  for(int task = 0; task < D->NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < D->NTopleaves; n++)
    recvcounts[D->TaskOfLeaf[n]]++;

  for(int task = 0; task < D->NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(leafnode_data);

  recvoffset[0] = 0;
  byteoffset[0] = 0;

  for(int task = 1; task < D->NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  loc_leaf_node_data = (leafnode_data *)Mem.mymalloc("loc_leaf_node_data", recvcounts[D->ThisTask] * sizeof(leafnode_data));

  int idx = 0;

  for(int n = 0; n < D->NTopleaves; n++)
    {
      if(D->TaskOfLeaf[n] == D->ThisTask)
        {
          int no        = NodeIndex[n];
          gravnode *nop = &TopNodes[no];

          leafnode_data *locp = &loc_leaf_node_data[idx];

          /* read out the multipole moments from the local base cells */
          locp->s    = nop->s;
          locp->mass = nop->mass;
#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
          locp->Q2Tensor = nop->Q2Tensor;
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
          locp->Q3Tensor = nop->Q3Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
          locp->Q4Tensor = nop->Q4Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
          locp->Q5Tensor = nop->Q5Tensor;
#endif
#if NSOFTCLASSES > 1
          locp->maxsofttype = nop->maxsofttype;
          locp->minsofttype = nop->minsofttype;
#endif
#ifdef FMM
          locp->MinOldAcc = nop->MinOldAcc;
#endif
          locp->not_empty = nop->not_empty;
          locp->level     = nop->level;
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
          locp->overlap_flag = nop->overlap_flag;
#endif
          idx++;
        }
    }

  // optimise this step - only need to update this once per shared memory node

  myMPI_Allgatherv(loc_leaf_node_data, bytecounts[D->ThisTask], MPI_BYTE, glob_leaf_node_data, bytecounts, byteoffset, MPI_BYTE,
                   D->Communicator);

  for(int task = 0; task < D->NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < D->NTopleaves; n++)
    {
      int task = D->TaskOfLeaf[n];
      if(task != D->ThisTask)
        {
          int no        = NodeIndex[n];
          gravnode *nop = &TopNodes[no];

          idx                  = recvoffset[task] + recvcounts[task]++;
          leafnode_data *globp = &glob_leaf_node_data[idx];

          nop->s    = globp->s;
          nop->mass = globp->mass;
#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
          nop->Q2Tensor = globp->Q2Tensor;
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
          nop->Q3Tensor = globp->Q3Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
          nop->Q4Tensor = globp->Q4Tensor;
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
          nop->Q5Tensor = globp->Q5Tensor;
#endif
#if NSOFTCLASSES > 1
          nop->maxsofttype = globp->maxsofttype;
          nop->minsofttype = globp->minsofttype;
#endif
          nop->not_empty = globp->not_empty;
#ifdef FMM
          nop->MinOldAcc = globp->MinOldAcc;
#endif
          nop->level = globp->level;
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
          nop->overlap_flag = globp->overlap_flag;
#endif
        }
    }

  Mem.myfree(loc_leaf_node_data);
  Mem.myfree(byteoffset);
  Mem.myfree(bytecounts);
  Mem.myfree(recvoffset);
  Mem.myfree(recvcounts);
  Mem.myfree(glob_leaf_node_data);
}

/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *  'no' is the node for which the moments shall be found
 *  'sib' is the sibling of this node
 */
template <typename partset>
void gravtree<partset>::update_node_recursive(int no, int sib, int mode)
{
  if(!(no >= MaxPart && no < MaxPart + MaxNodes)) /* are we an internal node? */
    Terminate("no internal node\n");

  gravnode *nop = get_nodep(no);

  if(mode == TREE_MODE_TOPLEVEL)
    {
      int p = nop->nextnode;

      /* if the next node is not a top-level node, we have reached a leaf node, and we need to do nothing */
      if(p < MaxPart || p >= FirstNonTopLevelNode)
        return;
    }

  MyReal mass = 0;
  vector<MyReal> s(0.0);

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
  symtensor2<MyReal> Q2Tensor(0.0); /**< quadrupole tensor */
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
  symtensor3<MyReal> Q3Tensor(0.0); /**< octupole tensor */
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
  symtensor4<MyReal> Q4Tensor(0.0); /**< hexadecupole tensor */
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
  symtensor5<MyReal> Q5Tensor(0.0); /**< triakontadipole tensor */
#endif
#if NSOFTCLASSES > 1
  unsigned char maxsofttype = NSOFTCLASSES + NSOFTCLASSES_HYDRO;
  unsigned char minsofttype = NSOFTCLASSES + NSOFTCLASSES_HYDRO + 1;
#endif
  unsigned char not_empty = 0;
#ifdef FMM
  float minOldAcc = MAX_FLOAT_NUMBER;
#endif

  int p = nop->nextnode;

  while(p != nop->sibling)
    {
      if(p >= 0)
        {
          if(p >= MaxPart && p < MaxPart + MaxNodes) /* we have an internal node */
            {
              int nextsib = get_nodep(p)->sibling;

              update_node_recursive(p, nextsib, mode);
            }

          if(p < MaxPart) /* a particle */
            {
              vector<MyReal> dxyz;
              Tp->nearest_image_intpos_to_pos(Tp->P[p].IntPos, nop->center.da, dxyz.da);

              MyReal m = Tp->P[p].getMass();

              mass += m;
              s += m * dxyz;

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
              vector<MyReal> mdxyz = m * dxyz;
              symtensor2<MyReal> mdxyz2(mdxyz, dxyz);
              Q2Tensor += mdxyz2;
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
              symtensor3<MyReal> mdxyz3(dxyz, mdxyz2);
              Q3Tensor += mdxyz3;
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
              symtensor4<MyReal> mdxyz4(dxyz, mdxyz3);
              Q4Tensor += mdxyz4;
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
              symtensor5<MyReal> mdxyz5(dxyz, mdxyz4);
              Q5Tensor += mdxyz5;
#endif
              not_empty = 1;
#ifndef HIERARCHICAL_GRAVITY
              if(Tp->getTimeBinSynchronized(Tp->P[p].getTimeBinGrav()))
#endif
                {
#ifdef FMM
                  if(minOldAcc > Tp->P[p].getOldAcc())
                    minOldAcc = Tp->P[p].getOldAcc();
#endif
                }

#if NSOFTCLASSES > 1
              if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[Tp->P[p].getSofteningClass()])
                maxsofttype = Tp->P[p].getSofteningClass();
              if(All.ForceSoftening[minsofttype] > All.ForceSoftening[Tp->P[p].getSofteningClass()])
                minsofttype = Tp->P[p].getSofteningClass();
#endif
              p = Nextnode[p];
            }
          else if(p < MaxPart + MaxNodes) /* an internal node  */
            {
              gravnode *noptr = get_nodep(p);

              vector<MyReal> dxyz;
              Tp->nearest_image_intpos_to_pos(noptr->s.da, nop->center.da, dxyz.da);

              MyReal m = noptr->mass;

              mass += m;
              s += m * dxyz;

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
              Q2Tensor += noptr->Q2Tensor;

              vector<MyReal> mdxyz = m * dxyz;
              symtensor2<MyReal> mdxyz2(mdxyz, dxyz);
              Q2Tensor += mdxyz2;
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
              Q3Tensor += noptr->Q3Tensor;

              symtensor3<MyReal> mdxyz3(dxyz, mdxyz2);
              Q3Tensor += mdxyz3;

              Q3Tensor += outer_prod_sum(noptr->Q2Tensor, dxyz);
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
              Q4Tensor += noptr->Q4Tensor;

              symtensor4<MyReal> mdxyz4(dxyz, mdxyz3);
              Q4Tensor += mdxyz4;

              Q4Tensor += outer_prod_sum(noptr->Q3Tensor, dxyz);
              symtensor2<MyReal> dxyz2(dxyz, dxyz);
              Q4Tensor += outer_prod_sum(noptr->Q2Tensor, dxyz2);
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
              Q5Tensor += noptr->Q5Tensor;

              symtensor5<MyReal> mdxyz5(dxyz, mdxyz4);
              Q5Tensor += mdxyz5;

              Q5Tensor += outer_prod_sum(noptr->Q4Tensor, dxyz);
              Q5Tensor += outer_prod_sum(noptr->Q3Tensor, dxyz2);
              symtensor3<MyReal> dxyz3(dxyz, dxyz2);
              Q5Tensor += outer_prod_sum(dxyz3, noptr->Q2Tensor);
#endif

#if NSOFTCLASSES > 1
              if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[noptr->maxsofttype])
                maxsofttype = noptr->maxsofttype;
              if(All.ForceSoftening[minsofttype] > All.ForceSoftening[noptr->minsofttype])
                minsofttype = noptr->minsofttype;
#endif
              not_empty |= noptr->not_empty;
#ifdef FMM
              if(minOldAcc > noptr->MinOldAcc)
                minOldAcc = noptr->MinOldAcc;
#endif

              p = noptr->sibling;
            }
          else if(p < MaxPart + MaxNodes + D->NTopleaves) /* a pseudo particle */
            {
              /* we are processing a local leaf-node which does not have any particles.
               * can continue to the next element, which should end the work.
               */
              p = Nextnode[p - MaxNodes];
            }
          else
            {
              /* an imported point */
              int n = p - ImportedNodeOffset;

              if(n >= NumPartImported)
                Terminate("n=%d >= NumPartImported=%d   MaxPart=%d MaxNodes=%d  D->NTopleaves=%d", n, NumPartImported, MaxPart,
                          MaxNodes, D->NTopleaves);

              vector<MyReal> dxyz;
              Tp->nearest_image_intpos_to_pos(Points[n].IntPos, nop->center.da, dxyz.da);

              MyReal m = Points[n].Mass;

              mass += m;
              s += m * dxyz;
#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
              vector<MyReal> mdxyz = m * dxyz;
              symtensor2<MyReal> mdxyz2(mdxyz, dxyz);
              Q2Tensor += mdxyz2;
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
              symtensor3<MyReal> mdxyz3(dxyz, mdxyz2);
              Q3Tensor += mdxyz3;
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
              symtensor4<MyReal> mdxyz4(dxyz, mdxyz3);
              Q4Tensor += mdxyz4;
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
              symtensor5<MyReal> mdxyz5(dxyz, mdxyz4);
              Q5Tensor += mdxyz5;
#endif
              not_empty = 1;
#ifndef HIERARCHICAL_GRAVITY
              if(Points[n].ActiveFlag)
#endif
                {
#ifdef FMM
                  if(minOldAcc > Points[n].OldAcc)
                    minOldAcc = Points[n].OldAcc;
#endif
                }
#if NSOFTCLASSES > 1
              if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[Points[n].SofteningClass])
                maxsofttype = Points[n].SofteningClass;
              if(All.ForceSoftening[minsofttype] > All.ForceSoftening[Points[n].SofteningClass])
                minsofttype = Points[n].SofteningClass;
#endif
              p = Nextnode[p - MaxNodes];
            }
        }
    }

  if(mass)
    {
      s *= (1 / mass);
    }

  nop->mass = mass;

  vector<MySignedIntPosType> off;
  Tp->pos_to_signedintpos(s.da, off.da);

  nop->s[0] = off[0] + nop->center[0];
  nop->s[1] = off[1] + nop->center[1];
  nop->s[2] = off[2] + nop->center[2];

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
  vector<MyReal> ms = mass * s;
  symtensor2<MyReal> ms2(ms, s);
  Q2Tensor -= ms2;
  nop->Q2Tensor = Q2Tensor;
#endif

#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
  symtensor3<MyReal> ms3(s, ms2);
  Q3Tensor -= ms3;
  Q3Tensor -= outer_prod_sum(Q2Tensor, s);
  nop->Q3Tensor = Q3Tensor;
#endif

#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
  symtensor4<MyReal> ms4(s, ms3);
  Q4Tensor -= ms4;
  Q4Tensor -= outer_prod_sum(Q3Tensor, s);
  symtensor2<MyReal> s2(s, s);
  Q4Tensor -= outer_prod_sum(Q2Tensor, s2);
  nop->Q4Tensor = Q4Tensor;
#endif

#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
  symtensor5<MyReal> ms5(s, ms4);
  Q5Tensor -= ms5;
  Q5Tensor -= outer_prod_sum(Q4Tensor, s);
  Q5Tensor -= outer_prod_sum(Q3Tensor, s2);
  symtensor3<MyReal> s3(s, s2);
  Q5Tensor -= outer_prod_sum(s3, Q2Tensor);
  nop->Q5Tensor = Q5Tensor;
#endif

#if NSOFTCLASSES > 1
  nop->maxsofttype = maxsofttype;
  nop->minsofttype = minsofttype;
#endif
  nop->cannot_be_opened_locally = 0;
  nop->not_empty                = not_empty;
#ifdef FMM
  nop->MinOldAcc = minOldAcc;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  MyIntPosType halflen = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - nop->level);
  nop->overlap_flag    = Tp->check_high_res_overlap(nop->center.da, halflen);
#endif
}

/*! \brief This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].
 *
 *  A check is performed that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
template <typename partset>
void gravtree<partset>::set_softenings(void)
{
  if(All.ComovingIntegrationOn)
    {
      for(int i = 0; i < NSOFTCLASSES; i++)
        if(All.SofteningComoving[i] * All.Time > All.SofteningMaxPhys[i])
          All.SofteningTable[i] = All.SofteningMaxPhys[i] / All.Time;
        else
          All.SofteningTable[i] = All.SofteningComoving[i];
    }
  else
    {
      for(int i = 0; i < NSOFTCLASSES; i++)
        All.SofteningTable[i] = All.SofteningComoving[i];
    }

#ifdef ADAPTIVE_HYDRO_SOFTENING
  for(int i = 0; i < NSOFTCLASSES_HYDRO; i++)
    All.SofteningTable[i + NSOFTCLASSES] = All.MinimumComovingHydroSoftening * pow(All.AdaptiveHydroSofteningSpacing, i);

  if(All.AdaptiveHydroSofteningSpacing < 1)
    Terminate("All.AdaptiveHydroSofteningSpacing < 1");

  /* we check that type=0 has its own slot 0 in the softening types, so that only gas masses are stored there */
  if(All.SofteningClassOfPartType[0] != 0)
    Terminate("All.SofteningClassOfPartType[0] != 0");

  for(int i = 1; i < NTYPES; i++)
    if(All.SofteningClassOfPartType[i] == All.SofteningClassOfPartType[0])
      Terminate("i=%d: All.SofteningClassOfPartType[i] == All.SofteningClassOfPartType[0]", i);
#endif

  for(int i = 0; i < NSOFTCLASSES + NSOFTCLASSES_HYDRO; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.ForceSoftening[NSOFTCLASSES + NSOFTCLASSES_HYDRO] =
      0; /* important - this entry is actually used in the tree construction for the search of the maximum softening in a node */
  All.ForceSoftening[NSOFTCLASSES + NSOFTCLASSES_HYDRO + 1] =
      MAX_FLOAT_NUMBER; /* important - this entry is actually used in the tree construction for the search of the maximum softening in
                           a node */
}

/* make sure that we instantiate the template */
#include "../data/simparticles.h"
template class gravtree<simparticles>;

/* make sure that we instantiate the template */
#if defined(LIGHTCONE) && (defined(LIGHTCONE_PARTICLES_GROUPS) || defined(LIGHTCONE_IMAGE_COMP_HSML_VELDISP))
#include "../data/lcparticles.h"
template class gravtree<lcparticles>;
#endif
