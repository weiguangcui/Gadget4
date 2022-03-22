/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  ngbtree_build.cc
 *
 *  \brief contains the code for the neighbor tree construction
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../gravtree/gravtree.h"
#include "../io/io.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../ngbtree/ngbtree.h"
#include "../sort/peano.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"
#include "../time_integration/timestep.h"

void ngbtree::report_log_message(void)
{
  double numnodes = NumNodes, tot_numnodes;
  MPI_Reduce(&numnodes, &tot_numnodes, 1, MPI_DOUBLE, MPI_SUM, 0, D->Communicator);

  D->mpi_printf("NGBTREE: Ngb-tree construction done. took %g sec  <numnodes>=%g  NTopnodes=%d NTopleaves=%d\n", Buildtime,
                tot_numnodes / D->NTask, D->NTopnodes, D->NTopleaves);
}

void ngbtree::fill_in_export_points(ngbpoint_data *exp_point, int i, int no) { Terminate("we don't expect to get here"); }

void ngbtree::exchange_topleafdata(void)
{
  struct leafnode_data
  {
    MySignedIntPosType center_offset_min[3];
    MySignedIntPosType center_offset_max[3];
    MyNgbTreeFloat vmin[3];
    MyNgbTreeFloat vmax[3];
    MyNgbTreeFloat MaxCsnd;
    MyNgbTreeFloat MaxHsml;
    MyNgbTreeFloat MaxDtHsml;
    unsigned char not_empty;
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
          int no       = NodeIndex[n];
          ngbnode *nop = &TopNodes[no];

          leafnode_data *locp = &loc_leaf_node_data[idx];

          locp->MaxCsnd   = nop->MaxCsnd;
          locp->MaxHsml   = nop->MaxHsml;
          locp->MaxDtHsml = nop->MaxDtHsml;
          locp->not_empty = nop->not_empty;

          for(int k = 0; k < 3; k++)
            {
              locp->center_offset_min[k] = nop->center_offset_min[k];
              locp->center_offset_max[k] = nop->center_offset_max[k];
              locp->vmin[k]              = nop->vmin[k];
              locp->vmax[k]              = nop->vmax[k];
            }

          idx++;
        }
    }

  myMPI_Allgatherv(loc_leaf_node_data, bytecounts[D->ThisTask], MPI_BYTE, glob_leaf_node_data, bytecounts, byteoffset, MPI_BYTE,
                   D->Communicator);

  for(int task = 0; task < D->NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < D->NTopleaves; n++)
    {
      int task = D->TaskOfLeaf[n];
      if(task != D->ThisTask)
        {
          int no       = NodeIndex[n];
          ngbnode *nop = &TopNodes[no];

          int idx              = recvoffset[task] + recvcounts[task]++;
          leafnode_data *globp = &glob_leaf_node_data[idx];

          nop->MaxCsnd    = globp->MaxCsnd;
          nop->MaxHsml    = globp->MaxHsml;
          nop->MaxDtHsml  = globp->MaxDtHsml;
          nop->Ti_Current = All.Ti_Current;
          nop->not_empty  = globp->not_empty;

          for(int k = 0; k < 3; k++)
            {
              nop->center_offset_min[k] = globp->center_offset_min[k];
              nop->center_offset_max[k] = globp->center_offset_max[k];
              nop->vmin[k]              = globp->vmin[k];
              nop->vmax[k]              = globp->vmax[k];
            }

          nop->Ti_Current = All.Ti_Current;
        }
    }

  Mem.myfree(loc_leaf_node_data);
  Mem.myfree(byteoffset);
  Mem.myfree(bytecounts);
  Mem.myfree(recvoffset);
  Mem.myfree(recvcounts);
  Mem.myfree(glob_leaf_node_data);
}

void ngbtree::check_bounds(void)
{
  for(int i = 0; i < Ninsert; i++)
    {
      if(Tp->P[i].get_Ti_Current() != All.Ti_Current)
        Tp->drift_particle(&Tp->P[i], &Tp->SphP[i], All.Ti_Current);  // this function avoids race conditions

      int no = Father[i];

      while(no >= 0)
        {
          ngbnode *nop = get_nodep(no);

          if(nop->level <= LEVEL_ALWAYS_OPEN)  // don't test the root node
            break;

          if(nop->Ti_Current != All.Ti_Current)
            nop->drift_node(All.Ti_Current, Tp);

          int errflag = 0;

          MyIntPosType left[3], right[3];

          left[0]  = Tp->nearest_image_intpos_to_intpos_X(nop->center_offset_min[0] + nop->center[0], Tp->P[i].IntPos[0]);
          right[0] = Tp->nearest_image_intpos_to_intpos_X(nop->center_offset_max[0] + nop->center[0], Tp->P[i].IntPos[0]);

          /* check whether we can stop walking along this branch */
          if(left[0] > 0 && right[0] > left[0])
            errflag |= 1;

          left[1]  = Tp->nearest_image_intpos_to_intpos_Y(nop->center_offset_min[1] + nop->center[1], Tp->P[i].IntPos[1]);
          right[1] = Tp->nearest_image_intpos_to_intpos_Y(nop->center_offset_max[1] + nop->center[1], Tp->P[i].IntPos[1]);

          /* check whether we can stop walking along this branch */
          if(left[1] > 0 && right[1] > left[1])
            errflag |= 2;

          left[2]  = Tp->nearest_image_intpos_to_intpos_Z(nop->center_offset_min[2] + nop->center[2], Tp->P[i].IntPos[2]);
          right[2] = Tp->nearest_image_intpos_to_intpos_Z(nop->center_offset_max[2] + nop->center[2], Tp->P[i].IntPos[2]);

          /* check whether we can stop walking along this branch */
          if(left[2] > 0 && right[2] > left[2])
            errflag |= 4;

          if(errflag)
            {
              MyIntPosType range_min[3], range_max[3];
              for(int k = 0; k < 3; k++)
                {
                  range_min[k] = nop->center_offset_min[k] + nop->center[k];
                  range_max[k] = nop->center_offset_max[k] + nop->center[k];
                }

              double pos[3], min[3], max[3];

              Tp->intpos_to_pos(Tp->P[i].IntPos, pos);
              Tp->intpos_to_pos(range_min, min);
              Tp->intpos_to_pos(range_max, max);

              Terminate(
                  "level=%d  errflag=%d  pos=%g %g %g  vel=%g %g %g    min=%g %g %g   max=%g %g %g   vmin=%g %g %g  vmax=%g %g %g    "
                  "\n",
                  nop->level, errflag, pos[0], pos[1], pos[2], Tp->P[i].Vel[0], Tp->P[i].Vel[1], Tp->P[i].Vel[2], min[0], min[1],
                  min[2], max[0], max[1], max[2], nop->vmin[0], nop->vmin[1], nop->vmin[2], nop->vmax[0], nop->vmax[1], nop->vmax[2]);
            }

          errflag = 0;

          for(int k = 0; k < 3; k++)
            {
              if(nop->vmin[k] > Tp->P[i].Vel[k])
                errflag = 1;

              if(nop->vmax[k] < Tp->P[i].Vel[k])
                errflag = 1;
            }

          if(errflag)
            {
              Terminate("vel=%g %g %g   min=%g %g %g   max=%g %g %g\n", Tp->P[i].Vel[0], Tp->P[i].Vel[1], Tp->P[i].Vel[2],
                        nop->vmin[0], nop->vmin[1], nop->vmin[2], nop->vmax[0], nop->vmax[1], nop->vmax[2]);
            }

          no = nop->father;
        }
    }
}

/*! this routine determines the node ranges a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *  mode = 0: process a leaf branch, mode = 1: process top-level nodes
 */

void ngbtree::update_node_recursive(int no, int sib, int mode)
{
  if(!(no >= MaxPart && no < MaxPart + MaxNodes)) /* are we an internal node? */
    Terminate("no internal node\n");

  ngbnode *nop = get_nodep(no);

  if(mode == TREE_MODE_TOPLEVEL)
    {
      int p = nop->nextnode;

      /* if the next node is not a top-level node, we have reached a leaf node, and we need to do nothing */
      if(p < MaxPart || p >= FirstNonTopLevelNode)
        return;
    }

  MyNgbTreeFloat maxcsnd   = 0;
  MyNgbTreeFloat maxhsml   = 0;
  MyNgbTreeFloat maxDtHsml = 0;

  MySignedIntPosType center_offset_min[3];
  MySignedIntPosType center_offset_max[3];
  MyNgbTreeFloat vmin[3], vmax[3];

  unsigned char not_empty = 0;

  MyIntPosType halflen = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - nop->level);

  for(int k = 0; k < 3; k++)
    {
      center_offset_min[k] = (halflen - 1);
      center_offset_max[k] = -halflen;

      vmin[k] = MAX_FLOAT_NUMBER;
      vmax[k] = -MAX_FLOAT_NUMBER;
    }

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
              if(maxcsnd < Tp->get_Csnd(p))
                maxcsnd = Tp->get_Csnd(p);

              if(maxhsml < Tp->SphP[p].get_Hsml())
                maxhsml = Tp->SphP[p].get_Hsml();

              if(maxDtHsml < Tp->get_DtHsml(p))
                maxDtHsml = Tp->get_DtHsml(p);

              MySignedIntPosType offset[3];

              for(int k = 0; k < 3; k++)
                {
                  offset[k] = Tp->P[p].IntPos[k] - nop->center[k];

                  if(offset[k] < center_offset_min[k])
                    center_offset_min[k] = offset[k];

                  if(offset[k] > center_offset_max[k])
                    center_offset_max[k] = offset[k];

                  if(vmin[k] > Tp->P[p].Vel[k])
                    vmin[k] = Tp->P[p].Vel[k];

                  if(vmax[k] < Tp->P[p].Vel[k])
                    vmax[k] = Tp->P[p].Vel[k];
                }

              not_empty = 1;

              p = Nextnode[p];
            }
          else if(p < MaxPart + MaxNodes) /* an internal node  */
            {
              ngbnode *noptr = get_nodep(p);

              if(maxcsnd < noptr->MaxCsnd)
                maxcsnd = noptr->MaxCsnd;

              if(maxhsml < noptr->MaxHsml)
                maxhsml = noptr->MaxHsml;

              if(maxDtHsml < noptr->MaxDtHsml)
                maxDtHsml = noptr->MaxDtHsml;

              MySignedIntPosType offset_min[3], offset_max[3];

              for(int k = 0; k < 3; k++)
                {
                  offset_min[k] = noptr->center_offset_min[k] + (MySignedIntPosType)(noptr->center[k] - nop->center[k]);
                  offset_max[k] = noptr->center_offset_max[k] + (MySignedIntPosType)(noptr->center[k] - nop->center[k]);

                  if(offset_min[k] < center_offset_min[k])
                    center_offset_min[k] = offset_min[k];

                  if(offset_max[k] > center_offset_max[k])
                    center_offset_max[k] = offset_max[k];

                  if(vmin[k] > noptr->vmin[k])
                    vmin[k] = noptr->vmin[k];

                  if(vmax[k] < noptr->vmax[k])
                    vmax[k] = noptr->vmax[k];
                }

              not_empty |= noptr->not_empty;

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

              Terminate("Ups!");
              p = Nextnode[p - MaxNodes];
            }
        }
    }

  nop->MaxCsnd   = maxcsnd;
  nop->MaxHsml   = maxhsml;
  nop->MaxDtHsml = maxDtHsml;

  nop->cannot_be_opened_locally = 0;
  nop->not_empty                = not_empty;

  for(int k = 0; k < 3; k++)
    {
      nop->center_offset_min[k] = center_offset_min[k];
      nop->center_offset_max[k] = center_offset_max[k];
      nop->vmin[k]              = vmin[k];
      nop->vmax[k]              = vmax[k];
    }

  nop->Ti_Current = All.Ti_Current;
}

void ngbtree::update_vbounds(int i, int *nchanged, int *nodelist, char *flag_changed)
{
  int no = Father[i];

  while(no >= 0)
    {
      ngbnode *nop = get_nodep(no);

      if(nop->Ti_Current != All.Ti_Current)
        nop->drift_node(All.Ti_Current, Tp);

      int has_changed = 0;

      for(int j = 0; j < 3; j++)
        {
          if(nop->vmin[j] > Tp->P[i].Vel[j])
            {
              nop->vmin[j] = Tp->P[i].Vel[j];
              has_changed  = 1;
            }

          if(nop->vmax[j] < Tp->P[i].Vel[j])
            {
              nop->vmax[j] = Tp->P[i].Vel[j];
              has_changed  = 1;
            }
        }

      if(has_changed == 0)
        break;

      if(no < FirstNonTopLevelNode) /* top-level tree-node reached */
        {
          int top_no = no - MaxPart;

          if(flag_changed[top_no] == 0)
            {
              flag_changed[top_no] = 1;

              nodelist[*nchanged] = no;
              *nchanged           = *nchanged + 1;
            }
          break;
        }

      no = nop->father;
    }
}

void ngbtree::finish_vounds_update(int nchanged, int *nodelist)
{
  int i, j, no, task, tot_nchanged;
  int *recvcounts, *recvoffset, *bytecounts, *byteoffset;
  int *tot_nodelist;
  struct leafnode_data
  {
    MyNgbTreeFloat vmin[3];
    MyNgbTreeFloat vmax[3];
  };
  leafnode_data *glob_leaf_node_data, *loc_leaf_node_data;

  /* share the pseudo-particle data accross CPUs */
  recvcounts = (int *)Mem.mymalloc("recvcounts", sizeof(int) * D->NTask);
  recvoffset = (int *)Mem.mymalloc("recvoffset", sizeof(int) * D->NTask);
  bytecounts = (int *)Mem.mymalloc("bytecounts", sizeof(int) * D->NTask);
  byteoffset = (int *)Mem.mymalloc("byteoffset", sizeof(int) * D->NTask);

  MPI_Allgather(&nchanged, 1, MPI_INT, recvcounts, 1, MPI_INT, D->Communicator);

  for(task = 0; task < D->NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(leafnode_data);

  for(task = 1, recvoffset[0] = 0, byteoffset[0] = 0; task < D->NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  loc_leaf_node_data = (leafnode_data *)Mem.mymalloc("loc_leaf_node_data", recvcounts[D->ThisTask] * sizeof(leafnode_data));

  for(i = 0; i < nchanged; i++)
    {
      for(j = 0; j < 3; j++)
        {
          loc_leaf_node_data[i].vmin[j] = get_nodep(nodelist[i])->vmin[j];
          loc_leaf_node_data[i].vmax[j] = get_nodep(nodelist[i])->vmax[j];
        }
    }

  for(task = 0, tot_nchanged = 0; task < D->NTask; task++)
    tot_nchanged += recvcounts[task];

  tot_nodelist        = (int *)Mem.mymalloc("tot_nodelist", tot_nchanged * sizeof(int));
  glob_leaf_node_data = (leafnode_data *)Mem.mymalloc("glob_leaf_node_data", tot_nchanged * sizeof(leafnode_data));

  myMPI_Allgatherv(nodelist, nchanged, MPI_INT, tot_nodelist, recvcounts, recvoffset, MPI_INT, D->Communicator);
  myMPI_Allgatherv(loc_leaf_node_data, bytecounts[D->ThisTask], MPI_BYTE, glob_leaf_node_data, bytecounts, byteoffset, MPI_BYTE,
                   D->Communicator);

  if(TreeSharedMem_ThisTask == 0) /* only one of the shared memory threads needs to update the toplevel tree */
    {
      for(i = 0; i < tot_nchanged; i++)
        {
          no = tot_nodelist[i];

          ngbnode *nop = get_nodep(no);

          if(nop->Ti_Current != All.Ti_Current)
            nop->drift_node(All.Ti_Current, Tp);

          for(j = 0; j < 3; j++)
            {
              nop->vmin[j] = glob_leaf_node_data[i].vmin[j];
              nop->vmax[j] = glob_leaf_node_data[i].vmax[j];
            }

          no = nop->father;

          while(no >= 0)
            {
              ngbnode *nop = get_nodep(no);

              if(nop->Ti_Current != All.Ti_Current)
                nop->drift_node(All.Ti_Current, Tp);

              int flag_changed = 0;

              for(j = 0; j < 3; j++)
                {
                  if(nop->vmin[j] > glob_leaf_node_data[i].vmin[j])
                    {
                      nop->vmin[j] = glob_leaf_node_data[i].vmin[j];
                      flag_changed = 1;
                    }

                  if(nop->vmax[j] < glob_leaf_node_data[i].vmax[j])
                    {
                      nop->vmax[j] = glob_leaf_node_data[i].vmax[j];
                      flag_changed = 1;
                    }
                }

              if(flag_changed == 0)
                break;

              no = nop->father;
            }
        }
    }

  Mem.myfree(glob_leaf_node_data);
  Mem.myfree(tot_nodelist);
  Mem.myfree(loc_leaf_node_data);
  Mem.myfree(byteoffset);
  Mem.myfree(bytecounts);
  Mem.myfree(recvoffset);
  Mem.myfree(recvcounts);
}

void ngbtree::finish_maxhsml_update(int nchanged, int *nodelist)
{
  int i, no, task, tot_nchanged;
  int *recvcounts, *recvoffset, *bytecounts, *byteoffset;
  int *tot_nodelist;
  struct leafnode_data
  {
    MyNgbTreeFloat MaxHsml;
    MyNgbTreeFloat MaxDtHsml;
  };
  leafnode_data *glob_leaf_node_data, *loc_leaf_node_data;

  /* share the pseudo-particle data accross CPUs */
  recvcounts = (int *)Mem.mymalloc("recvcounts", sizeof(int) * D->NTask);
  recvoffset = (int *)Mem.mymalloc("recvoffset", sizeof(int) * D->NTask);
  bytecounts = (int *)Mem.mymalloc("bytecounts", sizeof(int) * D->NTask);
  byteoffset = (int *)Mem.mymalloc("byteoffset", sizeof(int) * D->NTask);

  MPI_Allgather(&nchanged, 1, MPI_INT, recvcounts, 1, MPI_INT, D->Communicator);

  for(task = 0; task < D->NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(leafnode_data);

  for(task = 1, recvoffset[0] = 0, byteoffset[0] = 0; task < D->NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  loc_leaf_node_data = (leafnode_data *)Mem.mymalloc("loc_leaf_node_data", recvcounts[D->ThisTask] * sizeof(leafnode_data));

  for(i = 0; i < nchanged; i++)
    {
      loc_leaf_node_data[i].MaxHsml   = get_nodep(nodelist[i])->MaxHsml;
      loc_leaf_node_data[i].MaxDtHsml = get_nodep(nodelist[i])->MaxDtHsml;
    }

  for(task = 0, tot_nchanged = 0; task < D->NTask; task++)
    tot_nchanged += recvcounts[task];

  tot_nodelist        = (int *)Mem.mymalloc("tot_nodelist", tot_nchanged * sizeof(int));
  glob_leaf_node_data = (leafnode_data *)Mem.mymalloc("glob_leaf_node_data", tot_nchanged * sizeof(leafnode_data));

  myMPI_Allgatherv(nodelist, nchanged, MPI_INT, tot_nodelist, recvcounts, recvoffset, MPI_INT, D->Communicator);
  myMPI_Allgatherv(loc_leaf_node_data, bytecounts[D->ThisTask], MPI_BYTE, glob_leaf_node_data, bytecounts, byteoffset, MPI_BYTE,
                   D->Communicator);

  if(TreeSharedMem_ThisTask == 0) /* only one of the shared memory threads needs to update the toplevel tree */
    {
      for(i = 0; i < tot_nchanged; i++)
        {
          no = tot_nodelist[i];

          ngbnode *nop = get_nodep(no);

          if(nop->Ti_Current != All.Ti_Current)
            nop->drift_node(All.Ti_Current, Tp);

          nop->MaxHsml   = glob_leaf_node_data[i].MaxHsml;
          nop->MaxDtHsml = glob_leaf_node_data[i].MaxDtHsml;

          no = nop->father;

          while(no >= 0)
            {
              ngbnode *nop = get_nodep(no);

              if(nop->Ti_Current != All.Ti_Current)
                nop->drift_node(All.Ti_Current, Tp);

              if(glob_leaf_node_data[i].MaxHsml <= nop->MaxHsml && glob_leaf_node_data[i].MaxDtHsml <= nop->MaxDtHsml)
                break;
              else
                {
                  if(glob_leaf_node_data[i].MaxHsml > nop->MaxHsml)
                    nop->MaxHsml = glob_leaf_node_data[i].MaxHsml;

                  if(glob_leaf_node_data[i].MaxDtHsml > nop->MaxDtHsml)
                    nop->MaxDtHsml = glob_leaf_node_data[i].MaxDtHsml;
                }

              no = nop->father;
            }
        }
    }

  Mem.myfree(glob_leaf_node_data);
  Mem.myfree(tot_nodelist);
  Mem.myfree(loc_leaf_node_data);
  Mem.myfree(byteoffset);
  Mem.myfree(bytecounts);
  Mem.myfree(recvoffset);
  Mem.myfree(recvcounts);
}

void ngbtree::update_velocities(void)
{
  TIMER_START(CPU_NGBTREEUPDATEVEL);

  int nchanged       = 0;
  int *nodelist      = (int *)Mem.mymalloc("nodelist", D->NTopleaves * sizeof(int));
  char *flag_changed = (char *)Mem.mymalloc_clear("flag_changed", D->NTopnodes * sizeof(char));

  for(int i = 0; i < Tp->TimeBinsHydro.NActiveParticles; i++)
    {
      int target = Tp->TimeBinsHydro.ActiveParticleList[i];

      if(Tp->P[target].getType() == 0)
        update_vbounds(target, &nchanged, nodelist, flag_changed);
    }

  for(int timebin = All.HighestSynchronizedTimeBin; timebin >= 0; timebin--)
    {
      for(int target = Tp->TimeBinsHydro.FirstInTimeBin[timebin]; target >= 0; target = Tp->TimeBinsHydro.NextInTimeBin[target])
        if(Tp->P[target].getType() == 0)
          {
            update_vbounds(target, &nchanged, nodelist, flag_changed);
          }
    }

  finish_vounds_update(nchanged, nodelist);

  Mem.myfree(flag_changed);
  Mem.myfree(nodelist);

  TIMER_STOP(CPU_NGBTREEUPDATEVEL);
}

void ngbtree::update_maxhsml(void)
{
  TIMER_START(CPU_NGBTREEUPDATEMAXHSML);

  int nchanged       = 0;
  int *nodelist      = (int *)Mem.mymalloc("nodelist", D->NTopleaves * sizeof(int));
  char *flag_changed = (char *)Mem.mymalloc_clear("flag_changed", D->NTopnodes * sizeof(char));

  for(int i = 0; i < Tp->TimeBinsHydro.NActiveParticles; i++)
    {
      int target = Tp->TimeBinsHydro.ActiveParticleList[i];
      if(Tp->P[target].getType() == 0)
        {
          int no = Father[target];

          while(no >= 0)
            {
              ngbnode *nop = get_nodep(no);

              if(nop->Ti_Current != All.Ti_Current)
                nop->drift_node(All.Ti_Current, Tp);

              if(Tp->SphP[target].Hsml <= nop->MaxHsml && Tp->SphP[target].DtHsml <= nop->MaxDtHsml)
                break;
              else
                {
                  if(Tp->SphP[target].Hsml > nop->MaxHsml)
                    nop->MaxHsml = Tp->SphP[target].Hsml;

                  if(Tp->SphP[target].DtHsml > nop->MaxDtHsml)
                    nop->MaxDtHsml = Tp->SphP[target].DtHsml;
                }

              if(no < FirstNonTopLevelNode) /* top-level tree-node reached */
                {
                  int top_no = no - MaxPart;

                  if(top_no < 0 || top_no >= D->NTopnodes)
                    Terminate("top_no=%d   D->NTopleaves=%d\n", top_no, D->NTopnodes);

                  if(flag_changed[top_no] == 0)
                    {
                      flag_changed[top_no] = 1;

                      nodelist[nchanged++] = no;
                    }
                  break;
                }

              no = nop->father;
            }
        }
    }

  finish_maxhsml_update(nchanged, nodelist);

  Mem.myfree(flag_changed);
  Mem.myfree(nodelist);

  TIMER_STOP(CPU_NGBTREEUPDATEMAXHSML);
}
