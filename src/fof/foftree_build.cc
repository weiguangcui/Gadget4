/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file foftree_build.cc
 *
 *  \brief routines needed for FOF neighbor tree construction
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fof/foftree.h"
#include "../gravtree/gravtree.h"
#include "../io/io.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../sort/peano.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"
#include "../time_integration/timestep.h"

template <typename partset>
void foftree<partset>::report_log_message(void)
{
  double numnodes = NumNodes, tot_numnodes;
  MPI_Reduce(&numnodes, &tot_numnodes, 1, MPI_DOUBLE, MPI_SUM, 0, D->Communicator);

  D->mpi_printf("FOFTREE: Ngb-tree construction done. took %g sec  <numnodes>=%g  NTopnodes=%d NTopleaves=%d\n", Buildtime,
                tot_numnodes / D->NTask, D->NTopnodes, D->NTopleaves);
}

template <typename partset>
void foftree<partset>::fill_in_export_points(fofpoint_data *exp_point, int i, int no)
{
  Terminate("we don't expect to get here");
}

template <typename partset>
void foftree<partset>::exchange_topleafdata(void)
{
  struct leafnode_data
  {
    MyIntPosType range_min[3];
    MyIntPosType range_max[3];
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
          fofnode *nop = &TopNodes[no];

          leafnode_data *locp = &loc_leaf_node_data[idx];

          for(int k = 0; k < 3; k++)
            {
              locp->range_min[k] = nop->range_min[k];
              locp->range_max[k] = nop->range_max[k];
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
          fofnode *nop = &TopNodes[no];

          int idx              = recvoffset[task] + recvcounts[task]++;
          leafnode_data *globp = &glob_leaf_node_data[idx];

          for(int k = 0; k < 3; k++)
            {
              nop->range_min[k] = globp->range_min[k];
              nop->range_max[k] = globp->range_max[k];
            }
        }
    }

  Mem.myfree(loc_leaf_node_data);
  Mem.myfree(byteoffset);
  Mem.myfree(bytecounts);
  Mem.myfree(recvoffset);
  Mem.myfree(recvcounts);
  Mem.myfree(glob_leaf_node_data);
}

/*! this routine determines the node ranges a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *  mode = 0: process a leaf branch, mode = 1: process top-level nodes
 */
template <typename partset>
void foftree<partset>::update_node_recursive(int no, int sib, int mode)
{
  MyIntPosType range_min[3];
  MyIntPosType range_max[3];

  if(!(no >= MaxPart && no < MaxPart + MaxNodes)) /* are we an internal node? */
    Terminate("no internal node\n");

  fofnode *nop = get_nodep(no);

  if(mode == TREE_MODE_TOPLEVEL)
    {
      int p = nop->nextnode;

      /* if the next node is not a top-level node, we have reached a leaf node, and we need to do nothing */
      if(p < MaxPart || p >= FirstNonTopLevelNode)
        return;
    }

  for(int k = 0; k < 3; k++)
    {
      range_min[k] = ~((MyIntPosType)0);
      range_max[k] = 0;
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
              for(int k = 0; k < 3; k++)
                {
                  if(range_min[k] > Tp->P[p].IntPos[k])
                    range_min[k] = Tp->P[p].IntPos[k];

                  if(range_max[k] < Tp->P[p].IntPos[k])
                    range_max[k] = Tp->P[p].IntPos[k];
                }

              p = Nextnode[p];
            }
          else if(p < MaxPart + MaxNodes) /* an internal node  */
            {
              fofnode *nop = get_nodep(p);

              for(int k = 0; k < 3; k++)
                {
                  if(range_min[k] > nop->range_min[k])
                    range_min[k] = nop->range_min[k];

                  if(range_max[k] < nop->range_max[k])
                    range_max[k] = nop->range_max[k];
                }

              p = nop->sibling;
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

  for(int k = 0; k < 3; k++)
    {
      nop->range_min[k] = range_min[k];
      nop->range_max[k] = range_max[k];
    }
}

/* make sure that we instantiate the template */
#include "../data/simparticles.h"
template class foftree<simparticles>;

#if defined(LIGHTCONE) && (defined(LIGHTCONE_PARTICLES_GROUPS) || defined(LIGHTCONE_IMAGE_COMP_HSML_VELDISP))
/* make sure that we instantiate the template */
#include "../data/lcparticles.h"
template class foftree<lcparticles>;
#endif
