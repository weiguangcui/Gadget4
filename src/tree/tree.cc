/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file tree.cc
 *
 *  \brief basic routines for oct-tree building
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
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/cxxsort.h"
#include "../sort/peano.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"
#include "../tree/tree.h"

/*! This file contains the construction of the tree used for calculating the gravitational force
 *  and the neighbor tree for SPH.
 *  The type of tree implemented is a geometrical oct-tree, starting from a cube encompassing
 *  all particles. This cube is automatically found in the domain decomposition, which also
 *  splits up the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. In the present version of the code, the tree
 *  construction may be repeated without a renewed domain decomposition. In this case,
 *  if particles are on the "wrong" processor because a new domain decomposition has not been
 *  carried out, they are sent as temporary points to the right insertion processor according
 *  to the layout of the top-level nodes.
 */

/*! This function is a driver routine for constructing the oct-tree.
 *
 *  \return number of local nodes (including top level nodes) of the constructed tree
 */
template <typename node, typename partset, typename point_data, typename foreign_point_data>
int tree<node, partset, point_data, foreign_point_data>::treebuild(int ninsert, int *indexlist)
{
  if(MaxPart == 0 && ninsert == 0)
    return 0;  // nothing to be done

  if(MaxPart == 0 && ninsert > 0)
    Terminate("Strange, we'll try to construct a tree for %d particles, but it appears not be allocated\n", ninsert);

  if(ninsert == Tp->NumPart)
    D->mpi_printf("TREE: Full tree construction for all particles. (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());

  Ninsert   = ninsert;
  IndexList = indexlist;

  TIMER_START(CPU_TREEBUILD);

  double t0 = Logs.second();

  int flag, iter = 0;
  do /* try constructing tree until successful in terms of storage allocation */
    {
      TIMER_START(CPU_TREEBUILD_INSERT);

      int flag_single = treebuild_construct();

      TIMER_STOP(CPU_TREEBUILD_INSERT);

      MPI_Allreduce(&flag_single, &flag, 1, MPI_INT, MPI_MIN, D->Communicator);

      if(flag < 0)
        {
          /* tree construction was not successful and needs to be repeated */
          treefree();
          All.TreeAllocFactor *= 1.15;
          // D->mpi_printf("TREE: Increasing TreeAllocFactor, new value=%g\n", All.TreeAllocFactor);
          treeallocate(Tp->NumPart, Tp, D);
        }
      else
        {
          /* treebuild was successfull, but let's check if we allocated clearly too much storage, and if so repeat it */
          int max_numnodes;
          MPI_Allreduce(&NumNodes, &max_numnodes, 1, MPI_INT, MPI_MAX, D->Communicator);

          if((MaxNodes - D->NTopnodes) > 1.5 * (max_numnodes - D->NTopnodes))
            {
              double oldvalue = All.TreeAllocFactor;
              double newvalue = std::max<double>(1.1 * (max_numnodes - D->NTopnodes) / (MaxPart + BASENUMBER), 0.02);

              if(newvalue < oldvalue)
                {
                  treefree();

                  All.TreeAllocFactor = newvalue;
                  flag                = -1;
                  treeallocate(Tp->NumPart, Tp, D);
                }
            }
        }
      iter++;

      if(iter > TREE_MAX_ITER)
        Terminate("tree construction failed\n");
    }
  while(flag < 0);

  TIMER_STOPSTART(CPU_TREEBUILD, CPU_TREEBUILD_BRANCHES);

  /* first, construct the properties of the tree branches below the top leaves */

  int ntopleaves = D->NumTopleafOfTask[D->ThisTask];
  int *list      = D->ListOfTopleaves + D->FirstTopleafOfTask[D->ThisTask];

  for(int i = 0; i < ntopleaves; i++)
    {
      int no  = NodeIndex[list[i]];
      int sib = NodeSibling[list[i]];

      update_node_recursive(no, sib, TREE_MODE_BRANCH);
    }

  TIMER_STOPSTART(CPU_TREEBUILD_BRANCHES, CPU_TREEBUILD_TOPLEVEL);

  exchange_topleafdata();

  /* now update the top-level tree nodes */
  if(TreeSharedMem_ThisTask == 0)
    update_node_recursive(MaxPart, -1, TREE_MODE_TOPLEVEL);

  for(int k = D->NTopnodes; k < NumNodes; k++)
    {
      int index               = MaxPart + k;
      Nodes[index].OriginTask = D->ThisTask;
      Nodes[index].OriginNode = index;
      Nodes[index].access.clear();
      Nodes[index].flag_already_fetched = 0;
    }

  if(TreeSharedMem_ThisTask == 0)
    {
      for(int k = 0; k < D->NTopnodes; k++)
        {
          int index                  = MaxPart + k;
          TopNodes[index].OriginTask = D->ThisTask;
          TopNodes[index].OriginNode = index;
          TopNodes[index].access.clear();
          TopNodes[index].flag_already_fetched = 0;
        }
    }

  tree_initialize_leaf_node_access_info();

  double t1 = Logs.second();
  Buildtime = Logs.timediff(t0, t1);

  report_log_message();

  TIMER_STOP(CPU_TREEBUILD_TOPLEVEL);

  return NumNodes;
}

template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::tree_initialize_leaf_node_access_info(void)
{
  if(TreeSharedMem_ThisTask == 0)
    {
      for(int k = 0; k < D->NTopleaves; k++)
        {
          int index = NodeIndex[k];

          if(Shmem.GetNodeIDForSimulCommRank[D->TaskOfLeaf[k]] == Shmem.GetNodeIDForSimulCommRank[D->ThisTask])
            TopNodes[index].cannot_be_opened_locally = 0;
          else
            {
              TopNodes[index].cannot_be_opened_locally = 1;
              TopNodes[index].flag_already_fetched     = 0;
              TopNodes[index].nextnode                 = MaxPart + MaxNodes + k;
              TopNodes[index].nextnode_shmrank         = TreeSharedMem_ThisTask;
            }
          TopNodes[index].OriginTask = D->TaskOfLeaf[k];
          TopNodes[index].OriginNode = index;
        }
    }
}

/*! Constructs the gravitational oct-tree.
 *
 *  The index convention for accessing tree nodes is the following: \n
 *  node index \n
 *  [0...            MaxPart-1]                references single particles, the indices \n
 *  [MaxPart... MaxPart+MaxNodes-1]  references tree nodes \n
 *  [MaxPart+MaxNodes...  MaxPart+MaxNodes+NTopleaves-1]                                  references "pseudo particles", i.e. markers
 * for branches on foreign CPUs \n [MaxPart+MaxNodes+NTopleaves... MaxPart+MaxNodes+NTopleaves+NumPartImported-1]   references imported
 * points \n
 *
 *  the pointer `Nodes' is shifted such that Nodes[MaxPart] gives the first tree node (i.e. the root node).
 *
 *  \return if successful returns the number of local+top nodes of the constructed tree \n
 *          -1 if the number of allocated tree nodes is too small \n
 *          -2 if the number of allocated tree nodes is even too small to fit the top nodes \n
 *          -3 if a particle out of domain box condition was encountered
 */
template <typename node, typename partset, typename point_data, typename foreign_point_data>
int tree<node, partset, point_data, foreign_point_data>::treebuild_construct(void)
{
  if(TreeSharedMem_ThisTask == 0)
    {
      /* create an empty root node  */
      NextFreeNode = MaxPart; /* index of first free node */

      node *nfreep             = &TopNodes[NextFreeNode]; /* select first node        */
      nfreep->nextnode         = -1;
      nfreep->sibling          = -1;
      nfreep->father           = -1;
      nfreep->level            = 0;
      nfreep->sibling_shmrank  = TreeSharedMem_ThisTask;
      nfreep->nextnode_shmrank = TreeSharedMem_ThisTask;

      for(int j = 0; j < 3; j++)
        nfreep->center[j] = ((MyIntPosType)1) << (BITS_FOR_POSITIONS - 1);

      NodeIndex[0]   = NextFreeNode;
      NodeLevel[0]   = 0;
      NodeSibling[0] = -1;
      NumNodes       = 1;
      NextFreeNode++;

      /* create a set of empty nodes corresponding to the top-level domain
       * grid. We need to generate these nodes first to make sure that we have a
       * complete top-level tree which allows the easy insertion of the
       * pseudo-particles at the right place.
       */

      if(create_empty_nodes(MaxPart, 1, 0, 1, -1, 0, 0, 0) < 0)
        return -2;
    }

  MPI_Bcast(&NumNodes, 1, MPI_INT, 0, D->Communicator);
  MPI_Bcast(&NextFreeNode, 1, MPI_INT, 0, D->Communicator);

  if(NumNodes != D->NTopnodes)
    Terminate("NumNodes=%d != D->NTopnodes=%d", NumNodes, D->NTopnodes);

  FirstNonTopLevelNode = NextFreeNode;

  if(MaxPart + D->NTopnodes != FirstNonTopLevelNode)
    Terminate("unexpected");

  for(int j = 0; j < D->NTask; j++)
    Send_count[j] = 0;

#ifndef LEAN
  int *node_list = (int *)Mem.mymalloc_movable(&node_list, "node_list", Tp->NumPart * sizeof(int));
  int *task_list = (int *)Mem.mymalloc_movable(&task_list, "task_list", Tp->NumPart * sizeof(int));
#endif

  /* now we determine for each point the insertion top-level node, and the task on which this lies */
  for(int idx = 0; idx < Ninsert; idx++)
    {
      int i;
      if(IndexList)
        i = IndexList[idx];
      else
        i = idx;

      if(Tp->P[i].get_Ti_Current() != All.Ti_Current)
        Tp->drift_particle(&Tp->P[i], &Tp->SphP[i], All.Ti_Current);

      int no;
      int task;
      tree_get_node_and_task(i, no, task);

#ifndef LEAN
      node_list[i] = no;
      task_list[i] = task;
#endif

      if(task != D->ThisTask)
        Send_count[task]++;
    }

  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, D->Communicator);

  NumPartImported = 0;
  NumPartExported = 0;
  Recv_offset[0]  = 0;
  Send_offset[0]  = 0;

  for(int j = 0; j < D->NTask; j++)
    {
      NumPartImported += Recv_count[j];
      NumPartExported += Send_count[j];
      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  Points   = (point_data *)Mem.mymalloc_movable(&Points, "Points", NumPartImported * sizeof(point_data));
  Nextnode = (int *)Mem.mymalloc_movable(&Nextnode, "Nextnode", (MaxPart + D->NTopleaves + NumPartImported) * sizeof(int));
  Father   = (int *)Mem.mymalloc_movable(&Father, "Father", (MaxPart + NumPartImported) * sizeof(int));

  /* now put in markers ("pseudo" particles) in top-leaf nodes to indicate on which task the branch lies */
  for(int i = 0; i < D->NTopleaves; i++)
    {
      int index = NodeIndex[i];

      if(TreeSharedMem_ThisTask == 0)
        TopNodes[index].nextnode = MaxPart + MaxNodes + i;

      /* set nextnode for pseudo-particle (Nextnode exists on all ranks) */
      Nextnode[MaxPart + i] = TopNodes[index].sibling;
    }

  point_data *export_Points = (point_data *)Mem.mymalloc("export_Points", NumPartExported * sizeof(point_data));

  for(int j = 0; j < D->NTask; j++)
    Send_count[j] = 0;

  for(int idx = 0; idx < Ninsert; idx++) /* prepare particle data to be copied to other tasks */
    {
      int i;
      if(IndexList)
        i = IndexList[idx];
      else
        i = idx;

#ifdef LEAN
      int no, task;
      tree_get_node_and_task(i, no, task);
#else
      int task = task_list[i];
      int no   = node_list[i];
#endif

      if(task != D->ThisTask)
        {
          int n = Send_offset[task] + Send_count[task]++;

          fill_in_export_points(&export_Points[n], i, no);
        }
    }

  /* exchange data */
  for(int ngrp = 1; ngrp < (1 << D->PTask); ngrp++)
    {
      int recvTask = D->ThisTask ^ ngrp;
      if(recvTask < D->NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_Points[Send_offset[recvTask]], Send_count[recvTask] * sizeof(point_data), MPI_BYTE, recvTask,
                         TAG_DENS_A, &Points[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(point_data), MPI_BYTE, recvTask,
                         TAG_DENS_A, D->Communicator, MPI_STATUS_IGNORE);
    }

  Mem.myfree(export_Points);

  ImportedNodeOffset = MaxPart + MaxNodes + D->NTopleaves;

  MPI_Allreduce(&NumPartImported, &EndOfTreePoints, 1, MPI_INT, MPI_MAX, D->Communicator);
  EndOfTreePoints += ImportedNodeOffset;
  EndOfForeignNodes = EndOfTreePoints + (INT_MAX - EndOfTreePoints) / 2;

  /* make a list that holds the particles that belong to a certain node */
  index_data *index_list =
      (index_data *)Mem.mymalloc_movable(&index_list, "index_list", (Ninsert + NumPartImported) * sizeof(index_data));
  int count = 0;

  for(int idx = 0; idx < Ninsert; idx++)
    {
      int i;
      if(IndexList)
        i = IndexList[idx];
      else
        i = idx;

#ifdef LEAN
      int no, task;
      tree_get_node_and_task(i, no, task);
#else
      int task = task_list[i];
      int no   = node_list[i];
#endif

      if(task == D->ThisTask)
        {
          index_list[count].p       = i;
          index_list[count].subnode = no;
          count++;
        }
    }

  for(int i = 0; i < NumPartImported; i++)
    {
      index_list[count].p       = i + ImportedNodeOffset;
      index_list[count].subnode = Points[i].no;
      count++;
    }

#ifndef LEAN
  Mem.myfree_movable(task_list);
  Mem.myfree_movable(node_list);
#endif

  /* sort according to node so that particles indices in the same node are grouped together */
  mycxxsort(index_list, index_list + count, compare_index_data_subnode);

  int full_flag  = 0;
  int ntopleaves = D->NumTopleafOfTask[D->ThisTask];
  int start      = 0;

  for(int n = 0; n < ntopleaves; n++)
    {
      int no              = D->ListOfTopleaves[D->FirstTopleafOfTask[D->ThisTask] + n];
      int th              = NodeIndex[no];
      unsigned char level = NodeLevel[no];
      int sibling         = NodeSibling[no];

      while(start < count && index_list[start].subnode < no)
        start++;

      if(start < count)
        {
          int last = start;
          while(last < count && index_list[last].subnode == no)
            last++;

          int num = last - start;

          if(treebuild_insert_group_of_points(num, &index_list[start], th, level, sibling))
            {
              full_flag = 1;
              break;
            }

          start += num;
        }
    }

  if((NumNodes = NextFreeNode - MaxPart) >= MaxNodes)
    {
      if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
        {
          Tp->dump_particles();
          Terminate(
              "task %d: looks like a serious problem, stopping with particle dump.  NumNodes=%d MaxNodes=%d  NumPartImported=%d "
              "NumPart=%d\n",
              D->ThisTask, NumNodes, MaxNodes, NumPartImported, Tp->NumPart);
        }
    }

  Mem.myfree(index_list);

  if(full_flag)
    return -1;

  return NumNodes;
}

/*! inserts a group of particles into the gravitational tree
 *
 *  level - level of target node
 *  th    - target node
 *
 *  \return 0 if successful \n
 *          1 if too few nodes have been allocated in the Nodes array
 */
template <typename node, typename partset, typename point_data, typename foreign_point_data>
int tree<node, partset, point_data, foreign_point_data>::treebuild_insert_group_of_points(int num, index_data *index_list, int th,
                                                                                          unsigned char level, int sibling)
{
  if(level >= BITS_FOR_POSITIONS)
    Terminate(
        "It appears we have reached the bottom of the tree because there are more than TREE_NUM_BEFORE_NODESPLIT=%d particles in the "
        "smallest tree node representable for  BITS_FOR_POSITIONS=%d.\n"
        "Either eliminate the particles at (nearly) indentical coordinates, increase the setting for TREE_NUM_BEFORE_NODESPLIT, or "
        "possibly enlarge BITS_FOR_POSITIONS if you have really not enough dynamic range\n",
        (int)TREE_NUM_BEFORE_NODESPLIT, (int)BITS_FOR_POSITIONS);

  MyIntPosType mask       = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - 1 - level));
  unsigned char shiftx    = (BITS_FOR_POSITIONS - 1 - level);
  unsigned char shifty    = (BITS_FOR_POSITIONS - 2 - level);
  unsigned char shiftz    = (BITS_FOR_POSITIONS - 3 - level);
  MyIntPosType centermask = ~(~((MyIntPosType)0) >> level);

  int subcount[8] = {0, 0, 0, 0, 0, 0, 0, 0}, subnode[8];
  MyIntPosType *subintpos[8];

  for(int i = 0; i < num; i++)
    {
      MyIntPosType *intpos;

      int p = index_list[i].p;
      if(p < MaxPart)
        {
          intpos = Tp->P[p].IntPos;
        }
      else
        {
          int n  = p - ImportedNodeOffset;
          intpos = Points[n].IntPos;
        }

      unsigned char subnode = (((unsigned char)((intpos[0] & mask) >> shiftx)) | ((unsigned char)((intpos[1] & mask) >> shifty)) |
                               ((unsigned char)((intpos[2] & mask) >> shiftz)));
      if(subnode > 7)
        Terminate("stop: subnode > 7");

      subcount[subnode]++;

      subintpos[subnode] = intpos;

      index_list[i].subnode = subnode;
    }

  /* sort */
  mycxxsort(index_list, index_list + num, compare_index_data_subnode);

  centermask >>= 1;
  centermask |= ~(~((MyIntPosType)0) >> 1); /* this sets the MSB */
  mask >>= 1;

  node *nfreep_last = NULL;

  /* create the daughter nodes */
  for(int i = 0; i < 8; i++)
    {
      if(subcount[i] > TREE_NUM_BEFORE_NODESPLIT)
        {
          int thnew;

          thnew = NextFreeNode++;

          if(thnew - MaxPart >= MaxNodes)
            return 1; /* we are out of space */

          subnode[i] = thnew;

          if(nfreep_last)
            {
              nfreep_last->sibling  = thnew;
              nfreep_last->nextnode = thnew;

              nfreep_last->sibling_shmrank  = TreeSharedMem_ThisTask;
              nfreep_last->nextnode_shmrank = TreeSharedMem_ThisTask;
            }
          else
            {
              get_nodep(th)->nextnode         = thnew;
              get_nodep(th)->nextnode_shmrank = TreeSharedMem_ThisTask;
            }

          if(thnew < MaxPart + D->NTopnodes)
            Terminate("thnew = %d   <  MaxPart=%d + D->NTopnodes=%d", thnew, MaxPart, D->NTopnodes);

          node *nfreep = get_nodep(thnew);

          nfreep->father = th;
          nfreep->level  = level + 1;

          nfreep->center[0] = ((subintpos[i][0] & centermask) | mask);
          nfreep->center[1] = ((subintpos[i][1] & centermask) | mask);
          nfreep->center[2] = ((subintpos[i][2] & centermask) | mask);

          nfreep->nextnode_shmrank = TreeSharedMem_ThisTask;
          nfreep->sibling_shmrank  = TreeSharedMem_ThisTask;

          nfreep_last = nfreep;
        }
    }

  /* now insert the particles that are chained to the node */

  int p_last = -1, p_first = -1;

  /* now insert the particle groups into the created daughter nodes, or chain them to the node */
  for(int i = 0, off = 0; i < 8; i++)
    {
      if(subcount[i] <= TREE_NUM_BEFORE_NODESPLIT)
        {
          /* put the particles into the node as a chain */
          for(int j = 0; j < subcount[i]; j++)
            {
              int p = index_list[off + j].p;

              if(nfreep_last == NULL && p_first == -1)
                {
                  get_nodep(th)->nextnode         = p;
                  get_nodep(th)->nextnode_shmrank = TreeSharedMem_ThisTask;
                }

              if(p < MaxPart)
                Father[p] = th;
              else
                Father[p - MaxNodes - D->NTopleaves] = th;

              if(p_last >= 0)
                {
                  if(p_last < MaxPart)
                    Nextnode[p_last] = p;
                  else
                    Nextnode[p_last - MaxNodes] = p; /* imported point */
                }

              p_last = p;

              if(p_first == -1)
                p_first = p;
            }
        }

      off += subcount[i];
    }

  if(p_last >= 0)
    {
      if(p_last < MaxPart)
        Nextnode[p_last] = sibling;
      else
        Nextnode[p_last - MaxNodes] = sibling; /* imported point */
    }

  if(nfreep_last)
    {
      if(p_first >= 0)
        {
          nfreep_last->sibling  = p_first;
          nfreep_last->nextnode = p_first;

          nfreep_last->sibling_shmrank  = TreeSharedMem_ThisTask;
          nfreep_last->nextnode_shmrank = TreeSharedMem_ThisTask;
        }
      else
        {
          nfreep_last->sibling  = sibling;
          nfreep_last->nextnode = sibling;

          nfreep_last->sibling_shmrank  = TreeSharedMem_ThisTask;
          nfreep_last->nextnode_shmrank = TreeSharedMem_ThisTask;
        }
    }

  for(int i = 0, off = 0; i < 8; i++)
    {
      if(subcount[i] > TREE_NUM_BEFORE_NODESPLIT)
        {
          if(subnode[i] < MaxPart + D->NTopnodes)
            Terminate("subnode[i]=%d < MaxPart=%d + D->NTopnodes=%d", subnode[i], MaxPart, D->NTopnodes);

          int out_of_space =
              treebuild_insert_group_of_points(subcount[i], &index_list[off], subnode[i], level + 1, get_nodep(subnode[i])->sibling);
          if(out_of_space)
            return out_of_space;
        }

      off += subcount[i];
    }

  return 0; /* success */
}

/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 *
 * \return 0 if successful \n
 *         -1 if number of allocated tree nodes is too small to fit the newly created nodes
 */
template <typename node, typename partset, typename point_data, typename foreign_point_data>
int tree<node, partset, point_data, foreign_point_data>::create_empty_nodes(
    int no_parent /*!< parent node for which daughter nodes shall be created */, int level /*!< level of new nodes */,
    int topnode /*!< index of the parent node in the TopNodes array */,
    int bits /*!< 2^bits is the number of nodes per dimension at the level of the daughter nodes */, int sibling,
    MyIntPosType x /*!< position of the parent node in the x direction, in the range [0,2^(bits-1) - 1] */,
    MyIntPosType y /*!< position of the parent node in the y direction, in the range [0,2^(bits-1) - 1] */,
    MyIntPosType z /*!< position of the parent node in the z direction, in the range [0,2^(bits-1) - 1] */)
{
  if(D->TopNodes[topnode].Daughter >= 0)
    {
      int firstflag = 0;
      int nostart   = NextFreeNode;

      /* loop over daughter nodes */
      for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
          for(int k = 0; k < 2; k++)
            {
              if(NumNodes >= MaxNodes)
                {
                  if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                    {
                      char buf[MAXLEN_PATH_EXTRA];
                      snprintf(buf, MAXLEN_PATH_EXTRA,
                               "task %d: looks like a serious problem (NTopnodes=%d), stopping with particle dump.\n", D->ThisTask,
                               D->NTopnodes);
                      Tp->dump_particles();
                      Terminate(buf);
                    }
                  return -1;
                }

              int no = NextFreeNode++;
              NumNodes++;

              if(firstflag == 0)
                {
                  TopNodes[no_parent].nextnode = no;
                  firstflag                    = 1;
                }

              TopNodes[no].father = no_parent;

              if(i + 2 * j + 4 * k == 7)
                TopNodes[no].sibling = sibling;
              else
                TopNodes[no].sibling = no + 1;

              TopNodes[no].nextnode = TopNodes[no].sibling;

              TopNodes[no].level = level;

              MyIntPosType lenhalf   = ((MyIntPosType)1) << (BITS_FOR_POSITIONS - level - 1);
              TopNodes[no].center[0] = TopNodes[no_parent].center[0] + (2 * i - 1) * lenhalf;
              TopNodes[no].center[1] = TopNodes[no_parent].center[1] + (2 * j - 1) * lenhalf;
              TopNodes[no].center[2] = TopNodes[no_parent].center[2] + (2 * k - 1) * lenhalf;

              TopNodes[no].sibling_shmrank  = TreeSharedMem_ThisTask;
              TopNodes[no].nextnode_shmrank = TreeSharedMem_ThisTask;
            }

      /* loop over daughter nodes */
      for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
          for(int k = 0; k < 2; k++)
            {
              int no = nostart++;

              peanokey key = peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);
              int sub      = 7 & key.ls;

              if(D->TopNodes[D->TopNodes[topnode].Daughter + sub].Daughter == -1)
                {
                  NodeIndex[D->TopNodes[D->TopNodes[topnode].Daughter + sub].Leaf]   = no;
                  NodeLevel[D->TopNodes[D->TopNodes[topnode].Daughter + sub].Leaf]   = level;
                  NodeSibling[D->TopNodes[D->TopNodes[topnode].Daughter + sub].Leaf] = TopNodes[no].sibling;
                }

              /* create grand daughter nodes for current daughter node */
              if(create_empty_nodes(no, level + 1, D->TopNodes[topnode].Daughter + sub, bits + 1, TopNodes[no].sibling, 2 * x + i,
                                    2 * y + j, 2 * z + k) < 0)
                return -1;
            }
    }

  return 0;
}

/*! This function allocates the memory used for storage of the tree nodes. Usually,
 *  the number of required nodes is of order 0.7*maxpart, but if this is insufficient,
 *  the code will try to allocated more space by increasing TreeAllocFactor.
 */
template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::treeallocate(int max_partindex, partset *Tp_ptr, domain<partset> *Dptr)
{
  D  = Dptr;
  Tp = Tp_ptr;

  /* split up the communicator into pieces overlap with different shared memory regions */
  if(max_partindex != -1)
    MPI_Comm_split(D->Communicator, Shmem.Island_Smallest_WorldTask, 0, &TreeSharedMemComm);

  if(max_partindex != -1)
    MPI_Allreduce(&max_partindex, &MaxPart, 1, MPI_INT, MPI_MAX, D->Communicator);

  if(MaxPart == 0)
    return;  // nothing to be done

  if(Nodes)
    Terminate("Nodes already allocated");

  Send_count  = (int *)Mem.mymalloc_movable(&Send_count, "Send_count", sizeof(int) * D->NTask);
  Send_offset = (int *)Mem.mymalloc_movable(&Send_offset, "Send_offset", sizeof(int) * D->NTask);
  Recv_count  = (int *)Mem.mymalloc_movable(&Recv_count, "Recv_count", sizeof(int) * D->NTask);
  Recv_offset = (int *)Mem.mymalloc_movable(&Recv_offset, "Recv_offset", sizeof(int) * D->NTask);

  if(max_partindex != -1)
    {
      MPI_Allreduce(MPI_IN_PLACE, &All.TreeAllocFactor, 1, MPI_DOUBLE, MPI_MAX, D->Communicator);

      MaxNodes = (int)(All.TreeAllocFactor * (MaxPart + BASENUMBER)) + D->NTopnodes;

      int max_nodes;
      MPI_Allreduce(&MaxNodes, &max_nodes, 1, MPI_INT, MPI_MAX, D->Communicator);

      if(max_nodes != MaxNodes)
        Terminate("Strange: different maxnodes detected: %d %d", max_nodes, MaxNodes);

      int max_leaves;
      MPI_Allreduce(&D->NTopleaves, &max_leaves, 1, MPI_INT, MPI_MAX, D->Communicator);
      if(max_leaves != D->NTopleaves)
        Terminate("Strange: different maxnodes detected: %d %d", max_leaves, D->NTopleaves);
    }
  else
    {
      max_partindex = MaxPart;
    }

  MPI_Comm_rank(TreeSharedMemComm, &TreeSharedMem_ThisTask);
  MPI_Comm_size(TreeSharedMemComm, &TreeSharedMem_NTask);

  TreeNodes_offsets          = (ptrdiff_t *)Mem.mymalloc("TreeNodes_offsets", TreeSharedMem_NTask * sizeof(ptrdiff_t));
  TreePoints_offsets         = (ptrdiff_t *)Mem.mymalloc("TreePoints_offsets", TreeSharedMem_NTask * sizeof(ptrdiff_t));
  TreeNextnode_offsets       = (ptrdiff_t *)Mem.mymalloc("TreeNextnode_offsets", TreeSharedMem_NTask * sizeof(ptrdiff_t));
  TreeForeign_Nodes_offsets  = (ptrdiff_t *)Mem.mymalloc("TreeForeign_Nodes_offsets", TreeSharedMem_NTask * sizeof(ptrdiff_t));
  TreeForeign_Points_offsets = (ptrdiff_t *)Mem.mymalloc("TreeForeign_Points_offsets", TreeSharedMem_NTask * sizeof(ptrdiff_t));
  TreeP_offsets              = (ptrdiff_t *)Mem.mymalloc("TreeP_offsets", TreeSharedMem_NTask * sizeof(ptrdiff_t));
  TreeSphP_offsets           = (ptrdiff_t *)Mem.mymalloc("TreeSphP_offsets", TreeSharedMem_NTask * sizeof(ptrdiff_t));
  TreePS_offsets             = (ptrdiff_t *)Mem.mymalloc("TreePS_offsets", TreeSharedMem_NTask * sizeof(ptrdiff_t));

  TreeSharedMemBaseAddr = (void **)Mem.mymalloc("TreeSharedMemBaseAddr", TreeSharedMem_NTask * sizeof(void *));

  for(int i = 0; i < TreeSharedMem_NTask; i++)
    {
      int island_rank = Shmem.Island_ThisTask + (i - TreeSharedMem_ThisTask);

      if(island_rank < 0 || island_rank >= Shmem.Island_NTask)
        Terminate("island_rank=%d  < 0 || island_rank >= Shmem.Island_NTask=%d", island_rank, Shmem.Island_NTask);

      TreeSharedMemBaseAddr[i] = Shmem.SharedMemBaseAddr[island_rank];
    }

  /* allocate the Top-Level tree only once per shared-memory section in the communicator */
  if(TreeSharedMem_ThisTask == 0)
    {
      /* if we have a single shared memory node, or we are building a tree on a subcommunicator, allocate locally */
      if(TreeSharedMem_NTask == Shmem.World_NTask || D->NTask < Shmem.Sim_NTask)
        {
          NodeLevel   = (unsigned char *)Mem.mymalloc("NodeLevel", D->NTopleaves * sizeof(unsigned char));
          NodeSibling = (int *)Mem.mymalloc("NodeSibling", D->NTopleaves * sizeof(int));
          NodeIndex   = (int *)Mem.mymalloc("NodeIndex", D->NTopleaves * sizeof(int));

          TopNodes = (node *)Mem.mymalloc("TopNodes", D->NTopnodes * sizeof(node));
          TopNodes -= MaxPart;
        }
      else /* otherwise, allocate the storage on the ghost processor, and get the address in our space from him */
        {
          int ghost_rank = Shmem.GetGhostRankForSimulCommRank[Shmem.Sim_ThisTask];

          /* request the storage from the responsible ghost rank, and map it into the local address space */
          size_t tab_len[4] = {D->NTopleaves * sizeof(unsigned char), D->NTopleaves * sizeof(int), D->NTopleaves * sizeof(int),
                               D->NTopnodes * sizeof(node)};

          MPI_Send(tab_len, 4 * sizeof(tab_len), MPI_BYTE, ghost_rank, TAG_TOPNODE_ALLOC, MPI_COMM_WORLD);

          ptrdiff_t off[4];
          MPI_Recv(off, 4 * sizeof(ptrdiff_t), MPI_BYTE, ghost_rank, TAG_TOPNODE_OFFSET, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          NodeLevel   = (unsigned char *)((char *)Shmem.SharedMemBaseAddr[Shmem.Island_NTask - 1] + off[0]);
          NodeSibling = (int *)((char *)Shmem.SharedMemBaseAddr[Shmem.Island_NTask - 1] + off[1]);
          NodeIndex   = (int *)((char *)Shmem.SharedMemBaseAddr[Shmem.Island_NTask - 1] + off[2]);

          TopNodes = (node *)((char *)Shmem.SharedMemBaseAddr[Shmem.Island_NTask - 1] + off[3]);
          TopNodes -= MaxPart;

          MPI_Recv(&TreeInfoHandle, 1, MPI_INT, ghost_rank, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

  Nodes = (node *)Mem.mymalloc_movable(&Nodes, "Nodes", (MaxNodes - D->NTopnodes + 1) * sizeof(node));
  Nodes -= (MaxPart + D->NTopnodes);

  if(max_partindex != -1)
    treeallocate_share_topnode_addresses();
}

template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::treeallocate_share_topnode_addresses(void)
{
  MPI_Bcast(&TreeInfoHandle, 1, MPI_INT, 0, TreeSharedMemComm);

  ptrdiff_t off[4] = {((char *)NodeLevel - Mem.Base), ((char *)NodeSibling - Mem.Base), ((char *)NodeIndex - Mem.Base),
                      ((char *)TopNodes - Mem.Base)};

  MPI_Bcast(off, 4 * sizeof(ptrdiff_t), MPI_BYTE, 0, TreeSharedMemComm);

  int shmrank = Shmem.GetShmRankForSimulCommRank[Shmem.Sim_ThisTask];
  MPI_Bcast(&shmrank, 1, MPI_INT, 0, TreeSharedMemComm);

  NodeLevel   = (unsigned char *)((char *)Shmem.SharedMemBaseAddr[shmrank] + off[0]);
  NodeSibling = (int *)((char *)Shmem.SharedMemBaseAddr[shmrank] + off[1]);
  NodeIndex   = (int *)((char *)Shmem.SharedMemBaseAddr[shmrank] + off[2]);
  TopNodes    = (node *)((char *)Shmem.SharedMemBaseAddr[shmrank] + off[3]);
}

template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::prepare_shared_memory_access(void)
{
  if(D->NTask != Shmem.Sim_NTask)
    Terminate("This version of the tree communication algorithm only works with the full simulation partition");

  ptrdiff_t off;
  off = (char *)Nodes - Mem.Base;
  MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, TreeNodes_offsets, sizeof(ptrdiff_t), MPI_BYTE, TreeSharedMemComm);
  off = (char *)Points - Mem.Base;
  MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, TreePoints_offsets, sizeof(ptrdiff_t), MPI_BYTE, TreeSharedMemComm);
  off = (char *)Nextnode - Mem.Base;
  MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, TreeNextnode_offsets, sizeof(ptrdiff_t), MPI_BYTE, TreeSharedMemComm);
  off = (char *)Foreign_Nodes - Mem.Base;
  MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, TreeForeign_Nodes_offsets, sizeof(ptrdiff_t), MPI_BYTE, TreeSharedMemComm);
  off = (char *)Foreign_Points - Mem.Base;
  MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, TreeForeign_Points_offsets, sizeof(ptrdiff_t), MPI_BYTE, TreeSharedMemComm);
  off = (char *)Tp->P - Mem.Base;
  MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, TreeP_offsets, sizeof(ptrdiff_t), MPI_BYTE, TreeSharedMemComm);
  off = (char *)Tp->SphP - Mem.Base;
  MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, TreeSphP_offsets, sizeof(ptrdiff_t), MPI_BYTE, TreeSharedMemComm);

  if(Shmem.Sim_NTask < Shmem.World_NTask)
    {
      Shmem.tree_info[TreeInfoHandle].Bd.MaxPart            = MaxPart;
      Shmem.tree_info[TreeInfoHandle].Bd.MaxNodes           = MaxNodes;
      Shmem.tree_info[TreeInfoHandle].Bd.NTopnodes          = D->NTopnodes;
      Shmem.tree_info[TreeInfoHandle].Bd.ImportedNodeOffset = ImportedNodeOffset;
      Shmem.tree_info[TreeInfoHandle].Bd.EndOfTreePoints    = EndOfTreePoints;
      Shmem.tree_info[TreeInfoHandle].Bd.EndOfForeignNodes  = EndOfForeignNodes;

      // need to inform also our shared shared memory processor
      if(TreeSharedMem_ThisTask == 0)
        {
          MPI_Send(&TreeInfoHandle, 1, MPI_INT, Shmem.MyShmRankInGlobal, TAG_METDATA, MPI_COMM_WORLD);
          MPI_Send(&All.Ti_Current, sizeof(All.Ti_Current), MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_METDATA + 1, MPI_COMM_WORLD);
          MPI_Send(&Shmem.tree_info[TreeInfoHandle].Bd, sizeof(Shmem.tree_info[TreeInfoHandle].Bd), MPI_BYTE, Shmem.MyShmRankInGlobal,
                   TAG_METDATA + 2, MPI_COMM_WORLD);

          intposconvert *convfac = Tp;
          MPI_Send(convfac, sizeof(intposconvert), MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_METDATA + 3, MPI_COMM_WORLD);
        }

      Shmem.inform_offset_table(TopNodes);
      Shmem.inform_offset_table(Nodes);
      Shmem.inform_offset_table(Nextnode);
      Shmem.inform_offset_table(Points);
      Shmem.inform_offset_table(Tp->P);
      Shmem.inform_offset_table(Tp->SphP);
      Shmem.inform_offset_table(Foreign_Nodes);
      Shmem.inform_offset_table(Foreign_Points);

      MPI_Barrier(Shmem.SharedMemComm);  // this barrier is in principle superfluous, but on some systems,
                                         // the MPI_Gather in prepare_offset_table() can return prematurely
                                         // on the target rank before all data has arrived

      /* the following is needed to make sure that the shared memory handler on different nodes is already properly initialized */
      MPI_Barrier(D->Communicator);
    }
}

template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::cleanup_shared_memory_access(void)
{
  if(Shmem.Sim_NTask < Shmem.World_NTask)
    {
      if(TreeSharedMem_ThisTask == 0)
        {
          // need to send this flag to the correct processor rank (our shared memoin the global communicator
          MPI_Send(&TreeInfoHandle, 1, MPI_INT, Shmem.MyShmRankInGlobal, TAG_HEADER, MPI_COMM_WORLD);
        }
    }
}

template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::tree_fetch_foreign_nodes(enum ftype fetch_type)
{
  // find out how many we have to fetch from each node
  int *CountFetch  = (int *)Mem.mymalloc_movable_clear(&CountFetch, "CountFetch", Shmem.World_NTask * sizeof(int));
  int *OffsetFetch = (int *)Mem.mymalloc_movable(&OffsetFetch, "OffsetFetch", Shmem.World_NTask * sizeof(int));

  for(int i = 0; i < NumOnFetchStack; i++)
    CountFetch[StackToFetch[i].GhostRank]++;

  OffsetFetch[0] = 0;
  for(int i = 1; i < Shmem.World_NTask; i++)
    OffsetFetch[i] = OffsetFetch[i - 1] + CountFetch[i - 1];

  mycxxsort(StackToFetch, StackToFetch + NumOnFetchStack, compare_ghostrank);

  /* now go through each node in turn, and import from them the requested nodes
   */
  int CommPTask;
  for(CommPTask = 0; Shmem.World_NTask > (1 << CommPTask); CommPTask++)
    ;

  for(int ngrp = 1; ngrp < (1 << CommPTask); ngrp++)
    {
      int ghost_rank = Shmem.World_ThisTask ^ ngrp;  // this is already the responsible ghost rank

      if(ghost_rank < Shmem.World_NTask)
        {
          if(CountFetch[ghost_rank] > 0)
            {
              node_req *node_req_send =
                  (node_req *)Mem.mymalloc_movable(&node_req_send, "node_req_send", CountFetch[ghost_rank] * sizeof(node_req));

              for(int i = 0; i < CountFetch[ghost_rank]; i++)
                {
                  int k = OffsetFetch[ghost_rank] + i;

                  node *nop = get_nodep(StackToFetch[k].NodeToOpen, StackToFetch[k].ShmRank);

                  node_req_send[i].foreigntask = nop->OriginTask;
                  node_req_send[i].foreignnode = nop->OriginNode;
                }

              /* now make the node request */
              int tag;
              if(fetch_type == FETCH_GRAVTREE)
                tag = TAG_FETCH_GRAVTREE + TreeInfoHandle;
              else if(fetch_type == FETCH_SPH_DENSITY)
                tag = TAG_FETCH_SPH_DENSITY + TreeInfoHandle;
              else if(fetch_type == FETCH_SPH_HYDRO)
                tag = TAG_FETCH_SPH_HYDRO + TreeInfoHandle;
              else if(fetch_type == FETCH_SPH_TREETIMESTEP)
                tag = TAG_FETCH_SPH_TREETIMESTEP + TreeInfoHandle;
              else
                {
                  tag = 0;
                  Terminate("tag undefined");
                }

              MPI_Send(node_req_send, CountFetch[ghost_rank] * sizeof(node_req), MPI_BYTE, ghost_rank, tag, MPI_COMM_WORLD);
              Mem.myfree(node_req_send);

              // get the information about how many nodes and particles hang below each of the nodes
              node_count_info *node_info_send = (node_count_info *)Mem.mymalloc_movable(
                  &node_info_send, "node_info_send", CountFetch[ghost_rank] * sizeof(node_count_info));

              MPI_Recv(node_info_send, CountFetch[ghost_rank] * sizeof(node_count_info), MPI_BYTE, ghost_rank, TAG_N, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);

              /* now find out how many nodes and points we want to import in total */
              int n_sendpoints = 0;
              int n_sendnodes  = 0;

              for(int i = 0; i < CountFetch[ghost_rank]; i++)
                {
                  n_sendpoints += node_info_send[i].count_parts;
                  n_sendnodes += node_info_send[i].count_nodes;
                }

              foreign_point_data *buf_foreignpoints =
                  (foreign_point_data *)Mem.mymalloc("buf_foreignpoints", n_sendpoints * sizeof(foreign_point_data));

              node *buf_foreignnodes = (node *)Mem.mymalloc("buf_foreignnodes", n_sendnodes * sizeof(node));

              /* now receive the points and nodes */
              if(n_sendpoints > 0)
                MPI_Recv(buf_foreignpoints, n_sendpoints * sizeof(foreign_point_data), MPI_BYTE, ghost_rank, TAG_PDATA, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

              if(n_sendnodes > 0)
                MPI_Recv(buf_foreignnodes, n_sendnodes * sizeof(node), MPI_BYTE, ghost_rank, TAG_SPHDATA, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

              /* now we have to link the nodes and particles into the tree */

              int n_nodes_used = 0;
              int n_parts_used = 0;

              for(int i = 0; i < CountFetch[ghost_rank]; i++)
                {
                  int k = OffsetFetch[ghost_rank] + i;

                  int no      = StackToFetch[k].NodeToOpen;
                  int shmrank = StackToFetch[k].ShmRank;

                  node *nop = get_nodep(no, shmrank);

                  int n_nodes = node_info_send[i].count_nodes;
                  int n_parts = node_info_send[i].count_parts;

                  while(nop->access.test_and_set(std::memory_order_acquire))
                    {
                      // acquire spin lock
                    }

                  if(nop->cannot_be_opened_locally)  // make sure the node hasn't been inserted by another task yet
                    {
                      int pfirst = -1;
                      int plast  = -1;

                      for(int j = 0; j < n_nodes; j++)
                        {
                          if(NumForeignNodes >= MaxForeignNodes)
                            Terminate(
                                "We are out of storage for foreign nodes: NumForeignNodes=%d MaxForeignNodes=%d  j=%d n_parts=%d",
                                NumForeignNodes, MaxForeignNodes, j, n_parts);

                          /* tree index of this node */
                          int p = EndOfTreePoints + NumForeignNodes++;

                          //  cannot do a Foreign_Nodes[p - EndOfTreePoints] = buf_foreignnodes[n_nodes_used++]; because out
                          //  std::atomic_flag has a deleted default copy operator
                          memcpy(static_cast<void *>(&Foreign_Nodes[p - EndOfTreePoints]),
                                 static_cast<void *>(&buf_foreignnodes[n_nodes_used++]), sizeof(node));

                          Foreign_Nodes[p - EndOfTreePoints].access.clear();
                          Foreign_Nodes[p - EndOfTreePoints].flag_already_fetched = 0;

                          Foreign_Nodes[p - EndOfTreePoints].nextnode = -1;

                          if(plast >= 0) /* this is here still the previous one */
                            {
                              if(plast < EndOfForeignNodes) /* have a foreign node */
                                {
                                  Foreign_Nodes[plast - EndOfTreePoints].sibling         = p;
                                  Foreign_Nodes[plast - EndOfTreePoints].sibling_shmrank = Shmem.Island_ThisTask;
                                }
                              else
                                {
                                  Foreign_Points[plast - EndOfForeignNodes].Nextnode         = p;
                                  Foreign_Points[plast - EndOfForeignNodes].Nextnode_shmrank = Shmem.Island_ThisTask;
                                }
                            }

                          if(pfirst < 0)
                            pfirst = p;

                          plast = p;
                        }

                      for(int j = 0; j < n_parts; j++)
                        {
                          if(NumForeignPoints >= MaxForeignPoints)
                            Terminate(
                                "We are out of storage for foreign points: NumForeignPoints=%d MaxForeignPoints=%d  j=%d n_parts=%d",
                                NumForeignPoints, MaxForeignPoints, j, n_parts);

                          /* global tree index of this foreign point */
                          int p                                 = EndOfForeignNodes + NumForeignPoints++;
                          Foreign_Points[p - EndOfForeignNodes] = buf_foreignpoints[n_parts_used++];

                          if(plast >= 0) /* this is here still the previous one */
                            {
                              if(plast < EndOfForeignNodes) /* have a foreign node */
                                {
                                  Foreign_Nodes[plast - EndOfTreePoints].sibling         = p;
                                  Foreign_Nodes[plast - EndOfTreePoints].sibling_shmrank = Shmem.Island_ThisTask;
                                }
                              else
                                {
                                  Foreign_Points[plast - EndOfForeignNodes].Nextnode         = p;
                                  Foreign_Points[plast - EndOfForeignNodes].Nextnode_shmrank = Shmem.Island_ThisTask;
                                }
                            }

                          if(pfirst < 0)
                            pfirst = p;

                          plast = p;
                        }

                      if(plast < 0 || pfirst < 0)
                        Terminate("plast=%d < 0 || pfirst=%d < 0   n_nodes=%d n_parts=%d", plast, pfirst, n_nodes, n_parts);

                      if(plast >= 0) /* this is here still the previous one */
                        {
                          if(plast < EndOfForeignNodes) /* have a foreign node */
                            {
                              Foreign_Nodes[plast - EndOfTreePoints].sibling         = nop->sibling;
                              Foreign_Nodes[plast - EndOfTreePoints].sibling_shmrank = nop->sibling_shmrank;

                              if(Foreign_Nodes[plast - EndOfTreePoints].cannot_be_opened_locally == 0)
                                if(D->ThisTask == 0)
                                  Terminate("what?    plast - EndOfTreePoints=%d", plast - EndOfTreePoints);
                            }
                          else
                            {
                              Foreign_Points[plast - EndOfForeignNodes].Nextnode         = nop->sibling;
                              Foreign_Points[plast - EndOfForeignNodes].Nextnode_shmrank = nop->sibling_shmrank;
                            }
                        }

                      nop->nextnode         = pfirst;
                      nop->nextnode_shmrank = Shmem.Island_ThisTask;

                      nop->cannot_be_opened_locally.store(0, std::memory_order_release);

                      sum_NumForeignNodes += n_nodes;
                      sum_NumForeignPoints += n_parts;
                    }
                  else
                    {
                      // skip this node, because apparently it was fetched in the meantime by another mpi task
                      n_nodes_used += n_nodes;
                      n_parts_used += n_parts;
                    }

                  nop->access.clear(std::memory_order_release);
                }

              if(n_sendpoints != n_parts_used || n_sendnodes != n_nodes_used)
                Terminate("n_sendpoints != n_parts_used || n_sendnodes != n_nodes_used");

              Mem.myfree(buf_foreignnodes);
              Mem.myfree(buf_foreignpoints);
              Mem.myfree(node_info_send);
            }
        }
    }

  NumOnFetchStack = 0;

  Mem.myfree(OffsetFetch);
  Mem.myfree(CountFetch);
}

/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::treefree(void)
{
  MPI_Comm_free(&TreeSharedMemComm);

  if(MaxPart == 0)
    return;  // nothing to be done

  if(Nodes)
    {
      if(Father)
        {
          Mem.myfree_movable(Father);
          Mem.myfree_movable(Nextnode);
        }
      if(Points)
        {
          Mem.myfree_movable(Points);
        }

      Mem.myfree_movable(Nodes + MaxPart + D->NTopnodes);

      if(TreeSharedMem_ThisTask == 0)
        {
          if(TreeSharedMem_NTask == Shmem.World_NTask || D->NTask < Shmem.Sim_NTask)
            {
              Mem.myfree_movable(TopNodes + MaxPart);
              Mem.myfree_movable(NodeIndex);
              Mem.myfree_movable(NodeSibling);
              Mem.myfree_movable(NodeLevel);
            }
          else
            {
              int ghost_rank = Shmem.GetGhostRankForSimulCommRank[Shmem.Sim_ThisTask];

              // tell the ghost rank to free the storage
              MPI_Send(&TreeInfoHandle, 1, MPI_INT, ghost_rank, TAG_TOPNODE_FREE, MPI_COMM_WORLD);
            }
        }

      Mem.myfree(TreeSharedMemBaseAddr);

      Mem.myfree(TreePS_offsets);
      Mem.myfree(TreeSphP_offsets);
      Mem.myfree(TreeP_offsets);
      Mem.myfree(TreeForeign_Points_offsets);
      Mem.myfree(TreeForeign_Nodes_offsets);
      Mem.myfree(TreeNextnode_offsets);
      Mem.myfree(TreePoints_offsets);
      Mem.myfree(TreeNodes_offsets);

      Mem.myfree_movable(Recv_offset);
      Mem.myfree_movable(Recv_count);
      Mem.myfree_movable(Send_offset);
      Mem.myfree_movable(Send_count);

      Nodes       = NULL;
      TopNodes    = NULL;
      NodeIndex   = NULL;
      NodeSibling = NULL;
      NodeLevel   = NULL;
      Points      = NULL;
      Nextnode    = NULL;
      Father      = NULL;
    }
  else
    Terminate("trying to free the tree even though it's not allocated");
}

template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::tree_export_node_threads(int no, int i, thread_data *thread,
                                                                                   offset_tuple off)
{
  int task      = D->TaskOfLeaf[no - (MaxPart + MaxNodes)];
  int nodeindex = NodeIndex[no - (MaxPart + MaxNodes)];

  tree_export_node_threads_by_task_and_node(task, nodeindex, i, thread, off);
}

template <typename node, typename partset, typename point_data, typename foreign_point_data>
void tree<node, partset, point_data, foreign_point_data>::tree_export_node_threads_by_task_and_node(int task, int nodeindex, int i,
                                                                                                    thread_data *thread,
                                                                                                    offset_tuple off)
{
  if(task < 0 || task >= D->NTask)
    Terminate("task < 0 || task >= D->NTask");

  if(task != D->ThisTask)
    {
      if(thread->Exportflag[task] != i)
        {
          thread->Exportflag[task]     = i;
          int nexp                     = thread->Nexport++;
          thread->PartList[nexp].Task  = task;
          thread->PartList[nexp].Index = i;
          thread->ExportSpace -= thread->ItemSize;
        }

      int nexp                     = thread->NexportNodes++;
      nexp                         = -1 - nexp;
      data_nodelist *nodelist      = (data_nodelist *)(((char *)thread->PartList) + thread->InitialSpace);
      nodelist[nexp].Task          = task;
      nodelist[nexp].Index         = i;
      nodelist[nexp].NodeInfo.Node = nodeindex;
      thread->ExportSpace -= (sizeof(data_nodelist) + sizeof(int));
    }
}

#include "../ngbtree/ngbtree.h"
template class tree<ngbnode, simparticles, ngbpoint_data, foreign_sphpoint_data>;

#include "../gravtree/gravtree.h"
template class tree<gravnode, simparticles, gravpoint_data, foreign_gravpoint_data>;

#ifdef FOF
#include "../fof/foftree.h"
template class tree<fofnode, simparticles, fofpoint_data, foreign_fofpoint_data>;
#if defined(LIGHTCONE) && (defined(LIGHTCONE_PARTICLES_GROUPS) || defined(LIGHTCONE_IMAGE_COMP_HSML_VELDISP))
/* make sure that we instantiate the template */
#include "../data/lcparticles.h"
template class tree<fofnode, lcparticles, fofpoint_data, foreign_fofpoint_data>;
template class tree<gravnode, lcparticles, gravpoint_data, foreign_gravpoint_data>;
#endif
#endif
