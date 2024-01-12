/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  shared_mem_handler.cc
 *
 *  \brief implements code for the shared-memory fetching of remote date through  designated MPI handler ranks
 */

#include "gadgetconfig.h"

#include <hdf5.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstring>

#include "../gravtree/gravtree.h"
#include "../ngbtree/ngbtree.h"
#include "../time_integration/driftfac.h"

#include "../mpi_utils/shared_mem_handler.h"

typedef gravtree<simparticles> gtree;
typedef ngbtree ntree;

void shmem::prepare_offset_table(void *p, ptrdiff_t *&offset_tab)  // called by ghost task
{
  ptrdiff_t off = ((char *)p - Mem.Base);

  offset_tab = (ptrdiff_t *)Mem.mymalloc("offset_tab", Island_NTask * sizeof(ptrdiff_t));

  MPI_Gather(&off, sizeof(ptrdiff_t), MPI_BYTE, offset_tab, sizeof(ptrdiff_t), MPI_BYTE, Island_NTask - 1, SharedMemComm);
}

void shmem::inform_offset_table(void *p)  // called by worked tasks
{
  ptrdiff_t off = ((char *)p - Mem.Base);

  MPI_Gather(&off, sizeof(ptrdiff_t), MPI_BYTE, NULL, sizeof(ptrdiff_t), MPI_BYTE, Island_NTask - 1, SharedMemComm);
}

void shmem::free_offset_table(ptrdiff_t *&offset_tab) { Mem.myfree(offset_tab); }

void shmem::shared_memory_handler(void)
{
  simparticles Dp{MPI_COMM_WORLD}; /* dummy needed to access drift functions for ngbtree */

  /* first, we wait for the parameter All.MaxMemSize, so that we can initialize the memory handler */
  MPI_Bcast(All.get_data_ptr(), All.get_data_size(), MPI_BYTE, 0, MPI_COMM_WORLD);
  Mem.mymalloc_init(All.MaxMemSize, RST_BEGIN);

  while(true)
    {
      /* wait for an incoming message */
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      int source = status.MPI_SOURCE;
      int tag    = status.MPI_TAG;

      int length;
      MPI_Get_count(&status, MPI_BYTE, &length);

      /* now pick it up */
      char *message = (char *)Mem.mymalloc("message", length);
      MPI_Recv(message, length, MPI_BYTE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if(tag == TAG_METDATA)  // signals that we are synchronizing addresses and values for tree access
        {
          int handle = *((int *)message);
          Mem.myfree(message);

          MPI_Recv(&All.Ti_Current, sizeof(All.Ti_Current), MPI_BYTE, source, TAG_METDATA + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(&tree_info[handle].Bd, sizeof(bookkeeping_data), MPI_BYTE, source, TAG_METDATA + 2, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);

          intposconvert *convfac = &Dp;
          MPI_Recv(convfac, sizeof(intposconvert), MPI_BYTE, source, TAG_METDATA + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          prepare_offset_table(NULL, tree_info[handle].TopNodes_offsets);
          prepare_offset_table(NULL, tree_info[handle].Nodes_offsets);
          prepare_offset_table(NULL, tree_info[handle].Nextnode_offsets);
          prepare_offset_table(NULL, tree_info[handle].Points_offsets);
          prepare_offset_table(NULL, tree_info[handle].P_offsets);
          prepare_offset_table(NULL, tree_info[handle].SphP_offsets);
          prepare_offset_table(NULL, tree_info[handle].Foreign_Nodes_offsets);
          prepare_offset_table(NULL, tree_info[handle].Foreign_Points_offsets);

          MPI_Barrier(SharedMemComm);  // this barrier is in principle superfluous, but on some systems,
                                       // the MPI_Gather in prepare_offset_table() can return prematurely before all data has arrived
        }
      else if(tag == TAG_HEADER)  // signals that we are freeing addresses we stored for tree access
        {
          int handle = *((int *)message);
          Mem.myfree(message);

          free_offset_table(tree_info[handle].Foreign_Points_offsets);
          free_offset_table(tree_info[handle].Foreign_Nodes_offsets);
          free_offset_table(tree_info[handle].SphP_offsets);
          free_offset_table(tree_info[handle].P_offsets);
          free_offset_table(tree_info[handle].Points_offsets);
          free_offset_table(tree_info[handle].Nextnode_offsets);
          free_offset_table(tree_info[handle].Nodes_offsets);
          free_offset_table(tree_info[handle].TopNodes_offsets);
        }
      else if(tag >= TAG_FETCH_SPH_DENSITY && tag < TAG_FETCH_SPH_DENSITY + MAX_TREE_INFOS)  // fetch from SPH density tree
        {
          int handle = tag - TAG_FETCH_SPH_DENSITY;
          deal_with_sph_node_request(message, length, source, handle, &Dp);
          Mem.myfree(message);
        }
      else if(tag >= TAG_FETCH_SPH_HYDRO && tag < TAG_FETCH_SPH_HYDRO + MAX_TREE_INFOS)  // fetch from SPH hydro tree
        {
          int handle = tag - TAG_FETCH_SPH_HYDRO;
          deal_with_sph_node_request(message, length, source, handle, &Dp);
          Mem.myfree(message);
        }
      else if(tag >= TAG_FETCH_SPH_TREETIMESTEP && tag < TAG_FETCH_SPH_TREETIMESTEP + MAX_TREE_INFOS)  // fetch from SPH timesteptree
        {
          int handle = tag - TAG_FETCH_SPH_TREETIMESTEP;
          deal_with_sph_node_request(message, length, source, handle, &Dp);
          Mem.myfree(message);
        }
      else if(tag >= TAG_FETCH_GRAVTREE && tag < TAG_FETCH_GRAVTREE + MAX_TREE_INFOS)  // fetch from gravity tree
        {
          int handle = tag - TAG_FETCH_GRAVTREE;
          deal_with_gravity_node_request(message, length, source, handle);
          Mem.myfree(message);
        }
      else if(tag == TAG_KEY)  // request to terminate gracefully
        {
          H5Eset_auto(H5E_DEFAULT, NULL, NULL);
          MPI_Finalize();
          exit(0);
        }
      else if(tag == TAG_TABLE_ALLOC)  // take over storage for TableData
        {
          size_t tab_len = *((size_t *)message);
          Mem.myfree(message);

          TableData = (char *)Mem.mymalloc("table", tab_len);
          MPI_Recv(TableData, tab_len, MPI_BYTE, source, TAG_DMOM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ptrdiff_t off = ((char *)TableData - Mem.Base);
          MPI_Bcast(&off, sizeof(ptrdiff_t), MPI_BYTE, Island_ThisTask, SharedMemComm);
        }
      else if(tag == TAG_TABLE_FREE)
        {
          Mem.myfree(message);
          Mem.myfree(TableData);
        }
      else if(tag == TAG_EWALD_ALLOC)  // take over storage for EwaldData
        {
          size_t tab_len = *((size_t *)message);
          Mem.myfree(message);

          EwaldData = (char *)Mem.mymalloc("table", tab_len);
          MPI_Recv(EwaldData, tab_len, MPI_BYTE, source, TAG_DMOM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ptrdiff_t off = ((char *)EwaldData - Mem.Base);
          MPI_Bcast(&off, sizeof(ptrdiff_t), MPI_BYTE, Island_ThisTask, SharedMemComm);
        }
      else if(tag == TAG_TOPNODE_ALLOC)  // take over the storage of common top-level tree data
        {
          int handle = num_tree_info++;
          if(num_tree_info > MAX_TREE_INFOS)
            Terminate("num_tree_info > MAX_TREE_INFOS");

          size_t *sizep   = ((size_t *)message);
          size_t sizes[4] = {sizep[0], sizep[1], sizep[2], sizep[3]};
          Mem.myfree(message);

          tree_info[handle].NodeLevel_storage   = (char *)Mem.mymalloc("NodeLevel_storage", sizes[0]);
          tree_info[handle].NodeSibling_storage = (char *)Mem.mymalloc("NodeSibling_storage", sizes[1]);
          tree_info[handle].NodeIndex_storage   = (char *)Mem.mymalloc("NodeIndex_storage", sizes[2]);
          tree_info[handle].TopNodes_storage    = (char *)Mem.mymalloc("TopNodes_storage", sizes[3]);

          ptrdiff_t off[4] = {
              ((char *)tree_info[handle].NodeLevel_storage - Mem.Base), ((char *)tree_info[handle].NodeSibling_storage - Mem.Base),
              ((char *)tree_info[handle].NodeIndex_storage - Mem.Base), ((char *)tree_info[handle].TopNodes_storage - Mem.Base)};

          MPI_Send(off, 4 * sizeof(ptrdiff_t), MPI_BYTE, source, TAG_TOPNODE_OFFSET, MPI_COMM_WORLD);

          MPI_Send(&handle, 1, MPI_INT, source, TAG_N, MPI_COMM_WORLD);
        }
      else if(tag == TAG_TOPNODE_FREE)  // free the top-level storage for a tree again
        {
          int handle = *((int *)message);
          Mem.myfree(message);

          num_tree_info--;
          if(handle != num_tree_info)
            Terminate("unexpected handle");

          Mem.myfree(tree_info[handle].TopNodes_storage);
          Mem.myfree(tree_info[handle].NodeIndex_storage);
          Mem.myfree(tree_info[handle].NodeSibling_storage);
          Mem.myfree(tree_info[handle].NodeLevel_storage);
        }
      else if(tag == TAG_DRIFT_INIT)  // make the shared memory handler update All and init the local drift tables
        {
          memcpy(All.get_data_ptr(), message, All.get_data_size());
          Mem.myfree(message);
          Driftfac.init_drift_table();
        }
      else if(tag == TAG_ALL_UPDATE)  // make the shared memory handler update the contents of the All structure
        {
          memcpy(All.get_data_ptr(), message, All.get_data_size());
          Mem.myfree(message);
        }
    }
}

void shmem::deal_with_sph_node_request(char *message, int length, int source, int handle, simparticles *Sp)
{
#ifndef LEAN
  bookkeeping_data &Bdat = tree_info[handle].Bd;

  // we got the list of requested nodes
  ntree::node_req *node_req_recv = (ntree::node_req *)message;
  int nrecv                      = length / sizeof(ntree::node_req);

  /* as part of this message, we get the real actually targeted rank in
   * the simulation communicator. We will translate this to the rank in our shared memory
   * block
   */

  /* now prepare answer message by reading from shared memory */
  /******************* prepare tree answer ********************/

  ntree::node_count_info *node_info_recv =
      (ntree::node_count_info *)Mem.mymalloc("node_info_recv", nrecv * sizeof(ntree::node_count_info));

  /* first let's count how many nodes and particles are hanging below in each case */
  int n_recvpoints = 0;
  int n_recvnodes  = 0;

  for(int i = 0; i < nrecv; i++)
    {
      node_info_recv[i].count_nodes = 0;
      node_info_recv[i].count_parts = 0;

      int no      = node_req_recv[i].foreignnode;
      int task    = node_req_recv[i].foreigntask;
      int shmrank = GetShmRankForSimulCommRank[task];  // corresponding local shared memory rank, stays fixed below

      if(no < Bdat.MaxPart || no >= Bdat.MaxPart + Bdat.MaxNodes)
        Terminate("not an internal node");

      ngbnode *nop = ((ngbnode *)get_basenodep(no, shmrank, handle)) + no;

      int p = nop->nextnode;

      while(p != nop->sibling)
        {
          if(p < 0)
            Terminate("p=%d < 0", p);

          if(p < Bdat.MaxPart) /* a local particle */
            {
              node_info_recv[i].count_parts++;
              p = get_nextnodep(shmrank, handle)[p];
            }
          else if(p < Bdat.MaxPart + Bdat.MaxNodes) /* an internal node  */
            {
              node_info_recv[i].count_nodes++;
              p = (((ngbnode *)get_basenodep(p, shmrank, handle)) + p)->sibling;
            }
          else if(p >= Bdat.ImportedNodeOffset && p < Bdat.EndOfTreePoints) /* an imported tree point */
            {
              node_info_recv[i].count_parts++;
              p = get_nextnodep(shmrank, handle)[p - Bdat.MaxNodes];
            }
          else
            Terminate("p=%d MaxPart=%d MaxNodes=%d", p, Bdat.MaxPart, Bdat.MaxNodes);
        }

      if(node_info_recv[i].count_parts == 0 && node_info_recv[i].count_nodes == 0)
        Terminate("strange: we have we requested an empty node?\n");

      n_recvpoints += node_info_recv[i].count_parts;
      n_recvnodes += node_info_recv[i].count_nodes;
    }

  foreign_sphpoint_data *exportbuf_points = (foreign_sphpoint_data *)Mem.mymalloc_movable(
      &exportbuf_points, "exportbuf_points", n_recvpoints * sizeof(foreign_sphpoint_data));
  ngbnode *exportbuf_nodes = (ngbnode *)Mem.mymalloc_movable(&exportbuf_nodes, "exportbuf_nodes", n_recvnodes * sizeof(ngbnode));

  n_recvpoints = 0;
  n_recvnodes  = 0;

  for(int i = 0; i < nrecv; i++)
    {
      int no      = node_req_recv[i].foreignnode;
      int task    = node_req_recv[i].foreigntask;
      int shmrank = GetShmRankForSimulCommRank[task];  // corresponding local shared memory rank, stays fixed below

      ngbnode *nop = ((ngbnode *)get_basenodep(no, shmrank, handle)) + no;

      int p = nop->nextnode;

      while(p != nop->sibling)
        {
          if(p < Bdat.MaxPart) /* a local particle */
            {
              int off = n_recvpoints++;

              foreign_sphpoint_data *expoints = &exportbuf_points[off];

              particle_data *ptr         = (particle_data *)get_Pp(shmrank, handle) + p;
              sph_particle_data *sph_ptr = (sph_particle_data *)get_SphPp(shmrank, handle) + p;

              particle_data p_copy;
              sph_particle_data sphp_copy;

              if(ptr->get_Ti_Current() != All.Ti_Current)
                {
                  /* Because of possible lightcone output, shared memory fetch not allowed to drift original,
                   * because this rank doesn't have a lightcone buffer. The original node needs to do it.
                   * We thus drift a copy of the particle without allowing lightcone access
                   */
#ifndef LEAN
                  while(ptr->access.test_and_set(std::memory_order_acquire))
                    ;  // acquire spin lock
#endif
                  p_copy    = *ptr;
                  sphp_copy = *sph_ptr;

#ifndef LEAN
                  ptr->access.clear(std::memory_order_release);  // release spin lock
#endif
                  p_copy.access.clear(std::memory_order_release);  // clear spin lock in copy

                  /* use the copy from now on */

                  ptr     = &p_copy;
                  sph_ptr = &sphp_copy;

                  // the final flag tells the drift to not consider lightcone crossings
                  Sp->drift_particle(ptr, sph_ptr, All.Ti_Current, true);
                }

              expoints->IntPos[0]    = ptr->IntPos[0];
              expoints->IntPos[1]    = ptr->IntPos[1];
              expoints->IntPos[2]    = ptr->IntPos[2];
              expoints->Mass         = ptr->getMass();
              expoints->TimeBinHydro = ptr->TimeBinHydro;
              expoints->SphCore      = *sph_ptr;

              expoints->Nextnode = -1;

              p = get_nextnodep(shmrank, handle)[p];
            }
          else if(p < Bdat.MaxPart + Bdat.MaxNodes) /* an internal node  */
            {
              int off = n_recvnodes++;

              ngbnode *sourcep = (((ngbnode *)get_basenodep(p, shmrank, handle)) + p);
              ngbnode *exnodes = &exportbuf_nodes[off];

              if(sourcep->Ti_Current != All.Ti_Current)
                sourcep->drift_node(All.Ti_Current, Sp);

              *exnodes = *sourcep;

              exnodes->cannot_be_opened_locally = 1;
              exnodes->nextnode                 = -1;
              exnodes->sibling                  = -1;
              exnodes->OriginTask               = task;
              exnodes->OriginNode               = p;

              p = sourcep->sibling;
            }
          else if(p >= Bdat.ImportedNodeOffset) /* an imported treepoint particle  */
            {
              Terminate("not expected here");
            }
        }
    }

  /************************************************************/

  MPI_Send(node_info_recv, nrecv * sizeof(ntree::node_count_info), MPI_BYTE, source, TAG_N, MPI_COMM_WORLD);

  /* now transfer the points and nodes */
  if(n_recvpoints > 0)
    MPI_Send(exportbuf_points, n_recvpoints * sizeof(foreign_sphpoint_data), MPI_BYTE, source, TAG_PDATA, MPI_COMM_WORLD);

  if(n_recvnodes > 0)
    MPI_Send(exportbuf_nodes, n_recvnodes * sizeof(ngbnode), MPI_BYTE, source, TAG_SPHDATA, MPI_COMM_WORLD);

  Mem.myfree(exportbuf_nodes);
  Mem.myfree(exportbuf_points);
  Mem.myfree(node_info_recv);
#endif
}

void shmem::deal_with_gravity_node_request(char *message, int length, int source, int handle)
{
  bookkeeping_data &Bdat = tree_info[handle].Bd;

  // we got the list of requested nodes
  gtree::node_req *node_req_recv = (gtree::node_req *)message;
  int nrecv                      = length / sizeof(gtree::node_req);

  /* as part of this message, we get the real actually targeted rank in
   * the simulation communicator. We will translate this to the rank in our shared memory
   * block
   */

  /* now prepare answer message by reading from shared memory */
  /******************* prepare tree answer ********************/

  gtree::node_count_info *node_info_recv =
      (gtree::node_count_info *)Mem.mymalloc("node_info_recv", nrecv * sizeof(gtree::node_count_info));

  /* first let's count how many nodes and particles are hanging below in each case */
  int n_recvpoints = 0;
  int n_recvnodes  = 0;

  for(int i = 0; i < nrecv; i++)
    {
      node_info_recv[i].count_nodes = 0;
      node_info_recv[i].count_parts = 0;

      int no      = node_req_recv[i].foreignnode;
      int task    = node_req_recv[i].foreigntask;
      int shmrank = GetShmRankForSimulCommRank[task];  // corresponding local shared memory rank, stays fixed below

      if(no < Bdat.MaxPart || no >= Bdat.MaxPart + Bdat.MaxNodes)
        Terminate("not an internal node");

      gravnode *nop = ((gravnode *)get_basenodep(no, shmrank, handle)) + no;

      int p = nop->nextnode;

      while(p != nop->sibling)
        {
          if(p < 0)
            Terminate("p=%d < 0", p);

          if(p < Bdat.MaxPart) /* a local particle */
            {
              node_info_recv[i].count_parts++;
              p = get_nextnodep(shmrank, handle)[p];
            }
          else if(p < Bdat.MaxPart + Bdat.MaxNodes) /* an internal node  */
            {
              node_info_recv[i].count_nodes++;
              p = (((gravnode *)get_basenodep(p, shmrank, handle)) + p)->sibling;
            }
          else if(p >= Bdat.ImportedNodeOffset && p < Bdat.EndOfTreePoints) /* an imported tree point */
            {
              node_info_recv[i].count_parts++;
              p = get_nextnodep(shmrank, handle)[p - Bdat.MaxNodes];
            }
          else
            Terminate("p=%d MaxPart=%d MaxNodes=%d", p, Bdat.MaxPart, Bdat.MaxNodes);
        }

      if(node_info_recv[i].count_parts == 0 && node_info_recv[i].count_nodes == 0)
        Terminate("strange: we have we requested an empty node?\n");

      n_recvpoints += node_info_recv[i].count_parts;
      n_recvnodes += node_info_recv[i].count_nodes;
    }

  foreign_gravpoint_data *exportbuf_points = (foreign_gravpoint_data *)Mem.mymalloc_movable(
      &exportbuf_points, "exportbuf_points", n_recvpoints * sizeof(foreign_gravpoint_data));
  gravnode *exportbuf_nodes = (gravnode *)Mem.mymalloc_movable(&exportbuf_nodes, "exportbuf_nodes", n_recvnodes * sizeof(gravnode));

  n_recvpoints = 0;
  n_recvnodes  = 0;

  for(int i = 0; i < nrecv; i++)
    {
      int no      = node_req_recv[i].foreignnode;
      int task    = node_req_recv[i].foreigntask;
      int shmrank = GetShmRankForSimulCommRank[task];  // corresponding local shared memory rank, stays fixed below

      gravnode *nop = ((gravnode *)get_basenodep(no, shmrank, handle)) + no;

      int p = nop->nextnode;

      while(p != nop->sibling)
        {
          if(p < Bdat.MaxPart) /* a local particle */
            {
              int off = n_recvpoints++;

              foreign_gravpoint_data *expoints = &exportbuf_points[off];

              particle_data *ptr = (particle_data *)get_Pp(shmrank, handle) + p;

              expoints->IntPos[0] = ptr->IntPos[0];
              expoints->IntPos[1] = ptr->IntPos[1];
              expoints->IntPos[2] = ptr->IntPos[2];
              expoints->Mass      = ptr->getMass();
              expoints->Type      = ptr->getType();
              expoints->OldAcc    = ptr->OldAcc;
#if NSOFTCLASSES > 1
              expoints->SofteningClass = ptr->getSofteningClass();
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
              expoints->InsideOutsideFlag = ptr->InsideOutsideFlag;
#endif
              expoints->Nextnode = -1;

              p = get_nextnodep(shmrank, handle)[p];
            }
          else if(p < Bdat.MaxPart + Bdat.MaxNodes) /* an internal node  */
            {
              int off = n_recvnodes++;

              gravnode *exnodes = &exportbuf_nodes[off];
              gravnode *sourcep = (((gravnode *)get_basenodep(p, shmrank, handle)) + p);

              memcpy(static_cast<void *>(exnodes), static_cast<void *>(sourcep),
                     sizeof(gravnode));  //  cannot do a  *exnodes = *sourcep; because out std::atomic_flag
                                         //  has a deleted default copy operator

              exnodes->cannot_be_opened_locally = 1;
              exnodes->nextnode                 = -1;
              exnodes->sibling                  = -1;
              exnodes->OriginTask               = task;
              exnodes->OriginNode               = p;

              p = sourcep->sibling;
            }
          else if(p >= Bdat.ImportedNodeOffset) /* an imported Treepoint particle  */
            {
              int off = n_recvpoints++;

              foreign_gravpoint_data *expoints = &exportbuf_points[off];
              int n                            = p - Bdat.ImportedNodeOffset;
              gravpoint_data *pointsp          = ((gravpoint_data *)get_pointsp(shmrank, handle)) + n;

              expoints->IntPos[0] = pointsp->IntPos[0];
              expoints->IntPos[1] = pointsp->IntPos[1];
              expoints->IntPos[2] = pointsp->IntPos[2];
              expoints->Mass      = pointsp->Mass;
              expoints->Type      = pointsp->Type;
              expoints->OldAcc    = pointsp->OldAcc;
#if NSOFTCLASSES > 1
              expoints->SofteningClass = pointsp->SofteningClass;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
              expoints->InsideOutsideFlag = pointsp->InsideOutsideFlag;
#endif
              expoints->Nextnode = -1;

              p = get_nextnodep(shmrank, handle)[p - Bdat.MaxNodes];
            }
        }
    }

  /************************************************************/

  MPI_Send(node_info_recv, nrecv * sizeof(gtree::node_count_info), MPI_BYTE, source, TAG_N, MPI_COMM_WORLD);

  /* now transfer the points and nodes */
  if(n_recvpoints > 0)
    MPI_Send(exportbuf_points, n_recvpoints * sizeof(foreign_gravpoint_data), MPI_BYTE, source, TAG_PDATA, MPI_COMM_WORLD);

  if(n_recvnodes > 0)
    MPI_Send(exportbuf_nodes, n_recvnodes * sizeof(gravnode), MPI_BYTE, source, TAG_SPHDATA, MPI_COMM_WORLD);

  Mem.myfree(exportbuf_nodes);
  Mem.myfree(exportbuf_points);
  Mem.myfree(node_info_recv);
}
