/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  shared_mem_handler.h
 *
 *  \brief provides a class for accessing data of other MPI ranks via shared memory and designated MPI handler ranks
 */

#ifndef SHAREDMEM_H
#define SHAREDMEM_H

#include "gadgetconfig.h"

#include <hdf5.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <atomic>
#include <cstring>

#include "../data/simparticles.h"

#define MAX_TREE_INFOS 10

class shmem
{
 public:
  MPI_Comm SharedMemComm;   // the communicator linking the processes that have mutual shared memory access in the same node
  MPI_Comm SimulationComm;  // the communicator containing all the compute processors (or all the ghost processors)

  int World_ThisTask;  // rank
  int World_NTask;     // total number of MPI processes

  int Island_ThisTask;  // rank in current shared memory region
  int Island_NTask;     // number of MPI tasks in shared memory region

  int Sim_ThisTask;  // rank in simulation partition
  int Sim_NTask;     // size of MPI tasks in simulation partition

  int GhostRank;  // equal to 1 if we are a ghost rank, otherwise zero

  int Island_Smallest_WorldTask;  // this is the smallest global rank in the shared memory node

  // we need a table that maps the rank of a destination processor in the
  // simulation communicator to the rank of the responsible ghost processor in
  // the global communicator
  int *GetGhostRankForSimulCommRank;

  // we need a table that maps the rank of a destination processor in the
  // simulation communicator to the rank in the shared memory communicator
  int *GetShmRankForSimulCommRank;

  // we need a table that maps the rank of a simulation processor to the
  // smallest world rank in its shared memory node. With this we can decide
  // whether two ranks are on the same node
  int *GetNodeIDForSimulCommRank;

  // the rank in the global communicator that a processor should turn to for a shared memory request
  int MyShmRankInGlobal;

  MPI_Win SharedMemWin;

  void **SharedMemBaseAddr;

#ifdef ALLOCATE_SHARED_MEMORY_VIA_POSIX
  char **SharedMemBaseAddrRaw;
#endif

  char *TableData;
  char *EwaldData;

  struct bookkeeping_data
  {
    int MaxPart;
    int MaxNodes;
    int NTopnodes;
    int ImportedNodeOffset;
    int EndOfTreePoints;
    int EndOfForeignNodes;
  };

  struct tree_storage_info
  {
    bookkeeping_data Bd;

    ptrdiff_t *TopNodes_offsets;
    ptrdiff_t *Nodes_offsets;
    ptrdiff_t *Nextnode_offsets;
    ptrdiff_t *Points_offsets;
    ptrdiff_t *P_offsets;
    ptrdiff_t *SphP_offsets;
    ptrdiff_t *Foreign_Points_offsets;
    ptrdiff_t *Foreign_Nodes_offsets;

    char *TopNodes_storage;
    char *NodeLevel_storage;
    char *NodeSibling_storage;
    char *NodeIndex_storage;
  };

  tree_storage_info tree_info[MAX_TREE_INFOS];
  int num_tree_info = 0;

  inline char *get_basenodep(int no, unsigned char shmrank, int handle)
  {
    if(no < tree_info[handle].Bd.MaxPart + tree_info[handle].Bd.NTopnodes)
      return ((char *)SharedMemBaseAddr[shmrank] + tree_info[handle].TopNodes_offsets[shmrank]);
    else
      return ((char *)SharedMemBaseAddr[shmrank] + tree_info[handle].Nodes_offsets[shmrank]);
  }

  inline int *get_nextnodep(unsigned char shmrank, int handle)
  {
    return (int *)((char *)SharedMemBaseAddr[shmrank] + tree_info[handle].Nextnode_offsets[shmrank]);
  }

  inline char *get_pointsp(unsigned char shmrank, int handle)
  {
    return ((char *)SharedMemBaseAddr[shmrank] + tree_info[handle].Points_offsets[shmrank]);
  }

  inline char *get_Pp(unsigned char shmrank, int handle)
  {
    return ((char *)SharedMemBaseAddr[shmrank] + tree_info[handle].P_offsets[shmrank]);
  }

  inline char *get_SphPp(unsigned char shmrank, int handle)
  {
    return ((char *)SharedMemBaseAddr[shmrank] + tree_info[handle].SphP_offsets[shmrank]);
  }

  void deal_with_gravity_node_request(char *message, int length, int source, int handle);
  void deal_with_sph_node_request(char *message, int length, int source, int handle, simparticles *Sp);

  void prepare_offset_table(void *p, ptrdiff_t *&offset_tab);
  void inform_offset_table(void *p);
  void free_offset_table(ptrdiff_t *&offset_tab);

  void shared_memory_handler(void);
};

extern shmem Shmem;

#endif
