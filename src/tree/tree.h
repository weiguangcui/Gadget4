/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  tree.h
 *
 *  \brief declaration of the base class for building oct-trees
 */

#ifndef TREE_H
#define TREE_H

#ifndef TREE_NUM_BEFORE_NODESPLIT
#define TREE_NUM_BEFORE_NODESPLIT 3  // daughter nodes are only created if there are more than this number of particles in a node
#endif

#define TREE_MODE_BRANCH 0
#define TREE_MODE_TOPLEVEL 1

#define MAX_TREE_ALLOC_FACTOR 30.0

#define TREE_MAX_ITER 100

#include "gadgetconfig.h"

#include <mpi.h>

#include "../domain/domain.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../sort/peano.h"
#include "../system/system.h"

class sim;

/** The tree node data structure. Nodes points to the actual memory
    allocated for the internal nodes, but is shifted such that
    Nodes[Sp.MaxPart] gives the first allocated node. Note that node
    numbers less than Sp.MaxPart are the leaf nodes that contain a
    single particle, and node numbers >= MaxPart+MaxNodes are "pseudo
    particles" that hang off the toplevel leaf nodes belonging to
    other tasks. These are not represented by this structure. Instead,
    the tree traversal for these are saved in the Nextnode, Prevnode
    and Father arrays, indexed with the node number in the case of
    real particles and by nodenumber-MaxNodes for pseudo
    particles.  */

struct basenode
{
  std::atomic<node_bit_field> flag_already_fetched;

  vector<MyIntPosType> center; /**< geometrical center of node */

  int sibling;
  /** The next node in case the current node needs to be
      opened. Applying nextnode repeatedly results in a pure
      depth-first traversal of the tree. */
  int nextnode;
  /** The parent node of the node. (Is -1 for the root node.) */
  int father;

  int OriginTask; /* MPI rank (in full compute communicator) on which this node and its daughter nodes are natively stored */
  int OriginNode; /* Index of the node on the MPI rank that stores it and its daughter nodes natively */

  unsigned char level; /**< hold the tree level, used to store the side length of node in space efficient way */
  unsigned char sibling_shmrank;
  unsigned char nextnode_shmrank;

  std::atomic_flag access;

  std::atomic<unsigned char> cannot_be_opened_locally;

  // unsigned char cannot_be_opened_locally : 1;
  unsigned char not_empty : 1;
};

struct node_info
{
  int Node;
};

struct data_nodelist
{
  int Task;           /** target process */
  int Index;          /** local index that wants to open this node */
  node_info NodeInfo; /** info about node to be opened on foreign process, as well as perdiodic box offset (needed for Ewald summation
                      algorithm for periodic gravity */
};

template <typename node, typename partset, typename point_data, typename foreign_point_data>
class tree
{
  /* class for oct tree */

 public:
  struct index_data
  {
    int p;
    int subnode;
  };

  domain<partset> *D;
  partset *Tp;

  int *Father;
  int *Nextnode;
  int *NodeSibling;
  int *NodeIndex;

  node *TopNodes;
  node *Nodes;
  node *Foreign_Nodes;
  foreign_point_data *Foreign_Points;

  ptrdiff_t *TreeNodes_offsets;
  ptrdiff_t *TreePoints_offsets;
  ptrdiff_t *TreeNextnode_offsets;
  ptrdiff_t *TreeForeign_Nodes_offsets;
  ptrdiff_t *TreeForeign_Points_offsets;
  ptrdiff_t *TreeP_offsets;
  ptrdiff_t *TreeSphP_offsets;
  ptrdiff_t *TreePS_offsets;

  void **TreeSharedMemBaseAddr;

  unsigned char *NodeLevel;

  point_data *Points;
  int *IndexList;
  int *ResultIndexList;

  int *Send_offset;
  int *Send_count;
  int *Recv_count;
  int *Recv_offset;

  int MaxPart;
  int MaxNodes;
  int NumNodes;
  int NumPartImported;
  int NumPartExported;

  int NumForeignNodes;  // number of imported foreign tree nodes
  int MaxForeignNodes;

  int NumForeignPoints;  // number of imported foreign particles to allow completion of local tree walks
  int MaxForeignPoints;

  // for some statistics about the number of imported nodes and points
  long long sum_NumForeignNodes;
  long long sum_NumForeignPoints;

  int FirstNonTopLevelNode;

  int EndOfTreePoints;
  int EndOfForeignNodes;

  int ImportedNodeOffset;
  int Ninsert;
  int NextFreeNode;

  MPI_Comm TreeSharedMemComm;
  int TreeSharedMem_ThisTask = 0;
  int TreeSharedMem_NTask    = 0;

  int TreeInfoHandle;

  double Buildtime;

  int NumOnFetchStack;
  int MaxOnFetchStack;

  struct fetch_data
  {
    int NodeToOpen;
    int ShmRank;
    int GhostRank;
  };
  fetch_data *StackToFetch;

  static bool compare_ghostrank(const fetch_data &a, const fetch_data &b) { return a.GhostRank < b.GhostRank; }

  int NumOnWorkStack;
  int MaxOnWorkStack;
  int NewOnWorkStack;
  int AllocWorkStackBaseLow;
  int AllocWorkStackBaseHigh;

  struct workstack_data
  {
    int Target;
    int Node;
    int ShmRank;
    int MinTopLeafNode;
  };
  workstack_data *WorkStack;

  static bool compare_workstack(const workstack_data &a, const workstack_data &b)
  {
    if(a.MinTopLeafNode < b.MinTopLeafNode)
      return true;
    if(a.MinTopLeafNode > b.MinTopLeafNode)
      return false;

    return a.Target < b.Target;
  }

  void tree_add_to_fetch_stack(node *nop, int nodetoopen, unsigned char shmrank)
  {
    if(NumOnFetchStack >= MaxOnFetchStack)
      {
        Terminate("we shouldn't get here");
        MaxOnFetchStack *= 1.1;
        StackToFetch = (fetch_data *)Mem.myrealloc_movable(StackToFetch, MaxOnFetchStack * sizeof(fetch_data));
      }

    node_bit_field mybit = (((node_bit_field)1) << Shmem.Island_ThisTask);

    node_bit_field oldval = nop->flag_already_fetched.fetch_or(mybit);

    if((oldval & mybit) == 0)  // it wasn't fetched by me yet
      {
        int ghostrank = Shmem.GetGhostRankForSimulCommRank[nop->OriginTask];

        StackToFetch[NumOnFetchStack].NodeToOpen = nodetoopen;
        StackToFetch[NumOnFetchStack].ShmRank    = shmrank;
        StackToFetch[NumOnFetchStack].GhostRank  = ghostrank;

        NumOnFetchStack++;
      }
  }

  void tree_add_to_work_stack(int target, int no, unsigned char shmrank, int mintopleafnode)
  {
    if(NumOnWorkStack + NewOnWorkStack >= MaxOnWorkStack)
      {
        Terminate("we shouldn't get here");
        MaxOnWorkStack *= 1.1;
        WorkStack = (workstack_data *)Mem.myrealloc_movable(WorkStack, MaxOnWorkStack * sizeof(workstack_data));
      }

    WorkStack[NumOnWorkStack + NewOnWorkStack].Target         = target;
    WorkStack[NumOnWorkStack + NewOnWorkStack].Node           = no;
    WorkStack[NumOnWorkStack + NewOnWorkStack].ShmRank        = shmrank;
    WorkStack[NumOnWorkStack + NewOnWorkStack].MinTopLeafNode = mintopleafnode;

    NewOnWorkStack++;
  }

  struct node_count_info
  {
    int count_nodes;
    int count_parts;
  };

  struct node_req
  {
    int foreigntask;
    int foreignnode;
  };

  void prepare_shared_memory_access(void);
  void cleanup_shared_memory_access(void);

  enum ftype
  {
    FETCH_GRAVTREE,
    FETCH_SPH_DENSITY,
    FETCH_SPH_HYDRO,
    FETCH_SPH_TREETIMESTEP,
  };

  void tree_fetch_foreign_nodes(enum ftype fetch_type);
  void tree_initialize_leaf_node_access_info(void);

  inline foreign_point_data *get_foreignpointsp(int n, unsigned char shmrank)
  {
    return (foreign_point_data *)((char *)TreeSharedMemBaseAddr[shmrank] + TreeForeign_Points_offsets[shmrank]) + n;
  }

  inline subfind_data *get_PSp(int n, unsigned char shmrank)
  {
    return (subfind_data *)((char *)TreeSharedMemBaseAddr[shmrank] + TreePS_offsets[shmrank]) + n;
  }

  typedef decltype(Tp->P) pdata;

  inline pdata get_Pp(int n, unsigned char shmrank)
  {
    return (pdata)((char *)TreeSharedMemBaseAddr[shmrank] + TreeP_offsets[shmrank]) + n;
  }

  inline sph_particle_data *get_SphPp(int n, unsigned char shmrank)
  {
    return (sph_particle_data *)((char *)TreeSharedMemBaseAddr[shmrank] + TreeSphP_offsets[shmrank]) + n;
  }

 private:
  /** Gives next node in tree walk for the "particle" nodes. Entries 0
         -- MaxPart-1 are the real particles, and the "pseudoparticles" are
           indexed by the node number-MaxNodes. */

  /** Gives previous node in tree walk for the leaf (particle)
      nodes. Entries 0 -- MaxPart-1 are the real particles, and the
      "pseudoparticles" are indexed by the node number-MaxNodes. */

 public:
  tree() /* constructor */
  {
    TopNodes    = NULL;
    Nodes       = NULL;
    NodeIndex   = NULL;
    NodeSibling = NULL;
    NodeLevel   = NULL;
    Points      = NULL;
    Nextnode    = NULL;
    Father      = NULL;
    D           = NULL;
  }

  /** public functions */
  int treebuild(int ninsert, int *indexlist);
  void treefree(void);
  void treeallocate(int max_partindex, partset *Pptr, domain<partset> *Dptr);
  void treeallocate_share_topnode_addresses(void);

  void tree_export_node_threads(int no, int i, thread_data *thread, offset_tuple off = 0);
  void tree_export_node_threads_by_task_and_node(int task, int nodeindex, int i, thread_data *thread, offset_tuple off = 0);

  virtual void update_node_recursive(int no, int sib, int mode)            = 0;
  virtual void exchange_topleafdata(void)                                  = 0;
  virtual void report_log_message(void)                                    = 0;
  virtual void fill_in_export_points(point_data *exp_point, int i, int no) = 0;

  inline node *get_nodep(int no)
  {
    node *nop;

    if(no < MaxPart + D->NTopnodes)
      nop = &TopNodes[no];
    else if(no < MaxPart + MaxNodes)
      nop = &Nodes[no];
    else
      Terminate("illegale node index");

    return nop;
  }

  inline node *get_nodep(int no, unsigned char shmrank)
  {
    node *nop;

    if(no < MaxPart + D->NTopnodes)
      nop = &TopNodes[no];
    else if(no < MaxPart + MaxNodes) /* an internal node */
      {
        node *Nodes_shmrank = (node *)((char *)TreeSharedMemBaseAddr[shmrank] + TreeNodes_offsets[shmrank]);
        nop                 = &Nodes_shmrank[no];
      }
    else if(no >= EndOfTreePoints && no < EndOfForeignNodes) /* an imported tree node */
      {
        node *Foreign_Nodes_shmrank = (node *)((char *)TreeSharedMemBaseAddr[shmrank] + TreeForeign_Nodes_offsets[shmrank]);

        nop = &Foreign_Nodes_shmrank[no - EndOfTreePoints];
      }
    else
      Terminate("illegale node index");

    return nop;
  }

  inline int *get_nextnodep(unsigned char shmrank)
  {
    return (int *)((char *)TreeSharedMemBaseAddr[shmrank] + TreeNextnode_offsets[shmrank]);
  }

  inline point_data *get_pointsp(int no, unsigned char shmrank)
  {
    return (point_data *)((char *)TreeSharedMemBaseAddr[shmrank] + TreePoints_offsets[shmrank]) + no;
  }

  inline void tree_get_node_and_task(int i, int &no, int &task)
  {
    MyIntPosType xxb       = Tp->P[i].IntPos[0];
    MyIntPosType yyb       = Tp->P[i].IntPos[1];
    MyIntPosType zzb       = Tp->P[i].IntPos[2];
    MyIntPosType mask      = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - 1));
    unsigned char shiftx   = (BITS_FOR_POSITIONS - 3);
    unsigned char shifty   = (BITS_FOR_POSITIONS - 2);
    unsigned char shiftz   = (BITS_FOR_POSITIONS - 1);
    unsigned char level    = 0;
    unsigned char rotation = 0;

#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
    Tp->P[i].InsideOutsideFlag = Tp->check_high_res_point_location(Tp->P[i].IntPos);
#endif

    no = 0;
    while(D->TopNodes[no].Daughter >= 0)  // walk down top tree to find correct leaf
      {
        unsigned char pix     = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                             ((unsigned char)((zzb & mask) >> (shiftz--))));
        unsigned char subnode = peano_incremental_key(pix, &rotation);

        mask >>= 1;
        level++;

        no = D->TopNodes[no].Daughter + subnode;
      }

    no   = D->TopNodes[no].Leaf;
    task = D->TaskOfLeaf[no];
  }

 private:
  /* private member functions */

  int treebuild_construct(void);
  int treebuild_insert_group_of_points(int num, index_data *index_list, int th, unsigned char level, int sibling);
  int create_empty_nodes(int no, int level, int topnode, int bits, int sibling, MyIntPosType x, MyIntPosType y, MyIntPosType z);

 private:
  /* sort kernel */
  static inline bool compare_index_data_subnode(const index_data &a, const index_data &b) { return a.subnode < b.subnode; }
};

#endif
