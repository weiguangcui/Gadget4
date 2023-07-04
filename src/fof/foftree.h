/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file foftree.h
 *
 *  \brief declares structures and classes needed for the FOF neighbor tree
 */

#ifndef FOFTREE_H
#define FOFTREE_H

#include "gadgetconfig.h"

#include "../data/simparticles.h"
#include "../tree/tree.h"

struct fofpoint_data
{
  MyIntPosType IntPos[3];
  int no;
};

struct foreign_fofpoint_data
{
  int Nextnode;
  unsigned char Nextnode_shmrank;
};

struct fofnode : public basenode
{
  MyIntPosType range_min[3];
  MyIntPosType range_max[3];
};

template <typename partset>
class foftree : public tree<fofnode, partset, fofpoint_data, foreign_fofpoint_data>
{
 public:
  typedef tree<fofnode, partset, fofpoint_data, foreign_fofpoint_data> basetree;
  using basetree::Buildtime;
  using basetree::D;
  using basetree::Father;
  using basetree::FirstNonTopLevelNode;
  using basetree::get_nextnodep;
  using basetree::get_nodep;
  using basetree::get_Pp;
  using basetree::get_PSp;
  using basetree::ImportedNodeOffset;
  using basetree::MaxNodes;
  using basetree::MaxPart;
  using basetree::Nextnode;
  using basetree::NodeIndex;
  using basetree::Nodes;
  using basetree::NumNodes;
  using basetree::Recv_count;
  using basetree::Recv_offset;
  using basetree::Send_count;
  using basetree::Send_offset;
  using basetree::TopNodes;
  using basetree::Tp;  // This makes sure that we can access Tp from the base class without having to use "this"
  using basetree::tree_export_node_threads;
  using basetree::tree_export_node_threads_by_task_and_node;
  using basetree::TreeP_offsets;
  using basetree::TreePS_offsets;
  using basetree::TreeSharedMem_ThisTask;
  using basetree::TreeSharedMemBaseAddr;
  using basetree::TreeSphP_offsets;

  int *FullyLinkedNodePIndex;  // If a tree node is fully linked, this gives one particle of the encompassing group

  void update_node_recursive(int no, int sib, int mode) override;
  void exchange_topleafdata(void) override;
  void fill_in_export_points(fofpoint_data *exp_point, int i, int no) override;
  void report_log_message(void) override;

  int treefind_fof_primary(MyIntPosType *searchcenter, MyNgbTreeFloat hsml, int target, int mode, thread_data *thread, int numnodes,
                           node_info *firstnode, int *ngblist, MyIDStorage target_MinID, int target_MinIDTask,
                           double target_DistanceOrigin);
  int treefind_fof_return_a_particle_in_cell_recursive(int no);

  void fof_link_particles_in_cell_recursive(int no, int q);

  int treefind_fof_check_single_node_for_full_linking(int no);
};

#endif
