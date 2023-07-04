/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  ngbtree.h
 *
 *  \brief declaration of the class providing the neighbor tree routines
 */

#ifndef NGBTREE_H_
#define NGBTREE_H_

#include "gadgetconfig.h"

#include "../data/simparticles.h"
#include "../time_integration/driftfac.h"
#include "../tree/tree.h"

/** The ngb-tree data structure
 */

struct ngbpoint_data
{
  MyIntPosType IntPos[3];
  int no;
};

struct ngbnode : public basenode
{
  // we do this ugly trick of using memcpy for our own copy constructor and assignment operator
  // because the atomic_flag in particle_data has an implicitly deleted copy operator... so that the implicit functions
  // for this are unavailable. But we know what we are doing here, and surrounding this with an ugly hack
  // is the easiest way at the moment to work around this in our case unnecessary protection

  ngbnode() {}

  // declare our own copy constructor
  ngbnode(ngbnode &other) { memcpy(static_cast<void *>(this), static_cast<void *>(&other), sizeof(ngbnode)); }

  // declare our own assignment operator
  ngbnode &operator=(ngbnode &other)
  {
    memcpy(static_cast<void *>(this), static_cast<void *>(&other), sizeof(ngbnode));
    return *this;
  }

  MySignedIntPosType center_offset_min[3];
  MySignedIntPosType center_offset_max[3];

  MyNgbTreeFloat vmin[3];
  MyNgbTreeFloat vmax[3];
  MyNgbTreeFloat MaxHsml;
  MyNgbTreeFloat MaxDtHsml;
  MyNgbTreeFloat MaxCsnd;

  std::atomic<integertime> Ti_Current;

  void drift_node(integertime time1, simparticles *Sp)
  {
    while(access.test_and_set(std::memory_order_acquire))
      ;  // acquire spin lock

    if(Ti_Current != time1)
      {
        double dt_drift;

        if(All.ComovingIntegrationOn)
          dt_drift = Driftfac.get_drift_factor(Ti_Current, time1);
        else
          dt_drift = (time1 - Ti_Current) * All.Timebase_interval;

        /* get the shift in the enclosing box in coordinate space */
        double posdiff_min[3];
        double posdiff_max[3];
        for(int j = 0; j < 3; j++)
          {
            posdiff_min[j] = vmin[j] * dt_drift;
            posdiff_max[j] = vmax[j] * dt_drift;
          }

        /* convert to integer coordinate shifts */
        MySignedIntPosType delta_min[3];
        MySignedIntPosType delta_max[3];
        Sp->pos_to_signedintpos(posdiff_min, delta_min);
        Sp->pos_to_signedintpos(posdiff_max, delta_max);

        /* not adjust the bounding box in integer coordinates */
        for(int j = 0; j < 3; j++)
          {
            center_offset_min[j] += delta_min[j];
            center_offset_max[j] += delta_max[j];
          }

        MaxHsml += MaxDtHsml * dt_drift;

        Ti_Current = time1;
      }

    access.clear(std::memory_order_release);
  }
};

struct foreign_sphpoint_data
{
  sph_particle_data_hydrocore SphCore;

  MyDouble Mass;
  MyIntPosType IntPos[3];

  int Nextnode;
  unsigned char Nextnode_shmrank;

  signed char TimeBinHydro;
};

class ngbtree : public tree<ngbnode, simparticles, ngbpoint_data, foreign_sphpoint_data>
{
 public:
  typedef tree<ngbnode, simparticles, ngbpoint_data, foreign_sphpoint_data> basetree;

  // The various using statements  make sure that we can access the elements from the base class without having to use "this"
  using basetree::Buildtime;
  using basetree::D;
  using basetree::Father;
  using basetree::FirstNonTopLevelNode;
  using basetree::get_nextnodep;
  using basetree::get_nodep;
  using basetree::get_Pp;
  using basetree::get_SphPp;
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
  using basetree::Tp;
  using basetree::tree_export_node_threads;
  using basetree::TreeSharedMem_ThisTask;

  void update_velocities(void);
  void update_maxhsml(void);
  void update_vbounds(int i, int *nchanged, int *nodelist, char *flag_changed);
  void check_bounds(void);

  void update_node_recursive(int no, int sib, int mode) override;
  void exchange_topleafdata(void) override;
  void fill_in_export_points(ngbpoint_data *exp_point, int i, int no) override;
  void report_log_message(void) override;

 private: /* private member functions */
  void finish_vounds_update(int nchanged, int *nodelist);
  void finish_maxhsml_update(int nchanged, int *nodelist);
};

#endif
