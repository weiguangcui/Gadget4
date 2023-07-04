/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  mergertree.h
 *
 *  \brief definition of a class that provides the merger tree functionality
 */

#ifndef MERGERTREE_H
#define MERGERTREE_H

#ifdef MERGERTREE

#include "gadgetconfig.h"

#include <hdf5.h>

#include "../data/simparticles.h"
#include "../fof/fof.h"
#include "../subfind/subfind.h"

#define NUM_MOST_BOUND_PARTICLES_USED_FOR_TRACKING 16  // may not be larger than 255 due to storage limits

class mergertree : public setcomm
{
 private:
  simparticles *Sp;

 public:
  mergertree(MPI_Comm comm, simparticles *Sp_ptr) : setcomm(comm)
  {
    Sp               = Sp_ptr;
    PrevTotNsubhalos = 0;
    PrevNsubhalos    = 0;
    MtrP_NumPart     = 0;
  }

  void mergertree_determine_descendants_on_the_fly(int num);
  void descendants_in_postprocessing(int num);
  void halotrees_construct(int lastsnapnr);
  void get_previous_size_of_subhalo_for_each_particle(int num);

  unsigned long long CurrTotNsubhalos;
  unsigned long long PrevTotNsubhalos;

  int CurrNsubhalos;
  int PrevNsubhalos;

  int MtrP_NumPart;
  int PrevMtrP_NumPart;

  int Nhalos;
  long long TotNhalos;

  int Ntrees;
  long long TotNtrees;

  int LargestHaloCount;
  int LastSnapShotNr;  // we will construct trees composed of dumps 0 to LastSnapShotNr

  typedef fof<simparticles>::treehalo_t treehalo_type;
  treehalo_type *Halos; /* This will contain the subhalos making up the merger tree information */

  struct treehalo_ids_type
  {
    MyIDType SubMostBoundID;
    long long TreeID;
  };
  treehalo_ids_type *HaloIDdata;

  static bool compare_HaloIDdata_ID(const treehalo_ids_type &a, const treehalo_ids_type &b)
  {
    return a.SubMostBoundID < b.SubMostBoundID;
  }

  struct desc_list
  {
    double MaxScore;
    long long PrevSubhaloNr;
    long long DescSubhaloNr;

    long long FirstDescSubhaloNr;
    long long NextProgSubhaloNr;

    long long FileOffset;
    char FirstProgFlag;
  };
  desc_list *Descendants;

  struct prog_list
  {
    double MaxScore;
    long long SubhaloNr;

    long long ProgSubhaloNr;

    long long FirstProgSubhaloNr;
    long long NextDescSubhaloNr;

    long long MainProgSubhaloNr;

    long long FileOffset;

    char FirstDescFlag;
  };
  prog_list *Progenitors;

  struct mergertree_particle_data
  {
    int Type;
    MyIDType ID;
    long long FileOffset;

    long long PrevGroupNr;
    long long PrevSubhaloNr;

    long long GroupNr;
    long long SubhaloNr;

    MyLenType SubhaloLen;
    MyLenType PrevSubhaloLen;

    unsigned short RankInSubhalo;
    unsigned short PrevRankInSubhalo;
  };
  mergertree_particle_data *MtrP, *PrevMtrP;

  halotrees_table *TreeTable;

  struct treelink_data
  {
    long long TreeID;
    int TreeIndex;
  };
  treelink_data *TreeLink;

  times_catalogue *CatTimes;

 private:
  void mergertree_determine_descendants(int num);
  void mergertree_determine_descendants_postproc(int num);
  void mergertree_save_progenitors(int num);
  void mergertree_read_progenitors(int num);
  void desc_subfind_init_io_fields(void);
  void prog_subfind_init_io_fields(void);
  void treelink_subfind_init_io_fields(void);
  void mergertree_read_snap_ids(int num);
  void mergertree_match_ids_of_previous_snap(void);
  void mrgtr_init_io_fields(void);
  int halotrees_join_via_descendants(int num);
  int halotrees_join_via_progenitors(int num);
  void halotrees_assign_global_subhalonr_and_groupnr(void);
  void halotrees_load_catalogues(fof<simparticles> *FoF);
  void halotrees_initial_treeassignment(void);
  void halotrees_assign_new_treeindices(void);
  void halotrees_save_trees(void);
  void halotrees_save_subhalo_treelinks(void);
  void halotrees_collect_treehalos(void);
  void halotrees_link_trees(void);
  void halotrees_propagate_max_branch_length_descendants(int num);
  void halotrees_propagate_max_branch_length_progenitors(int num);
  void halotrees_determine_mainprogenitor(void);
  void halotrees_reshuffle(char **ptr, size_t len, int ncurrent, int ntarget);
  void halotrees_remap_treepointers(void);
  void mergertree_assign_group_numbers(fof<simparticles> *FoF);
  void trees_subfind_init_io_fields(void);
  long long halotrees_join_trees_via_fof_or_mostboundid_bridges(int mode);
  void mergertree_match_ids_of_current_snap(void);

  void mergertree_chain_up_progenitors_with_same_descendant(void);
  void mergertree_set_first_progenitor_with_same_descendant(void);
  void mergertree_select_maximum_score_progenitors(int nmatch);
  void mergertree_select_maximum_score_descendants(int nmatch);
  int mergertree_find_matching_segments_and_scores(void);
  void mergertree_chain_up_descendants_with_same_progenitor(void);
  void mergertree_set_first_descendant_with_same_progenitor(void);

  struct subhalo_extension
  {
    int TreeDescendant;
    int TreeFirstProgenitor;
    int TreeNextProgenitor;
    int TreeFirstHaloInFOFgroup;
    int TreeNextHaloInFOFgroup;
    int TreeProgenitor;
    int TreeFirstDescendant;
    int TreeNextDescendant;
    int TreeMainProgenitor;

    long long MaxLenProgBranch;
  };

  /* PrevTotNsubhalos is the total number of subhalos in the previous catalogue. This list is split onto
   * the processors, with 'PrevNsubhalos' being the number of these subhalos assigned to the local MPI task.
   */

  static bool compare_MtrP_FileOffset(const mergertree_particle_data &a, const mergertree_particle_data &b)
  {
    return a.FileOffset < b.FileOffset;
  }

  static bool compare_MtrP_SubhaloNr(const mergertree_particle_data &a, const mergertree_particle_data &b)
  {
    if(a.PrevSubhaloNr < b.PrevSubhaloNr)
      return true;
    if(a.PrevSubhaloNr > b.PrevSubhaloNr)
      return false;

    if(a.PrevRankInSubhalo < b.PrevRankInSubhalo)
      return true;
    if(a.PrevRankInSubhalo > b.PrevRankInSubhalo)
      return false;

    return a.ID < b.ID;
  }

  static bool compare_Group_FileOffset(const fof<simparticles>::group_properties &a, const fof<simparticles>::group_properties &b)
  {
    return a.FileOffset < b.FileOffset;
  }

  static bool compare_Subhalo_FileOffset(const fof<simparticles>::subhalo_properties &a,
                                         const fof<simparticles>::subhalo_properties &b)
  {
    return a.FileOffset < b.FileOffset;
  }

  static bool compare_MtrP_ID(const mergertree_particle_data &a, const mergertree_particle_data &b) { return a.ID < b.ID; }

  struct assign_particle_data
  {
    MyIDType ID;
    long long PrevSubhaloNr;
    MyLenType PrevSubhaloLen;
    int OriginTask;
    int OriginIndex;
  };

  static bool compare_AssignP_ID(const assign_particle_data &a, const assign_particle_data &b) { return a.ID < b.ID; }

  static bool compare_AssignP_Origin(const assign_particle_data &a, const assign_particle_data &b)
  {
    if(a.OriginTask < b.OriginTask)
      return true;
    if(a.OriginTask > b.OriginTask)
      return false;

    return a.OriginIndex < b.OriginIndex;
  }

  /* This structure is used to do the matching/identification of the descendants. First, every particle is listed here
   * with its subhalo number and rank in the previous subhalo catalogue, while GroupNr/SubRankInGr (or equivalently the combined
   * SubhaloNr) give the number of the subhalo they are in the new catalogue. By first sorting according to PrevSubhaloNr
   * and then by SubhaloNr, we get all descendant candidates a previous subhalo can have. We then need to figure out
   * which one to pick, which we do via the Score variable.
   */
  struct desc_partdata
  {
    MyHaloNrType PrevSubhaloNr;  // Global subhalo number of the particle in the previous group catalogue
    MyHaloNrType CurrSubhaloNr;  // New global subhalo number in the current group catalogue

    float DescScore;  // auxiliary variable to score the different possible descendants
    float ProgScore;  // auxiliary variable to score the different possible progenitors

    compactrank_t PrevRankInSubhalo;  // Rank of particle within previous subhalo
    compactrank_t CurrRankInSubhalo;  // Rank of particle within current subhalo
  };
  desc_partdata *desc;

  /* sort kernel, first by PrevSubhaloNr, then by SubhaloNr */
  static bool mergertree_compare_PrevSubNr_NewSubNr(const desc_partdata &a, const desc_partdata &b)
  {
    if(a.PrevSubhaloNr < b.PrevSubhaloNr)
      return true;
    if(a.PrevSubhaloNr > b.PrevSubhaloNr)
      return false;

    return a.CurrSubhaloNr < b.CurrSubhaloNr;
  }

  static bool mergertree_compare_NewSubNr_PrevSubNr(const desc_partdata &a, const desc_partdata &b)
  {
    if(a.CurrSubhaloNr < b.CurrSubhaloNr)
      return true;
    if(a.CurrSubhaloNr > b.CurrSubhaloNr)
      return false;

    return a.PrevSubhaloNr < b.PrevSubhaloNr;
  }

  /* sort kernel */
  static bool mergertree_compare_DescSubhaloNr(const desc_list &a, const desc_list &b)
  {
    if(a.DescSubhaloNr < b.DescSubhaloNr)
      return true;
    if(a.DescSubhaloNr > b.DescSubhaloNr)
      return false;

    return a.PrevSubhaloNr < b.PrevSubhaloNr;
  }

  /* sort kernel */
  static bool mergertree_compare_ProgSubhaloNr(const prog_list &a, const prog_list &b)
  {
    if(a.ProgSubhaloNr < b.ProgSubhaloNr)
      return true;
    if(a.ProgSubhaloNr > b.ProgSubhaloNr)
      return false;

    return a.SubhaloNr < b.SubhaloNr;
  }

  static bool mergertree_compare_SubhaloNr(const prog_list &a, const prog_list &b) { return a.SubhaloNr < b.SubhaloNr; }

  /* sort kernel */
  static bool mergertree_compare_PrevSubhaloNr(const desc_list &a, const desc_list &b) { return a.PrevSubhaloNr < b.PrevSubhaloNr; }

  /* This structure is needed to organize the information of all the group and subhalo catalogs that are read in.
   * The array Cats[] holds a table with all the group catalogs, indexes by the snapshot number.
   */
  struct halo_catalogue
  {
    fof<simparticles>::group_properties *Group;  // table of local groups
    long long FirstGroup;                        // first group number on this processor
    long long TotNgroups;                        // total number of groups in this catalog
    int Ngroups;                                 // number of groups stored on local processor
    int *TabNgroups;                             // used to store a table with the group numbers on all processors

    fof<simparticles>::subhalo_properties *Subhalo;  // table of local subhalos
    long long FirstSubhalo;                          // first subhalo number on this processor
    long long TotNsubhalos;                          // total number of subhalos in catalog
    int Nsubhalos;                                   // local number of subhalos on this processor
    int *TabNsubhalos;                               // used to store a table with the subhalo numbers on all processors

    subhalo_extension *SubExt;  // additional subhalo properties that we determine as part of the tree building

    desc_list *Descendants;  // stores descendant information
    prog_list *Progenitors;  // stores progenitor information
  };
  halo_catalogue *Cats;

  struct tlink
  {
    long long UniqueGroupNr;
    long long OrderIndex;

    long long TreeID;
    long long NewTreeID;

    int TreeTask;
    int NewTreeTask;

    int OrigTask;
  };

  static bool compare_tlink_GroupNr(const tlink &a, const tlink &b) { return a.UniqueGroupNr < b.UniqueGroupNr; }

  static bool compare_tlink_TreeID(const tlink &a, const tlink &b) { return a.TreeID < b.TreeID; }

  static bool compare_tlink_OrigTask_OrderIndex(const tlink &a, const tlink &b)
  {
    if(a.OrigTask < b.OrigTask)
      return true;
    if(a.OrigTask > b.OrigTask)
      return false;

    return a.OrderIndex < b.OrderIndex;
  }

  static bool compare_Halos_TreeID_TreeIndex(const treehalo_type &a, const treehalo_type &b)
  {
    if(a.TreeID < b.TreeID)
      return true;
    if(a.TreeID > b.TreeID)
      return false;

    return a.TreeIndex < b.TreeIndex;
  }

  static bool compare_Halos_UniqueGroupNr_SubhaloNr(const treehalo_type &a, const treehalo_type &b)
  {
    if(a.UniqueGroupNr < b.UniqueGroupNr)
      return true;
    if(a.UniqueGroupNr > b.UniqueGroupNr)
      return false;

    return a.SubhaloNr < b.SubhaloNr;
  }

  static bool compare_Desc_FileOffset(const desc_list &a, const desc_list &b) { return a.FileOffset < b.FileOffset; }

  static bool compare_Prog_FileOffset(const prog_list &a, const prog_list &b) { return a.FileOffset < b.FileOffset; }

  /*--------------------------------------------------------------------------------------------------------------*/

  struct data_list
  {
    long long targetsubhalonr;
    long long intreeid;
    long long originsubhalonr;
    int origin;
  };

  static bool compare_data_list_subhalonnr(const data_list &a, const data_list &b) { return a.targetsubhalonr < b.targetsubhalonr; }

  struct remap_data
  {
    int loc_index;
    int new_treeindexptr;
    long long targetsubhalonr;
    long long originsubhalonr;
    long long treeid;
    long long intreeid;
    int orig_index;
  };

  static bool compare_remap_data_subhalonr(const remap_data &a, const remap_data &b) { return a.targetsubhalonr < b.targetsubhalonr; }

  static bool compare_remap_data_orig_index(const remap_data &a, const remap_data &b) { return a.orig_index < b.orig_index; }

  /*--------------------------------------------------------------------------------------------------------------*/

  /* the following one is used to assign consecutive indices within each tree */
  struct assign_data
  {
    int origin_task;
    int origin_num;
    int origin_index;

    long long treeid;
    long long newtreeid;
    long long treeindex;
  };

  /* some sort kernels */
  static bool compare_assign_data_treeid_origin_num_origin_task_origin_index(const assign_data &a, const assign_data &b)
  {
    if(a.treeid < b.treeid)
      return true;
    if(a.treeid > b.treeid)
      return false;

    if(a.origin_num > b.origin_num)
      return true;
    if(a.origin_num < b.origin_num)
      return false;

    if(a.origin_task < b.origin_task)
      return true;
    if(a.origin_task > b.origin_task)
      return false;

    return a.origin_index < b.origin_index;
  }

  static bool compare_assign_data_origin_task_origin_num_origin_index(const assign_data &a, const assign_data &b)
  {
    if(a.origin_task < b.origin_task)
      return true;
    if(a.origin_task > b.origin_task)
      return false;

    if(a.origin_num < b.origin_num)
      return true;
    if(a.origin_num > b.origin_num)
      return false;

    return a.origin_index < b.origin_index;
  }

  /*--------------------------------------------------------------------------------------------------------------*/

  struct halotrees_data
  {
    long long descendantnr;
    long long progenitornr;
    int loc_index;
    int orig_order;

    long long treeid;
    int treetask;
  };

  static bool compare_halotrees_data_descendantnr(const halotrees_data &a, const halotrees_data &b)
  {
    return a.descendantnr < b.descendantnr;
  }

  static bool compare_halotrees_data_progenitornr(const halotrees_data &a, const halotrees_data &b)
  {
    return a.progenitornr < b.progenitornr;
  }

  static bool compare_halotrees_data_orig_order(const halotrees_data &a, const halotrees_data &b)
  {
    return a.orig_order < b.orig_order;
  }

  struct halotrees_propagate_data
  {
    long long DescSubhaloNr;
    long long ProgSubhaloNr;
    long long SubhaloNr;
    long long MaxLenProgBranch;
    int index;
    int orig_order;
  };

  static bool compare_halotrees_propagate_data_orig_order(const halotrees_propagate_data &a, const halotrees_propagate_data &b)
  {
    return a.orig_order < b.orig_order;
  }

  static bool compare_halotrees_propagate_data_DescSubhaloNr(const halotrees_propagate_data &a, const halotrees_propagate_data &b)
  {
    return a.DescSubhaloNr < b.DescSubhaloNr;
  }

  static bool compare_halotrees_propagate_data_ProgSubhaloNr(const halotrees_propagate_data &a, const halotrees_propagate_data &b)
  {
    return a.ProgSubhaloNr < b.ProgSubhaloNr;
  }

  struct halotrees_firstprog_data
  {
    long long DescSubhaloNr;
    long long SubhaloNr;
  };

  static bool compare_halotrees_firstprog_data_DescSubhaloNr(const halotrees_firstprog_data &a, const halotrees_firstprog_data &b)
  {
    return a.DescSubhaloNr < b.DescSubhaloNr;
  }

  struct descnr_data
  {
    long long DescSubhaloNr;
    long long SubhaloNr;
    long long CumulLen;
    long long TreeID;
    int TreeTask;
    int orig_index;
  };

  static bool compare_sorted_list_descsubhalonr(const descnr_data &a, const descnr_data &b)
  {
    return a.DescSubhaloNr < b.DescSubhaloNr;
  }

  struct prognr_data
  {
    long long ProgSubhaloNr;
    long long SubhaloNr;
    long long TreeID;
    int TreeTask;
    int orig_index;
  };

  static bool compare_sorted_list_progsubhalonr(const prognr_data &a, const prognr_data &b)
  {
    return a.ProgSubhaloNr < b.ProgSubhaloNr;
  }

  void halotrees_select_interior_min_newtreeid(int mode, tlink *treehalos, long long totnsubs);
};

#endif
#endif
