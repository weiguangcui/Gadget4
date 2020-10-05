/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file fof.h
 *
 *  \brief declaration of class for the FOF functionality
 */

#ifndef FOF_H
#define FOF_H

#include "gadgetconfig.h"

#ifdef FOF

#include "../data/simparticles.h"
#include "../domain/domain.h"
#include "../fof/foftree.h"
#include "../gravtree/gravtree.h"

/*! \brief FOF routines
 *
 * The template type 'partset' will be given either by 'simparticles' (normal particle set) or by 'lightconeparticles' (for group
 * finding on the lightcone).
 */
template <typename partset>
class fof : public setcomm
{
 public:
  fof(MPI_Comm comm, partset *Tp_ptr, domain<partset> *D_ptr) : setcomm(comm)
  {
    Tp        = Tp_ptr;
    FoFDomain = D_ptr;
  }

  partset *Tp;

  int Ngroups;
  int MaxNgroups;
  long long TotNgroups;

  int Nsubhalos;
  int MaxNsubhalos = 0;
  long long TotNsubhalos;

  long long Nids;
  long long TotNids;

  double Time;
  double Redshift;

  struct group_properties
  {
    long long GroupNr;
    long long OffsetType[NTYPES];
    MyLenType Len;
    MyLenType LenType[NTYPES];
    MyIDType MinID;
    int MinIDTask;
    int OriginTask;

    MyDouble MassType[NTYPES];
    MyDouble Mass;
    MyDouble Ascale;
    MyDouble CM[3];
    MyDouble Pos[3];
    MyIntPosType IntPos[3];
    MyIntPosType FirstIntPos[3];

    MyFloat Vel[3];
#ifdef STARFORMATION
    MyFloat Sfr;
#endif

#ifdef SUBFIND
    int TargetTask; /* primary CPU responsible for running subfind on this group */
    int Nsubs;
    long long FirstSub;
    MyFloat M_Mean200, R_Mean200;
    MyFloat M_Crit200, R_Crit200;
    MyFloat M_Crit500, R_Crit500;
    MyFloat M_TopHat200, R_TopHat200;
#endif
#ifdef MERGERTREE
    long long FileOffset;
#endif
#if defined(SUBFIND_ORPHAN_TREATMENT)
    int LenPrevMostBnd;
#endif
  };
  group_properties *Group;

#ifdef SUBFIND
  struct subhalo_properties
  {
    MyLenType Len;
    MyLenType LenType[NTYPES];
    long long OffsetType[NTYPES]; /*!< gives first particle index of subhalo relative to beginning of file (simplifies work with the
                                       catalogs, as one doesn't need group lengths to know about fuzz particles) */
    long long GroupNr;            /*!< global parent FOF group number */
    int SubRankInGr;              /* local subhalo index within a given FOF group */
    int SubParentRank;
    MyIDType SubMostBoundID;
    MyFloat Mass;
    MyFloat MassType[NTYPES];
    MyFloat SubVelDisp;
    MyFloat SubVmax;
    MyFloat SubVmaxRad;
    MyFloat SubHalfMassRad;
    MyFloat SubHalfMassRadType[NTYPES];
    MyFloat Pos[3];
    MyIntPosType IntPos[3];
    MyFloat CM[3];
    MyFloat Vel[3];
    MyFloat Spin[3];

#ifdef STARFORMATION
    MyFloat Sfr;
    MyFloat GasMassSfr;
#endif
#ifdef MERGERTREE
    long long SubhaloNr;     /* global subhalo number within a given snapshot */
    long long UniqueGroupNr; /* global group number within whole merger tree  */
    long long FileOffset;
    long long TreeID;
    int TreeTask;
    int TreeIndex;
    MyFloat M_Crit200; /* will only be set for main subhalos in halos */
#endif
#if defined(SUBFIND_ORPHAN_TREATMENT)
    int SubhaloLenPrevMostBnd;
#endif
  };
  subhalo_properties *Subhalo;

#if defined(MERGERTREE)
  /* Subhalos as they appear in the merger tree */
  struct treehalo_t
  {
    /* the following pointers are always meant to be relative to the halos in the same tree */
    int TreeDescendant;
    int TreeFirstProgenitor;
    int TreeNextProgenitor;
    int TreeFirstHaloInFOFgroup;
    int TreeNextHaloInFOFgroup;
    int TreeProgenitor;
    int TreeFirstDescendant;
    int TreeNextDescendant;
    int TreeMainProgenitor;

    long long TreeID;
    int TreeIndex;

    /* original position of the subhalo in the subhalo group catalogue output */
    int SnapNum;
    long long SubhaloNr; /* number of subhalo in the full group catalogue of this snapshot */
    long long GroupNr;
    long long UniqueGroupNr;

    /* properties of subhalo */
    subhalo_properties SubProp;
  };
#endif

#endif  // SUBFIND

  void fof_fof(int num, const char *grpcat_basename, const char *grpcat_dirbasename, double inner_distance);
  double fof_get_comoving_linking_length(void);
  void fof_subfind_exchange(MPI_Comm Communicator);
  void fof_subfind_write_file(char *fname, int writeTask, int lastTask, void *CommBuffer);
  double fof_find_nearest_dmparticle(void);
  double fof_find_groups(void);
  void fof_subfind_save_groups(int num, const char *basename, const char *grpcat_dirbasename);
  void fof_subfind_init_io_fields(void);
  void fof_subfind_load_groups(int num);
  void fof_assign_group_offset(void);
  void fof_reorder_PS(int *Id, int Nstart, int N);

 private:
  gravtree<partset> FoFGravTree; /*!< an instance of a gravitational tree */
  foftree<partset> FoFNgbTree;   /*!< an instance of a neighbour search tree */
  domain<partset> *FoFDomain;    /*!< a pointer to the parent domain decomposition */

  int NgroupsExt;

  void fof_compile_catalogue(double inner_distance);
  void fof_compute_group_properties(int gr, int start, int len);
  void fof_prepare_output_order(void);
  void fof_add_in_properties_of_group_segments(void);
  void fof_finish_group_properties(void);
  void fof_assign_group_numbers(void);
  void fof_get_halo_position(MyIntPosType *intpos, double *pos);

#if defined(LIGHTCONE_PARTICLES_GROUPS)
  double fof_distance_to_origin(int i);
#endif

 public:
  struct fof_particle_list
  {
    MyIDStorage MinID;
    int MinIDTask;
    int Pindex;
#if defined(LIGHTCONE_PARTICLES_GROUPS)
    double DistanceOrigin;
#endif
  };
  fof_particle_list *FOF_PList;

  struct fof_group_list
  {
    MyIDStorage MinID;
    int MinIDTask;
    MyLenType Count;
#if defined(LIGHTCONE_PARTICLES_GROUPS)
    double DistanceOrigin;
#endif
  };
  fof_group_list *FOF_GList;

 public:
#ifndef LEAN
  static bool fof_compare_subfind_data_Type(const subfind_data &a, const subfind_data &b) { return a.Type < b.Type; }
#endif

  static bool fof_compare_subfind_data_GroupNr_SubNr_Egy_Key(const subfind_data &a, const subfind_data &b)
  {
    if(a.GroupNr < b.GroupNr)
      return true;
    if(a.GroupNr > b.GroupNr)
      return false;

#ifdef SUBFIND
    if(a.SubRankInGr < b.SubRankInGr)
      return true;
    if(a.SubRankInGr > b.SubRankInGr)
      return false;

    if(a.v.DM_BindingEnergy < b.v.DM_BindingEnergy)
      return true;
    if(a.v.DM_BindingEnergy > b.v.DM_BindingEnergy)
      return false;
#endif

    return a.u.Key < b.u.Key;
  }

#if defined(MERGERTREE) && defined(SUBFIND)

  static bool fof_compare_subfind_data_GroupNr_SubRankInNr_BndEgy(const subfind_data &a, const subfind_data &b)
  {
    if(a.GroupNr < b.GroupNr)
      return true;
    if(a.GroupNr > b.GroupNr)
      return false;

    if(a.SubRankInGr < b.SubRankInGr)
      return true;
    if(a.SubRankInGr > b.SubRankInGr)
      return false;

    return a.v.DM_BindingEnergy < b.v.DM_BindingEnergy;
  }

#endif

  static bool fof_compare_subfind_data_OriginTask_OriginIndex(const subfind_data &a, const subfind_data &b)
  {
    if(a.OriginTask < b.OriginTask)
      return true;
    if(a.OriginTask > b.OriginTask)
      return false;

    return a.OriginIndex < b.OriginIndex;
  }

  static bool fof_compare_FOF_PList_MinID(const fof_particle_list &a, const fof_particle_list &b)
  {
    return a.MinID.get() < b.MinID.get();
  }

  static bool fof_compare_FOF_GList_MinID(const fof_group_list &a, const fof_group_list &b) { return a.MinID.get() < b.MinID.get(); }

  static bool fof_compare_FOF_GList_MinIDTask(const fof_group_list &a, const fof_group_list &b) { return a.MinIDTask < b.MinIDTask; }

  static bool fof_compare_Group_Len_MinID_DiffOriginTaskMinIDTask(const group_properties &a, const group_properties &b)
  {
    if(a.Len > b.Len)
      return true;
    if(a.Len < b.Len)
      return false;

    if(a.MinID < b.MinID)
      return true;
    if(a.MinID > b.MinID)
      return false;

    return labs(a.OriginTask - a.MinIDTask) < labs(b.OriginTask - b.MinIDTask);
  }

  static bool fof_compare_Group_OriginTask_MinID(const group_properties &a, const group_properties &b)
  {
    if(a.OriginTask < b.OriginTask)
      return true;
    if(a.OriginTask > b.OriginTask)
      return false;

    return a.MinID < b.MinID;
  }

  /* now comes the group catalogue created with execution of FOF */

 public:
  static inline bool fof_compare_Group_GroupNr(const group_properties &a, const group_properties &b) { return a.GroupNr < b.GroupNr; }

  static bool fof_compare_Group_MinID(const group_properties &a, const group_properties &b) { return a.MinID < b.MinID; }

  static bool fof_compare_Group_MinIDTask(const group_properties &a, const group_properties &b) { return a.MinIDTask < b.MinIDTask; }

  /********************************************  SubFind part ****************************************/
#ifdef SUBFIND
 public:
  unsigned char *ProcessedFlag;

  unsigned long long GroupNr;
  double Ascale;

  MPI_Comm SubComm;
  int CommSplitColor;
  int SubNTask, SubThisTask;

  int Ncollective;
  int NprocsCollective;
  int MaxSerialGroupLen;

  void subfind_density_hsml_guess(void);
  void subfind_find_subhalos(int num, const char *basename, const char *grpcat_dirbasename);

 public:
  struct sort_r2list
  {
    double r;
    double mass;
  };

  static bool subfind_compare_dist_rotcurve(const sort_r2list &a, const sort_r2list &b) { return a.r < b.r; }

  double subfind_density(void);
  double subfind_overdensity(void);
  double subfind_get_overdensity_value(int type, double ascale);
  void subfind_save_final(int num, const char *basename, const char *grpcat_dirbasename);

  void subfind_processing(domain<partset> *SubDomain, domain_options mode);
  void subfind_potential_compute(domain<partset> *SubDomain, int num, int *d);
  void subfind_find_linkngb(domain<partset> *SubDomain, int num, int *list);
  void subfind_find_nearesttwo(domain<partset> *SubDomain, int num, int *list);
  void subfind_redetermine_groupnr(void);

  void subfind_process_groups_serially(void);
  void subfind_distribute_particles(MPI_Comm Communicator);
  void subfind_distribute_groups(void);
  double subfind_get_particle_balance(void);
  void subfind_assign_subhalo_offsettype(void);
  void subfind_match_ids_of_previously_most_bound_ids(partset *Tp);

  double subfind_locngb_treefind(MyDouble xyz[3], int desngb, double hguess);
  int subfind_unbind(domain<partset> *D, MPI_Comm Communicator, int *unbind_list, int len);

  int subfind_determine_sub_halo_properties(int *d, int num, subhalo_properties *subhalo, MPI_Comm Communicator);

  void subfind_hbt_single_group(domain<partset> *SubDomain, domain<partset> *SingleDomain, domain_options mode, int gr);

  struct proc_assign_data
  {
    long long GroupNr;
    MyLenType Len;
    int FirstTask;
    int NTask;
  };
  proc_assign_data *ProcAssign;

  struct submp_data
  {
    long long GroupNr;
    int index;

#ifndef SUBFIND_HBT
    MyFloat DM_Density;
#endif
  };
  submp_data *submp;

  struct cand_dat
  {
    location head;
    MyLenType len;
    MyLenType bound_length;

    int nsub;
    int rank, subnr, parent;
  };

  struct coll_cand_dat
  {
    location head;
    MyLenType rank;
    MyLenType len;
    MyLenType bound_length;

    int nsub;
    int subnr, parent;
  };
  coll_cand_dat *coll_candidates;

  struct SubDMData
  {
    double rho;
    double vx, vy, vz;
    double v2;
  };

  static inline bool subfind_compare_submp_GroupNr_DM_Density(const submp_data &a, const submp_data &b)
  {
#ifndef SUBFIND_HBT
    if(a.GroupNr < b.GroupNr)
      return true;
    if(a.GroupNr > b.GroupNr)
      return false;

    return (a.DM_Density > b.DM_Density);
#else
    return a.GroupNr < b.GroupNr;
#endif
  }

  static inline bool subfind_compare_binding_energy(const double &a, const double &b) { return a > b; }

  static inline bool subfind_compare_potential(const double &a, const double &b) { return a < b; }

  static inline bool subfind_compare_Subhalo_GroupNr_SubRankInGr(const subhalo_properties &a, const subhalo_properties &b)
  {
    if(a.GroupNr < b.GroupNr)
      return true;
    if(a.GroupNr > b.GroupNr)
      return false;

    return a.SubRankInGr < b.SubRankInGr;
  }

  static inline bool subfind_compare_procassign_GroupNr(const proc_assign_data &a, const proc_assign_data &b)
  {
    return a.GroupNr < b.GroupNr;
  }

 private:
  long long count_decisions;
  long long count_different_decisions;

  struct sort_density_data
  {
    MyFloat density;
    int ngbcount;
    location index; /* this will store the task in the upper word */
    location ngb_index1;
    location ngb_index2;
    approxlen PrevSizeOfSubhalo;
  };
  sort_density_data *sd;

  struct PPS_data
  {
    int index;
    int submark;
  };
  PPS_data *PPS;

  void subfind_col_find_coll_candidates(long long totgrouplen);

  void subfind_poll_for_requests(void);

  int subfind_distlinklist_get_tail_set_tail_increaselen(location index, location &tail, location newtail, approxlen prevlen);

  void subfind_distlinklist_get_two_heads(location ngb_index1, location ngb_index2, location &head, location &head_attach);
  void subfind_distlinklist_set_next(location index, location next);
  void subfind_distlinklist_add_particle(location index);
  void subfind_distlinklist_add_bound_particles(location index, int nsub);
  void subfind_distlinklist_mark_particle(location index, int target, int submark);
  void subfind_distlinklist_set_headandnext(location index, location head, location next);
  void subfind_distlinklist_set_tailandlen(location index, location tail, MyLenType len, double prevlen);
  void subfind_distlinklist_get_tailandlen(location index, location &tail, MyLenType &len, double &prevlen);
  void subfind_distlinklist_set_all(location index, location head, location tail, MyLenType len, location next, approxlen prevlen);
  location subfind_distlinklist_get_next(location index);
  location subfind_distlinklist_get_head(location index);
  location subfind_distlinklist_setrank_and_get_next(location index, MyLenType &rank);
  location subfind_distlinklist_set_head_get_next(location index, location head);
  MyLenType subfind_distlinklist_get_rank(location index);

  void subfind_get_factors(double &fac_vel_to_phys, double &fac_hubbleflow, double &fac_comov_to_phys);

  void subfind_process_single_group(domain<partset> *SubDomain, domain<partset> *SingleDomain, domain_options mode, int gr);
  void subfind_unbind_independent_ones(domain<partset> *SingleDomain, int count);
  double subfind_vel_to_phys_factor(void);

  void subfind_collective_printf(const char *fmt, ...)
  {
    if(SubNTask > 1 && SubThisTask == 0)
      {
        va_list l;
        va_start(l, fmt);
        vprintf(fmt, l);
        va_end(l);
      }
  }

  static bool subfind_compare_densities(const sort_density_data &a, const sort_density_data &b) /* largest density first */
  {
    return a.density > b.density;
  }

  static bool subfind_PS_compare_DM_density(const subfind_data &a, const subfind_data &b) /* largest density first */
  {
    return a.u.s.u.DM_Density > b.u.s.u.DM_Density;
  }

  static bool subfind_PS_compare_origintask_originindex(const subfind_data &a, const subfind_data &b)
  {
    if(a.u.s.origintask < b.u.s.origintask)
      return true;
    if(a.u.s.origintask > b.u.s.origintask)
      return false;

    return a.u.s.originindex < b.u.s.originindex;
  }

  static bool subfind_compare_coll_candidates_subnr(const coll_cand_dat &a, const coll_cand_dat &b) { return a.subnr < b.subnr; }

  static bool subfind_compare_coll_candidates_nsubs(const coll_cand_dat &a, const coll_cand_dat &b) { return a.nsub < b.nsub; }

  static bool subfind_compare_coll_candidates_boundlength(const coll_cand_dat &a, const coll_cand_dat &b)
  {
    if(a.bound_length > b.bound_length)
      return true;
    if(a.bound_length < b.bound_length)
      return false;

    return a.rank < b.rank;
  }

  static bool subfind_compare_coll_candidates_rank(const coll_cand_dat &a, const coll_cand_dat &b)
  {
    if(a.rank < b.rank)
      return true;
    if(a.rank > b.rank)
      return false;
    return a.len > b.len;
  }

  static bool subfind_compare_PPS(const PPS_data &a, const PPS_data &b) { return a.submark < b.submark; }

  MyLenType *SFLen;
  location *SFHead;
  location *SFNext;
  location *SFTail;
  double *SFPrevLen;

  int count_cand, max_coll_candidates;
  int *unbind_list;

  int NumPartGroup;
  int *IndexList;
  int LocalLen;

  struct sort_as_data
  {
    double density;
    int targettask;
    long long origin;
  };

  static bool subfind_compare_as_density(const sort_as_data &a, const sort_as_data &b) /* largest density first */
  {
    return a.density > b.density;
  }

  static bool subfind_compare_as_origin(const sort_as_data &a, const sort_as_data &b) /* largest density first */
  {
    return a.origin < b.origin;
  }

  struct hbt_pcand_t
  {
    MyHaloNrType SubhaloNr;
    approxlen PrevSizeOfSubhalo;
    int index;
  };

  static bool subfind_hbt_compare_pcand_subhalonr(const hbt_pcand_t &a, const hbt_pcand_t &b) { return a.SubhaloNr < b.SubhaloNr; }

  struct hbt_subcand_t
  {
    MyHaloNrType SubhaloNr;
    MyLenType len;
    bool DoIt;

    int TargetTask;
    int TargetIndex;
    long long summedprevlen;

    /*
    location head;
    MyLenType rank;

    MyLenType bound_length;

    int nsub;
    int subnr, parent;
    */
  };

  static bool subfind_hbt_compare_subcand_subhalonr(const hbt_subcand_t &a, const hbt_subcand_t &b)
  {
    return a.SubhaloNr < b.SubhaloNr;
  }

  static bool subfind_hbt_compare_subcand_len(const hbt_subcand_t &a, const hbt_subcand_t &b) { return a.len < b.len; }

  static bool subfind_hbt_compare_subcand_summedprevlen(const hbt_subcand_t &a, const hbt_subcand_t &b)
  {
    return a.summedprevlen < b.summedprevlen;
  }

  struct hbt_subhalo_t
  {
    MyLenType Len;
    int SubRankInGr;
    int ThisTask;
    int ThisIndex;
    long long SubhaloNr;
  };

  static bool subfind_hbt_compare_subhalolist_len(const hbt_subhalo_t &a, const hbt_subhalo_t &b) { return a.Len > b.Len; }

  static bool subfind_hbt_compare_subhalolist_thistask_thisindex(const hbt_subhalo_t &a, const hbt_subhalo_t &b)

  {
    if(a.ThisTask < b.ThisTask)
      return true;
    if(a.ThisTask > b.ThisTask)
      return false;

    return a.ThisIndex < b.ThisIndex;
  }

  static bool subfind_hbt_compare_subhalolist_prevsubhalonr(const hbt_subhalo_t &a, const hbt_subhalo_t &b)
  {
    return a.SubhaloNr < b.SubhaloNr;
  }

#endif
};

inline bool is_type_primary_link_type(int type)
{
  if((1 << type) & (FOF_PRIMARY_LINK_TYPES))
    return true;
  else
    return false;
}

inline bool is_type_secondary_link_type(int type)
{
  if((1 << type) & (FOF_SECONDARY_LINK_TYPES))
    return true;
  else
    return false;
}
#endif  // end of FOF

#endif
