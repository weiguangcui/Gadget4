/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file gravtree.h
 *
 *  \brief declares the gravitational tree node structure and a class for gravity tree calculation
 */

#ifndef FORCETREE_H
#define FORCETREE_H

#include "gadgetconfig.h"

#ifndef NTAB
#define NTAB 256 /* size of short-range look up */
#endif

#define TREE_DO_BASE_PM 1
#define TREE_DO_HIGHRES_PM 2
#define TREE_ACTIVE_CUTTOFF_BASE_PM 4
#define TREE_ACTIVE_CUTTOFF_HIGHRES_PM 8

#define TREE_MIN_WORKSTACK_SIZE 100000
#define TREE_EXPECTED_CYCLES 80

#define NODE_TYPE_LOCAL_PARTICLE 0
#define NODE_TYPE_TREEPOINT_PARTICLE 1
#define NODE_TYPE_FETCHED_PARTICLE 2
#define NODE_TYPE_LOCAL_NODE 3
#define NODE_TYPE_FETCHED_NODE 4

#define NODE_USE 0
#define NODE_OPEN 1
#define NODE_DISCARD 2

#define MAX_TREE_ALLOC_FACTOR 30.0

#define TAKE_NSLOTS_IN_ONE_GO 32

#include <string.h>

#include "../data/simparticles.h"
#include "../data/symtensors.h"
#include "../domain/domain.h"
#include "../mpi_utils/mpi_utils.h"
#include "../tree/tree.h"

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

struct gravnode : public basenode
{
  MyDouble mass;          /**< mass of node */
  vector<MyIntPosType> s; /**< center of mass of node (in integer coordinates!) */
#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
  symtensor2<MyDouble> Q2Tensor; /**< quadrupole tensor */
#endif
#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM))
  symtensor3<MyDouble> Q3Tensor; /**< octupole tensor */
#endif
#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM))
  symtensor4<MyDouble> Q4Tensor; /**< hexadecupole tensor */
#endif
#if(MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM))
  symtensor5<MyDouble> Q5Tensor; /**< triakontadipole tensor */
#endif
#ifdef FMM
  float MinOldAcc; /**< minimum magnitude of old gravitational force. Used in relative opening criterion */
#endif
#if(NSOFTCLASSES > 1)
  unsigned char maxsofttype; /**< hold the maximum gravitational softening of particles */
  unsigned char minsofttype; /**< hold the minimum gravitational softening of particles */
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  unsigned char overlap_flag : 2;
#endif

  inline int getSofteningClass(void)
  {
#if NSOFTCLASSES > 1
    return maxsofttype;
#else
    return 0;
#endif
  }
};

struct gravpoint_data
{
  MyIntPosType IntPos[3];
  MyDouble Mass;
  float OldAcc;
  int index;
  int no;
  unsigned char Type;
#if NSOFTCLASSES > 1
  unsigned char SofteningClass : 7;
#endif
#ifndef HIERARCHICAL_GRAVITY
  unsigned char ActiveFlag : 1; /* we don't need this for hierarchical gravity as then the particles are always active */
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  unsigned char InsideOutsideFlag : 1;
#endif

  inline unsigned char getSofteningClass(void)
  {
#if NSOFTCLASSES > 1
    return SofteningClass;
#else
    return 0;
#endif
  }
};

struct foreign_gravpoint_data
{
  MyIntPosType IntPos[3];
  MyDouble Mass;
  int Nextnode;
  unsigned char Nextnode_shmrank;
  unsigned char Type;
  float OldAcc;
#if NSOFTCLASSES > 1
  unsigned char SofteningClass;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  unsigned char InsideOutsideFlag : 1;
#endif

  inline unsigned char getSofteningClass(void)
  {
#if NSOFTCLASSES > 1
    return SofteningClass;
#else
    return 0;
#endif
  }
};

#define HIGHEST_NEEDEDORDER_EWALD_DPHI 1

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM) || (MULTIPOLE_ORDER >= 2 && defined(FMM)))
#undef HIGHEST_NEEDEDORDER_EWALD_DPHI
#define HIGHEST_NEEDEDORDER_EWALD_DPHI 2
#endif

#if(MULTIPOLE_ORDER >= 4) || (MULTIPOLE_ORDER >= 3 && defined(EXTRAPOTTERM) || (MULTIPOLE_ORDER >= 3 && defined(FMM)))
#undef HIGHEST_NEEDEDORDER_EWALD_DPHI
#define HIGHEST_NEEDEDORDER_EWALD_DPHI 3
#endif

#if(MULTIPOLE_ORDER >= 5) || (MULTIPOLE_ORDER >= 4 && defined(EXTRAPOTTERM) || (MULTIPOLE_ORDER >= 4 && defined(FMM)))
#undef HIGHEST_NEEDEDORDER_EWALD_DPHI
#define HIGHEST_NEEDEDORDER_EWALD_DPHI 4
#endif

#if((MULTIPOLE_ORDER >= 5 && defined(EXTRAPOTTERM) && defined(EVALPOTENTIAL)) || (MULTIPOLE_ORDER >= 5 && defined(FMM)) || \
    defined(EWALD_TEST))
#undef HIGHEST_NEEDEDORDER_EWALD_DPHI
#define HIGHEST_NEEDEDORDER_EWALD_DPHI 5
#endif

#ifdef EXTRA_HIGH_EWALD_ACCURACY
#define EWALD_TAYLOR_ORDER 3
#else
#define EWALD_TAYLOR_ORDER 2
#endif

/* variables for Ewald correction lookup table */
struct ewald_data
{
  MyReal D0phi;
  vector<MyReal> D1phi;
  symtensor2<MyReal> D2phi;
  symtensor3<MyReal> D3phi;
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 4
  symtensor4<MyReal> D4phi;
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 5
  symtensor5<MyReal> D5phi;
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 6
  symtensor6<MyReal> D6phi;
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 7
  symtensor7<MyReal> D7phi;
#endif
};

template <typename partset> /* partset will either be 'simparticles' or 'lightconeparticles'  */
class gravtree : public tree<gravnode, partset, gravpoint_data, foreign_gravpoint_data>
{
 public:
  typedef tree<gravnode, partset, gravpoint_data, foreign_gravpoint_data> basetree;
  using basetree::Buildtime;
  using basetree::D;  // this avoids that we have to use a "this->" when accessing these variables from the base class
  using basetree::EndOfForeignNodes;
  using basetree::EndOfTreePoints;
  using basetree::Father;
  using basetree::FirstNonTopLevelNode;
  using basetree::Foreign_Nodes;
  using basetree::get_nodep;
  using basetree::ImportedNodeOffset;
  using basetree::IndexList;
  using basetree::MaxForeignNodes;
  using basetree::MaxNodes;
  using basetree::MaxPart;
  using basetree::Nextnode;
  using basetree::Ninsert;
  using basetree::NodeIndex;
  using basetree::Nodes;
  using basetree::NumForeignNodes;
  using basetree::NumForeignPoints;
  using basetree::NumNodes;
  using basetree::NumPartExported;
  using basetree::NumPartImported;
  using basetree::Points;
  using basetree::Recv_count;
  using basetree::Recv_offset;
  using basetree::ResultIndexList;
  using basetree::Send_count;
  using basetree::Send_offset;
  using basetree::TopNodes;
  using basetree::Tp;
  using basetree::TreeForeign_Nodes_offsets;
  using basetree::TreeForeign_Points_offsets;
  using basetree::TreeNextnode_offsets;
  using basetree::TreeNodes_offsets;
  using basetree::TreeP_offsets;
  using basetree::TreePoints_offsets;
  using basetree::TreePS_offsets;
  using basetree::TreeSharedMem_ThisTask;
  using basetree::TreeSharedMemBaseAddr;
  using basetree::TreeSharedMemComm;
  using basetree::TreeSphP_offsets;

  struct mesh_factors
  {
    double rcut2;
    double asmthfac;

    double asmthinv1;
    double asmthinv2;
#if(MULTIPOLE_ORDER >= 2)
    double asmthinv3;
#endif
#if(MULTIPOLE_ORDER >= 3)
    double asmthinv4;
#endif
#if(MULTIPOLE_ORDER >= 4)
    double asmthinv5;
#endif
#if(MULTIPOLE_ORDER >= 5)
    double asmthinv6;
#endif
    MyIntPosType intrcut[3];
  };
  mesh_factors mf[2];

  struct resultsactiveimported_data
  {
    vector<MyFloat> GravAccel;
#ifdef EVALPOTENTIAL
    MyFloat Potential;
#endif
    int GravCost;
    int index;
  };

  struct index_data
  {
    int p;
    int subnode;
  };

  /** Variables for gravitational tree */

  char MeasureCostFlag;

  resultsactiveimported_data *ResultsActiveImported;

  ewald_data PotTaylor;
  int num_layers = 0;

 private:
  /** Gives next node in tree walk for the "particle" nodes. Entries 0
   -- MaxPart-1 are the real particles, and the "pseudoparticles" are
   indexed by the node number-MaxNodes. */

  /** Gives previous node in tree walk for the leaf (particle)
   nodes. Entries 0 -- MaxPart-1 are the real particles, and the
   "pseudoparticles" are indexed by the node number-MaxNodes. */

 private:
  /* private member functions */

  void update_node_recursive(int no, int sib, int mode) override;
  void exchange_topleafdata(void) override;
  void fill_in_export_points(gravpoint_data *exp_point, int i, int no) override;
  void report_log_message(void) override;

 public:
  void set_softenings(void);

#ifdef ALLOW_DIRECT_SUMMATION
  void gravity_direct(partset *Pptr, domain<partset> *Dptr, int timebin);
#endif

  void gravity_exchange_forces(void);

  /** public functions */
 public:
#ifdef PMGRID
  void short_range_init(void);
  void set_mesh_factors(void);
#endif

 public:
#ifdef PMGRID
  /*! \brief variable for short-range lookup table
   *
   *  contains the factor needed for the short range
   *  contribution of the tree to the gravity force
   */
  struct
  {
    float fac0;
    float fac1;
#if(MULTIPOLE_ORDER >= 2)
    float fac2;
#endif
#if(MULTIPOLE_ORDER >= 3)
    float fac3;
#endif
#if(MULTIPOLE_ORDER >= 4)
    float fac4;
#endif
#if(MULTIPOLE_ORDER >= 5)
    float fac5;
#endif
  } shortrange_factors[NTAB + 1];
#endif

 public:
#if defined(PMGRID)
  char DoPM;
#endif
  char DoEwald;

  struct gfactors
  {
    MyReal rinv2;
#if(MULTIPOLE_ORDER >= 2)
    MyReal rinv3;
#endif

    MyReal fac0;  //    g
    MyReal fac1;  //    g'
#if(MULTIPOLE_ORDER >= 2)
    MyReal fac2;  //    (g'' - g'/r) = (g'' - g'/r)
#endif
#if(MULTIPOLE_ORDER >= 3)
    MyReal fac3;  //    (g''' - 3*g''/r + 3 * g'/r^2) = (g''' - (3/r) * fac2)
#endif
#if(MULTIPOLE_ORDER >= 4)
    MyReal fac4;  //    (g'''' - 6*g'''/r + 15*g''/r^2 - 15*g'/r^3)  =  (g'''' - (6/r) * fac3 - (3/r^2) * fac2)
#endif
#if(MULTIPOLE_ORDER >= 5)
    MyReal fac5;  //    (g''''' - 10*g''''/r + 45*g'''/r^2 - 105*g''/r^3 +  105 g'/r^4)  =  (g''''' - (10/r) * fac4 - (15/r^2) * fac3)
#endif

    gfactors() /* constructor, initialize the factors to zero */
    {
      fac0 = 0;
      fac1 = 0;
#if(MULTIPOLE_ORDER >= 2)
      fac2 = 0;
#endif
#if(MULTIPOLE_ORDER >= 3)
      fac3 = 0;
#endif
#if(MULTIPOLE_ORDER >= 4)
      fac4 = 0;
#endif
#if(MULTIPOLE_ORDER >= 5)
      fac5 = 0;
#endif
    }
  };

  template <typename T>
  inline void get_gfactors_multipole(gfactors &res, const T r, const T h_max, const T rinv)
  {
    res.rinv2 = rinv * rinv;

#if(MULTIPOLE_ORDER >= 2)
    res.rinv3 = res.rinv2 * rinv;
#endif

    if(r >= h_max)
      {
#ifdef EVALPOTENTIAL
        res.fac0 += rinv;
#endif
        res.fac1 -= res.rinv2;
#if(MULTIPOLE_ORDER >= 2)
        res.fac2 += 3 * res.rinv3;
#endif
#if(MULTIPOLE_ORDER >= 3)
        res.fac3 -= 15 * res.rinv3 * rinv;
#endif
#if(MULTIPOLE_ORDER >= 4)
        res.fac4 += 105 * res.rinv3 * res.rinv2;
#endif
#if(MULTIPOLE_ORDER >= 5)
        res.fac5 -= 945 * res.rinv3 * res.rinv3;
#endif
      }
    else
      {
        T h_inv  = 1 / h_max;
        T h2_inv = h_inv * h_inv;
        T u      = r * h_inv;
        T fac1;
#if(MULTIPOLE_ORDER >= 2)
        T h3_inv = h_inv * h2_inv;
        T gpp;
#endif
#if(MULTIPOLE_ORDER >= 3)
        T gppp;
#endif
#if(MULTIPOLE_ORDER >= 4)
        T gpppp;
#endif
#if(MULTIPOLE_ORDER >= 5)
        T gppppp;
#endif
        if(u < static_cast<T>(0.5))
          {
            T u2 = u * u;
#ifdef EVALPOTENTIAL
            res.fac0 -= h_inv * (static_cast<T>(SOFTFAC4) +
                                 u2 * (static_cast<T>(SOFTFAC5) + u2 * (static_cast<T>(SOFTFAC6) * u + static_cast<T>(SOFTFAC7))));
#endif
            fac1 = -h2_inv * u * (static_cast<T>(SOFTFAC1) + u2 * (static_cast<T>(SOFTFAC2) * u + static_cast<T>(SOFTFAC3)));
            res.fac1 += fac1;
#if(MULTIPOLE_ORDER >= 2)
            gpp = -h3_inv * (static_cast<T>(SOFTFAC30) + (static_cast<T>(SOFTFAC31) + static_cast<T>(SOFTFAC32) * u) * u2);
#endif
#if(MULTIPOLE_ORDER >= 3)
            gppp = -h3_inv * h_inv * (static_cast<T>(SOFTFAC33) * u + static_cast<T>(SOFTFAC34) * u2);
#endif
#if(MULTIPOLE_ORDER >= 4)
            gpppp = -h3_inv * h2_inv * (static_cast<T>(SOFTFAC33) + static_cast<T>(SOFTFAC35) * u);
#endif
#if(MULTIPOLE_ORDER >= 5)
            gppppp = -h3_inv * h3_inv * static_cast<T>(SOFTFAC35);
#endif
          }
        else
          {
            T u2    = u * u;
            T u3    = u2 * u;
            T u3inv = 1 / u3;
#ifdef EVALPOTENTIAL
            res.fac0 -=
                h_inv * (static_cast<T>(SOFTFAC13) + static_cast<T>(SOFTFAC14) / u +
                         u2 * (static_cast<T>(SOFTFAC1) +
                               u * (static_cast<T>(SOFTFAC15) + u * (static_cast<T>(SOFTFAC16) + static_cast<T>(SOFTFAC17) * u))));
#endif
            fac1 = -h2_inv * u *
                   (static_cast<T>(SOFTFAC8) + static_cast<T>(SOFTFAC9) * u + static_cast<T>(SOFTFAC10) * u2 +
                    static_cast<T>(SOFTFAC11) * u3 + static_cast<T>(SOFTFAC12) * u3inv);
            res.fac1 += fac1;
#if(MULTIPOLE_ORDER >= 2)
            gpp = -h3_inv * (static_cast<T>(SOFTFAC40) + static_cast<T>(SOFTFAC41) / (u * u2) + static_cast<T>(SOFTFAC42) * u +
                             static_cast<T>(SOFTFAC43) * u2 + static_cast<T>(SOFTFAC44) * u2 * u);
#endif
#if(MULTIPOLE_ORDER >= 3)
            gppp = -h3_inv * h_inv *
                   (static_cast<T>(SOFTFAC45) + static_cast<T>(SOFTFAC46) / (u2 * u2) + static_cast<T>(SOFTFAC47) * u +
                    static_cast<T>(SOFTFAC48) * u2);
#endif
#if(MULTIPOLE_ORDER >= 4)
            gpppp =
                -h3_inv * h2_inv * (static_cast<T>(SOFTFAC49) / (u2 * u3) + static_cast<T>(SOFTFAC47) + static_cast<T>(SOFTFAC50) * u);
#endif
#if(MULTIPOLE_ORDER >= 5)
            gppppp = -h3_inv * h3_inv * (static_cast<T>(SOFTFAC51) / (u3 * u3) + static_cast<T>(SOFTFAC50));
#endif
          }
#if(MULTIPOLE_ORDER >= 2)
        T fac2 = (gpp - rinv * fac1);
        res.fac2 += fac2;
#endif
#if(MULTIPOLE_ORDER >= 3)
        T fac3 = (gppp - 3 * rinv * fac2);
        res.fac3 += fac3;
#endif
#if(MULTIPOLE_ORDER >= 4)
        T fac4 = (gpppp - 6 * rinv * fac3 - 3 * rinv * rinv * fac2);
        res.fac4 += fac4;
#endif
#if(MULTIPOLE_ORDER >= 5)
        T fac5 = (gppppp - 10 * rinv * fac4 - 15 * rinv * rinv * fac3);
        res.fac5 += fac5;
#endif
      }
  }

  template <typename T>
  inline void get_gfactors_monopole(gfactors &res, const T r, const T h_max, const T rinv)
  {
    if(r > h_max)
      {
        res.fac1 -= rinv * rinv;
#ifdef EVALPOTENTIAL
        res.fac0 += rinv;
#endif
      }
    else
      {
        T h_inv  = 1 / h_max;
        T h2_inv = h_inv * h_inv;
        T u      = r * h_inv;

        if(u < 0.5f)
          {
            T u2 = u * u;
            res.fac1 -= h2_inv * u * (static_cast<T>(SOFTFAC1) + u2 * (static_cast<T>(SOFTFAC2) * u + static_cast<T>(SOFTFAC3)));
#ifdef EVALPOTENTIAL
            res.fac0 -= h_inv * (static_cast<T>(SOFTFAC4) +
                                 u2 * (static_cast<T>(SOFTFAC5) + u2 * (static_cast<T>(SOFTFAC6) * u + static_cast<T>(SOFTFAC7))));
#endif
          }
        else
          {
            T u2 = u * u;
            T u3 = u2 * u;
            res.fac1 -= h2_inv * u *
                        (static_cast<T>(SOFTFAC8) + static_cast<T>(SOFTFAC9) * u + static_cast<T>(SOFTFAC10) * u2 +
                         static_cast<T>(SOFTFAC11) * u3 + static_cast<T>(SOFTFAC12) / u3);
#ifdef EVALPOTENTIAL
            res.fac0 -=
                h_inv * (static_cast<T>(SOFTFAC13) + static_cast<T>(SOFTFAC14) / u +
                         u2 * (static_cast<T>(SOFTFAC1) +
                               u * (static_cast<T>(SOFTFAC15) + u * (static_cast<T>(SOFTFAC16) + static_cast<T>(SOFTFAC17) * u))));
#endif
          }
      }
  }

  template <typename T>
  inline void get_gfactors_potential(gfactors &res, const T r, const T hmax, const T rinv)
  {
    if(r >= hmax)
      {
        res.fac0 = rinv;
#if(MULTIPOLE_ORDER >= 3)
        T rinv3  = rinv * rinv * rinv;
        res.fac1 = -rinv * rinv;
        res.fac2 = 3 * rinv3;
#endif
      }
    else
      {
        T h_inv = 1 / hmax;
#if(MULTIPOLE_ORDER >= 3)
        T h2_inv = h_inv * h_inv;
        T h3_inv = h2_inv * h_inv;
#endif
        T u = r * h_inv;
#if(MULTIPOLE_ORDER >= 3)
        T fac1;
#endif
        if(u < 0.5)
          {
            T u2 = u * u;
#if(MULTIPOLE_ORDER >= 3)
            fac1 = -h2_inv * u * (static_cast<T>(SOFTFAC1) + u2 * (static_cast<T>(SOFTFAC2) * u + static_cast<T>(SOFTFAC3)));
            res.fac1 += fac1;
#endif
            res.fac0 = -h_inv * (static_cast<T>(SOFTFAC4) +
                                 u2 * (static_cast<T>(SOFTFAC5) + u2 * (static_cast<T>(SOFTFAC6) * u + static_cast<T>(SOFTFAC7))));
          }
        else
          {
            T u2 = u * u;
#if(MULTIPOLE_ORDER >= 3)
            T u3 = u2 * u;
            fac1 = -h2_inv * u *
                   (static_cast<T>(SOFTFAC8) + static_cast<T>(SOFTFAC9) * u + static_cast<T>(SOFTFAC10) * u2 +
                    static_cast<T>(SOFTFAC11) * u3 + static_cast<T>(SOFTFAC12) / u3);
            res.fac1 += fac1;
#endif
            res.fac0 =
                -h_inv * (static_cast<T>(SOFTFAC13) + static_cast<T>(SOFTFAC14) / u +
                          u2 * (static_cast<T>(SOFTFAC1) +
                                u * (static_cast<T>(SOFTFAC15) + u * (static_cast<T>(SOFTFAC16) + static_cast<T>(SOFTFAC17) * u))));
          }

#if(MULTIPOLE_ORDER >= 3)
        T gpp;
        T u2 = u * u;

        if(u < 0.5)
          gpp = -h3_inv * (static_cast<T>(SOFTFAC30) + (static_cast<T>(SOFTFAC31) + static_cast<T>(SOFTFAC32) * u) * u2);
        else
          gpp = -h3_inv * (static_cast<T>(SOFTFAC40) + static_cast<T>(SOFTFAC41) / (u * u2) + static_cast<T>(SOFTFAC42) * u +
                           static_cast<T>(SOFTFAC43) * u2 + static_cast<T>(SOFTFAC44) * u2 * u);

        res.fac2 = (gpp - rinv * fac1);
#endif
      }
  }

#ifdef PMGRID
  inline bool modify_gfactors_pm_multipole(gfactors &res, const double r, const double rinv, const mesh_factors *mfp)
  {
    double tabentry = mfp->asmthfac * r;
    int tabindex    = (int)tabentry;

    if(tabindex < NTAB)
      {
        double w1 = tabentry - tabindex;
        double w0 = 1 - w1;

#ifdef EVALPOTENTIAL
        res.fac0 += mfp->asmthinv1 * (w0 * shortrange_factors[tabindex].fac0 + w1 * shortrange_factors[tabindex + 1].fac0);
#endif
        res.fac1 += mfp->asmthinv2 * (w0 * shortrange_factors[tabindex].fac1 + w1 * shortrange_factors[tabindex + 1].fac1);

#if(MULTIPOLE_ORDER >= 2)
        res.fac2 += mfp->asmthinv3 * (w0 * shortrange_factors[tabindex].fac2 + w1 * shortrange_factors[tabindex + 1].fac2);
#endif
#if(MULTIPOLE_ORDER >= 3)
        res.fac3 += mfp->asmthinv4 * (w0 * shortrange_factors[tabindex].fac3 + w1 * shortrange_factors[tabindex + 1].fac3);
#endif
#if(MULTIPOLE_ORDER >= 4)
        res.fac4 += mfp->asmthinv5 * (w0 * shortrange_factors[tabindex].fac4 + w1 * shortrange_factors[tabindex + 1].fac4);
#endif
#if(MULTIPOLE_ORDER >= 5)
        res.fac5 += mfp->asmthinv6 * (w0 * shortrange_factors[tabindex].fac5 + w1 * shortrange_factors[tabindex + 1].fac5);
#endif
        return false;
      }
    else
      return true;
  }

  inline bool modify_gfactors_pm_monopole(gfactors &res, const double r, const double rinv, const mesh_factors *mfp)
  {
    double tabentry = mfp->asmthfac * r;
    int tabindex    = (int)tabentry;

    if(tabindex < NTAB)
      {
        double w1 = tabentry - tabindex;
        double w0 = 1 - w1;

#ifdef EVALPOTENTIAL
        res.fac0 += mfp->asmthinv1 * (w0 * shortrange_factors[tabindex].fac0 + w1 * shortrange_factors[tabindex + 1].fac0);
#endif
        res.fac1 += mfp->asmthinv2 * (w0 * shortrange_factors[tabindex].fac1 + w1 * shortrange_factors[tabindex + 1].fac1);

        return false;
      }
    else
      return true;
  }

#endif
};

#endif
