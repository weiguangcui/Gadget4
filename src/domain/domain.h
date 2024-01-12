/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file domain.h
 *
 *  \brief declares the class used for the domain decomposition
 */
#ifndef DOMAIN_H
#define DOMAIN_H

#include "gadgetconfig.h"

#ifndef ALLVARS_H
#include "../data/allvars.h"
#endif
#include "../data/dtypes.h"
#include "../mpi_utils/setcomm.h"

enum domain_options
{
  STANDARD,
  COLL_SUBFIND,
  SERIAL_SUBFIND
};

template <typename partset> /* partset will either be the 'simparticles' or the 'lightconeparticles' class is the matching
                               particle_data struct */
class domain : public setcomm
{
 private:
  partset *Tp;

 public:
  domain(MPI_Comm comm, partset *Tp_ptr) : setcomm(comm) /* constructor */
  {
    Tp = Tp_ptr;

    ListOfTopleaves    = NULL;
    TaskOfLeaf         = NULL;
    TopNodes           = NULL;
    NumTopleafOfTask   = NULL;
    FirstTopleafOfTask = NULL;
  }

  typedef typename partset::pdata pdata;

 public:
  domain_options Mode;

  int NTopnodes;
  int NTopleaves;

  int MultipleDomains = 1;

  int MaxTopNodes; /**< Maximum number of nodes in the top-level tree used for domain decomposition */

  int *ListOfTopleaves;
  int *NumTopleafOfTask;
  int *FirstTopleafOfTask;

  int *TaskOfLeaf;

  /** The top node structure is an octree used for encoding the domain
      decomposition. Its leaf nodes are the units into which the domain
      is decomposed. */
  struct topnode_data
  {
    int Daughter; /*!< index of first daughter cell (out of 8) of top-level node */
    int Leaf;     /*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
  };

  topnode_data *TopNodes;

  /** Array of task numbers holding the respective top-level nodes. For
      the topnodes entries, it is indexed by the Leaf member, for
      pseudoparticles it is indexed by the node
      number-MaxPart-MaxNodes.  */

  void domain_decomposition(domain_options mode);
  void domain_allocate(void);
  void domain_allocate(int maxtopnodes);
  void domain_free(void);
  void domain_resize_storage(int count_get, int count_get_sph, int option_flag);

  size_t domain_sizeof_topnode_data(void) { return sizeof(topnode_data); }

 private:
  struct local_topnode_data
  {
    double Cost;
    long long CountTot; /*!< counts the global number of particles in this top-level node */
    peanokey StartKey;  /*!< first Peano-Hilbert key in top-level node */
    int Count;          /*!< counts the local number of particles in this top-level node */
    int Level;          /*!< encodes side-length of the node, level 0 is the root node */
    int PIndex;         /*!< first particle in node */
  };

  struct domain_peano_hilbert_data
  {
    peanokey key;
    float cost;
  } * mp;

  struct domain_cost_data
  {
    double Cost;
  };

  struct domain_count_data
  {
    int task;
    int count;
    int origintask;
  };

  struct domain_segments_data
  {
    int task, start, end, used;
    double bin_GravCost[TIMEBINS];
    double bin_HydroCost[TIMEBINS];
    double load;
    double loadsph;
  };

  // domain_segments_data *domainAssign;

  domain_cost_data *domain_leaf_cost;

  local_topnode_data *topNodes;

  peanokey *domain_key;

  int NumTimeBinsToBeBalanced;
  int ListOfTimeBinsToBeBalanced[TIMEBINS];
  double GravCostPerListedTimeBin[TIMEBINS];
  double MaxGravCostPerListedTimeBin[TIMEBINS];
  double GravCostNormFactors[TIMEBINS];
  double HydroCostPerListedTimeBin[TIMEBINS];
  double HydroCostNormFactors[TIMEBINS];
  double NormFactorLoad;
  double NormFactorLoadSph;
  double TotalCost;

  int domain_grav_weight[TIMEBINS];
  int domain_hydro_weight[TIMEBINS];
  int domain_to_be_balanced[TIMEBINS];

  void domain_find_total_cost(void);
  void domain_report_balance(void);
  void domain_combine_multipledomains(void);
  void domain_init_sum_cost(void);
  void domain_countToGo(int *toGoDM, int *toGoSph);
  void domain_exchange(void);
  void domain_determineTopTree(void);

  void domain_printf(char *buf);
  void domain_rearrange_particle_sequence(void);

  void do_box_wrapping(void);
  void domain_coll_subfind_prepare_exchange(void);
  void domain_do_local_refine(int n, int *list);
  void domain_walktoptree(int no);
  double domain_get_cost_summed_over_timebins(int i);

  /* note: static here needed to suppress the hidden 'this' argument so that we can use this in a call of sort */
  static bool domain_compare_count(const domain_count_data &a, const domain_count_data &b) { return a.count > b.count; }

  static bool domain_compare_key(const domain_peano_hilbert_data &a, const domain_peano_hilbert_data &b) { return a.key < b.key; }

  static bool domain_sort_candidates(const int &a, const int &b) { return a < b; }

  struct cost_queue_data
  {
    double value;
#ifdef SIMPLE_DOMAIN_AGGREGATION
    double aggregated_value;
#endif
    int index;
  };

  static bool domain_sort_cost_queue_data(const cost_queue_data &a, const cost_queue_data &b)
  {
#ifdef SIMPLE_DOMAIN_AGGREGATION
    if(a.aggregated_value > b.aggregated_value)
      return true;
    else if(a.aggregated_value < b.aggregated_value)
      return false;
#endif
    return a.value > b.value;
  }

  void domain_determinate_aggregated_value(cost_queue_data *data, int ndomains);

  inline int n_to_no(int n)
  {
    int no       = 0;
    peanokey key = domain_key[n];

    while(TopNodes[no].Daughter >= 0)
      {
        unsigned int off = ((key.hs & (~((~((MyIntPosType)0)) >> 3))) >> (BITS_FOR_POSITIONS - 3));

        no = TopNodes[no].Daughter + off;

        key.hs <<= 3;
        key.hs |= (key.is & (~((~((MyIntPosType)0)) >> 3))) >> (BITS_FOR_POSITIONS - 3);

        key.is <<= 3;
        key.is |= (key.ls & (~((~((MyIntPosType)0)) >> 3))) >> (BITS_FOR_POSITIONS - 3);

        key.ls <<= 3;
      }

    return TopNodes[no].Leaf;
  }

  struct peano_hilbert_data
  {
    peanokey key;
    int index;
  };

  static bool compare_peano_hilbert_data(const peano_hilbert_data &a, const peano_hilbert_data &b) { return a.key < b.key; }

  void peano_hilbert_order(peanokey *key);

  struct balance_try_data
  {
    int nextra;
    double try_balance;
  };

  static bool domain_compare_trybalance(const balance_try_data &a, const balance_try_data &b) { return a.try_balance < b.try_balance; }

 public:
  void reorder_particles(int *Id, int Nstart, int N);
  void reorder_gas(int *Id);
  void reorder_PS(int *Id, int Nstart, int N);
  void reorder_P_and_PS(int *Id);
  void reorder_P_PS(int NumGas, int NumPart);

  void particle_exchange_based_on_PS(MPI_Comm Communicator);

 private:
  struct local_sort_data
  {
    int targetindex;
    int index;
  };

  static inline bool compare_local_sort_data_targetindex(const local_sort_data &a, const local_sort_data &b)
  {
    return a.targetindex < b.targetindex;
  }

#if defined(RANDOMIZE_DOMAINCENTER_TYPES) || defined(RANDOMIZE_DOMAINCENTER)
  MyIntPosType domainInnersize;
  MyIntPosType domainReferenceIntPos[3];
  MySignedIntPosType domainXmintot[3], domainXmaxtot[3];

  void domain_find_type_extension(void);
  int domain_type_extension_overlap(int j);
#endif

#ifdef DOMAIN_SPECIAL_CHECK
  void domain_special_check(int mode, int ndomains);
#endif

  void domain_printf(const char *fmt, ...)
  {
    if((Mode == STANDARD && ThisTask == 0))  // || (Mode == COLL_SUBFIND))
      {
        va_list l;
        va_start(l, fmt);
        vprintf(fmt, l);
        //        myflush(stdout);
        va_end(l);
      }
  }
};

#endif
