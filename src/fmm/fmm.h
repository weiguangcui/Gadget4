/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file fmm.h
 *
 *  \brief declares the class used for the fast multipole method (FMM)
 */

#ifndef FMM_H
#define FMM_H

#include "gadgetconfig.h"

#ifdef FMM

#include "../data/symtensors.h"
#include "../gravtree/gravtree.h"

class fmm : public gravtree<simparticles>
{
 public:
  void gravity_fmm(int timebin);

 private:
  long long interactioncountPP;
  long long interactioncountPN;
  long long interactioncountNN;
  long long sum_NumForeignNodes;
  long long sum_NumForeignPoints;
  long long interactioncountEffective;

  char *Topnode_depends_on_local_mass;

  MyReal errTolTheta2;
  MyReal errTolThetaMax2;
  MyReal errTolForceAcc;

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  bool skip_actual_force_computation;
#endif

  struct fieldcoeff
  {
    MyReal phi;
    vector<MyReal> dphi; /* first potential derivative at center of mass (i.e. minus the acceleration) */
#if(MULTIPOLE_ORDER >= 2)
    symtensor2<MyReal> d2phi; /* second derivatives  */
#endif
#if(MULTIPOLE_ORDER >= 3)
    symtensor3<MyReal> d3phi; /* third derivatives */
#endif
#if(MULTIPOLE_ORDER >= 4)
    symtensor4<MyReal> d4phi; /* fourth derivatives */
#endif
#if(MULTIPOLE_ORDER >= 5)
    symtensor5<MyReal> d5phi; /* fifth derivatives */
#endif

    MyReal interactions;
  };

  struct taylor_data
  {
    fieldcoeff coeff;
  };
  taylor_data *TaylorCoeff;

  struct fmm_workstack_data
  {
    int Node1;
    int Node2;
    unsigned char ShmRank1;
    unsigned char ShmRank2;
    int MinTopLeafNode;
  };

  fmm_workstack_data *FMM_WorkStack;

  static bool compare_fmm_workstack(const fmm_workstack_data &a, const fmm_workstack_data &b)
  {
    if(a.MinTopLeafNode < b.MinTopLeafNode)
      return true;
    if(a.MinTopLeafNode > b.MinTopLeafNode)
      return false;

    return a.Node1 < b.Node1;
  }

  inline void fmm_add_to_work_stack(int node1, int node2, unsigned char shmrank1, unsigned char shmrank2, int mintopleafnode)
  {
    if(NumOnWorkStack + NewOnWorkStack >= MaxOnWorkStack)
      {
        Terminate("we have run out of space:  NumOnWorkStack=%d + NewOnWorkStack=%d >= MaxOnWorkStack=%d", NumOnWorkStack,
                  NewOnWorkStack, MaxOnWorkStack);
      }

    FMM_WorkStack[NumOnWorkStack + NewOnWorkStack].Node1          = node1;
    FMM_WorkStack[NumOnWorkStack + NewOnWorkStack].Node2          = node2;
    FMM_WorkStack[NumOnWorkStack + NewOnWorkStack].ShmRank1       = shmrank1;
    FMM_WorkStack[NumOnWorkStack + NewOnWorkStack].ShmRank2       = shmrank2;
    FMM_WorkStack[NumOnWorkStack + NewOnWorkStack].MinTopLeafNode = mintopleafnode;

    NewOnWorkStack++;
  }

  inline bool fmm_depends_on_local_mass(int no, unsigned char shmrank)
  {
    if(no >= MaxPart && no < FirstNonTopLevelNode)
      {
        if(Topnode_depends_on_local_mass[no])
          return true;
        else
          return false;
      }
    else if(no >= FirstNonTopLevelNode && no < MaxPart + MaxNodes)
      {
        if(shmrank == Shmem.Island_ThisTask)
          return true;
        else
          return false;
      }
    else
      return false;
  }

  void fmm_force_interact(int no_sink, int no_source, char type_sink, char type_source, unsigned char shmrank_sink,
                          unsigned char shmrank_source, int mintopleafnode, int committed);

  void fmm_force_passdown(int no, unsigned char shmrank, taylor_data taylor_current);
  void fmm_open_both(gravnode *noptr_sink, gravnode *noptr_source, int mintopleafnode, int committed);

  void fmm_open_node(int no_particle, gravnode *nop, char type_particle, unsigned char shmrank_particle, int mintopleafnode,
                     int committed);
  void fmm_particle_particle_interaction(int no_sink, int no_source, int type_sink, int type_source, unsigned char shmrank_sink,
                                         unsigned char shmrank_source);

  void fmm_particle_node_interaction(int no_sink, int no_source, int type_sink, int type_source, unsigned char shmrank_sink,
                                     unsigned char shmrank_source, gravnode *noptr_source, vector<MyReal> &dxyz, MyReal &r2);

  void fmm_node_node_interaction(int no_sink, int no_source, int type_sink, int type_source, unsigned char shmrank_sink,
                                 unsigned char shmrank_source, gravnode *noptr_sink, gravnode *noptr_source, vector<MyReal> &dxyz,
                                 MyReal &r2);

  int fmm_evaluate_node_node_opening_criterion(gravnode *noptr_sink, gravnode *noptr_source, vector<MyReal> &dxyz, MyReal &r2);

  int fmm_evaluate_particle_node_opening_criterion(int no_sink, char type_sink, unsigned char shmrank_sink, gravnode *nop_source,
                                                   vector<MyReal> &dxyz, MyReal &r2);

  void fmm_determine_nodes_with_local_mass(int no, int sib);
};

#endif
#endif
