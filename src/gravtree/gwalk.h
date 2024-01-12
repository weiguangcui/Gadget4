/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  gwalk.h
 *
 *  \brief defines a class for walking the gravitational tree
 */

#ifndef GRAVTREE_WALK_H
#define GRAVTREE_WALK_H

#include "gadgetconfig.h"

#include "../mpi_utils/shared_mem_handler.h"

class gwalk : public gravtree<simparticles>
{
 public:
  void gravity_tree(int timebin);

 private:
  long long interactioncountPP;
  long long interactioncountPN;

  MyReal theta2;
  MyReal thetamax2;
  MyReal errTolForceAcc;

#ifdef PRESERVE_SHMEM_BINARY_INVARIANCE
  bool skip_actual_force_computation;
#endif

  struct pinfo
  {
    MyIntPosType *intpos;
    MyReal aold;
    MyReal h_i;
    int Type;
#if NSOFTCLASSES > 1
    int SofteningClass;
#endif
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
    int InsideOutsideFlag;
#endif

    vector<MyFloat> *acc;
    MyFloat *pot;
    int *GravCost;
  };

  inline int get_pinfo(int i, pinfo &pdat)
  {
    int ptype;

    if(i < Tp->NumPart)
      {
        ptype = NODE_TYPE_LOCAL_PARTICLE;

        pdat.intpos = Tp->P[i].IntPos;

        pdat.Type = Tp->P[i].getType();
#if NSOFTCLASSES > 1
        pdat.SofteningClass = Tp->P[i].getSofteningClass();
#endif
        pdat.aold = Tp->P[i].OldAcc;
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
        pdat.InsideOutsideFlag = Tp->P[i].InsideOutsideFlag;
#endif

        pdat.acc = &Tp->P[i].GravAccel;
#ifdef EVALPOTENTIAL
        pdat.pot = &Tp->P[i].Potential;
#endif
        pdat.GravCost = &Tp->P[i].GravCost;
      }
    else
      {
        ptype = NODE_TYPE_TREEPOINT_PARTICLE;

        int n = i - ImportedNodeOffset;

        pdat.intpos = Points[n].IntPos;

        pdat.Type = Points[n].Type;
#if NSOFTCLASSES > 1
        pdat.SofteningClass = Points[n].SofteningClass;
#endif
        pdat.aold = Points[n].OldAcc;
#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
        pdat.InsideOutsideFlag = Points[n].InsideOutsideFlag;
#endif

        int idx  = ResultIndexList[n];
        pdat.acc = &ResultsActiveImported[idx].GravAccel;
#ifdef EVALPOTENTIAL
        pdat.pot = &ResultsActiveImported[idx].Potential;
#endif
        pdat.GravCost = &ResultsActiveImported[idx].GravCost;
      }

#if NSOFTCLASSES > 1
    pdat.h_i = All.ForceSoftening[pdat.SofteningClass];
#else
    pdat.h_i = All.ForceSoftening[0];
#endif

    return ptype;
  }

  inline void gwalk_open_node(const pinfo &pdat, int i, char ptype, gravnode *nop, int mintopleafnode, int committed);
  void gravity_force_interact(const pinfo &pdat, int i, int no, char ptype, char no_type, unsigned char shmrank, int mintopleafnode,
                              int committed);

  inline int evaluate_particle_node_opening_criterion_and_interaction(const pinfo &pdat, gravnode *nop);
  inline void evaluate_particle_particle_interaction(const pinfo &pdat, const int no, const char jtype, int no_task);
};

#endif
