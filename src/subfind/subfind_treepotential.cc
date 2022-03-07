/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_treepotential.cc
 *
 *  \brief routines to compute the gravitational potential of the particles making up a group
 */

#include "gadgetconfig.h"

#ifdef SUBFIND

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fof/fof.h"
#include "../gravtree/gravtree.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/generic_comm.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/peano.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
struct potdata_in : data_in_generic
{
  MyIntPosType IntPos[3];

  unsigned char Type;
#if NSOFTCLASSES > 1
  unsigned char SofteningClass;
#endif
};

struct potdata_out
{
  MyFloat Potential;
};

template <typename T_tree, typename T_domain, typename T_partset>
class potdata_comm : public generic_comm<potdata_in, potdata_out, T_tree, T_domain, T_partset>
{
 public:
  typedef generic_comm<potdata_in, potdata_out, T_tree, T_domain, T_partset> gcomm;
  typedef gravtree<T_partset> gtree;
  using gcomm::D;
  using gcomm::Thread;
  using gcomm::Tp;  // This makes sure that we can access Tp from the base class without having to use "this->Tp"
  using gcomm::Tree;

  /* need to call the base class constructor explicitly */
  potdata_comm(T_domain *dptr, T_tree *tptr, T_partset *pptr) : gcomm(dptr, tptr, pptr) {}

  /* routine that fills the relevant particle/cell data into the input structure defined above */
  void particle2in(potdata_in *in, int i) override
  {
    for(int k = 0; k < 3; k++)
      in->IntPos[k] = Tp->P[i].IntPos[k];

    in->Type = Tp->P[i].getType();
#if NSOFTCLASSES > 1
    in->SofteningClass = Tp->P[i].getSofteningClass();
#endif
  }

  /* routine to store or combine result data */
  void out2particle(potdata_out *out, int i, int mode) override
  {
    if(mode == MODE_LOCAL_PARTICLES) /* initial store */
      Tp->PS[i].u.s.u.DM_Potential = out->Potential;
    else /* combine */
      Tp->PS[i].u.s.u.DM_Potential += out->Potential;
  }

  int evaluate(int target, int mode, int thread_id, int action, potdata_in *in, int numnodes, node_info *firstnode,
               potdata_out &out) override
  {
    gravnode *nop          = NULL;
    double pot             = 0.0;
    int ninteractions      = 0;
    MyIntPosType intpos[3] = {in->IntPos[0], in->IntPos[1], in->IntPos[2]};

    double theta2 = All.ErrTolTheta * All.ErrTolTheta;

#if NSOFTCLASSES > 1
    double hmax, h_i = All.ForceSoftening[in->SofteningClass];
#else
    double hmax = All.ForceSoftening[0];
#endif

    for(int k = 0; k < numnodes; k++)
      {
        int no;

        if(mode == 0)
          no = Tree->MaxPart; /* root node */
        else
          {
            no = firstnode[k].Node;

            if(numnodes > Tree->MaxNodes)
              Terminate("numnodes=%d Tree->MaxNodes=%d  Tree->NumNodes=%d no=%d", numnodes, Tree->MaxNodes, Tree->NumNodes, no);

            no = Tree->get_nodep(no)->nextnode; /* open it */
          }

        unsigned int shmrank = Tree->TreeSharedMem_ThisTask;

        while(no >= 0)
          {
            vector<double> dxyz;
            double r2, mass;

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
            char flag_Q2Tensor = 0;
#endif

            if(no < Tree->MaxPart) /* single particle */
              {
                auto P = Tree->get_Pp(no, shmrank);

                /* convert the integer distance to floating point */
                Tp->nearest_image_intpos_to_pos(P->IntPos, intpos, dxyz.da);

                r2 = dxyz.r2();

                mass = P->getMass();

#if NSOFTCLASSES > 1
                double h_j = All.ForceSoftening[Tp->P[no].getSofteningClass()];
                hmax       = (h_j > h_i) ? h_j : h_i;
#endif

                no = Tree->get_nextnodep(shmrank)[no]; /* note: here shmrank cannot change */
              }
            else if(no < Tree->MaxPart + Tree->MaxNodes) /* we have an internal node */
              {
                if(mode == 1)
                  {
                    if(no < Tree->FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the
                                                           branch */
                      {
                        no = -1;
                        continue;
                      }
                  }

                nop = Tree->get_nodep(no, shmrank);

                if(nop->level <= LEVEL_ALWAYS_OPEN)
                  {
                    /* we always open the root node (its full node length couldn't be stored in the integer type */
                    no      = nop->nextnode;
                    shmrank = nop->nextnode_shmrank;
                    continue;
                  }

                MyIntPosType halflen = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - nop->level);
                MyIntPosType intlen  = halflen << 1;

                {
                  /* check whether we lie very close to the cell, and if yes, open it */
                  MyIntPosType dist[3] = {nop->center[0] - intpos[0], nop->center[1] - intpos[1], nop->center[2] - intpos[2]};

                  dist[0] = (((MySignedIntPosType)dist[0]) >= 0) ? dist[0] : -dist[0];
                  dist[1] = (((MySignedIntPosType)dist[1]) >= 0) ? dist[1] : -dist[1];
                  dist[2] = (((MySignedIntPosType)dist[2]) >= 0) ? dist[2] : -dist[2];

                  if(dist[0] < intlen && dist[1] < intlen && dist[2] < intlen)
                    {
                      /* open cell */
                      no      = nop->nextnode;
                      shmrank = nop->nextnode_shmrank;
                      continue;
                    }
                }

                /* converts the integer distance to floating point */
                Tp->nearest_image_intpos_to_pos(nop->s.da, intpos, dxyz.da);

                r2 = dxyz.r2();

                mass = nop->mass;

                double len  = intlen * Tp->FacIntToCoord;
                double len2 = len * len;

                /* check Barnes-Hut opening criterion */
                if(len2 > r2 * theta2)
                  {
                    /* open cell */
                    no      = nop->nextnode;
                    shmrank = nop->nextnode_shmrank;
                    continue;
                  }

#if NSOFTCLASSES > 1
                double h_j = All.ForceSoftening[nop->maxsofttype];

                if(h_j > h_i)
                  {
                    if(r2 < h_j * h_j)
                      {
                        if(All.ForceSoftening[nop->minsofttype] < All.ForceSoftening[nop->maxsofttype])
                          {
                            /* open cell */
                            no      = nop->nextnode;
                            shmrank = nop->nextnode_shmrank;
                            continue;
                          }
                      }
                    hmax = h_j;
                  }
                else
                  hmax = h_i;
#endif

                  /* ok, node can be used */

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
                /* will need to account for quadrupole tensor of node */
                flag_Q2Tensor = 1;
#endif

                no      = nop->sibling;
                shmrank = nop->sibling_shmrank;
              }
            else if(no >= Tree->ImportedNodeOffset) /* point from imported nodelist */
              {
                Terminate("We don't expect TreePoints here");
              }
            else /* pseudo particle */
              {
                if(mode == MODE_LOCAL_PARTICLES)
                  Tree->tree_export_node_threads(no, target, &Thread);

                no = Tree->Nextnode[no - Tree->MaxNodes];
                continue;
              }

            /* now evaluate the multipole moment */
            if(mass)
              {
                double r = sqrt(r2);

                double rinv = (r > 0) ? 1.0 / r : 0;

                typename gtree::gfactors gfac;

                Tree->get_gfactors_potential(gfac, r, hmax, rinv);

                pot -= mass * gfac.fac0;

                ninteractions++;

#if(MULTIPOLE_ORDER >= 3) || (MULTIPOLE_ORDER >= 2 && defined(EXTRAPOTTERM))
                if(flag_Q2Tensor)
                  {
                    double g1               = gfac.fac1 * rinv;
                    double g2               = gfac.fac2 * rinv * rinv;
                    vector<MyDouble> Q2dxyz = nop->Q2Tensor * dxyz;
                    double Q2dxyz2          = Q2dxyz * dxyz;
                    double Q2trace          = nop->Q2Tensor.trace();

                    pot -= 0.5 * (g1 * Q2trace + g2 * Q2dxyz2);  //  quadrupole potential
                  }
#endif
              }
          }
      }

    /* now store the result */

    out.Potential = pot;

    return ninteractions;
  }
};

template <typename partset>
void fof<partset>::subfind_potential_compute(domain<partset> *SubDomain, int num, int *darg)
{
  /* create an object for handling the communication */
  potdata_comm<gravtree<partset>, domain<partset>, partset> commpattern{SubDomain, &FoFGravTree, Tp};

  commpattern.execute(num, darg, MODE_DEFAULT);
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
