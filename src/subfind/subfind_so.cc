/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_so.cc
 *
 *  \brief spherical overdensity virial radius determination and property calculation
 */

#include "gadgetconfig.h"

#ifdef SUBFIND

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fof/fof.h"
#include "../gravtree/gravtree.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/generic_comm.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

static double *R200, *M200;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
struct sodata_in : data_in_generic
{
  MyIntPosType IntPos[3];
  double R200;
};

struct sodata_out
{
  double Mass;
};

template <typename T_tree, typename T_domain, typename T_partset>
class sodata_comm : public generic_comm<sodata_in, sodata_out, T_tree, T_domain, T_partset>
{
 public:
  typedef generic_comm<sodata_in, sodata_out, T_tree, T_domain, T_partset> gcomm;
  using gcomm::D;
  using gcomm::Thread;
  using gcomm::Tp;  // This makes sure that we can access Tp from the base class without having to use "this->Tp"
  using gcomm::Tree;

  typename fof<T_partset>::group_properties *Group;

  /* need to call the base class constructor explicitly */
  sodata_comm(T_domain *dptr, T_tree *tptr, T_partset *pptr, typename fof<T_partset>::group_properties *grpptr)
      : gcomm(dptr, tptr, pptr)
  {
    Group = grpptr;
  }

  /* routine that fills the relevant particle/cell data into the input structure defined above */
  void particle2in(sodata_in *in, int i) override
  {
    in->IntPos[0] = Group[i].IntPos[0];
    in->IntPos[1] = Group[i].IntPos[1];
    in->IntPos[2] = Group[i].IntPos[2];

    in->R200 = R200[i];
  }

  /* routine to store or combine result data */
  void out2particle(sodata_out *out, int i, int mode) override
  {
    if(mode == MODE_LOCAL_PARTICLES) /* initial store */
      M200[i] = out->Mass;
    else /* combine */
      M200[i] += out->Mass;
  }

  /*! This function represents the core of the SPH density computation. The
   *  target particle may either be local, or reside in the communication
   *  buffer.
   */
  int evaluate(int target, int mode, int thread_id, int action, sodata_in *in, int numnodes, node_info *firstnode,
               sodata_out &out) override
  {
    MyIntPosType *intpos = in->IntPos;
    double hsml          = in->R200;

    double mass = 0;

    for(int k = 0; k < numnodes; k++)
      {
        int no;

        if(mode == MODE_LOCAL_PARTICLES)
          {
            no = Tree->MaxPart; /* root node */
          }
        else
          {
            no = firstnode[k].Node;
            no = Tree->get_nodep(no)->nextnode; /* open it */
          }

        unsigned int shmrank = Tree->TreeSharedMem_ThisTask;

        while(no >= 0)
          {
            if(no < Tree->MaxPart) /* single particle */
              {
                auto P = Tree->get_Pp(no, shmrank);

                no = Tree->get_nextnodep(shmrank)[no]; /* note: here shmrank cannot change */

                double dxyz[3];
                Tp->nearest_image_intpos_to_pos(P->IntPos, intpos, dxyz); /* converts the integer distance to floating point */

                double h2 = hsml * hsml;

                double r2 = dxyz[0] * dxyz[0];
                if(r2 > h2)
                  continue;

                r2 += dxyz[1] * dxyz[1];
                if(r2 > h2)
                  continue;

                r2 += dxyz[2] * dxyz[2];
                if(r2 > h2)
                  continue;

                mass += P->getMass();
              }
            else if(no < Tree->MaxPart + Tree->MaxNodes) /* internal node */
              {
                if(mode == MODE_IMPORTED_PARTICLES)
                  {
                    /* we reached a top-level node again, which means that we are done with the branch */
                    if(no < Tree->FirstNonTopLevelNode)
                      break;
                  }

                gravnode *current = Tree->get_nodep(no, shmrank);

                no      = current->sibling; /* in case the node can be discarded */
                shmrank = current->sibling_shmrank;

                /* converts the integer distance to floating point */
                double dxyz[3];
                Tp->nearest_image_intpos_to_pos(current->center.da, intpos, dxyz);

                double lenhalf = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - 1 - current->level)) * Tp->FacIntToCoord;

                double dist = hsml + lenhalf;

                if(fabs(dxyz[0]) > dist)
                  continue;
                if(fabs(dxyz[1]) > dist)
                  continue;
                if(fabs(dxyz[2]) > dist)
                  continue;

                /* now test against the minimal sphere enclosing everything */
                dist += FACT1 * 2.0 * lenhalf;

                double r2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];

                if(r2 > dist * dist)
                  continue;

                if(no >= Tree->FirstNonTopLevelNode) /* only do this for fully local nodes */
                  {
                    double lenhalf = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - 1 - current->level)) * Tp->FacIntToCoord;

                    /* test whether the node is contained within the sphere */
                    dist = hsml - FACTSQRT3 * lenhalf;
                    if(dist > 0)
                      if(r2 < dist * dist)
                        {
                          mass += current->mass;
                          continue;
                        }
                  }

                no      = current->nextnode; /* ok, we need to open the node */
                shmrank = current->nextnode_shmrank;
              }
            else if(no >= Tree->ImportedNodeOffset) /* point from imported nodelist */
              {
                int n = no - Tree->ImportedNodeOffset;
                no    = Tree->Nextnode[no - Tree->MaxNodes];
                /* note: here shmrank cannot change */

                double dxyz[3];
                Tp->nearest_image_intpos_to_pos(Tree->Points[n].IntPos, intpos,
                                                dxyz); /* converts the integer distance to floating point */

                double h2 = hsml * hsml;

                double r2 = dxyz[0] * dxyz[0];
                if(r2 > h2)
                  continue;

                r2 += dxyz[1] * dxyz[1];
                if(r2 > h2)
                  continue;

                r2 += dxyz[2] * dxyz[2];
                if(r2 > h2)
                  continue;

                mass += Tree->Points[n].Mass;
              }
            else /* pseudo particle */
              {
                if(mode == MODE_LOCAL_PARTICLES)
                  Tree->tree_export_node_threads(no, target, &Thread);

                no = Tree->Nextnode[no - Tree->MaxNodes];
              }
          }
      }

    out.Mass = mass;

    return 0;
  }
};

template <typename partset>
double fof<partset>::subfind_get_overdensity_value(int type, double ascale)
{
  double z = 1 / ascale - 1;
  double omegaz =
      All.Omega0 * pow(1 + z, 3) / (All.Omega0 * pow(1 + z, 3) + (1 - All.Omega0 - All.OmegaLambda) * pow(1 + z, 2) + All.OmegaLambda);
  double x = omegaz - 1;

  if(type == 0)
    {
      return 200.0;  // Mean200
    }
  else if(type == 1)
    {  // Generalized Tophat overdensity
      return (18 * M_PI * M_PI + 82 * x - 39 * x * x) / omegaz;
    }
  else if(type == 2)
    {
      return 200.0 / omegaz;  // DeltaCrit200
    }
  else if(type == 3)
    {
      return 500.0 / omegaz;  // DeltaCrit500
    }
  else
    Terminate("can't be");

  return 0;
}

template <typename partset>
double fof<partset>::subfind_overdensity(void)
{
  int *TargetList = (int *)Mem.mymalloc("TargetList", Ngroups * sizeof(int));

  double *Left  = (double *)Mem.mymalloc("Left", sizeof(double) * Ngroups);
  double *Right = (double *)Mem.mymalloc("Right", sizeof(double) * Ngroups);

  R200 = (double *)Mem.mymalloc("R200", sizeof(double) * Ngroups);
  M200 = (double *)Mem.mymalloc("M200", sizeof(double) * Ngroups);

  double rhoback = 3 * All.Omega0 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  double tstart = Logs.second();

  sodata_comm<gravtree<partset>, domain<partset>, partset> commpattern{FoFDomain, &FoFGravTree, Tp, Group};

  for(int rep = 0; rep < 4; rep++) /* repeat for all four overdensity values */
    {
      int Nso = 0;

      for(int i = 0; i < Ngroups; i++)
        {
          if(Group[i].Nsubs > 0)
            {
              double rguess = pow(All.G * Group[i].Mass / (100 * All.Hubble * All.Hubble), 1.0 / 3);

              TargetList[Nso++] = i;

              Right[i] = 3 * rguess;
              Left[i]  = 0;

              R200[i] = 0.5 * (Left[i] + Right[i]);
            }
        }

      int iter = 0;
      long long ntot;

      /* we will repeat the whole thing for those groups where we didn't converge to a SO radius yet */
      do
        {
          double t0 = Logs.second();

          commpattern.execute(Nso, TargetList, MODE_DEFAULT);

          /* do final operations on results */
          int npleft = 0;
          for(int n = 0; n < Nso; n++)
            {
              int i = TargetList[n];

              double overdensity = M200[i] / (4.0 * M_PI / 3.0 * R200[i] * R200[i] * R200[i]) / rhoback;

              if((Right[i] - Left[i]) > 1.0e-4 * Left[i])
                {
                  /* need to redo this group */
                  TargetList[npleft++] = i;

                  double delta = subfind_get_overdensity_value(rep, Group[i].Ascale);
                  if(overdensity > delta)
                    Left[i] = R200[i];
                  else
                    Right[i] = R200[i];

                  R200[i] = 0.5 * (Left[i] + Right[i]);

                  if(iter >= MAXITER - 10)
                    {
                      printf("gr=%d task=%d  R200=%g Left=%g Right=%g Menclosed=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", i, ThisTask,
                             R200[i], Left[i], Right[i], M200[i], Right[i] - Left[i], Group[i].Pos[0], Group[i].Pos[1],
                             Group[i].Pos[2]);
                      myflush(stdout);
                    }
                }
            }

          Nso = npleft;

          sumup_large_ints(1, &npleft, &ntot, Communicator);

          double t1 = Logs.second();

          if(ntot > 0)
            {
              iter++;

              if(iter > 0)
                mpi_printf("SUBFIND: SO iteration %2d: need to repeat for %12lld halo centers. (took %g sec)\n", iter, ntot,
                           Logs.timediff(t0, t1));

              if(iter > MAXITER)
                Terminate("failed to converge in SO iteration");
            }
        }
      while(ntot > 0);

      for(int i = 0; i < Ngroups; i++)
        {
          if(Group[i].Nsubs > 0)
            {
              double overdensity = M200[i] / (4.0 * M_PI / 3.0 * R200[i] * R200[i] * R200[i]) / rhoback;

              double delta = subfind_get_overdensity_value(rep, Group[i].Ascale);

              if((overdensity - delta) > 0.1 * delta)
                {
                  R200[i] = M200[i] = 0;
                }
              else if(M200[i] < 5 * Group[i].Mass / Group[i].Len)
                {
                  R200[i] = M200[i] = 0;
                }
            }
          else
            R200[i] = M200[i] = 0;

          switch(rep)
            {
              case 0:
                Group[i].M_Mean200 = M200[i];
                Group[i].R_Mean200 = R200[i];
                break;
              case 1:
                Group[i].M_TopHat200 = M200[i];
                Group[i].R_TopHat200 = R200[i];
                break;
              case 2:
                Group[i].M_Crit200 = M200[i];
                Group[i].R_Crit200 = R200[i];
                break;
              case 3:
                Group[i].M_Crit500 = M200[i];
                Group[i].R_Crit500 = R200[i];
                break;
            }
        }
    }

  Mem.myfree(M200);
  Mem.myfree(R200);
  Mem.myfree(Right);
  Mem.myfree(Left);
  Mem.myfree(TargetList);

  double tend = Logs.second();
  return Logs.timediff(tstart, tend);
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
