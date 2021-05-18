/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_findlinkngb.cc
 *
 *  \brief find the nearest linking neighbors used in looking for saddle points
 */

#include "gadgetconfig.h"

#ifdef SUBFIND
#ifndef SUBFIND_HBT

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <cstdio>

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
#include "../mpi_utils/mpi_utils.h"
#include "../sort/cxxsort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

struct r2type
{
  MyFloat r2;
  int index;
};

static bool subfind_ngb_compare_dist(const r2type &a, const r2type &b) { return a.r2 < b.r2; }

static int *DM_NumNgb;
static double *Dist2list;
static MyFloat *Left, *Right;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
struct nearest_in : data_in_generic
{
  MyIntPosType IntPos[3];
  MyFloat DM_Hsml;
};

struct nearest_out
{
  int Ngb;
};

template <typename T_tree, typename T_domain, typename T_partset>
class nearest_comm : public generic_comm<nearest_in, nearest_out, T_tree, T_domain, T_partset>
{
 public:
  typedef generic_comm<nearest_in, nearest_out, T_tree, T_domain, T_partset> gcomm;
  using gcomm::D;
  using gcomm::Thread;
  using gcomm::Tp;  // This makes sure that we can access Tp from the base class without having to use "this->Tp"
  using gcomm::Tree;

  /* need to call the base class constructor explicitly */
  nearest_comm(T_domain *dptr, T_tree *tptr, T_partset *pptr) : gcomm(dptr, tptr, pptr) {}

  /* routine that fills the relevant particle/cell data into the input structure defined above */
  void particle2in(nearest_in *in, int i)
  {
    in->IntPos[0] = Tp->P[i].IntPos[0];
    in->IntPos[1] = Tp->P[i].IntPos[1];
    in->IntPos[2] = Tp->P[i].IntPos[2];
    in->DM_Hsml   = Tp->PS[i].v.DM_Hsml;
  }

  /* routine to store or combine result data */
  void out2particle(nearest_out *out, int i, int mode)
  {
    if(mode == MODE_LOCAL_PARTICLES) /* initial store */
      DM_NumNgb[i] = out->Ngb;
    else /* combine */
      DM_NumNgb[i] += out->Ngb;
  }

  /*! This function represents the core of the SPH density computation. The
   *  target particle may either be local, or reside in the communication
   *  buffer.
   */
  int evaluate(int target, int mode, int thread_id, int action, nearest_in *in, int numnodes, node_info *firstnode, nearest_out &out)
  {
    MyIntPosType *intpos = in->IntPos;
    double hsml          = in->DM_Hsml;
    int numngb           = 0;
    int exported         = 0;

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
                auto *P = Tree->get_Pp(no, shmrank);

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

                if(numngb >= Tp->NumPart)
                  Terminate("numngb >= Tp->NumPart");

                Dist2list[numngb++] = r2;
              }
            else if(no < Tree->MaxPart + Tree->MaxNodes) /* internal node */
              {
                if(mode == 1)
                  {
                    if(no < Tree->FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the
                                                           branch */
                      break;
                  }

                gravnode *current = Tree->get_nodep(no, shmrank);

                no      = current->sibling; /* in case the node can be discarded */
                shmrank = current->sibling_shmrank;

                double dxyz[3];
                Tp->nearest_image_intpos_to_pos(current->center.da, intpos,
                                                dxyz); /* converts the integer distance to floating point */

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
                if(dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2] > dist * dist)
                  continue;

                no      = current->nextnode; /* ok, we need to open the node */
                shmrank = current->nextnode_shmrank;
              }
            else
              {
                /* pseudo particle */

                if(mode == MODE_LOCAL_PARTICLES)
                  if(target >= 0) /* if no target is given, export will not occur */
                    {
                      exported = 1;

                      if(mode == MODE_LOCAL_PARTICLES)
                        Tree->tree_export_node_threads(no, target, &Thread);
                    }

                no = Tree->Nextnode[no - Tree->MaxNodes];
              }
          }
      }

    if(mode == MODE_LOCAL_PARTICLES) /* local particle */
      if(exported == 0)              /* completely local */
        if(numngb >= All.DesLinkNgb)
          {
            r2type *R2list = (r2type *)Mem.mymalloc("R2list", sizeof(r2type) * numngb);
            for(int i = 0; i < numngb; i++)
              {
                R2list[i].r2 = Dist2list[i];
              }

            mycxxsort(R2list, R2list + numngb, subfind_ngb_compare_dist);

            Tp->PS[target].v.DM_Hsml = sqrt(R2list[All.DesLinkNgb - 1].r2);
            numngb                   = All.DesLinkNgb;

            for(int i = 0; i < numngb; i++)
              {
                Dist2list[i] = R2list[i].r2;
              }

            Mem.myfree(R2list);
          }

    out.Ngb = numngb;

    return 0;
  }
};

template <typename partset>
void fof<partset>::subfind_find_linkngb(domain<partset> *SubDomain, int num, int *list)
{
  subfind_collective_printf("SUBFIND: root-task=%d: Start find_linkngb. (%d particles on root-task)\n", ThisTask, num);

  Dist2list = (double *)Mem.mymalloc("Dist2list", Tp->NumPart * sizeof(double));
  Left      = (MyFloat *)Mem.mymalloc("Left", sizeof(MyFloat) * Tp->NumPart);
  Right     = (MyFloat *)Mem.mymalloc("Right", sizeof(MyFloat) * Tp->NumPart);
  DM_NumNgb = (int *)Mem.mymalloc_movable(&DM_NumNgb, "DM_NumNgb", sizeof(int) * Tp->NumPart);

  int *targetlist = (int *)Mem.mymalloc("targetlist", num * sizeof(int));

  for(int idx = 0; idx < num; idx++)
    {
      targetlist[idx] = list[idx]; /* to preserve the input list, we make a copy */

      int i   = list[idx];
      Left[i] = Right[i] = 0;
    }

  nearest_comm<gravtree<partset>, domain<partset>, partset> commpattern{SubDomain, &FoFGravTree, Tp};

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  long long ntot;
  int iter = 0;

  do
    {
      double t0 = Logs.second();

      commpattern.execute(num, targetlist, MODE_DEFAULT);

      /* do final operations on results */
      int npleft = 0;
      for(int idx = 0; idx < num; idx++)
        {
          int i = targetlist[idx];

          /* now check whether we had enough neighbours */

          if(DM_NumNgb[i] != All.DesLinkNgb && ((Right[i] - Left[i]) > 1.0e-6 * Left[i] || Left[i] == 0 || Right[i] == 0))
            {
              /* need to redo this particle */
              targetlist[npleft++] = i;

              if(DM_NumNgb[i] < All.DesLinkNgb)
                Left[i] = std::max<double>(Tp->PS[i].v.DM_Hsml, Left[i]);
              else
                {
                  if(Right[i] != 0)
                    {
                      if(Tp->PS[i].v.DM_Hsml < Right[i])
                        Right[i] = Tp->PS[i].v.DM_Hsml;
                    }
                  else
                    Right[i] = Tp->PS[i].v.DM_Hsml;
                }

              if(iter >= MAXITER - 10)
                {
                  double pos[3];
                  Tp->intpos_to_pos(Tp->P[i].IntPos, pos);

                  printf("i=%d task=%d ID=%d DM_Hsml=%g Left=%g Right=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", i, ThisTask,
                         (int)Tp->P[i].ID.get(), Tp->PS[i].v.DM_Hsml, Left[i], Right[i], (double)(Right[i] - Left[i]), pos[0], pos[1],
                         pos[2]);
                  myflush(stdout);
                }

              if(Right[i] > 0 && Left[i] > 0)
                Tp->PS[i].v.DM_Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
              else
                {
                  if(Right[i] == 0 && Left[i] == 0)
                    Terminate("can't occur");

                  if(Right[i] == 0 && Left[i] > 0)
                    Tp->PS[i].v.DM_Hsml *= 1.26;

                  if(Right[i] > 0 && Left[i] == 0)
                    Tp->PS[i].v.DM_Hsml /= 1.26;
                }
            }
        }

      num = npleft;

      sumup_large_ints(1, &npleft, &ntot, SubComm);

      double t1 = Logs.second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            subfind_collective_printf(
                "SUBFIND: root-task=%d: find linkngb iteration %d, need to repeat for %lld particles. (took %g sec)\n", ThisTask, iter,
                ntot, Logs.timediff(t0, t1));

          if(iter > MAXITER)
            Terminate("failed to converge in neighbour iteration in density_findlinkngb()\n");
        }
    }
  while(ntot > 0);

  Mem.myfree(targetlist);
  Mem.myfree(DM_NumNgb);
  Mem.myfree(Right);
  Mem.myfree(Left);

  Mem.myfree(Dist2list);

  subfind_collective_printf("SUBFIND: root-task=%d: Done with find_linkngb\n", ThisTask);
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
#endif
