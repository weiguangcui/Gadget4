/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file fof_nearest.cc
 *
 *  \brief routines to find nearest neighbors
 */

#include "gadgetconfig.h"

#ifdef FOF

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
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/generic_comm.h"
#include "../ngbtree/ngbtree.h"
#include "../system/system.h"

/* local data structure for collecting particle/cell data that is sent to other processors if needed */

struct fofdata_in : data_in_generic
{
  MyIntPosType IntPos[3];
  MyFloat Hsml;
};

/* local data structure that holds results acquired on remote processors */
struct fofdata_out
{
  MyFloat Distance;
  MyIDStorage MinID;
  int MinIDTask;
#if defined(SUBFIND)
  MyFloat DM_Hsml;
#endif
};

template <typename T_tree, typename T_domain, typename T_partset>
class fofdata_comm : public generic_comm<fofdata_in, fofdata_out, T_tree, T_domain, T_partset>
{
 public:
  typedef generic_comm<fofdata_in, fofdata_out, T_tree, T_domain, T_partset> gcomm;
  using gcomm::D;
  using gcomm::Thread;
  using gcomm::Tp;  // This makes sure that we can access Tp from the base class without having to use "this->Tp"
  using gcomm::Tree;

  /* need to call the base class constructor explicitly */
  fofdata_comm(T_domain *dptr, T_tree *tptr, T_partset *pptr)
      : generic_comm<fofdata_in, fofdata_out, T_tree, T_domain, T_partset>(dptr, tptr, pptr)
  {
  }

  /* routine that fills the relevant particle/cell data into the input structure defined above */
  void particle2in(fofdata_in *in, int i)
  {
    for(int k = 0; k < 3; k++)
      in->IntPos[k] = Tp->P[i].IntPos[k];

    in->Hsml = Tp->fof_nearest_hsml[i];
  }

  /* routine to store or combine result data */
  void out2particle(fofdata_out *out, int i, int mode)
  {
    if(out->Distance < Tp->fof_nearest_distance[i])
      {
        Tp->fof_nearest_distance[i] = out->Distance;
        Tp->MinID[i]                = out->MinID;
        Tp->MinIDTask[i]            = out->MinIDTask;
#if defined(SUBFIND)
        Tp->PS[i].v.DM_Hsml = out->DM_Hsml;
#endif
      }
  }

  int evaluate(int target, int mode, int thread_id, int action, fofdata_in *in, int numnodes, node_info *firstnode, fofdata_out &out)
  {
    MyIntPosType *intpos = in->IntPos;
    double h             = in->Hsml;
    int index            = -1;
    double r2max         = MAX_REAL_NUMBER;

    /* Now start the actual tree-walk computation for this particle */

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

        int shmrank = Tree->TreeSharedMem_ThisTask;

        while(no >= 0)
          {
            if(no < Tree->MaxPart) /* single particle */
              {
                int p  = no;
                auto P = Tree->get_Pp(no, shmrank);

                no = Tree->get_nextnodep(shmrank)[no]; /* note: here shmrank cannot change */

                if(shmrank != Tree->TreeSharedMem_ThisTask)
                  Terminate("this routine may not consider shared memory particles");

                double dxyz[3];
                Tp->nearest_image_intpos_to_pos(P->IntPos, intpos, dxyz); /* converts the integer distance to floating point */

                double dist = h;

                if(fabs(dxyz[0]) > dist)
                  continue;
                if(fabs(dxyz[1]) > dist)
                  continue;
                if(fabs(dxyz[2]) > dist)
                  continue;

                double r2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];

                if(r2 < r2max && r2 < h * h)
                  {
                    index = p;
                    r2max = r2;
                  }
              }
            else if(no < Tree->MaxPart + Tree->MaxNodes) /* internal node */
              {
                if(mode == MODE_IMPORTED_PARTICLES)
                  {
                    if(no < Tree->FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the
                                                           branch */
                      break;
                  }

                fofnode *current = Tree->get_nodep(no, shmrank);

                if(current->level <= LEVEL_ALWAYS_OPEN)
                  {
                    /* we always open the root node (its full node length can't be stored in the integer type */
                    no      = current->nextnode; /* no change in shmrank expected here */
                    shmrank = current->nextnode_shmrank;
                    continue;
                  }

                int nosaved = no;

                no      = current->sibling; /* in case the node can be discarded */
                shmrank = current->sibling_shmrank;

                double dxyz[3];
                Tp->nearest_image_intpos_to_pos(current->center.da, intpos,
                                                dxyz); /* converts the integer distance to floating point */

                double len = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - current->level)) * Tp->FacIntToCoord;

                double dist = h + 0.5 * len;

                if(fabs(dxyz[0]) > dist)
                  continue;
                if(fabs(dxyz[1]) > dist)
                  continue;
                if(fabs(dxyz[2]) > dist)
                  continue;

                /* now test against the minimal sphere enclosing everything */
                dist += FACT1 * len;
                if(dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2] > dist * dist)
                  continue;

                int p = current->nextnode;

                /* in case the next node after opening is not a top-level node, we have either reached a leaf node or are in a local
                 * branch we need to do nothing if we would end up on different shared memory thread */
                if(p < Tree->MaxPart || (p >= Tree->FirstNonTopLevelNode && p < Tree->MaxPart + Tree->MaxNodes))
                  {
                    if(current->nextnode_shmrank != Tree->TreeSharedMem_ThisTask)
                      {
                        int task = D->ThisTask + current->nextnode_shmrank - Tree->TreeSharedMem_ThisTask;

                        if(target >= 0) /* export */
                          Tree->tree_export_node_threads_by_task_and_node(task, nosaved, target, &Thread);

                        no      = current->sibling; /* in case the node can be discarded */
                        shmrank = current->sibling_shmrank;
                        continue;
                      }
                  }

                no      = current->nextnode; /* ok, we need to open the node */
                shmrank = current->nextnode_shmrank;
              }
            else if(no >= Tree->ImportedNodeOffset) /* point from imported nodelist */
              {
                Terminate("do not expect imported points here");
              }
            else /* pseudo particle */
              {
                if(mode == MODE_LOCAL_PARTICLES)
                  if(target >= 0)
                    Tree->tree_export_node_threads(no, target, &Thread);

                no = Tree->get_nextnodep(shmrank)[no - Tree->MaxNodes];
                /* note: here shmrank does not need to change */
              }
          }
      }

    if(index >= 0)
      {
        out.Distance  = sqrt(r2max);
        out.MinID     = Tp->MinID[Tp->Head[index]];
        out.MinIDTask = Tp->MinIDTask[Tp->Head[index]];
#if defined(SUBFIND)
        out.DM_Hsml = Tp->PS[index].v.DM_Hsml;
#endif
      }
    else
      {
        out.Distance = MAX_REAL_NUMBER;
        out.MinID.set(0);
        out.MinIDTask = 0;
#if defined(SUBFIND)
        out.DM_Hsml = 0;
#endif
      }

    return 0;
  }
};

template <typename partset>
double fof<partset>::fof_find_nearest_dmparticle(void)
{
#ifdef LEAN
  return 0;
#endif

  double tstart = Logs.second();

  mpi_printf("FOF: Start finding nearest dm-particle (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());

  Tp->fof_nearest_distance = (MyFloat *)Mem.mymalloc("fof_nearest_distance", sizeof(MyFloat) * Tp->NumPart);
  Tp->fof_nearest_hsml     = (MyFloat *)Mem.mymalloc("fof_nearest_hsml", sizeof(MyFloat) * Tp->NumPart);

  int *TargetList = (int *)Mem.mymalloc("TargetList", Tp->NumPart * sizeof(int));

  int Nsearch = 0;

  for(int n = 0; n < Tp->NumPart; n++)
    {
      if(is_type_secondary_link_type(Tp->P[n].getType()))
        {
          Tp->fof_nearest_distance[n] = MAX_REAL_NUMBER;
          if(Tp->P[n].getType() == 0 && Tp->SphP[n].Hsml > 0)
            Tp->fof_nearest_hsml[n] = Tp->SphP[n].Hsml;
          else
            Tp->fof_nearest_hsml[n] = 0.1 * Tp->LinkL;

          TargetList[Nsearch++] = n;
        }
    }

  fofdata_comm<foftree<partset>, domain<partset>, partset> commpattern{FoFDomain, &FoFNgbTree, Tp};

  int iter = 0;
  long long ntot;
  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      double t0 = Logs.second();

      commpattern.execute(Nsearch, TargetList, MODE_DEFAULT);

      /* do final operations on results */
      int npleft = 0;
      for(int n = 0; n < Nsearch; n++)
        {
          int i = TargetList[n];

          if(Tp->fof_nearest_distance[i] > 1.0e29)
            {
              if(Tp->fof_nearest_hsml[i] < 4 * Tp->LinkL) /* we only search out to a maximum distance */
                {
                  /* need to redo this particle */
                  TargetList[npleft++] = i;
                  Tp->fof_nearest_hsml[i] *= 2.0;
                  if(iter >= MAXITER - 10)
                    {
                      double pos[3];
                      Tp->intpos_to_pos(Tp->P[i].IntPos, pos);

                      printf("FOF: i=%d task=%d ID=%d P[i].Type=%d Hsml=%g LinkL=%g nearest=%g pos=(%g|%g|%g)\n", i, ThisTask,
                             (int)Tp->P[i].ID.get(), Tp->P[i].getType(), Tp->fof_nearest_hsml[i], Tp->LinkL,
                             Tp->fof_nearest_distance[i], pos[0], pos[1], pos[2]);
                      myflush(stdout);
                    }
                }
              else
                {
                  Tp->fof_nearest_distance[i] = 0; /* we do not continue to search for this particle */
                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot, Communicator);

      Nsearch = npleft;

      double t1 = Logs.second();
      if(ntot > 0)
        {
          iter++;
          if(iter > 0)
            mpi_printf("FOF: fof-nearest iteration %d: need to repeat for %lld particles. (took = %g sec)\n", iter, ntot,
                       Logs.timediff(t0, t1));

          if(iter > MAXITER)
            Terminate("FOF: failed to converge in fof-nearest\n");
        }
    }
  while(ntot > 0);

  Mem.myfree(TargetList);
  Mem.myfree(Tp->fof_nearest_hsml);
  Mem.myfree(Tp->fof_nearest_distance);

  mpi_printf("FOF: done finding nearest dm-particle\n");

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

#endif /* of FOF */
