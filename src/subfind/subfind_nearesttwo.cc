/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_nearesttwo.cc
 *
 *  \brief determine the nearest two denser neighbours for linking them in excursion set formalism of SUBFIND
 */

#include "gadgetconfig.h"

#ifdef SUBFIND
#ifndef SUBFIND_HBT

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
#include "../mpi_utils/mpi_utils.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
struct sngb_in : data_in_generic
{
  MyIntPosType IntPos[3];
  MyIDType ID;
  MyFloat Hsml;
  MyFloat Density;
  MyFloat Dist[2];
  int Count;
  location Index[2];
};

struct sngb_out
{
  MyFloat Dist[2];
  location Index[2];
  int Count;
};

template <typename T_tree, typename T_domain, typename T_partset>
class sngb_comm : public generic_comm<sngb_in, sngb_out, T_tree, T_domain, T_partset>
{
 public:
  typedef generic_comm<sngb_in, sngb_out, T_tree, T_domain, T_partset> gcomm;
  using gcomm::D;
  using gcomm::Thread;
  using gcomm::Tp;  // This makes sure that we can access Tp from the base class without having to use "this->Tp"
  using gcomm::Tree;

  /* need to call the base class constructor explicitly */
  sngb_comm(T_domain *dptr, T_tree *tptr, T_partset *pptr) : gcomm(dptr, tptr, pptr) {}

  /* routine that fills the relevant particle/cell data into the input structure defined above */
  void particle2in(sngb_in *in, int i) override
  {
    in->IntPos[0] = Tp->P[i].IntPos[0];
    in->IntPos[1] = Tp->P[i].IntPos[1];
    in->IntPos[2] = Tp->P[i].IntPos[2];

    in->Hsml    = Tp->PS[i].v.DM_Hsml;
    in->ID      = Tp->P[i].ID.get();
    in->Density = Tp->PS[i].u.s.u.DM_Density;
    in->Count   = Tp->PS[i].nearest.count;
    for(int k = 0; k < Tp->PS[i].nearest.count; k++)
      {
        in->Dist[k]  = Tp->R2Loc[i].dist[k];
        in->Index[k] = Tp->PS[i].nearest.index[k];
      }
  }

  /* routine to store or combine result data */
  void out2particle(sngb_out *out, int i, int mode) override
  {
    if(mode == MODE_LOCAL_PARTICLES) /* initial store */
      {
        Tp->PS[i].nearest.count = out->Count;

        for(int k = 0; k < out->Count; k++)
          {
            Tp->R2Loc[i].dist[k]       = out->Dist[k];
            Tp->PS[i].nearest.index[k] = out->Index[k];
          }
      }
    else /* combine */
      {
        for(int k = 0; k < out->Count; k++)
          {
            if(Tp->PS[i].nearest.count >= 1)
              if(Tp->PS[i].nearest.index[0] == out->Index[k])
                continue;

            if(Tp->PS[i].nearest.count == 2)
              if(Tp->PS[i].nearest.index[1] == out->Index[k])
                continue;

            int l;

            if(Tp->PS[i].nearest.count < 2)
              {
                l = Tp->PS[i].nearest.count;
                Tp->PS[i].nearest.count++;
              }
            else
              {
                l = (Tp->R2Loc[i].dist[0] > Tp->R2Loc[i].dist[1]) ? 0 : 1;

                if(out->Dist[k] >= Tp->R2Loc[i].dist[l])
                  continue;
              }

            Tp->R2Loc[i].dist[l]       = out->Dist[k];
            Tp->PS[i].nearest.index[l] = out->Index[k];

            if(Tp->PS[i].nearest.count == 2)
              if(Tp->PS[i].nearest.index[0] == Tp->PS[i].nearest.index[1])
                Terminate("this is not supposed to happen");
          }
      }
  }

  /*! This function represents the core of the neighbor search. The target particle may either be local, or reside in the communication
   *  buffer.
   */
  int evaluate(int target, int mode, int thread_id, int action, sngb_in *in, int numnodes, node_info *firstnode,
               sngb_out &out) override
  {
    MyIntPosType *intpos = in->IntPos;
    MyIDType ID          = in->ID;
    double density       = in->Density;
    double hsml          = in->Hsml;
    int count            = in->Count;

    location index[2];
    double dist[2];
    for(int k = 0; k < count; k++)
      {
        dist[k]  = in->Dist[k];
        index[k] = in->Index[k];
      }

    if(count == 2)
      if(index[0] == index[1])
        {
          Terminate("target=%d mode=%d\n", target, mode);
        }

    count = 0;

    hsml *= 1.00001; /* prevents that the most distant neighbour on the edge of the search region may not be found.
                      * (needed for consistency with serial algorithm)
                      */

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
                auto *P          = Tree->get_Pp(no, shmrank);
                subfind_data *PS = Tree->get_PSp(no, shmrank);

                no = Tree->get_nextnodep(shmrank)[no]; /* note: here shmrank cannot change */

                if(P->ID.get() != ID) /* exclude the target particle itself */
                  {
                    if(PS->u.s.u.DM_Density > density) /* we only need to look at neighbors that are denser */
                      {
                        /* converts the integer distance to floating point */
                        double dxyz[3];
                        Tp->nearest_image_intpos_to_pos(P->IntPos, intpos, dxyz);

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

                        // ok, we found a particle. Only store up to two closest ones.

                        int task = D->ThisTask + (shmrank - Tree->TreeSharedMem_ThisTask);

                        if(task < 0 || task >= D->NTask)
                          Terminate("illegal task=%d  D->NTask=%d", task, D->NTask);

                        if(count < 2)
                          {
                            dist[count]  = r2;
                            index[count] = {task, PS->InvIndex}; /* note: ThisTask refers here to the Subdomain */
                            count++;
                          }
                        else
                          {
                            int k = (dist[0] > dist[1]) ? 0 : 1;

                            if(r2 < dist[k])
                              {
                                dist[k]  = r2;
                                index[k] = {task, PS->InvIndex}; /* note: ThisTask refers here to the Subdomain */
                              }
                          }
                      }
                  }
              }
            else if(no < Tree->MaxPart + Tree->MaxNodes) /* internal node */
              {
                if(mode == 1)
                  {
                    if(no < Tree->FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the
                                                           branch */
                      {
                        break;
                      }
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
            else if(no >= Tree->ImportedNodeOffset) /* point from imported nodelist */
              {
                Terminate("do not expect imported points here");
              }
            else /* pseudo particle */
              {
                if(mode == MODE_LOCAL_PARTICLES)
                  if(target >= 0) /* note: if no target is given, export will not occur */
                    Tree->tree_export_node_threads(no, target, &Thread);

                no = Tree->Nextnode[no - Tree->MaxNodes];
              }
          }
      }

    out.Count = count;

    for(int k = 0; k < count; k++)
      {
        out.Dist[k]  = dist[k];
        out.Index[k] = index[k];
      }

    return 0;
  }
};

template <typename partset>
void fof<partset>::subfind_find_nearesttwo(domain<partset> *SubDomain, int num, int *list)
{
  subfind_collective_printf("SUBFIND: root-task=%d: Start finding nearest two.\n", ThisTask);

  for(int i = 0; i < num; i++)
    Tp->PS[list[i]].nearest.count = 0;

  /* create an object for handling the communication */
  sngb_comm<gravtree<partset>, domain<partset>, partset> commpattern{SubDomain, &FoFGravTree, Tp};

  commpattern.execute(num, list, MODE_DEFAULT);

  subfind_collective_printf("SUBFIND: root-task=%d: Done with nearest two.\n", ThisTask);
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
