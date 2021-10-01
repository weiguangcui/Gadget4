/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_density.cc
 *
 *  \brief local matter density calculation for Subfind algorithm
 */

#include "gadgetconfig.h"

#ifdef SUBFIND
#ifndef SUBFIND_HBT

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

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

// Structure for communication during the density computation. Holds data that is sent to other processors.
struct subdens_in : data_in_generic
{
  MyIntPosType IntPos[3];
  MyFloat Hsml;
};

/* local data structure that holds results acquired on remote processors */
struct subdens_out
{
  int Ngb;
  MyFloat Rho;
#ifdef SUBFIND_STORE_LOCAL_DENSITY
  MyFloat VelDisp, Vx, Vy, Vz;
#endif
};

static int *DM_NumNgb;
#ifdef SUBFIND_STORE_LOCAL_DENSITY
static MyFloat *Vx, *Vy, *Vz;
#endif

template <typename T_tree, typename T_domain, typename T_partset>
class subdens_comm : public generic_comm<subdens_in, subdens_out, T_tree, T_domain, T_partset>
{
 public:
  typedef generic_comm<subdens_in, subdens_out, T_tree, T_domain, T_partset> gcomm;
  using gcomm::D;
  using gcomm::Thread;
  using gcomm::Tp;  // This makes sure that we can access Tp from the base class without having to use "this->Tp"
  using gcomm::Tree;

  /* need to call the base class constructor explicitly */
  subdens_comm(T_domain *dptr, T_tree *tptr, T_partset *pptr) : gcomm(dptr, tptr, pptr) {}

  /* routine that fills the relevant particle/cell data into the input structure defined above */
  void particle2in(subdens_in *in, int i) override
  {
    in->IntPos[0] = Tp->P[i].IntPos[0];
    in->IntPos[1] = Tp->P[i].IntPos[1];
    in->IntPos[2] = Tp->P[i].IntPos[2];
    in->Hsml      = Tp->PS[i].v.DM_Hsml;
  }

  /* routine to store or combine result data */
  void out2particle(subdens_out *out, int i, int mode) override
  {
    if(mode == MODE_LOCAL_PARTICLES) /* initial store */
      {
        DM_NumNgb[i]               = out->Ngb;
        Tp->PS[i].u.s.u.DM_Density = out->Rho;
#ifdef SUBFIND_STORE_LOCAL_DENSITY
        Vx[i]                    = out->Vx;
        Vy[i]                    = out->Vy;
        Vz[i]                    = out->Vz;
        Tp->PS[i].SubfindVelDisp = out->VelDisp;
#endif
      }
    else /* combine */
      {
        DM_NumNgb[i] += out->Ngb;
        Tp->PS[i].u.s.u.DM_Density += out->Rho;
#ifdef SUBFIND_STORE_LOCAL_DENSITY
        Vx[i] += out->Vx;
        Vy[i] += out->Vy;
        Vz[i] += out->Vz;
        Tp->PS[i].SubfindVelDisp += out->VelDisp;
#endif
      }
  }

  /*! This function represents the core of the SPH density computation. The
   *  target particle may either be local, or reside in the communication
   *  buffer.
   */
  int evaluate(int target, int mode, int thread_id, int action, subdens_in *in, int numnodes, node_info *firstnode,
               subdens_out &out) override
  {
    MyIntPosType *intpos = in->IntPos;
    double hsml          = in->Hsml;

    double h2    = hsml * hsml;
    double hinv  = 1.0 / hsml;
    double hinv3 = hinv * hinv * hinv;

    int numngb    = 0;
    double rhosum = 0;
#ifdef SUBFIND_STORE_LOCAL_DENSITY
    double vxsum = 0;
    double vysum = 0;
    double vzsum = 0;
    double v2sum = 0;
#endif

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
            int p, type;
            double mass, r2;
            typename T_partset::pdata *P;

            if(no < Tree->MaxPart) /* single particle */
              {
                p = no;
                P = Tree->get_Pp(no, shmrank);

                no = Tree->get_nextnodep(shmrank)[no]; /* note: here shmrank cannot change */

                double dxyz[3];
                Tp->nearest_image_intpos_to_pos(P->IntPos, intpos, dxyz); /* converts the integer distance to floating point */

                r2 = dxyz[0] * dxyz[0];
                if(r2 > h2)
                  continue;

                r2 += dxyz[1] * dxyz[1];
                if(r2 > h2)
                  continue;

                r2 += dxyz[2] * dxyz[2];
                if(r2 > h2)
                  continue;

                mass = P->getMass();
                type = P->getType();
              }
            else if(no < Tree->MaxPart + Tree->MaxNodes) /* internal node */
              {
                if(mode == MODE_IMPORTED_PARTICLES)
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

                continue;
              }
            else if(no >= Tree->ImportedNodeOffset) /* point from imported nodelist */
              {
                int n = no - Tree->ImportedNodeOffset;
                no    = Tree->Nextnode[no - Tree->MaxNodes];
                /* note: here shmrank cannot change */

                double dxyz[3];
                Tp->nearest_image_intpos_to_pos(Tree->Points[n].IntPos, intpos,
                                                dxyz); /* converts the integer distance to floating point */

                r2 = dxyz[0] * dxyz[0];
                if(r2 > h2)
                  continue;

                r2 += dxyz[1] * dxyz[1];
                if(r2 > h2)
                  continue;

                r2 += dxyz[2] * dxyz[2];
                if(r2 > h2)
                  continue;

                mass = Tree->Points[n].Mass;
                type = Tree->Points[n].Type;

                p = -1;
              }
            else /* pseudo particle */
              {
                if(mode == MODE_LOCAL_PARTICLES)
                  if(target >= 0) /* if no target is given, export will not occur */
                    Tree->tree_export_node_threads(no, target, &Thread);

                no = Tree->Nextnode[no - Tree->MaxNodes];
                continue;
              }

            if(r2 < h2)
              {
                if(is_type_primary_link_type(type))
                  {
                    numngb++;

                    if(p < 0)
                      Terminate("this should not occur");

#ifdef SUBFIND_STORE_LOCAL_DENSITY
                    vxsum += P->Vel[0];
                    vysum += P->Vel[1];
                    vzsum += P->Vel[2];
                    v2sum += P->Vel[0] * P->Vel[0] + P->Vel[1] * P->Vel[1] + P->Vel[2] * P->Vel[2];
#endif
                  }

                if(is_type_primary_link_type(type) || is_type_secondary_link_type(type))
                  {
                    double r = sqrt(r2);

                    double u = r * hinv, wk;

                    if(u < 0.5)
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    else
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                    rhosum += mass * wk;
                  }
              }
          }
      }

    out.Ngb = numngb;
    out.Rho = rhosum;
#ifdef SUBFIND_STORE_LOCAL_DENSITY
    out.Vx      = vxsum;
    out.Vy      = vysum;
    out.Vz      = vzsum;
    out.VelDisp = v2sum;
#endif
    return 0;
  }
};

template <typename partset>
double fof<partset>::subfind_density(void)
{
  mpi_printf("SUBFIND: finding total densities around all particles\n");

  double tstart = Logs.second();

  DM_NumNgb      = (int *)Mem.mymalloc_movable(&DM_NumNgb, "DM_NumNgb", sizeof(int) * Tp->NumPart);
  MyFloat *Left  = (MyFloat *)Mem.mymalloc_movable(&Left, "Left", sizeof(MyFloat) * Tp->NumPart);
  MyFloat *Right = (MyFloat *)Mem.mymalloc_movable(&Right, "Right", sizeof(MyFloat) * Tp->NumPart);

  int *targetlist = (int *)Mem.mymalloc_movable(&targetlist, "targetlist", sizeof(int) * Tp->NumPart);

#ifdef SUBFIND_STORE_LOCAL_DENSITY
  Vx = (MyFloat *)Mem.mymalloc("Vx", sizeof(MyFloat) * Tp->NumPart);
  Vy = (MyFloat *)Mem.mymalloc("Vy", sizeof(MyFloat) * Tp->NumPart);
  Vz = (MyFloat *)Mem.mymalloc("Vz", sizeof(MyFloat) * Tp->NumPart);
#endif

  int ntodo = 0;

  for(int i = 0; i < Tp->NumPart; i++)
    {
      Left[i] = Right[i] = 0;
      DM_NumNgb[i]       = 0;

      Tp->PS[i].u.s.u.DM_Density = 0;

#ifdef SUBFIND_STORE_LOCAL_DENSITY
      Tp->PS[i].SubfindHsml    = 0;
      Tp->PS[i].SubfindDensity = 0;
      Tp->PS[i].SubfindVelDisp = 0;
      if(is_type_primary_link_type(Tp->P[i].getType()) || is_type_secondary_link_type(Tp->P[i].getType()))
        targetlist[ntodo++] = i;
#else
      if(Tp->PS[i].GroupNr.get() != HALONR_MAX) /* do it only for particles is in a group */
        targetlist[ntodo++] = i;
#endif
    }

  subdens_comm<gravtree<partset>, domain<partset>, partset> commpattern{FoFDomain, &FoFGravTree, Tp};

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  long long ntot;
  int iter = 0;
  do
    {
      double t0 = Logs.second();

      commpattern.execute(ntodo, targetlist, MODE_DEFAULT);

      /* do final operations on results */
      int npleft = 0;
      for(int n = 0; n < ntodo; n++)
        {
          int i = targetlist[n];

          /* now check whether we had enough neighbours */
          if(abs(DM_NumNgb[i] - All.DesNumNgb) > All.MaxNumNgbDeviation &&
             ((Right[i] - Left[i]) > 1.0e-4 * Left[i] || Left[i] == 0 || Right[i] == 0))
            {
              /* need to redo this particle */
              targetlist[npleft++] = i;

              if(DM_NumNgb[i] < All.DesNumNgb)
                Left[i] = (MyFloat)std::max<double>(Tp->PS[i].v.DM_Hsml, Left[i]);
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

                  printf("SUBFIND: i=%d task=%d ID=%lld Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", i,
                         ThisTask, (long long)Tp->P[i].ID.get(), Tp->PS[i].v.DM_Hsml, Left[i], Right[i], (double)DM_NumNgb[i],
                         Right[i] - Left[i], pos[0], pos[1], pos[2]);
                  myflush(stdout);
                }

              if(Right[i] > 0 && Left[i] > 0)
                Tp->PS[i].v.DM_Hsml = (MyFloat)pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3.0);
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

      ntodo = npleft;

      sumup_large_ints(1, &npleft, &ntot, Communicator);

      double t1 = Logs.second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("SUBFIND: ngb iteration %2d: need to repeat for %15lld particles. (took %g sec)\n", iter, ntot,
                       Logs.timediff(t0, t1));

          if(iter > MAXITER)
            Terminate("failed to converge in neighbour iteration in subfind_density()\n");
        }
    }
  while(ntot > 0);

#ifdef SUBFIND_STORE_LOCAL_DENSITY

  double vel_to_phys = subfind_vel_to_phys_factor();

  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(DM_NumNgb[i] > 0)
        {
          Vx[i] /= DM_NumNgb[i];
          Vy[i] /= DM_NumNgb[i];
          Vz[i] /= DM_NumNgb[i];
          Tp->PS[i].SubfindVelDisp /= DM_NumNgb[i];
          Tp->PS[i].SubfindVelDisp = vel_to_phys * sqrt(Tp->PS[i].SubfindVelDisp - Vx[i] * Vx[i] - Vy[i] * Vy[i] - Vz[i] * Vz[i]);
        }
      else
        Tp->PS[i].SubfindVelDisp = 0;
    }

  Mem.myfree_movable(Vz);
  Mem.myfree_movable(Vy);
  Mem.myfree_movable(Vx);
#endif

  Mem.myfree_movable(targetlist);
  Mem.myfree_movable(Right);
  Mem.myfree_movable(Left);
  Mem.myfree_movable(DM_NumNgb);

#ifdef SUBFIND_STORE_LOCAL_DENSITY
  for(int i = 0; i < Tp->NumPart; i++)
    {
      Tp->PS[i].SubfindHsml    = Tp->PS[i].v.DM_Hsml;
      Tp->PS[i].SubfindDensity = Tp->PS[i].u.s.u.DM_Density;
    }
#endif

  double tend = Logs.second();
  return Logs.timediff(tstart, tend);
}

template <>
double fof<simparticles>::subfind_vel_to_phys_factor(void)
{
  if(All.ComovingIntegrationOn)
    return 1.0 / All.Time;
  else
    return 1.0;
}

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
template <>
double fof<lcparticles>::subfind_vel_to_phys_factor(void)
{
  return 1.0;
}
#endif

template <typename partset>
void fof<partset>::subfind_density_hsml_guess(void) /* set the initial guess for the smoothing length */
{
  double hsml_prev = 0;

  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(is_type_primary_link_type(Tp->P[i].getType()))
        {
          int no = FoFGravTree.Father[i];

          while(8.0 * All.DesNumNgb * Tp->P[i].getMass() > FoFGravTree.get_nodep(no)->mass)
            {
              int p = FoFGravTree.get_nodep(no)->father;

              if(p < 0)
                break;

              no = p;
            }

          double len = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - FoFGravTree.get_nodep(no)->level)) * Tp->FacIntToCoord;

          Tp->PS[i].v.DM_Hsml = hsml_prev =
              (pow(3.0 / (4.0 * M_PI) * All.DesNumNgb * Tp->P[i].getMass() / FoFGravTree.get_nodep(no)->mass, 1.0 / 3.0) * len);

          if(Tp->PS[i].v.DM_Hsml == 0)
            {
              double pos[3];
              Tp->intpos_to_pos(Tp->P[i].IntPos, pos);

              Terminate(
                  "zero hsml guess: Hsml=0 task=%d i=%d no=%d Nodes[no].len=%g Nodes[no].mass=%g P[i].Mass=%g type=%d ID=%llu  "
                  "pos=(%g|%g|%g)\n",
                  ThisTask, i, no, len, FoFGravTree.get_nodep(no)->mass, Tp->P[i].getMass(), Tp->P[i].getType(),
                  (long long)Tp->P[i].ID.get(), pos[0], pos[1], pos[2]);
            }
        }
      else
        {
          if(hsml_prev)
            Tp->PS[i].v.DM_Hsml = hsml_prev;
          else
            Tp->PS[i].v.DM_Hsml = All.SofteningTable[All.SofteningClassOfPartType[Tp->P[i].getType()]];
        }
    }
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
