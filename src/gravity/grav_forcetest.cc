/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file grav_forcetest.cc
 *
 *  \brief routines for testing the force accuracy though comparison with direct summation
 */

#include "gadgetconfig.h"

#ifdef FORCETEST

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
#include "../gravity/ewald.h"
#include "../gravity/grav_forcetest.h"
#include "../gravtree/gravtree.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mpi_utils/generic_comm.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../pm/pm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/*  This function computes the gravitational forces for all active particles.
 *  A new tree is constructed, if the number of force computations since
 *  it's last construction exceeds some fraction of the total
 *  particle number, otherwise tree nodes are dynamically updated if needed.
 */

/*  FORCETEST_TESTFORCELAW=1   for special test to check force law for TreePM
 *  FORCETEST_TESTFORCELAW=2   for special test to check force law for TreePM+PLACEHIGHRESREGION
 */

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
struct frctest_in : data_in_generic
{
  MyIntPosType IntPos[3];
  unsigned char Type;
#if NSOFTCLASSES > 1
  unsigned char SofteningClass;
#endif
};

/* local data structure that holds results acquired on remote processors */
struct frctest_out
{
  double Acc[3];
  double Pot;
  double DistToID1;
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
  double AccShortRange[3];
  double PotShortRange;
#ifdef PLACEHIGHRESREGION
  double AccVeryShortRange[3];
  double PotVeryShortRange;
#endif
#endif
};

typedef generic_comm<frctest_in, frctest_out, gravtree<simparticles>, domain<simparticles>, simparticles> my_comm;

class frctest_comm : public my_comm
{
 public:
  // need to implement a constructor that calls the constructor of the base class
  frctest_comm(domain<simparticles> *dptr, gravtree<simparticles> *tptr, simparticles *pptr, gravtest *gtptr)
      : my_comm(dptr, tptr, pptr)
  {
  }

  using my_comm::D;
  using my_comm::Thread;
  using my_comm::Tp;
  using my_comm::Tree;

  void particle2in(frctest_in *in, int i) override
  {
    for(int j = 0; j < 3; j++)
      in->IntPos[j] = Tp->P[i].IntPos[j];

    in->Type = Tp->P[i].getType();
#if NSOFTCLASSES > 1
    in->SofteningClass = Tp->P[i].getSofteningClass();
#endif
  }

  void out2particle(frctest_out *out, int i, int mode) override
  {
    if(mode == MODE_LOCAL_PARTICLES) /* initial store */
      {
        Tp->P[i].GravAccelDirect[0] = out->Acc[0];
        Tp->P[i].GravAccelDirect[1] = out->Acc[1];
        Tp->P[i].GravAccelDirect[2] = out->Acc[2];
        Tp->P[i].PotentialDirect    = out->Pot;

        Tp->P[i].DistToID1 = out->DistToID1;
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
        Tp->P[i].GravAccelShortRange[0] = out->AccShortRange[0];
        Tp->P[i].GravAccelShortRange[1] = out->AccShortRange[1];
        Tp->P[i].GravAccelShortRange[2] = out->AccShortRange[2];
        Tp->P[i].PotentialShortRange    = out->PotShortRange;
#ifdef PLACEHIGHRESREGION
        Tp->P[i].GravAccelVeryShortRange[0] = out->AccVeryShortRange[0];
        Tp->P[i].GravAccelVeryShortRange[1] = out->AccVeryShortRange[1];
        Tp->P[i].GravAccelVeryShortRange[2] = out->AccVeryShortRange[2];
        Tp->P[i].PotentialVeryShortRange    = out->PotVeryShortRange;
#endif
#endif
      }
    else /* combine */
      {
        Tp->P[i].GravAccelDirect[0] += out->Acc[0];
        Tp->P[i].GravAccelDirect[1] += out->Acc[1];
        Tp->P[i].GravAccelDirect[2] += out->Acc[2];
        Tp->P[i].PotentialDirect += out->Pot;
        if(out->DistToID1 > 0)
          Tp->P[i].DistToID1 = out->DistToID1;
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
        Tp->P[i].GravAccelShortRange[0] += out->AccShortRange[0];
        Tp->P[i].GravAccelShortRange[1] += out->AccShortRange[1];
        Tp->P[i].GravAccelShortRange[2] += out->AccShortRange[2];
        Tp->P[i].PotentialShortRange += out->PotShortRange;
#ifdef PLACEHIGHRESREGION
        Tp->P[i].GravAccelVeryShortRange[0] += out->AccVeryShortRange[0];
        Tp->P[i].GravAccelVeryShortRange[1] += out->AccVeryShortRange[1];
        Tp->P[i].GravAccelVeryShortRange[2] += out->AccVeryShortRange[2];
        Tp->P[i].PotentialVeryShortRange += out->PotVeryShortRange;
#endif
#endif
      }
  }

  int evaluate(int target, int mode, int thread_id, int action, frctest_in *in, int numnodes, node_info *firstnode,
               frctest_out &out) override
  {
    /* make sure that the particle is exported to all other tasks exactly once */
    if(mode == MODE_LOCAL_PARTICLES)
      {
        for(int n = 0; n < this->D->NTopleaves; n++)
          {
            int task = D->TaskOfLeaf[n];

            if(task == D->ThisTask)
              continue;

            if(Thread.Exportflag[task] == target)
              continue;

            int no = n + this->Tree->MaxPart + this->Tree->MaxNodes; /* a pseudo node for this task */

            this->Tree->tree_export_node_threads(no, target, &Thread);
          }
      }

    MyIntPosType *intpos = in->IntPos;

    out.Pot = 0;

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
    out.PotShortRange = 0;
#ifdef PLACEHIGHRESREGION
    out.PotVeryShortRange = 0;
#endif
#endif
    for(int i = 0; i < 3; i++)
      {
        out.Acc[i] = 0;
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
        out.AccShortRange[i] = 0;
#ifdef PLACEHIGHRESREGION
        out.AccVeryShortRange[i] = 0;
#endif
#endif
      }

    double disttoid1 = 0;

    for(int idx = 0; idx < Tp->nsource; idx++)
      {
        int j = Tp->indexlist[idx];

        double hmax;
#if NSOFTCLASSES > 1
        double h_i = All.ForceSoftening[in->SofteningClass];
        double h_j = All.ForceSoftening[Tp->P[j].getSofteningClass()];

        if(h_j > h_i)
          hmax = h_j;
        else
          hmax = h_i;
#else
        hmax = All.ForceSoftening[0];
#endif

        double dxyz[3];
        Tp->nearest_image_intpos_to_pos(Tp->P[j].IntPos, intpos, dxyz); /* converts the integer distance to floating point */

        double r2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];

        double mass = Tp->P[j].getMass();

        /* now evaluate the multipole moment */

        double r = sqrt(r2);

        if(Tp->P[j].ID.get() == 1)
          disttoid1 = r;

        /* we compute 3 different forces:
         * (1) The correct direct summation force, if needed with Ewald correction: ftrue
         * In the case of PM:
         * (2) The short range direct summation force with only the erfc cut-off (this is what the tree can at best deliver): fsr
         * (3) The expected PM force based on the long-range part of the Ewald sum. This is equal to ftrue - fsr - fsfr_periodic_images
         * */

        double wp_newton, fac_newton;

        if(r > 0)
          {
            fac_newton = mass / (r2 * r);
            wp_newton  = -mass / r;
          }
        else
          {
            fac_newton = 0;
            wp_newton  = 0;
          }

        double fac, wp;

        if(r >= hmax)
          {
            fac = fac_newton;
            wp  = wp_newton;
          }
        else
          {
            double h_inv  = 1.0 / hmax;
            double h3_inv = h_inv * h_inv * h_inv;
            double u      = r * h_inv;

            if(u < 0.5)
              {
                double u2 = u * u;
                fac       = mass * h3_inv * (SOFTFAC1 + u2 * (SOFTFAC2 * u + SOFTFAC3));
                wp        = mass * h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
              }
            else
              {
                double u2 = u * u, u3 = u2 * u;
                fac = mass * h3_inv * (SOFTFAC8 + SOFTFAC9 * u + SOFTFAC10 * u2 + SOFTFAC11 * u3 + SOFTFAC12 / u3);
                wp  = mass * h_inv * (SOFTFAC13 + SOFTFAC14 / u + u2 * (SOFTFAC1 + u * (SOFTFAC15 + u * (SOFTFAC16 + SOFTFAC17 * u))));
              }
          }

          // The Newtonian force is:      fac * dxyz
          // The Newtonian potential is:  wp

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
        {
          double asmth = Tp->Asmth[0];
          double u     = 0.5 / asmth * r;

          double factor_force = (erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u) - 1.0);
          double factor_pot   = erfc(u);

          double facs = fac + fac_newton * factor_force;
          double wps  = wp + (r > 0 ? wp_newton * (factor_pot - 1.0) : mass / (asmth * sqrt(M_PI)));

          double acc_short_x = dxyz[0] * facs;
          double acc_short_y = dxyz[1] * facs;
          double acc_short_z = dxyz[2] * facs;

#ifndef GRAVITY_TALLBOX
          double alpha = 0.5 / asmth;
          double pot_short =
              wps + mass * M_PI / (alpha * alpha * All.BoxSize * All.BoxSize * All.BoxSize) * (LONG_X * LONG_Y * LONG_Z);
#else
          double pot_short = wps; /* the constant potential term is here computed as part of the long-range force and not the
                                     short-range force, unlike in the ordinary periodic case */
#endif
          out.AccShortRange[0] += acc_short_x;
          out.AccShortRange[1] += acc_short_y;
          out.AccShortRange[2] += acc_short_z;
          out.PotShortRange += pot_short;

#ifdef PLACEHIGHRESREGION
          if(Tp->check_high_res_point_location(Tp->P[j].IntPos) == FLAG_INSIDE &&
             Tp->check_high_res_point_location(intpos) == FLAG_INSIDE)
            {
              double asmth = Tp->Asmth[1];
              double u     = 0.5 / asmth * r;

              double factor_force = (erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u) - 1.0);
              double factor_pot   = erfc(u);

              double facs = fac + fac_newton * factor_force;
              double wps  = wp + (r > 0 ? wp_newton * (factor_pot - 1.0) : mass / (asmth * sqrt(M_PI)));

              double alpha = 0.5 / asmth;

              acc_short_x = dxyz[0] * facs;
              acc_short_y = dxyz[1] * facs;
              acc_short_z = dxyz[2] * facs;

              pot_short = wps + mass * M_PI / (alpha * alpha * All.BoxSize * All.BoxSize * All.BoxSize) * (LONG_X * LONG_Y * LONG_Z);
            }

          out.AccVeryShortRange[0] += acc_short_x;
          out.AccVeryShortRange[1] += acc_short_y;
          out.AccVeryShortRange[2] += acc_short_z;
          out.PotVeryShortRange += pot_short;
#endif
        }
#endif

        // the direct force is the ordinary Newtonian force
        out.Acc[0] += fac * dxyz[0];
        out.Acc[1] += fac * dxyz[1];
        out.Acc[2] += fac * dxyz[2];
        out.Pot += wp;

#ifdef PERIODIC
        // and in the periodic case, we add in the correction potential/force
        ewald_data ew;
        Ewald.ewald_gridlookup(Tp->P[j].IntPos, intpos, ewald::POINTMASS, ew);

        out.Pot += mass * ew.D0phi;
        out.Acc[0] += mass * ew.D1phi[0];
        out.Acc[1] += mass * ew.D1phi[1];
        out.Acc[2] += mass * ew.D1phi[2];
#endif
      }

    out.DistToID1 = disttoid1;

    return 0;
  }
};

void gravtest::gravity_forcetest(int timebin)
{
  TIMER_START(CPU_FORCETEST);

  int *TargetList = (int *)Mem.mymalloc("TargetList", Sp->NumPart * sizeof(int));
  int nloc        = 0;

  particle_data *P = Sp->P;

  // create a random selection of target particles for which we compute direct summation forces
  for(int idx = 0; idx < Sp->TimeBinsGravity.NActiveParticles; idx++)
    {
      int target = Sp->TimeBinsGravity.ActiveParticleList[idx];

      if(target < 0)
        continue;

#ifdef FORCETEST_FIXEDPARTICLESET
      if(All.NumCurrentTiStep == 0)
        {
          if(get_random_number() < FORCETEST)
            P[target].SelectedFlag = true;
          else
            P[target].SelectedFlag = false;
        }
      if(P[target].SelectedFlag)
        TargetList[nloc++] = target;
#else
      if(get_random_number() < FORCETEST)
        TargetList[nloc++] = target;
#endif
    }

  long long ntot;
  sumup_large_ints(1, &nloc, &ntot, D->Communicator);

  /* we pull put a separate list of the particles with non-zero masses to accelerate some of our special tests
   * where there are potentially many particles with zero masses
   */

  Sp->nsource   = 0;
  Sp->indexlist = (int *)Mem.mymalloc("indexlist", Sp->NumPart * sizeof(int));

#ifdef HIERARCHICAL_GRAVITY
  int numidx = Sp->TimeBinsGravity.NActiveParticles;
#else
  int numidx = Sp->NumPart;
#endif

  for(int idx = 0; idx < numidx; idx++)
    {
#ifdef HIERARCHICAL_GRAVITY
      int target = Sp->TimeBinsGravity.ActiveParticleList[idx];

      if(target < 0)
        continue;
#else
      int target = idx;
#endif

      if(P[target].getMass() == 0)
        continue;

      Sp->indexlist[Sp->nsource++] = target;
    }

  D->mpi_printf("FORCETEST: Testing forces for %lld particles out of %lld active ones.\n", ntot,
                Sp->TimeBinsGravity.GlobalNActiveParticles);

  frctest_comm commpattern(this->D, this->GravTree, this->Sp, this);

  double t0 = Logs.second();

  commpattern.execute(nloc, TargetList, MODE_DEFAULT);

  double t1   = Logs.second();
  double maxt = Logs.timediff(t0, t1);

  D->mpi_printf("FORCETEST: Testing forces took %g sec.\n", maxt);

  Mem.myfree(Sp->indexlist);
  Sp->indexlist = NULL;

  /* now add things for comoving integration */

  if(All.ComovingIntegrationOn)
    {
#ifndef PERIODIC
      double fac1 = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(int idx = 0; idx < nloc; idx++)
        {
          int i = TargetList[idx];

          double pos[3];
          Sp->intpos_to_pos(Sp->P[i].IntPos, pos);

          for(int j = 0; j < 3; j++)
            Sp->P[i].GravAccelDirect[j] += fac1 * pos[j];
        }
#endif
    }

  /*  muliply by G */
  for(int idx = 0; idx < nloc; idx++)
    {
      int i = TargetList[idx];

      for(int j = 0; j < 3; j++)
        {
          Sp->P[i].GravAccelDirect[j] *= All.G;
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
          Sp->P[i].GravAccelShortRange[j] *= All.G;
#ifdef PLACEHIGHRESREGION
          Sp->P[i].GravAccelVeryShortRange[j] *= All.G;
#endif
#endif
#if(FORCETEST_TESTFORCELAW == 3)
          Sp->P[i].GravAccelDirectTest[j] *= All.G;
#endif
        }

#if NSOFTCLASSES > 1
      double selfpot = Sp->P[i].getMass() / (All.ForceSoftening[Sp->P[i].getSofteningClass()] / 2.8);
#else
      double selfpot = Sp->P[i].getMass() / (All.ForceSoftening[0] / 2.8);
#endif

      Sp->P[i].PotentialDirect += selfpot; /* remove self-potential */
      Sp->P[i].PotentialDirect *= All.G;

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
      Sp->P[i].PotentialShortRange += selfpot; /* remove self-potential */
      Sp->P[i].PotentialShortRange *= All.G;

#ifdef PLACEHIGHRESREGION
      Sp->P[i].PotentialVeryShortRange += selfpot; /* remove self-potential */
      Sp->P[i].PotentialVeryShortRange *= All.G;
#endif
#endif
#if(FORCETEST_TESTFORCELAW == 3)
      Sp->P[i].PotentialDirectTest += selfpot; /* remove self-potential */
      Sp->P[i].PotentialDirectTest *= All.G;
#endif
    }

  /* Finally, the following factor allows a computation of cosmological simulation
     with vacuum energy in physical coordinates */

  if(All.ComovingIntegrationOn == 0)
    {
      double fac1 = All.OmegaLambda * All.Hubble * All.Hubble;

      for(int idx = 0; idx < nloc; idx++)
        {
          int i = TargetList[idx];

          double pos[3];
          Sp->intpos_to_pos(Sp->P[i].IntPos, pos);

          for(int j = 0; j < 3; j++)
            Sp->P[i].GravAccelDirect[j] += fac1 * pos[j];
        }
    }

  /* now output the forces to a file */

  int *nloc_tab = (int *)Mem.mymalloc("nloc_tab", D->NTask * sizeof(int));
  MPI_Allgather(&nloc, 1, MPI_INT, nloc_tab, 1, MPI_INT, D->Communicator);

  for(int nthis = 0; nthis < D->NTask; nthis++)
    {
      if(nloc_tab[nthis] > 0)
        {
          if(nthis == D->ThisTask)
            {
              char buf[MAXLEN_PATH_EXTRA];
              snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "forcetest.txt");

              if(!(Logs.FdForceTest = fopen(buf, "a")))
                Terminate("error in opening file '%s'\n", buf);

              for(int idx = 0; idx < nloc; idx++)
                {
                  int i = TargetList[idx];

                  double pos[3];
                  Sp->intpos_to_pos(Sp->P[i].IntPos, pos);

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)

#ifdef PLACEHIGHRESREGION
                  int flaginside = Sp->check_high_res_point_location(P[i].IntPos);
                  fprintf(
                      Logs.FdForceTest,
                      "%d %d %lld  %g  %17.12g %17.12g %17.12g  %17.12g  %12.12g %17.12g %17.12g  %17.12g %17.12g %17.12g  %17.12g "
                      "%17.12g %17.12g  %17.12g %17.12g %17.12g  %17.12g %17.12g %17.12g  %17.12g "
                      "%17.12g %17.12g     %17.12g %17.12g %17.12g %17.12g %17.12g %17.12g   %17.12g\n",
                      P[i].getType(), flaginside, (long long)P[i].ID.get(), All.Time, pos[0], pos[1], pos[2], P[i].DistToID1,
                      P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2], P[i].GravAccelShortRange[0],
                      P[i].GravAccelShortRange[1], P[i].GravAccelShortRange[2], P[i].GravAccelVeryShortRange[0],
                      P[i].GravAccelVeryShortRange[1], P[i].GravAccelVeryShortRange[2], P[i].GravAccel[0], P[i].GravAccel[1],
                      P[i].GravAccel[2], P[i].GravPM[0], P[i].GravPM[1], P[i].GravPM[2], P[i].GravAccelHPM[0], P[i].GravAccelHPM[1],
                      P[i].GravAccelHPM[2], P[i].PotentialDirect, P[i].PotentialShortRange, P[i].PotentialVeryShortRange,
                      P[i].Potential, P[i].PM_Potential, P[i].PotentialHPM, Sp->Asmth[1]);
#else
                  fprintf(
                      Logs.FdForceTest,
                      "%d %d %lld  %g  %17.12g %17.12g %17.12g  %17.12g  %17.12g %17.12g %17.12g  %17.12g %17.12g %17.12g  %17.12g "
                      "%17.12g "
                      "%17.12g  %17.12g %17.12g %17.12g  %17.12g %17.12g %17.12g  %17.12g\n",
                      P[i].getType(), timebin, (long long)P[i].ID.get(), All.Time, pos[0], pos[1], pos[2], P[i].DistToID1,
                      P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2], P[i].GravAccelShortRange[0],
                      P[i].GravAccelShortRange[1], P[i].GravAccelShortRange[2], P[i].GravAccel[0], P[i].GravAccel[1],
                      P[i].GravAccel[2], P[i].GravPM[0], P[i].GravPM[1], P[i].GravPM[2], P[i].PotentialDirect,
                      P[i].PotentialShortRange, P[i].Potential, P[i].PM_Potential);
#endif
#else
                  fprintf(Logs.FdForceTest,
                          "%d %d %lld %g  %17.12g %17.12g %17.12g  %17.12g  %17.12g %17.12g %17.12g  %17.12g %17.12g %17.12g  %17.12g "
                          "%17.12g\n",
                          P[i].getType(), timebin, (long long)P[i].ID.get(), All.Time, pos[0], pos[1], pos[2], P[i].DistToID1,
                          P[i].GravAccelDirect[0], P[i].GravAccelDirect[1], P[i].GravAccelDirect[2], Sp->P[i].GravAccel[0],
                          Sp->P[i].GravAccel[1], Sp->P[i].GravAccel[2], P[i].Potential, P[i].PotentialDirect);
#endif
                }

              fclose(Logs.FdForceTest);
            }

          MPI_Barrier(D->Communicator);
        }
    }
  Mem.myfree(nloc_tab);

  /* Now the force computation is finished */
  if(D->ThisTask == 0)
    {
      double costtotal = Sp->NumPart * ntot;

      fprintf(Logs.FdTimings, "DIRECT Nf= %lld  step=%d  timebin=%d  part/sec=%g | %g  ia/part=%g\n\n", ntot, All.NumCurrentTiStep,
              timebin, ((double)ntot) / (Sp->NTask * maxt + 1.0e-20), ntot / ((maxt + 1.0e-20) * Sp->NTask),
              ((double)(costtotal)) / (ntot + 1.0e-20));

      myflush(Logs.FdTimings);
    }

  Mem.myfree(TargetList);

  TIMER_STOP(CPU_FORCETEST);
}

#ifdef FORCETEST_TESTFORCELAW /* in this option we assume NSOFTCLASSES >= 2 for FORCETEST_TESTFORCELAW=1, and  NSOFTCLASSES >= 3 for \
                                 FORCETEST_TESTFORCELAW=2*/
void sim::gravity_forcetest_testforcelaw(void)
{
  int Ncycles = 40;
  double xyz[6], eps;

  NgbTree.treefree();
  Sp.mark_active_timebins();

  for(int cycle = 0; cycle < Ncycles; cycle++)
    {
      Domain.mpi_printf("\nTEST-FORCE-LAW: cycle=%d|%d ----------------------------------\n\n", cycle, Ncycles);

      double epsloc = 0, xyzloc[6] = {0, 0, 0, 0, 0, 0};

      /* set particle with ID=1 to new random coordinate in box */
      for(int n = 0; n < Sp.NumPart; n++)
        {
          Sp.P[n].setType(1);
          Sp.P[n].setSofteningClass(1);

          if(Sp.P[n].ID.get() == 1)
            {
              xyzloc[0] = All.BoxSize / LONG_X * get_random_number();
              xyzloc[1] = All.BoxSize / LONG_Y * get_random_number();
              xyzloc[2] = All.BoxSize / LONG_Z * get_random_number();

              for(int i = 0; i < 3; i++)
                xyzloc[3 + i] = xyzloc[i];

#if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 2)
              if(get_random_number() < 0.5)
                {
                  xyzloc[3] = All.BoxSize / LONG_X * get_random_number();
                  xyzloc[4] = All.BoxSize / LONG_Y * get_random_number();
                  xyzloc[5] = All.BoxSize / LONG_Z * get_random_number();
                }
#endif
              Sp.pos_to_intpos(&xyzloc[3], Sp.P[n].IntPos);

              epsloc = All.ForceSoftening[Sp.P[n].getSofteningClass()];
            }
        }

      MPI_Allreduce(xyzloc, xyz, 6, MPI_DOUBLE, MPI_SUM, Communicator);
      MPI_Allreduce(&epsloc, &eps, 1, MPI_DOUBLE, MPI_SUM, Communicator);

      double rmin = 0.01 * eps;
      double rmax = sqrt(pow(0.5 * All.BoxSize / LONG_X, 2) + pow(0.5 * All.BoxSize / LONG_Y, 2) + pow(0.5 * All.BoxSize / LONG_Z, 2));

      for(int n = 0; n < Sp.NumPart; n++)
        {
          if(Sp.P[n].ID.get() != 1)
            {
              double r     = exp(log(rmin) + (log(rmax) - log(rmin)) * get_random_number());
              double theta = acos(2 * get_random_number() - 1);
              double phi   = 2 * M_PI * get_random_number();

#if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 2)
              if(get_random_number() < 0.5)
                {
                  r = exp(log(rmin) + (log(rmax / 100) - log(rmin)) * get_random_number());
                  Sp.P[n].setType(2);
                }
              else
                {
                  Sp.P[n].setType(3);
                }
#endif

              double dx = r * sin(theta) * cos(phi);
              double dy = r * sin(theta) * sin(phi);
              double dz = r * cos(theta);

              double pos[3];
              pos[0] = xyz[0] + dx;
              pos[1] = xyz[1] + dy;
              pos[2] = xyz[2] + dz;

#if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 2)
              if(Sp.P[n].getType() == 3)
                {
                  pos[0] = xyz[3] + dx;
                  pos[1] = xyz[4] + dy;
                  pos[2] = xyz[5] + dz;
                }
#endif
              Sp.pos_to_intpos(pos, Sp.P[n].IntPos);
              Sp.constrain_intpos(Sp.P[n].IntPos);
            }
        }

      Domain.domain_free();
      Domain.domain_decomposition(STANDARD); /* do domain decomposition if needed */

      /* allocate space for gravity accelerations */

      Sp.TimeBinsGravity.NActiveParticles = 0;
      for(int timebin = All.HighestSynchronizedTimeBin; timebin >= 0; timebin--)
        {
          for(int i = Sp.TimeBinsGravity.FirstInTimeBin[timebin]; i >= 0; i = Sp.TimeBinsGravity.NextInTimeBin[i])
            Sp.TimeBinsGravity.ActiveParticleList[Sp.TimeBinsGravity.NActiveParticles++] = i;
        }

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
      gravity_long_range_force();
#endif

      compute_grav_accelerations(All.HighestActiveTimeBin);
    }

  endrun();
}
#endif

#endif
