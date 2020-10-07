/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file gravity.cc
 *
 * \brief main driver routines for computing the gravitational accelerations for all active particles
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fmm/fmm.h"
#include "../gravtree/gravtree.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/*! \brief This routine computes the gravitational accelerations for all active particles.
 *
 * If the particle mesh is used and the current time step
 * requires a PM force computation, new long range forces are
 * computed by long_range_force(). Then the shortrange tree forces
 * are computed by gravity(). The force tree is rebuild every time step.
 */
void sim::compute_grav_accelerations(int timebin)
{
  sumup_large_ints(1, &Sp.TimeBinsGravity.NActiveParticles, &Sp.TimeBinsGravity.GlobalNActiveParticles, Communicator);

  mpi_printf("ACCEL: Start tree gravity force computation... (%lld particles)\n", Sp.TimeBinsGravity.GlobalNActiveParticles);

  if(Sp.TimeBinsGravity.GlobalNActiveParticles > 0)
    {
      GravTree.DoEwald = 0;

#ifdef PMGRID
      GravTree.DoPM = 0; /* default value */
#endif

#if defined(PMGRID) && defined(PERIODIC) && \
    !defined(TREEPM_NOTIMESPLIT) /* classic TreePM in periodic box with time integration split */
      GravTree.DoEwald = 0;
      GravTree.DoPM    = TREE_ACTIVE_CUTTOFF_BASE_PM;

#ifdef PLACEHIGHRESREGION
      if(Sp.TimeBinsGravity.GlobalNActiveParticles > All.ActivePartFracForPMinsteadOfEwald * Sp.TotNumPart)
        GravTree.DoPM += TREE_DO_HIGHRES_PM + TREE_ACTIVE_CUTTOFF_HIGHRES_PM;
#endif

#else /* everything else */

#if defined(PMGRID) /* with PM acceleration (we force here TREEPM_NOTIMESPLIT to be set) */

      if(Sp.TimeBinsGravity.GlobalNActiveParticles > All.ActivePartFracForPMinsteadOfEwald * Sp.TotNumPart)
        {
          GravTree.DoEwald = 0;
#ifdef PLACEHIGHRESREGION
          GravTree.DoPM    = TREE_DO_BASE_PM + TREE_DO_HIGHRES_PM + TREE_ACTIVE_CUTTOFF_BASE_PM + TREE_ACTIVE_CUTTOFF_HIGHRES_PM;
#else
          GravTree.DoPM = TREE_DO_BASE_PM + TREE_ACTIVE_CUTTOFF_BASE_PM;
#endif
        }
      else
        {
          GravTree.DoPM    = 0;
#if defined(PERIODIC)
          GravTree.DoEwald = 1;
#endif
        }

#else /* here no PM acceleration is used */

#if defined(PERIODIC)
      GravTree.DoEwald = 1; /* periodic boundaries with Ewald summation */
#endif

#endif

#endif

#ifdef SECOND_ORDER_LPT_ICS
      if(All.Ti_Current == 0)
        second_order_ic_correction();
#endif

      if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0 && All.RelOpeningCriterionInUse == 0)
        {
          /* For the first timestep, we do one gravity calculation up front
           * with the Barnes & Hut Criterion to allow usage of relative opening
           * criterion with consistent accuracy.
           */
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
          gravity_long_range_force();
#endif
          gravity(timebin);

          gravity_set_oldacc(timebin);

          /* now may switch on relative opening criterion since we have an old acceleration */
          if(All.TypeOfOpeningCriterion == 1)
            All.RelOpeningCriterionInUse = 1;
        }

      gravity(timebin); /* computes gravity acceleration */

#ifdef FORCETEST
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
      if(timebin == All.HighestOccupiedGravTimeBin)
#endif
        {
          GravTree.treeallocate(Sp.NumPart, &Sp, &Domain);
          GravTest.gravity_forcetest(timebin);
          GravTree.treefree();

          if(FORCETEST >= 2.0) /* this is for a special test where we repeat the calculation */
            {
              for(int i = 0; i < (int)(FORCETEST)-1; i++)
                {
                  NgbTree.treefree();
                  Domain.domain_free();
                  Domain.domain_decomposition(STANDARD);
                  NgbTree.treeallocate(Sp.NumGas, &Sp, &Domain);
                  NgbTree.treebuild(Sp.NumGas, NULL);

                  gravity(timebin);

                  GravTree.treeallocate(Sp.NumPart, &Sp, &Domain);
                  GravTest.gravity_forcetest(timebin);
                  GravTree.treefree();
                }
            }

          if(FORCETEST > 1.0)
            endrun();
        }
#endif
    }

  mpi_printf("ACCEL: tree force computation done.\n");

  if(All.TimeLimitCPU == 0)
    endrun();
}

/*! \brief main driver routine of gravity tree/fmm force calculation
 *
 *  This routine handles the whole tree force calculation. First it
 *  builds a new force tree force_treebuild() every timestep. This tree is then
 *  used to calculate a new tree force for every active particle ( gravity_tree() ).
 *  The passed variable 'timebin' is only used to inform about the largest timebin
 *  of the particles in the list of active particles.
 *
 *  The tree will be constructed for the NActiveGravity particles listed in the
 *  array ActiveActiveGravityParticles[]
 */
void sim::gravity(int timebin)
{
  if(timebin == All.HighestOccupiedGravTimeBin)
    GravTree.MeasureCostFlag = 1;
  else
    GravTree.MeasureCostFlag = 0;

  /* let's first initialize the results */
  for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
    {
      int target = Sp.TimeBinsGravity.ActiveParticleList[i];
#ifdef EVALPOTENTIAL
      Sp.P[target].Potential = 0;
#endif
      for(int j = 0; j < 3; j++)
        Sp.P[target].GravAccel[j] = 0;

      if(GravTree.MeasureCostFlag)
        Sp.P[target].GravCost = 0;
    }

#ifdef SELFGRAVITY
  /* set new softening lengths on global steps to take into account possible cosmological time variation */
  if(timebin == All.HighestOccupiedGravTimeBin)
    GravTree.set_softenings();

#if defined(PMGRID) && (defined(TREEPM_NOTIMESPLIT) || defined(PLACEHIGHRESREGION))
  if((GravTree.DoPM & (TREE_DO_BASE_PM + TREE_DO_HIGHRES_PM)))
    gravity_pm(timebin);

#if defined(FORCETEST) && defined(PLACEHIGHRESREGION)
  for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
    {
      int target                = Sp.TimeBinsGravity.ActiveParticleList[i];
      Sp.P[target].PotentialHPM = All.G * Sp.P[target].Potential;
      for(int j = 0; j < 3; j++)
        Sp.P[target].GravAccelHPM[j] = All.G * Sp.P[target].GravAccel[j];
    }
#endif

#endif

#ifdef ALLOW_DIRECT_SUMMATION
  if(Sp.TimeBinsGravity.GlobalNActiveParticles < DIRECT_SUMMATION_THRESHOLD)
    {
      GravTree.gravity_direct(&Sp, &Domain, timebin);
    }
  else
#endif
    {
      GravTree.treeallocate(Sp.NumPart, &Sp, &Domain);

#ifdef HIERARCHICAL_GRAVITY
      GravTree.treebuild(Sp.TimeBinsGravity.NActiveParticles, Sp.TimeBinsGravity.ActiveParticleList);
#else
    GravTree.treebuild(Sp.NumPart, NULL);
#endif

#ifdef FMM
      GravTree.gravity_fmm(timebin);
#else
    GravTree.gravity_tree(timebin);
#endif

      GravTree.treefree();
    }

  /* now multplify with G and add things for comoving integration */
  gravity_comoving_factors(timebin);

#endif

#ifdef EXTERNALGRAVITY
  gravity_external();
#endif
}

void sim::gravity_set_oldacc(int timebin)
{
#ifdef HIERARCHICAL_GRAVITY
  if(timebin == All.HighestOccupiedGravTimeBin)
#endif
    {
      mpi_printf("GRAVTREE/FMM: Setting OldAcc!\n");

      particle_data *P = Sp.P;

      double ginv = 1 / All.G;

      for(int idx = 0; idx < Sp.TimeBinsGravity.NActiveParticles; idx++)
        {
          int target = Sp.TimeBinsGravity.ActiveParticleList[idx];
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
          double ax = P[target].GravAccel[0] + P[target].GravPM[0];
          double ay = P[target].GravAccel[1] + P[target].GravPM[1];
          double az = P[target].GravAccel[2] + P[target].GravPM[2];
#else
        double ax = P[target].GravAccel[0];
        double ay = P[target].GravAccel[1];
        double az = P[target].GravAccel[2];
#endif
          P[target].OldAcc = sqrt(ax * ax + ay * ay + az * az) * ginv;
        }
    }
}

void sim::gravity_comoving_factors(int timebin)
{
  particle_data *P = Sp.P;

#ifndef PERIODIC
  if(All.ComovingIntegrationOn)
    {
      /* here we carry out an integration in comoving coordinates but in a non-periodic space (i.e. the 'big sphere setup') */
      double fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
        {
          int target = Sp.TimeBinsGravity.ActiveParticleList[i];

          double pos[3];
          Sp.intpos_to_pos(P[target].IntPos, pos); /* converts the integer distance to floating point */

          for(int j = 0; j < 3; j++)
            P[target].GravAccel[j] += fac * pos[j];

#ifdef EVALPOTENTIAL
          double r2 = 0;
          for(int k = 0; k < 3; k++)
            r2 += pos[k] * pos[k];
          P[target].Potential -= fac * r2;
#endif
        }
    }
#endif

  /*  muliply by G */
  for(int idx = 0; idx < Sp.TimeBinsGravity.NActiveParticles; idx++)
    {
      int target = Sp.TimeBinsGravity.ActiveParticleList[idx];

      for(int j = 0; j < 3; j++)
        P[target].GravAccel[j] *= All.G;

#if defined(EVALPOTENTIAL) && defined(PMGRID) && defined(PERIODIC)
        /* To get correct zero point in potential for TreePM calculation, need to add this term,
         * because we cannot include the pi/(V*alpha^2) term in the correction potential in real space
         * since we only touch a restricted set of particles in the tree calculation
         */
#ifndef GRAVITY_TALLBOX
      /* note: for the tallbox,  the constant potential term is comuted as part of the long-range force and not the short-range force,
       * unlike in the ordinary periodic case */
      if(GravTree.DoPM)
        {
          double alpha = 0.5 / Sp.Asmth[0];
          P[target].Potential +=
              All.TotalMass * M_PI / (alpha * alpha * All.BoxSize * All.BoxSize * All.BoxSize) * (LONG_X * LONG_Y * LONG_Z);
        }
#endif
#endif

#if defined(EVALPOTENTIAL)
#ifndef FMM
        /* remove self-interaction */
#if NSOFTCLASSES > 1
      P[target].Potential += P[target].getMass() / (All.ForceSoftening[P[target].getSofteningClass()] / 2.8);
#else
      P[target].Potential += P[target].getMass() / (All.ForceSoftening[0] / 2.8);
#endif

#endif

#if defined(FMM) && defined(PERIODIC) && !defined(PMGRID)
      /* in FMM case, add in interaction with other images in periodic grid */
      P[target].Potential += P[target].getMass() * Ewald.ewald_gridlookup_origin_D0();
#endif

      P[target].Potential *= All.G;

#if defined(PMGRID) && !defined(FORCETEST_TESTFORCELAW)
      P[target].Potential += P[target].PM_Potential; /* add in long-range potential */
#endif
#endif

      if(All.ComovingIntegrationOn == 0 && All.OmegaLambda != 0)
        {
#ifdef PERIODIC
          Terminate(
              "You specified a periodic simulation in physical coordinates but with a non-zero cosmological constant - this can't be "
              "run");
#endif
          /* Finally, the following factor allows a computation of a cosmological simulation
               with vacuum energy in physical coordinates */

          double pos[3];
          Sp.intpos_to_pos(P[target].IntPos, pos); /* converts the integer distance to floating point */

          double fac = All.OmegaLambda * All.Hubble * All.Hubble;

          for(int j = 0; j < 3; j++)
            Sp.P[target].GravAccel[j] += fac * pos[j];

#ifdef EVALPOTENTIAL
          double r2 = 0;
          for(int k = 0; k < 3; k++)
            r2 += pos[k] * pos[k];
          P[target].Potential -= 0.5 * fac * r2;
#endif
        }
    }
}

#if defined(PMGRID) && (defined(TREEPM_NOTIMESPLIT) || defined(PLACEHIGHRESREGION))
void sim::gravity_pm(int timebin)
{
  double tstart = Logs.second();
  TIMER_START(CPU_PM_GRAVITY);
  mpi_printf("TREEPM: Starting PM part of force calculation. (timebin=%d)\n", timebin);

#if !defined(PERIODIC) || defined(PLACEHIGHRESREGION)
  PM.pm_init_regionsize();
#endif

  if((GravTree.DoPM & TREE_DO_BASE_PM))
    {
#ifdef PERIODIC
      PM.pmforce_periodic(LOW_MESH, NULL);
#else
      /* non periodic PM mesh */
      PM.pmforce_nonperiodic(LOW_MESH);
#endif
    }

#ifdef PLACEHIGHRESREGION
  if((GravTree.DoPM & TREE_DO_HIGHRES_PM))
    PM.pmforce_nonperiodic(HIGH_MESH);
#endif

  mpi_printf("TREEPM: Finished PM part of force calculation.\n");
  TIMER_STOP(CPU_PM_GRAVITY);
  double tend               = Logs.second();
  All.CPUForLastPMExecution = Logs.timediff(tstart, tend);
}
#endif

/*! \brief This function computes the long-range PM force for all the particles.
 *
 */
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
void sim::gravity_long_range_force(void)
{
#ifndef SELFGRAVITY
  return;
#endif

  double tstart = Logs.second();
  TIMER_START(CPU_PM_GRAVITY);

  for(int i = 0; i < Sp.NumPart; i++)
    {
      Sp.P[i].GravPM[0] = Sp.P[i].GravPM[1] = Sp.P[i].GravPM[2] = 0;
#ifdef EVALPOTENTIAL
      Sp.P[i].PM_Potential = 0;
#endif
    }

  PM.pmforce_periodic(0, NULL);

  /* multiply with the gravitational constant */
  for(int i = 0; i < Sp.NumPart; i++)
    {
      for(int j = 0; j < 3; j++)
        Sp.P[i].GravPM[j] *= All.G;
#ifdef EVALPOTENTIAL
      Sp.P[i].PM_Potential *= All.G;
#endif
    }

  TIMER_STOP(CPU_PM_GRAVITY);
  double tend               = Logs.second();
  All.CPUForLastPMExecution = Logs.timediff(tstart, tend);

  Sp.find_long_range_step_constraint();
}
#endif
