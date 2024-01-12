/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  kicks.cc
 *
 *  \brief drives gravitational and hydrodynamical force calculations and applies corresponding kicks to particles
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../domain/domain.h"
#include "../gravtree/gravtree.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"
#include "../time_integration/timestep.h"

/*! \brief performs the first half step kick operator for the gravity
 *
 * This function applies a half step kick similar to do_gravity_step_second_half().
 * If we are on a PM step the kick due to the particle mesh's long range gravity
 * is applied first. Afterwards the short range kick due to the tree force is added.
 */
void sim::find_timesteps_and_do_gravity_step_first_half(void)
{
  TIMER_START(CPU_DRIFTS);

  All.set_cosmo_factors_for_current_time();

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  if(All.PM_Ti_endstep == All.Ti_Current) /* need to do long-range kick */
    {
      /* first, determine the new PM timestep */
      integertime ti_step = Sp.get_timestep_pm();

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      integertime tstart = All.PM_Ti_begstep;
      integertime tend   = tstart + ti_step / 2;

      double dt_gravkick;

      if(All.ComovingIntegrationOn)
        dt_gravkick = Driftfac.get_gravkick_factor(tstart, tend);
      else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;

      for(int i = 0; i < Sp.NumPart; i++)
        {
          for(int j = 0; j < 3; j++)
            Sp.P[i].Vel[j] += Sp.P[i].GravPM[j] * dt_gravkick;
        }
    }
#endif

  Sp.TimeBinsGravity.timebin_make_list_of_active_particles_up_to_timebin(All.HighestSynchronizedTimeBin);
  sumup_large_ints(1, &Sp.TimeBinsGravity.NActiveParticles, &Sp.TimeBinsGravity.GlobalNActiveParticles, Communicator);

#ifdef FORCE_EQUAL_TIMESTEPS
  find_global_timesteps();
#endif

#ifdef HIERARCHICAL_GRAVITY
  /* First, move all active particles to the highest allowed timestep for this synchronization time.
   * They will then cascade down to smaller timesteps as needed.
   */
  for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
    {
      int target = Sp.TimeBinsGravity.ActiveParticleList[i];
      int bin    = All.HighestSynchronizedTimeBin;
      int binold = Sp.P[target].TimeBinGrav;

      Sp.TimeBinsGravity.timebin_move_particle(target, binold, bin);
      Sp.P[target].TimeBinGrav = bin;
    }

  long long Previous_GlobalNActiveGravity = Sp.TimeBinsGravity.GlobalNActiveParticles;

  double dt_gravsum = 0;

  int bin_highest_occupied = 0;

  /* go over all timebins */
  for(int timebin = All.HighestSynchronizedTimeBin; timebin >= 0; timebin--)
    {
      Sp.TimeBinsGravity.NActiveParticles = 0;

      Sp.TimeBinsGravity.timebin_add_particles_of_timebin_to_list_of_active_particles(timebin);
      sumup_large_ints(1, &Sp.TimeBinsGravity.NActiveParticles, &Sp.TimeBinsGravity.GlobalNActiveParticles, Communicator);

      if(Sp.TimeBinsGravity.GlobalNActiveParticles == 0) /* we are done at this point */
        break;

      /* calculate gravity for all active particles */
      if(Sp.TimeBinsGravity.GlobalNActiveParticles != Previous_GlobalNActiveGravity)
        {
          TIMER_STOP(CPU_DRIFTS);

          compute_grav_accelerations(timebin);

          TIMER_START(CPU_DRIFTS);
        }

      /* now check whether the current timestep should be reduced */

      int nfine = 0;
      for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
        {
          int target = Sp.TimeBinsGravity.ActiveParticleList[i];
          int binold = Sp.P[target].TimeBinGrav;

          if(Sp.test_if_grav_timestep_is_too_large(target, binold))
            nfine++;
        }

      long long nfine_tot;
      sumup_large_ints(1, &nfine, &nfine_tot, Communicator);

      int push_down_flag = 0;
      if(Sp.TimeBinsGravity.GlobalNActiveParticles == Sp.TotNumPart && nfine_tot < Sp.TimeBinsGravity.GlobalNActiveParticles &&
         nfine_tot > 0.33 * Sp.TimeBinsGravity.GlobalNActiveParticles)
        {
          mpi_printf(
              "KICKS: We reduce the highest occupied timestep by pushing %lld particles on timebin=%d down in timestep (fraction "
              "wanting lower bin is: %g)\n",
              Sp.TimeBinsGravity.GlobalNActiveParticles - nfine_tot, timebin,
              ((double)nfine_tot) / Sp.TimeBinsGravity.GlobalNActiveParticles);

          push_down_flag = 1;
        }

      for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
        {
          int target = Sp.TimeBinsGravity.ActiveParticleList[i];
          int binold = Sp.P[target].TimeBinGrav;

          if(push_down_flag || Sp.test_if_grav_timestep_is_too_large(target, binold))
            {
              int bin = binold - 1;

              if(bin == 0)
                {
                  Sp.print_particle_info(target);
                  Terminate("timestep too small");
                }

              Sp.TimeBinsGravity.timebin_move_particle(target, binold, bin);
              Sp.P[target].TimeBinGrav = bin;
            }
          else if(binold > bin_highest_occupied)
            bin_highest_occupied = binold;
        }

      if(All.HighestOccupiedGravTimeBin == 0) /* this will only be the case in the very first step */
        {
          MPI_Allreduce(&bin_highest_occupied, &All.HighestOccupiedGravTimeBin, 1, MPI_INT, MPI_MAX, Communicator);
          if(All.HighestOccupiedTimeBin > 0)
            mpi_printf("KICKS: Special Start-up: All.HighestOccupiedGravTimeBin=%d\n", All.HighestOccupiedGravTimeBin);
        }

      if(Sp.TimeBinsGravity.GlobalNActiveParticles)
        {
          integertime ti_step = timebin ? (((integertime)1) << timebin) : 0;

          integertime tstart = All.Ti_begstep[timebin]; /* beginning of step */
          integertime tend   = tstart + ti_step / 2;    /* midpoint of step */

          double dt_gravkick;

          if(All.ComovingIntegrationOn)
            dt_gravkick = Driftfac.get_gravkick_factor(tstart, tend);
          else
            dt_gravkick = (tend - tstart) * All.Timebase_interval;

          double dt_save = dt_gravkick;

          if(timebin < All.HighestSynchronizedTimeBin)
            {
              ti_step = (timebin + 1) ? (((integertime)1) << (timebin + 1)) : 0;

              tstart = All.Ti_begstep[timebin + 1]; /* beginning of step */
              tend   = tstart + ti_step / 2;        /* midpoint of step */

              if(All.ComovingIntegrationOn)
                dt_gravkick -= Driftfac.get_gravkick_factor(tstart, tend);
              else
                dt_gravkick -= (tend - tstart) * All.Timebase_interval;
            }

          dt_gravsum += dt_gravkick;

          mpi_printf("KICKS: 1st gravity for hierarchical timebin=%d:  %lld particles   dt_gravkick=%g  %g %g\n", timebin,
                     Sp.TimeBinsGravity.GlobalNActiveParticles, dt_gravkick, dt_gravsum, dt_save);

          for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
            {
              int target = Sp.TimeBinsGravity.ActiveParticleList[i];

              for(int j = 0; j < 3; j++)
                Sp.P[target].Vel[j] += Sp.P[target].GravAccel[j] * dt_gravkick;
            }
        }

      Previous_GlobalNActiveGravity = Sp.TimeBinsGravity.GlobalNActiveParticles;
    }

#else

  mpi_printf("KICKS: 1st gravity for highest active timebin=%d:  particles %lld\n", All.HighestActiveTimeBin,
             Sp.TimeBinsGravity.GlobalNActiveParticles);

  for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
    {
      int target          = Sp.TimeBinsGravity.ActiveParticleList[i];

#ifndef FORCE_EQUAL_TIMESTEPS
      integertime ti_step = Sp.get_timestep_grav(target);
      int timebin;

      Sp.timebins_get_bin_and_do_validity_checks(ti_step, &timebin, Sp.P[target].TimeBinGrav);

      ti_step = timebin ? (((integertime)1) << timebin) : 0;

      Sp.TimeBinsGravity.timebin_move_particle(target, Sp.P[target].TimeBinGrav, timebin);
      Sp.P[target].TimeBinGrav = timebin;
#else
      int timebin         = Sp.P[target].TimeBinGrav;
      integertime ti_step = timebin ? (((integertime)1) << timebin) : 0;
#endif

      integertime tstart = All.Ti_begstep[timebin]; /* beginning of step */
      integertime tend   = tstart + ti_step / 2;    /* midpoint of step */

      double dt_gravkick;

      if(All.ComovingIntegrationOn)
        dt_gravkick = Driftfac.get_gravkick_factor(tstart, tend);
      else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;

      for(int j = 0; j < 3; j++)
        Sp.P[target].Vel[j] += Sp.P[target].GravAccel[j] * dt_gravkick;
    }

#endif

  TIMER_STOP(CPU_DRIFTS);
}

/*! \brief performs the second gravity half step kick operator
 *
 * This function applies a half step kick similar to do_gravity_step_first_half().
 * First the short range kick due to the tree force is added. If we are on a PM step the kick
 * due to the particle mesh's long range gravity is applied too. In both cases
 * the momentum and energy for Sph particles is updated.
 */
void sim::do_gravity_step_second_half(void)
{
  TIMER_START(CPU_DRIFTS);

  char fullmark[8];

  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    snprintf(fullmark, 8, "(*)");
  else
    fullmark[0] = 0;

  if(ThisTask == 0)
    {
      fprintf(Logs.FdTimings, "\nStep%s: %d, t: %g, dt: %g, highest active timebin: %d  (lowest active: %d, highest occupied: %d)\n",
              fullmark, All.NumCurrentTiStep, All.Time, All.TimeStep, All.HighestActiveTimeBin, All.LowestActiveTimeBin,
              All.HighestOccupiedTimeBin);

      fprintf(Logs.FdDensity, "\nStep%s: %d, t: %g, dt: %g, highest active timebin: %d  (lowest active: %d, highest occupied: %d)\n",
              fullmark, All.NumCurrentTiStep, All.Time, All.TimeStep, All.HighestActiveTimeBin, All.LowestActiveTimeBin,
              All.HighestOccupiedTimeBin);

      fprintf(Logs.FdHydro, "\nStep%s: %d, t: %g, dt: %g, highest active timebin: %d  (lowest active: %d, highest occupied: %d)\n",
              fullmark, All.NumCurrentTiStep, All.Time, All.TimeStep, All.HighestActiveTimeBin, All.LowestActiveTimeBin,
              All.HighestOccupiedTimeBin);
    }

  double dt_gravkick;

#ifdef HIERARCHICAL_GRAVITY
  /* go over all timebins, in inverse sequence so that we end up getting the cumulative force at the end */
  for(int timebin = 0; timebin <= All.HighestActiveTimeBin; timebin++)
    {
      if(Sp.TimeBinSynchronized[timebin])
        {
          /* need to make all timebins below the current one active */
          Sp.TimeBinsGravity.timebin_make_list_of_active_particles_up_to_timebin(timebin);
          sumup_large_ints(1, &Sp.TimeBinsGravity.NActiveParticles, &Sp.TimeBinsGravity.GlobalNActiveParticles, Communicator);

          if(Sp.TimeBinsGravity.GlobalNActiveParticles)
            {
              /* calculate gravity for all active particles */

              TIMER_STOP(CPU_DRIFTS);

              compute_grav_accelerations(timebin);

              TIMER_START(CPU_DRIFTS);

              mpi_printf("KICKS: 2nd gravity for hierarchical timebin=%d:  particles %lld\n", timebin,
                         Sp.TimeBinsGravity.GlobalNActiveParticles);

              integertime ti_step = timebin ? (((integertime)1) << timebin) : 0;

              integertime tend = All.Ti_begstep[timebin]; /* end of step (Note: All.Ti_begstep[] has already been advanced for the next
                                                             step at this point)   */
              integertime tstart = tend - ti_step / 2;    /* midpoint of step */

              if(All.ComovingIntegrationOn)
                dt_gravkick = Driftfac.get_gravkick_factor(tstart, tend);
              else
                dt_gravkick = (tend - tstart) * All.Timebase_interval;

              if(timebin < All.HighestActiveTimeBin)
                {
                  ti_step = (timebin + 1) ? (((integertime)1) << (timebin + 1)) : 0;

                  tend = All.Ti_begstep[timebin + 1]; /* end of step (Note: All.Ti_begstep[] has already been advanced for the next
                                                         step at this point)   */
                  tstart = tend - ti_step / 2;        /* midpoint of step */

                  if(All.ComovingIntegrationOn)
                    dt_gravkick -= Driftfac.get_gravkick_factor(tstart, tend);
                  else
                    dt_gravkick -= (tend - tstart) * All.Timebase_interval;
                }

              for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
                {
                  int target = Sp.TimeBinsGravity.ActiveParticleList[i];

                  for(int j = 0; j < 3; j++)
                    Sp.P[target].Vel[j] += Sp.P[target].GravAccel[j] * dt_gravkick;

                  if(Sp.P[target].getType() == 0 && All.HighestOccupiedGravTimeBin == timebin)
                    {
                      for(int j = 0; j < 3; j++)
                        {
                          Sp.SphP[target].VelPred[j]       = Sp.P[target].Vel[j];
                          Sp.SphP[target].FullGravAccel[j] = Sp.P[target].GravAccel[j];
                        }
                    }
                }
            }
        }
    }
#else

  Sp.TimeBinsGravity.timebin_make_list_of_active_particles_up_to_timebin(All.HighestActiveTimeBin);
  sumup_large_ints(1, &Sp.TimeBinsGravity.NActiveParticles, &Sp.TimeBinsGravity.GlobalNActiveParticles, Communicator);

  if(Sp.TimeBinsGravity.GlobalNActiveParticles)
    {
      TIMER_STOP(CPU_DRIFTS);

      /* calculate gravity for all active particles */
      compute_grav_accelerations(All.HighestActiveTimeBin);

      TIMER_START(CPU_DRIFTS);

      mpi_printf("KICKS: 2nd gravity for highest active timebin=%d:  particles %lld\n", All.HighestActiveTimeBin,
                 Sp.TimeBinsGravity.GlobalNActiveParticles);

      for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
        {
          int target  = Sp.TimeBinsGravity.ActiveParticleList[i];
          int timebin = Sp.P[target].TimeBinGrav;

          integertime ti_step = (timebin) ? (((integertime)1) << (timebin)) : 0;
          integertime tend    = All.Ti_Current;
          integertime tstart  = tend - ti_step / 2; /* midpoint of step */

          if(All.ComovingIntegrationOn)
            dt_gravkick = Driftfac.get_gravkick_factor(tstart, tend);
          else
            dt_gravkick = (tend - tstart) * All.Timebase_interval;

          for(int j = 0; j < 3; j++)
            Sp.P[target].Vel[j] += Sp.P[target].GravAccel[j] * dt_gravkick;

          if(Sp.P[target].getType() == 0)
            {
              for(int j = 0; j < 3; j++)
                Sp.SphP[target].VelPred[j] = Sp.P[target].Vel[j];
            }
        }
    }

#endif

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
  if(All.PM_Ti_endstep == All.Ti_Current) /* need to do long-range kick */
    {
      TIMER_STOP(CPU_DRIFTS);

      gravity_long_range_force();

      TIMER_START(CPU_DRIFTS);

      integertime ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
      integertime tstart  = All.PM_Ti_begstep + ti_step / 2;
      integertime tend    = tstart + ti_step / 2;

      if(All.ComovingIntegrationOn)
        dt_gravkick = Driftfac.get_gravkick_factor(tstart, tend);
      else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;

      for(int i = 0; i < Sp.NumPart; i++)
        for(int j = 0; j < 3; j++)
          Sp.P[i].Vel[j] += Sp.P[i].GravPM[j] * dt_gravkick;

      for(int i = 0; i < Sp.NumGas; i++)
        if(Sp.P[i].getType() == 0)
          for(int j = 0; j < 3; j++)
            Sp.SphP[i].VelPred[j] = Sp.P[i].Vel[j];

      gravity_set_oldacc(All.HighestActiveTimeBin);
    }
#else
  gravity_set_oldacc(All.HighestActiveTimeBin);
#endif

  TIMER_STOP(CPU_DRIFTS);
}

void sim::do_hydro_step_first_half(void)
{
  if(All.NumCurrentTiStep == 0) /* special domain decomposition now that we know the timebins of both gravity and hydro */
    {
      Sp.mark_active_timebins();

      NgbTree.treefree();

      Domain.domain_free();
      Domain.domain_decomposition(STANDARD);

      NgbTree.treeallocate(Sp.NumGas, &Sp, &Domain);
      NgbTree.treebuild(Sp.NumGas, NULL);
    }

  /* now we can calculate the hydro forces */
  hydro_force(FIRST_HALF_STEP); /* computes hydrodynamical accelerations and rate of chnange of entropy,
                                   and applies this where appropriate directly (half step kicks)  */
}

void sim::do_hydro_step_second_half(void)
{
  /* now we can calculate the hydro forces */
  hydro_force(SECOND_HALF_STEP); /* computes hydrodynamical accelerations and change rate of entropy,
                                    and applies this where appropriate directly (half step kicks)  */
}

/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 */
void sim::hydro_force(int step_indicator)
{
  if(Sp.TimeBinsHydro.GlobalNActiveParticles == 0)
    return;

  /* Create list of targets. */
  int *targetlist = (int *)Mem.mymalloc("targetlist", Sp.NumGas * sizeof(int));

  struct old_hydro_accel
  {
    MyFloat HydroAccel[3];
  };
  old_hydro_accel *Old = NULL;

  if(step_indicator == SECOND_HALF_STEP)
    Old = (old_hydro_accel *)Mem.mymalloc("Old", Sp.TimeBinsHydro.NActiveParticles * sizeof(old_hydro_accel));

  int Nhydroforces = 0;

  for(int i = 0; i < Sp.TimeBinsHydro.NActiveParticles; i++)
    {
      int target = Sp.TimeBinsHydro.ActiveParticleList[i];

      if(target < 0 || target >= Sp.NumGas)
        Terminate("target=%d i=%d\n", target, i);

      if(step_indicator == SECOND_HALF_STEP)
        {
          Old[i].HydroAccel[0] = Sp.SphP[target].HydroAccel[0];
          Old[i].HydroAccel[1] = Sp.SphP[target].HydroAccel[1];
          Old[i].HydroAccel[2] = Sp.SphP[target].HydroAccel[2];
        }

      targetlist[Nhydroforces++] = target;
    }

#ifdef REUSE_HYDRO_ACCELERATIONS_FROM_PREVIOUS_STEP
  /* in this case, the forces for the first hydro step are simply kept */
  if(step_indicator != FIRST_HALF_STEP)
#endif
    {
      NgbTree.hydro_forces_determine(Nhydroforces, targetlist);
    }

  /* let's now do the hydrodynamical kicks */
  for(int i = 0; i < Nhydroforces; i++)
    {
      int target = Sp.TimeBinsHydro.ActiveParticleList[i];

      int timebin = Sp.P[target].getTimeBinHydro();

      integertime tstart, tend, ti_step = timebin ? (((integertime)1) << timebin) : 0;

      if(step_indicator == SECOND_HALF_STEP)
        {
          tend   = All.Ti_Current;
          tstart = tend - ti_step / 2; /* midpoint of step */
        }
      else
        {
          tstart = All.Ti_Current;
          tend   = tstart + ti_step / 2; /* midpoint of step */
        }

      double dt_hydrokick, dt_entr = (tend - tstart) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dt_hydrokick = Driftfac.get_hydrokick_factor(tstart, tend);
      else
        dt_hydrokick = dt_entr;

      Sp.SphP[target].Entropy += Sp.SphP[target].DtEntropy * dt_entr;

      Sp.P[target].Vel[0] += Sp.SphP[target].HydroAccel[0] * dt_hydrokick;
      Sp.P[target].Vel[1] += Sp.SphP[target].HydroAccel[1] * dt_hydrokick;
      Sp.P[target].Vel[2] += Sp.SphP[target].HydroAccel[2] * dt_hydrokick;

      if(step_indicator == SECOND_HALF_STEP)
        {
          Sp.SphP[target].EntropyPred = Sp.SphP[target].Entropy;
#ifdef PRESSURE_ENTROPY_SPH
          Sp.SphP[target].EntropyToInvGammaPred = pow(Sp.SphP[target].EntropyPred, 1.0 / GAMMA);
#endif
          Sp.SphP[target].set_thermodynamic_variables();

          Sp.SphP[target].VelPred[0] += (Sp.SphP[target].HydroAccel[0] - Old[i].HydroAccel[0]) * dt_hydrokick;
          Sp.SphP[target].VelPred[1] += (Sp.SphP[target].HydroAccel[1] - Old[i].HydroAccel[1]) * dt_hydrokick;
          Sp.SphP[target].VelPred[2] += (Sp.SphP[target].HydroAccel[2] - Old[i].HydroAccel[2]) * dt_hydrokick;

          /* note: if there is no gravity, we should instead set VelPred = Vel (if this is not done anymore in the gravity
           * routine)
           */
        }
    }

  if(step_indicator == SECOND_HALF_STEP)
    Mem.myfree(Old);

  Mem.myfree(targetlist);
}
