/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  timestep.cc
 *
 *  \brief routines for determining the timesteps of particles
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../data/simparticles.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"
#include "../time_integration/timestep.h"

/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */

void sim::find_hydro_timesteps(void)
{
#ifndef FORCE_EQUAL_TIMESTEPS

  All.set_cosmo_factors_for_current_time();

  NgbTree.tree_based_timesteps();

  TIMER_START(CPU_TIMELINE);

  Sp.assign_hydro_timesteps();

  TIMER_STOP(CPU_TIMELINE);

#endif
}

void simparticles::assign_hydro_timesteps(void)
{
  /* Now assign new timesteps for the hydro particles that are synchronized */
  for(int i = 0; i < TimeBinsHydro.NActiveParticles; i++)
    {
      int target = TimeBinsHydro.ActiveParticleList[i];
      if(P[target].getType() != 0)
        continue;

      if(TimeBinSynchronized[P[target].getTimeBinHydro()])
        {
          integertime ti_step = get_timestep_hydro(target);

          int bin;
          timebins_get_bin_and_do_validity_checks(ti_step, &bin, P[target].getTimeBinHydro());

#ifdef SELFGRAVITY
          /* we enforce that the hydro timestep is nested inside the gravity step */
          if(bin > P[target].TimeBinGrav)
            bin = P[target].TimeBinGrav;
#endif

          TimeBinsHydro.timebin_move_particle(target, P[target].getTimeBinHydro(), bin);

          P[target].setTimeBinHydro(bin);
        }
    }
}

#ifdef FORCE_EQUAL_TIMESTEPS
void sim::find_global_timesteps(void)
{
  All.set_cosmo_factors_for_current_time();

  NgbTree.tree_based_timesteps();

  TIMER_START(CPU_TIMELINE);

  integertime globTimeStep = TIMEBASE;

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  globTimeStep = Sp.get_timestep_pm();
#endif

#if defined(SELFGRAVITY) || defined(EXTERNALGRAVITY)
  for(int idx = 0; idx < Sp.TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = Sp.TimeBinsGravity.ActiveParticleList[idx];

      integertime ti_step = Sp.get_timestep_grav(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }
#endif

  for(int idx = 0; idx < Sp.TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = Sp.TimeBinsHydro.ActiveParticleList[idx];
      if(Sp.P[i].getType() != 0)
        continue;

      integertime ti_step = Sp.get_timestep_hydro(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &globTimeStep, &All.GlobalTimeStep, Communicator);
#else
  MPI_Allreduce(&globTimeStep, &All.GlobalTimeStep, 1, MPI_INT, MPI_MIN, Communicator);
#endif

  for(int idx = 0; idx < Sp.TimeBinsGravity.NActiveParticles; idx++)
    {
      int target = Sp.TimeBinsGravity.ActiveParticleList[idx];

      int bin;

      Sp.timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, Sp.P[target].TimeBinGrav);
      Sp.TimeBinsGravity.timebin_move_particle(target, Sp.P[target].TimeBinGrav, bin);
      Sp.P[target].TimeBinGrav = bin;
    }

  for(int idx = 0; idx < Sp.TimeBinsHydro.NActiveParticles; idx++)
    {
      int target = Sp.TimeBinsHydro.ActiveParticleList[idx];
      if(Sp.P[target].getType() != 0)
        continue;

      int bin;

      Sp.timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, Sp.P[target].getTimeBinHydro());
      Sp.TimeBinsHydro.timebin_move_particle(target, Sp.P[target].getTimeBinHydro(), bin);
#ifndef LEAN
      Sp.P[target].TimeBinHydro = bin;
#endif
    }

  TIMER_STOP(CPU_TIMELINE);
}
#endif

int simparticles::test_if_grav_timestep_is_too_large(int p, int bin)
{
  integertime ti_step_bin = bin ? (((integertime)1) << bin) : 0;

  integertime ti_step = get_timestep_grav(p);

  if(ti_step < ti_step_bin)
    return 1;
  else
    return 0;
}

/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
integertime simparticles::get_timestep_grav(int p /*!< particle index */)
{
  double ax = All.cf_a2inv * P[p].GravAccel[0];
  double ay = All.cf_a2inv * P[p].GravAccel[1];
  double az = All.cf_a2inv * P[p].GravAccel[2];

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  ax += All.cf_a2inv * P[p].GravPM[0];
  ay += All.cf_a2inv * P[p].GravPM[1];
  az += All.cf_a2inv * P[p].GravPM[2];
#endif

  double ac = sqrt(ax * ax + ay * ay + az * az); /* this is now the physical acceleration */

  if(ac == 0)
    ac = MIN_FLOAT_NUMBER;

    /* determine the "kinematic" timestep dt_grav, in physical units */
#if NSOFTCLASSES > 1
  double dt_grav = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].getSofteningClass()] / ac);
#else
  double dt_grav = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[0] / ac);
#endif

  double dt = dt_grav;

  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     All.cf_hubble_a=1.
   */
  dt *= All.cf_hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

#if defined(SELFGRAVITY) && defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  if(dt >= All.DtDisplacement)
    dt = All.DtDisplacement;
#endif

  if(dt < All.MinSizeTimestep)
    {
      dt = All.MinSizeTimestep;
#ifndef NO_STOP_BELOW_MINTIMESTEP
      Terminate(
          "Timestep wants to be below the limit MinSizeTimestep=%g\n"
          "Part-ID=%lld task=%d type=%d dtgrav=%g ac=%g soft=%g\n",
          All.MinSizeTimestep, (long long)P[p].ID.get(), ThisTask, P[p].getType(), dt, ac,
          All.SofteningTable[P[p].getSofteningClass()]);
#endif
    }

  integertime ti_step = (integertime)(dt / All.Timebase_interval);

  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      double pos[3];
      intpos_to_pos(P[p].IntPos, pos); /* converts the integer coordinates to floating point */

      Terminate(
          "\nError: A timestep of size zero was assigned on the integer timeline!\n"
          "We better stop.\n"
          "Task=%d Part-ID=%lld type=%d dt_grav=%g dt=%g tibase=%g ac=%g xyz=(%g|%g|%g) vel=(%g|%g|%g) tree=(%g|%g|%g) mass=%g\n\n",
          ThisTask, (long long)P[p].ID.get(), P[p].getType(), dt_grav, dt, All.Timebase_interval, ac, pos[0], pos[1], pos[2],
          P[p].Vel[0], P[p].Vel[1], P[p].Vel[2], P[p].GravAccel[0], P[p].GravAccel[1], P[p].GravAccel[2], P[p].getMass());

      myflush(stdout);
      Terminate("integer timestep reached zero");
    }

  return ti_step;
}

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
integertime simparticles::get_timestep_pm(void)
{
  integertime ti_step = TIMEBASE;
  while(ti_step > (All.DtDisplacement / All.Timebase_interval))
    ti_step >>= 1;

  if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep)) /* PM-timestep wants to increase */
    {
      int bin    = get_timestep_bin(ti_step);
      int binold = get_timestep_bin(All.PM_Ti_endstep - All.PM_Ti_begstep);

      while(TimeBinSynchronized[bin] == 0 && bin > binold) /* make sure the new step is synchronized */
        bin--;

      ti_step = bin ? (((integertime)1) << bin) : 0;
    }

  if(All.Ti_Current == TIMEBASE) /* we here finish the last timestep. */
    ti_step = 0;

  return ti_step;
}
#endif

integertime simparticles::get_timestep_hydro(int p /*!< particle index */)
{
  if(P[p].getType() != 0)
    Terminate("P[p].getType() != 0");

  double ax = All.cf_afac2 * SphP[p].HydroAccel[0];
  double ay = All.cf_afac2 * SphP[p].HydroAccel[1];
  double az = All.cf_afac2 * SphP[p].HydroAccel[2];

  ax += All.cf_a2inv * P[p].GravAccel[0];
  ay += All.cf_a2inv * P[p].GravAccel[1];
  az += All.cf_a2inv * P[p].GravAccel[2];

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  ax += All.cf_a2inv * P[p].GravPM[0];
  ay += All.cf_a2inv * P[p].GravPM[1];
  az += All.cf_a2inv * P[p].GravPM[2];
#endif

  double ac = sqrt(ax * ax + ay * ay + az * az); /* this is now the physical acceleration */

  if(ac == 0)
    ac = MIN_FLOAT_NUMBER;

  /* determine the "kinematic" timestep dt_grav, in physical units */
  double dt_kin = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * SphP[p].Hsml / ac);

  /* calculate local Courant timestep and treebased maximum timestep in physical units */

  double dt_courant = (All.cf_atime / All.cf_afac3) * All.CourantFac * 2.0 * SphP[p].Hsml / (SphP[p].MaxSignalVel + MIN_FLOAT_NUMBER);

  double dt_treebased = (All.cf_atime / All.cf_afac3) * SphP[p].CurrentMaxTiStep;

  /* calculate a timestep that restricts the rate at which the smoothing length may change,
   * in physical units
   */
  double dt_hsml = All.cf_atime2 * All.CourantFac * fabs(SphP[p].Hsml / (SphP[p].DtHsml + MIN_FLOAT_NUMBER));

  /* now take the smallest of these four criteria */
  double dt = dt_kin;
  if(dt > dt_courant)
    dt = dt_courant;
  if(dt > dt_treebased)
    dt = dt_treebased;
  if(dt > dt_hsml)
    dt = dt_hsml;

#ifdef STARFORMATION
  if(P[p].getType() == 0) /* to protect using a particle that has been turned into a star */
    {
      if(SphP[p].Sfr > 0)
        {
          double dt_sfr =
              0.1 * P[p].getMass() / (SphP[p].Sfr / ((All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR)));
          if(dt_sfr < dt)
            dt = dt_sfr;
        }
    }
#endif

  /* convert the physical timestep to dloga in the cosmological case.
   * Note: If comoving integration has not been selected, All.cf_hubble_a = 1.0.
   */
  dt *= All.cf_hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  if(dt >= All.DtDisplacement)
    dt = All.DtDisplacement;
#endif

  if(dt < All.MinSizeTimestep)
    {
      if(P[p].getType() == 0)
        Terminate(
            "Timestep wants to be below the limit MinSizeTimestep=%g\n"
            "Part-ID=%lld task=%d dtkin=%g dtcourant=%g ac=%g\n",
            All.MinSizeTimestep, (long long)P[p].ID.get(), ThisTask, dt_kin * All.cf_hubble_a, dt_courant * All.cf_hubble_a, ac);
      dt = All.MinSizeTimestep;
    }

  integertime ti_step = (integertime)(dt / All.Timebase_interval);

  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      double pos[3];
      intpos_to_pos(P[p].IntPos, pos); /* converts the integer coordinates to floating point */

      Terminate(
          "\nError: A timestep of size zero was assigned on the integer timeline!\n"
          "We better stop.\n"
          "Task=%d Part-ID=%lld type=%d dt=%g dtc=%g dt_kin=%g dt_treebased=%g dt_hsml=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) "
          "vel=(%g|%g|%g) "
          "tree=(%g|%g|%g) mass=%g  All.cf_hubble_a=%g\n\n",
          ThisTask, (long long)P[p].ID.get(), P[p].getType(), dt, dt_courant, dt_kin, dt_treebased, dt_hsml, All.Timebase_interval,
          (int)ti_step, ac, pos[0], pos[1], pos[2], P[p].Vel[0], P[p].Vel[1], P[p].Vel[2], P[p].GravAccel[0], P[p].GravAccel[1],
          P[p].GravAccel[2], P[p].getMass(), All.cf_hubble_a);
    }

  return ti_step;
}

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
void simparticles::find_long_range_step_constraint(void)
{
  double dtmin = MAX_DOUBLE_NUMBER;

  for(int p = 0; p < NumPart; p++)
    {
      if(P[p].getType() == 0)
        continue;

      /* calculate acceleration */
      double ax = All.cf_a2inv * P[p].GravPM[0];
      double ay = All.cf_a2inv * P[p].GravPM[1];
      double az = All.cf_a2inv * P[p].GravPM[2];

      double ac = sqrt(ax * ax + ay * ay + az * az); /* this is now the physical acceleration */

      if(ac < 1.0e-30)
        ac = 1.0e-30;

#if NSOFTCLASSES > 1
      double dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.ForceSoftening[P[p].getSofteningClass()] / 2.8 / ac);
#else
      double dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.ForceSoftening[0] / 2.8 / ac);
#endif
      dt *= All.cf_hubble_a;

      if(dt < dtmin)
        dtmin = dt;
    }

  dtmin *= 2.0; /* move it one timebin higher to prevent being too conservative */

  MPI_Allreduce(&dtmin, &All.DtDisplacement, 1, MPI_DOUBLE, MPI_MIN, Communicator);

  mpi_printf("TIMESTEPS: displacement time constraint: %g  (%g)\n", All.DtDisplacement, All.MaxSizeTimestep);

  if(All.DtDisplacement > All.MaxSizeTimestep)
    All.DtDisplacement = All.MaxSizeTimestep;
}
#endif

int simparticles::get_timestep_bin(integertime ti_step)
{
  int bin = -1;

  if(ti_step == 0)
    return 0;

  if(ti_step == 1)
    Terminate("time-step of integer size 1 not allowed\n");

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  return bin;
}

void TimeBinData::timebins_init(const char *name, int *maxPart)
{
  NActiveParticles   = 0;
  ActiveParticleList = 0;

  for(int i = 0; i < TIMEBINS; i++)
    {
      FirstInTimeBin[i] = -1;
      LastInTimeBin[i]  = -1;
    }

  NextInTimeBin = 0;
  PrevInTimeBin = 0;

  strncpy(Name, name, 99);
  Name[99] = 0;
  MaxPart  = maxPart;
}

void TimeBinData::timebins_allocate(void)
{
  char Identifier[200];
  Identifier[199] = 0;

  snprintf(Identifier, 199, "NextActiveParticle%s", Name);
  ActiveParticleList = (int *)Mem.mymalloc_movable(&ActiveParticleList, Identifier, *(MaxPart) * sizeof(int));

  snprintf(Identifier, 199, "NextInTimeBin%s", Name);
  NextInTimeBin = (int *)Mem.mymalloc_movable(&NextInTimeBin, Identifier, *(MaxPart) * sizeof(int));

  snprintf(Identifier, 199, "PrevInTimeBin%s", Name);
  PrevInTimeBin = (int *)Mem.mymalloc_movable(&PrevInTimeBin, Identifier, *(MaxPart) * sizeof(int));
}

void TimeBinData::timebins_free(void)
{
  Mem.myfree_movable(PrevInTimeBin);
  Mem.myfree_movable(NextInTimeBin);
  Mem.myfree_movable(ActiveParticleList);

  PrevInTimeBin      = NULL;
  NextInTimeBin      = NULL;
  ActiveParticleList = NULL;
}

void TimeBinData::timebins_reallocate(void)
{
  if(ActiveParticleList != NULL)
    {
      ActiveParticleList = (int *)Mem.myrealloc_movable(ActiveParticleList, *(MaxPart) * sizeof(int));
      NextInTimeBin      = (int *)Mem.myrealloc_movable(NextInTimeBin, *(MaxPart) * sizeof(int));
      PrevInTimeBin      = (int *)Mem.myrealloc_movable(PrevInTimeBin, *(MaxPart) * sizeof(int));
    }
}

void simparticles::timebins_get_bin_and_do_validity_checks(integertime ti_step, int *bin_new, int bin_old)
{
  /* make it a power 2 subdivision */
  integertime ti_min = TIMEBASE;
  while(ti_min > ti_step)
    ti_min >>= 1;
  ti_step = ti_min;

  /* get timestep bin */
  int bin = -1;

  if(ti_step == 0)
    bin = 0;

  if(ti_step == 1)
    Terminate("time-step of integer size 1 not allowed\n");

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  if(bin > bin_old) /* timestep wants to increase */
    {
      while(TimeBinSynchronized[bin] == 0 && bin > bin_old) /* make sure the new step is synchronized */
        bin--;

      ti_step = bin ? (((integertime)1) << bin) : 0;
    }

  if(All.Ti_Current >= TIMEBASE) /* we here finish the last timestep. */
    {
      ti_step = 0;
      bin     = 0;
    }

  if((TIMEBASE - All.Ti_Current) < ti_step) /* check that we don't run beyond the end */
    {
      Terminate("we are beyond the end of the timeline"); /* should not happen */
    }

  *bin_new = bin;
}

void TimeBinData::timebin_move_particle(int p, int timeBin_old, int timeBin_new)
{
  if(timeBin_old == timeBin_new)
    return;

  TimeBinCount[timeBin_old]--;

  int prev = PrevInTimeBin[p];
  int next = NextInTimeBin[p];

  if(FirstInTimeBin[timeBin_old] == p)
    FirstInTimeBin[timeBin_old] = next;
  if(LastInTimeBin[timeBin_old] == p)
    LastInTimeBin[timeBin_old] = prev;
  if(prev >= 0)
    NextInTimeBin[prev] = next;
  if(next >= 0)
    PrevInTimeBin[next] = prev;

  if(TimeBinCount[timeBin_new] > 0)
    {
      PrevInTimeBin[p]                          = LastInTimeBin[timeBin_new];
      NextInTimeBin[LastInTimeBin[timeBin_new]] = p;
      NextInTimeBin[p]                          = -1;
      LastInTimeBin[timeBin_new]                = p;
    }
  else
    {
      FirstInTimeBin[timeBin_new] = LastInTimeBin[timeBin_new] = p;
      PrevInTimeBin[p] = NextInTimeBin[p] = -1;
    }

  TimeBinCount[timeBin_new]++;
}

void TimeBinData::timebin_remove_particle(int idx, int bin)
{
  int p                   = ActiveParticleList[idx];
  ActiveParticleList[idx] = -1;

  TimeBinCount[bin]--;

  if(p >= 0)
    {
      int prev = PrevInTimeBin[p];
      int next = NextInTimeBin[p];

      if(prev >= 0)
        NextInTimeBin[prev] = next;
      if(next >= 0)
        PrevInTimeBin[next] = prev;

      if(FirstInTimeBin[bin] == p)
        FirstInTimeBin[bin] = next;
      if(LastInTimeBin[bin] == p)
        LastInTimeBin[bin] = prev;
    }
}

/* insert a particle into the timebin struct behind another already existing particle */
void TimeBinData::timebin_add_particle(int i_new, int i_old, int timeBin, int addToListOfActiveParticles)
{
  TimeBinCount[timeBin]++;

  if(i_old < 0)
    {
      /* if we don't have an existing particle to add if after, let's take the last one in this timebin */
      i_old = LastInTimeBin[timeBin];

      if(i_old < 0)
        {
          /* the timebin is empty at the moment, so just add the new particle */
          FirstInTimeBin[timeBin] = i_new;
          LastInTimeBin[timeBin]  = i_new;
          NextInTimeBin[i_new]    = -1;
          PrevInTimeBin[i_new]    = -1;
        }
    }

  if(i_old >= 0)
    {
      /* otherwise we added it already */
      PrevInTimeBin[i_new] = i_old;
      NextInTimeBin[i_new] = NextInTimeBin[i_old];
      if(NextInTimeBin[i_old] >= 0)
        PrevInTimeBin[NextInTimeBin[i_old]] = i_new;
      NextInTimeBin[i_old] = i_new;
      if(LastInTimeBin[timeBin] == i_old)
        LastInTimeBin[timeBin] = i_new;
    }

  if(addToListOfActiveParticles)
    {
      ActiveParticleList[NActiveParticles] = i_new;
      NActiveParticles++;
    }
}

void simparticles::timebin_cleanup_list_of_active_particles(void)
{
  for(int idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].ID.get() == 0 && P[i].getMass() == 0)
        {
          TimeBinsGravity.timebin_remove_particle(idx, P[i].TimeBinGrav);
        }
    }

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].ID.get() == 0 && P[i].getMass() == 0 && P[i].getType() == 0)
        {
          TimeBinsHydro.timebin_remove_particle(idx, P[i].getTimeBinHydro());
        }
    }
}

void TimeBinData::timebin_make_list_of_active_particles_up_to_timebin(int timebin)
{
  NActiveParticles = 0;
  for(int tbin = timebin; tbin >= 0; tbin--)
    timebin_add_particles_of_timebin_to_list_of_active_particles(tbin);
}

void TimeBinData::timebin_add_particles_of_timebin_to_list_of_active_particles(int timebin)
{
  for(int i = FirstInTimeBin[timebin]; i >= 0; i = NextInTimeBin[i])
    {
      ActiveParticleList[NActiveParticles] = i;
      NActiveParticles++;
    }
}
