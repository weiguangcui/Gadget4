/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  predict.cc
 *
 *  \brief find the next sync point, drift particles forward in time, and (re)build the timebin lists
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/simparticles.h"
#include "../lightcone/lightcone.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"
#include "../time_integration/timestep.h"

/*
 * It counts the number of particles in each timebin and updates the
 * linked lists containing the particles of each time bin. Afterwards the
 * linked list of active particles is updated by make_list_of_active_particles().
 *
 * The linked lists for each timebin are stored in #FirstInTimeBin[], #LastInTimeBin[],
 * #PrevInTimeBin[] and #NextInTimeBin[]. The counters of particles per timebin are
 * #TimeBinCount and #TimeBinCountSph.
 */
void simparticles::reconstruct_timebins(void)
{
  TIMER_START(CPU_TIMELINE);

  for(int bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsHydro.TimeBinCount[bin]   = 0;
      TimeBinsHydro.FirstInTimeBin[bin] = -1;
      TimeBinsHydro.LastInTimeBin[bin]  = -1;

      TimeBinsGravity.TimeBinCount[bin]   = 0;
      TimeBinsGravity.FirstInTimeBin[bin] = -1;
      TimeBinsGravity.LastInTimeBin[bin]  = -1;

#ifdef STARFORMATION
      TimeBinSfr[bin] = 0;
#endif
    }

  for(int i = 0; i < NumPart; i++)
    {
      int bin = P[i].TimeBinGrav;

      if(TimeBinsGravity.TimeBinCount[bin] > 0)
        {
          TimeBinsGravity.PrevInTimeBin[i]                                  = TimeBinsGravity.LastInTimeBin[bin];
          TimeBinsGravity.NextInTimeBin[i]                                  = -1;
          TimeBinsGravity.NextInTimeBin[TimeBinsGravity.LastInTimeBin[bin]] = i;
          TimeBinsGravity.LastInTimeBin[bin]                                = i;
        }
      else
        {
          TimeBinsGravity.FirstInTimeBin[bin] = TimeBinsGravity.LastInTimeBin[bin] = i;
          TimeBinsGravity.PrevInTimeBin[i] = TimeBinsGravity.NextInTimeBin[i] = -1;
        }

      TimeBinsGravity.TimeBinCount[bin]++;

      if(P[i].getType() == 0)
        {
          bin = P[i].getTimeBinHydro();

          if(TimeBinsHydro.TimeBinCount[bin] > 0)
            {
              TimeBinsHydro.PrevInTimeBin[i]                                = TimeBinsHydro.LastInTimeBin[bin];
              TimeBinsHydro.NextInTimeBin[i]                                = -1;
              TimeBinsHydro.NextInTimeBin[TimeBinsHydro.LastInTimeBin[bin]] = i;
              TimeBinsHydro.LastInTimeBin[bin]                              = i;
            }
          else
            {
              TimeBinsHydro.FirstInTimeBin[bin] = TimeBinsHydro.LastInTimeBin[bin] = i;
              TimeBinsHydro.PrevInTimeBin[i] = TimeBinsHydro.NextInTimeBin[i] = -1;
            }

          TimeBinsHydro.TimeBinCount[bin]++;

#ifdef STARFORMATION
          TimeBinSfr[bin] += SphP[i].Sfr;
#endif
        }
    }

  make_list_of_active_particles();

  TIMER_STOP(CPU_TIMELINE);
}

/*! \brief This function finds the next synchronization point of the system.
 * (i.e. the earliest point of time any of the particles needs a force
 * computation).
 *
 * This function drifts all particles, including inactive particles to the
 * next sync point. This is done by drift_particles(). Afterwards the linked
 * list of active particles is updated to the new sync point by
 * make_list_of_active_particles(). Particles become active/inactive here.
 */
integertime simparticles::find_next_sync_point(void)
{
  /* find the next kick time */
  integertime ti_next_kick = TIMEBASE;

  for(int n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinsGravity.TimeBinCount[n] || TimeBinsHydro.TimeBinCount[n])
        {
          integertime ti_next_for_bin;
          if(n > 0)
            {
              integertime dt_bin = (((integertime)1) << n);
              ti_next_for_bin    = (All.Ti_Current / dt_bin) * dt_bin + dt_bin; /* next kick time for this timebin */
            }
          else
            {
              ti_next_for_bin = All.Ti_Current;
            }

          if(ti_next_for_bin < ti_next_kick)
            ti_next_kick = ti_next_for_bin;
        }
    }

  integertime ti_next_kick_global;
#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &ti_next_kick, &ti_next_kick_global, Communicator);
#else
  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, Communicator);
#endif

  return ti_next_kick_global;
}

void simparticles::mark_active_timebins(void)
{
  int lowest_active_bin = TIMEBINS, highest_active_bin = 0;
  int lowest_occupied_bin = TIMEBINS, highest_occupied_bin = 0;
  int lowest_occupied_gravity_bin = TIMEBINS, highest_occupied_gravity_bin = 0;
  int highest_synchronized_bin = 0;
  int nsynchronized_gravity = 0, nsynchronized_hydro = 0;

  /* mark the bins that will be synchronized/active */

  for(int n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinsGravity.TimeBinCount[n])
        {
          if(highest_occupied_gravity_bin < n)
            highest_occupied_gravity_bin = n;

          if(lowest_occupied_gravity_bin > n)
            lowest_occupied_gravity_bin = n;
        }

      int active = TimeBinsHydro.TimeBinCount[n] + TimeBinsGravity.TimeBinCount[n];

      if(active)
        {
          if(highest_occupied_bin < n)
            highest_occupied_bin = n;

          if(lowest_occupied_bin > n)
            lowest_occupied_bin = n;
        }

      integertime dt_bin = (((integertime)1) << n);

      if((All.Ti_Current % dt_bin) == 0)
        {
          TimeBinSynchronized[n] = 1;
          All.Ti_begstep[n]      = All.Ti_Current;

          nsynchronized_gravity += TimeBinsGravity.TimeBinCount[n];
          nsynchronized_hydro += TimeBinsHydro.TimeBinCount[n];

          if(highest_synchronized_bin < n)
            highest_synchronized_bin = n;

          if(active)
            {
              if(highest_active_bin < n)
                highest_active_bin = n;

              if(lowest_active_bin > n)
                lowest_active_bin = n;
            }
        }
      else
        TimeBinSynchronized[n] = 0;
    }

  int lowest_in[3], lowest_out[3];
  lowest_in[0] = lowest_occupied_bin;
  lowest_in[1] = lowest_occupied_gravity_bin;
  lowest_in[2] = lowest_active_bin;
  MPI_Allreduce(lowest_in, lowest_out, 3, MPI_INT, MPI_MIN, Communicator);
  All.LowestOccupiedTimeBin     = lowest_out[0];
  All.LowestOccupiedGravTimeBin = lowest_out[1];
  All.LowestActiveTimeBin       = lowest_out[2];

  int highest_in[4], highest_out[4];
  highest_in[0] = highest_occupied_bin;
  highest_in[1] = highest_occupied_gravity_bin;
  highest_in[2] = highest_active_bin;
  highest_in[3] = highest_synchronized_bin;
  MPI_Allreduce(highest_in, highest_out, 4, MPI_INT, MPI_MAX, Communicator);
  All.HighestOccupiedTimeBin     = highest_out[0];
  All.HighestOccupiedGravTimeBin = highest_out[1];
  All.HighestActiveTimeBin       = highest_out[2];
  All.HighestSynchronizedTimeBin = highest_out[3];

  /* note: the lowest synchronized bin is always 1 */

  int input_ints[2 + 2 * TIMEBINS];
  long long output_longs[2 + 2 * TIMEBINS];

  input_ints[0] = nsynchronized_hydro;
  input_ints[1] = nsynchronized_gravity;
  memcpy(input_ints + 2, TimeBinsGravity.TimeBinCount, TIMEBINS * sizeof(int));
  memcpy(input_ints + 2 + TIMEBINS, TimeBinsHydro.TimeBinCount, TIMEBINS * sizeof(int));

  sumup_large_ints(2 + 2 * TIMEBINS, input_ints, output_longs, Communicator);

  All.GlobalNSynchronizedHydro   = output_longs[0];
  All.GlobalNSynchronizedGravity = output_longs[1];
  long long *tot_count_grav      = output_longs + 2;
  long long *tot_count_sph       = output_longs + 2 + TIMEBINS;

  long long tot_grav = 0, tot_sph = 0;

  for(int n = 0; n < TIMEBINS; n++)
    {
      tot_grav += tot_count_grav[n];
      tot_sph += tot_count_sph[n];

      if(n > 0)
        {
          tot_count_grav[n] += tot_count_grav[n - 1];
          tot_count_sph[n] += tot_count_sph[n - 1];
        }
    }

  All.SmallestTimeBinWithDomainDecomposition = All.HighestOccupiedTimeBin;

  for(int n = All.HighestOccupiedTimeBin; n >= All.LowestOccupiedTimeBin; n--)
    {
      if(tot_count_grav[n] > All.ActivePartFracForNewDomainDecomp * tot_grav ||
         tot_count_sph[n] > All.ActivePartFracForNewDomainDecomp * tot_sph)
        All.SmallestTimeBinWithDomainDecomposition = n;
    }
}

void simparticles::drift_all_particles(void)
{
  TIMER_START(CPU_DRIFTS);

  for(int i = 0; i < NumPart; i++)
    {
#ifdef LIGHTCONE_MASSMAPS
      int flag = drift_particle(&P[i], &SphP[i], All.Ti_Current);

      if(flag)
        {
          MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, Communicator);
          LightCone->lightcone_massmap_flush(0);
        }
#else
      drift_particle(&P[i], &SphP[i], All.Ti_Current);
#endif
    }

#ifdef LIGHTCONE_MASSMAPS
  int flag = 0;
  do
    {
      flag = 0;
      MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, Communicator);
      if(flag)
        LightCone->lightcone_massmap_flush(0);
    }
  while(flag);
#endif

  TIMER_STOP(CPU_DRIFTS);
}

/*! \brief This function drifts a particle i to time1
 *
 * @param time_previous current time
 * @param time1 time to which particles get drifted
 */
int simparticles::drift_particle(particle_data *P, sph_particle_data *SphP, integertime time1, bool ignore_light_cone)
{
  int buffer_full_flag = 0;
#ifndef LEAN
  while(P->access.test_and_set(std::memory_order_acquire))
    {
      // acquire spin lock
    }
#endif

  integertime time0 = P->Ti_Current.load(std::memory_order_acquire);

  if(time1 == time0)
    {
#ifndef LEAN
      P->access.clear(std::memory_order_release);
#endif
      return buffer_full_flag;
    }

  if(time1 < time0)
    Terminate("no prediction into past allowed: time0=%lld time1=%lld\n", (long long)time0, (long long)time1);

  double dt_drift;

  if(All.ComovingIntegrationOn)
    dt_drift = Driftfac.get_drift_factor(time0, time1);
  else
    dt_drift = (time1 - time0) * All.Timebase_interval;

#ifdef LIGHTCONE
  if(ignore_light_cone == false)
    buffer_full_flag = LightCone->lightcone_test_for_particle_addition(P, time0, time1, dt_drift);
#endif

  double posdiff[3];
  for(int j = 0; j < 3; j++)
    posdiff[j] = P->Vel[j] * dt_drift;

  MyIntPosType delta[3];
  pos_to_signedintpos(posdiff, (MySignedIntPosType *)delta);

  for(int j = 0; j < 3; j++)
    P->IntPos[j] += delta[j];

  constrain_intpos(P->IntPos); /* will only do something if we have a stretched box */

  if(P->getType() == 0)
    {
      double dt_hydrokick, dt_entr, dt_gravkick;

      if(All.ComovingIntegrationOn)
        {
          dt_entr      = (time1 - time0) * All.Timebase_interval;
          dt_hydrokick = Driftfac.get_hydrokick_factor(time0, time1);
          dt_gravkick  = Driftfac.get_gravkick_factor(time0, time1);
        }
      else
        dt_gravkick = dt_entr = dt_hydrokick = dt_drift;

      for(int j = 0; j < 3; j++)
        {
          SphP->VelPred[j] += SphP->HydroAccel[j] * dt_hydrokick;
#ifdef HIERARCHICAL_GRAVITY
          SphP->VelPred[j] += SphP->FullGravAccel[j] * dt_gravkick;
#else
          SphP->VelPred[j] += P->GravAccel[j] * dt_gravkick;
#endif
#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
          SphP->VelPred[j] += P->GravPM[j] * dt_gravkick;
#endif
        }

      SphP->EntropyPred += SphP->DtEntropy * dt_entr;

      SphP->Density += SphP->DtDensity * dt_drift;

      SphP->Hsml += SphP->DtHsml * dt_drift;

#ifdef PRESSURE_ENTROPY_SPH
      SphP->PressureSphDensity += SphP->DtPressureSphDensity * dt_drift;

      SphP->EntropyToInvGammaPred = pow(SphP->EntropyPred, 1.0 / GAMMA);
#endif

      SphP->set_thermodynamic_variables();
    }

  P->Ti_Current = time1;

#ifndef LEAN
  P->access.clear(std::memory_order_release);
#endif

  return buffer_full_flag;
}

void simparticles::make_list_of_active_particles(void)
{
  TIMER_START(CPU_DRIFTS);

  TimeBinsHydro.NActiveParticles = 0;

  for(int n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(int i = TimeBinsHydro.FirstInTimeBin[n]; i >= 0; i = TimeBinsHydro.NextInTimeBin[i])
            {
              if(P[i].getType() == 0)
                {
                  if(P[i].getTimeBinHydro() != n)
                    Terminate("P[i].TimeBinHydro=%d != timebin=%d", P[i].getTimeBinHydro(), n);

                  if(P[i].Ti_Current.load(std::memory_order_acquire) != All.Ti_Current)
                    drift_particle(&P[i], &SphP[i], All.Ti_Current);

                  TimeBinsHydro.ActiveParticleList[TimeBinsHydro.NActiveParticles++] = i;
                }
            }
        }
    }

  TimeBinsGravity.NActiveParticles = 0;

  for(int n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(int i = TimeBinsGravity.FirstInTimeBin[n]; i >= 0; i = TimeBinsGravity.NextInTimeBin[i])
            {
              if(P[i].Ti_Current.load(std::memory_order_acquire) != All.Ti_Current)
                drift_particle(&P[i], &SphP[i], All.Ti_Current);

              TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles++] = i;
            }
        }
    }

  int in[2] = {TimeBinsGravity.NActiveParticles, TimeBinsHydro.NActiveParticles};
  long long out[2];

  sumup_large_ints(2, in, out, Communicator);

  TimeBinsGravity.GlobalNActiveParticles = out[0];
  TimeBinsHydro.GlobalNActiveParticles   = out[1];

  TIMER_STOP(CPU_DRIFTS);
}
