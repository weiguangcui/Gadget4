/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  run.cc
 *
 *  \brief contains the basic simulation loop that iterates over timesteps
 */

// clang-format off
#include "gadgetconfig.h"
// clang-format on

#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../gravtree/gravtree.h"
#include "../io/io.h"
#include "../io/snap_io.h"
#include "../lightcone/lightcone_massmap_io.h"
#include "../lightcone/lightcone_particle_io.h"
#include "../logs/logs.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../ngbtree/ngbtree.h"
#include "../sort/parallel_sort.h"
#include "../system/system.h"

/*!
 * Main driver routine for advancing the simulation forward in time.
 * The loop terminates when the cpu-time limit is reached, when a `stop' file
 * is found in the output directory, or when the simulation ends because we
 * arrived at TimeMax.
 *
 * If the simulation is started from initial conditions, an initial domain
 * decomposition is performed, the gravitational forces are computed and
 * initial hydro forces are calculated.
 */
void sim::run(void)
{
#if defined(NGENIC_TEST) && defined(PERIODIC) && defined(PMGRID)
  snap_io Snap(&Sp, Communicator, All.SnapFormat);             /* get an I/O object */
  Snap.write_snapshot(All.SnapshotFileCount, NORMAL_SNAPSHOT); /* write snapshot file */
#if defined(POWERSPEC_ON_OUTPUT)
  PM.calculate_power_spectra(All.SnapshotFileCount);
#endif
  return;
#endif

  while(1) /* main loop over synchronization points */
    {
      /* store old time for logging purposes */
      All.TimeOld = All.Time;

      /* determine the next synchronization time at which we have active particles */
      integertime ti_next_kick_global = Sp.find_next_sync_point();

#ifdef OUTPUT_NON_SYNCHRONIZED_ALLOWED
      while(ti_next_kick_global > All.Ti_nextoutput && All.Ti_nextoutput >= 0)
        {
          All.Ti_Current = All.Ti_nextoutput;
          All.Time       = All.get_absolutetime_from_integertime(All.Ti_Current);
          All.set_cosmo_factors_for_current_time();

          Sp.drift_all_particles();
          create_snapshot_if_desired();
        }
#endif

      All.Ti_Current = ti_next_kick_global;
      All.Time       = All.get_absolutetime_from_integertime(All.Ti_Current);
      All.set_cosmo_factors_for_current_time();
      All.TimeStep = All.Time - All.TimeOld;

#ifdef LIGHTCONE
#ifdef LIGHTCONE_PARTICLES
      mpi_printf("LIGHTCONE_PARTICLES: Lp.NumPart=%d\n", Lp.NumPart);
#endif
#ifdef LIGHTCONE_MASSMAPS
      mpi_printf("LIGHTCONE_MASSMAPS:  Mp.NumPart=%d \n", Mp.NumPart);
#endif
#if defined(LIGHTCONE_MASSMAPS) || defined(LIGHTCONE_PARTICLES_GROUPS)
      Sp.drift_all_particles();  // we do this here to be able to avoid large buffer sizes, if needed multiple binning operations are
                                 // done
#endif
#endif

      /* mark the timebins that are active on this step */
      Sp.mark_active_timebins();

      /* create lists with the particles that are synchronized on this step */
      Sp.make_list_of_active_particles();

      /* produce some further log messages */
      Logs.output_log_messages();

      /* call functions that update certain 'extra' physics settings to new current time */
      set_non_standard_physics_for_current_time();

      /* for sufficiently large steps, carry out a new domain decomposition */
      if(All.HighestActiveTimeBin >= All.SmallestTimeBinWithDomainDecomposition)
        {
          NgbTree.treefree();
          Domain.domain_free();

          Sp.drift_all_particles();

#ifdef LIGHTCONE_PARTICLES
          LightCone.lightcone_clear_boxlist(All.Time);
#endif

#ifdef DEBUG_MD5
          Logs.log_debug_md5("C");
#endif
          Domain.domain_decomposition(STANDARD);

#ifdef DEBUG_MD5
          Logs.log_debug_md5("D");
#endif

          NgbTree.treeallocate(Sp.NumGas, &Sp, &Domain);
          NgbTree.treebuild(Sp.NumGas, NULL);
        }

      /* compute SPH densities and smoothing lengths for active SPH particles, and optionally those
       * accessed passively. This also creates the list of active hydro particles at this
       * synchronization point, which is stored in the list TimeBinsHydro.ActiveParticleList[].
       * This list is reused for the subsequent second and first hydro step. */
      NgbTree.compute_densities();

      /* if particles have increased their smoothing lengths, this is recorded in parent tree nodes */
      NgbTree.update_maxhsml();

      /* hydro-forces, second half-step. This will also update the predicted velocities/entropies with the new current ones */
      do_hydro_step_second_half();

      /* this does the closing gravity half-step for the timebins that end at the current synchronization point */
      do_gravity_step_second_half();

      /* do any extra physics, in a Strang-split way, for the timesteps that are finished */
      calculate_non_standard_physics_end_of_step();

#ifdef DEBUG_MD5
      Logs.log_debug_md5("A");
#endif

      Logs.compute_statistics();

      Logs.flush_everything();

#ifdef DEBUG_MD5
      Logs.log_debug_md5("BEFORE SNAP");
#endif
      create_snapshot_if_desired();

#ifdef DEBUG_MD5
      Logs.log_debug_md5("AFTER SNAP");
#endif

      if(All.Ti_Current >= TIMEBASE || All.Time > All.TimeMax) /* did we reached the final time? */
        {
          mpi_printf("\nFinal time=%g reached. Simulation ends.\n", All.TimeMax);

          /* make a snapshot at the final time in case none has been produced at this time yet */
          if(All.Ti_lastoutput != All.Ti_Current)
            {
              All.Ti_nextoutput = All.Ti_Current;
              create_snapshot_if_desired();
            }

          break;
        }

      /* kicks particles by half a gravity step */
      find_timesteps_and_do_gravity_step_first_half();

#ifdef DEBUG_MD5
      Logs.log_debug_md5("B");
#endif

      /* Find new hydro timesteps. This will not change the set of active hydro particles at this synchronization point,
       * but it can change how they are distributed over timebins. */
      find_hydro_timesteps();

      /* compute hydro-forces and apply momentum changes to interacting particle pairs for first half-steps */
      do_hydro_step_first_half();

      /* update the neighbor tree with the new velocities */
      NgbTree.update_velocities();

      /* output some CPU usage log-info (accounts for everything needed up to complete the previous timestep) */
      Logs.write_cpu_log();

#ifdef STOP_AFTER_STEP
      if(All.NumCurrentTiStep == STOP_AFTER_STEP)
        {
          mpi_printf("RUN: We have reached the timestep specified with STOP_AFTER_STEP and therefore stop.");
          endrun();
        }
#endif

      All.NumCurrentTiStep++;

      /* Check whether we should write a restart file */
      if(check_for_interruption_of_run())
        return;
    }

  restart Restart{Communicator};
  Restart.write(this); /* write a restart file at final time - can be used to continue simulation beyond final time */

  Logs.write_cpu_log(); /* output final cpu measurements */
}

/*! \brief calls extra modules after drift operator
 *
 * This routine is called after a new synchronization time has been determined.
 */
void sim::set_non_standard_physics_for_current_time(void)
{
#ifdef COOLING
  CoolSfr.IonizeParams(); /* set UV background for the current time */
#endif
}

/*! \brief calls extra modules at the end of the run loop
 *
 * The second gravitational half kick has already been applied to the
 * particles at this time, i.e. the particles at the sync-point have finished their regular timestep.
 *
 */
void sim::calculate_non_standard_physics_end_of_step(void)
{
#ifdef COOLING
#ifdef STARFORMATION
  CoolSfr.sfr_create_star_particles(&Sp);
  CoolSfr.cooling_and_starformation(&Sp);
#else
  CoolSfr.cooling_only(&Sp);
#endif
#endif

#ifdef MEASURE_TOTAL_MOMENTUM
  Logs.compute_total_momentum();
#endif
}

/*! \brief checks whether the run must interrupted
 *
 * The run is interrupted either if the stop file is present or,
 * if 85% of the CPU time are up. This routine also handles the
 * regular writing of restart files. The restart file is also
 * written if the restart file is present
 *
 * \return 1 if the run has to be interrupted, 0 otherwise
 */
int sim::check_for_interruption_of_run(void)
{
  /* Check whether we need to interrupt the run */
  int stopflag = 0;
  if(ThisTask == 0)
    {
      FILE *fd;
      char stopfname[MAXLEN_PATH_EXTRA];

      snprintf(stopfname, MAXLEN_PATH_EXTRA, "%sstop", All.OutputDir);
      if((fd = fopen(stopfname, "r"))) /* Is the stop-file present? If yes, interrupt the run. */
        {
          fclose(fd);
          printf("stop-file detected. stopping.\n");
          stopflag = 1;
          unlink(stopfname);
        }

      snprintf(stopfname, MAXLEN_PATH_EXTRA, "%srestart", All.OutputDir);
      if((fd = fopen(stopfname, "r"))) /* Is the restart-file present? If yes, write a user-requested restart file. */
        {
          fclose(fd);
          printf("restart-file detected. writing restart files.\n");
          stopflag = 3;
          unlink(stopfname);
        }

      if(Logs.CPUThisRun > 0.85 * All.TimeLimitCPU) /* are we running out of CPU-time ? If yes, interrupt run. */
        {
          printf("reaching time-limit. stopping.\n");
          stopflag = 2;
        }
    }

  MPI_Bcast(&stopflag, 1, MPI_INT, 0, Communicator);

  if(stopflag)
    {
      restart Restart{Communicator};
      Restart.write(this); /* write restart file */
      MPI_Barrier(Communicator);

      if(stopflag == 3)
        return 0;

      if(stopflag == 2 && ThisTask == 0)
        {
          FILE *fd;
          char contfname[MAXLEN_PATH_EXTRA];
          snprintf(contfname, MAXLEN_PATH_EXTRA, "%scont", All.OutputDir);
          if((fd = fopen(contfname, "w")))
            fclose(fd);
        }
      return 1;
    }

  /* is it time to write a regular restart-file? (for security) */
  if(ThisTask == 0)
    {
      if((Logs.CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
        {
          All.TimeLastRestartFile = Logs.CPUThisRun;
          stopflag                = 3;
        }
      else
        stopflag = 0;
    }

  MPI_Bcast(&stopflag, 1, MPI_INT, 0, Communicator);

  if(stopflag == 3)
    {
      restart Restart{Communicator};
      Restart.write(this); /* write an occasional restart file */
      stopflag = 0;
    }
  return 0;
}

/*! \brief Returns the next output time that is equal or larger than
 *  ti_curr
 *
 *  \param ti_curr current simulation time
 */
integertime sim::find_next_outputtime(integertime ti_curr)
{
  integertime ti;
  integertime ti_next = -1;
  double time;

  All.DumpFlag_nextoutput = 1;

  if(All.OutputListOn)
    {
      for(int i = 0; i < All.OutputListLength; i++)
        {
          time = All.OutputListTimes[i];

          if(time >= All.TimeBegin && time <= All.TimeMax)
            {
              if(All.ComovingIntegrationOn)
                ti = (integertime)(log(time / All.TimeBegin) / All.Timebase_interval);
              else
                ti = (integertime)((time - All.TimeBegin) / All.Timebase_interval);

#ifndef OUTPUT_NON_SYNCHRONIZED_ALLOWED
              /* We will now modify 'ti' to map it to the closest available output time according to the specified MaxSizeTimestep.
               * The real output time may hence deviate by  +/- 0.5*MaxSizeTimestep from the desired output time.
               */

              /* first, determine maximum output interval based on All.MaxSizeTimestep */
              integertime timax = (integertime)(All.MaxSizeTimestep / All.Timebase_interval);

              /* make it a power 2 subdivision */
              integertime ti_min = TIMEBASE;
              while(ti_min > timax)
                ti_min >>= 1;
              timax = ti_min;

              double multiplier = ti / ((double)timax);

              /* now round this to the nearest multiple of timax */
              ti = ((integertime)(multiplier + 0.5)) * timax;
#endif

              if(ti >= ti_curr)
                {
                  if(ti_next == -1)
                    {
                      ti_next                 = ti;
                      All.DumpFlag_nextoutput = All.OutputListFlag[i];
                    }

                  if(ti_next > ti)
                    {
                      ti_next                 = ti;
                      All.DumpFlag_nextoutput = All.OutputListFlag[i];
                    }
                }
            }
        }
    }
  else
    {
      if(All.ComovingIntegrationOn)
        {
          if(All.TimeBetSnapshot <= 1.0)
            Terminate("TimeBetSnapshot > 1.0 required for your simulation.\n");
        }
      else
        {
          if(All.TimeBetSnapshot <= 0.0)
            Terminate("TimeBetSnapshot > 0.0 required for your simulation.\n");
        }
      time     = All.TimeOfFirstSnapshot;
      int iter = 0;

      while(time < All.TimeBegin)
        {
          if(All.ComovingIntegrationOn)
            time *= All.TimeBetSnapshot;
          else
            time += All.TimeBetSnapshot;

          iter++;

          if(iter > 1000000)
            Terminate("Can't determine next output time.\n");
        }

      while(time <= All.TimeMax)
        {
          if(All.ComovingIntegrationOn)
            ti = (integertime)(log(time / All.TimeBegin) / All.Timebase_interval);
          else
            ti = (integertime)((time - All.TimeBegin) / All.Timebase_interval);

#ifndef OUTPUT_NON_SYNCHRONIZED_ALLOWED
          /* We will now modify 'ti' to map it to the closest available output time according to the specified MaxSizeTimestep.
           * The real output time may hence deviate by  +/- 0.5*MaxSizeTimestep from the desired output time.
           */

          /* first, determine maximum output interval based on All.MaxSizeTimestep */
          integertime timax = (integertime)(All.MaxSizeTimestep / All.Timebase_interval);

          /* make it a power 2 subdivision */
          integertime ti_min = TIMEBASE;
          while(ti_min > timax)
            ti_min >>= 1;
          timax = ti_min;

          double multiplier = ti / ((double)timax);

          /* now round this to the nearest multiple of timax */
          ti = ((integertime)(multiplier + 0.5)) * timax;
#endif

          if(ti >= ti_curr)
            {
              ti_next = ti;
              break;
            }

          if(All.ComovingIntegrationOn)
            time *= All.TimeBetSnapshot;
          else
            time += All.TimeBetSnapshot;

          iter++;

          if(iter > MAXITER)
            Terminate("Can't determine next output time.\n");
        }
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE; /* this will prevent any further output */

      mpi_printf("\nSNAPSHOT: There is no valid time for a further snapshot file.\n");
    }
  else
    {
      double next = All.get_absolutetime_from_integertime(ti_next);

      mpi_printf("\nSNAPSHOT: Setting next time for snapshot file to Time_next= %g  (DumpFlag=%d)\n\n", next, All.DumpFlag_nextoutput);
    }

  return ti_next;
}

/*! \brief Check if a snapshot should be saved
 *
 * This function checks whether a snapshot file or other kinds of output files,
 * such as a projection, should be saved at the current time-step.
 * If that is the case, the appropriate functions to produce the
 * desired file are called and the parameter controlling the output are updated
 * accordingly.
 */
void sim::create_snapshot_if_desired(void)
{
#if defined(LIGHTCONE_MASSMAPS)
  /* we may do this on partial timesteps since for massmaps we always drift all particles, i.e. the lightcone is complete up to
   * All.Time */
  LightCone.lightcone_massmap_flush(1);
#endif

#ifndef OUTPUT_NON_SYNCHRONIZED_ALLOWED
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) /* allow only top-level synchronization points */
#endif
    if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
      {
        for(int i = 0; i < Sp.NumPart; i++)
          if(Sp.P[i].Ti_Current != All.Ti_Current)
            Terminate("P[i].Ti_Current != All.Ti_Current");

#if defined(STARFORMATION) && defined(FOF)
        // do an extra domain decomposition here to make sure that there are no new stars among the block of gas particles
        NgbTree.treefree();
        Domain.domain_free();
        Domain.domain_decomposition(STANDARD);
        NgbTree.treeallocate(Sp.NumGas, &Sp, &Domain);
        NgbTree.treebuild(Sp.NumGas, NULL);
#endif

#ifndef OUTPUT_NON_SYNCHRONIZED_ALLOWED
        NgbTree.treefree();
        Sp.TimeBinsGravity.timebins_free();
        Sp.TimeBinsHydro.timebins_free();
#endif

#ifdef FOF
        mpi_printf("\nFOF: We shall first compute a group catalog for this snapshot file\n");

        /* this structure will hold auxiliary information for each particle, needed only during group finding */
        Sp.PS = (subfind_data *)Mem.mymalloc_movable(&Sp.PS, "PS", Sp.MaxPart * sizeof(subfind_data));
        memset(Sp.PS, 0, Sp.MaxPart * sizeof(subfind_data));

        /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
        for(int i = 0; i < Sp.NumPart; i++)
          {
            Sp.PS[i].OriginTask  = ThisTask;
            Sp.PS[i].OriginIndex = i;
          }

        fof<simparticles> FoF{Communicator, &Sp, &Domain};

        FoF.fof_fof(All.SnapshotFileCount, "fof", "groups", 0);

#if defined(MERGERTREE) && defined(SUBFIND)
        MergerTree.CurrTotNsubhalos = FoF.TotNsubhalos;
        MergerTree.CurrNsubhalos    = FoF.Nsubhalos;

        MergerTree.mergertree_determine_descendants_on_the_fly(All.SnapshotFileCount);

        MergerTree.PrevTotNsubhalos = FoF.TotNsubhalos;
        MergerTree.PrevNsubhalos    = FoF.Nsubhalos;

        for(int n = 0; n < Sp.NumPart; n++)
          {
            Sp.P[n].PrevSubhaloNr     = Sp.PS[n].SubhaloNr;
            Sp.P[n].PrevSizeOfSubhalo = Sp.PS[n].SizeOfSubhalo;
            Sp.P[n].PrevRankInSubhalo = Sp.PS[n].RankInSubhalo;

            if(Sp.P[n].PrevSubhaloNr.get() >= MergerTree.PrevTotNsubhalos && Sp.P[n].PrevSubhaloNr.get() != HALONR_MAX)
              Terminate("Sp.P[n].PrevSubhaloNr=%lld  MergerTree.PrevTotNsubhalos=%lld\n", (long long)Sp.P[n].PrevSubhaloNr.get(),
                        (long long)MergerTree.PrevTotNsubhalos);

            if(Sp.P[n].PrevSizeOfSubhalo.get() > 0 && Sp.P[n].PrevSubhaloNr.get() == HALONR_MAX)
              Terminate("Sp.P[n].PrevSizeOfSubhalo=%d  Sp.P[n].PrevSubhaloNr=%lld\n", (int)Sp.P[n].PrevSizeOfSubhalo.get(),
                        (long long)Sp.P[n].PrevSubhaloNr.get());
          }
#endif
#endif

        if(All.DumpFlag_nextoutput)
          {
            snap_io Snap(&Sp, Communicator, All.SnapFormat);             /* get an I/O object */
            Snap.write_snapshot(All.SnapshotFileCount, NORMAL_SNAPSHOT); /* write snapshot file */
          }

#ifdef SUBFIND_ORPHAN_TREATMENT
        {
          snap_io Snap(&Sp, Communicator, All.SnapFormat);
          Snap.write_snapshot(All.SnapshotFileCount, MOST_BOUND_PARTICLE_SNAPHOT); /* write special snapshot file */
        }
#endif

#ifdef FOF
        /* now revert from output order to the original order */
        for(int n = 0; n < Sp.NumPart; n++)
          {
            Sp.PS[n].TargetTask  = Sp.PS[n].OriginTask;
            Sp.PS[n].TargetIndex = Sp.PS[n].OriginIndex;
          }

        TIMER_START(CPU_FOF);

        Domain.particle_exchange_based_on_PS(Communicator);

        TIMER_STOP(CPU_FOF);

        Mem.myfree(Sp.PS);
#endif

#if defined(POWERSPEC_ON_OUTPUT) && defined(PERIODIC) && defined(PMGRID)
        PM.calculate_power_spectra(All.SnapshotFileCount);
#endif

        All.SnapshotFileCount++;
        All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);

#ifndef OUTPUT_NON_SYNCHRONIZED_ALLOWED
        Sp.TimeBinsHydro.timebins_allocate();
        Sp.TimeBinsGravity.timebins_allocate();

        /* we need to reconstruct the timebins here. Even though the particles are in the same place again,
         * it could have happened that Sp.P was reduced in size temporarily below NumPart on a certain task,
         * in which case the timebin data may have become invalid.
         */
        Sp.reconstruct_timebins();

        NgbTree.treeallocate(Sp.NumGas, &Sp, &Domain);
        NgbTree.treebuild(Sp.NumGas, NULL);
#endif
      }

#if defined(LIGHTCONE_PARTICLES)
  if(Lp.TestIfAboveFillFactor(std::min<int>(Lp.MaxPart, Sp.MaxPart)))
    {
#if defined(LIGHTCONE_PARTICLES_GROUPS) && defined(FOF)
      /* do this only on full timesteps if groups are calculated on lightcone */
      if(All.Ti_Current >= TIMEBASE || All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
        {
          mpi_printf("\nLIGHTCONE_PARTICLES_GROUPS: We shall first compute a group catalogue for the lightcone particles\n");

          /* assign unique IDs to Lp particles */

          int *numlist = (int *)Mem.mymalloc("numlist", Lp.NumPart * sizeof(int));

          MPI_Allgather(&Lp.NumPart, 1, MPI_INT, numlist, 1, MPI_INT, Communicator);

          long long newID = 1;
          for(int i = 0; i < ThisTask; i++)
            newID += numlist[i];

          for(int i = 0; i < Lp.NumPart; i++)
            Lp.P[i].ID.set(newID++);

          Mem.myfree(numlist);

          domain<lcparticles> LcDomain(Communicator, &Lp);

          LcDomain.domain_decomposition(STANDARD);

          /* this structure will hold auxiliary information for each particle, needed only during group finding */
          Lp.PS = (subfind_data *)Mem.mymalloc_movable(&Lp.PS, "PS", Lp.MaxPart * sizeof(subfind_data));
          memset(Lp.PS, 0, Lp.MaxPart * sizeof(subfind_data));

          /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
          for(int i = 0; i < Lp.NumPart; i++)
            {
              Lp.PS[i].OriginTask  = ThisTask;
              Lp.PS[i].OriginIndex = i;
            }

          fof<lcparticles> FoF{Communicator, &Lp, &LcDomain};

          double inner_distance = Driftfac.get_comoving_distance_for_scalefactor(All.Time);

          FoF.fof_fof(All.LightconeFileCount, "lc_fof", "lc_groups", inner_distance);

#endif

          {
#ifdef MERGERTREE
            MergerTree.Ntrees = 0;
            lightcone_particle_io Lcone(&Lp, &LightCone, &MergerTree, Communicator, All.SnapFormat); /* get an I/O object */
#else
        lightcone_particle_io Lcone(&Lp, &LightCone, Communicator, All.SnapFormat); /* get an I/O object */
#endif
            long long NumLP_tot = Lp.NumPart;
            MPI_Allreduce(MPI_IN_PLACE, &NumLP_tot, 1, MPI_LONG_LONG, MPI_SUM, Communicator);
            mpi_printf("\nLIGHTCONE: writing particle lightcone conesnap files #%d ... (NumLP_tot = %lld)\n", All.LightconeFileCount,
                       NumLP_tot);

            for(int i = 0; i < Lp.NumPart; i++)
              {
                double pos[3];
                Lp.signedintpos_to_pos((MySignedIntPosType *)Lp.P[i].IntPos, pos);
                vec2pix_ring(LIGHTCONE_ORDER_NSIDE, pos, &Lp.P[i].ipnest);
              }

#if !defined(LIGHTCONE_PARTICLES_GROUPS)
            /* let's now sort the lightcone_particle_data according to healpix pixel number */
            mycxxsort_parallel(Lp.P, Lp.P + Lp.NumPart, Lp.compare_ipnest, Communicator);
#endif

#if !defined(LIGHTCONE_PARTICLES_SKIP_SAVING)
            for(int conenr = 0; conenr < LightCone.Nlightcones; conenr++)
              Lcone.lightcone_save(All.LightconeFileCount, conenr, false);
#endif

            mpi_printf("LIGHTCONE: done with writing files.\n");

            All.LightconeFileCount++;
            /* let's put this in brackets so that the object's destructor will be called already here */
          }

#if defined(LIGHTCONE_PARTICLES_GROUPS) && defined(FOF)
          /* now revert from output order to the original order */
          for(int n = 0; n < Lp.NumPart; n++)
            {
              Lp.PS[n].TargetTask  = Lp.PS[n].OriginTask;
              Lp.PS[n].TargetIndex = Lp.PS[n].OriginIndex;
            }

          TIMER_START(CPU_FOF);

          LcDomain.particle_exchange_based_on_PS(Communicator);

          Mem.myfree(Lp.PS);

          LcDomain.domain_free();

          TIMER_STOP(CPU_FOF);

          int ncount[2] = {0, 0};
          long long nsum[2];
          for(int n = 0; n < Lp.NumPart; n++)
            {
              if(Lp.P[n].getFlagSaveDistance())
                {
                  Lp.P[n--] = Lp.P[--Lp.NumPart];
                  ncount[0]++;
                }
              else
                {
                  ncount[1]++;
                }
            }

          sumup_large_ints(2, ncount, nsum, Communicator);
          mpi_printf("LIGHTCONE_PARTICLES_GROUPS: We could store %lld particles from the buffer, but had to keep %lld\n", nsum[0],
                     nsum[1]);
        }
#else
      Lp.NumPart = 0;
#endif

      if(Lp.MaxPart > LIGHTCONE_ALLOC_FAC * Sp.MaxPart + 1 && Lp.NumPart < LIGHTCONE_ALLOC_FAC * Sp.MaxPart)
        Lp.reallocate_memory_maxpart(LIGHTCONE_ALLOC_FAC * Sp.MaxPart);
    }
#endif
}
