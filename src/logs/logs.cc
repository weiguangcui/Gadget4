/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file logs.cc
 *
 *  \brief routines for log-file handling
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../logs/logs.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/*! \brief Open files for logging.
 *
 *   This function opens various log-files that report on the status and
 *   performance of the simulation. Upon restart, the code will append to
 *   these files.
 */
void logs::open_logfiles(void)
{
  char mode[2], buf[MAXLEN_PATH_EXTRA];

  if(All.RestartFlag == RST_BEGIN)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  if(ThisTask != 0) /* only the root processors writes to the log files */
    return;

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "cpu.txt");
  if(!(FdCPU = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "cpu.csv");
  if(!(FdCPUCSV = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "info.txt");
  if(!(FdInfo = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "energy.txt");
  if(!(FdEnergy = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "timings.txt");
  if(!(FdTimings = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "density.txt");
  if(!(FdDensity = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "hydro.txt");
  if(!(FdHydro = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "timebins.txt");
  if(!(FdTimebin = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "domain.txt");
  if(!(FdDomain = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

#ifdef MEASURE_TOTAL_MOMENTUM
  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "momentum.txt");
  if(!(FdMomentum = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);
#endif

#ifdef FORCETEST
  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "forcetest.txt");
  if(!(FdForceTest = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  fclose(FdForceTest);
#endif

#ifdef DEBUG_MD5
  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "debug_md5.txt");
  if(!(FdDebug = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);
#endif

  fprintf(FdBalance, "\n");

  fprintf(FdCPUCSV, "STEP, TIME, CPUS, MULTIPLEDOMAIN, HIGHESTTIMEBIN, ");

  for(int i = 0; i < CPU_LAST; i++)
    {
      if(Timer_data[i].symb != 0 && Timer_data[i].symbImbal != 0)
        fprintf(FdBalance, "%-20s = '%c' / '%c'\n", Timer_data[i].longname, Timer_data[i].symb, Timer_data[i].symbImbal);

      fprintf(FdCPUCSV, "%s1, %s2, %s3, ", Timer_data[i].shortname, Timer_data[i].shortname, Timer_data[i].shortname);
    }
  fprintf(FdBalance, "\n");

  fprintf(FdCPUCSV, "\n");

#ifdef STARFORMATION
  snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);
#endif
}

int logs::flush_everything(void)
{
#ifndef REDUCE_FLUSH
  return 0;
#else
  if(ThisTask == 0)
    {
      if((CPUThisRun - All.FlushLast) < All.FlushCpuTimeDiff)
        {
          return 0;
        }
      else
        {
          All.FlushLast = CPUThisRun;
        }
    }
  else
    {
      return 0;
    }
#endif

  mpi_printf("Flushing...\n");

  fflush(FdDomain);
  fflush(FdTimings);
  fflush(FdInfo);
  fflush(FdTimebin);
  fflush(FdBalance);
  fflush(FdCPU);
  fflush(FdEnergy);
  fflush(FdCPUCSV);

#ifdef STARFORMATION
  fflush(FdSfr);
#endif

  return 1;
}

/*! \brief Write the FdInfo and FdTimeBin files.
 *
 * At each time step this function writes on to two log-files.
 * In FdInfo, it just lists the timesteps that have been done, while in
 * FdTimeBin it outputs information about the active and occupied time-bins.
 */
void logs::output_log_messages(void)
{
  TIMER_START(CPU_LOGS);

  long long tot_count_grav[TIMEBINS], tot_count_sph[TIMEBINS];
  sumup_large_ints(TIMEBINS, Sp->TimeBinsGravity.TimeBinCount, tot_count_grav, Communicator);
  sumup_large_ints(TIMEBINS, Sp->TimeBinsHydro.TimeBinCount, tot_count_sph, Communicator);

  Mem.report_detailed_memory_usage_of_largest_task();

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
        {
          double z = 1.0 / (All.Time) - 1;
          fprintf(FdInfo, "\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n",
                  All.NumCurrentTiStep, All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep),
                  All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);
          printf("\n\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n",
                 All.NumCurrentTiStep, All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep),
                 All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);
          fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep, All.Time, z,
                  All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
          myflush(FdInfo);
        }
      else
        {
          fprintf(FdInfo, "\nSync-Point %d, Time: %g, Systemstep: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n", All.NumCurrentTiStep,
                  All.Time, All.TimeStep, All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);
          printf("\n\nSync-Point %d, Time: %g, Systemstep: %g, Nsync-grv: %10llu, Nsync-hyd: %10llu\n", All.NumCurrentTiStep, All.Time,
                 All.TimeStep, All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro);
          fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
          myflush(FdInfo);
        }

      long long tot_cumulative_grav[TIMEBINS], tot_cumulative_sph[TIMEBINS];
      tot_cumulative_grav[0] = tot_count_grav[0];
      tot_cumulative_sph[0]  = tot_count_sph[0];

      for(int i = 1; i < TIMEBINS; i++)
        {
          tot_cumulative_grav[i] = tot_count_grav[i] + tot_cumulative_grav[i - 1];
          tot_cumulative_sph[i]  = tot_count_sph[i] + tot_cumulative_sph[i - 1];
        }

      double avg_CPU_TimeBin[TIMEBINS];

      for(int i = 0; i < TIMEBINS; i++)
        {
          double sum = 0;
          if(tot_count_sph[i] > 0 || tot_count_grav[i] > 0)
            for(int j = 0; j < All.CPU_TimeBinCountMeasurements[i]; j++)
              sum += All.CPU_TimeBinMeasurements[i][j];

          if(All.CPU_TimeBinCountMeasurements[i] && (tot_count_sph[i] > 0 || tot_count_grav[i] > 0))
            avg_CPU_TimeBin[i] = sum / All.CPU_TimeBinCountMeasurements[i];
          else
            avg_CPU_TimeBin[i] = 0;
        }

      int weight = 1;
      double sum = 0;
      double frac_CPU_TimeBin[TIMEBINS];

      for(int i = All.HighestOccupiedTimeBin; i >= 0; i--, weight *= 2)
        {
          int corr_weight;

          if(weight > 1)
            corr_weight = weight / 2;
          else
            corr_weight = weight;

          frac_CPU_TimeBin[i] = corr_weight * avg_CPU_TimeBin[i];
          sum += frac_CPU_TimeBin[i];
        }

      for(int i = All.HighestOccupiedTimeBin; i >= 0; i--)
        {
          if(sum)
            frac_CPU_TimeBin[i] /= sum;
        }

      fprintf(FdTimebin,
              "Occupied timebins: gravity         sph          dt              cumul-grav   cumul-sph A D    avg-time  cpu-frac\n");

      long long tot_grav = 0, tot_sph = 0;
      for(int i = TIMEBINS - 1; i >= 0; i--)
        if(tot_count_sph[i] > 0 || tot_count_grav[i] > 0)
          {
            fprintf(
                FdTimebin, " %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu  %10llu %c %c  %10.2f    %5.1f%%\n",
                Sp->TimeBinSynchronized[i] ? 'X' : ' ', i, tot_count_grav[i], tot_count_sph[i],
                i > 0 ? (((integertime)1) << i) * All.Timebase_interval : 0.0, tot_cumulative_grav[i], tot_cumulative_sph[i],
                (i == All.HighestActiveTimeBin) ? '<' : ' ',
                (All.HighestActiveTimeBin >= All.SmallestTimeBinWithDomainDecomposition && i == All.HighestActiveTimeBin) ? '*' : ' ',
                avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);

            if(Sp->TimeBinSynchronized[i])
              {
                tot_grav += tot_count_grav[i];
                tot_sph += tot_count_sph[i];
              }
          }
      fprintf(FdTimebin, "               ------------------------\n");
#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
      if(All.PM_Ti_endstep == All.Ti_Current)
        {
          fprintf(FdTimebin, "PM-Step. Total: %10llu  %10llu\n", tot_grav, tot_sph);
        }
      else
#endif
        {
          fprintf(FdTimebin, "Total active:   %10llu  %10llu\n", tot_grav, tot_sph);
        }
      fprintf(FdTimebin, "\n");
      myflush(FdTimebin);
    }

  TIMER_STOP(CPU_LOGS);
}

void logs::init_cpu_log(simparticles *Sp_ptr)
{
  Sp = Sp_ptr;

  for(int i = 0; i < CPU_LAST; i++)
    {
      if(Timer_data[i].parent >= 0)
        Timer_data[i].depth = Timer_data[Timer_data[i].parent].depth + 1;
      else
        Timer_data[i].depth = 0;
    }

  for(int i = 0; i < CPU_LAST; i++)
    {
      CPU_Sum[i]  = 0.;
      CPU_Step[i] = 0.;
    }
  CPUThisRun = 0.;

  TimerStackPos = 0;
  TimerStack[0] = CPU_MISC;

  WallclockTime = Logs.second();
  StartOfRun    = Logs.second();
}

/*! \brief Write the FdBalance and FdCPU files.
 *
 * At each time step this function writes on to two log-files.
 * In FdBalance, it outputs in a graphical way the amount of
 * time spent in the various parts of the code, while
 * in FdCPU it writes information about the cpu-time consumption
 * of the various modules.
 */
void logs::write_cpu_log(void)
{
  TIMER_START(CPU_LOGS);

  double local_total = 0;
  for(int i = 0; i < CPU_LAST; i++)
    local_total += CPU_Step[i];

  double max_total = 0;
  MPI_Reduce(&local_total, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, Communicator);

  double max_CPU_Step[CPU_LAST], avg_CPU_Step[CPU_LAST];
  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_LAST, MPI_DOUBLE, MPI_MAX, 0, Communicator);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_LAST, MPI_DOUBLE, MPI_SUM, 0, Communicator);

  if(ThisTask == 0)
    {
      double summed_CPU_Step[CPU_LAST];

      /* sum up cpu items into groups */
      for(int i = 0; i < CPU_LAST; i++)
        summed_CPU_Step[i] = avg_CPU_Step[i];

      for(int i = CPU_LAST - 1; i > CPU_ALL; i--)
        if(Timer_data[i].parent >= 0)
          summed_CPU_Step[Timer_data[i].parent] += summed_CPU_Step[i];

      /* calc averages, update CPU_Sum */
      double avg_total = 0;
      for(int i = 0; i < CPU_LAST; i++)
        {
          avg_CPU_Step[i] /= NTask;
          avg_total += avg_CPU_Step[i];

          summed_CPU_Step[i] /= NTask;
          CPU_Sum[i] += summed_CPU_Step[i];
        }

      /* create balance.txt string */
      char cpu_String[CPU_STRING_LEN + 1];
      put_symbol(cpu_String, 0., 1.0, '-');

      double tsum = 0.0;
      for(int i = 1; i < CPU_LAST; i++)
        {
          if(max_CPU_Step[i] > 0 && Timer_data[i].symb != 0 && Timer_data[i].symbImbal != 0)
            {
              double t0 = tsum;
              double t1 = tsum + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
              put_symbol(cpu_String, t0 / avg_total, t1 / avg_total, Timer_data[i].symb);
              tsum += t1 - t0;

              t0 = tsum;
              t1 = tsum + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
              put_symbol(cpu_String, t0 / avg_total, t1 / avg_total, Timer_data[i].symbImbal);
              tsum += t1 - t0;
            }
        }

      // put_symbol(cpu_String, tsum / max_total, 1.0, '-');
      fprintf(FdBalance, "Step=%7d  sec=%10.3f Nsync-grv=%10llu Nsync-hyd=%10llu  %s\n", All.NumCurrentTiStep, max_total,
              All.GlobalNSynchronizedGravity, All.GlobalNSynchronizedHydro, cpu_String);
      myflush(FdBalance);

      if(All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin] == NUMBER_OF_MEASUREMENTS_TO_RECORD)
        {
          All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]--;
          memmove(&All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][0], &All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][1],
                  (NUMBER_OF_MEASUREMENTS_TO_RECORD - 1) * sizeof(double));
        }

      All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]++] = max_total;

      fprintf(FdCPUCSV, "%d, %g, %d, %d, ", All.NumCurrentTiStep, All.Time, NTask, All.HighestActiveTimeBin);

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d, HighestActiveTimeBin: %d\n", All.NumCurrentTiStep, All.Time, NTask,
              All.HighestActiveTimeBin);
      fprintf(FdCPU, "                          diff                 cumulative\n");

      for(int i = 0; i < CPU_LAST; i++)
        {
          fprintf(FdCPU, "%*s%*s%10.2f  %5.1f%% %10.2f  %5.1f%%\n", 2 * Timer_data[i].depth, "", -20 + 2 * Timer_data[i].depth,
                  Timer_data[i].longname, summed_CPU_Step[i], summed_CPU_Step[i] / summed_CPU_Step[CPU_ALL] * 100., CPU_Sum[i],
                  CPU_Sum[i] / CPU_Sum[CPU_ALL] * 100.);
          fprintf(FdCPUCSV, "%f, %f, %f, ", summed_CPU_Step[i], CPU_Sum[i], CPU_Sum[i] / CPU_Sum[CPU_ALL] * 100.);
        }

      fprintf(FdCPU, "\n");
      myflush(FdCPU);

      fprintf(FdCPUCSV, "\n");
      myflush(FdCPUCSV);
    }

  CPUThisRun = Logs.timediff(StartOfRun, Logs.second());

  for(int i = 0; i < CPU_LAST; i++)
    CPU_Step[i] = 0.;

  TIMER_STOP(CPU_LOGS);
}

/*! \brief Fill the cpu balance string representing the cpu usage in a graphical way
 *
 * This function fills a fraction, specified by the parameters t0 and t1,
 * of the array string with the debug symbol given by c.
 *
 * \param string string to fill
 * \param t0 initial position of the symbol in the array as a fraction of its maximum dimension
 * \param t1 final position of the symbol in the array as a fraction of its maximum dimension
 * \param c symbol to be put on string
 */
void logs::put_symbol(char *string, double t0, double t1, char c)
{
  int i = (int)(t0 * CPU_STRING_LEN + 0.5);
  int j = (int)(t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    string[i++] = c;

  string[CPU_STRING_LEN] = 0;
}

double logs::measure_time(void) /* strategy: call this at end of functions to account for time in this function, and before another
                                   (nontrivial) function is called */
{
  double t      = Logs.second();
  double dt     = t - WallclockTime;
  WallclockTime = t;

  return dt;
}

/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double logs::second(void)
{
  return MPI_Wtime();

  /*
   * possible alternative:
   *
   * return ((double) clock()) / CLOCKS_PER_SEC;
   *
   * but note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}

/* returns the time difference between two measurements
 * obtained with Logs.second().
 */
double logs::timediff(double t0, double t1)
{
  double dt = t1 - t0;

  if(dt < 0)
    dt = 0;

  return dt;
}

/*! \brief Computes new global statistics if needed (done by energy_statistics())
 *
 */
void logs::compute_statistics(void)
{
  /* check whether we want a full energy statistics */
  if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics &&
     All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) /* allow only top-level synchronization points */
    {
      compute_global_quantities_of_system();

      if(ThisTask == 0)
        {
          fprintf(FdEnergy, "%14.8g %14.8g %14.8g %14.8g", All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin);

          for(int i = 0; i < NTYPES; i++)
            fprintf(FdEnergy, " %14.8g %14.8g %14.8g", SysState.EnergyIntComp[i], SysState.EnergyPotComp[i],
                    SysState.EnergyKinComp[i]);

          for(int i = 0; i < NTYPES; i++)
            fprintf(FdEnergy, " %14.8g", SysState.MassComp[i]);

          fprintf(FdEnergy, "\n");

          myflush(FdEnergy);
        }

      All.TimeLastStatistics += All.TimeBetStatistics;
    }
}

#ifdef MEASURE_TOTAL_MOMENTUM
void logs::compute_total_momentum(void)
{
  double mom[3] = {0, 0, 0};

  for(int i = 0; i < Sp->NumPart; i++)
    {
      for(int j = 0; j < 3; j++)
        mom[j] += Sp->P[i].getMass() * Sp->P[i].Vel[j];
    }

  // some the stuff over all processors
  double momsum[3] = {0, 0, 0};
  MPI_Reduce(mom, momsum, 3, MPI_DOUBLE, MPI_SUM, 0, Communicator);

  if(ThisTask == 0)
    {
      fprintf(FdMomentum, "%14.8g   %25.15g %25.15g %25.15g\n", All.Time, momsum[0], momsum[1], momsum[2]);

      myflush(FdMomentum);
    }
}
#endif

/*! \brief This routine computes various global properties of the particle
 * distribution and stores the result in the struct `SysState'.
 *
 * Currently, not all the information that's computed here is
 * actually used (e.g. momentum is not really used anywhere),
 * just the energies are written to a log-file every once in a while.
 */
void logs::compute_global_quantities_of_system(void)
{
  state_of_system sys;
  double egyspec;

  particle_data *P = Sp->P;
  // sph_particle_data *SphP = Sp->SphP;

  All.set_cosmo_factors_for_current_time();

  for(int n = 0; n < NTYPES; n++)
    {
      sys.MassComp[n] = sys.EnergyKinComp[n] = sys.EnergyPotComp[n] = sys.EnergyIntComp[n] = 0;

      for(int j = 0; j < 4; j++)
        sys.CenterOfMassComp[n][j] = sys.MomentumComp[n][j] = sys.AngMomentumComp[n][j] = 0;
    }

  for(int i = 0; i < Sp->NumPart; i++)
    {
      sys.MassComp[P[i].getType()] += Sp->P[i].getMass();

#if defined(SELFGRAVITY) && defined(EVALPOTENTIAL)
      sys.EnergyPotComp[P[i].getType()] += 0.5 * Sp->P[i].getMass() * P[i].Potential / All.cf_atime;
#endif

#if defined(EXTERNALGRAVITY) && defined(EVALPOTENTIAL)
      sys.EnergyPotComp[P[i].getType()] += P[i].getMass() * P[i].ExtPotential;
#endif

      double vel[3] = {0, 0, 0};

      if(P[i].getType() == 0)
        {
          for(int j = 0; j < 3; j++)
            vel[j] = P[i].Vel[j];

          sys.EnergyKinComp[0] += 0.5 * Sp->P[i].getMass() * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

          egyspec = Sp->get_utherm_from_entropy(i);

          sys.EnergyIntComp[0] += Sp->P[i].getMass() * egyspec;
        }
#if(NTYPES > 1)
      else
        {
          for(int j = 0; j < 3; j++)
            vel[j] = P[i].Vel[j];

          sys.EnergyKinComp[P[i].getType()] +=
              0.5 * Sp->P[i].getMass() * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) * All.cf_a2inv;
        }
#endif

      double pos[3];
      Sp->intpos_to_pos(P[i].IntPos, pos);  // converts the integer coordinates to floating point

      for(int j = 0; j < 3; j++)
        {
          sys.MomentumComp[P[i].getType()][j] += Sp->P[i].getMass() * vel[j];
          sys.CenterOfMassComp[P[i].getType()][j] += Sp->P[i].getMass() * pos[j];
        }

      sys.AngMomentumComp[P[i].getType()][0] += Sp->P[i].getMass() * (pos[1] * vel[2] - pos[2] * vel[1]);
      sys.AngMomentumComp[P[i].getType()][1] += Sp->P[i].getMass() * (pos[2] * vel[0] - pos[0] * vel[2]);
      sys.AngMomentumComp[P[i].getType()][2] += Sp->P[i].getMass() * (pos[0] * vel[1] - pos[1] * vel[0]);
    }

  // some the stuff over all processors
  MPI_Reduce(&sys.MassComp[0], &SysState.MassComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, Communicator);
  MPI_Reduce(&sys.EnergyPotComp[0], &SysState.EnergyPotComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, Communicator);
  MPI_Reduce(&sys.EnergyIntComp[0], &SysState.EnergyIntComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, Communicator);
  MPI_Reduce(&sys.EnergyKinComp[0], &SysState.EnergyKinComp[0], NTYPES, MPI_DOUBLE, MPI_SUM, 0, Communicator);
  MPI_Reduce(&sys.MomentumComp[0][0], &SysState.MomentumComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, Communicator);
  MPI_Reduce(&sys.AngMomentumComp[0][0], &SysState.AngMomentumComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, Communicator);
  MPI_Reduce(&sys.CenterOfMassComp[0][0], &SysState.CenterOfMassComp[0][0], NTYPES * 4, MPI_DOUBLE, MPI_SUM, 0, Communicator);

  if(ThisTask == 0)
    {
      for(int i = 0; i < NTYPES; i++)
        SysState.EnergyTotComp[i] = SysState.EnergyKinComp[i] + SysState.EnergyPotComp[i] + SysState.EnergyIntComp[i];

      SysState.Mass = SysState.EnergyKin = SysState.EnergyPot = SysState.EnergyInt = SysState.EnergyTot = 0;

      for(int j = 0; j < 3; j++)
        SysState.Momentum[j] = SysState.AngMomentum[j] = SysState.CenterOfMass[j] = 0;

      for(int i = 0; i < NTYPES; i++)
        {
          SysState.Mass += SysState.MassComp[i];
          SysState.EnergyKin += SysState.EnergyKinComp[i];
          SysState.EnergyPot += SysState.EnergyPotComp[i];
          SysState.EnergyInt += SysState.EnergyIntComp[i];
          SysState.EnergyTot += SysState.EnergyTotComp[i];

          for(int j = 0; j < 3; j++)
            {
              SysState.Momentum[j] += SysState.MomentumComp[i][j];
              SysState.AngMomentum[j] += SysState.AngMomentumComp[i][j];
              SysState.CenterOfMass[j] += SysState.CenterOfMassComp[i][j];
            }
        }

      for(int i = 0; i < NTYPES; i++)
        for(int j = 0; j < 3; j++)
          if(SysState.MassComp[i] > 0)
            SysState.CenterOfMassComp[i][j] /= SysState.MassComp[i];

      for(int j = 0; j < 3; j++)
        if(SysState.Mass > 0)
          SysState.CenterOfMass[j] /= SysState.Mass;

      for(int i = 0; i < NTYPES; i++)
        {
          SysState.CenterOfMassComp[i][3] = SysState.MomentumComp[i][3] = SysState.AngMomentumComp[i][3] = 0;
          for(int j = 0; j < 3; j++)
            {
              SysState.CenterOfMassComp[i][3] += SysState.CenterOfMassComp[i][j] * SysState.CenterOfMassComp[i][j];
              SysState.MomentumComp[i][3] += SysState.MomentumComp[i][j] * SysState.MomentumComp[i][j];
              SysState.AngMomentumComp[i][3] += SysState.AngMomentumComp[i][j] * SysState.AngMomentumComp[i][j];
            }
          SysState.CenterOfMassComp[i][3] = sqrt(SysState.CenterOfMassComp[i][3]);
          SysState.MomentumComp[i][3]     = sqrt(SysState.MomentumComp[i][3]);
          SysState.AngMomentumComp[i][3]  = sqrt(SysState.AngMomentumComp[i][3]);
        }

      SysState.CenterOfMass[3] = SysState.Momentum[3] = SysState.AngMomentum[3] = 0;

      for(int j = 0; j < 3; j++)
        {
          SysState.CenterOfMass[3] += SysState.CenterOfMass[j] * SysState.CenterOfMass[j];
          SysState.Momentum[3] += SysState.Momentum[j] * SysState.Momentum[j];
          SysState.AngMomentum[3] += SysState.AngMomentum[j] * SysState.AngMomentum[j];
        }

      SysState.CenterOfMass[3] = sqrt(SysState.CenterOfMass[3]);
      SysState.Momentum[3]     = sqrt(SysState.Momentum[3]);
      SysState.AngMomentum[3]  = sqrt(SysState.AngMomentum[3]);
    }

  // give everyone the result, maybe the want to do something with it
  MPI_Bcast(&SysState, sizeof(state_of_system), MPI_BYTE, 0, Communicator);
}
