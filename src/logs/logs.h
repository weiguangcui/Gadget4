/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file logs.h
 *
 *  \brief declares constants and some classes used for log-file handling
 */

#ifndef LOGS_H
#define LOGS_H

#include "gadgetconfig.h"

#include <stdio.h>

#include "../main/simulation.h"
#include "../mpi_utils/setcomm.h"

#define CPU_STRING_LEN 120

#define TIMER_STACK_DEPTH 30

class simparticles;
class lcparticles;

class logs : public setcomm
{
 public:
  logs() : setcomm("delayed init")  // constructor
  {

#define TIMER_STRUCT
#include "../logs/timer.h"
#undef TIMER_STRUCT

  }



  simparticles *Sp;

  double CPUThisRun; /*!< Sums CPU time of current process */

  FILE *FdInfo,   /**< file handle for info.txt log-file. */
      *FdEnergy,  /**< file handle for energy.txt log-file. */
      *FdTimings, /**< file handle for timings.txt log-file. */
      *FdDensity, /**< file handle for timings.txt log-file. */
      *FdHydro,   /**< file handle for timings.txt log-file. */
      *FdBalance, /**< file handle for balance.txt log-file. */
      *FdTimebin, /**< file handle for timebins.txt log-file. */
      *FdDomain,  /**< file handle for domain.txt log-file. */
      *FdCPU;     /**< file handle for cpu.txt log-file. */

  FILE *FdCPUCSV; /**< file handle for cpu.csv log-file. Used if the cpu log is printed in csv format as well. */

#ifdef MEASURE_TOTAL_MOMENTUM
  FILE *FdMomentum;
#endif

#ifdef STARFORMATION
  FILE *FdSfr; /**< file handle for sfr.txt log-file. */
#endif

#ifdef DEBUG_MD5
  FILE *FdDebug; /**< file handle for debug_md5.txt log-file. */
#endif

#ifdef FORCETEST
  FILE *FdForceTest; /*!< file handle for forcetest.txt log-file. */
#endif

  void init_cpu_log(simparticles *Sp_ptr);
  void open_logfiles(void);
  void write_cpu_log(void);
  void output_log_messages(void);
  void compute_statistics(void);
  void print_particle_info_from_ID(MyIDType ID);
  void print_particle_info(int i);
  void log_debug_md5(const char *msg);
  void calc_memory_checksum(const char *msg, void *base, size_t bytes);
  void block_checksum(void *base, size_t bytes, int res[4]);
  int flush_everything(void);
  void compute_total_momentum(void);

  double measure_time(void);
  double timediff(double t0, double t1);
  double second(void);

  enum timers
  {
    CPU_NONE = -2, /*!< used for counters without a parent */
    CPU_ROOT = -1, /*!< root node of the tree */

#define TIMER_ENUM
#include "../logs/timer.h"
#undef TIMER_ENUM

  };

  /*! \brief struct containing the information of a CPU timer
   *
   */
  struct timer_d
  {
    int parent;         /*!< id of the parent timer */
    char shortname[30]; /*!< string containing the internal name of the timer */
    char longname[30];  /*!< name of the timer */
    char symb;          /*!< symbol used in balance.txt for the active part */
    char symbImbal;     /*!< symbol used in balance.txt for imbalance */
    char depth;         /*!< depth in the tree-like structure of this timer */
  };
  timer_d Timer_data[CPU_LAST + 1];

  double CPU_Step[CPU_LAST];
  double CPU_Step_Stored[CPU_LAST];
  double CPU_Sum[CPU_LAST]; /**< sums wallclock time/CPU consumption in whole run */

  enum timers TimerStack[TIMER_STACK_DEPTH];
  int TimerStackPos = 0;

 private:
  double StartOfRun; /*!< This stores the time of the start of the run for evaluating the elapsed time */

  double WallclockTime; /*!< This holds the last wallclock time measurement for timings measurements */

  void put_symbol(char *string, double t0, double t1, char c);

  /* global state of system */
  struct state_of_system
  {
    double Mass, EnergyKin, EnergyPot, EnergyInt, EnergyTot, Momentum[4], AngMomentum[4], CenterOfMass[4], MassComp[NTYPES],
        EnergyKinComp[NTYPES], EnergyPotComp[NTYPES], EnergyIntComp[NTYPES], EnergyTotComp[NTYPES], MomentumComp[NTYPES][4],
        AngMomentumComp[NTYPES][4], CenterOfMassComp[NTYPES][4];
  };
  state_of_system SysState;

  void compute_global_quantities_of_system(void);

 public:
  void timer_stop(enum timers counter)
  {
    if(TimerStack[TimerStackPos] != counter)
      Terminate("Wrong use of timer_stop(), you must stop the timer started last");

    CPU_Step[TimerStack[TimerStackPos--]] += measure_time();

    if(TimerStackPos < 0)
      Terminate("Do not stop the out CPU_MISC timer");
  }

  void timer_start(enum timers counter)
  {
    CPU_Step[TimerStack[TimerStackPos]] += measure_time();

    for(int itimer = 0; itimer <= TimerStackPos; itimer++)
      if(counter == TimerStack[itimer])
        Terminate("Try to start timer %d, but it is already running.\n", counter);

    if(++TimerStackPos >= TIMER_STACK_DEPTH)
      {
        Terminate("Run out of timer stack space, increase TIMER_STACK_DEPTH");
      }
    else
      TimerStack[TimerStackPos] = counter;
  }
};

#include "../data/simparticles.h"

#define TIMER_START_INTERNAL(counter)                                                      \
  {                                                                                        \
    Logs.CPU_Step[Logs.TimerStack[Logs.TimerStackPos]] += Logs.measure_time();             \
    for(int itimer = 0; itimer <= Logs.TimerStackPos; itimer++)                            \
      if(logs::counter == Logs.TimerStack[itimer])                                         \
        {                                                                                  \
          Terminate("Try to start timer %d, but it is already running.\n", logs::counter); \
        };                                                                                 \
    if(++Logs.TimerStackPos >= TIMER_STACK_DEPTH)                                          \
      {                                                                                    \
        Terminate("Run out of timer stack space, increase TIMER_STACK_DEPTH");             \
      }                                                                                    \
    else                                                                                   \
      {                                                                                    \
        Logs.TimerStack[Logs.TimerStackPos] = logs::counter;                               \
      }                                                                                    \
  }

/*! \def  TIMER_START(counter)
 * \brief Starts the timer counter
 *
 * Use this macro instead of directly accessing the CPU_Step array,
 * so manual instrumentation APIs can be attached.
 *
 * \param counter Name of the timer to start
 */
#define TIMER_START(counter) TIMER_START_INTERNAL(counter)

#define TIMER_STOP_INTERNAL(counter)                                                \
  {                                                                                 \
    if(Logs.TimerStack[Logs.TimerStackPos] != logs::counter)                        \
      {                                                                             \
        Terminate("Wrong use of TIMER_STOP, you must stop the timer started last"); \
      }                                                                             \
    Logs.CPU_Step[Logs.TimerStack[Logs.TimerStackPos--]] += Logs.measure_time();    \
    if(Logs.TimerStackPos < 0)                                                      \
      {                                                                             \
        Terminate("Do not stop the out CPU_MISC timer");                            \
      }                                                                             \
  }

/*! \def TIMER_STOP(counter)
 * \brief Stops the timer counter
 *
 * Use this macro instead of directly accessing the CPU_Step array,
 * so manual instrumentation APIs can be attached.
 *
 * \param counter Name of the timer to stop
 */
#define TIMER_STOP(counter) TIMER_STOP_INTERNAL(counter)

/*! \def TIMER_STOPSTART(stop, start)
 * \brief Stops the timer 'stop' and starts the timer 'start'
 *
 * Use this macro instead of directly accessing the CPU_Step array,
 * so manual instrumentation APIs can be attached.
 *
 * \param stop Name of the timer to stop
 * \param start Name of the timer to start
 */
#define TIMER_STOPSTART(stop, start) \
  {                                  \
    TIMER_STOP_INTERNAL(stop);       \
    TIMER_START_INTERNAL(start);     \
  }

/*! \def TIMER_ADD(counter, amount)
 * \brief Adds amount to the timer counter

 * \param counter Name of the timer to add to
 * \param amount amount to add to timer counter
 */
#define TIMER_ADD(counter, amount) Logs.CPU_Step[counter] += (amount);

/*! \def TIMER_DIFF(counter)
 * \brief Returns amount elapsed for the timer since last save with TIMER_STORE

 * \param counter Name of the timer to add to
 */
#define TIMER_DIFF(counter) (Logs.CPU_Step[logs::counter] - Logs.CPU_Step_Stored[logs::counter])

/*! \def TIMER_STORE
 * \brief Copies the current value of CPU times to a stored variable, such that differences with respect to this reference can be
 * calculated
 *
 */
#define TIMER_STORE memcpy(Logs.CPU_Step_Stored, Logs.CPU_Step, sizeof(Logs.CPU_Step));

extern logs Logs;

#endif
