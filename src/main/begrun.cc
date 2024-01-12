/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file begrun.cc
 *
 *  \brief initial set-up of a simulation run
 */

#include "gadgetconfig.h"
#include "compiler-command-line-args.h"

#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../data/simparticles.h"
#include "../fmm/fmm.h"
#include "../fof/fof.h"
#include "../gitversion/version.h"
#include "../gravity/ewald.h"
#include "../gravtree/gravtree.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../io/parameters.h"
#include "../lightcone/lightcone.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../ngbtree/ngbtree.h"
#include "../pm/pm.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"
#include "../time_integration/timestep.h"

/*!
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameter file is read in and parsed and global variables
 *  are initialized to their proper values.
 */

void sim::hello(void)
{
  mpi_printf(
      "\n  ___    __    ____    ___  ____  ____       __\n / __)  /__\\  (  _ \\  / __)( ___)(_  _)___  /. |\n"
      "( (_-. /(__)\\  )(_) )( (_-. )__)   )( (___)(_  _)\n \\___/(__)(__)(____/  \\___/(____) (__)       (_)\n\n");

  mpi_printf("This is Gadget, version %s.\nGit commit %s, %s\n\n", GADGET_VERSION, GIT_COMMIT, GIT_DATE);

  mpi_printf("Code was compiled with the following compiler and flags:\n%s\n\n\n", compiler_flags);

  mpi_printf("Code was compiled with the following settings:\n");

  if(ThisTask == 0)
    {
      output_compile_time_options();
    }

  mpi_printf("\n\nRunning on %d MPI tasks.\n\n", NTask);

#ifdef EXPLICIT_VECTORIZATION

  int instrset = instrset_detect();
  mpi_printf("CPU supports instruction sets: ");
  if(instrset >= 1)
    mpi_printf("SSE ");
  if(instrset >= 2)
    mpi_printf("SSE2 ");
  if(instrset >= 3)
    mpi_printf("SSE3 ");
  if(instrset >= 4)
    mpi_printf("SSSE3 ");
  if(instrset >= 5)
    mpi_printf("SSE4.1 ");
  if(instrset >= 6)
    mpi_printf("SSE4.2 ");
  if(instrset >= 7)
    mpi_printf("AVX ");
  if(instrset >= 8)
    mpi_printf("AVX2 ");
  if(instrset >= 9)
    mpi_printf("AVX512F ");
  if(instrset >= 10)
    mpi_printf("AVX512VL AVX512BW AVX512DQ ");
  mpi_printf("\n\n");

  if(instrset < 7)
    {
      mpi_printf(
          "You compiled with explicit vectorization in terms of AVX instructions, but this CPU does not support AVX or higher.\n\n");
      endrun();
    }
#endif

  mpi_printf("\n");
  mpi_printf("BEGRUN: Size of particle structure       %4d  [bytes]\n", (int)sizeof(particle_data));
  mpi_printf("BEGRUN: Size of sph particle structure   %4d  [bytes]\n", (int)sizeof(sph_particle_data));
  mpi_printf("BEGRUN: Size of gravity tree node        %4d  [bytes]\n", (int)sizeof(gravnode));
  mpi_printf("BEGRUN: Size of neighbour tree node      %4d  [bytes]\n", (int)sizeof(ngbnode));
  mpi_printf("BEGRUN: Size of subfind auxiliary data   %4d  [bytes]\n", (int)sizeof(subfind_data));

  mpi_printf("\n");
}

/*! \brief This function performs the initial set-up of the simulation.
 *
 *  First, the parameter file is read by read_parameter_file(),
 *  then routines for setting units, etc are called. This function only does
 *  the setup necessary to load the IC file. After the IC file has been loaded
 *  and prepared by init(), setup continues with begrun2(). This splitting is
 *  done so that we can return cleanly from operations that don't actually
 *  start the simulation (converting snapshots, making projected images, etc.)
 */
void sim::begrun1(const char *parameterFile)
{
  All.register_parameters();

  int errorFlag = All.read_parameter_file(parameterFile); /* ... read in parameters for this run on task 0*/

  if(ThisTask == 0)
    {
      int n = strlen(All.OutputDir);
      if(n > 0)
        if(All.OutputDir[n - 1] != '/')
          strcat(All.OutputDir, "/");

      mkdir(All.OutputDir, 02755);
    }

  /* now communicate the relevant parameters to the other processes, *including* the shared memory handler */
  /* this also tells the shared memory handler how much memory it may allocate */
  MPI_Bcast(All.get_data_ptr(), All.get_data_size(), MPI_BYTE, 0, MPI_COMM_WORLD);

#ifdef HOST_MEMORY_REPORTING
  Mem.check_maxmemsize_setting(All.MaxMemSize);
#endif

  Mem.mymalloc_init(All.MaxMemSize, All.RestartFlag);

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, Communicator);

  if(errorFlag)
    {
      if(Shmem.Island_ThisTask == 0 && Shmem.Island_NTask != Shmem.World_NTask)
        {
          char c = 0;
          // need to send this flag to our shared memory rank so that it also ends itself
          MPI_Send(&c, 1, MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_KEY, MPI_COMM_WORLD);
        }

      mpi_printf("\nWe stop because of an error in the parameterfile.\n\n");
      MPI_Finalize();
      exit(0);
    }

  if(All.OutputListOn)
    All.read_outputlist(All.OutputListFilename);
  else
    All.OutputListLength = 0;

  All.write_used_parameters(All.OutputDir, "parameters-usedvalues");

  All.some_parameter_checks();

#ifdef ENABLE_HEALTHTEST
  healthtest();
#endif

  set_units();

  my_create_HDF5_halfprec_handler();
  my_create_HDF5_special_integer_types();

  my_mpi_types_init();

#ifdef LIGHTCONE_PARTICLES
  LightCone.lightcone_init_geometry(All.LightConeDefinitionFile);
#endif

  /* disable automatic error printing of HDF5 - we check for errors ourselves */
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

#ifdef DEBUG
  enable_core_dumps_and_fpu_exceptions();
#endif

#ifdef PMGRID
  GravTree.short_range_init();
#endif

  Sp.TimeBinsHydro.timebins_init("Hydro", &Sp.MaxPartSph);
  Sp.TimeBinsGravity.timebins_init("Gravity", &Sp.MaxPart);

#ifdef COOLING
  All.Time = All.TimeBegin;
  All.set_cosmo_factors_for_current_time();
  CoolSfr.InitCool();
#endif

#ifdef STARFORMATION
  CoolSfr.init_clouds();
#endif

#if((!defined(PMGRID) || (defined(PMGRID) && defined(TREEPM_NOTIMESPLIT))) && defined(SELFGRAVITY) && defined(PERIODIC)) || \
    defined(FORCETEST)
  Ewald.ewald_init();
#endif

  init_rng(ThisTask);

#ifdef DEBUG_SYMTENSORS
  symtensor_test();
#endif

#ifdef PMGRID
  if(All.RestartFlag == RST_BEGIN || All.RestartFlag == RST_RESUME || All.RestartFlag == RST_STARTFROMSNAP ||
     All.RestartFlag == RST_POWERSPEC)
    {
#ifdef PERIODIC
      PM.pm_init_periodic(&Sp);
#ifdef PLACEHIGHRESREGION
      PM.pm_init_nonperiodic(&Sp);
#endif
#else
      PM.pm_init_nonperiodic(&Sp);
#endif
    }
#endif

  Logs.open_logfiles();

  All.TimeLastRestartFile = Logs.CPUThisRun;

#ifdef DEBUG_MD5
  Logs.log_debug_md5("START");
#endif
}

/*! \brief This function does late setup, after the IC file has been loaded
 *  but before run() is called.
 *
 *  The output files are opened and various modules are initialized. The next output
 *  time is determined by find_next_outputtime() and various timers are set.
 *
 */
void sim::begrun2(void)
{
  char contfname[MAXLEN_PATH_EXTRA];
  snprintf(contfname, MAXLEN_PATH_EXTRA, "%scont", All.OutputDir);
  unlink(contfname);

  if(All.RestartFlag != RST_BEGIN && All.RestartFlag != RST_RESUME && All.RestartFlag != RST_STARTFROMSNAP)
    Logs.open_logfiles();

  if(All.ComovingIntegrationOn)
    {
      Driftfac.init_drift_table();

      if(Shmem.Island_NTask != Shmem.World_NTask && Shmem.Island_ThisTask == 0)
        {
          // We actually have multiple shared memory nodes in which we set aside one MPI rank for shared memory communication.
          // Tell the shared memory node to initialize its drift table as well, so that it can be used for particle preditions on the
          // fly
          if(Shmem.Island_ThisTask == 0)
            {
              // update All on shared memory handler, to be sure that be get the correct times
              MPI_Send(All.get_data_ptr(), All.get_data_size(), MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_DRIFT_INIT, MPI_COMM_WORLD);
            }
        }
    }

#ifdef LIGHTCONE
#ifdef LIGHTCONE_PARTICLES
  if(LightCone.lightcone_init_times())
    endrun();
#endif
#ifdef LIGHTCONE_MASSMAPS
  LightCone.lightcone_init_massmaps();
  if(LightCone.lightcone_massmap_report_boundaries())
    endrun();
#endif
  if(LightCone.lightcone_init_boxlist())
    endrun();

  double linklength = 0;

#ifdef FOF
  fof<simparticles> FoF{Communicator, &Sp, &Domain};
  linklength = FoF.fof_get_comoving_linking_length();
#endif

  LightCone.lightcone_init_intposconverter(linklength);
#endif

  if(All.RestartFlag == RST_STARTFROMSNAP)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
  else if(All.RestartFlag == RST_RESUME)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);

  All.TimeLastRestartFile = Logs.CPUThisRun;

#ifdef REDUCE_FLUSH
  All.FlushLast = Logs.CPUThisRun;
#endif

  // update All on shared memory handler, just to allow it to access its elements if needed
  if(Shmem.Island_NTask != Shmem.World_NTask && Shmem.Island_ThisTask == 0)
    MPI_Send(All.get_data_ptr(), All.get_data_size(), MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_ALL_UPDATE, MPI_COMM_WORLD);

#if defined(FORCETEST) && defined(FORCETEST_TESTFORCELAW)
  gravity_forcetest_testforcelaw();
#endif
}

/*! \brief Computes conversion factors between internal code units and the
 *  cgs-system.
 *
 *  In addition constants like the gravitation constant are set.
 */
void sim::set_units(void)
{
  All.UnitTime_in_s         = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;
  All.UnitTime_in_years     = All.UnitTime_in_s / SEC_PER_YEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs     = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs    = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs      = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  if(All.ComovingIntegrationOn)
    {
      All.OmegaCurvature = 1.0 - (All.Omega0 + All.OmegaLambda);

      /* check whether the supplied value of All.Hubble makes sense */
      if(All.HubbleParam != 1.0)
        {
          double Hubble_expected = HUBBLE * All.UnitTime_in_s;
          if(fabs(Hubble_expected - All.Hubble) > 1.0e-3 * All.Hubble)
            Terminate(
                "You have supplied All.Hubble=%g/All.HubbleParam=%g. For this choice, we would expect All.Hubble=%g. We better stop\n",
                All.Hubble, All.HubbleParam, Hubble_expected);
        }
      else
        {
          double Hubble_expected = 0.7 * HUBBLE * All.UnitTime_in_s;
          if(All.Hubble < 0.5 * Hubble_expected || All.Hubble > 2.0 * Hubble_expected)
            Terminate(
                "You have supplied All.Hubble=%g/All.HubbleParam=%g. For this choice, we would expect All.Hubble somewhere in the "
                "ballpark of %g. We better stop\n",
                All.Hubble, All.HubbleParam, Hubble_expected);
        }
    }

  mpi_printf("BEGRUN: Hubble (internal units)  = %g\n", All.Hubble);
  mpi_printf("BEGRUN: h                        = %g\n", All.HubbleParam);
  mpi_printf("BEGRUN: G (internal units)       = %g\n", All.G);
  mpi_printf("BEGRUN: UnitMass_in_g            = %g\n", All.UnitMass_in_g);
  mpi_printf("BEGRUN: UnitLenth_in_cm          = %g\n", All.UnitLength_in_cm);
  mpi_printf("BEGRUN: UnitTime_in_s            = %g\n", All.UnitTime_in_s);
  mpi_printf("BEGRUN: UnitVelocity_in_cm_per_s = %g\n", All.UnitVelocity_in_cm_per_s);
  mpi_printf("BEGRUN: UnitDensity_in_cgs       = %g\n", All.UnitDensity_in_cgs);
  mpi_printf("BEGRUN: UnitEnergy_in_cgs        = %g\n", All.UnitEnergy_in_cgs);
  mpi_printf("\n");

#ifdef STARFORMATION
  CoolSfr.set_units_sfr();
#endif
}

/** \brief This function aborts the simulations.
 *
 * This method has to be called by all processes. It should be used only
 * if the simulation ends without a errors or a an error message is already printed.
 * Otherwise Terminate() should be used instead.
 */
void sim::endrun(void)
{
  mpi_printf("endrun called, calling MPI_Finalize()\nbye!\n\n");
  fflush(stdout);

  if(Shmem.Island_ThisTask == 0 && Shmem.Island_NTask != Shmem.World_NTask)
    {
      char c = 0;
      // need to send this flag to our shared memory rank so that it also ends itself
      MPI_Send(&c, 1, MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_KEY, MPI_COMM_WORLD);
    }

  /* The hdf5 library will sometimes register an atexit() handler that calls its error handler.
   * This is set to my_hdf_error_handler, which calls MPI_Abort.
   * Calling MPI_Abort after MPI_Finalize is not allowed.
   * Hence unset the HDF error handler here*/
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  MPI_Finalize();
  exit(0);
}
