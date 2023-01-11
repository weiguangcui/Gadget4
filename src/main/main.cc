/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  main.cc
 *
 *  \brief start of the program
 */

#include "gadgetconfig.h"

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../half/half.hpp"
#include "../io/io.h"
#include "../io/restart.h"
#include "../io/snap_io.h"
#include "../logs/logs.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../ngenic/ngenic.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"

/* create instances of global objects */

global_data_all_processes All;
driftfac Driftfac;
ewald Ewald; /* get an instance of the Ewald correction tables */
logs Logs;
memory Mem; /* our instance of the memory object */
shmem Shmem;

/*!
 *  This function initializes the MPI communication packages, and sets
 *  cpu-time counters to 0.  Then begrun1() is called, which sets up
 *  the simulation. Then either IC's or restart files are loaded. In
 *  case of IC's init() is called which prepares the IC's for the run.
 *  A call to begrun2() finishes the initialization. Finally, run() is
 *  started, the main simulation loop, which iterates over the timesteps.
 */
int main(int argc, char **argv)
{
  /* find out how many core we have per CPU and which ones we can use */
  pinning Pin;
  Pin.detect_topology();
  Pin.get_core_set();

  /* initialize MPI, this may already impose some pinning */
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &Shmem.World_ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &Shmem.World_NTask);

#if NUMBER_OF_MPI_LISTENERS_PER_NODE > 1
  MPI_Comm fullsharedmemnode;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &fullsharedmemnode);

  int fullsharedmemnode_ThisTask, fullsharedmemnode_NTask;
  MPI_Comm_rank(fullsharedmemnode, &fullsharedmemnode_ThisTask);
  MPI_Comm_size(fullsharedmemnode, &fullsharedmemnode_NTask);

  int bin;
  subdivide_evenly_get_bin(fullsharedmemnode_NTask, NUMBER_OF_MPI_LISTENERS_PER_NODE, fullsharedmemnode_ThisTask, &bin);

  MPI_Comm_split(fullsharedmemnode, bin, fullsharedmemnode_ThisTask, &Shmem.SharedMemComm);
#else
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &Shmem.SharedMemComm);
#endif

  MPI_Comm_rank(Shmem.SharedMemComm, &Shmem.Island_ThisTask);
  MPI_Comm_size(Shmem.SharedMemComm, &Shmem.Island_NTask);

  int min_ntask, max_ntask;
  MPI_Allreduce(&Shmem.Island_NTask, &max_ntask, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&Shmem.Island_NTask, &min_ntask, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&Shmem.World_ThisTask, &Shmem.Island_Smallest_WorldTask, 1, MPI_INT, MPI_MIN, Shmem.SharedMemComm);

  if(Shmem.World_ThisTask == 0)
    printf("Shared memory islands host a minimum of %d and a maximum of %d MPI ranks.\n", min_ntask, max_ntask);

  Shmem.GhostRank = 0;

  if(max_ntask < Shmem.World_NTask)
    {
      if(Shmem.Island_ThisTask == Shmem.Island_NTask - 1)  // selected the ghost MPI ranks
        Shmem.GhostRank = 1;

      if(min_ntask > 1)
        {
          int comm_ranks;
          MPI_Allreduce(&Shmem.GhostRank, &comm_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

          if(Shmem.World_ThisTask == 0)
            printf("We shall use %d MPI ranks in total for assisting one-sided communication (%d per shared memory node).\n",
                   comm_ranks, NUMBER_OF_MPI_LISTENERS_PER_NODE);
        }
      else
        {
          if(Shmem.World_ThisTask == 0)
            Terminate("We have shared memory islands with just one MPI rank -- can't put aside this one just for communication");
        }

      if(max_ntask > MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY)
        {
          if(Shmem.World_ThisTask == 0)
            Terminate(
                "We have shared memory islands with %d MPI ranks, which is larger than  MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY=%d\n"
                "You may consider increasing NUMBER_OF_MPI_LISTENERS_PER_NODE, current value is %d\n",
                max_ntask, MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY, NUMBER_OF_MPI_LISTENERS_PER_NODE);
        }
    }

  /* we can now split the communicator into the processing ones, and the ones reserved for communication */
  MPI_Comm_split(MPI_COMM_WORLD, Shmem.GhostRank, Shmem.World_ThisTask, &Shmem.SimulationComm);

  MPI_Comm_rank(Shmem.SimulationComm, &Shmem.Sim_ThisTask);
  MPI_Comm_size(Shmem.SimulationComm, &Shmem.Sim_NTask);

  /* Let's now find out for everyone the global rank of the responsible shared memory ghost */
  Shmem.MyShmRankInGlobal = Shmem.World_ThisTask;

  MPI_Bcast(&Shmem.MyShmRankInGlobal, 1, MPI_INT, Shmem.Island_NTask - 1, Shmem.SharedMemComm);

  if(Shmem.GhostRank == 1)
    {
      Mem.initcomm(Shmem.SimulationComm);
      Shmem.shared_memory_handler();  // note: this call will not return
    }

  /* creating main simulation object, holding all the data and objects of our simulation */
  sim Sim{Shmem.SimulationComm};

  Sim.determine_compute_nodes();

  /* initialize the communicator structures of our global objects */
  Ewald.initcomm(Sim.Communicator);  // a global class
  Logs.initcomm(Sim.Communicator);   // a global class
  All.initcomm(Sim.Communicator);    // a global class
  Mem.initcomm(Sim.Communicator);    // a global class
  Mem.determine_compute_nodes();

  /* output a welcome message */
  Sim.hello();

  /* initialize CPU-time/Wallclock-time measurement */
  Logs.init_cpu_log(&Sim.Sp);

  /* pin the MPI ranks to the available core set */
  Pin.pin_to_core_set(&Sim);

#ifdef HOST_MEMORY_REPORTING
  Sim.mpi_report_comittable_memory();
  Mem.MemoryOnNode       = Sim.MemoryOnNode;
  Mem.SharedMemoryOnNode = Sim.SharedMemoryOnNode;
#endif

  if(argc < 2)
    {
      if(Sim.ThisTask == 0)
        {
          printf("\nParameters are missing.\n");
          printf("Start as ./Gadget4 <ParameterFile> [<RestartFlag>] [<SpecialOptions>]\n");
          printf("\n");
          printf("   RestartFlag    Action\n");
          printf("      %2d          Read initial conditions and start simulation\n", RST_BEGIN);
          printf("      %2d          Read restart files and resume simulation\n", RST_RESUME);
          printf("      %2d          Restart from specified snapshot dump and resume simulation\n", RST_STARTFROMSNAP);
          printf("      %2d          Run FOF and optionally SUBFIND\n", RST_FOF);
          printf("      %2d          Calculate a matter power spectrum\n", RST_POWERSPEC);
          printf("      %2d          Convert snapshot file to different format [input=ICFormat  output=SnapFormat\n", RST_CONVERTSNAP);
          printf("      %2d          Create cosmological initial conditions\n", RST_CREATEICS);
          printf("      %2d          Calculate descendants/progenitors [connecting group catalogues SnapNum-1 and SnapNum]\n",
                 RST_CALCDESC);
          printf("      %2d          Arrange halos in merger trees [using group catalogues up to SnapNum]\n", RST_MAKETREES);
          printf("      %2d          Carry out I/O bandwidth test to determine best setting for number of concurrent reads/writes\n",
                 RST_IOBANDWIDTH);
          printf("      %2d          Rearrange particle-lightcone data in merger tree order <conenr>  <firstnum>  <lastnum>\n",
                 RST_LCREARRANGE);
          printf("      %2d          Rearrange most-bound snapshot data in merger tree order <firstnum>  <lastnum>\n",
                 RST_SNPREARRANGE);
          printf("\n");
        }
      Sim.endrun();
    }

  /*  argv[1]  holds the parameterfile */

  if(argc >= 3)
    All.RestartFlag = (enum restart_options)atoi(argv[2]);
  else
    All.RestartFlag = RST_BEGIN;

  int restartSnapNum = -1;
  if(argc >= 4)
    restartSnapNum = atoi(argv[3]);

  /* Do minimal validation of arguments here rather than in random places in the code */
  if((All.RestartFlag == RST_FOF || All.RestartFlag == RST_POWERSPEC || All.RestartFlag == RST_CONVERTSNAP ||
      All.RestartFlag == RST_CALCDESC || All.RestartFlag == RST_MAKETREES) &&
     restartSnapNum < 0)
    {
      Terminate("Need to give the snapshot number for the chosen option.\n");
    }

  /* set-up run based on parameterfile */
  Sim.begrun1(argv[1]);

  /* see if we are loading a restart file or an IC file */
  if(All.RestartFlag == RST_RESUME)
    {
      restart Restart{Sim.Communicator};
      Restart.load(&Sim);
      All.RestartFlag = RST_RESUME;  // prevent that this is overwritten by All.RestartFlag in restart set
    }
  else
    {
      /* We're reading an IC file. Is it a snapshot or really an IC? */
      char fname[MAXLEN_PATH_EXTRA];

      if(All.RestartFlag != RST_BEGIN && All.RestartFlag != RST_RESUME && restartSnapNum >= 0)
        {
          if(All.NumFilesPerSnapshot > 1)
            snprintf(fname, MAXLEN_PATH_EXTRA, "%s/snapdir_%03d/%s_%03d", All.OutputDir, restartSnapNum, All.SnapshotFileBase,
                     restartSnapNum);
          else
            snprintf(fname, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, restartSnapNum);
        }
      else
        {
          strcpy(fname, All.InitCondFile);
        }

      if(All.RestartFlag == RST_STARTFROMSNAP)
        {
          All.ICFormat = All.SnapFormat;
        }

      if(All.RestartFlag == RST_CALCDESC)
        {
#ifdef MERGERTREE
          Sim.MergerTree.descendants_in_postprocessing(restartSnapNum);
          Sim.endrun();
#else
          Terminate("Compile with option MERGERTREE for this option");
#endif
        }

      if(All.RestartFlag == RST_MAKETREES)
        {
#ifdef MERGERTREE
          Sim.MergerTree.halotrees_construct(restartSnapNum);
          Sim.endrun();
#else
          Terminate("Compile with option MERGERTREE for this option");
#endif
        }

      if(All.RestartFlag == RST_IOBANDWIDTH)
        {
          Sim.measure_io_bandwidth();
          Sim.endrun();
        }

      if(All.RestartFlag == RST_LCREARRANGE)
        {
#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES) && defined(REARRANGE_OPTION)
          Sim.rearrange_lightcone(argc, argv);
#else
          Terminate("need to compile with REARRANGE_OPTION for this option to work\n");
#endif
          Sim.endrun();
        }

      if(All.RestartFlag == RST_SNPREARRANGE)
        {
#if defined(MERGERTREE) && defined(REARRANGE_OPTION)
          Sim.rearrange_snapshot(argc, argv);
#else
          Terminate("need to compile with REARRANGE_OPTION for this option to work\n");
#endif
          Sim.endrun();
        }

#ifdef CREATE_GRID
      if(All.RestartFlag == RST_BEGIN || All.RestartFlag == RST_CREATEICS)
        Sim.Ngenic.create_grid();
      else
#endif
        {
          snap_io Snap(&Sim.Sp, Sim.Communicator, All.ICFormat); /* get an I/O object */

          Snap.read_ic(fname);
        }

      /* If we are supposed to just convert the file, write and exit here. */
      if(All.RestartFlag == RST_CONVERTSNAP)
        {
#ifdef COOLING
          Sim.CoolSfr.InitCool();
#endif
          strcat(All.SnapshotFileBase, "_converted");
          Sim.mpi_printf("Start writing file %s\nSnapNum %d\n", All.SnapshotFileBase, restartSnapNum);

          snap_io Snap(&Sim.Sp, Sim.Communicator, All.SnapFormat); /* get an I/O object */
          Snap.write_snapshot(restartSnapNum, NORMAL_SNAPSHOT);

          Sim.endrun();
        }

      Sim.init(restartSnapNum);
    }

  Sim.begrun2();

  Sim.run(); /* main simulation loop */

  Sim.endrun(); /* clean up & finalize MPI */

  return 0;
}
