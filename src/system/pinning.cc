/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pinning.cc
 *
 *  \brief routines to report and modify the pinning of processes/threads to CPU cores
 */

#include "gadgetconfig.h"

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../system/pinning.h"

#define MAX_CORES 4096

void pinning::get_core_set(void)
{
#ifdef IMPOSE_PINNING
  cpuset = hwloc_bitmap_alloc();
  hwloc_get_proc_cpubind(topology, getpid(), cpuset, 0);
#endif
}

void pinning::detect_topology(void)
{
#ifdef IMPOSE_PINNING
  /* Allocate and initialize topology object. */
  hwloc_topology_init(&topology);

  /* Perform the topology detection. */
  hwloc_topology_load(topology);

  /* Get some additional topology information
     in case we need the topology depth later. */
  topodepth = hwloc_topology_get_depth(topology);

  int depth = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);

  if(depth == HWLOC_TYPE_DEPTH_UNKNOWN)
    sockets = -1;
  else
    sockets = hwloc_get_nbobjs_by_depth(topology, depth);

  depth = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);

  if(depth == HWLOC_TYPE_DEPTH_UNKNOWN)
    cores = -1;
  else
    cores = hwloc_get_nbobjs_by_depth(topology, depth);

  depth = hwloc_get_type_depth(topology, HWLOC_OBJ_PU);

  if(depth == HWLOC_TYPE_DEPTH_UNKNOWN)
    pus = -1;
  else
    pus = hwloc_get_nbobjs_by_depth(topology, depth);
#endif
}

void pinning::pin_to_core_set(setcomm *sc)
{
#ifdef IMPOSE_PINNING
  sc->mpi_printf("PINNING: We have %d sockets, %d physical cores and %d logical cores on the first MPI-task's node.\n", sockets, cores,
                 pus);
  if(cores <= 0 || sockets <= 0 || pus <= 0)
    {
      sc->mpi_printf("PINNING: The topology cannot be recognized. We refrain from any pinning attempt.\n");
      flag_pinning_error = 1;
      return;
    }

  hyperthreads_per_core = pus / cores;

  if(hyperthreads_per_core < 1)
    Terminate("Need at least one logical thread per physical core\n");

  if(pus > cores)
    sc->mpi_printf("PINNING: Looks like %d hyperthreads per physical core are in principle possible.\n", hyperthreads_per_core);

  cpuset_after_MPI_init = hwloc_bitmap_alloc();
  hwloc_get_proc_cpubind(topology, getpid(), cpuset_after_MPI_init, 0);

  if(!hwloc_bitmap_isequal(cpuset, cpuset_after_MPI_init))
    sc->mpi_printf("PINNING: Apparently, the MPI library set some pinning itself. We'll override this.\n");

  int available_pus = 0;

  for(int id = hwloc_bitmap_first(cpuset); id != -1; id = hwloc_bitmap_next(cpuset, id))
    available_pus++;

  sc->mpi_printf("PINNING: Looks like %d logical cores are available.\n", available_pus);

  if(available_pus == pus)
    {
      sc->mpi_printf("PINNING: Looks like all available logical cores are at our disposal.\n");
    }
  else
    {
      if(available_pus >= 1)
        {
          sc->mpi_printf("PINNING: Looks like already before start of the code, a tight binding was imposed.\n");
#ifdef IMPOSE_PINNING_OVERRIDE_MODE
          for(int id = 0; id < pus; id++)
            hwloc_bitmap_set(cpuset, id);
          available_pus = pus;
          sc->mpi_printf("PINNING: We are overriding this and make all %d available to us.\n", available_pus);
#else
          sc->mpi_printf(
              "PINNING: We refrain from any pinning attempt ourselves. (This can be changed by setting the compile flag "
              "IMPOSE_PINNING_OVERRIDE_MODE.)\n");
          flag_pinning_error = 1;
          return;
#endif
        }
    }

  char buf[MAX_CORES + 1];

  for(int i = 0; i < pus && i < MAX_CORES; i++)
    if(hwloc_bitmap_isset(cpuset, i))
      buf[i] = '1';
    else
      buf[i] = '-';
  buf[pus] = 0;

  sc->mpi_printf("PINNING: Available logical cores on first node:  %s\n", buf);

  int pus_per_task = available_pus / sc->TasksInThisNode;

  sc->mpi_printf("PINNING: %d logical cores are available per MPI Task.\n", pus_per_task);

  if(pus_per_task <= 0)
    Terminate("Need at least one logical core per MPI task for pinning to make sense.\n");

  /* go through all logical cores in sequence of proximity */
  int depth        = hwloc_get_type_depth(topology, HWLOC_OBJ_PU);
  int cores_before = 0;
  int cid;

  for(cid = 0; cores_before < sc->RankInThisNode * pus_per_task && cid < pus; cid++)
    {
      hwloc_obj_t obj = hwloc_get_obj_by_depth(topology, depth, cid);

      hwloc_cpuset_t cpuset_core = hwloc_bitmap_dup(obj->cpuset);
      if(hwloc_bitmap_isincluded(cpuset_core, cpuset))
        {
          cores_before++;
        }
      hwloc_bitmap_free(cpuset_core);
    }

  /* cid should now be the logical index of the first PU for this MPI task */

  hwloc_obj_t obj            = hwloc_get_obj_by_depth(topology, depth, cid);
  hwloc_cpuset_t current_cpu = hwloc_bitmap_dup(obj->cpuset);

  hwloc_set_proc_cpubind(topology, getpid(), current_cpu, HWLOC_CPUBIND_PROCESS);
#endif
}

void pinning::report_pinning(setcomm *sc)
{
#ifdef IMPOSE_PINNING
  if(flag_pinning_error)
    return;

  hwloc_get_cpubind(topology, cpuset, 0);

  char buf[MAX_CORES + 1];

  for(int i = 0; i < pus && i < MAX_CORES; i++)
    if(hwloc_bitmap_isset(cpuset, i))
      buf[i] = '1';
    else
      buf[i] = '-';
  buf[pus] = 0;

  for(int i = 0; i < sc->NTask; i++)
    {
      if(sc->ThisTask == i && sc->ThisNode == 0)
        printf("PINNING: Node=%4d: Task=%04d:                   %s\n", sc->ThisNode, sc->ThisTask, buf);
      fflush(stdout);
      MPI_Barrier(sc->Communicator);
    }
#endif
}
