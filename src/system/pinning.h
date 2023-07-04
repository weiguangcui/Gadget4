/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pinning.h
 *
 *  \brief declares a class for pinning related work
 */

#ifndef PINNING_H
#define PINNING_H

#include "gadgetconfig.h"

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../main/main.h"
#include "../mpi_utils/setcomm.h"
#include "../system/system.h"

/*! \file pinning.c
 *  \brief examines cpu topology and binds processes and threads to cores
 */

#ifdef IMPOSE_PINNING
#include <hwloc.h>
#endif

class pinning
{
 private:
#ifdef IMPOSE_PINNING
  int flag_pinning_error = 0;

  hwloc_cpuset_t cpuset, cpuset_after_MPI_init;
  hwloc_topology_t topology;
  int topodepth;
  int sockets;
  int cores;
  int pus;
  int hyperthreads_per_core;

#endif
 public:
  void get_core_set(void);
  void detect_topology(void);
  void pin_to_core_set(setcomm *sc);
  void report_pinning(setcomm *sc);
};

#endif
