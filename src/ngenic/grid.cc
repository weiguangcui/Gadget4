/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  grid.cc
 *
 *  \brief routines for setting up an unperturbed particle load
 */

#include "gadgetconfig.h"

#ifdef CREATE_GRID

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngenic/ngenic.h"
#include "../system/system.h"

void ngenic::create_grid(void)
{
  long long gridSize    = All.GridSize;
  long long partTotal   = gridSize * gridSize * gridSize;
  long long partPerTask = partTotal / NTask;

  Sp->RegionLen     = All.BoxSize;
  Sp->FacCoordToInt = pow(2.0, BITS_FOR_POSITIONS) / Sp->RegionLen;
  Sp->FacIntToCoord = Sp->RegionLen / pow(2.0, BITS_FOR_POSITIONS);

  All.Time = All.TimeBegin;

  double masstot = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G) * All.BoxSize * All.BoxSize * All.BoxSize;

  double m = masstot / (partTotal);

  for(int i = 0; i < NTYPES; i++)
    All.MassTable[i] = 0.;

  All.MassTable[1] = m;

  Sp->NumGas  = 0;
  Sp->NumPart = partPerTask;

  if(ThisTask == NTask - 1)
    {
      Sp->NumPart = partTotal - Sp->NumPart * (NTask - 1);
    }

  int max_load, max_sphload;
  MPI_Allreduce(&Sp->NumPart, &max_load, 1, MPI_INT, MPI_MAX, Communicator);
  MPI_Allreduce(&Sp->NumGas, &max_sphload, 1, MPI_INT, MPI_MAX, Communicator);

#ifdef GENERATE_GAS_IN_ICS
  Sp->TotNumGas  = partTotal;
  Sp->TotNumPart = 2 * partTotal;
  max_sphload    = max_load;
  max_load *= 2;
#else
  Sp->TotNumPart = partTotal;
  Sp->TotNumGas  = 0;
#endif

  Sp->MaxPart    = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
  Sp->MaxPartSph = max_sphload / (1.0 - 2 * ALLOC_TOLERANCE);

  Sp->allocate_memory();

  for(int i = 0; i < Sp->NumPart; i++)
    {
      long long ipcell = ThisTask * partPerTask + i;
      int x            = ipcell / (All.GridSize * All.GridSize);
      int xr           = ipcell % (All.GridSize * All.GridSize);
      int y            = xr / All.GridSize;
      int z            = xr % All.GridSize;

      double xyz[3];
      xyz[0] = x * All.BoxSize / All.GridSize;
      xyz[1] = y * All.BoxSize / All.GridSize;
      xyz[2] = z * All.BoxSize / All.GridSize;

      Sp->pos_to_intpos(xyz, Sp->P[i].IntPos);

      Sp->P[i].Vel[0] = 0.;
      Sp->P[i].Vel[1] = 0.;
      Sp->P[i].Vel[2] = 0.;

      Sp->P[i].ID.set(ipcell + 1);

      Sp->P[i].setType(1);
    }

  mpi_printf("NGENIC: generated grid of size %d\n", All.GridSize);
}

#endif
