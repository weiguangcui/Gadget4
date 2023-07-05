/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file mmparticles.h
 *
 *  \brief defines class that holds particles for projection on lightcone massmaps
 */

#ifndef MMPART_H
#define MMPART_H

/* mass map particle */

#if defined(LIGHTCONE) && defined(LIGHTCONE_MASSMAPS)

#include "gadgetconfig.h"

#include <math.h>

#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/lightcone_massmap_data.h"
#include "../data/macros.h"
#include "../data/mymalloc.h"
#include "../mpi_utils/setcomm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

class mmparticles : public setcomm
{
 public:
  mmparticles(MPI_Comm comm) : setcomm(comm) {}

  int NumPart; /**< number of particles on the LOCAL processor */
  int MaxPart; /**< This gives the maxmimum number of particles that can be stored on one processor. */

  int Npix; /* total number of pixels of Healpix tessellation of lightcone */
  int FirstPix;
  int NpixLoc;

  lightcone_massmap_data *P; /*!< holds the mass point particle data on local processor */

  void allocate_memory(void)
  {
    P = (lightcone_massmap_data *)Mem.mymalloc_movable_clear(&P, "P", MaxPart * sizeof(lightcone_massmap_data));
  }

  void reallocate_memory_maxpart(int maxpartNew)
  {
    MaxPart = maxpartNew;

    P = (lightcone_massmap_data *)Mem.myrealloc_movable(P, MaxPart * sizeof(lightcone_massmap_data));

    if(NumPart > MaxPart)
      Terminate("NumPart=%d > MaxPart=%d", NumPart, MaxPart);
  }
};

#endif

#endif
