/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pm.h
 *
 *  \brief definition of a class to bundle the PM-force calculation algorithms
 */

#ifndef PM_H
#define PM_H

#if defined(PMGRID) || defined(NGENIC)

#include "gadgetconfig.h"

#include <fftw3.h>

typedef ptrdiff_t fft_ptrdiff_t;

#ifdef DOUBLEPRECISION_FFTW
typedef double fft_real;
typedef fftw_complex fft_complex;
#else
typedef float fft_real;
typedef fftwf_complex fft_complex;
#endif

#include "../mpi_utils/setcomm.h"
#include "../pm/pm_nonperiodic.h"
#include "../pm/pm_periodic.h"

class pm : public pm_periodic, public pm_nonperiodic
{
 public:
  pm(MPI_Comm comm) : setcomm(comm), pm_periodic(comm), pm_nonperiodic(comm) {}
};

#endif

#endif
