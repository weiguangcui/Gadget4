/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  mpi_vars.cc
 *
 *  \brief contains the global variables defined in the MPI helper functions
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../mpi_utils/mpi_utils.h"

MPI_Datatype MPI_MyIntPosType;

MPI_Op MPI_MIN_MyIntPosType;
MPI_Op MPI_MAX_MyIntPosType;
MPI_Op MPI_MIN_MySignedIntPosType;
MPI_Op MPI_MAX_MySignedIntPosType;
