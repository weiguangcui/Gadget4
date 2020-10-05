/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  sums_and_minmax.cc
 *
 *  \brief some simple extensions of MPI-collectives
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
#include "../mpi_utils/mpi_utils.h"

void minimum_large_ints(int n, long long *src, long long *res, MPI_Comm comm)
{
  if(src == res)
    MPI_Allreduce(MPI_IN_PLACE, res, n, MPI_LONG_LONG, MPI_MIN, comm);
  else
    MPI_Allreduce(src, res, n, MPI_LONG_LONG, MPI_MIN, comm);
}

void sumup_large_ints(int n, int *src, long long *res, MPI_Comm comm)
{
  long long *numlist = (long long *)Mem.mymalloc("numlist", n * sizeof(long long));

  for(int j = 0; j < n; j++)
    numlist[j] = src[j];
  MPI_Allreduce(numlist, res, n, MPI_LONG_LONG, MPI_SUM, comm);

  Mem.myfree(numlist);
}

void sumup_longs(int n, long long *src, long long *res, MPI_Comm comm)
{
  if(src == res)
    MPI_Allreduce(MPI_IN_PLACE, res, n, MPI_LONG_LONG, MPI_SUM, comm);
  else
    MPI_Allreduce(src, res, n, MPI_LONG_LONG, MPI_SUM, comm);
}
