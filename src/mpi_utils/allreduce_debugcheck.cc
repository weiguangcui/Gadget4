/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  allreduce_debugcheck.cc
 *
 *  \brief some routines for cross-checking the use of collective MPI routines
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

int myMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  int mi, ma;

  MPI_Allreduce(&count, &mi, 1, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce(&count, &ma, 1, MPI_INT, MPI_MAX, comm);

  if(mi != ma)
    {
      int thistask, ntask;
      MPI_Comm_rank(comm, &thistask);
      MPI_Comm_size(comm, &ntask);

      printf("Error in MPI_Allreduce:  task=%d out of %d  has size = %d\n", thistask, ntask, count);
      fflush(stdout);
      MPI_Barrier(comm);

      Terminate("mi=%d ma=%d\n", mi, ma);
    }

  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}
