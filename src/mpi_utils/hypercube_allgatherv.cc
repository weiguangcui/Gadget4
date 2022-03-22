/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file hypercube_allgatherv.cc
 *
 *  \brief a simple version of an Allgatherv implemented with a hypercube communication model
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

#define TAG 100

int myMPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int *recvcount, int *displs,
                     MPI_Datatype recvtype, MPI_Comm comm)
{
#ifndef MPI_HYPERCUBE_ALLGATHERV
  return MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcount, displs, recvtype, comm);
#else
  int ntask, thistask, ptask, ngrp, size_sendtype, size_recvtype;
  MPI_Status status;

  MPI_Comm_rank(comm, &thistask);
  MPI_Comm_size(comm, &ntask);

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  for(ptask = 0; ntask > (1 << ptask); ptask++)
    ;

  for(ngrp = 1; ngrp < (1 << ptask); ngrp++)
    {
      int recvtask = thistask ^ ngrp;

      if(recvtask < ntask)
        {
          if(sendbuf == MPI_IN_PLACE)
            myMPI_Sendrecv((char *)recvbuf + displs[thistask] * size_recvtype, recvcount[thistask], sendtype, recvtask, TAG,
                           (char *)recvbuf + displs[recvtask] * size_recvtype, recvcount[recvtask], recvtype, recvtask, TAG, comm,
                           &status);
          else
            myMPI_Sendrecv(sendbuf, sendcount, sendtype, recvtask, TAG, (char *)recvbuf + displs[recvtask] * size_recvtype,
                           recvcount[recvtask], recvtype, recvtask, TAG, comm, &status);
        }
    }

  if(sendbuf != MPI_IN_PLACE)
    if((char *)sendbuf != (char *)recvbuf + displs[thistask] * size_recvtype)
      memcpy((char *)recvbuf + displs[thistask] * size_recvtype, sendbuf, sendcount * size_sendtype);

  return 0;

#endif
}
