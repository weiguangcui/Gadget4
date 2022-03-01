/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  sizelimited_sendrecv.cc
 *
 *  \brief implements a wrapper around myMPI_Sendrecv that if needed transmits the data in smaller pieces than a prescribed maximum
 * size
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

int myMPI_Sendrecv(void *sendb, size_t sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvb, size_t recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status)
{
  int iter      = 0, size_sendtype, size_recvtype, send_now, recv_now;
  char *sendbuf = (char *)sendb;
  char *recvbuf = (char *)recvb;

  if(dest != source)
    Terminate("dest != source");

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  int thistask;
  MPI_Comm_rank(comm, &thistask);

  if(dest == thistask)
    {
      memcpy(recvbuf, sendbuf, recvcount * size_recvtype);
      return 0;
    }

  size_t count_limit = MPI_MESSAGE_SIZELIMIT_IN_BYTES / size_sendtype;

  while(sendcount > 0 || recvcount > 0)
    {
      if(sendcount > count_limit)
        {
          send_now = count_limit;

          iter++;
        }
      else
        send_now = sendcount;

      if(recvcount > count_limit)
        recv_now = count_limit;
      else
        recv_now = recvcount;

      MPI_Sendrecv(sendbuf, send_now, sendtype, dest, sendtag, recvbuf, recv_now, recvtype, source, recvtag, comm, status);

      sendcount -= send_now;
      recvcount -= recv_now;

      sendbuf += send_now * size_sendtype;
      recvbuf += recv_now * size_recvtype;
    }

  return 0;
}
