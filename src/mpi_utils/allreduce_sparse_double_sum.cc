/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  allreduce_sparse_double_sum.cc
 *
 *  \brief implementation of a reduction operation for sparsely populated data
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

void allreduce_sparse_double_sum(double *loc, double *glob, int N, MPI_Comm Communicator)
{
  int ntask, thistask, ptask;
  MPI_Comm_size(Communicator, &ntask);
  MPI_Comm_rank(Communicator, &thistask);

  for(ptask = 0; ntask > (1 << ptask); ptask++)
    ;

  int *send_count  = (int *)Mem.mymalloc("send_count", sizeof(int) * ntask);
  int *recv_count  = (int *)Mem.mymalloc("recv_count", sizeof(int) * ntask);
  int *send_offset = (int *)Mem.mymalloc("send_offset", sizeof(int) * ntask);
  int *recv_offset = (int *)Mem.mymalloc("recv_offset", sizeof(int) * ntask);
  int *blocksize   = (int *)Mem.mymalloc("blocksize", sizeof(int) * ntask);

  int blk     = N / ntask;
  int rmd     = N - blk * ntask; /* remainder */
  int pivot_n = rmd * (blk + 1);

  int loc_first_n = 0;
  for(int task = 0; task < ntask; task++)
    {
      if(task < rmd)
        blocksize[task] = blk + 1;
      else
        blocksize[task] = blk;

      if(task < thistask)
        loc_first_n += blocksize[task];
    }

  double *loc_data = (double *)Mem.mymalloc("loc_data", blocksize[thistask] * sizeof(double));
  memset(loc_data, 0, blocksize[thistask] * sizeof(double));

  for(int j = 0; j < ntask; j++)
    send_count[j] = 0;

  /* find for each non-zero element the processor where it should go for being summed */
  for(int n = 0; n < N; n++)
    {
      if(loc[n] != 0)
        {
          int task;
          if(n < pivot_n)
            task = n / (blk + 1);
          else
            task = rmd + (n - pivot_n) / blk; /* note: if blk=0, then this case can not occur */

          send_count[task]++;
        }
    }

  myMPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, Communicator);

  int nimport = 0, nexport = 0;

  recv_offset[0] = 0, send_offset[0] = 0;

  for(int j = 0; j < ntask; j++)
    {
      nexport += send_count[j];
      nimport += recv_count[j];
      if(j > 0)
        {
          send_offset[j] = send_offset[j - 1] + send_count[j - 1];
          recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
        }
    }

  struct ind_data
  {
    int n;
    double val;
  };
  ind_data *export_data, *import_data;

  export_data = (ind_data *)Mem.mymalloc("export_data", nexport * sizeof(ind_data));
  import_data = (ind_data *)Mem.mymalloc("import_data", nimport * sizeof(ind_data));

  for(int j = 0; j < ntask; j++)
    send_count[j] = 0;

  for(int n = 0; n < N; n++)
    {
      if(loc[n] != 0)
        {
          int task;

          if(n < pivot_n)
            task = n / (blk + 1);
          else
            task = rmd + (n - pivot_n) / blk; /* note: if blk=0, then this case can not occur */

          int index              = send_offset[task] + send_count[task]++;
          export_data[index].n   = n;
          export_data[index].val = loc[n];
        }
    }

  for(int ngrp = 0; ngrp < (1 << ptask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = thistask ^ ngrp;
      if(recvTask < ntask)
        if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_data[send_offset[recvTask]], send_count[recvTask] * sizeof(ind_data), MPI_BYTE, recvTask, TAG_DENS_B,
                         &import_data[recv_offset[recvTask]], recv_count[recvTask] * sizeof(ind_data), MPI_BYTE, recvTask, TAG_DENS_B,
                         Communicator, MPI_STATUS_IGNORE);
    }

  for(int i = 0; i < nimport; i++)
    {
      int j = import_data[i].n - loc_first_n;

      if(j < 0 || j >= blocksize[thistask])
        Terminate("j=%d < 0 || j>= blocksize[thistask]=%d", j, blocksize[thistask]);

      loc_data[j] += import_data[i].val;
    }

  Mem.myfree(import_data);
  Mem.myfree(export_data);

  /* now share the cost data across all processors */
  int *bytecounts = (int *)Mem.mymalloc("bytecounts", sizeof(int) * ntask);
  int *byteoffset = (int *)Mem.mymalloc("byteoffset", sizeof(int) * ntask);

  for(int task = 0; task < ntask; task++)
    bytecounts[task] = blocksize[task] * sizeof(double);

  byteoffset[0] = 0;
  for(int task = 1; task < ntask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  myMPI_Allgatherv(loc_data, bytecounts[thistask], MPI_BYTE, glob, bytecounts, byteoffset, MPI_BYTE, Communicator);

  Mem.myfree(byteoffset);
  Mem.myfree(bytecounts);

  Mem.myfree(loc_data);
  Mem.myfree(blocksize);
  Mem.myfree(recv_offset);
  Mem.myfree(send_offset);
  Mem.myfree(recv_count);
  Mem.myfree(send_count);
}
