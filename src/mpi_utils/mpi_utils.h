/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  mpi_utils.h
 *  \brief declares some numerical values for MPI tags and function prototypes for MPI helper functions
 */

#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include "gadgetconfig.h"

#include <mpi.h>
#include "../data/dtypes.h"
#include "../data/mymalloc.h"

/*!< Various tags used for labeling MPI messages */
#define TAG_TOPNODE_FREE 4
#define TAG_TOPNODE_OFFSET 5
#define TAG_TOPNODE_ALLOC 6
#define TAG_EWALD_ALLOC 7
#define TAG_TABLE_ALLOC 8
#define TAG_TABLE_FREE 9
#define TAG_N 10
#define TAG_HEADER 11
#define TAG_PDATA 12
#define TAG_SPHDATA 13
#define TAG_KEY 14
#define TAG_DMOM 15
#define TAG_NODELEN 16
#define TAG_HMAX 17
#define TAG_GRAV_A 18
#define TAG_GRAV_B 19
#define TAG_DIRECT_A 20
#define TAG_DIRECT_B 21
#define TAG_HYDRO_A 22
#define TAG_HYDRO_B 23
#define TAG_NFORTHISTASK 24
#define TAG_PERIODIC_A 25
#define TAG_PERIODIC_B 26
#define TAG_PERIODIC_C 27
#define TAG_PERIODIC_D 28
#define TAG_NONPERIOD_A 29
#define TAG_NONPERIOD_B 30
#define TAG_NONPERIOD_C 31
#define TAG_NONPERIOD_D 32
#define TAG_POTENTIAL_A 33
#define TAG_POTENTIAL_B 34
#define TAG_DENS_A 35
#define TAG_DENS_B 36
#define TAG_LOCALN 37
#define TAG_BH_A 38
#define TAG_BH_B 39
#define TAG_SMOOTH_A 40
#define TAG_SMOOTH_B 41
#define TAG_ENRICH_A 42
#define TAG_CONDUCT_A 43
#define TAG_CONDUCT_B 44
#define TAG_FOF_A 45
#define TAG_FOF_B 46
#define TAG_FOF_C 47
#define TAG_FOF_D 48
#define TAG_FOF_E 49
#define TAG_FOF_F 50
#define TAG_FOF_G 51
#define TAG_HOTNGB_A 52
#define TAG_HOTNGB_B 53
#define TAG_GRAD_A 54
#define TAG_GRAD_B 55

#define TAG_SE 56

#define TAG_SEARCH_A 58
#define TAG_SEARCH_B 59

#define TAG_INJECT_A 61

#define TAG_PDATA_SPH 70
#define TAG_KEY_SPH 71

#define TAG_PDATA_STAR 72
#define TAG_STARDATA 73
#define TAG_KEY_STAR 74

#define TAG_PDATA_BH 75
#define TAG_BHDATA 76
#define TAG_KEY_BH 77

#define TAG_GRAVCOST_A 79
#define TAG_GRAVCOST_B 80

#define TAG_PM_FOLD 81

#define TAG_BARRIER 82
#define TAG_PART_DATA 83
#define TAG_NODE_DATA 84
#define TAG_RESULTS 85
#define TAG_DRIFT_INIT 86
#define TAG_ALL_UPDATE 87
#define TAG_METDATA 500
#define TAG_FETCH_GRAVTREE 1000
#define TAG_FETCH_SPH_DENSITY 2000
#define TAG_FETCH_SPH_HYDRO 3000
#define TAG_FETCH_SPH_TREETIMESTEP 4000

void my_mpi_types_init(void);

int myMPI_Sendrecv(void *sendbuf, size_t sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, size_t recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);

int myMPI_Alltoallv_new_prep(int *sendcnt, int *recvcnt, int *rdispls, MPI_Comm comm, int method);

void myMPI_Alltoallv_new(void *sendb, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void *recvb, int *recvcounts, int *rdispls,
                         MPI_Datatype recvtype, MPI_Comm comm, int method);

void myMPI_Alltoallv(void *sendbuf, size_t *sendcounts, size_t *sdispls, void *recvbuf, size_t *recvcounts, size_t *rdispls, int len,
                     int big_flag, MPI_Comm comm);

void my_int_MPI_Alltoallv(void *sendb, int *sendcounts, int *sdispls, void *recvb, int *recvcounts, int *rdispls, int len,
                          int big_flag, MPI_Comm comm);

int myMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int myMPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int *recvcount, int *displs,
                     MPI_Datatype recvtype, MPI_Comm comm);

int myMPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm);

void allreduce_sparse_double_sum(double *loc, double *glob, int N, MPI_Comm comm);

void minimum_large_ints(int n, long long *src, long long *res, MPI_Comm comm);
void sumup_longs(int n, long long *src, long long *res, MPI_Comm comm);
void sumup_large_ints(int n, int *src, long long *res, MPI_Comm comm);

extern MPI_Datatype MPI_MyIntPosType;

extern MPI_Op MPI_MIN_MyIntPosType;
extern MPI_Op MPI_MAX_MyIntPosType;
extern MPI_Op MPI_MIN_MySignedIntPosType;
extern MPI_Op MPI_MAX_MySignedIntPosType;

template <typename T>
void allreduce_sum(T *glob, int N, MPI_Comm Communicator)
{
  int ntask, thistask, ptask;
  MPI_Comm_size(Communicator, &ntask);
  MPI_Comm_rank(Communicator, &thistask);

  for(ptask = 0; ntask > (1 << ptask); ptask++)
    ;

  // we are responsible for a certain stretch of the result, namely the one starting at loc_first_n, of length blocksize[thistask]

  int *blocksize  = (int *)Mem.mymalloc("blocksize", sizeof(int) * ntask);
  int *blockstart = (int *)Mem.mymalloc("blockstart", sizeof(int) * ntask);

  int blk     = N / ntask;
  int rmd     = N - blk * ntask; /* remainder */
  int pivot_n = rmd * (blk + 1);

  int loc_first_n = 0;
  blockstart[0]   = 0;

  for(int task = 0; task < ntask; task++)
    {
      if(task < rmd)
        blocksize[task] = blk + 1;
      else
        blocksize[task] = blk;

      if(task < thistask)
        loc_first_n += blocksize[task];

      if(task > 0)
        blockstart[task] = blockstart[task - 1] + blocksize[task - 1];
    }

  /* here we store the local result */
  T *loc_data = (T *)Mem.mymalloc_clear("loc_data", blocksize[thistask] * sizeof(T));

  int *send_count = (int *)Mem.mymalloc("send_count", sizeof(int) * ntask);
  int *recv_count = (int *)Mem.mymalloc("recv_count", sizeof(int) * ntask);

  int *send_offset = (int *)Mem.mymalloc("send_offset", sizeof(int) * ntask);

  struct ind_data
  {
    int n;
    T val;
  };

  ind_data *export_data = NULL;
  int nexport           = 0;

  for(int rep = 0; rep < 2; rep++)
    {
      for(int j = 0; j < ntask; j++)
        send_count[j] = 0;

      /* find for each non-zero element the processor where it should go for being summed */
      for(int n = 0; n < N; n++)
        {
          if(glob[n] != 0)
            {
              int task;
              if(n < pivot_n)
                task = n / (blk + 1);
              else
                task = rmd + (n - pivot_n) / blk; /* note: if blk=0, then this case can not occur */

              if(rep == 0)
                send_count[task]++;
              else
                {
                  int index              = send_offset[task] + send_count[task]++;
                  export_data[index].n   = n;
                  export_data[index].val = glob[n];
                }
            }
        }

      if(rep == 0)
        {
          myMPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, Communicator);

          send_offset[0] = 0;

          for(int j = 0; j < ntask; j++)
            {
              nexport += send_count[j];

              if(j > 0)
                send_offset[j] = send_offset[j - 1] + send_count[j - 1];
            }

          export_data = (ind_data *)Mem.mymalloc("export_data", nexport * sizeof(ind_data));
        }
      else
        {
          for(int ngrp = 0; ngrp < (1 << ptask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
            {
              int recvTask = thistask ^ ngrp;
              if(recvTask < ntask)
                if(send_count[recvTask] > 0 || recv_count[recvTask] > 0)
                  {
                    int nimport = recv_count[recvTask];

                    ind_data *import_data = (ind_data *)Mem.mymalloc("import_data", nimport * sizeof(ind_data));

                    myMPI_Sendrecv(&export_data[send_offset[recvTask]], send_count[recvTask] * sizeof(ind_data), MPI_BYTE, recvTask,
                                   TAG_DENS_B, import_data, recv_count[recvTask] * sizeof(ind_data), MPI_BYTE, recvTask, TAG_DENS_B,
                                   Communicator, MPI_STATUS_IGNORE);

                    for(int i = 0; i < nimport; i++)
                      {
                        int j = import_data[i].n - loc_first_n;

                        if(j < 0 || j >= blocksize[thistask])
                          Terminate("j=%d < 0 || j>= blocksize[thistask]=%d", j, blocksize[thistask]);

                        loc_data[j] += import_data[i].val;
                      }

                    Mem.myfree(import_data);
                  }
            }

          Mem.myfree(export_data);
        }
    }

  Mem.myfree(send_offset);
  Mem.myfree(recv_count);
  Mem.myfree(send_count);

  /* now share the result across all processors */
  for(int ngrp = 0; ngrp < (1 << ptask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = thistask ^ ngrp;
      if(recvTask < ntask)
        if(blocksize[thistask] > 0 || blocksize[recvTask] > 0)
          myMPI_Sendrecv(loc_data, blocksize[thistask] * sizeof(T), MPI_BYTE, recvTask, TAG_DENS_A, &glob[blockstart[recvTask]],
                         blocksize[recvTask] * sizeof(T), MPI_BYTE, recvTask, TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);
    }

  Mem.myfree(loc_data);
  Mem.myfree(blockstart);
  Mem.myfree(blocksize);
}

#endif
