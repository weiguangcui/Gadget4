/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pm_mpi_fft.cc
 *
 *  \brief code for doing different variants of parallel FFT transforms
 */

#include "gadgetconfig.h"

#if defined(PMGRID) || defined(NGENIC)

#include <fftw3.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../system/system.h"

/* We only use the one-dimensional FFTW3 routines, because the MPI versions of FFTW3
 * allocated memory for themselves during the transforms (which we want to strictly avoid),
 * and because we want to allow transforms that are so big that more than 2GB may be
 * transferred betweeen processors.
 */

#ifndef FFT_COLUMN_BASED

void pm_mpi_fft::my_slab_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ)
{
  subdivide_evenly(NgridX, NTask, ThisTask, &plan->slabstart_x, &plan->nslab_x);
  subdivide_evenly(NgridY, NTask, ThisTask, &plan->slabstart_y, &plan->nslab_y);

  plan->slab_to_task = (int *)Mem.mymalloc_movable(&plan->slab_to_task, "slab_to_task", NgridX * sizeof(int));

  for(int task = 0; task < NTask; task++)
    {
      int start, n;

      subdivide_evenly(NgridX, NTask, task, &start, &n);

      for(int i = start; i < start + n; i++)
        plan->slab_to_task[i] = task;
    }

  MPI_Allreduce(&plan->nslab_x, &plan->largest_x_slab, 1, MPI_INT, MPI_MAX, Communicator);
  MPI_Allreduce(&plan->nslab_y, &plan->largest_y_slab, 1, MPI_INT, MPI_MAX, Communicator);

  plan->slabs_x_per_task = (int *)Mem.mymalloc_movable(&plan->slabs_x_per_task, "slabs_x_per_task", NTask * sizeof(int));
  MPI_Allgather(&plan->nslab_x, 1, MPI_INT, plan->slabs_x_per_task, 1, MPI_INT, Communicator);

  plan->first_slab_x_of_task = (int *)Mem.mymalloc_movable(&plan->first_slab_x_of_task, "first_slab_x_of_task", NTask * sizeof(int));
  MPI_Allgather(&plan->slabstart_x, 1, MPI_INT, plan->first_slab_x_of_task, 1, MPI_INT, Communicator);

  plan->slabs_y_per_task = (int *)Mem.mymalloc_movable(&plan->slabs_y_per_task, "slabs_y_per_task", NTask * sizeof(int));
  MPI_Allgather(&plan->nslab_y, 1, MPI_INT, plan->slabs_y_per_task, 1, MPI_INT, Communicator);

  plan->first_slab_y_of_task = (int *)Mem.mymalloc_movable(&plan->first_slab_y_of_task, "first_slab_y_of_task", NTask * sizeof(int));
  MPI_Allgather(&plan->slabstart_y, 1, MPI_INT, plan->first_slab_y_of_task, 1, MPI_INT, Communicator);

  plan->NgridX = NgridX;
  plan->NgridY = NgridY;
  plan->NgridZ = NgridZ;

  int Ngridz = NgridZ / 2 + 1; /* dimension needed in complex space */

  plan->Ngridz = Ngridz;
  plan->Ngrid2 = 2 * Ngridz;
}

void pm_mpi_fft::my_slab_based_fft_free(fft_plan *plan)
{
  Mem.myfree(plan->first_slab_y_of_task);
  Mem.myfree(plan->slabs_y_per_task);
  Mem.myfree(plan->first_slab_x_of_task);
  Mem.myfree(plan->slabs_x_per_task);
  Mem.myfree(plan->slab_to_task);
}

/*! \brief Transposes the array field
 *
 * The array field is transposed such that the data in x direction is local to only one task.
 * This is done, so the force in x-direction can be obtained by finite differencing.
 * However the array is not fully transposed, i.e. the x-direction is not the fastest
 * running array index
 *
 * \param field The array to transpose
 * \param scratch scratch space used during communication (same size as field)
 */
void pm_mpi_fft::my_slab_transposeA(fft_plan *plan, fft_real *field, fft_real *scratch)
{
  int n, prod, task, flag_big = 0, flag_big_all = 0;

  prod = NTask * plan->nslab_x;

  for(n = 0; n < prod; n++)
    {
      int x    = n / NTask;
      int task = n % NTask;

      int y;

      for(y = plan->first_slab_y_of_task[task]; y < plan->first_slab_y_of_task[task] + plan->slabs_y_per_task[task]; y++)
        memcpy(scratch + ((size_t)plan->NgridZ) * (plan->first_slab_y_of_task[task] * plan->nslab_x +
                                                   x * plan->slabs_y_per_task[task] + (y - plan->first_slab_y_of_task[task])),
               field + ((size_t)plan->Ngrid2) * (plan->NgridY * x + y), plan->NgridZ * sizeof(fft_real));
    }

  size_t *scount = (size_t *)Mem.mymalloc("scount", NTask * sizeof(size_t));
  size_t *rcount = (size_t *)Mem.mymalloc("rcount", NTask * sizeof(size_t));
  size_t *soff   = (size_t *)Mem.mymalloc("soff", NTask * sizeof(size_t));
  size_t *roff   = (size_t *)Mem.mymalloc("roff", NTask * sizeof(size_t));

  for(task = 0; task < NTask; task++)
    {
      scount[task] = plan->nslab_x * plan->slabs_y_per_task[task] * (plan->NgridZ * sizeof(fft_real));
      rcount[task] = plan->nslab_y * plan->slabs_x_per_task[task] * (plan->NgridZ * sizeof(fft_real));

      soff[task] = plan->first_slab_y_of_task[task] * plan->nslab_x * (plan->NgridZ * sizeof(fft_real));
      roff[task] = plan->first_slab_x_of_task[task] * plan->nslab_y * (plan->NgridZ * sizeof(fft_real));

      if(scount[task] > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
        flag_big = 1;
    }

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  myMPI_Alltoallv(scratch, scount, soff, field, rcount, roff, 1, flag_big_all, Communicator);

  Mem.myfree(roff);
  Mem.myfree(soff);
  Mem.myfree(rcount);
  Mem.myfree(scount);
}

/*! \brief Undo the transposition of the array field
 *
 * The transposition of the array field is undone such that the data in
 * x direction is distributed among all tasks again. Thus the result of
 * force computation in x-direction is sent back to the original task.
 *
 * \param field The array to transpose
 * \param scratch scratch space used during communication (same size as field)
 */
void pm_mpi_fft::my_slab_transposeB(fft_plan *plan, fft_real *field, fft_real *scratch)
{
  int n, prod, task, flag_big = 0, flag_big_all = 0;

  size_t *scount = (size_t *)Mem.mymalloc("scount", NTask * sizeof(size_t));
  size_t *rcount = (size_t *)Mem.mymalloc("rcount", NTask * sizeof(size_t));
  size_t *soff   = (size_t *)Mem.mymalloc("soff", NTask * sizeof(size_t));
  size_t *roff   = (size_t *)Mem.mymalloc("roff", NTask * sizeof(size_t));

  for(task = 0; task < NTask; task++)
    {
      rcount[task] = plan->nslab_x * plan->slabs_y_per_task[task] * (plan->NgridZ * sizeof(fft_real));
      scount[task] = plan->nslab_y * plan->slabs_x_per_task[task] * (plan->NgridZ * sizeof(fft_real));

      roff[task] = plan->first_slab_y_of_task[task] * plan->nslab_x * (plan->NgridZ * sizeof(fft_real));
      soff[task] = plan->first_slab_x_of_task[task] * plan->nslab_y * (plan->NgridZ * sizeof(fft_real));

      if(scount[task] > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
        flag_big = 1;
    }

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  myMPI_Alltoallv(field, scount, soff, scratch, rcount, roff, 1, flag_big_all, Communicator);

  Mem.myfree(roff);
  Mem.myfree(soff);
  Mem.myfree(rcount);
  Mem.myfree(scount);

  prod = NTask * plan->nslab_x;

  for(n = 0; n < prod; n++)
    {
      int x    = n / NTask;
      int task = n % NTask;

      int y;
      for(y = plan->first_slab_y_of_task[task]; y < plan->first_slab_y_of_task[task] + plan->slabs_y_per_task[task]; y++)
        memcpy(field + ((size_t)plan->Ngrid2) * (plan->NgridY * x + y),
               scratch + ((size_t)plan->NgridZ) * (plan->first_slab_y_of_task[task] * plan->nslab_x +
                                                   x * plan->slabs_y_per_task[task] + (y - plan->first_slab_y_of_task[task])),
               plan->NgridZ * sizeof(fft_real));
    }
}

/* Given a slab-decomposed 3D field a[...] with total dimension [nx x ny x nz], whose first dimension is
 * split across the processors, this routine outputs in b[] the transpose where then the second dimension is split
 * across the processors. sx[] gives for each MPI task how many slabs it has, and firstx[] is the first
 * slab for a given task. Likewise, sy[]/firsty[] gives the same thing for the transposed order. Note, the
 * contents of the array a[] will be destroyed by the routine.
 *
 * An element (x,y,z) is accessed in a[] with index [([x - firstx] * ny + y) * nz + z]
 * and in b[] as [((y - firsty) * nx + x) * nz + z]
 *
 * if mode = 1, the reverse operation is carried out.
 */
void pm_mpi_fft::my_slab_transpose(void *av, void *bv, int *sx, int *firstx, int *sy, int *firsty, int nx, int ny, int nz, int mode)
{
  char *a = (char *)av;
  char *b = (char *)bv;

  size_t *scount = (size_t *)Mem.mymalloc("scount", NTask * sizeof(size_t));
  size_t *rcount = (size_t *)Mem.mymalloc("rcount", NTask * sizeof(size_t));
  size_t *soff   = (size_t *)Mem.mymalloc("soff", NTask * sizeof(size_t));
  size_t *roff   = (size_t *)Mem.mymalloc("roff", NTask * sizeof(size_t));
  int i, n, prod, flag_big = 0, flag_big_all = 0;

  for(i = 0; i < NTask; i++)
    {
      scount[i] = sy[i] * sx[ThisTask] * ((size_t)nz);
      rcount[i] = sy[ThisTask] * sx[i] * ((size_t)nz);
      soff[i]   = firsty[i] * sx[ThisTask] * ((size_t)nz);
      roff[i]   = sy[ThisTask] * firstx[i] * ((size_t)nz);

      if(scount[i] * sizeof(fft_complex) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
        flag_big = 1;
    }

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  if(mode == 0)
    {
      /* first pack the data into contiguous blocks */
      prod = NTask * sx[ThisTask];
      for(n = 0; n < prod; n++)
        {
          int k = n / NTask;
          int i = n % NTask;
          int j;

          for(j = 0; j < sy[i]; j++)
            memcpy(b + (k * sy[i] + j + firsty[i] * sx[ThisTask]) * (nz * sizeof(fft_complex)),
                   a + (k * ny + (firsty[i] + j)) * (nz * sizeof(fft_complex)), nz * sizeof(fft_complex));
        }

      /* tranfer the data */
      myMPI_Alltoallv(b, scount, soff, a, rcount, roff, sizeof(fft_complex), flag_big_all, Communicator);

      /* unpack the data into the right order */
      prod = NTask * sy[ThisTask];
      for(n = 0; n < prod; n++)
        {
          int j = n / NTask;
          int i = n % NTask;
          int k;

          for(k = 0; k < sx[i]; k++)
            memcpy(b + (j * nx + k + firstx[i]) * (nz * sizeof(fft_complex)),
                   a + ((k + firstx[i]) * sy[ThisTask] + j) * (nz * sizeof(fft_complex)), nz * sizeof(fft_complex));
        }
    }
  else
    {
      /* first pack the data into contiguous blocks */
      prod = NTask * sy[ThisTask];
      for(n = 0; n < prod; n++)
        {
          int j = n / NTask;
          int i = n % NTask;
          int k;

          for(k = 0; k < sx[i]; k++)
            memcpy(b + ((k + firstx[i]) * sy[ThisTask] + j) * (nz * sizeof(fft_complex)),
                   a + (j * nx + k + firstx[i]) * (nz * sizeof(fft_complex)), nz * sizeof(fft_complex));
        }

      /* tranfer the data */
      myMPI_Alltoallv(b, rcount, roff, a, scount, soff, sizeof(fft_complex), flag_big_all, Communicator);

      /* unpack the data into the right order */
      prod = NTask * sx[ThisTask];
      for(n = 0; n < prod; n++)
        {
          int k = n / NTask;
          int i = n % NTask;
          int j;

          for(j = 0; j < sy[i]; j++)
            memcpy(b + (k * ny + (firsty[i] + j)) * (nz * sizeof(fft_complex)),
                   a + (k * sy[i] + j + firsty[i] * sx[ThisTask]) * (nz * sizeof(fft_complex)), nz * sizeof(fft_complex));
        }
    }

  /* now the result is in b[] */

  Mem.myfree(roff);
  Mem.myfree(soff);
  Mem.myfree(rcount);
  Mem.myfree(scount);
}

void pm_mpi_fft::my_slab_based_fft(fft_plan *plan, void *data, void *workspace, int forward)
{
  int n, prod;
  int slabsx = plan->slabs_x_per_task[ThisTask];
  int slabsy = plan->slabs_y_per_task[ThisTask];

  int ngridx  = plan->NgridX;
  int ngridy  = plan->NgridY;
  int ngridz  = plan->Ngridz;
  int ngridz2 = 2 * ngridz;

  size_t ngridx_long  = ngridx;
  size_t ngridy_long  = ngridy;
  size_t ngridz_long  = ngridz;
  size_t ngridz2_long = ngridz2;

  fft_real *data_real       = (fft_real *)data;
  fft_complex *data_complex = (fft_complex *)data, *workspace_complex = (fft_complex *)workspace;

  if(forward == 1)
    {
      /* do the z-direction FFT, real to complex */
      prod = slabsx * ngridy;
      for(n = 0; n < prod; n++)
        {
          FFTW(execute_dft_r2c)(plan->forward_plan_zdir, data_real + n * ngridz2_long, workspace_complex + n * ngridz_long);
        }

      /* do the y-direction FFT, complex to complex */
      prod = slabsx * ngridz;
      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->forward_plan_ydir, workspace_complex + i * ngridz * ngridy_long + j, data_complex + i * ngridz * ngridy_long + j);
        }

      /* now our data resides in data_complex[] */

      /* do the transpose */
      my_slab_transpose(data_complex, workspace_complex, plan->slabs_x_per_task, plan->first_slab_x_of_task, plan->slabs_y_per_task,
                        plan->first_slab_y_of_task, ngridx, ngridy, ngridz, 0);

      /* now the data is in workspace_complex[] */

      /* finally, do the transform along the x-direction (we are in transposed order, x and y have interchanged */
      prod = slabsy * ngridz;
      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->forward_plan_xdir, workspace_complex + i * ngridz * ngridx_long + j, data_complex + i * ngridz * ngridx_long + j);
        }

      /* now the result is in data_complex[] */
    }
  else
    {
      prod = slabsy * ngridz;

      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->backward_plan_xdir, data_complex + i * ngridz * ngridx_long + j, workspace_complex + i * ngridz * ngridx_long + j);
        }

      my_slab_transpose(workspace_complex, data_complex, plan->slabs_x_per_task, plan->first_slab_x_of_task, plan->slabs_y_per_task,
                        plan->first_slab_y_of_task, ngridx, ngridy, ngridz, 1);

      prod = slabsx * ngridz;

      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->backward_plan_ydir, data_complex + i * ngridz * ngridy_long + j, workspace_complex + i * ngridz * ngridy_long + j);
        }

      prod = slabsx * ngridy;

      for(n = 0; n < prod; n++)
        {
          FFTW(execute_dft_c2r)(plan->backward_plan_zdir, workspace_complex + n * ngridz_long, data_real + n * ngridz2_long);
        }

      /* now the result is in data[] */
    }
}

#else

void pm_mpi_fft::my_column_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ)
{
  plan->NgridX = NgridX;
  plan->NgridY = NgridY;
  plan->NgridZ = NgridZ;

  int Ngridz = NgridZ / 2 + 1;

  plan->Ngridz = Ngridz;
  plan->Ngrid2 = 2 * Ngridz;

  subdivide_evenly(plan->NgridX * plan->NgridY, NTask, ThisTask, &plan->firstcol_XY, &plan->ncol_XY);
  subdivide_evenly(plan->NgridX * plan->Ngrid2, NTask, ThisTask, &plan->firstcol_XZ, &plan->ncol_XZ);
  subdivide_evenly(plan->Ngrid2 * plan->NgridY, NTask, ThisTask, &plan->firstcol_ZY, &plan->ncol_ZY);

  plan->lastcol_XY = plan->firstcol_XY + plan->ncol_XY - 1;
  plan->lastcol_XZ = plan->firstcol_XZ + plan->ncol_XZ - 1;
  plan->lastcol_ZY = plan->firstcol_ZY + plan->ncol_ZY - 1;

  subdivide_evenly(NgridX * Ngridz, NTask, ThisTask, &plan->transposed_firstcol, &plan->transposed_ncol);
  subdivide_evenly(NgridY * Ngridz, NTask, ThisTask, &plan->second_transposed_firstcol, &plan->second_transposed_ncol);

  plan->second_transposed_ncells = ((size_t)plan->NgridX) * plan->second_transposed_ncol;

  plan->max_datasize = ((size_t)plan->Ngrid2) * plan->ncol_XY;
  plan->max_datasize = std::max<size_t>(plan->max_datasize, 2 * ((size_t)plan->NgridY) * plan->transposed_ncol);
  plan->max_datasize = std::max<size_t>(plan->max_datasize, 2 * ((size_t)plan->NgridX) * plan->second_transposed_ncol);
  plan->max_datasize = std::max<size_t>(plan->max_datasize, ((size_t)plan->ncol_XZ) * plan->NgridY);
  plan->max_datasize = std::max<size_t>(plan->max_datasize, ((size_t)plan->ncol_ZY) * plan->NgridX);

  plan->fftsize = plan->max_datasize;

  plan->offsets_send_A = (size_t *)Mem.mymalloc_movable_clear(&plan->offsets_send_A, "offsets_send_A", NTask * sizeof(size_t));
  plan->offsets_recv_A = (size_t *)Mem.mymalloc_movable_clear(&plan->offsets_recv_A, "offsets_recv_A", NTask * sizeof(size_t));
  plan->offsets_send_B = (size_t *)Mem.mymalloc_movable_clear(&plan->offsets_send_B, "offsets_send_B", NTask * sizeof(size_t));
  plan->offsets_recv_B = (size_t *)Mem.mymalloc_movable_clear(&plan->offsets_recv_B, "offsets_recv_B", NTask * sizeof(size_t));
  plan->offsets_send_C = (size_t *)Mem.mymalloc_movable_clear(&plan->offsets_send_C, "offsets_send_C", NTask * sizeof(size_t));
  plan->offsets_recv_C = (size_t *)Mem.mymalloc_movable_clear(&plan->offsets_recv_C, "offsets_recv_C", NTask * sizeof(size_t));
  plan->offsets_send_D = (size_t *)Mem.mymalloc_movable_clear(&plan->offsets_send_D, "offsets_send_D", NTask * sizeof(size_t));
  plan->offsets_recv_D = (size_t *)Mem.mymalloc_movable_clear(&plan->offsets_recv_D, "offsets_recv_D", NTask * sizeof(size_t));

  plan->count_send_A  = (size_t *)Mem.mymalloc_movable_clear(&plan->count_send_A, "count_send_A", NTask * sizeof(size_t));
  plan->count_recv_A  = (size_t *)Mem.mymalloc_movable_clear(&plan->count_recv_A, "count_recv_A", NTask * sizeof(size_t));
  plan->count_send_B  = (size_t *)Mem.mymalloc_movable_clear(&plan->count_send_B, "count_send_B", NTask * sizeof(size_t));
  plan->count_recv_B  = (size_t *)Mem.mymalloc_movable_clear(&plan->count_recv_B, "count_recv_B", NTask * sizeof(size_t));
  plan->count_send_C  = (size_t *)Mem.mymalloc_movable_clear(&plan->count_send_C, "count_send_C", NTask * sizeof(size_t));
  plan->count_recv_C  = (size_t *)Mem.mymalloc_movable_clear(&plan->count_recv_C, "count_recv_C", NTask * sizeof(size_t));
  plan->count_send_D  = (size_t *)Mem.mymalloc_movable_clear(&plan->count_send_D, "count_send_D", NTask * sizeof(size_t));
  plan->count_recv_D  = (size_t *)Mem.mymalloc_movable_clear(&plan->count_recv_D, "count_recv_D", NTask * sizeof(size_t));
  plan->count_send_13 = (size_t *)Mem.mymalloc_movable_clear(&plan->count_send_13, "count_send_13", NTask * sizeof(size_t));
  plan->count_recv_13 = (size_t *)Mem.mymalloc_movable_clear(&plan->count_recv_13, "count_recv_13", NTask * sizeof(size_t));
  plan->count_send_23 = (size_t *)Mem.mymalloc_movable_clear(&plan->count_send_23, "count_send_23", NTask * sizeof(size_t));
  plan->count_recv_23 = (size_t *)Mem.mymalloc_movable_clear(&plan->count_recv_23, "count_recv_23", NTask * sizeof(size_t));
  plan->count_send_13back =
      (size_t *)Mem.mymalloc_movable_clear(&plan->count_send_13back, "count_send_13back", NTask * sizeof(size_t));
  plan->count_recv_13back =
      (size_t *)Mem.mymalloc_movable_clear(&plan->count_recv_13back, "count_recv_13back", NTask * sizeof(size_t));
  plan->count_send_23back =
      (size_t *)Mem.mymalloc_movable_clear(&plan->count_send_23back, "count_send_23back", NTask * sizeof(size_t));
  plan->count_recv_23back =
      (size_t *)Mem.mymalloc_movable_clear(&plan->count_recv_23back, "count_recv_23back", NTask * sizeof(size_t));

  int dimA[3]  = {plan->NgridX, plan->NgridY, plan->Ngridz};
  int permA[3] = {0, 2, 1};

  my_fft_column_remap(NULL, dimA, plan->firstcol_XY, plan->ncol_XY, NULL, permA, plan->transposed_firstcol, plan->transposed_ncol,
                      plan->offsets_send_A, plan->offsets_recv_A, plan->count_send_A, plan->count_recv_A, 1);

  int dimB[3]  = {plan->NgridX, plan->Ngridz, plan->NgridY};
  int permB[3] = {2, 1, 0};

  my_fft_column_remap(NULL, dimB, plan->transposed_firstcol, plan->transposed_ncol, NULL, permB, plan->second_transposed_firstcol,
                      plan->second_transposed_ncol, plan->offsets_send_B, plan->offsets_recv_B, plan->count_send_B, plan->count_recv_B,
                      1);

  int dimC[3]  = {plan->NgridY, plan->Ngridz, plan->NgridX};
  int permC[3] = {2, 1, 0};

  my_fft_column_remap(NULL, dimC, plan->second_transposed_firstcol, plan->second_transposed_ncol, NULL, permC,
                      plan->transposed_firstcol, plan->transposed_ncol, plan->offsets_send_C, plan->offsets_recv_C, plan->count_send_C,
                      plan->count_recv_C, 1);

  int dimD[3]  = {plan->NgridX, plan->Ngridz, plan->NgridY};
  int permD[3] = {0, 2, 1};

  my_fft_column_remap(NULL, dimD, plan->transposed_firstcol, plan->transposed_ncol, NULL, permD, plan->firstcol_XY, plan->ncol_XY,
                      plan->offsets_send_D, plan->offsets_recv_D, plan->count_send_D, plan->count_recv_D, 1);

  int dim23[3]  = {plan->NgridX, plan->NgridY, plan->Ngrid2};
  int perm23[3] = {0, 2, 1};

  my_fft_column_transpose(NULL, dim23, plan->firstcol_XY, plan->ncol_XY, NULL, perm23, plan->firstcol_XZ, plan->ncol_XZ,
                          plan->count_send_23, plan->count_recv_23, 1);

  int dim23back[3]  = {plan->NgridX, plan->Ngrid2, plan->NgridY};
  int perm23back[3] = {0, 2, 1};

  my_fft_column_transpose(NULL, dim23back, plan->firstcol_XZ, plan->ncol_XZ, NULL, perm23back, plan->firstcol_XY, plan->ncol_XY,
                          plan->count_send_23back, plan->count_recv_23back, 1);

  int dim13[3]  = {plan->NgridX, plan->NgridY, plan->Ngrid2};
  int perm13[3] = {2, 1, 0};

  my_fft_column_transpose(NULL, dim13, plan->firstcol_XY, plan->ncol_XY, NULL, perm13, plan->firstcol_ZY, plan->ncol_ZY,
                          plan->count_send_13, plan->count_recv_13, 1);

  int dim13back[3]  = {plan->Ngrid2, plan->NgridY, plan->NgridX};
  int perm13back[3] = {2, 1, 0};

  my_fft_column_transpose(NULL, dim13back, plan->firstcol_ZY, plan->ncol_ZY, NULL, perm13back, plan->firstcol_XY, plan->ncol_XY,
                          plan->count_send_13back, plan->count_recv_13back, 1);
}

void pm_mpi_fft::my_column_based_fft_free(fft_plan *plan)
{
  Mem.myfree(plan->count_recv_23back);
  Mem.myfree(plan->count_send_23back);
  Mem.myfree(plan->count_recv_13back);
  Mem.myfree(plan->count_send_13back);
  Mem.myfree(plan->count_recv_23);
  Mem.myfree(plan->count_send_23);
  Mem.myfree(plan->count_recv_13);
  Mem.myfree(plan->count_send_13);
  Mem.myfree(plan->count_recv_D);
  Mem.myfree(plan->count_send_D);
  Mem.myfree(plan->count_recv_C);
  Mem.myfree(plan->count_send_C);
  Mem.myfree(plan->count_recv_B);
  Mem.myfree(plan->count_send_B);
  Mem.myfree(plan->count_recv_A);
  Mem.myfree(plan->count_send_A);

  Mem.myfree(plan->offsets_recv_D);
  Mem.myfree(plan->offsets_send_D);
  Mem.myfree(plan->offsets_recv_C);
  Mem.myfree(plan->offsets_send_C);
  Mem.myfree(plan->offsets_recv_B);
  Mem.myfree(plan->offsets_send_B);
  Mem.myfree(plan->offsets_recv_A);
  Mem.myfree(plan->offsets_send_A);
}

void pm_mpi_fft::my_fft_swap23(fft_plan *plan, fft_real *data, fft_real *out)
{
  int dim23[3]  = {plan->NgridX, plan->NgridY, plan->Ngrid2};
  int perm23[3] = {0, 2, 1};

  my_fft_column_transpose(data, dim23, plan->firstcol_XY, plan->ncol_XY, out, perm23, plan->firstcol_XZ, plan->ncol_XZ,
                          plan->count_send_23, plan->count_recv_23, 0);
}

void pm_mpi_fft::my_fft_swap23back(fft_plan *plan, fft_real *data, fft_real *out)
{
  int dim23back[3]  = {plan->NgridX, plan->Ngrid2, plan->NgridY};
  int perm23back[3] = {0, 2, 1};

  my_fft_column_transpose(data, dim23back, plan->firstcol_XZ, plan->ncol_XZ, out, perm23back, plan->firstcol_XY, plan->ncol_XY,
                          plan->count_send_23back, plan->count_recv_23back, 0);
}

void pm_mpi_fft::my_fft_swap13(fft_plan *plan, fft_real *data, fft_real *out)
{
  int dim13[3]  = {plan->NgridX, plan->NgridY, plan->Ngrid2};
  int perm13[3] = {2, 1, 0};

  my_fft_column_transpose(data, dim13, plan->firstcol_XY, plan->ncol_XY, out, perm13, plan->firstcol_ZY, plan->ncol_ZY,
                          plan->count_send_13, plan->count_recv_13, 0);
}

void pm_mpi_fft::my_fft_swap13back(fft_plan *plan, fft_real *data, fft_real *out)
{
  int dim13back[3]  = {plan->Ngrid2, plan->NgridY, plan->NgridX};
  int perm13back[3] = {2, 1, 0};

  my_fft_column_transpose(data, dim13back, plan->firstcol_ZY, plan->ncol_ZY, out, perm13back, plan->firstcol_XY, plan->ncol_XY,
                          plan->count_send_13back, plan->count_recv_13back, 0);
}

void pm_mpi_fft::my_column_based_fft(fft_plan *plan, void *data, void *workspace, int forward)
{
  long long n;
  fft_real *data_real = (fft_real *)data, *workspace_real = (fft_real *)workspace;
  fft_complex *data_complex = (fft_complex *)data, *workspace_complex = (fft_complex *)workspace;

  if(forward == 1)
    {
      /* do the z-direction FFT, real to complex */
      for(n = 0; n < plan->ncol_XY; n++)
        FFTW(execute_dft_r2c)(plan->forward_plan_zdir, data_real + n * plan->Ngrid2, workspace_complex + n * plan->Ngridz);

      int dimA[3]  = {plan->NgridX, plan->NgridY, plan->Ngridz};
      int permA[3] = {0, 2, 1};

      my_fft_column_remap(workspace_complex, dimA, plan->firstcol_XY, plan->ncol_XY, data_complex, permA, plan->transposed_firstcol,
                          plan->transposed_ncol, plan->offsets_send_A, plan->offsets_recv_A, plan->count_send_A, plan->count_recv_A,
                          0);

      /* do the y-direction FFT in 'data', complex to complex */
      for(n = 0; n < plan->transposed_ncol; n++)
        FFTW(execute_dft)(plan->forward_plan_ydir, data_complex + n * plan->NgridY, workspace_complex + n * plan->NgridY);

      int dimB[3]  = {plan->NgridX, plan->Ngridz, plan->NgridY};
      int permB[3] = {2, 1, 0};

      my_fft_column_remap(workspace_complex, dimB, plan->transposed_firstcol, plan->transposed_ncol, data_complex, permB,
                          plan->second_transposed_firstcol, plan->second_transposed_ncol, plan->offsets_send_B, plan->offsets_recv_B,
                          plan->count_send_B, plan->count_recv_B, 0);

      /* do the x-direction FFT in 'data', complex to complex */
      for(n = 0; n < plan->second_transposed_ncol; n++)
        FFTW(execute_dft)(plan->forward_plan_xdir, data_complex + n * plan->NgridX, workspace_complex + n * plan->NgridX);

      /* result is now in workspace */
    }
  else
    {
      /* do inverse FFT in 'data' */
      for(n = 0; n < plan->second_transposed_ncol; n++)
        FFTW(execute_dft)(plan->backward_plan_xdir, data_complex + n * plan->NgridX, workspace_complex + n * plan->NgridX);

      int dimC[3]  = {plan->NgridY, plan->Ngridz, plan->NgridX};
      int permC[3] = {2, 1, 0};

      my_fft_column_remap(workspace_complex, dimC, plan->second_transposed_firstcol, plan->second_transposed_ncol, data_complex, permC,
                          plan->transposed_firstcol, plan->transposed_ncol, plan->offsets_send_C, plan->offsets_recv_C,
                          plan->count_send_C, plan->count_recv_C, 0);

      /* do inverse FFT in 'data' */
      for(n = 0; n < plan->transposed_ncol; n++)
        FFTW(execute_dft)(plan->backward_plan_ydir, data_complex + n * plan->NgridY, workspace_complex + n * plan->NgridY);

      int dimD[3]  = {plan->NgridX, plan->Ngridz, plan->NgridY};
      int permD[3] = {0, 2, 1};

      my_fft_column_remap(workspace_complex, dimD, plan->transposed_firstcol, plan->transposed_ncol, data_complex, permD,
                          plan->firstcol_XY, plan->ncol_XY, plan->offsets_send_D, plan->offsets_recv_D, plan->count_send_D,
                          plan->count_recv_D, 0);

      /* do complex-to-real inverse transform on z-coordinates */
      for(n = 0; n < plan->ncol_XY; n++)
        FFTW(execute_dft_c2r)(plan->backward_plan_zdir, data_complex + n * plan->Ngridz, workspace_real + n * plan->Ngrid2);
    }
}

void pm_mpi_fft::my_fft_column_remap(fft_complex *data, int Ndims[3], /* global dimensions of data cube */
                                     int in_firstcol, int in_ncol,    /* first column and number of columns */
                                     fft_complex *out, int perm[3], int out_firstcol, int out_ncol, size_t *offset_send,
                                     size_t *offset_recv, size_t *count_send, size_t *count_recv, size_t just_count_flag)
{
  int j, target, origin, ngrp, recvTask, perm_rev[3], xyz[3], uvw[3];
  size_t nimport, nexport;

  /* determine the inverse permutation */
  for(j = 0; j < 3; j++)
    perm_rev[j] = perm[j];

  if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2)) /* not yet the inverse */
    {
      for(j = 0; j < 3; j++)
        perm_rev[j] = perm[perm[j]];

      if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2))
        Terminate("bummer");
    }

  int in_colums          = Ndims[0] * Ndims[1];
  int in_avg             = (in_colums - 1) / NTask + 1;
  int in_exc             = NTask * in_avg - in_colums;
  int in_tasklastsection = NTask - in_exc;
  int in_pivotcol        = in_tasklastsection * in_avg;

  int out_colums          = Ndims[perm[0]] * Ndims[perm[1]];
  int out_avg             = (out_colums - 1) / NTask + 1;
  int out_exc             = NTask * out_avg - out_colums;
  int out_tasklastsection = NTask - out_exc;
  int out_pivotcol        = out_tasklastsection * out_avg;

  size_t i, ncells = ((size_t)in_ncol) * Ndims[2];

  xyz[0] = in_firstcol / Ndims[1];
  xyz[1] = in_firstcol % Ndims[1];
  xyz[2] = 0;

  memset(count_send, 0, NTask * sizeof(size_t));

  /* loop over all cells in input array and determine target processor */
  for(i = 0; i < ncells; i++)
    {
      /* determine target task */
      uvw[0] = xyz[perm[0]];
      uvw[1] = xyz[perm[1]];
      uvw[2] = xyz[perm[2]];

      int newcol = Ndims[perm[1]] * uvw[0] + uvw[1];
      if(newcol < out_pivotcol)
        target = newcol / out_avg;
      else
        target = (newcol - out_pivotcol) / (out_avg - 1) + out_tasklastsection;

      /* move data element to targettask */

      if(just_count_flag)
        count_send[target]++;
      else
        {
          size_t off  = offset_send[target] + count_send[target]++;
          out[off][0] = data[i][0];
          out[off][1] = data[i][1];
        }
      xyz[2]++;
      if(xyz[2] == Ndims[2])
        {
          xyz[2] = 0;
          xyz[1]++;
          if(xyz[1] == Ndims[1])
            {
              xyz[1] = 0;
              xyz[0]++;
            }
        }
    }

  if(just_count_flag)
    {
      myMPI_Alltoall(count_send, sizeof(size_t), MPI_BYTE, count_recv, sizeof(size_t), MPI_BYTE, Communicator);

      for(j = 0, nimport = 0, nexport = 0, offset_send[0] = 0, offset_recv[0] = 0; j < NTask; j++)
        {
          nexport += count_send[j];
          nimport += count_recv[j];

          if(j > 0)
            {
              offset_send[j] = offset_send[j - 1] + count_send[j - 1];
              offset_recv[j] = offset_recv[j - 1] + count_recv[j - 1];
            }
        }

      if(nexport != ncells)
        Terminate("nexport=%lld != ncells=%lld", (long long)nexport, (long long)ncells);
    }
  else
    {
      nimport = 0;

      /* exchange all the data */
      for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(count_send[recvTask] > 0 || count_recv[recvTask] > 0)
                myMPI_Sendrecv(&out[offset_send[recvTask]], count_send[recvTask] * sizeof(fft_complex), MPI_BYTE, recvTask, TAG_DENS_A,
                               &data[offset_recv[recvTask]], count_recv[recvTask] * sizeof(fft_complex), MPI_BYTE, recvTask,
                               TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);

              nimport += count_recv[recvTask];
            }
        }

      /* now loop over the new cell layout */
      /* find enclosing rectangle around columns in new plane */

      int first[3], last[3];

      first[0] = out_firstcol / Ndims[perm[1]];
      first[1] = out_firstcol % Ndims[perm[1]];
      first[2] = 0;

      last[0] = (out_firstcol + out_ncol - 1) / Ndims[perm[1]];
      last[1] = (out_firstcol + out_ncol - 1) % Ndims[perm[1]];
      last[2] = Ndims[perm[2]] - 1;

      if(first[1] + out_ncol >= Ndims[perm[1]])
        {
          first[1] = 0;
          last[1]  = Ndims[perm[1]] - 1;
        }

      /* now need to map this back to the old coordinates */

      int xyz_first[3], xyz_last[3];

      for(j = 0; j < 3; j++)
        {
          xyz_first[j] = first[perm_rev[j]];
          xyz_last[j]  = last[perm_rev[j]];
        }

      memset(count_recv, 0, NTask * sizeof(size_t));

      size_t count = 0;

      /* traverse an enclosing box around the new cell layout in the old order */
      for(xyz[0] = xyz_first[0]; xyz[0] <= xyz_last[0]; xyz[0]++)
        for(xyz[1] = xyz_first[1]; xyz[1] <= xyz_last[1]; xyz[1]++)
          for(xyz[2] = xyz_first[2]; xyz[2] <= xyz_last[2]; xyz[2]++)
            {
              /* check that the point is actually part of a column */
              uvw[0] = xyz[perm[0]];
              uvw[1] = xyz[perm[1]];
              uvw[2] = xyz[perm[2]];

              int col = uvw[0] * Ndims[perm[1]] + uvw[1];

              if(col >= out_firstcol && col < out_firstcol + out_ncol)
                {
                  /* determine origin task */
                  int newcol = Ndims[1] * xyz[0] + xyz[1];
                  if(newcol < in_pivotcol)
                    origin = newcol / in_avg;
                  else
                    origin = (newcol - in_pivotcol) / (in_avg - 1) + in_tasklastsection;

                  size_t index = ((size_t)Ndims[perm[2]]) * (col - out_firstcol) + uvw[2];

                  /* move data element from origin task */
                  size_t off    = offset_recv[origin] + count_recv[origin]++;
                  out[index][0] = data[off][0];
                  out[index][1] = data[off][1];

                  count++;
                }
            }

      if(count != nimport)
        {
          int fi = out_firstcol % Ndims[perm[1]];
          int la = (out_firstcol + out_ncol - 1) % Ndims[perm[1]];

          Terminate("count=%lld nimport=%lld   ncol=%d fi=%d la=%d first=%d last=%d\n", (long long)count, (long long)nimport, out_ncol,
                    fi, la, first[1], last[1]);
        }
    }
}

void pm_mpi_fft::my_fft_column_transpose(fft_real *data, int Ndims[3], /* global dimensions of data cube */
                                         int in_firstcol, int in_ncol, /* first column and number of columns */
                                         fft_real *out, int perm[3], int out_firstcol, int out_ncol, size_t *count_send,
                                         size_t *count_recv, size_t just_count_flag)
{
  /* determine the inverse permutation */
  int perm_rev[3];
  for(int j = 0; j < 3; j++)
    perm_rev[j] = perm[j];

  if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2)) /* not yet the inverse */
    {
      for(int j = 0; j < 3; j++)
        perm_rev[j] = perm[perm[j]];

      if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2))
        Terminate("bummer");
    }

  int in_colums          = Ndims[0] * Ndims[1];
  int in_avg             = (in_colums - 1) / NTask + 1;
  int in_exc             = NTask * in_avg - in_colums;
  int in_tasklastsection = NTask - in_exc;
  int in_pivotcol        = in_tasklastsection * in_avg;

  int out_colums          = Ndims[perm[0]] * Ndims[perm[1]];
  int out_avg             = (out_colums - 1) / NTask + 1;
  int out_exc             = NTask * out_avg - out_colums;
  int out_tasklastsection = NTask - out_exc;
  int out_pivotcol        = out_tasklastsection * out_avg;

  if(just_count_flag)
    memset(count_send, 0, NTask * sizeof(size_t));

  /* exchange all the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          // check whether we have anything to do
          if(count_send[target] == 0 && count_recv[target] == 0 && just_count_flag == 0)
            continue;

          /* determine enclosing rectangle of current region */
          int source_first[3];
          source_first[0] = in_firstcol / Ndims[1];
          source_first[1] = in_firstcol % Ndims[1];
          source_first[2] = 0;

          int source_last[3];
          source_last[0] = (in_firstcol + in_ncol - 1) / Ndims[1];
          source_last[1] = (in_firstcol + in_ncol - 1) % Ndims[1];
          source_last[2] = Ndims[2] - 1;

          if(source_first[1] + in_ncol >= Ndims[1])
            {
              source_first[1] = 0;
              source_last[1]  = Ndims[1] - 1;
            }

          /* determine target columns */

          int target_first_col    = 0;
          int long target_num_col = 0;

          if(target < out_tasklastsection)
            {
              target_first_col = target * out_avg;
              target_num_col   = out_avg;
            }
          else
            {
              target_first_col = (target - out_tasklastsection) * (out_avg - 1) + out_pivotcol;
              target_num_col   = (out_avg - 1);
            }

          /* find enclosing rectangle around columns in new plane */
          int first[3], last[3];

          first[0] = target_first_col / Ndims[perm[1]];
          first[1] = target_first_col % Ndims[perm[1]];
          first[2] = 0;

          last[0] = (target_first_col + target_num_col - 1) / Ndims[perm[1]];
          last[1] = (target_first_col + target_num_col - 1) % Ndims[perm[1]];
          last[2] = Ndims[perm[2]] - 1;

          if(first[1] + target_num_col >= Ndims[perm[1]])
            {
              first[1] = 0;
              last[1]  = Ndims[perm[1]] - 1;
            }

          /* now we map this back to the old coordinates */
          int xyz_first[3], xyz_last[3];

          for(int j = 0; j < 3; j++)
            {
              xyz_first[j] = first[perm_rev[j]];
              xyz_last[j]  = last[perm_rev[j]];
            }

          /* determine common box */
          int xyz_start[3], xyz_end[3];
          for(int j = 0; j < 3; j++)
            {
              xyz_start[j] = std::max<int>(xyz_first[j], source_first[j]);
              xyz_end[j]   = std::min<int>(xyz_last[j], source_last[j]);
            }

          int xyz[3];
          for(int j = 0; j < 3; j++)
            xyz[j] = xyz_start[j];

          /* now do the same determination for the flipped situation on the target side */

          int flip_in_firstcol = 0;
          int flip_in_ncol     = 0;

          if(target < in_tasklastsection)
            {
              flip_in_firstcol = target * in_avg;
              flip_in_ncol     = in_avg;
            }
          else
            {
              flip_in_firstcol = (target - in_tasklastsection) * (in_avg - 1) + in_pivotcol;
              flip_in_ncol     = (in_avg - 1);
            }

          /* determine enclosing rectangle of current region */
          int flip_source_first[3];
          flip_source_first[0] = flip_in_firstcol / Ndims[1];
          flip_source_first[1] = flip_in_firstcol % Ndims[1];
          flip_source_first[2] = 0;

          int flip_source_last[3];
          flip_source_last[0] = (flip_in_firstcol + flip_in_ncol - 1) / Ndims[1];
          flip_source_last[1] = (flip_in_firstcol + flip_in_ncol - 1) % Ndims[1];
          flip_source_last[2] = Ndims[2] - 1;

          if(flip_source_first[1] + flip_in_ncol >= Ndims[1])
            {
              flip_source_first[1] = 0;
              flip_source_last[1]  = Ndims[1] - 1;
            }

          /* determine target columns */

          int flip_first_col = 0;
          int flip_num_col   = 0;

          if(ThisTask < out_tasklastsection)
            {
              flip_first_col = ThisTask * out_avg;
              flip_num_col   = out_avg;
            }
          else
            {
              flip_first_col = (ThisTask - out_tasklastsection) * (out_avg - 1) + out_pivotcol;
              flip_num_col   = (out_avg - 1);
            }

          /* find enclosing rectangle around columns in new plane */
          int flip_first[3], flip_last[3];

          flip_first[0] = flip_first_col / Ndims[perm[1]];
          flip_first[1] = flip_first_col % Ndims[perm[1]];
          flip_first[2] = 0;

          flip_last[0] = (flip_first_col + flip_num_col - 1) / Ndims[perm[1]];
          flip_last[1] = (flip_first_col + flip_num_col - 1) % Ndims[perm[1]];
          flip_last[2] = Ndims[perm[2]] - 1;

          if(flip_first[1] + flip_num_col >= Ndims[perm[1]])
            {
              flip_first[1] = 0;
              flip_last[1]  = Ndims[perm[1]] - 1;
            }

          /* now we map this back to the old coordinates */
          int abc_first[3], abc_last[3];

          for(int j = 0; j < 3; j++)
            {
              abc_first[j] = flip_first[perm_rev[j]];
              abc_last[j]  = flip_last[perm_rev[j]];
            }

          /* determine common box */
          int abc_start[3], abc_end[3];
          for(int j = 0; j < 3; j++)
            {
              abc_start[j] = std::max<int>(abc_first[j], flip_source_first[j]);
              abc_end[j]   = std::min<int>(abc_last[j], flip_source_last[j]);
            }

          int abc[3];

          for(int j = 0; j < 3; j++)
            abc[j] = abc_start[j];

          size_t tot_count_send = 0;
          size_t tot_count_recv = 0;

          /* now check how much free memory there is on the two partners, use at most half of it */
          size_t parnter_freebytes;
          myMPI_Sendrecv(&Mem.FreeBytes, sizeof(size_t), MPI_BYTE, target, TAG_DENS_B, &parnter_freebytes, sizeof(size_t), MPI_BYTE,
                         target, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);

          size_t freeb = std::min<size_t>(parnter_freebytes, Mem.FreeBytes);

          size_t limit = 0.5 * freeb / (sizeof(fft_real) + sizeof(fft_real));

          if(just_count_flag)
            limit = SIZE_MAX;

          int iter = 0;
          do
            {
              size_t limit_send = count_send[target] - tot_count_send;
              size_t limit_recv = count_recv[target] - tot_count_recv;

              if(just_count_flag)
                {
                  limit_send = SIZE_MAX;
                  limit_recv = SIZE_MAX;
                }
              else
                {
                  if(limit_send > limit)
                    limit_send = limit;

                  if(limit_recv > limit)
                    limit_recv = limit;
                }

              fft_real *buffer_send = NULL;
              fft_real *buffer_recv = NULL;

              if(just_count_flag == 0)
                {
                  buffer_send = (fft_real *)Mem.mymalloc("buffer_send", limit_send * sizeof(fft_real));
                  buffer_recv = (fft_real *)Mem.mymalloc("buffer_recv", limit_recv * sizeof(fft_real));
                }

              /* traverse the common box between the new and old layout  */
              size_t count = 0;

              while(count < limit_send && xyz[0] <= xyz_end[0] && xyz[1] <= xyz_end[1] && xyz[2] <= xyz_end[2])
                {
                  /* check that the point is actually part of a column in the old layout */
                  int col_old = xyz[0] * Ndims[1] + xyz[1];

                  if(col_old >= in_firstcol && col_old < in_firstcol + in_ncol)
                    {
                      /* check that the point is actually part of a column in the new layout */
                      int uvw[3];
                      uvw[0] = xyz[perm[0]];
                      uvw[1] = xyz[perm[1]];
                      uvw[2] = xyz[perm[2]];

                      int col_new = uvw[0] * Ndims[perm[1]] + uvw[1];

                      if(col_new >= target_first_col && col_new < target_first_col + target_num_col)
                        {
                          // ok, we found a match

                          if(just_count_flag)
                            count_send[target]++;
                          else
                            {
                              long long source_cell = (Ndims[1] * xyz[0] + xyz[1] - in_firstcol) * Ndims[2] + xyz[2];

                              buffer_send[count++] = data[source_cell];
                              tot_count_send++;
                            }
                        }
                    }

                  xyz[2]++;
                  if(xyz[2] > xyz_end[2])
                    {
                      xyz[2] = xyz_start[2];
                      xyz[1]++;
                      if(xyz[1] > xyz_end[1])
                        {
                          xyz[1] = xyz_start[1];
                          xyz[0]++;
                        }
                    }
                }

              if(just_count_flag == 0)
                {
                  myMPI_Sendrecv(buffer_send, limit_send * sizeof(fft_real), MPI_BYTE, target, TAG_DENS_A, buffer_recv,
                                 limit_recv * sizeof(fft_real), MPI_BYTE, target, TAG_DENS_A, Communicator, MPI_STATUS_IGNORE);

                  size_t count = 0;
                  while(count < limit_recv && abc[0] <= abc_end[0] && abc[1] <= abc_end[1] && abc[2] <= abc_end[2])
                    {
                      /* check that the point is actually part of a column in the old layout */
                      int col_old = abc[0] * Ndims[1] + abc[1];

                      if(col_old >= flip_in_firstcol && col_old < flip_in_firstcol + flip_in_ncol)
                        {
                          /* check that the point is actually part of a column in the new layout */
                          int uvw[3];
                          uvw[0] = abc[perm[0]];
                          uvw[1] = abc[perm[1]];
                          uvw[2] = abc[perm[2]];

                          int col_new = uvw[0] * Ndims[perm[1]] + uvw[1];

                          if(col_new >= flip_first_col && col_new < flip_first_col + flip_num_col)
                            {
                              // ok, we found a match

                              long long target_cell = (Ndims[perm[1]] * uvw[0] + uvw[1] - flip_first_col) * Ndims[perm[2]] + uvw[2];

                              out[target_cell] = buffer_recv[count++];
                              tot_count_recv++;
                            }
                        }

                      abc[2]++;
                      if(abc[2] > abc_end[2])
                        {
                          abc[2] = abc_start[2];
                          abc[1]++;
                          if(abc[1] > abc_end[1])
                            {
                              abc[1] = abc_start[1];
                              abc[0]++;
                            }
                        }
                    }

                  Mem.myfree(buffer_recv);
                  Mem.myfree(buffer_send);
                }
              else
                break;

              iter++;

              if(iter > 20)
                Terminate("high number of iterations: limit=%lld", (long long)limit);
            }
          while(tot_count_send < count_send[target] || tot_count_recv < count_recv[target]);
        }
    }
  if(just_count_flag)
    myMPI_Alltoall(count_send, sizeof(size_t), MPI_BYTE, count_recv, sizeof(size_t), MPI_BYTE, Communicator);
}

#endif

#endif
