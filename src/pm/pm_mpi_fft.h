/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pm_mpi_fft.h
 *
 *  \brief declaration of a class for carrying out different variants of parallel FFT transforms
 */

#ifndef PM_MPI_FFT_H
#define PM_MPI_FFT_H

#include "gadgetconfig.h"

#include "../mpi_utils/setcomm.h"

#ifndef FFTW
#define CONCAT(prefix, name) prefix##name
#ifdef DOUBLEPRECISION_FFTW
#define FFTW(x) CONCAT(fftw_, x)
#else
#define FFTW(x) CONCAT(fftwf_, x)
#endif
#endif

class pm_mpi_fft : public virtual setcomm
{
 public:
  pm_mpi_fft(MPI_Comm comm) : setcomm(comm) {}

  struct fft_plan
  {
    int NgridX, NgridY, NgridZ;
    int Ngridz, Ngrid2;

    FFTW(plan) forward_plan_zdir;
    FFTW(plan) forward_plan_xdir;
    FFTW(plan) forward_plan_ydir;

    FFTW(plan) backward_plan_zdir;
    FFTW(plan) backward_plan_ydir;
    FFTW(plan) backward_plan_xdir;

#ifndef FFT_COLUMN_BASED

    int *slab_to_task; /*!< Maps a slab index to the task responsible for the slab */
    int *slabs_x_per_task;
    int *first_slab_x_of_task; /*!< Array containing the index of the first slab of each task */
    int *slabs_y_per_task;     /*!< Array containing the number of slabs each task is responsible for */
    int *first_slab_y_of_task; /*!< Array containing the index of the first slab of each task */

    int nslab_x, slabstart_x, nslab_y, slabstart_y;
    int largest_x_slab; /*!< size of the largest slab in x direction */
    int largest_y_slab; /*!< size of the largest slab in y direction */

#else
    size_t max_datasize;
    size_t fftsize;

    int firstcol_XY, ncol_XY, lastcol_XY;
    int firstcol_XZ, ncol_XZ, lastcol_XZ;
    int firstcol_ZY, ncol_ZY, lastcol_ZY;

    int transposed_firstcol, transposed_ncol;
    int second_transposed_firstcol, second_transposed_ncol;
    long long second_transposed_ncells;

    size_t *offsets_send_A;
    size_t *offsets_recv_A;
    size_t *offsets_send_B;
    size_t *offsets_recv_B;
    size_t *offsets_send_C;
    size_t *offsets_recv_C;
    size_t *offsets_send_D;
    size_t *offsets_recv_D;

    size_t *count_send_A;
    size_t *count_recv_A;
    size_t *count_send_B;
    size_t *count_recv_B;
    size_t *count_send_C;
    size_t *count_recv_C;
    size_t *count_send_D;
    size_t *count_recv_D;
    size_t *count_send_13;
    size_t *count_recv_13;
    size_t *count_send_23;
    size_t *count_recv_23;
    size_t *count_send_13back;
    size_t *count_recv_13back;
    size_t *count_send_23back;
    size_t *count_recv_23back;
#endif
  };

  void my_slab_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ);
  void my_slab_based_fft(fft_plan *plan, void *data, void *workspace, int forward);
  void my_slab_based_fft_free(fft_plan *plan);

  void my_column_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ);
  void my_column_based_fft(fft_plan *plan, void *data, void *workspace, int forward);
  void my_column_based_fft_free(fft_plan *plan);

  void my_slab_transposeA(fft_plan *plan, fft_real *field, fft_real *scratch);
  void my_slab_transposeB(fft_plan *plan, fft_real *field, fft_real *scratch);

  void my_fft_swap23(fft_plan *plan, fft_real *data, fft_real *out);
  void my_fft_swap13(fft_plan *plan, fft_real *data, fft_real *out);
  void my_fft_swap23back(fft_plan *plan, fft_real *data, fft_real *out);
  void my_fft_swap13back(fft_plan *plan, fft_real *data, fft_real *out);

 private:
#ifndef FFT_COLUMN_BASED

  void my_slab_transpose(void *av, void *bv, int *sx, int *firstx, int *sy, int *firsty, int nx, int ny, int nz, int mode);

#else
  void my_fft_column_remap(fft_complex *data, int Ndims[3], int in_firstcol, int in_ncol, fft_complex *out, int perm[3],
                           int out_firstcol, int out_ncol, size_t *offset_send, size_t *offset_recv, size_t *count_send,
                           size_t *count_recv, size_t just_count_flag);

  void my_fft_column_transpose(fft_real *data, int Ndims[3], /* global dimensions of data cube */
                               int in_firstcol, int in_ncol, /* first column and number of columns */
                               fft_real *out, int perm[3], int out_firstcol, int out_ncol, size_t *count_send, size_t *count_recv,
                               size_t just_count_flag);

#endif
};

#endif
