/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pm_periodic.cc
 *
 *  \brief routines for periodic PM-force calculation
 */

#include "gadgetconfig.h"

#if defined(PMGRID) && defined(PERIODIC)

#include <fftw3.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <algorithm>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../pm/pm_periodic.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/*!
 * These routines support two different strategies for doing the particle data exchange to assemble the density field
 * and to read out the forces and potentials:
 *
 * The default scheme sends the particle positions to the target slabs, and bins them there. This works usually well for
 * homogenuously loaded boxes, but can be problematic for zoom-in runs. In the latter case,  PM_ZOOM_OPTIMIZED can be
 * activated, where the data is binned on the originating processor followed by assembly of the binned density field.
 *
 * In addition, the routines can be either used with a slab-based FFT (as is traditionally done in FFTW), or with a
 * column-based FFT. The latter requires more communication and is hence usually slower than the slab-based one.
 * But if the number of MPI ranks exceeds the number of cells per dimension, then the column-based one can still scale
 * and offers a balanced memory consumption, whereas this is not the case for the slab-based approach. To select the
 * column-based FFT, the switch FFT_COLUMN_BASED can be activated.
 *
 * The switches PM_ZOOM_OPTIMIZED and FFT_COLUMN_BASED may also be combined, such that there are 4 main modes of how the
 * PM routines may operate.
 *
 * It is also possible to use non-cubical boxes, by means of setting one or several of the LONG_X, LONG_Y, and LONG_Z
 * options in the config file. The values need to be integers, and then BoxSize is stretched by that factor in the
 * corresponding dimension.
 *
 * Finally, one may also use the TreePM routine for simulations where gravity is perdiodic only in two spatial dimensions.
 * The non-periodic dimension is selected via the GRAVITY_TALLBOX flag. Also in this case, arbitrarily stretched boxes can
 * be used, and one can use PM_ZOOM_OPTIMIZED and/or FFT_COLUMN_BASED if desired.
 *
 * If eight times the particle load per processor exceeds 2^31 ~ 2 billion, one should activate NUMPART_PER_TASK_LARGE.
 * The code will check this condition and terminate if this is violated, so there should hopefully be no severe risk
 * to accidentally forget this.
 */

#define GRIDz (GRIDZ / 2 + 1)
#define GRID2 (2 * GRIDz)

/* short-cut macros for accessing different 3D arrays */
#define FI(x, y, z) (((large_array_offset)GRID2) * (GRIDY * (x) + (y)) + (z))
#define FCxy(c, z) (((large_array_offset)GRID2) * ((c)-myplan.firstcol_XY) + (z))
#define FCxz(c, y) (((large_array_offset)GRIDY) * ((c)-myplan.firstcol_XZ) + (y))
#define FCzy(c, x) (((large_array_offset)GRIDX) * ((c)-myplan.firstcol_ZY) + (x))

#ifndef FFT_COLUMN_BASED
#define NI(x, y, z) (((large_array_offset)GRIDZ) * ((y) + (x)*myplan.nslab_y) + (z))
#endif

/*! \brief This routine generates the FFTW-plans to carry out the FFTs later on.
 *
 *  Some auxiliary variables for bookkeeping are also initialized.
 */
void pm_periodic::pm_init_periodic(simparticles *Sp_ptr)
{
  Sp = Sp_ptr;

  Sp->Asmth[0] = ASMTH * All.BoxSize / PMGRID;
  Sp->Rcut[0]  = RCUT * Sp->Asmth[0];

  /* Set up the FFTW-3 plan files. */
  int ndimx[1] = {GRIDX}; /* dimension of the 1D transforms */
  int ndimy[1] = {GRIDY}; /* dimension of the 1D transforms */
  int ndimz[1] = {GRIDZ}; /* dimension of the 1D transforms */

  int max_GRID2 = 2 * (std::max<int>(std::max<int>(GRIDX, GRIDY), GRIDZ) / 2 + 1);

  /* temporarily allocate some arrays to make sure that out-of-place plans are created */
  rhogrid   = (fft_real *)Mem.mymalloc("rhogrid", max_GRID2 * sizeof(fft_real));
  forcegrid = (fft_real *)Mem.mymalloc("forcegrid", max_GRID2 * sizeof(fft_real));

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  int alignflag = FFTW_UNALIGNED;
#endif

  myplan.forward_plan_zdir = FFTW(plan_many_dft_r2c)(1, ndimz, 1, rhogrid, 0, 1, GRID2, (fft_complex *)forcegrid, 0, 1, GRIDz,
                                                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

#ifndef FFT_COLUMN_BASED
  int stride = GRIDz;
#else
  int stride    = 1;
#endif

  myplan.forward_plan_ydir =
      FFTW(plan_many_dft)(1, ndimy, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRIDY, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRIDY, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_xdir =
      FFTW(plan_many_dft)(1, ndimx, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRIDX, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRIDX, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_xdir =
      FFTW(plan_many_dft)(1, ndimx, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRIDX, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRIDX, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_ydir =
      FFTW(plan_many_dft)(1, ndimy, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRIDY, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRIDY, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_zdir = FFTW(plan_many_dft_c2r)(1, ndimz, 1, (fft_complex *)rhogrid, 0, 1, GRIDz, forcegrid, 0, 1, GRID2,
                                                      FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  Mem.myfree(forcegrid);
  Mem.myfree(rhogrid);

#ifndef FFT_COLUMN_BASED

  my_slab_based_fft_init(&myplan, GRIDX, GRIDY, GRIDZ);

  maxfftsize = std::max<int>(myplan.largest_x_slab * GRIDY, myplan.largest_y_slab * GRIDX) * ((size_t)GRID2);

#else

  my_column_based_fft_init(&myplan, GRIDX, GRIDY, GRIDZ);

  maxfftsize = myplan.max_datasize;

#endif

#if defined(GRAVITY_TALLBOX)
  kernel        = (fft_real *)Mem.mymalloc("kernel", maxfftsize * sizeof(fft_real));
  fft_of_kernel = (fft_complex *)kernel;

  pmforce_setup_tallbox_kernel();
#endif
}

/* Below, the two functions
 *
 *           pmforce_ ...... _prepare_density()
 * and
 *           pmforce_ ...... _readout_forces_or_potential(int dim)
 *
 * are defined in two different versions, one that works better for uniform
 * simulations, the other for zoom-in runs. Only one of the two sets is used,
 * depending on the setting of PM_ZOOM_OPTIMIZED.
 */

#ifdef PM_ZOOM_OPTIMIZED

void pm_periodic::pmforce_zoom_optimized_prepare_density(int mode, int *typelist)
{
  int level, recvTask;
  MPI_Status status;

  particle_data *P = Sp->P;

  part = (part_slab_data *)Mem.mymalloc("part", 8 * (NSource * sizeof(part_slab_data)));
  large_numpart_type *part_sortindex =
      (large_numpart_type *)Mem.mymalloc("part_sortindex", 8 * (NSource * sizeof(large_numpart_type)));

#ifdef FFT_COLUMN_BASED
  int columns         = GRIDX * GRIDY;
  int avg             = (columns - 1) / NTask + 1;
  int exc             = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol        = tasklastsection * avg;
#endif

  /* determine the cells each particle accesses */
  for(int idx = 0; idx < NSource; idx++)
    {
      int i = Sp->get_active_index(idx);

      if(P[i].Ti_Current != All.Ti_Current)
        Sp->drift_particle(&P[i], &Sp->SphP[i], All.Ti_Current);

      int slab_x, slab_y, slab_z;
      if(mode == 2)
        {
          slab_x = (P[i].IntPos[0] * POWERSPEC_FOLDFAC) / INTCELL;
          slab_y = (P[i].IntPos[1] * POWERSPEC_FOLDFAC) / INTCELL;
          slab_z = (P[i].IntPos[2] * POWERSPEC_FOLDFAC) / INTCELL;
        }
      else if(mode == 3)
        {
          slab_x = (P[i].IntPos[0] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          slab_y = (P[i].IntPos[1] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          slab_z = (P[i].IntPos[2] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
        }
      else
        {
          slab_x = P[i].IntPos[0] / INTCELL;
          slab_y = P[i].IntPos[1] / INTCELL;
          slab_z = P[i].IntPos[2] / INTCELL;
        }

      large_numpart_type index_on_grid = ((large_numpart_type)idx) << 3;

      for(int xx = 0; xx < 2; xx++)
        for(int yy = 0; yy < 2; yy++)
          for(int zz = 0; zz < 2; zz++)
            {
              int slab_xx = slab_x + xx;
              int slab_yy = slab_y + yy;
              int slab_zz = slab_z + zz;

              if(slab_xx >= GRIDX)
                slab_xx = 0;
              if(slab_yy >= GRIDY)
                slab_yy = 0;
              if(slab_zz >= GRIDZ)
                slab_zz = 0;

              large_array_offset offset = FI(slab_xx, slab_yy, slab_zz);

              part[index_on_grid].partindex   = (i << 3) + (xx << 2) + (yy << 1) + zz;
              part[index_on_grid].globalindex = offset;
              part_sortindex[index_on_grid]   = index_on_grid;
              index_on_grid++;
            }
    }

  /* note: num_on_grid will be  8 times larger than the particle number, but num_field_points will generally be much smaller */

  large_array_offset num_field_points;
  large_numpart_type num_on_grid = ((large_numpart_type)NSource) << 3;

  /* bring the part-field into the order of the accessed cells. This allows the removal of duplicates */

  mycxxsort(part_sortindex, part_sortindex + num_on_grid, pm_periodic_sortindex_comparator(part));

  if(num_on_grid > 0)
    num_field_points = 1;
  else
    num_field_points = 0;

  /* determine the number of unique field points */
  for(large_numpart_type i = 1; i < num_on_grid; i++)
    {
      if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
        num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex = (large_array_offset *)Mem.mymalloc_movable(&localfield_globalindex, "localfield_globalindex",
                                                                      num_field_points * sizeof(large_array_offset));
  localfield_data        = (fft_real *)Mem.mymalloc_movable(&localfield_data, "localfield_data", num_field_points * sizeof(fft_real));
  localfield_first       = (size_t *)Mem.mymalloc_movable(&localfield_first, "localfield_first", NTask * sizeof(size_t));
  localfield_sendcount   = (size_t *)Mem.mymalloc_movable(&localfield_sendcount, "localfield_sendcount", NTask * sizeof(size_t));
  localfield_offset      = (size_t *)Mem.mymalloc_movable(&localfield_offset, "localfield_offset", NTask * sizeof(size_t));
  localfield_recvcount   = (size_t *)Mem.mymalloc_movable(&localfield_recvcount, "localfield_recvcount", NTask * sizeof(size_t));

  for(int i = 0; i < NTask; i++)
    {
      localfield_first[i]     = 0;
      localfield_sendcount[i] = 0;
    }

  /* establish the cross link between the part[ ]-array and the local list of
   * mesh points. Also, count on which CPU the needed field points are stored.
   */
  num_field_points = 0;
  for(large_numpart_type i = 0; i < num_on_grid; i++)
    {
      if(i > 0)
        if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
          num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
        if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
          continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

#ifndef FFT_COLUMN_BASED
      int slab = part[part_sortindex[i]].globalindex / (GRIDY * GRID2);
      int task = myplan.slab_to_task[slab];
#else
      int task, column = part[part_sortindex[i]].globalindex / (GRID2);

      if(column < pivotcol)
        task = column / avg;
      else
        task = (column - pivotcol) / (avg - 1) + tasklastsection;
#endif

      if(localfield_sendcount[task] == 0)
        localfield_first[task] = num_field_points;

      localfield_sendcount[task]++;
    }
  num_field_points++;

  localfield_offset[0] = 0;
  for(int i = 1; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_sendcount[i - 1];

  Mem.myfree_movable(part_sortindex);
  part_sortindex = NULL;

  /* now bin the local particle data onto the mesh list */
  for(large_numpart_type i = 0; i < num_field_points; i++)
    localfield_data[i] = 0;

  for(large_numpart_type i = 0; i < num_on_grid; i += 8)
    {
      int pindex = (part[i].partindex >> 3);
      MyIntPosType rmd_x, rmd_y, rmd_z;

      if(mode == 2)
        {
          rmd_x = (P[pindex].IntPos[0] * POWERSPEC_FOLDFAC) % INTCELL;
          rmd_y = (P[pindex].IntPos[1] * POWERSPEC_FOLDFAC) % INTCELL;
          rmd_z = (P[pindex].IntPos[2] * POWERSPEC_FOLDFAC) % INTCELL;
        }
      else if(mode == 3)
        {
          rmd_x = (P[pindex].IntPos[0] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
          rmd_y = (P[pindex].IntPos[1] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
          rmd_z = (P[pindex].IntPos[2] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
        }
      else
        {
          rmd_x = P[pindex].IntPos[0] % INTCELL;
          rmd_y = P[pindex].IntPos[1] % INTCELL;
          rmd_z = P[pindex].IntPos[2] % INTCELL;
        }

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);
      double dz = rmd_z * (1.0 / INTCELL);

      double weight = P[pindex].getMass();

      if(mode) /* only for power spectrum calculation */
        if(typelist[P[pindex].getType()] == 0)
          continue;

      localfield_data[part[i + 0].localindex] += weight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 1].localindex] += weight * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 2].localindex] += weight * (1.0 - dx) * dy * (1.0 - dz);
      localfield_data[part[i + 3].localindex] += weight * (1.0 - dx) * dy * dz;
      localfield_data[part[i + 4].localindex] += weight * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 5].localindex] += weight * (dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 6].localindex] += weight * (dx)*dy * (1.0 - dz);
      localfield_data[part[i + 7].localindex] += weight * (dx)*dy * dz;
    }

  rhogrid = (fft_real *)Mem.mymalloc_clear("rhogrid", maxfftsize * sizeof(fft_real));

  /* exchange data and add contributions to the local mesh-path */
  myMPI_Alltoall(localfield_sendcount, sizeof(size_t), MPI_BYTE, localfield_recvcount, sizeof(size_t), MPI_BYTE, Communicator);

  for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_data        = (fft_real *)Mem.mymalloc("import_data", localfield_recvcount[recvTask] * sizeof(fft_real));
              import_globalindex = (large_array_offset *)Mem.mymalloc("import_globalindex",
                                                                      localfield_recvcount[recvTask] * sizeof(large_array_offset));

              if(localfield_sendcount[recvTask] > 0 || localfield_recvcount[recvTask] > 0)
                {
                  myMPI_Sendrecv(localfield_data + localfield_offset[recvTask], localfield_sendcount[recvTask] * sizeof(fft_real),
                                 MPI_BYTE, recvTask, TAG_NONPERIOD_A, import_data, localfield_recvcount[recvTask] * sizeof(fft_real),
                                 MPI_BYTE, recvTask, TAG_NONPERIOD_A, Communicator, &status);

                  myMPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_B,
                                 import_globalindex, localfield_recvcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                 TAG_NONPERIOD_B, Communicator, &status);
                }
            }
          else
            {
              import_data        = localfield_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

          /* note: here every element in rhogrid is only accessed once, so there should be no race condition */
          for(size_t i = 0; i < localfield_recvcount[recvTask]; i++)
            {
              /* determine offset in local FFT slab */
#ifndef FFT_COLUMN_BASED
              large_array_offset offset =
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRIDY * ((large_array_offset)GRID2);
#else
              large_array_offset offset = import_globalindex[i] - myplan.firstcol_XY * ((large_array_offset)GRID2);
#endif
              rhogrid[offset] += import_data[i];
            }

          if(level > 0)
            {
              Mem.myfree(import_globalindex);
              Mem.myfree(import_data);
            }
        }
    }
}

/* Function to read out the force component corresponding to spatial dimension 'dim'.
 * If dim is negative, potential values are read out and assigned to particles.
 */
void pm_periodic::pmforce_zoom_optimized_readout_forces_or_potential(fft_real *grid, int dim)
{
  particle_data *P = Sp->P;

#ifdef EVALPOTENTIAL
#ifdef GRAVITY_TALLBOX
  double fac = 1.0 / (((double)GRIDX) * GRIDY * GRIDZ); /* to get potential  */
#else
  double fac = 4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / pow(All.BoxSize, 3); /* to get potential  */
#endif
#endif

  for(int level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      int recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_data        = (fft_real *)Mem.mymalloc("import_data", localfield_recvcount[recvTask] * sizeof(fft_real));
              import_globalindex = (large_array_offset *)Mem.mymalloc("import_globalindex",
                                                                      localfield_recvcount[recvTask] * sizeof(large_array_offset));

              if(localfield_sendcount[recvTask] > 0 || localfield_recvcount[recvTask] > 0)
                {
                  MPI_Status status;
                  myMPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_C,
                                 import_globalindex, localfield_recvcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                 TAG_NONPERIOD_C, Communicator, &status);
                }
            }
          else
            {
              import_data        = localfield_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

          for(size_t i = 0; i < localfield_recvcount[recvTask]; i++)
            {
#ifndef FFT_COLUMN_BASED
              large_array_offset offset =
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRIDY * ((large_array_offset)GRID2);
#else
              large_array_offset offset = import_globalindex[i] - myplan.firstcol_XY * ((large_array_offset)GRID2);
#endif
              import_data[i] = grid[offset];
            }

          if(level > 0)
            {
              MPI_Status status;
              myMPI_Sendrecv(import_data, localfield_recvcount[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A,
                             localfield_data + localfield_offset[recvTask], localfield_sendcount[recvTask] * sizeof(fft_real),
                             MPI_BYTE, recvTask, TAG_NONPERIOD_A, Communicator, &status);

              Mem.myfree(import_globalindex);
              Mem.myfree(import_data);
            }
        }
    }

  /* read out the force/potential values, which all have been assembled in localfield_data */
  for(int idx = 0; idx < NSource; idx++)
    {
      int i = Sp->get_active_index(idx);

#if !defined(HIERARCHICAL_GRAVITY) && defined(TREEPM_NOTIMESPLIT)
      if(!Sp->TimeBinSynchronized[P[i].TimeBinGrav])
        continue;
#endif

      large_numpart_type j = (idx << 3);

      MyIntPosType rmd_x = P[i].IntPos[0] % INTCELL;
      MyIntPosType rmd_y = P[i].IntPos[1] % INTCELL;
      MyIntPosType rmd_z = P[i].IntPos[2] % INTCELL;

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);
      double dz = rmd_z * (1.0 / INTCELL);

      double value = localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                     localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz +
                     localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz) +
                     localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz +
                     localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz) +
                     localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz +
                     localfield_data[part[j + 6].localindex] * (dx)*dy * (1.0 - dz) +
                     localfield_data[part[j + 7].localindex] * (dx)*dy * dz;

      if(dim < 0)
        {
#ifdef EVALPOTENTIAL
#if defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
          P[i].PM_Potential += value * fac;
#else
          P[i].Potential += value * fac;
#endif
#endif
        }
      else
        {
#if defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
          Sp->P[i].GravPM[dim] += value;
#else
          Sp->P[i].GravAccel[dim] += value;
#endif
        }
    }
}

#else

/*
 *  Here come the routines for a different communication algorithm that is better suited for a homogeneously loaded boxes.
 */

void pm_periodic::pmforce_uniform_optimized_prepare_density(int mode, int *typelist)
{
  Sndpm_count = (size_t *)Mem.mymalloc("Sndpm_count", NTask * sizeof(size_t));
  Sndpm_offset = (size_t *)Mem.mymalloc("Sndpm_offset", NTask * sizeof(size_t));
  Rcvpm_count = (size_t *)Mem.mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *)Mem.mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

  particle_data *P = Sp->P;

  /* determine the slabs/columns each particles accesses */

#ifdef FFT_COLUMN_BASED
  int columns = GRIDX * GRIDY;
  int avg = (columns - 1) / NTask + 1;
  int exc = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol = tasklastsection * avg;
#endif

  for(int rep = 0; rep < 2; rep++)
    {
      /* each threads needs to do the loop to clear its send_count[] array */
      for(int j = 0; j < NTask; j++)
        Sndpm_count[j] = 0;

      for(int idx = 0; idx < NSource; idx++)
        {
          int i = Sp->get_active_index(idx);

          if(P[i].Ti_Current != All.Ti_Current)
            Sp->drift_particle(&Sp->P[i], &Sp->SphP[i], All.Ti_Current);

          if(mode) /* only for power spectrum calculation */
            if(typelist[P[i].getType()] == 0)
              continue;

          int slab_x;
          if(mode == 2)
            slab_x = (P[i].IntPos[0] * POWERSPEC_FOLDFAC) / INTCELL;
          else if(mode == 3)
            slab_x = (P[i].IntPos[0] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          else
            slab_x = P[i].IntPos[0] / INTCELL;

          int slab_xx = slab_x + 1;

          if(slab_xx >= GRIDX)
            slab_xx = 0;

#ifndef FFT_COLUMN_BASED
          if(rep == 0)
            {
              int task0 = myplan.slab_to_task[slab_x];
              int task1 = myplan.slab_to_task[slab_xx];

              Sndpm_count[task0]++;

              if(task0 != task1)
                Sndpm_count[task1]++;
            }
          else
            {
              int task0 = myplan.slab_to_task[slab_x];
              int task1 = myplan.slab_to_task[slab_xx];

              size_t ind0 = Sndpm_offset[task0] + Sndpm_count[task0]++;
#ifndef LEAN
              partout[ind0].Mass = P[i].getMass();
#endif
              for(int j = 0; j < 3; j++)
                partout[ind0].IntPos[j] = P[i].IntPos[j];

              if(task0 != task1)
                {
                  size_t ind1 = Sndpm_offset[task1] + Sndpm_count[task1]++;
#ifndef LEAN
                  partout[ind1].Mass = P[i].getMass();
#endif
                  for(int j = 0; j < 3; j++)
                    partout[ind1].IntPos[j] = P[i].IntPos[j];
                }
            }

#else
          int slab_y;
          if(mode == 2)
            slab_y = (P[i].IntPos[1] * POWERSPEC_FOLDFAC) / INTCELL;
          else if(mode == 3)
            slab_y = (P[i].IntPos[1] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          else
            slab_y = P[i].IntPos[1] / INTCELL;

          int slab_yy = slab_y + 1;

          if(slab_yy >= GRIDY)
            slab_yy = 0;

          int column0 = slab_x * GRIDY + slab_y;
          int column1 = slab_x * GRIDY + slab_yy;
          int column2 = slab_xx * GRIDY + slab_y;
          int column3 = slab_xx * GRIDY + slab_yy;

          int task0, task1, task2, task3;

          if(column0 < pivotcol)
            task0 = column0 / avg;
          else
            task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

          if(column1 < pivotcol)
            task1 = column1 / avg;
          else
            task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

          if(column2 < pivotcol)
            task2 = column2 / avg;
          else
            task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

          if(column3 < pivotcol)
            task3 = column3 / avg;
          else
            task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

          if(rep == 0)
            {
              Sndpm_count[task0]++;
              if(task1 != task0)
                Sndpm_count[task1]++;
              if(task2 != task1 && task2 != task0)
                Sndpm_count[task2]++;
              if(task3 != task0 && task3 != task1 && task3 != task2)
                Sndpm_count[task3]++;
            }
          else
            {
              size_t ind0        = Sndpm_offset[task0] + Sndpm_count[task0]++;
#ifndef LEAN
              partout[ind0].Mass = P[i].getMass();
#endif
              for(int j = 0; j < 3; j++)
                partout[ind0].IntPos[j] = P[i].IntPos[j];

              if(task1 != task0)
                {
                  size_t ind1        = Sndpm_offset[task1] + Sndpm_count[task1]++;
#ifndef LEAN
                  partout[ind1].Mass = P[i].getMass();
#endif
                  for(int j = 0; j < 3; j++)
                    partout[ind1].IntPos[j] = P[i].IntPos[j];
                }
              if(task2 != task1 && task2 != task0)
                {
                  size_t ind2        = Sndpm_offset[task2] + Sndpm_count[task2]++;
#ifndef LEAN
                  partout[ind2].Mass = P[i].getMass();
#endif
                  for(int j = 0; j < 3; j++)
                    partout[ind2].IntPos[j] = P[i].IntPos[j];
                }
              if(task3 != task0 && task3 != task1 && task3 != task2)
                {
                  size_t ind3        = Sndpm_offset[task3] + Sndpm_count[task3]++;
#ifndef LEAN
                  partout[ind3].Mass = P[i].getMass();
#endif
                  for(int j = 0; j < 3; j++)
                    partout[ind3].IntPos[j] = P[i].IntPos[j];
                }
            }
#endif
        }

      if(rep == 0)
        {
          myMPI_Alltoall(Sndpm_count, sizeof(size_t), MPI_BYTE, Rcvpm_count, sizeof(size_t), MPI_BYTE, Communicator);

          nimport = 0, nexport = 0, Rcvpm_offset[0] = 0, Sndpm_offset[0] = 0;
          for(int j = 0; j < NTask; j++)
            {
              nexport += Sndpm_count[j];
              nimport += Rcvpm_count[j];

              if(j > 0)
                {
                  Sndpm_offset[j] = Sndpm_offset[j - 1] + Sndpm_count[j - 1];
                  Rcvpm_offset[j] = Rcvpm_offset[j - 1] + Rcvpm_count[j - 1];
                }
            }

          /* allocate import and export buffer */
          partin = (partbuf *)Mem.mymalloc_movable(&partin, "partin", nimport * sizeof(partbuf));
          partout = (partbuf *)Mem.mymalloc("partout", nexport * sizeof(partbuf));
        }
    }

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange particle data */
  myMPI_Alltoallv(partout, Sndpm_count, Sndpm_offset, partin, Rcvpm_count, Rcvpm_offset, sizeof(partbuf), flag_big_all, Communicator);

  Mem.myfree(partout);

  /* allocate cleared density field */
  rhogrid = (fft_real *)Mem.mymalloc_movable_clear(&rhogrid, "rhogrid", maxfftsize * sizeof(fft_real));

#ifndef FFT_COLUMN_BASED
  /* bin particle data onto mesh, in multi-threaded fashion */

  for(size_t i = 0; i < nimport; i++)
    {
      int slab_x, slab_y, slab_z;
      MyIntPosType rmd_x, rmd_y, rmd_z;

      if(mode == 2)
        {
          slab_x = (partin[i].IntPos[0] * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_x = (partin[i].IntPos[0] * POWERSPEC_FOLDFAC) % INTCELL;
          slab_y = (partin[i].IntPos[1] * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_y = (partin[i].IntPos[1] * POWERSPEC_FOLDFAC) % INTCELL;
          slab_z = (partin[i].IntPos[2] * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_z = (partin[i].IntPos[2] * POWERSPEC_FOLDFAC) % INTCELL;
        }
      else if(mode == 3)
        {
          slab_x = (partin[i].IntPos[0] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_x = (partin[i].IntPos[0] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
          slab_y = (partin[i].IntPos[1] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_y = (partin[i].IntPos[1] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
          slab_z = (partin[i].IntPos[2] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_z = (partin[i].IntPos[2] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
        }
      else
        {
          slab_x = partin[i].IntPos[0] / INTCELL;
          rmd_x = partin[i].IntPos[0] % INTCELL;
          slab_y = partin[i].IntPos[1] / INTCELL;
          rmd_y = partin[i].IntPos[1] % INTCELL;
          slab_z = partin[i].IntPos[2] / INTCELL;
          rmd_z = partin[i].IntPos[2] % INTCELL;
        }

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);
      double dz = rmd_z * (1.0 / INTCELL);

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;
      if(slab_yy >= GRIDY)
        slab_yy = 0;
      if(slab_zz >= GRIDZ)
        slab_zz = 0;

#ifdef LEAN
      double mass = All.PartMass;
#else
      double mass = partin[i].Mass;
#endif

      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          rhogrid[FI(slab_x, slab_y, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
          rhogrid[FI(slab_x, slab_y, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));

          rhogrid[FI(slab_x, slab_yy, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
          rhogrid[FI(slab_x, slab_yy, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          rhogrid[FI(slab_xx, slab_y, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
          rhogrid[FI(slab_xx, slab_y, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));

          rhogrid[FI(slab_xx, slab_yy, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
          rhogrid[FI(slab_xx, slab_yy, slab_zz)] += (mass * (dx) * (dy) * (dz));
        }
    }

#else

  int first_col = myplan.firstcol_XY;
  int last_col = myplan.firstcol_XY + myplan.ncol_XY - 1;

  for(size_t i = 0; i < nimport; i++)
    {
      int slab_x, slab_y;
      MyIntPosType rmd_x, rmd_y;
      if(mode == 2)
        {
          slab_x = (partin[i].IntPos[0] * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_x = (partin[i].IntPos[0] * POWERSPEC_FOLDFAC) % INTCELL;
          slab_y = (partin[i].IntPos[1] * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_y = (partin[i].IntPos[1] * POWERSPEC_FOLDFAC) % INTCELL;
        }
      else if(mode == 3)
        {
          slab_x = (partin[i].IntPos[0] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_x = (partin[i].IntPos[0] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
          slab_y = (partin[i].IntPos[1] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_y = (partin[i].IntPos[1] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
        }
      else
        {
          slab_x = partin[i].IntPos[0] / INTCELL;
          rmd_x = partin[i].IntPos[0] % INTCELL;
          slab_y = partin[i].IntPos[1] / INTCELL;
          rmd_y = partin[i].IntPos[1] % INTCELL;
        }

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;

      if(slab_yy >= GRIDY)
        slab_yy = 0;

      int col0 = slab_x * GRIDY + slab_y;
      int col1 = slab_x * GRIDY + slab_yy;
      int col2 = slab_xx * GRIDY + slab_y;
      int col3 = slab_xx * GRIDY + slab_yy;

#ifdef LEAN
      double mass = All.PartMass;
#else
      double mass = partin[i].Mass;
#endif

      int slab_z;
      MyIntPosType rmd_z;
      if(mode == 2)
        {
          slab_z = (partin[i].IntPos[2] * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_z = (partin[i].IntPos[2] * POWERSPEC_FOLDFAC) % INTCELL;
        }
      else if(mode == 3)
        {
          slab_z = (partin[i].IntPos[2] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) / INTCELL;
          rmd_z = (partin[i].IntPos[2] * POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC) % INTCELL;
        }
      else
        {
          slab_z = partin[i].IntPos[2] / INTCELL;
          rmd_z = partin[i].IntPos[2] % INTCELL;
        }

      double dz = rmd_z * (1.0 / INTCELL);

      int slab_zz = slab_z + 1;

      if(slab_zz >= GRIDZ)
        slab_zz = 0;

      if(col0 >= first_col && col0 <= last_col)
        {
          rhogrid[FCxy(col0, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
          rhogrid[FCxy(col0, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
        }

      if(col1 >= first_col && col1 <= last_col)
        {
          rhogrid[FCxy(col1, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
          rhogrid[FCxy(col1, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
        }

      if(col2 >= first_col && col2 <= last_col)
        {
          rhogrid[FCxy(col2, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
          rhogrid[FCxy(col2, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
        }

      if(col3 >= first_col && col3 <= last_col)
        {
          rhogrid[FCxy(col3, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
          rhogrid[FCxy(col3, slab_zz)] += (mass * (dx) * (dy) * (dz));
        }
    }

#endif
}

/* If dim<0, this function reads out the potential, otherwise Cartesian force components.
 */
void pm_periodic::pmforce_uniform_optimized_readout_forces_or_potential_xy(fft_real *grid, int dim)
{
  particle_data *P = Sp->P;

#ifdef EVALPOTENTIAL
#ifdef GRAVITY_TALLBOX
  double fac = 1.0 / (((double)GRIDX) * GRIDY * GRIDZ); /* to get potential  */
#else
  double fac = 4 * M_PI * (LONG_X * LONG_Y * LONG_Z) / pow(All.BoxSize, 3); /* to get potential  */
#endif
#endif

  MyFloat *flistin = (MyFloat *)Mem.mymalloc("flistin", nimport * sizeof(MyFloat));
  MyFloat *flistout = (MyFloat *)Mem.mymalloc("flistout", nexport * sizeof(MyFloat));

#ifdef FFT_COLUMN_BASED
  int columns = GRIDX * GRIDY;
  int avg = (columns - 1) / NTask + 1;
  int exc = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol = tasklastsection * avg;
#endif

  for(size_t i = 0; i < nimport; i++)
    {
      flistin[i] = 0;

      int slab_x = partin[i].IntPos[0] / INTCELL;
      int slab_y = partin[i].IntPos[1] / INTCELL;
      int slab_z = partin[i].IntPos[2] / INTCELL;

      MyIntPosType rmd_x = partin[i].IntPos[0] % INTCELL;
      MyIntPosType rmd_y = partin[i].IntPos[1] % INTCELL;
      MyIntPosType rmd_z = partin[i].IntPos[2] % INTCELL;

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);
      double dz = rmd_z * (1.0 / INTCELL);

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;
      if(slab_yy >= GRIDY)
        slab_yy = 0;
      if(slab_zz >= GRIDZ)
        slab_zz = 0;

#ifndef FFT_COLUMN_BASED
      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_x, slab_y, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_y, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_x, slab_yy, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_yy, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_xx, slab_y, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_y, slab_zz)] * (dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_xx, slab_yy, slab_z)] * (dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_yy, slab_zz)] * (dx) * (dy) * (dz);
        }
#else
      int column0 = slab_x * GRIDY + slab_y;
      int column1 = slab_x * GRIDY + slab_yy;
      int column2 = slab_xx * GRIDY + slab_y;
      int column3 = slab_xx * GRIDY + slab_yy;

      if(column0 >= myplan.firstcol_XY && column0 <= myplan.lastcol_XY)
        {
          flistin[i] += grid[FCxy(column0, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FCxy(column0, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz);
        }
      if(column1 >= myplan.firstcol_XY && column1 <= myplan.lastcol_XY)
        {
          flistin[i] +=
              grid[FCxy(column1, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FCxy(column1, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.firstcol_XY && column2 <= myplan.lastcol_XY)
        {
          flistin[i] +=
              grid[FCxy(column2, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FCxy(column2, slab_zz)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.firstcol_XY && column3 <= myplan.lastcol_XY)
        {
          flistin[i] += grid[FCxy(column3, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FCxy(column3, slab_zz)] * (dx) * (dy) * (dz);
        }
#endif
    }

  /* exchange the potential component data */
  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(MyFloat) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange  data */
  myMPI_Alltoallv(flistin, Rcvpm_count, Rcvpm_offset, flistout, Sndpm_count, Sndpm_offset, sizeof(MyFloat), flag_big_all,
                  Communicator);

  /* each threads needs to do the loop to clear its send_count[] array */
  for(int j = 0; j < NTask; j++)
    Sndpm_count[j] = 0;

  for(int idx = 0; idx < NSource; idx++)
    {
      int i = Sp->get_active_index(idx);

      int slab_x = P[i].IntPos[0] / INTCELL;
      int slab_xx = slab_x + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;

#ifndef FFT_COLUMN_BASED
      int task0 = myplan.slab_to_task[slab_x];
      int task1 = myplan.slab_to_task[slab_xx];

      double value = flistout[Sndpm_offset[task0] + Sndpm_count[task0]++];

      if(task0 != task1)
        value += flistout[Sndpm_offset[task1] + Sndpm_count[task1]++];
#else
      int slab_y = P[i].IntPos[1] / INTCELL;
      int slab_yy = slab_y + 1;

      if(slab_yy >= GRIDY)
        slab_yy = 0;

      int column0 = slab_x * GRIDY + slab_y;
      int column1 = slab_x * GRIDY + slab_yy;
      int column2 = slab_xx * GRIDY + slab_y;
      int column3 = slab_xx * GRIDY + slab_yy;

      int task0, task1, task2, task3;

      if(column0 < pivotcol)
        task0 = column0 / avg;
      else
        task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

      if(column1 < pivotcol)
        task1 = column1 / avg;
      else
        task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

      if(column2 < pivotcol)
        task2 = column2 / avg;
      else
        task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

      if(column3 < pivotcol)
        task3 = column3 / avg;
      else
        task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

      double value = flistout[Sndpm_offset[task0] + Sndpm_count[task0]++];

      if(task1 != task0)
        value += flistout[Sndpm_offset[task1] + Sndpm_count[task1]++];

      if(task2 != task1 && task2 != task0)
        value += flistout[Sndpm_offset[task2] + Sndpm_count[task2]++];

      if(task3 != task0 && task3 != task1 && task3 != task2)
        value += flistout[Sndpm_offset[task3] + Sndpm_count[task3]++];
#endif

#if !defined(HIERARCHICAL_GRAVITY) && defined(TREEPM_NOTIMESPLIT)
      if(!Sp->TimeBinSynchronized[Sp->P[i].TimeBinGrav])
        continue;
#endif

      if(dim < 0)
        {
#ifdef EVALPOTENTIAL
#if defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
          Sp->P[i].PM_Potential += value * fac;
#else
          Sp->P[i].Potential += value * fac;
#endif
#endif
        }
      else
        {
#if defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
          Sp->P[i].GravPM[dim] += value;
#else
          Sp->P[i].GravAccel[dim] += value;
#endif
        }
    }

  Mem.myfree(flistout);
  Mem.myfree(flistin);
}

#ifdef FFT_COLUMN_BASED
void pm_periodic::pmforce_uniform_optimized_readout_forces_or_potential_xz(fft_real *grid, int dim)
{
  if(dim != 1)
    Terminate("bummer");

  size_t *send_count = (size_t *)Mem.mymalloc("send_count", NTask * sizeof(size_t));
  size_t *send_offset = (size_t *)Mem.mymalloc("send_offset", NTask * sizeof(size_t));
  size_t *recv_count = (size_t *)Mem.mymalloc("recv_count", NTask * sizeof(size_t));
  size_t *recv_offset = (size_t *)Mem.mymalloc("recv_offset", NTask * sizeof(size_t));

  struct partbuf
  {
    MyIntPosType IntPos[3];
  };

  partbuf *partin = NULL, *partout = NULL;
  size_t nimport = 0, nexport = 0;

  particle_data *P = Sp->P;

  int columns = GRIDX * GRID2;
  int avg = (columns - 1) / NTask + 1;
  int exc = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol = tasklastsection * avg;

  /* determine the slabs/columns each particles accesses */
  for(int rep = 0; rep < 2; rep++)
    {
      for(int j = 0; j < NTask; j++)
        send_count[j] = 0;

      for(int idx = 0; idx < NSource; idx++)
        {
          int i = Sp->get_active_index(idx);

          if(P[i].Ti_Current != All.Ti_Current)
            Sp->drift_particle(&Sp->P[i], &Sp->SphP[i], All.Ti_Current);

          int slab_x = P[i].IntPos[0] / INTCELL;
          int slab_xx = slab_x + 1;

          if(slab_xx >= GRIDX)
            slab_xx = 0;

          int slab_z = P[i].IntPos[2] / INTCELL;
          int slab_zz = slab_z + 1;

          if(slab_zz >= GRIDZ)
            slab_zz = 0;

          int column0 = slab_x * GRID2 + slab_z;
          int column1 = slab_x * GRID2 + slab_zz;
          int column2 = slab_xx * GRID2 + slab_z;
          int column3 = slab_xx * GRID2 + slab_zz;

          int task0, task1, task2, task3;

          if(column0 < pivotcol)
            task0 = column0 / avg;
          else
            task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

          if(column1 < pivotcol)
            task1 = column1 / avg;
          else
            task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

          if(column2 < pivotcol)
            task2 = column2 / avg;
          else
            task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

          if(column3 < pivotcol)
            task3 = column3 / avg;
          else
            task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

          if(rep == 0)
            {
              send_count[task0]++;
              if(task1 != task0)
                send_count[task1]++;
              if(task2 != task1 && task2 != task0)
                send_count[task2]++;
              if(task3 != task0 && task3 != task1 && task3 != task2)
                send_count[task3]++;
            }
          else
            {
              size_t ind0 = send_offset[task0] + send_count[task0]++;
              for(int j = 0; j < 3; j++)
                partout[ind0].IntPos[j] = P[i].IntPos[j];

              if(task1 != task0)
                {
                  size_t ind1 = send_offset[task1] + send_count[task1]++;
                  for(int j = 0; j < 3; j++)
                    partout[ind1].IntPos[j] = P[i].IntPos[j];
                }
              if(task2 != task1 && task2 != task0)
                {
                  size_t ind2 = send_offset[task2] + send_count[task2]++;

                  for(int j = 0; j < 3; j++)
                    partout[ind2].IntPos[j] = P[i].IntPos[j];
                }
              if(task3 != task0 && task3 != task1 && task3 != task2)
                {
                  size_t ind3 = send_offset[task3] + send_count[task3]++;

                  for(int j = 0; j < 3; j++)
                    partout[ind3].IntPos[j] = P[i].IntPos[j];
                }
            }
        }

      if(rep == 0)
        {
          myMPI_Alltoall(send_count, sizeof(size_t), MPI_BYTE, recv_count, sizeof(size_t), MPI_BYTE, Communicator);

          nimport = 0, nexport = 0;
          recv_offset[0] = send_offset[0] = 0;

          for(int j = 0; j < NTask; j++)
            {
              nexport += send_count[j];
              nimport += recv_count[j];

              if(j > 0)
                {
                  send_offset[j] = send_offset[j - 1] + send_count[j - 1];
                  recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
                }
            }

          /* allocate import and export buffer */
          partin = (partbuf *)Mem.mymalloc("partin", nimport * sizeof(partbuf));
          partout = (partbuf *)Mem.mymalloc("partout", nexport * sizeof(partbuf));
        }
    }

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(send_count[i] * sizeof(partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange particle data */
  myMPI_Alltoallv(partout, send_count, send_offset, partin, recv_count, recv_offset, sizeof(partbuf), flag_big_all, Communicator);

  Mem.myfree(partout);

  MyFloat *flistin = (MyFloat *)Mem.mymalloc("flistin", nimport * sizeof(MyFloat));
  MyFloat *flistout = (MyFloat *)Mem.mymalloc("flistout", nexport * sizeof(MyFloat));

  for(size_t i = 0; i < nimport; i++)
    {
      flistin[i] = 0;

      int slab_x = partin[i].IntPos[0] / INTCELL;
      int slab_y = partin[i].IntPos[1] / INTCELL;
      int slab_z = partin[i].IntPos[2] / INTCELL;

      MyIntPosType rmd_x = partin[i].IntPos[0] % INTCELL;
      MyIntPosType rmd_y = partin[i].IntPos[1] % INTCELL;
      MyIntPosType rmd_z = partin[i].IntPos[2] % INTCELL;

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);
      double dz = rmd_z * (1.0 / INTCELL);

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;
      if(slab_yy >= GRIDY)
        slab_yy = 0;
      if(slab_zz >= GRIDZ)
        slab_zz = 0;

      int column0 = slab_x * GRID2 + slab_z;
      int column1 = slab_x * GRID2 + slab_zz;
      int column2 = slab_xx * GRID2 + slab_z;
      int column3 = slab_xx * GRID2 + slab_zz;

      if(column0 >= myplan.firstcol_XZ && column0 <= myplan.lastcol_XZ)
        {
          flistin[i] += grid[FCxz(column0, slab_y)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FCxz(column0, slab_yy)] * (1.0 - dx) * (dy) * (1.0 - dz);
        }
      if(column1 >= myplan.firstcol_XZ && column1 <= myplan.lastcol_XZ)
        {
          flistin[i] +=
              grid[FCxz(column1, slab_y)] * (1.0 - dx) * (1.0 - dy) * (dz) + grid[FCxz(column1, slab_yy)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.firstcol_XZ && column2 <= myplan.lastcol_XZ)
        {
          flistin[i] +=
              grid[FCxz(column2, slab_y)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FCxz(column2, slab_yy)] * (dx) * (dy) * (1.0 - dz);
        }

      if(column3 >= myplan.firstcol_XZ && column3 <= myplan.lastcol_XZ)
        {
          flistin[i] += grid[FCxz(column3, slab_y)] * (dx) * (1.0 - dy) * (dz) + grid[FCxz(column3, slab_yy)] * (dx) * (dy) * (dz);
        }
    }

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  flag_big = 0;
  for(int i = 0; i < NTask; i++)
    if(send_count[i] * sizeof(MyFloat) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange data */
  myMPI_Alltoallv(flistin, recv_count, recv_offset, flistout, send_count, send_offset, sizeof(MyFloat), flag_big_all, Communicator);

  for(int j = 0; j < NTask; j++)
    send_count[j] = 0;

  /* now assign to original points */
  for(int idx = 0; idx < NSource; idx++)
    {
      int i = Sp->get_active_index(idx);

      int slab_x = P[i].IntPos[0] / INTCELL;
      int slab_xx = slab_x + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;

      int slab_z = P[i].IntPos[2] / INTCELL;
      int slab_zz = slab_z + 1;

      if(slab_zz >= GRIDZ)
        slab_zz = 0;

      int column0 = slab_x * GRID2 + slab_z;
      int column1 = slab_x * GRID2 + slab_zz;
      int column2 = slab_xx * GRID2 + slab_z;
      int column3 = slab_xx * GRID2 + slab_zz;

      int task0, task1, task2, task3;

      if(column0 < pivotcol)
        task0 = column0 / avg;
      else
        task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

      if(column1 < pivotcol)
        task1 = column1 / avg;
      else
        task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

      if(column2 < pivotcol)
        task2 = column2 / avg;
      else
        task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

      if(column3 < pivotcol)
        task3 = column3 / avg;
      else
        task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

      double value = flistout[send_offset[task0] + send_count[task0]++];

      if(task1 != task0)
        value += flistout[send_offset[task1] + send_count[task1]++];

      if(task2 != task1 && task2 != task0)
        value += flistout[send_offset[task2] + send_count[task2]++];

      if(task3 != task0 && task3 != task1 && task3 != task2)
        value += flistout[send_offset[task3] + send_count[task3]++];

#if !defined(HIERARCHICAL_GRAVITY) && defined(TREEPM_NOTIMESPLIT)
      if(!Sp->TimeBinSynchronized[Sp->P[i].TimeBinGrav])
        continue;
#endif

#if defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
      Sp->P[i].GravPM[dim] += value;
#else
      Sp->P[i].GravAccel[dim] += value;
#endif
    }

  Mem.myfree(flistout);
  Mem.myfree(flistin);
  Mem.myfree(partin);
  Mem.myfree(recv_offset);
  Mem.myfree(recv_count);
  Mem.myfree(send_offset);
  Mem.myfree(send_count);
}

void pm_periodic::pmforce_uniform_optimized_readout_forces_or_potential_zy(fft_real *grid, int dim)
{
  if(dim != 0)
    Terminate("bummer");

  size_t *send_count = (size_t *)Mem.mymalloc("send_count", NTask * sizeof(size_t));
  size_t *send_offset = (size_t *)Mem.mymalloc("send_offset", NTask * sizeof(size_t));
  size_t *recv_count = (size_t *)Mem.mymalloc("recv_count", NTask * sizeof(size_t));
  size_t *recv_offset = (size_t *)Mem.mymalloc("recv_offset", NTask * sizeof(size_t));

  struct partbuf
  {
    MyIntPosType IntPos[3];
  };

  partbuf *partin = NULL, *partout = NULL;
  size_t nimport = 0, nexport = 0;

  particle_data *P = Sp->P;

  int columns = GRID2 * GRIDY;
  int avg = (columns - 1) / NTask + 1;
  int exc = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol = tasklastsection * avg;

  /* determine the slabs/columns each particles accesses */
  for(int rep = 0; rep < 2; rep++)
    {
      for(int j = 0; j < NTask; j++)
        send_count[j] = 0;

      for(int idx = 0; idx < NSource; idx++)
        {
          int i = Sp->get_active_index(idx);

          if(P[i].Ti_Current != All.Ti_Current)
            Sp->drift_particle(&Sp->P[i], &Sp->SphP[i], All.Ti_Current);

          int slab_z = P[i].IntPos[2] / INTCELL;
          int slab_zz = slab_z + 1;

          if(slab_zz >= GRIDZ)
            slab_zz = 0;

          int slab_y = P[i].IntPos[1] / INTCELL;
          int slab_yy = slab_y + 1;

          if(slab_yy >= GRIDY)
            slab_yy = 0;

          int column0 = slab_z * GRIDY + slab_y;
          int column1 = slab_z * GRIDY + slab_yy;
          int column2 = slab_zz * GRIDY + slab_y;
          int column3 = slab_zz * GRIDY + slab_yy;

          int task0, task1, task2, task3;

          if(column0 < pivotcol)
            task0 = column0 / avg;
          else
            task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

          if(column1 < pivotcol)
            task1 = column1 / avg;
          else
            task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

          if(column2 < pivotcol)
            task2 = column2 / avg;
          else
            task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

          if(column3 < pivotcol)
            task3 = column3 / avg;
          else
            task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

          if(rep == 0)
            {
              send_count[task0]++;
              if(task1 != task0)
                send_count[task1]++;
              if(task2 != task1 && task2 != task0)
                send_count[task2]++;
              if(task3 != task0 && task3 != task1 && task3 != task2)
                send_count[task3]++;
            }
          else
            {
              size_t ind0 = send_offset[task0] + send_count[task0]++;
              for(int j = 0; j < 3; j++)
                partout[ind0].IntPos[j] = P[i].IntPos[j];

              if(task1 != task0)
                {
                  size_t ind1 = send_offset[task1] + send_count[task1]++;
                  for(int j = 0; j < 3; j++)
                    partout[ind1].IntPos[j] = P[i].IntPos[j];
                }
              if(task2 != task1 && task2 != task0)
                {
                  size_t ind2 = send_offset[task2] + send_count[task2]++;

                  for(int j = 0; j < 3; j++)
                    partout[ind2].IntPos[j] = P[i].IntPos[j];
                }
              if(task3 != task0 && task3 != task1 && task3 != task2)
                {
                  size_t ind3 = send_offset[task3] + send_count[task3]++;

                  for(int j = 0; j < 3; j++)
                    partout[ind3].IntPos[j] = P[i].IntPos[j];
                }
            }
        }

      if(rep == 0)
        {
          myMPI_Alltoall(send_count, sizeof(size_t), MPI_BYTE, recv_count, sizeof(size_t), MPI_BYTE, Communicator);

          nimport = 0, nexport = 0;
          recv_offset[0] = send_offset[0] = 0;

          for(int j = 0; j < NTask; j++)
            {
              nexport += send_count[j];
              nimport += recv_count[j];

              if(j > 0)
                {
                  send_offset[j] = send_offset[j - 1] + send_count[j - 1];
                  recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
                }
            }

          /* allocate import and export buffer */
          partin = (partbuf *)Mem.mymalloc("partin", nimport * sizeof(partbuf));
          partout = (partbuf *)Mem.mymalloc("partout", nexport * sizeof(partbuf));
        }
    }

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(send_count[i] * sizeof(partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange particle data */
  myMPI_Alltoallv(partout, send_count, send_offset, partin, recv_count, recv_offset, sizeof(partbuf), flag_big_all, Communicator);

  Mem.myfree(partout);

  MyFloat *flistin = (MyFloat *)Mem.mymalloc("flistin", nimport * sizeof(MyFloat));
  MyFloat *flistout = (MyFloat *)Mem.mymalloc("flistout", nexport * sizeof(MyFloat));

  for(size_t i = 0; i < nimport; i++)
    {
      flistin[i] = 0;

      int slab_x = partin[i].IntPos[0] / INTCELL;
      int slab_y = partin[i].IntPos[1] / INTCELL;
      int slab_z = partin[i].IntPos[2] / INTCELL;

      MyIntPosType rmd_x = partin[i].IntPos[0] % INTCELL;
      MyIntPosType rmd_y = partin[i].IntPos[1] % INTCELL;
      MyIntPosType rmd_z = partin[i].IntPos[2] % INTCELL;

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);
      double dz = rmd_z * (1.0 / INTCELL);

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;
      if(slab_yy >= GRIDY)
        slab_yy = 0;
      if(slab_zz >= GRIDZ)
        slab_zz = 0;

      int column0 = slab_z * GRIDY + slab_y;
      int column1 = slab_z * GRIDY + slab_yy;
      int column2 = slab_zz * GRIDY + slab_y;
      int column3 = slab_zz * GRIDY + slab_yy;

      if(column0 >= myplan.firstcol_ZY && column0 <= myplan.lastcol_ZY)
        {
          flistin[i] += grid[FCzy(column0, slab_x)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FCzy(column0, slab_xx)] * (dx) * (1.0 - dy) * (1.0 - dz);
        }
      if(column1 >= myplan.firstcol_ZY && column1 <= myplan.lastcol_ZY)
        {
          flistin[i] +=
              grid[FCzy(column1, slab_x)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FCzy(column1, slab_xx)] * (dx) * (dy) * (1.0 - dz);
        }

      if(column2 >= myplan.firstcol_ZY && column2 <= myplan.lastcol_ZY)
        {
          flistin[i] +=
              grid[FCzy(column2, slab_x)] * (1.0 - dx) * (1.0 - dy) * (dz) + grid[FCzy(column2, slab_xx)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.firstcol_ZY && column3 <= myplan.lastcol_ZY)
        {
          flistin[i] += grid[FCzy(column3, slab_x)] * (1.0 - dx) * (dy) * (dz) + grid[FCzy(column3, slab_xx)] * (dx) * (dy) * (dz);
        }
    }

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  flag_big = 0;
  for(int i = 0; i < NTask; i++)
    if(send_count[i] * sizeof(MyFloat) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange data */
  myMPI_Alltoallv(flistin, recv_count, recv_offset, flistout, send_count, send_offset, sizeof(MyFloat), flag_big_all, Communicator);

  for(int j = 0; j < NTask; j++)
    send_count[j] = 0;

  /* now assign to original points */
  for(int idx = 0; idx < NSource; idx++)
    {
      int i = Sp->get_active_index(idx);

      int slab_z = P[i].IntPos[2] / INTCELL;
      int slab_zz = slab_z + 1;

      if(slab_zz >= GRIDZ)
        slab_zz = 0;

      int slab_y = P[i].IntPos[1] / INTCELL;
      int slab_yy = slab_y + 1;

      if(slab_yy >= GRIDY)
        slab_yy = 0;

      int column0 = slab_z * GRIDY + slab_y;
      int column1 = slab_z * GRIDY + slab_yy;
      int column2 = slab_zz * GRIDY + slab_y;
      int column3 = slab_zz * GRIDY + slab_yy;

      int task0, task1, task2, task3;

      if(column0 < pivotcol)
        task0 = column0 / avg;
      else
        task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

      if(column1 < pivotcol)
        task1 = column1 / avg;
      else
        task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

      if(column2 < pivotcol)
        task2 = column2 / avg;
      else
        task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

      if(column3 < pivotcol)
        task3 = column3 / avg;
      else
        task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

      double value = flistout[send_offset[task0] + send_count[task0]++];

      if(task1 != task0)
        value += flistout[send_offset[task1] + send_count[task1]++];

      if(task2 != task1 && task2 != task0)
        value += flistout[send_offset[task2] + send_count[task2]++];

      if(task3 != task0 && task3 != task1 && task3 != task2)
        value += flistout[send_offset[task3] + send_count[task3]++];

#if !defined(HIERARCHICAL_GRAVITY) && defined(TREEPM_NOTIMESPLIT)
      if(!Sp->TimeBinSynchronized[Sp->P[i].TimeBinGrav])
        continue;
#endif

#if defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
      Sp->P[i].GravPM[dim] += value;
#else
      Sp->P[i].GravAccel[dim] += value;
#endif
    }

  Mem.myfree(flistout);
  Mem.myfree(flistin);
  Mem.myfree(partin);
  Mem.myfree(recv_offset);
  Mem.myfree(recv_count);
  Mem.myfree(send_offset);
  Mem.myfree(send_count);
}
#endif

#endif

/*! Calculates the long-range periodic force given the particle positions
 *  using the PM method.  The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potential by fast Fourier transform methods. The potential is finite-differenced
 *  using a 4-point finite differencing formula, and the forces are
 *  interpolated tri-linearly to the particle positions. The CIC kernel is
 *  deconvolved.
 *
 *  For mode=0, normal force calculation, mode=1, only density field construction
 *  for a power spectrum calculation. In the later case, typelist flags the particle
 *  types that should be included in the density field.
 */
void pm_periodic::pmforce_periodic(int mode, int *typelist)
{
  int x, y, z;

  double tstart = Logs.second();

  if(mode == 0)
    mpi_printf("PM-PERIODIC: Starting periodic PM calculation. (Rcut=%g)  presently allocated=%g MB\n", Sp->Rcut[0],
               Mem.getAllocatedBytesInMB());

#ifdef HIERARCHICAL_GRAVITY
  NSource = Sp->TimeBinsGravity.NActiveParticles;
#else
  NSource = Sp->NumPart;
#endif

#ifndef TREEPM_NOTIMESPLIT
  if(NSource != Sp->NumPart)
    Terminate("unexpected NSource != Sp->NumPart");
#endif

#ifndef NUMPART_PER_TASK_LARGE
  if((((long long)Sp->NumPart) << 3) >= (((long long)1) << 31))
    Terminate("We are dealing with a too large particle number per MPI rank - enabling NUMPART_PER_TASK_LARGE might help.");
#endif

  double asmth2 = Sp->Asmth[0] * Sp->Asmth[0];
  double d      = All.BoxSize / PMGRID;
  double dhalf  = 0.5 * d;

#ifdef GRAVITY_TALLBOX
  double fac = 1.0 / (((double)GRIDX) * GRIDY * GRIDZ); /* to get potential  */
#else
  double fac = 4 * M_PI * (LONG_X * LONG_Y * LONG_Z) / pow(All.BoxSize, 3); /* to get potential  */
#endif

  fac *= 1 / (2 * d); /* for finite differencing */

#ifdef PM_ZOOM_OPTIMIZED
  pmforce_zoom_optimized_prepare_density(mode, typelist);
#else
  pmforce_uniform_optimized_prepare_density(mode, typelist);
#endif

  /* note: after density, we still keep the field 'partin' from the density assignment,
   * as we can use this later on to return potential and z-force
   */

  /* allocate the memory to hold the FFT fields */

  forcegrid = (fft_real *)Mem.mymalloc_movable(&forcegrid, "forcegrid", maxfftsize * sizeof(fft_real));

  workspace = forcegrid;

#ifndef FFT_COLUMN_BASED
  fft_of_rhogrid = (fft_complex *)&rhogrid[0];
#else
  fft_of_rhogrid = (fft_complex *)&workspace[0];
#endif

  /* Do the FFT of the density field */
#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, &rhogrid[0], &workspace[0], 1);
#else
  my_column_based_fft(&myplan, rhogrid, workspace, 1); /* result is in workspace, not in rhogrid ! */
#endif

  if(mode != 0)
    {
      pmforce_measure_powerspec(mode - 1, typelist);

#if defined(FFT_COLUMN_BASED) && !defined(PM_ZOOM_OPTIMIZED)
      Mem.myfree_movable(partin);
      partin = NULL;
#endif
    }
  else
    {
      /* multiply with Green's function in order to obtain the potential (or forces for spectral diffencing) */

      double kfacx = 2.0 * M_PI / (GRIDX * d);
      double kfacy = 2.0 * M_PI / (GRIDY * d);
      double kfacz = 2.0 * M_PI / (GRIDZ * d);

#ifdef FFT_COLUMN_BASED
      for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
        {
          large_array_offset ipcell = ip + ((large_array_offset)myplan.second_transposed_firstcol) * GRIDX;
          y                         = ipcell / (GRIDX * GRIDz);
          int yr                    = ipcell % (GRIDX * GRIDz);
          z                         = yr / GRIDX;
          x                         = yr % GRIDX;
#else
      for(x = 0; x < GRIDX; x++)
        for(y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
          for(z = 0; z < GRIDz; z++)
            {
#endif
          int xx, yy, zz;

          if(x >= (GRIDX / 2))
            xx = x - GRIDX;
          else
            xx = x;
          if(y >= (GRIDY / 2))
            yy = y - GRIDY;
          else
            yy = y;
          if(z >= (GRIDZ / 2))
            zz = z - GRIDZ;
          else
            zz = z;

          double kx = kfacx * xx;
          double ky = kfacy * yy;
          double kz = kfacz * zz;

          double k2 = kx * kx + ky * ky + kz * kz;

          double smth = 1.0, deconv = 1.0;

          if(k2 > 0)
            {
              smth = -exp(-k2 * asmth2) / k2;

              /* do deconvolution */

              double fx = 1, fy = 1, fz = 1;

              if(xx != 0)
                {
                  fx = kx * dhalf;
                  fx = sin(fx) / fx;
                }
              if(yy != 0)
                {
                  fy = ky * dhalf;
                  fy = sin(fy) / fy;
                }
              if(zz != 0)
                {
                  fz = kz * dhalf;
                  fz = sin(fz) / fz;
                }

              double ff = 1 / (fx * fy * fz);
              deconv    = ff * ff * ff * ff;

              smth *= deconv; /* deconvolution */
            }

#ifndef FFT_COLUMN_BASED
          large_array_offset ip = ((large_array_offset)GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#endif

#ifdef GRAVITY_TALLBOX
          double re = fft_of_rhogrid[ip][0] * fft_of_kernel[ip][0] - fft_of_rhogrid[ip][1] * fft_of_kernel[ip][1];
          double im = fft_of_rhogrid[ip][0] * fft_of_kernel[ip][1] + fft_of_rhogrid[ip][1] * fft_of_kernel[ip][0];

          fft_of_rhogrid[ip][0] = re * deconv * exp(-k2 * asmth2);
          fft_of_rhogrid[ip][1] = im * deconv * exp(-k2 * asmth2);
#else
              fft_of_rhogrid[ip][0] *= smth;
              fft_of_rhogrid[ip][1] *= smth;
#endif
        }

#ifndef GRAVITY_TALLBOX
#ifdef FFT_COLUMN_BASED
      if(myplan.second_transposed_firstcol == 0)
        fft_of_rhogrid[0][0] = fft_of_rhogrid[0][1] = 0.0;
#else
      if(myplan.slabstart_y == 0)
        fft_of_rhogrid[0][0] = fft_of_rhogrid[0][1] = 0.0;
#endif
#endif

        /* Do the inverse FFT to get the potential/forces */

#ifndef FFT_COLUMN_BASED
      my_slab_based_fft(&myplan, &rhogrid[0], &workspace[0], -1);
#else
      my_column_based_fft(&myplan, workspace, rhogrid, -1);
#endif

      /* Now rhogrid holds the potential/forces */

#ifdef EVALPOTENTIAL
#ifdef PM_ZOOM_OPTIMIZED
      pmforce_zoom_optimized_readout_forces_or_potential(rhogrid, -1);
#else
      pmforce_uniform_optimized_readout_forces_or_potential_xy(rhogrid, -1);
#endif
#endif

      /* get the force components by finite differencing of the potential for each dimension,
       * and send the results back to the right CPUs
       */

      /* we do the x component last, because for differencing the potential in the x-direction, we need to construct the
       * transpose
       */

#ifndef FFT_COLUMN_BASED

      /* z-direction */
      for(y = 0; y < GRIDY; y++)
        for(x = 0; x < myplan.nslab_x; x++)
          for(z = 0; z < GRIDZ; z++)
            {
              int zr = z + 1, zl = z - 1, zrr = z + 2, zll = z - 2;
              if(zr >= GRIDZ)
                zr -= GRIDZ;
              if(zrr >= GRIDZ)
                zrr -= GRIDZ;
              if(zl < 0)
                zl += GRIDZ;
              if(zll < 0)
                zll += GRIDZ;

              forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, y, zl)] - rhogrid[FI(x, y, zr)]) -
                                              (1.0 / 6) * (rhogrid[FI(x, y, zll)] - rhogrid[FI(x, y, zrr)]));
            }

#ifdef PM_ZOOM_OPTIMIZED
      pmforce_zoom_optimized_readout_forces_or_potential(forcegrid, 2);
#else
      pmforce_uniform_optimized_readout_forces_or_potential_xy(forcegrid, 2);
#endif

      /* y-direction */
      for(y = 0; y < GRIDY; y++)
        for(x = 0; x < myplan.nslab_x; x++)
          for(z = 0; z < GRIDZ; z++)
            {
              int yr = y + 1, yl = y - 1, yrr = y + 2, yll = y - 2;
              if(yr >= GRIDY)
                yr -= GRIDY;
              if(yrr >= GRIDY)
                yrr -= GRIDY;
              if(yl < 0)
                yl += GRIDY;
              if(yll < 0)
                yll += GRIDY;

              forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, yl, z)] - rhogrid[FI(x, yr, z)]) -
                                              (1.0 / 6) * (rhogrid[FI(x, yll, z)] - rhogrid[FI(x, yrr, z)]));
            }

#ifdef PM_ZOOM_OPTIMIZED
      pmforce_zoom_optimized_readout_forces_or_potential(forcegrid, 1);
#else
      pmforce_uniform_optimized_readout_forces_or_potential_xy(forcegrid, 1);
#endif

      /* x-direction */
      my_slab_transposeA(&myplan, rhogrid, forcegrid); /* compute the transpose of the potential field for finite differencing */
                                                       /* note: for the x-direction, we difference the transposed field */

      for(x = 0; x < GRIDX; x++)
        for(y = 0; y < myplan.nslab_y; y++)
          for(z = 0; z < GRIDZ; z++)
            {
              int xrr = x + 2, xll = x - 2, xr = x + 1, xl = x - 1;
              if(xr >= GRIDX)
                xr -= GRIDX;
              if(xrr >= GRIDX)
                xrr -= GRIDX;
              if(xl < 0)
                xl += GRIDX;
              if(xll < 0)
                xll += GRIDX;

              forcegrid[NI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[NI(xl, y, z)] - rhogrid[NI(xr, y, z)]) -
                                              (1.0 / 6) * (rhogrid[NI(xll, y, z)] - rhogrid[NI(xrr, y, z)]));
            }

      my_slab_transposeB(&myplan, forcegrid, rhogrid); /* reverse the transpose from above */

#ifdef PM_ZOOM_OPTIMIZED
      pmforce_zoom_optimized_readout_forces_or_potential(forcegrid, 0);
#else
      pmforce_uniform_optimized_readout_forces_or_potential_xy(forcegrid, 0);
#endif

#else

      /* z-direction */
      for(large_array_offset i = 0; i < myplan.ncol_XY; i++)
        {
          fft_real *forcep = &forcegrid[GRID2 * i];
          fft_real *potp = &rhogrid[GRID2 * i];

          for(int z = 0; z < GRIDZ; z++)
            {
              int zr = z + 1;
              int zl = z - 1;
              int zrr = z + 2;
              int zll = z - 2;

              if(zr >= GRIDZ)
                zr -= GRIDZ;
              if(zrr >= GRIDZ)
                zrr -= GRIDZ;
              if(zl < 0)
                zl += GRIDZ;
              if(zll < 0)
                zll += GRIDZ;

              forcep[z] = fac * ((4.0 / 3) * (potp[zl] - potp[zr]) - (1.0 / 6) * (potp[zll] - potp[zrr]));
            }
        }

#ifdef PM_ZOOM_OPTIMIZED
      pmforce_zoom_optimized_readout_forces_or_potential(forcegrid, 2);
#else
      pmforce_uniform_optimized_readout_forces_or_potential_xy(forcegrid, 2);

      /* at this point we can free partin */
      Mem.myfree_movable(partin);
      partin = NULL;
#endif

      /* y-direction */
      my_fft_swap23(&myplan, rhogrid, forcegrid);  // rhogrid contains potential field, forcegrid the transposed field

      /* make an in-place computation */
      fft_real *column = (fft_real *)Mem.mymalloc("column", GRIDY * sizeof(fft_real));

      for(large_array_offset i = 0; i < myplan.ncol_XZ; i++)
        {
          memcpy(column, &forcegrid[GRIDY * i], GRIDY * sizeof(fft_real));

          fft_real *potp = column;
          fft_real *forcep = &forcegrid[GRIDY * i];

          for(int y = 0; y < GRIDY; y++)
            {
              int yr = y + 1;
              int yl = y - 1;
              int yrr = y + 2;
              int yll = y - 2;

              if(yr >= GRIDY)
                yr -= GRIDY;
              if(yrr >= GRIDY)
                yrr -= GRIDY;
              if(yl < 0)
                yl += GRIDY;
              if(yll < 0)
                yll += GRIDY;

              forcep[y] = fac * ((4.0 / 3) * (potp[yl] - potp[yr]) - (1.0 / 6) * (potp[yll] - potp[yrr]));
            }
        }

      Mem.myfree(column);

      /* now need to read out from forcegrid  in a non-standard way */

#ifdef PM_ZOOM_OPTIMIZED
      /* need a third field as scratch space */
      fft_real *scratch = (fft_real *)Mem.mymalloc("scratch", myplan.fftsize * sizeof(fft_real));

      my_fft_swap23back(&myplan, forcegrid, scratch);
      pmforce_zoom_optimized_readout_forces_or_potential(scratch, 1);

      Mem.myfree(scratch);
#else
      pmforce_uniform_optimized_readout_forces_or_potential_xz(forcegrid, 1);
#endif

      /* x-direction */
      my_fft_swap13(&myplan, rhogrid, forcegrid);  // rhogrid contains potential field

      for(large_array_offset i = 0; i < myplan.ncol_ZY; i++)
        {
          fft_real *forcep = &rhogrid[GRIDX * i];
          fft_real *potp = &forcegrid[GRIDX * i];

          for(int x = 0; x < GRIDX; x++)
            {
              int xr = x + 1;
              int xl = x - 1;
              int xrr = x + 2;
              int xll = x - 2;

              if(xr >= GRIDX)
                xr -= GRIDX;
              if(xrr >= GRIDX)
                xrr -= GRIDX;
              if(xl < 0)
                xl += GRIDX;
              if(xll < 0)
                xll += GRIDX;

              forcep[x] = fac * ((4.0 / 3) * (potp[xl] - potp[xr]) - (1.0 / 6) * (potp[xll] - potp[xrr]));
            }
        }

        /* now need to read out from forcegrid in a non-standard way */
#ifdef PM_ZOOM_OPTIMIZED
      my_fft_swap13back(&myplan, rhogrid, forcegrid);
      pmforce_zoom_optimized_readout_forces_or_potential(forcegrid, 0);
#else
      pmforce_uniform_optimized_readout_forces_or_potential_zy(rhogrid, 0);
#endif

#endif
    }

  /* free stuff */

  Mem.myfree(forcegrid);
  Mem.myfree(rhogrid);

#ifdef PM_ZOOM_OPTIMIZED
  Mem.myfree(localfield_recvcount);
  Mem.myfree(localfield_offset);
  Mem.myfree(localfield_sendcount);
  Mem.myfree(localfield_first);
  Mem.myfree(localfield_data);
  Mem.myfree(localfield_globalindex);
  Mem.myfree(part);
#else
#ifndef FFT_COLUMN_BASED
  Mem.myfree(partin);
#endif
  Mem.myfree(Rcvpm_offset);
  Mem.myfree(Rcvpm_count);
  Mem.myfree(Sndpm_offset);
  Mem.myfree(Sndpm_count);
#endif

  double tend = Logs.second();

  if(mode == 0)
    mpi_printf("PM-PERIODIC: done.  (took %g seconds)\n", Logs.timediff(tstart, tend));
}

#ifdef GRAVITY_TALLBOX

/*! This function sets-up the Greens function for calculating the tall-box potential
 *  in real space, with suitable zero padding in the direction of the tall box.
 */
void pm_periodic::pmforce_setup_tallbox_kernel(void)
{
  double d = All.BoxSize / PMGRID;

  mpi_printf("PM-PERIODIC: Setting up tallbox kernel (GRIDX=%d, GRIDY=%d, GRIDZ=%d)\n", GRIDX, GRIDY, GRIDZ);

  /* now set up kernel and its Fourier transform */

  for(int i = 0; i < maxfftsize; i++) /* clear local field */
    kernel[i] = 0;

#ifndef FFT_COLUMN_BASED
  for(int i = myplan.slabstart_x; i < (myplan.slabstart_x + myplan.nslab_x); i++)
    for(int j = 0; j < GRIDY; j++)
      {
#else
  for(int c = myplan.firstcol_XY; c < (myplan.firstcol_XY + myplan.ncol_XY); c++)
    {
      int i = c / GRIDY;
      int j = c % GRIDY;
#endif

        for(int k = 0; k < GRIDZ; k++)
          {
            int ii, jj, kk;

            if(i >= (GRIDX / 2))
              ii = i - GRIDX;
            else
              ii = i;
            if(j >= (GRIDY / 2))
              jj = j - GRIDY;
            else
              jj = j;
            if(k >= (GRIDZ / 2))
              kk = k - GRIDZ;
            else
              kk = k;

            double xx = ii * d;
            double yy = jj * d;
            double zz = kk * d;

            double pot = pmperiodic_tallbox_long_range_potential(xx, yy, zz);

#ifndef FFT_COLUMN_BASED
            size_t ip = FI(i - myplan.slabstart_x, j, k);
#else
          size_t ip = FCxy(c, k);
#endif
            kernel[ip] = pot / All.BoxSize;
          }

#ifndef FFT_COLUMN_BASED
      }
#else
    }
#endif

  /* Do the FFT of the kernel */

  fft_real *workspc = (fft_real *)Mem.mymalloc("workspc", maxfftsize * sizeof(fft_real));

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, kernel, workspc, 1);
#else
  my_column_based_fft(&myplan, kernel, workspc, 1); /* result is in workspace, not in kernel */
  memcpy(kernel, workspc, maxfftsize * sizeof(fft_real));
#endif

  Mem.myfree(workspc);

  mpi_printf("PM-PERIODIC: Done setting up tallbox kernel\n");
}

double pm_periodic::pmperiodic_tallbox_long_range_potential(double x, double y, double z)
{
  x /= All.BoxSize;
  y /= All.BoxSize;
  z /= All.BoxSize;

  double r = sqrt(x * x + y * y + z * z);

  double xx, yy, zz;
  switch(GRAVITY_TALLBOX)
    {
      case 0:
        xx = y;
        yy = z;
        zz = x;
        break;
      case 1:
        xx = x;
        yy = z;
        zz = y;
        break;
      case 2:
        xx = x;
        yy = y;
        zz = z;
        break;
    }
  x = xx;
  y = yy;
  z = zz;

  /* the third dimension, z, is now the non-periodic one */

  double leff  = sqrt(BOXX * BOXY);
  double alpha = 2.0 / leff;

  double sum1 = 0.0;

  int qxmax = (int)(10.0 / (BOXX * alpha) + 0.5);
  int qymax = (int)(10.0 / (BOXY * alpha) + 0.5);

  int nxmax = (int)(4.0 * alpha * BOXX + 0.5);
  int nymax = (int)(4.0 * alpha * BOXY + 0.5);

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double dx = x - nx * BOXX;
            double dy = y - ny * BOXY;
            double r  = sqrt(dx * dx + dy * dy + z * z);
            if(r > 0)
              sum1 += erfc(alpha * r) / r;
          }
        else
          {
            // in the nx/ny=0 case, correct for the short range force
            double alpha_star = 0.5 / (((double)ASMTH) / PMGRID);
            double u          = alpha_star * r;
            if(r > 0)
              sum1 += (erfc(alpha * r) - erfc(u)) / r;
            else
              sum1 += 2.0 / sqrt(M_PI) * (alpha_star - alpha);
          }
      }

  double alpha2 = alpha * alpha;

  double sum2 = 0.0;

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double kx = (2.0 * M_PI / BOXX) * nx;
            double ky = (2.0 * M_PI / BOXY) * ny;
            double k2 = kx * kx + ky * ky;
            double k  = sqrt(k2);

            if(k * z > 0)
              {
                double ex = exp(-k * z);
                if(ex > 0)
                  sum2 += cos(kx * x + ky * y) * (erfc(k / (2 * alpha) + alpha * z) / ex + ex * erfc(k / (2 * alpha) - alpha * z)) / k;
              }
            else
              {
                double ex = exp(k * z);
                if(ex > 0)
                  sum2 += cos(kx * x + ky * y) * (ex * erfc(k / (2 * alpha) + alpha * z) + erfc(k / (2 * alpha) - alpha * z) / ex) / k;
              }
          }
      }

  sum2 *= M_PI / (BOXX * BOXY);

  double psi = 2.0 * alpha / sqrt(M_PI) +
               (2 * sqrt(M_PI) / (BOXX * BOXY) * (exp(-alpha2 * z * z) / alpha + sqrt(M_PI) * z * erf(alpha * z))) - (sum1 + sum2);

  return psi;
}
#endif

/*----------------------------------------------------------------------------------------------------*/
/*           Here comes code for the power-spectrum computation                                       */
/*----------------------------------------------------------------------------------------------------*/

void pm_periodic::calculate_power_spectra(int num)
{
  int n_type[NTYPES];
  long long ntot_type_all[NTYPES];
  /* determine global and local particle numbers */
  for(int n = 0; n < NTYPES; n++)
    n_type[n] = 0;
  for(int n = 0; n < Sp->NumPart; n++)
    n_type[Sp->P[n].getType()]++;

  sumup_large_ints(NTYPES, n_type, ntot_type_all, Communicator);

  int typeflag[NTYPES];

  for(int i = 0; i < NTYPES; i++)
    typeflag[i] = 1;

#ifdef HIERARCHICAL_GRAVITY
  int flag_extra_allocate = 0;
  if(Sp->TimeBinsGravity.ActiveParticleList == NULL)
    {
      flag_extra_allocate = 1;
      Sp->TimeBinsGravity.timebins_allocate();
    }

  Sp->TimeBinsGravity.NActiveParticles = 0;
  for(int i = 0; i < Sp->NumPart; i++)
    Sp->TimeBinsGravity.ActiveParticleList[Sp->TimeBinsGravity.NActiveParticles++] = i;
#endif

  if(ThisTask == 0)
    {
      char buf[MAXLEN_PATH_EXTRA];
      snprintf(buf, MAXLEN_PATH_EXTRA, "%s/powerspecs", All.OutputDir);
      mkdir(buf, 02755);
    }

  snprintf(power_spec_fname, MAXLEN_PATH_EXTRA, "%s/powerspecs/powerspec_%03d.txt", All.OutputDir, num);

  pmforce_do_powerspec(typeflag); /* calculate power spectrum for all particle types */

  /* check whether whether more than one type is in use */
  int count_types = 0;
  for(int i = 0; i < NTYPES; i++)
    if(ntot_type_all[i] > 0)
      count_types++;

  if(count_types > 1)
    for(int i = 0; i < NTYPES; i++)
      {
        if(ntot_type_all[i] > 0)
          {
            for(int j = 0; j < NTYPES; j++)
              typeflag[j] = 0;

            typeflag[i] = 1;

            snprintf(power_spec_fname, MAXLEN_PATH_EXTRA, "%s/powerspecs/powerspec_type%d_%03d.txt", All.OutputDir, i, num);

            pmforce_do_powerspec(typeflag); /* calculate power spectrum for type i */
          }
      }

#ifdef HIERARCHICAL_GRAVITY
  if(flag_extra_allocate)
    Sp->TimeBinsGravity.timebins_free();
#endif
}

void pm_periodic::pmforce_do_powerspec(int *typeflag)
{
  mpi_printf("POWERSPEC: Begin power spectrum. (typeflag=[");
  for(int i = 0; i < NTYPES; i++)
    mpi_printf(" %d ", typeflag[i]);
  mpi_printf("])\n");

  double tstart = Logs.second();

  pmforce_periodic(1, typeflag); /* calculate regular power spectrum for selected particle types */

  pmforce_periodic(2, typeflag); /* calculate folded power spectrum for selected particle types  */

  pmforce_periodic(3, typeflag); /* calculate twice folded power spectrum for selected particle types  */

  double tend = Logs.second();

  mpi_printf("POWERSPEC: End power spectrum. took %g seconds\n", Logs.timediff(tstart, tend));
}

void pm_periodic::pmforce_measure_powerspec(int flag, int *typeflag)
{
  particle_data *P = Sp->P;

  long long CountModes[BINS_PS];
  double SumPowerUncorrected[BINS_PS]; /* without binning correction (as for shot noise) */
  double PowerUncorrected[BINS_PS];    /* without binning correction */
  double DeltaUncorrected[BINS_PS];    /* without binning correction */
  double ShotLimit[BINS_PS];
  double KWeightSum[BINS_PS];
  double Kbin[BINS_PS];

  double mass = 0, mass2 = 0, count = 0;
  for(int i = 0; i < Sp->NumPart; i++)
    if(typeflag[P[i].getType()])
      {
        mass += Sp->P[i].getMass();
        mass2 += Sp->P[i].getMass() * Sp->P[i].getMass();
        count += 1.0;
      }

  MPI_Allreduce(MPI_IN_PLACE, &mass, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(MPI_IN_PLACE, &mass2, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_DOUBLE, MPI_SUM, Communicator);

  double d     = All.BoxSize / PMGRID;
  double dhalf = 0.5 * d;

  double fac = 1.0 / mass;

  double K0     = 2 * M_PI / All.BoxSize;                                                        /* minimum k */
  double K1     = 2 * M_PI / All.BoxSize * (POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC * PMGRID / 2); /* maximum k that can be measured */
  double binfac = BINS_PS / (log(K1) - log(K0));

  double kfacx = 2.0 * M_PI * LONG_X / All.BoxSize;
  double kfacy = 2.0 * M_PI * LONG_Y / All.BoxSize;
  double kfacz = 2.0 * M_PI * LONG_Z / All.BoxSize;

  for(int i = 0; i < BINS_PS; i++)
    {
      SumPowerUncorrected[i] = 0;
      CountModes[i]          = 0;
      KWeightSum[i]          = 0;
    }

#ifdef FFT_COLUMN_BASED
  for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
    {
      large_array_offset ipcell = ip + ((large_array_offset)myplan.second_transposed_firstcol) * GRIDX;
      int y                     = ipcell / (GRIDX * GRIDz);
      int yr                    = ipcell % (GRIDX * GRIDz);
      int z                     = yr / GRIDX;
      int x                     = yr % GRIDX;
#else
  for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
    for(int x = 0; x < GRIDX; x++)
      for(int z = 0; z < GRIDz; z++)
        {
#endif
      int count_double;

      if(z >= 1 &&
         z < (GRIDZ + 1) / 2) /* these modes need to be counted twice due to the storage scheme for the FFT of a real field */
        count_double = 1;
      else
        count_double = 0;

      int xx, yy, zz;

      if(x >= (GRIDX / 2))
        xx = x - GRIDX;
      else
        xx = x;

      if(y >= (GRIDY / 2))
        yy = y - GRIDY;
      else
        yy = y;

      if(z >= (GRIDZ / 2))
        zz = z - GRIDZ;
      else
        zz = z;

      double kx = kfacx * xx;
      double ky = kfacy * yy;
      double kz = kfacz * zz;

      double k2 = kx * kx + ky * ky + kz * kz;

      if(k2 > 0)
        {
          /* do deconvolution */
          double fx = 1, fy = 1, fz = 1;

          if(xx != 0)
            {
              fx = kx * dhalf;
              fx = sin(fx) / fx;
            }
          if(yy != 0)
            {
              fy = ky * dhalf;
              fy = sin(fy) / fy;
            }
          if(zz != 0)
            {
              fz = kz * dhalf;
              fz = sin(fz) / fz;
            }
          double ff   = 1 / (fx * fy * fz);
          double smth = ff * ff * ff * ff;
          /*
           * Note: The Fourier-transform of the density field (rho_hat) must be multiplied with ff^2
           * in order to do the de-convolution. Thats why po = rho_hat^2 gains a factor of ff^4.
           */
          /* end deconvolution */

#ifndef FFT_COLUMN_BASED
          large_array_offset ip = ((large_array_offset)GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#endif

          double po = (fft_of_rhogrid[ip][0] * fft_of_rhogrid[ip][0] + fft_of_rhogrid[ip][1] * fft_of_rhogrid[ip][1]);

          po *= fac * fac * smth;

          double k = sqrt(k2);

          if(flag == 1)
            k *= POWERSPEC_FOLDFAC;
          else if(flag == 2)
            k *= POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC;

          if(k >= K0 && k < K1)
            {
              int bin = log(k / K0) * binfac;

              SumPowerUncorrected[bin] += po;
              CountModes[bin] += 1;
              KWeightSum[bin] += log(k);

              if(count_double)
                {
                  SumPowerUncorrected[bin] += po;
                  CountModes[bin] += 1;
                  KWeightSum[bin] += log(k);
                }
            }
        }
    }

  MPI_Allreduce(MPI_IN_PLACE, SumPowerUncorrected, BINS_PS, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(MPI_IN_PLACE, CountModes, BINS_PS, MPI_LONG_LONG, MPI_SUM, Communicator);
  MPI_Allreduce(MPI_IN_PLACE, KWeightSum, BINS_PS, MPI_DOUBLE, MPI_SUM, Communicator);

  int count_non_zero_bins = 0;
  for(int i = 0; i < BINS_PS; i++)
    {
      if(CountModes[i] > 0)
        {
          Kbin[i] = exp(KWeightSum[i] / CountModes[i]);
          count_non_zero_bins++;
        }
      else
        Kbin[i] = exp((i + 0.5) / binfac + log(K0));

      if(CountModes[i] > 0)
        PowerUncorrected[i] = SumPowerUncorrected[i] / CountModes[i];
      else
        PowerUncorrected[i] = 0;

      DeltaUncorrected[i] = 4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3) * PowerUncorrected[i];

      ShotLimit[i] = 4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3) * (mass2 / (mass * mass));
    }

  /* store the result */
  if(ThisTask == 0)
    {
      FILE *fd;

      if(flag == 0)
        {
          if(!(fd = fopen(power_spec_fname, "w"))) /* store the unfolded spectrum */
            Terminate("can't open file `%s`\n", power_spec_fname);
        }
      else if(flag == 1 || flag == 2)
        {
          if(!(fd = fopen(power_spec_fname, "a"))) /* append the file, store the folded spectrum */
            Terminate("can't open file `%s`\n", power_spec_fname);
        }
      else
        Terminate("Something wrong.\n");

      fprintf(fd, "%g\n", All.Time);
      fprintf(fd, "%d\n", count_non_zero_bins);
      fprintf(fd, "%g\n", All.BoxSize);
      fprintf(fd, "%d\n", (int)(PMGRID));
      if(All.ComovingIntegrationOn)
        fprintf(fd, "%g\n", All.ComovingIntegrationOn > 0 ? linear_growth_factor(All.Time, 1.0) : 1.0);

      for(int i = 0; i < BINS_PS; i++)
        if(CountModes[i] > 0)
          fprintf(fd, "%g %g %g %g %g\n", Kbin[i], DeltaUncorrected[i], PowerUncorrected[i], (double)CountModes[i], ShotLimit[i]);

      if(flag == 2)
        {
          fprintf(fd, "%g\n", mass);
          fprintf(fd, "%g\n", count);
          fprintf(fd, "%g\n", mass * mass / mass2);
        }

      fclose(fd);
    }
}

#endif
