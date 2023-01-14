/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pm_nonperiodic.cc
 *
 *  \brief code for non-periodic long-range PM force calculation
 */

#include "gadgetconfig.h"

#if defined(PMGRID) && (!defined(PERIODIC) || defined(PLACEHIGHRESREGION))

#include <fftw3.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../pm/pm_mpi_fft.h"
#include "../pm/pm_nonperiodic.h"
#include "../sort/cxxsort.h"
#include "../src/gravtree/gravtree.h"
#include "../src/time_integration/timestep.h"
#include "../system/system.h"

#define GRID (HRPMGRID)
#define GRIDz (GRID / 2 + 1)
#define GRID2 (2 * GRIDz)

#define FI(x, y, z) (((large_array_offset)GRID2) * (GRID * (x) + (y)) + (z))
#define FC(c, z) (((large_array_offset)GRID2) * ((c)-myplan.firstcol_XY) + (z))
#define TI(x, y, z) (((large_array_offset)GRID) * ((x) + (y)*myplan.nslab_x) + (z))

/*! This function determines the particle extension of all particles, and for
 *  those types selected with PLACEHIGHRESREGION if this is used, and then
 *  determines the boundaries of the non-periodic FFT-mesh that can be placed
 *  on this region. Note that a sufficient buffer region at the rim of the
 *  occupied part of the mesh needs to be reserved in order to allow a correct
 *  finite differencing using a 4-point formula. In addition, to allow
 *  non-periodic boundaries, the actual FFT mesh used is twice as large in
 *  each dimension compared with GRID.
 */
void pm_nonperiodic::pm_init_regionsize(void)
{
  /* first, find a reference coordinate by selecting an arbitrary particle in the respective regions. For definiteness, we choose the
   * first particle */

  particle_data *P  = Sp->P;
  int have_low_mesh = NTask, have_high_mesh = NTask; /* default is we don't have a particle */

  if(Sp->NumPart > 0)
    {
      for(int j = 0; j < 3; j++)
        Sp->ReferenceIntPos[LOW_MESH][j] = P[0].IntPos[j];

      have_low_mesh = ThisTask;
    }

  for(int i = 0; i < Sp->NumPart; i++)
    {
#ifdef PLACEHIGHRESREGION
      if(((1 << P[i].getType()) & (PLACEHIGHRESREGION)))
        {
          for(int j = 0; j < 3; j++)
            Sp->ReferenceIntPos[HIGH_MESH][j] = P[i].IntPos[j];

          have_high_mesh = ThisTask;
          break;
        }
#endif
    }

  int have_global[4] = {have_low_mesh, ThisTask, have_high_mesh, ThisTask};

  MPI_Allreduce(MPI_IN_PLACE, have_global, 2, MPI_2INT, MPI_MINLOC, Communicator);

  if(have_global[0] >= NTask)
    Terminate("have_global[0] >= NTask: Don't we have any particle?");

  MPI_Bcast(&Sp->ReferenceIntPos[LOW_MESH][0], 3 * sizeof(MyIntPosType), MPI_BYTE, have_global[1], Communicator);

#ifdef PLACEHIGHRESREGION
  if(have_global[2] >= NTask)
    Terminate("have_global[2] >= NTask: Apparently there are no particles in high res region");

  MPI_Bcast(&Sp->ReferenceIntPos[HIGH_MESH][0], 3 * sizeof(MyIntPosType), MPI_BYTE, have_global[3], Communicator);
#endif

  /* find enclosing rectangle */

  MySignedIntPosType xmin[2][3], xmax[2][3];

  for(int j = 0; j < 3; j++)
    {
      xmin[LOW_MESH][j] = xmin[HIGH_MESH][j] = 0;
      xmax[LOW_MESH][j] = xmax[HIGH_MESH][j] = 0;
    }

  for(int i = 0; i < Sp->NumPart; i++)
    {
      MyIntPosType diff[3] = {P[i].IntPos[0] - Sp->ReferenceIntPos[LOW_MESH][0], P[i].IntPos[1] - Sp->ReferenceIntPos[LOW_MESH][1],
                              P[i].IntPos[2] - Sp->ReferenceIntPos[LOW_MESH][2]};

      MySignedIntPosType *delta = (MySignedIntPosType *)diff;

      for(int j = 0; j < 3; j++)
        {
          if(delta[j] > xmax[LOW_MESH][j])
            xmax[LOW_MESH][j] = delta[j];
          if(delta[j] < xmin[LOW_MESH][j])
            xmin[LOW_MESH][j] = delta[j];
        }

#ifdef PLACEHIGHRESREGION
      if(((1 << P[i].getType()) & (PLACEHIGHRESREGION)))
        {
          MyIntPosType diff[3] = {P[i].IntPos[0] - Sp->ReferenceIntPos[HIGH_MESH][0],
                                  P[i].IntPos[1] - Sp->ReferenceIntPos[HIGH_MESH][1],
                                  P[i].IntPos[2] - Sp->ReferenceIntPos[HIGH_MESH][2]};

          MySignedIntPosType *delta = (MySignedIntPosType *)diff;

          for(int j = 0; j < 3; j++)
            {
              if(delta[j] > xmax[HIGH_MESH][j])
                xmax[HIGH_MESH][j] = delta[j];
              if(delta[j] < xmin[HIGH_MESH][j])
                xmin[HIGH_MESH][j] = delta[j];
            }
        }
#endif
    }

  MPI_Allreduce(xmin, Sp->Xmintot, 6, MPI_MyIntPosType, MPI_MIN_MySignedIntPosType, Communicator);
  MPI_Allreduce(xmax, Sp->Xmaxtot, 6, MPI_MyIntPosType, MPI_MAX_MySignedIntPosType, Communicator);

  for(int i = 0; i < 3; i++)
    {
      Sp->Xmaxtot[LOW_MESH][i] += 1; /* so that all particles fulfill   xmin <= pos < xmax instead of xmin <= pos <= xmax*/
      Sp->Xmaxtot[HIGH_MESH][i] += 1;
    }

  MyIntPosType inner_meshsize[2], enclosing_meshsize[2];
  int flag_recompute_kernel = 0;

  for(int mesh = 0; mesh < 2; mesh++)
    {
      inner_meshsize[mesh] = (MyIntPosType)(Sp->Xmaxtot[mesh][0] - Sp->Xmintot[mesh][0]);

      if((MyIntPosType)(Sp->Xmaxtot[mesh][1] - Sp->Xmintot[mesh][1]) > inner_meshsize[mesh])
        inner_meshsize[mesh] = Sp->Xmaxtot[mesh][1] - Sp->Xmintot[mesh][1];

      if((MyIntPosType)(Sp->Xmaxtot[mesh][2] - Sp->Xmintot[mesh][2]) > inner_meshsize[mesh])
        inner_meshsize[mesh] = Sp->Xmaxtot[mesh][2] - Sp->Xmintot[mesh][2];
    }

  for(int mesh = 0; mesh < 2; mesh++)
    {
#ifdef PERIODIC
      if(mesh == LOW_MESH)
        continue;
#endif
#ifndef PLACEHIGHRESREGION
      if(mesh == HIGH_MESH)
        continue;
#endif

      MyIntPosType blocksize = 1;
      MyIntPosType mask      = ~((MyIntPosType)0);

      if(mesh == LOW_MESH)
        {
          MyIntPosType refsize = inner_meshsize[mesh] >> 4; /* pick 1/8 as reference size */

          for(int i = 0; i < BITS_FOR_POSITIONS; i++)
            {
              if(blocksize >= refsize)
                break;

              blocksize <<= 1;
              mask <<= 1;
            }
        }
      else
        {
#if defined(RANDOMIZE_DOMAINCENTER_TYPES) || defined(PLACEHIGHRESREGION)
          blocksize = Sp->PlacingBlocksize;
#else
          Terminate("we should not get here");
#endif
        }

      mpi_printf(
          "PM-NONPERIODIC: Allowed region for isolated PM mesh (%s):  BEFORE  (%g|%g|%g) -> (%g|%g|%g)  inner meshsize=%g  "
          "blocksize=%g\n",
          mesh == LOW_MESH ? "coarse" : "fine", Sp->FacIntToCoord * Sp->Xmintot[mesh][0], Sp->FacIntToCoord * Sp->Xmintot[mesh][1],
          Sp->FacIntToCoord * Sp->Xmintot[mesh][2], Sp->FacIntToCoord * Sp->Xmaxtot[mesh][0], Sp->FacIntToCoord * Sp->Xmaxtot[mesh][1],
          Sp->FacIntToCoord * Sp->Xmaxtot[mesh][2], Sp->FacIntToCoord * inner_meshsize[mesh], Sp->FacIntToCoord * blocksize);

      enclosing_meshsize[mesh] = 0;

      /* expand the box so that it aligns with blocksize */
      for(int i = 0; i < 3; i++)
        {
          MyIntPosType left, right;
          if(mesh == LOW_MESH)
            {
              /* now round it down */
              left = ((Sp->Xmintot[mesh][i] + Sp->Xmaxtot[mesh][i]) / 2 - inner_meshsize[mesh] / 2) + Sp->ReferenceIntPos[mesh][i];
              left &= mask;

              /* now round it up */
              right = ((Sp->Xmintot[mesh][i] + Sp->Xmaxtot[mesh][i]) / 2 + inner_meshsize[mesh] / 2) + Sp->ReferenceIntPos[mesh][i];
              right &= mask;
              right += blocksize;
            }
          else
            {
#if defined(RANDOMIZE_DOMAINCENTER_TYPES) || defined(PLACEHIGHRESREGION)
              left  = (Sp->ReferenceIntPos[HIGH_MESH][i] + Sp->Xmintot[HIGH_MESH][i]) & Sp->PlacingMask;
              right = left + Sp->PlacingBlocksize;
#else
              Terminate("we should not get here");
#endif
            }

          Sp->Xmintot[mesh][i] = left - Sp->ReferenceIntPos[mesh][i];
          Sp->Xmaxtot[mesh][i] = right - Sp->ReferenceIntPos[mesh][i];

          Sp->Left[mesh][i] = left;

          Sp->MeshSize[mesh][i] = Sp->Xmaxtot[mesh][i] - Sp->Xmintot[mesh][i];

          if(Sp->MeshSize[mesh][i] > enclosing_meshsize[mesh])
            enclosing_meshsize[mesh] = Sp->MeshSize[mesh][i];
        }

      mpi_printf(
          "PM-NONPERIODIC: Allowed region for isolated PM mesh (%s):  AFTER   (%g|%g|%g) -> (%g|%g|%g)  enclosing_meshsize=%g   "
          "absolute space (%g|%g|%g) -> (%g|%g|%g)\n",
          mesh == LOW_MESH ? "coarse" : "fine", Sp->FacIntToCoord * Sp->Xmintot[mesh][0], Sp->FacIntToCoord * Sp->Xmintot[mesh][1],
          Sp->FacIntToCoord * Sp->Xmintot[mesh][2], Sp->FacIntToCoord * Sp->Xmaxtot[mesh][0], Sp->FacIntToCoord * Sp->Xmaxtot[mesh][1],
          Sp->FacIntToCoord * Sp->Xmaxtot[mesh][2], Sp->FacIntToCoord * enclosing_meshsize[mesh],
          Sp->FacIntToCoord * (Sp->Xmintot[mesh][0] + Sp->ReferenceIntPos[mesh][0]),
          Sp->FacIntToCoord * (Sp->Xmintot[mesh][1] + Sp->ReferenceIntPos[mesh][1]),
          Sp->FacIntToCoord * (Sp->Xmintot[mesh][2] + Sp->ReferenceIntPos[mesh][2]),
          Sp->FacIntToCoord * (Sp->Xmaxtot[mesh][0] + Sp->ReferenceIntPos[mesh][0]),
          Sp->FacIntToCoord * (Sp->Xmaxtot[mesh][1] + Sp->ReferenceIntPos[mesh][1]),
          Sp->FacIntToCoord * (Sp->Xmaxtot[mesh][2] + Sp->ReferenceIntPos[mesh][2]));

      if(enclosing_meshsize[mesh] != Sp->OldMeshSize[mesh])
        {
          flag_recompute_kernel = 1;
          Sp->OldMeshSize[mesh] = enclosing_meshsize[mesh];
        }

      /* this will produce enough room for zero-padding and buffer region to
       allow finite differencing of the potential  */
      Sp->TotalMeshSize[mesh] = 2.0 * enclosing_meshsize[mesh] * (GRID / ((double)(GRID - 10)));

      /* move lower left corner by two cells to allow finite differencing of the potential by a 4-point function */
      for(int i = 0; i < 3; i++)
        Sp->Corner[mesh][i] = Sp->Xmintot[mesh][i] - (2.0 * Sp->TotalMeshSize[mesh] / GRID);

      Sp->Asmth[mesh] = ASMTH * (Sp->FacIntToCoord * Sp->TotalMeshSize[mesh]) / GRID;
      Sp->Rcut[mesh]  = RCUT * Sp->Asmth[mesh];
    }

  static int first_init_done = 0; /* for detecting restart from restartfiles */
  if(flag_recompute_kernel || first_init_done == 0)
    {
#ifndef PERIODIC
      mpi_printf("PM-NONPERIODIC: Recompute kernel course mesh:  Asmth=%g   Rcut=%g     mesh cell size=%g\n", Sp->Asmth[LOW_MESH],
                 Sp->Rcut[LOW_MESH], Sp->FacIntToCoord * Sp->TotalMeshSize[LOW_MESH] / GRID);
#endif
#ifdef PLACEHIGHRESREGION
      mpi_printf("PM-NONPERIODIC: Recompute kernel fine mesh:    Asmth=%g   Rcut=%g     mesh cell size=%g\n", Sp->Asmth[HIGH_MESH],
                 Sp->Rcut[HIGH_MESH], Sp->FacIntToCoord * Sp->TotalMeshSize[HIGH_MESH] / GRID);
#endif

      pm_setup_nonperiodic_kernel();
      first_init_done = 1;
    }
}

/*! Initialization of the non-periodic PM routines. The plan-files for FFTW
 *  are created. Finally, the routine to set-up the non-periodic Greens
 *  function is called.
 */
void pm_nonperiodic::pm_init_nonperiodic(simparticles *Sp_ptr)
{
  Sp = Sp_ptr;

  /* Set up the FFTW-3 plan files. */
  int ndim[1] = {GRID}; /* dimension of the 1D transforms */

  /* temporarily allocate some arrays to make sure that out-of-place plans are created */
  rhogrid   = (fft_real *)Mem.mymalloc("rhogrid", GRID2 * sizeof(fft_real));
  forcegrid = (fft_real *)Mem.mymalloc("forcegrid", GRID2 * sizeof(fft_real));

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  int alignflag = FFTW_UNALIGNED;
#endif
#ifndef FFT_COLUMN_BASED
  int stride = GRIDz;
#else
  int stride = 1;
#endif

  myplan.forward_plan_zdir = FFTW(plan_many_dft_r2c)(1, ndim, 1, rhogrid, 0, 1, GRID2, (fft_complex *)forcegrid, 0, 1, GRIDz,
                                                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_xdir =
      FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRID, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_ydir =
      FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRID, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_zdir = FFTW(plan_many_dft_c2r)(1, ndim, 1, (fft_complex *)rhogrid, 0, 1, GRIDz, forcegrid, 0, 1, GRID2,
                                                      FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_xdir =
      FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRID, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRID, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_ydir =
      FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRID, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRID, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  Mem.myfree(forcegrid);
  Mem.myfree(rhogrid);

#ifndef FFT_COLUMN_BASED

  my_slab_based_fft_init(&myplan, GRID, GRID, GRID);

  maxfftsize = myplan.largest_x_slab * GRID * ((size_t)GRID2);

#else

  my_column_based_fft_init(&myplan, GRID, GRID, GRID);

  maxfftsize = myplan.max_datasize;

#endif

  /* now allocate memory to hold the FFT fields */

  size_t bytes, bytes_tot = 0;

#if !defined(PERIODIC)
  kernel[0] = (fft_real *)Mem.mymalloc("kernel[0]", bytes = maxfftsize * sizeof(fft_real));
  bytes_tot += bytes;
  fft_of_kernel[0] = (fft_complex *)kernel[0];
#endif

#if defined(PLACEHIGHRESREGION)
  kernel[1] = (fft_real *)Mem.mymalloc("kernel[1]", bytes = maxfftsize * sizeof(fft_real));
  bytes_tot += bytes;
  fft_of_kernel[1] = (fft_complex *)kernel[1];
#endif

  mpi_printf("\nPM-NONPERIODIC: Allocated %g MByte for FFT kernel(s).\n\n", bytes_tot / (1024.0 * 1024.0));
}

#ifdef PM_ZOOM_OPTIMIZED

void pm_nonperiodic::pmforce_nonperiodic_zoom_optimized_prepare_density(int grnr)
{
  MPI_Status status;

  particle_data *P = Sp->P;

  double to_slab_fac = GRID / ((double)Sp->TotalMeshSize[grnr]);

  part = (part_slab_data *)Mem.mymalloc("part", 8 * (NSource * sizeof(part_slab_data)));
  large_numpart_type *part_sortindex =
      (large_numpart_type *)Mem.mymalloc("part_sortindex", 8 * (NSource * sizeof(large_numpart_type)));

  int ngrid = 0;

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

      MyIntPosType diff[3] = {P[i].IntPos[0] - Sp->ReferenceIntPos[grnr][0], P[i].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                              P[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

      MySignedIntPosType *delta = (MySignedIntPosType *)diff;

      if(delta[0] < Sp->Xmintot[grnr][0] || delta[0] >= Sp->Xmaxtot[grnr][0])
        continue;

      if(delta[1] < Sp->Xmintot[grnr][1] || delta[1] >= Sp->Xmaxtot[grnr][1])
        continue;

      if(delta[2] < Sp->Xmintot[grnr][2] || delta[2] >= Sp->Xmaxtot[grnr][2])
        continue;

      int slab_x = (int)(to_slab_fac * (delta[0] - Sp->Corner[grnr][0]));
      int slab_y = (int)(to_slab_fac * (delta[1] - Sp->Corner[grnr][1]));
      int slab_z = (int)(to_slab_fac * (delta[2] - Sp->Corner[grnr][2]));
      int myngrid;

      myngrid = ngrid;
      ngrid += 1;

      large_numpart_type index_on_grid = ((large_numpart_type)myngrid) * 8;

      for(int xx = 0; xx < 2; xx++)
        for(int yy = 0; yy < 2; yy++)
          for(int zz = 0; zz < 2; zz++)
            {
              int slab_xx = slab_x + xx;
              int slab_yy = slab_y + yy;
              int slab_zz = slab_z + zz;

              if(slab_xx >= GRID)
                slab_xx -= GRID;
              if(slab_yy >= GRID)
                slab_yy -= GRID;
              if(slab_zz >= GRID)
                slab_zz -= GRID;

              large_array_offset offset = FI(slab_xx, slab_yy, slab_zz);

              part[index_on_grid].partindex   = (i << 3) + (xx << 2) + (yy << 1) + zz;
              part[index_on_grid].globalindex = offset;
              part_sortindex[index_on_grid]   = index_on_grid;
              index_on_grid++;
            }
    }

  /* note: num_on_grid will be  8 times larger than the particle number, but num_field_points will generally be much smaller */
  num_on_grid = ((large_numpart_type)ngrid) * 8;

  /* bring the part-field into the order of the accessed cells. This allows the removal of duplicates */
  mycxxsort(part_sortindex, part_sortindex + num_on_grid, pm_nonperiodic_sortindex_comparator(part));

  large_array_offset num_field_points;

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
      int slab = part[part_sortindex[i]].globalindex / (GRID * GRID2);
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

      MyIntPosType diff[3] = {P[pindex].IntPos[0] - Sp->ReferenceIntPos[grnr][0], P[pindex].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                              P[pindex].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

      MySignedIntPosType *delta = (MySignedIntPosType *)diff;

      double dx = to_slab_fac * (delta[0] - Sp->Corner[grnr][0]);
      double dy = to_slab_fac * (delta[1] - Sp->Corner[grnr][1]);
      double dz = to_slab_fac * (delta[2] - Sp->Corner[grnr][2]);

      int slab_x = (int)(dx);
      int slab_y = (int)(dy);
      int slab_z = (int)(dz);

      dx -= slab_x;
      dy -= slab_y;
      dz -= slab_z;

      double weight = P[pindex].getMass();

      localfield_data[part[i + 0].localindex] += weight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 1].localindex] += weight * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 2].localindex] += weight * (1.0 - dx) * dy * (1.0 - dz);
      localfield_data[part[i + 3].localindex] += weight * (1.0 - dx) * dy * dz;
      localfield_data[part[i + 4].localindex] += weight * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 5].localindex] += weight * (dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 6].localindex] += weight * (dx)*dy * (1.0 - dz);
      localfield_data[part[i + 7].localindex] += weight * (dx)*dy * dz;
    }

  rhogrid = (fft_real *)Mem.mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
  for(large_array_offset ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;

  /* exchange data and add contributions to the local mesh-path */
  myMPI_Alltoall(localfield_sendcount, sizeof(size_t), MPI_BYTE, localfield_recvcount, sizeof(size_t), MPI_BYTE, Communicator);

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
          for(large_numpart_type i = 0; i < localfield_recvcount[recvTask]; i++)
            {
              /* determine offset in local FFT slab */
#ifndef FFT_COLUMN_BASED
              large_array_offset offset =
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRID * ((large_array_offset)GRID2);
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
void pm_nonperiodic::pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(int grnr, int dim)
{
#ifdef EVALPOTENTIAL
  double fac = 1.0 / (Sp->FacIntToCoord * Sp->TotalMeshSize[grnr]) / pow(GRID, 3); /* to get potential */
#endif

  particle_data *P = Sp->P;

  MPI_Status status;

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

  double to_slab_fac = GRID / ((double)Sp->TotalMeshSize[grnr]);

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

          for(large_numpart_type i = 0; i < localfield_recvcount[recvTask]; i++)
            {
#ifndef FFT_COLUMN_BASED
              large_array_offset offset =
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRID * ((large_array_offset)GRID2);
#else
              large_array_offset offset = import_globalindex[i] - myplan.firstcol_XY * ((large_array_offset)GRID2);
#endif
              import_data[i] = grid[offset];
            }

          if(level > 0)
            {
              myMPI_Sendrecv(import_data, localfield_recvcount[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A,
                             localfield_data + localfield_offset[recvTask], localfield_sendcount[recvTask] * sizeof(fft_real),
                             MPI_BYTE, recvTask, TAG_NONPERIOD_A, Communicator, &status);

              Mem.myfree(import_globalindex);
              Mem.myfree(import_data);
            }
        }
    }

  /* read out the force/potential values, which all have been assembled in localfield_data */

  int ngrid = (num_on_grid >> 3);

  for(int k = 0; k < ngrid; k++)
    {
      large_numpart_type j = (((large_numpart_type)k) << 3);

      int i = (part[j].partindex >> 3);

#if !defined(HIERARCHICAL_GRAVITY) && defined(TREEPM_NOTIMESPLIT)
      if(!Sp->TimeBinSynchronized[P[i].TimeBinGrav])
        continue;
#endif

      MyIntPosType diff[3] = {P[i].IntPos[0] - Sp->ReferenceIntPos[grnr][0], P[i].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                              P[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

      MySignedIntPosType *delta = (MySignedIntPosType *)diff;

      double dx = to_slab_fac * (delta[0] - Sp->Corner[grnr][0]);
      double dy = to_slab_fac * (delta[1] - Sp->Corner[grnr][1]);
      double dz = to_slab_fac * (delta[2] - Sp->Corner[grnr][2]);

      int slab_x = (int)(dx);
      int slab_y = (int)(dy);
      int slab_z = (int)(dz);

      dx -= slab_x;
      dy -= slab_y;
      dz -= slab_z;

      double value = +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
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
          P[i].Potential += value * fac;
#endif
        }
      else
        {
          Sp->P[i].GravAccel[dim] += value;
        }
    }
}

#else

void pm_nonperiodic::pmforce_nonperiodic_uniform_optimized_prepare_density(int grnr)
{
  double to_slab_fac = GRID / ((double)Sp->TotalMeshSize[grnr]);

  Sndpm_count = (size_t *)Mem.mymalloc("Sndpm_count", NTask * sizeof(size_t));
  Sndpm_offset = (size_t *)Mem.mymalloc("Sndpm_offset", NTask * sizeof(size_t));
  Rcvpm_count = (size_t *)Mem.mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *)Mem.mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

#ifdef FFT_COLUMN_BASED
  int columns = GRIDX * GRIDY;
  int avg = (columns - 1) / NTask + 1;
  int exc = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol = tasklastsection * avg;
#endif

  /* determine the slabs/columns each particles accesses */
  {
    size_t *send_count = Sndpm_count;

    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    for(int idx = 0; idx < NSource; idx++)
      {
        int i = Sp->get_active_index(idx);

        if(Sp->P[i].Ti_Current != All.Ti_Current)
          Sp->drift_particle(&Sp->P[i], &Sp->SphP[i], All.Ti_Current);

        MyIntPosType diff[3] = {Sp->P[i].IntPos[0] - Sp->ReferenceIntPos[grnr][0], Sp->P[i].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                                Sp->P[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

        MySignedIntPosType *delta = (MySignedIntPosType *)diff;

        if(delta[0] < Sp->Xmintot[grnr][0] || delta[0] >= Sp->Xmaxtot[grnr][0])
          continue;

        if(delta[1] < Sp->Xmintot[grnr][1] || delta[1] >= Sp->Xmaxtot[grnr][1])
          continue;

        if(delta[2] < Sp->Xmintot[grnr][2] || delta[2] >= Sp->Xmaxtot[grnr][2])
          continue;

        int slab_x = (int)(to_slab_fac * (delta[0] - Sp->Corner[grnr][0]));
        int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        send_count[task0]++;
        if(task0 != task1)
          send_count[task1]++;
#else
        int slab_y  = (int)(to_slab_fac * (delta[1] - Sp->Corner[grnr][1]));
        int slab_yy = slab_y + 1;

        int column0 = slab_x * GRID + slab_y;
        int column1 = slab_x * GRID + slab_yy;
        int column2 = slab_xx * GRID + slab_y;
        int column3 = slab_xx * GRID + slab_yy;

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

        send_count[task0]++;
        if(task1 != task0)
          send_count[task1]++;
        if(task2 != task1 && task2 != task0)
          send_count[task2]++;
        if(task3 != task0 && task3 != task1 && task3 != task2)
          send_count[task3]++;
#endif
      }
  }

  /* collect thread-specific offset table and collect the results from the other threads */
  Sndpm_offset[0] = 0;
  for(int i = 1; i < NTask; i++)
    {
      int ind = i;
      int ind_prev = i - 1;

      Sndpm_offset[ind] = Sndpm_offset[ind_prev] + Sndpm_count[ind_prev];
    }

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
  partin = (partbuf *)Mem.mymalloc("partin", nimport * sizeof(partbuf));
  partout = (partbuf *)Mem.mymalloc("partout", nexport * sizeof(partbuf));

  {
    size_t *send_count = Sndpm_count;
    size_t *send_offset = Sndpm_offset;

    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    /* fill export buffer */
    for(int idx = 0; idx < NSource; idx++)
      {
        int i = Sp->get_active_index(idx);

        MyIntPosType diff[3] = {Sp->P[i].IntPos[0] - Sp->ReferenceIntPos[grnr][0], Sp->P[i].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                                Sp->P[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

        MySignedIntPosType *delta = (MySignedIntPosType *)diff;

        if(delta[0] < Sp->Xmintot[grnr][0] || delta[0] >= Sp->Xmaxtot[grnr][0])
          continue;

        if(delta[1] < Sp->Xmintot[grnr][1] || delta[1] >= Sp->Xmaxtot[grnr][1])
          continue;

        if(delta[2] < Sp->Xmintot[grnr][2] || delta[2] >= Sp->Xmaxtot[grnr][2])
          continue;

        int slab_x = (int)(to_slab_fac * (delta[0] - Sp->Corner[grnr][0]));
        int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        size_t ind0 = send_offset[task0] + send_count[task0]++;
        partout[ind0].Mass = Sp->P[i].getMass();
        for(int j = 0; j < 3; j++)
          partout[ind0].IntPos[j] = Sp->P[i].IntPos[j];

        if(task0 != task1)
          {
            size_t ind1 = send_offset[task1] + send_count[task1]++;
            partout[ind1].Mass = Sp->P[i].getMass();
            for(int j = 0; j < 3; j++)
              partout[ind1].IntPos[j] = Sp->P[i].IntPos[j];
          }
#else
        int slab_y  = (int)(to_slab_fac * (delta[1] - Sp->Corner[grnr][1]));
        int slab_yy = slab_y + 1;

        int column0 = slab_x * GRID + slab_y;
        int column1 = slab_x * GRID + slab_yy;
        int column2 = slab_xx * GRID + slab_y;
        int column3 = slab_xx * GRID + slab_yy;

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

        size_t ind0        = send_offset[task0] + send_count[task0]++;
        partout[ind0].Mass = Sp->P[i].getMass();
        for(int j = 0; j < 3; j++)
          partout[ind0].IntPos[j] = Sp->P[i].IntPos[j];

        if(task1 != task0)
          {
            size_t ind1        = send_offset[task1] + send_count[task1]++;
            partout[ind1].Mass = Sp->P[i].getMass();
            for(int j = 0; j < 3; j++)
              partout[ind1].IntPos[j] = Sp->P[i].IntPos[j];
          }
        if(task2 != task1 && task2 != task0)
          {
            size_t ind2        = send_offset[task2] + send_count[task2]++;
            partout[ind2].Mass = Sp->P[i].getMass();
            for(int j = 0; j < 3; j++)
              partout[ind2].IntPos[j] = Sp->P[i].IntPos[j];
          }
        if(task3 != task0 && task3 != task1 && task3 != task2)
          {
            size_t ind3        = send_offset[task3] + send_count[task3]++;
            partout[ind3].Mass = Sp->P[i].getMass();
            for(int j = 0; j < 3; j++)
              partout[ind3].IntPos[j] = Sp->P[i].IntPos[j];
          }
#endif
      }
  }

  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange particle data */
  myMPI_Alltoallv(partout, Sndpm_count, Sndpm_offset, partin, Rcvpm_count, Rcvpm_offset, sizeof(partbuf), flag_big_all, Communicator);

  Mem.myfree(partout);

  /* allocate density field */
  rhogrid = (fft_real *)Mem.mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
  for(size_t ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;

#ifndef FFT_COLUMN_BASED
  /* bin particle data onto mesh, in multi-threaded fashion */
  {
    int first_y, count_y;
    subdivide_evenly(GRID, 1, 0, &first_y, &count_y);
    int last_y = first_y + count_y - 1;

    for(size_t i = 0; i < nimport; i++)
      {
        MyIntPosType diff[3] = {partin[i].IntPos[0] - Sp->ReferenceIntPos[grnr][0], partin[i].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                                partin[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

        MySignedIntPosType *delta = (MySignedIntPosType *)diff;

        double dy = to_slab_fac * (delta[1] - Sp->Corner[grnr][1]);
        int slab_y = (int)(dy);
        dy -= slab_y;

        int slab_yy = slab_y + 1;
        int flag_slab_y, flag_slab_yy;

        if(slab_y >= first_y && slab_y <= last_y)
          flag_slab_y = 1;
        else
          flag_slab_y = 0;

        if(slab_yy >= first_y && slab_yy <= last_y)
          flag_slab_yy = 1;
        else
          flag_slab_yy = 0;

        if(flag_slab_y || flag_slab_yy)
          {
            double mass = partin[i].Mass;

            double dx = to_slab_fac * (delta[0] - Sp->Corner[grnr][0]);
            double dz = to_slab_fac * (delta[2] - Sp->Corner[grnr][2]);
            int slab_x = (int)(dx);
            int slab_z = (int)(dz);
            dx -= slab_x;
            dz -= slab_z;

            int slab_xx = slab_x + 1;
            int slab_zz = slab_z + 1;

            int flag_slab_x, flag_slab_xx;

            if(myplan.slab_to_task[slab_x] == ThisTask)
              {
                slab_x -= myplan.first_slab_x_of_task[ThisTask];
                flag_slab_x = 1;
              }
            else
              flag_slab_x = 0;

            if(myplan.slab_to_task[slab_xx] == ThisTask)
              {
                slab_xx -= myplan.first_slab_x_of_task[ThisTask];
                flag_slab_xx = 1;
              }
            else
              flag_slab_xx = 0;

            if(flag_slab_x)
              {
                if(flag_slab_y)
                  {
                    rhogrid[FI(slab_x, slab_y, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                    rhogrid[FI(slab_x, slab_y, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
                  }

                if(flag_slab_yy)
                  {
                    rhogrid[FI(slab_x, slab_yy, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                    rhogrid[FI(slab_x, slab_yy, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
                  }
              }

            if(flag_slab_xx)
              {
                if(flag_slab_y)
                  {
                    rhogrid[FI(slab_xx, slab_y, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                    rhogrid[FI(slab_xx, slab_y, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
                  }

                if(flag_slab_yy)
                  {
                    rhogrid[FI(slab_xx, slab_yy, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                    rhogrid[FI(slab_xx, slab_yy, slab_zz)] += (mass * (dx) * (dy) * (dz));
                  }
              }
          }
      }
  }

#else

  struct data_cols
  {
    int col0, col1, col2, col3;
    double dx, dy;
  };

  data_cols *aux = (data_cols *)Mem.mymalloc("aux", nimport * sizeof(data_cols));

  for(int i = 0; i < nimport; i++)
    {
      MyIntPosType diff[3] = {partin[i].IntPos[0] - Sp->ReferenceIntPos[grnr][0], partin[i].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                              partin[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

      MySignedIntPosType *delta = (MySignedIntPosType *)diff;

      aux[i].dx = to_slab_fac * (delta[0] - Sp->Corner[grnr][0]);
      aux[i].dy = to_slab_fac * (delta[1] - Sp->Corner[grnr][1]);

      int slab_x = (int)(aux[i].dx);
      int slab_y = (int)(aux[i].dy);

      aux[i].dx -= slab_x;
      aux[i].dy -= slab_y;

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;

      aux[i].col0 = slab_x * GRID + slab_y;
      aux[i].col1 = slab_x * GRID + slab_yy;
      aux[i].col2 = slab_xx * GRID + slab_y;
      aux[i].col3 = slab_xx * GRID + slab_yy;
    }

  {
    int first_col, last_col, count_col;
    subdivide_evenly(myplan.ncol_XY, 1, 0, &first_col, &count_col);
    last_col = first_col + count_col - 1;
    first_col += myplan.firstcol_XY;
    last_col += myplan.firstcol_XY;

    for(int i = 0; i < nimport; i++)
      {
        int flag0, flag1, flag2, flag3;
        int col0 = aux[i].col0;
        int col1 = aux[i].col1;
        int col2 = aux[i].col2;
        int col3 = aux[i].col3;

        if(col0 >= first_col && col0 <= last_col)
          flag0 = 1;
        else
          flag0 = 0;

        if(col1 >= first_col && col1 <= last_col)
          flag1 = 1;
        else
          flag1 = 0;

        if(col2 >= first_col && col2 <= last_col)
          flag2 = 1;
        else
          flag2 = 0;

        if(col3 >= first_col && col3 <= last_col)
          flag3 = 1;
        else
          flag3 = 0;

        if(flag0 || flag1 || flag2 || flag3)
          {
            double mass = partin[i].Mass;

            double dx = aux[i].dx;
            double dy = aux[i].dy;

            MySignedIntPosType deltaz = (MySignedIntPosType)(partin[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]);

            double dz = to_slab_fac * (deltaz - Sp->Corner[grnr][2]);
            int slab_z = (int)(dz);
            dz -= slab_z;
            int slab_zz = slab_z + 1;

            if(flag0)
              {
                rhogrid[FC(col0, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                rhogrid[FC(col0, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
              }

            if(flag1)
              {
                rhogrid[FC(col1, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                rhogrid[FC(col1, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
              }

            if(flag2)
              {
                rhogrid[FC(col2, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                rhogrid[FC(col2, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
              }

            if(flag3)
              {
                rhogrid[FC(col3, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                rhogrid[FC(col3, slab_zz)] += (mass * (dx) * (dy) * (dz));
              }
          }
      }
  }

  Mem.myfree(aux);

#endif
}

/* If dim<0, this function reads out the potential, otherwise Cartesian force components.
 */
void pm_nonperiodic::pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(int grnr, int dim)
{
#ifdef EVALPOTENTIAL
  double fac = 1.0 / (Sp->FacIntToCoord * Sp->TotalMeshSize[grnr]) / pow(GRID, 3); /* to get potential */
#endif

  double to_slab_fac = GRID / ((double)Sp->TotalMeshSize[grnr]);

  double *flistin = (double *)Mem.mymalloc("flistin", nimport * sizeof(double));
  double *flistout = (double *)Mem.mymalloc("flistout", nexport * sizeof(double));

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

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

      MyIntPosType diff[3] = {partin[i].IntPos[0] - Sp->ReferenceIntPos[grnr][0], partin[i].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                              partin[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

      MySignedIntPosType *delta = (MySignedIntPosType *)diff;

      double dx = to_slab_fac * (delta[0] - Sp->Corner[grnr][0]);
      double dy = to_slab_fac * (delta[1] - Sp->Corner[grnr][1]);
      double dz = to_slab_fac * (delta[2] - Sp->Corner[grnr][2]);

      int slab_x = (int)(dx);
      int slab_y = (int)(dy);
      int slab_z = (int)(dz);

      dx -= slab_x;
      dy -= slab_y;
      dz -= slab_z;

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

#ifndef FFT_COLUMN_BASED
      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += +grid[FI(slab_x, slab_y, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_y, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_x, slab_yy, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_yy, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += +grid[FI(slab_xx, slab_y, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_y, slab_zz)] * (dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_xx, slab_yy, slab_z)] * (dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_yy, slab_zz)] * (dx) * (dy) * (dz);
        }
#else
      int column0 = slab_x * GRID + slab_y;
      int column1 = slab_x * GRID + slab_yy;
      int column2 = slab_xx * GRID + slab_y;
      int column3 = slab_xx * GRID + slab_yy;

      if(column0 >= myplan.firstcol_XY && column0 <= myplan.lastcol_XY)
        {
          flistin[i] += +grid[FC(column0, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FC(column0, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz);
        }
      if(column1 >= myplan.firstcol_XY && column1 <= myplan.lastcol_XY)
        {
          flistin[i] +=
              +grid[FC(column1, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FC(column1, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.firstcol_XY && column2 <= myplan.lastcol_XY)
        {
          flistin[i] +=
              +grid[FC(column2, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FC(column2, slab_zz)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.firstcol_XY && column3 <= myplan.lastcol_XY)
        {
          flistin[i] += +grid[FC(column3, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FC(column3, slab_zz)] * (dx) * (dy) * (dz);
        }
#endif
    }

  /* exchange the potential component data */
  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(double) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange  data */
  myMPI_Alltoallv(flistin, Rcvpm_count, Rcvpm_offset, flistout, Sndpm_count, Sndpm_offset, sizeof(double), flag_big_all, Communicator);

  /* now assign them to the correct particles */

  size_t *send_count = Sndpm_count;
  size_t *send_offset = Sndpm_offset;

  for(int j = 0; j < NTask; j++)
    send_count[j] = 0;

  for(int idx = 0; idx < NSource; idx++)
    {
      int i = Sp->get_active_index(idx);

      MyIntPosType diff[3] = {Sp->P[i].IntPos[0] - Sp->ReferenceIntPos[grnr][0], Sp->P[i].IntPos[1] - Sp->ReferenceIntPos[grnr][1],
                              Sp->P[i].IntPos[2] - Sp->ReferenceIntPos[grnr][2]};

      MySignedIntPosType *delta = (MySignedIntPosType *)diff;

      if(delta[0] < Sp->Xmintot[grnr][0] || delta[0] >= Sp->Xmaxtot[grnr][0])
        continue;

      if(delta[1] < Sp->Xmintot[grnr][1] || delta[1] >= Sp->Xmaxtot[grnr][1])
        continue;

      if(delta[2] < Sp->Xmintot[grnr][2] || delta[2] >= Sp->Xmaxtot[grnr][2])
        continue;

      int slab_x = (int)(to_slab_fac * (delta[0] - Sp->Corner[grnr][0]));
      int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
      int task0 = myplan.slab_to_task[slab_x];
      int task1 = myplan.slab_to_task[slab_xx];

      double value = flistout[send_offset[task0] + send_count[task0]++];

      if(task0 != task1)
        value += flistout[send_offset[task1] + send_count[task1]++];
#else
      int slab_y = (int)(to_slab_fac * (delta[1] - Sp->Corner[grnr][1]));
      int slab_yy = slab_y + 1;

      int column0 = slab_x * GRID + slab_y;
      int column1 = slab_x * GRID + slab_yy;
      int column2 = slab_xx * GRID + slab_y;
      int column3 = slab_xx * GRID + slab_yy;

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
#endif

#if !defined(HIERARCHICAL_GRAVITY) && defined(TREEPM_NOTIMESPLIT)
      if(!Sp->TimeBinSynchronized[Sp->P[i].TimeBinGrav])
        continue;
#endif

      if(dim < 0)
        {
#ifdef EVALPOTENTIAL
          Sp->P[i].Potential += value * fac;
#endif
        }
      else
        {
          Sp->P[i].GravAccel[dim] += value;
        }
    }

  Mem.myfree(flistout);
  Mem.myfree(flistin);
}
#endif

/*! Calculates the long-range non-periodic forces using the PM method.  The
 *  potential is Gaussian filtered with Asmth, given in mesh-cell units. The
 *  potential is finite differenced using a 4-point finite differencing
 *  formula to obtain the force fields, which are then interpolated to the
 *  particle positions. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The CIC kernel is deconvolved.
 */
int pm_nonperiodic::pmforce_nonperiodic(int grnr)
{
  double tstart = Logs.second();

  mpi_printf("PM-NONPERIODIC: Starting non-periodic PM calculation (Rcut=%g, grid=%d)  presently allocated=%g MB.\n", Sp->Rcut[grnr],
             grnr, Mem.getAllocatedBytesInMB());

  if(grnr == 1 && Sp->Rcut[1] > Sp->Rcut[0])
    Terminate(
        "We have Sp->Rcut[1]=%g  >  Sp->Rcut[0]=%g, which means that the high-res cut-off is larger than the normal one... this is "
        "not good.",
        Sp->Rcut[1], Sp->Rcut[0]);

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

  double fac = 1.0 / (Sp->FacIntToCoord * Sp->TotalMeshSize[grnr]) / pow(GRID, 3); /* to get potential */
  fac *= 1 / (2 * (Sp->FacIntToCoord * Sp->TotalMeshSize[grnr]) / GRID);           /* for finite differencing */

#ifdef PM_ZOOM_OPTIMIZED
  pmforce_nonperiodic_zoom_optimized_prepare_density(grnr);
#else
  pmforce_nonperiodic_uniform_optimized_prepare_density(grnr);
#endif

  /* allocate the memory to hold the FFT fields */
  forcegrid = (fft_real *)Mem.mymalloc("forcegrid", maxfftsize * sizeof(fft_real));

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

  /* multiply with kernel in Fourier space */
  /* multiply with the Fourier transform of the Green's function (kernel) */

  /* multiply with Green's function in order to obtain the potential */

#ifdef FFT_COLUMN_BASED
  for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
    {
#else
  for(int x = 0; x < GRID; x++)
    for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
      for(int z = 0; z < GRIDz; z++)
        {
#endif

#ifndef FFT_COLUMN_BASED
      large_array_offset ip = ((large_array_offset)GRIDz) * (GRID * (y - myplan.slabstart_y) + x) + z;
#endif

      double re = fft_of_rhogrid[ip][0] * fft_of_kernel[grnr][ip][0] - fft_of_rhogrid[ip][1] * fft_of_kernel[grnr][ip][1];
      double im = fft_of_rhogrid[ip][0] * fft_of_kernel[grnr][ip][1] + fft_of_rhogrid[ip][1] * fft_of_kernel[grnr][ip][0];

      fft_of_rhogrid[ip][0] = re;
      fft_of_rhogrid[ip][1] = im;
    }

    /* Do the inverse FFT to get the potential */

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, rhogrid, workspace, -1);
#else
  my_column_based_fft(&myplan, workspace, rhogrid, -1);
#endif

  /* Now rhogrid holds the potential */

#ifdef EVALPOTENTIAL
#ifdef PM_ZOOM_OPTIMIZED
  pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(grnr, -1);
#else
  pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(grnr, -1);
#endif
#endif

  /* get the force components by finite differencing of the potential for each dimension,
   * and send the results back to the right CPUs
   */
  for(int dim = 2; dim >= 0; dim--) /* Calculate each component of the force. */
    {
      /* we do the x component last, because for differencing the potential in the x-direction, we need to construct the transpose */
#ifndef FFT_COLUMN_BASED
      if(dim == 0)
        my_slab_transposeA(&myplan, rhogrid, forcegrid); /* compute the transpose of the potential field for finite differencing */

      for(int y = 2; y < GRID / 2 - 2; y++)
        for(int x = 0; x < myplan.nslab_x; x++)
          if(x + myplan.slabstart_x >= 2 && x + myplan.slabstart_x < GRID / 2 - 2)
            for(int z = 2; z < GRID / 2 - 2; z++)
              {
                int yrr = y, yll = y, yr = y, yl = y;
                int zrr = z, zll = z, zr = z, zl = z;

                switch(dim)
                  {
                    case 0: /* note: for the x-direction, we difference the transposed direction (y) */
                    case 1:
                      yr  = y + 1;
                      yl  = y - 1;
                      yrr = y + 2;
                      yll = y - 2;

                      break;
                    case 2:
                      zr  = z + 1;
                      zl  = z - 1;
                      zrr = z + 2;
                      zll = z - 2;

                      break;
                  }

                if(dim == 0)
                  forcegrid[TI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[TI(x, yl, zl)] - rhogrid[TI(x, yr, zr)]) -
                                                  (1.0 / 6) * (rhogrid[TI(x, yll, zll)] - rhogrid[TI(x, yrr, zrr)]));
                else
                  forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, yl, zl)] - rhogrid[FI(x, yr, zr)]) -
                                                  (1.0 / 6) * (rhogrid[FI(x, yll, zll)] - rhogrid[FI(x, yrr, zrr)]));
              }

      if(dim == 0)
        my_slab_transposeB(&myplan, forcegrid, rhogrid); /* reverse the transpose from above */
#else
      fft_real *scratch;

      if(dim != 2)
        {
          scratch = (fft_real *)Mem.mymalloc("scratch", myplan.fftsize * sizeof(fft_real)); /* need a third field as scratch space */
          memcpy(scratch, rhogrid, myplan.fftsize * sizeof(fft_real));

          if(dim == 1)
            my_fft_swap23(&myplan, scratch, forcegrid);
          else
            my_fft_swap13(&myplan, scratch, forcegrid);
        }

      int ncols;
      if(dim == 2)
        ncols = myplan.ncol_XY;
      else if(dim == 1)
        ncols = myplan.ncol_XZ;
      else
        ncols = myplan.ncol_ZY;

      for(int i = 0; i < ncols; i++)
        {
          fft_real *forcep, *potp;

          if(dim != 2)
            {
              forcep = &scratch[GRID * i];
              potp = &forcegrid[GRID * i];
            }
          else
            {
              forcep = &forcegrid[GRID2 * i];
              potp = &rhogrid[GRID2 * i];
            }

          for(int z = 2; z < GRID / 2 - 2; z++)
            {
              int zr = z + 1;
              int zl = z - 1;
              int zrr = z + 2;
              int zll = z - 2;

              forcep[z] = fac * ((4.0 / 3) * (potp[zl] - potp[zr]) - (1.0 / 6) * (potp[zll] - potp[zrr]));
            }
        }

      if(dim != 2)
        {
          if(dim == 1)
            my_fft_swap23back(&myplan, scratch, forcegrid);
          else
            my_fft_swap13back(&myplan, scratch, forcegrid);

          Mem.myfree(scratch);
        }
#endif

#ifdef PM_ZOOM_OPTIMIZED
      pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(grnr, dim);
#else
      pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(grnr, dim);
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
  Mem.myfree(partin);
  Mem.myfree(Rcvpm_offset);
  Mem.myfree(Rcvpm_count);
  Mem.myfree(Sndpm_offset);
  Mem.myfree(Sndpm_count);
#endif

  double tend = Logs.second();

  mpi_printf("PM-NONPERIODIC: done.  (took %g seconds)\n", Logs.timediff(tstart, tend));

  return 0;
}

/*! This function sets-up the Greens function for the non-periodic potential
 *  in real space, and then converts it to Fourier space by means of a FFT.
 */
void pm_nonperiodic::pm_setup_nonperiodic_kernel(void)
{
  mpi_printf("PM-NONPERIODIC: Setting up non-periodic PM kernel(s) (GRID=%d)  presently allocated=%g MB).\n", (int)GRID,
             Mem.getAllocatedBytesInMB());

  /* now set up kernel and its Fourier transform */

#if !defined(PERIODIC)
  for(size_t i = 0; i < maxfftsize; i++) /* clear local field */
    kernel[0][i] = 0;

#ifndef FFT_COLUMN_BASED
  for(int i = myplan.slabstart_x; i < (myplan.slabstart_x + myplan.nslab_x); i++)
    for(int j = 0; j < GRID; j++)
      {
#else
  for(int c = myplan.firstcol_XY; c < (myplan.firstcol_XY + myplan.ncol_XY); c++)
    {
      int i = c / GRID;
      int j = c % GRID;
#endif
        for(int k = 0; k < GRID; k++)
          {
            double xx = ((double)i) / GRID;
            double yy = ((double)j) / GRID;
            double zz = ((double)k) / GRID;

            if(xx >= 0.5)
              xx -= 1.0;
            if(yy >= 0.5)
              yy -= 1.0;
            if(zz >= 0.5)
              zz -= 1.0;

            double r = sqrt(xx * xx + yy * yy + zz * zz);

            double u = 0.5 * r / (((double)ASMTH) / GRID);

            double fac = 1 - erfc(u);

#ifndef FFT_COLUMN_BASED
            size_t ip = FI(i - myplan.slabstart_x, j, k);
#else
          size_t ip = FC(c, k);
#endif
            if(r > 0)
              kernel[0][ip] = -fac / r;
            else
              kernel[0][ip] = -1 / (sqrt(M_PI) * (((double)ASMTH) / GRID));
          }
      }

  {
    fft_real *workspc = (fft_real *)Mem.mymalloc("workspc", maxfftsize * sizeof(fft_real));
    /* Do the FFT of the kernel */
#ifndef FFT_COLUMN_BASED
    my_slab_based_fft(&myplan, kernel[0], workspc, 1);
#else
    my_column_based_fft(&myplan, kernel[0], workspc, 1); /* result is in workspace, not in kernel */
    memcpy(kernel[0], workspc, maxfftsize * sizeof(fft_real));
#endif
    Mem.myfree(workspc);
  }

#endif

#if defined(PLACEHIGHRESREGION)

  for(int i = 0; i < maxfftsize; i++) /* clear local field */
    kernel[1][i] = 0;

#ifndef FFT_COLUMN_BASED
  for(int i = myplan.slabstart_x; i < (myplan.slabstart_x + myplan.nslab_x); i++)
    for(int j = 0; j < GRID; j++)
      {
#else
  for(int c = myplan.firstcol_XY; c < (myplan.firstcol_XY + myplan.ncol_XY); c++)
    {
      int i = c / GRID;
      int j = c % GRID;
#endif
        for(int k = 0; k < GRID; k++)
          {
            double xx = ((double)i) / GRID;
            double yy = ((double)j) / GRID;
            double zz = ((double)k) / GRID;

            if(xx >= 0.5)
              xx -= 1.0;
            if(yy >= 0.5)
              yy -= 1.0;
            if(zz >= 0.5)
              zz -= 1.0;

            double r = sqrt(xx * xx + yy * yy + zz * zz);

            double u = 0.5 * r / (((double)ASMTH) / GRID);

            double fac = erfc(u * Sp->Asmth[1] / Sp->Asmth[0]) - erfc(u);

#ifndef FFT_COLUMN_BASED
            size_t ip = FI(i - myplan.slabstart_x, j, k);
#else
          size_t ip = FC(c, k);
#endif

            if(r > 0)
              kernel[1][ip] = -fac / r;
            else
              {
                fac           = 1 - Sp->Asmth[1] / Sp->Asmth[0];
                kernel[1][ip] = -fac / (sqrt(M_PI) * (((double)ASMTH) / GRID));
              }
          }
#ifndef FFT_COLUMN_BASED
      }
#else
    }
#endif

  {
    fft_real *workspc = (fft_real *)Mem.mymalloc("workspc", maxfftsize * sizeof(fft_real));
    /* Do the FFT of the kernel */
#ifndef FFT_COLUMN_BASED
    my_slab_based_fft(&myplan, kernel[1], workspc, 1);
#else
    my_column_based_fft(&myplan, kernel[1], workspc, 1); /* result is in workspace, not in kernel */
    memcpy(kernel[1], workspc, maxfftsize * sizeof(fft_real));
#endif
    Mem.myfree(workspc);
  }

#endif

  /* deconvolve the Greens function twice with the CIC kernel */
#ifdef FFT_COLUMN_BASED

  for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
    {
      large_array_offset ipcell = ip + myplan.transposed_firstcol * GRID;
      int y                     = ipcell / (GRID * GRIDz);
      int yr                    = ipcell % (GRID * GRIDz);
      int z                     = yr / GRID;
      int x                     = yr % GRID;
#else
  for(int x = 0; x < GRID; x++)
    for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
      for(int z = 0; z < GRIDz; z++)
        {
#endif
      double kx, ky, kz;

      if(x > GRID / 2)
        kx = x - GRID;
      else
        kx = x;
      if(y > GRID / 2)
        ky = y - GRID;
      else
        ky = y;
      if(z > GRID / 2)
        kz = z - GRID;
      else
        kz = z;

      double k2 = kx * kx + ky * ky + kz * kz;

      if(k2 > 0)
        {
          double fx = 1, fy = 1, fz = 1;

          if(kx != 0)
            {
              fx = (M_PI * kx) / GRID;
              fx = sin(fx) / fx;
            }
          if(ky != 0)
            {
              fy = (M_PI * ky) / GRID;
              fy = sin(fy) / fy;
            }
          if(kz != 0)
            {
              fz = (M_PI * kz) / GRID;
              fz = sin(fz) / fz;
            }

          double ff = 1 / (fx * fy * fz);
          ff        = ff * ff * ff * ff;

#ifndef FFT_COLUMN_BASED
          large_array_offset ip = ((large_array_offset)GRIDz) * (GRID * (y - myplan.slabstart_y) + x) + z;
#endif
#if !defined(PERIODIC)
          fft_of_kernel[0][ip][0] *= ff;
          fft_of_kernel[0][ip][1] *= ff;
#endif
#if defined(PLACEHIGHRESREGION)
          fft_of_kernel[1][ip][0] *= ff;
          fft_of_kernel[1][ip][1] *= ff;
#endif
        }
    }

  /* end deconvolution */
}

#endif
