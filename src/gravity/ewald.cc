/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file ewald.cc
 *
 *  \brief Code for Ewald correction computations.
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
#include "../gravity/ewald.h"
#include "../gravity/ewaldtensors.h"
#include "../gravtree/gravtree.h"
#include "../io/io.h"
#include "../main/simulation.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"

/*!
 *  This file contains the computation of the Ewald correction table, and the corresponding lookup functions.
 *
 *   in D0phi we store the correction potential:
 *
 *      phi = 1/x + pi/alpha^2 - sum_q (erfc(alpha |x-q|)/|x-q|)  - 4pi/V sum_k exp(-k^2/(4alpha^2))/k^2 cos(k*x)
 *
 *   in D1phi we store the first derivative of correction potential
 *
 *      dphi/dx_i
 *
 *   in D2phi we store the correction tensor (second derivatives of correction potential)
 *
 *      d2phi/(dx_i dx_j)
 *
 *   in D3phi we store the third order correction tensor (third derivatives of correction potential)
 *
 *      d3phi/(dx_i dx_j dx_k)
 *
 *   and so on also for D4phi and D5phi
 */

/*! \brief This function initializes tables with the correction force and the
 *  correction potential due to the periodic images of a point mass located
 *  at the origin.
 *
 *  These corrections are obtained by Ewald summation. (See for example
 *  Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231) The correction fields
 *  are used to obtain the full periodic force if periodic boundaries
 *  combined with the pure tree algorithm are used. For the TreePM/FMM-PM
 *  algorithms, the Ewald correction is not used.
 *
 *  The correction fields are stored on disk once they are computed. If a
 *  corresponding file is found, they are loaded from disk to speed up the
 *  initialization. The Ewald summation issrc/gravtree_forcetest.c done in parallel, i.e. the
 *  processors share the work to compute the tables if needed.
 */
void ewald::ewald_init(void)
{
  mpi_printf("EWALD: initialize Ewald correction...\n");

  RegionLen     = All.BoxSize;
  FacCoordToInt = pow(2.0, BITS_FOR_POSITIONS) / RegionLen;
  FacIntToCoord = RegionLen / pow(2.0, BITS_FOR_POSITIONS);

  Ewd = (ewald_data *)Mem.mymalloc("Ewd", sizeof(ewald_data) * (ENX + 1) * (ENY + 1) * (ENZ + 1));

  char buf[MAXLEN_PATH_EXTRA];
  snprintf(buf, MAXLEN_PATH_EXTRA, "ewald_table_%d-%d-%d_%d-%d-%d_precision%d-order%d.dat", LONG_X, LONG_Y, LONG_Z, ENX, ENY, ENZ,
           (int)sizeof(MyReal), HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER);

  int recomputeflag = 0;

  if(ThisTask == 0)
    {
      FILE *fd;
      if((fd = fopen(buf, "r")))
        {
          mpi_printf("\nEWALD: reading Ewald tables from file `%s'\n", buf);

          ewald_header tabh;
          my_fread(&tabh, sizeof(ewald_header), 1, fd);

#ifndef GRAVITY_TALLBOX
          int ewaldtype = -1;
#else
          int ewaldtype = GRAVITY_TALLBOX + 1;
#endif
          if(tabh.resx != ENX || tabh.resy != ENY || tabh.resz != ENZ || tabh.varsize != sizeof(MyFloat) ||
             tabh.ewaldtype != ewaldtype)
            {
              mpi_printf("\nEWALD: something's wrong with this table file. Discarding it.\n");
              recomputeflag = 1;
            }
          else
            {
              my_fread(Ewd, sizeof(ewald_data), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);

              recomputeflag = 0;
            }
          fclose(fd);
        }
      else
        recomputeflag = 1;
    }

  MPI_Bcast(&recomputeflag, 1, MPI_INT, 0, Communicator);

  if(recomputeflag)
    {
      mpi_printf("\nEWALD: No usable Ewald tables in file `%s' found. Recomputing them...\n", buf);

      /* ok, let's recompute things. Actually, we do that in parallel. */

      int size = (ENX + 1) * (ENY + 1) * (ENZ + 1);
      int first, count;

      subdivide_evenly(size, NTask, ThisTask, &first, &count);

      for(int n = first; n < first + count; n++)
        {
          int i = n / ((ENY + 1) * (ENZ + 1));
          int j = (n - i * (ENY + 1) * (ENZ + 1)) / (ENZ + 1);
          int k = (n - i * (ENY + 1) * (ENZ + 1) - j * (ENZ + 1));

          if(ThisTask == 0)
            {
              if(count > 20)
                if(((n - first) % (count / 20)) == 0)
                  {
                    printf("%4.1f percent done\n", (n - first) / (count / 100.0));
                    myflush(stdout);
                  }
            }

          double xx = 0.5 * DBX * (1.0 / LONG_X) * ((double)i) / ENX;
          double yy = 0.5 * DBY * (1.0 / LONG_Y) * ((double)j) / ENY;
          double zz = 0.5 * DBZ * (1.0 / LONG_Z) * ((double)k) / ENZ;

          ewald_data *ewdp = Ewd + ewd_offset(i, j, k);

#ifndef GRAVITY_TALLBOX
          ewdp->D0phi = ewald_D0(xx, yy, zz);
          ewdp->D1phi = ewald_D1(xx, yy, zz);
          ewdp->D2phi = ewald_D2(xx, yy, zz);
          ewdp->D3phi = ewald_D3(xx, yy, zz);
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 4
          ewdp->D4phi = ewald_D4(xx, yy, zz);
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 5
          ewdp->D5phi = ewald_D5(xx, yy, zz);
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 6
          ewdp->D6phi = ewald_D6(xx, yy, zz);
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 7
          ewdp->D7phi = ewald_D7(xx, yy, zz);
#endif
#else

          switch(GRAVITY_TALLBOX)
            {
              case 0:
                {
                  ewdp->D0phi = ewald_D0(yy, zz, xx);
                  auto D1phi  = ewald_D1(yy, zz, xx);
                  auto D2phi  = ewald_D2(yy, zz, xx);
                  auto D3phi  = ewald_D3(yy, zz, xx);

                  ewdp->D1phi[vX] = D1phi[vZ];
                  ewdp->D1phi[vY] = D1phi[vX];
                  ewdp->D1phi[vZ] = D1phi[vY];

                  ewdp->D2phi[qXX] = D2phi[qZZ];
                  ewdp->D2phi[qXY] = D2phi[qZX];
                  ewdp->D2phi[qXZ] = D2phi[qZY];
                  ewdp->D2phi[qYY] = D2phi[qXX];
                  ewdp->D2phi[qYZ] = D2phi[qXY];
                  ewdp->D2phi[qZZ] = D2phi[qYY];

                  ewdp->D3phi[dXXX] = D3phi[dZZZ];
                  ewdp->D3phi[dXXY] = D3phi[dZZX];
                  ewdp->D3phi[dXXZ] = D3phi[dZZY];
                  ewdp->D3phi[dXYY] = D3phi[dZXX];
                  ewdp->D3phi[dXYZ] = D3phi[dZXY];
                  ewdp->D3phi[dXZZ] = D3phi[dZYY];
                  ewdp->D3phi[dYYY] = D3phi[dXXX];
                  ewdp->D3phi[dYYZ] = D3phi[dXXY];
                  ewdp->D3phi[dYZZ] = D3phi[dXYY];
                  ewdp->D3phi[dZZZ] = D3phi[dYYY];
                }
                break;

              case 1:
                {
                  ewdp->D0phi = ewald_D0(xx, zz, yy);
                  auto D1phi  = ewald_D1(xx, zz, yy);
                  auto D2phi  = ewald_D2(xx, zz, yy);
                  auto D3phi  = ewald_D3(xx, zz, yy);

                  ewdp->D1phi[vX] = D1phi[vX];
                  ewdp->D1phi[vY] = D1phi[vZ];
                  ewdp->D1phi[vZ] = D1phi[vY];

                  ewdp->D2phi[qXX] = D2phi[qXX];
                  ewdp->D2phi[qXY] = D2phi[qXZ];
                  ewdp->D2phi[qXZ] = D2phi[qXY];
                  ewdp->D2phi[qYY] = D2phi[qZZ];
                  ewdp->D2phi[qYZ] = D2phi[qZY];
                  ewdp->D2phi[qZZ] = D2phi[qYY];

                  ewdp->D3phi[dXXX] = D3phi[dXXX];
                  ewdp->D3phi[dXXY] = D3phi[dXXZ];
                  ewdp->D3phi[dXXZ] = D3phi[dXXY];
                  ewdp->D3phi[dXYY] = D3phi[dXZZ];
                  ewdp->D3phi[dXYZ] = D3phi[dXZY];
                  ewdp->D3phi[dXZZ] = D3phi[dXYY];
                  ewdp->D3phi[dYYY] = D3phi[dZZZ];
                  ewdp->D3phi[dYYZ] = D3phi[dZZY];
                  ewdp->D3phi[dYZZ] = D3phi[dZYY];
                  ewdp->D3phi[dZZZ] = D3phi[dYYY];
                }
                break;

              case 2:
                {
                  ewdp->D0phi = ewald_D0(xx, yy, zz);
                  ewdp->D1phi = ewald_D1(xx, yy, zz);
                  ewdp->D2phi = ewald_D2(xx, yy, zz);
                  ewdp->D3phi = ewald_D3(xx, yy, zz);
                }
                break;
            }
#endif
        }

      int *recvcnts = (int *)Mem.mymalloc("recvcnts", NTask * sizeof(int));
      int *recvoffs = (int *)Mem.mymalloc("recvoffs", NTask * sizeof(int));

      for(int i = 0; i < NTask; i++)
        {
          int off, cnt;
          subdivide_evenly(size, NTask, i, &off, &cnt);
          recvcnts[i] = cnt * sizeof(ewald_data);
          recvoffs[i] = off * sizeof(ewald_data);
        }

      myMPI_Allgatherv(MPI_IN_PLACE, size * sizeof(ewald_data), MPI_BYTE, Ewd, recvcnts, recvoffs, MPI_BYTE, Communicator);

      Mem.myfree(recvoffs);
      Mem.myfree(recvcnts);

      mpi_printf("\nEWALD: writing Ewald tables to file `%s'\n", buf);
      if(ThisTask == 0)
        {
          FILE *fd;
          if((fd = fopen(buf, "w")))
            {
              ewald_header tabh;
              tabh.resx    = ENX;
              tabh.resy    = ENY;
              tabh.resz    = ENZ;
              tabh.varsize = sizeof(MyFloat);
#ifndef GRAVITY_TALLBOX
              tabh.ewaldtype = -1;
#else
              tabh.ewaldtype = GRAVITY_TALLBOX + 1;
#endif
              my_fwrite(&tabh, sizeof(ewald_header), 1, fd);

              my_fwrite(Ewd, sizeof(ewald_data), (ENX + 1) * (ENY + 1) * (ENZ + 1), fd);
              fclose(fd);
            }
          else
            Terminate("can't write to file '%s'\n", buf);
        }
    }
  else
    {
      /* here we got them from disk */
      int len = (ENX + 1) * (ENY + 1) * (ENZ + 1) * sizeof(ewald_data);
      MPI_Bcast(Ewd, len, MPI_BYTE, 0, Communicator);
    }

  Ewd_fac_intp[0] = 2.0 * EN * LONG_X / All.BoxSize;
  Ewd_fac_intp[1] = 2.0 * EN * LONG_Y / All.BoxSize;
  Ewd_fac_intp[2] = 2.0 * EN * LONG_Z / All.BoxSize;

  /* now scale things to the boxsize that is actually used */
  for(int i = 0; i <= ENX; i++)
    for(int j = 0; j <= ENY; j++)
      for(int k = 0; k <= ENZ; k++)
        {
          ewald_data *ewdp = Ewd + ewd_offset(i, j, k);

          ewdp->D0phi *= 1 / All.BoxSize; /* potential */
          ewdp->D1phi *= 1 / pow(All.BoxSize, 2);
          ewdp->D2phi *= 1 / pow(All.BoxSize, 3);
          ewdp->D3phi *= 1 / pow(All.BoxSize, 4);
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 4
          ewdp->D4phi *= 1 / pow(All.BoxSize, 5);
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 5
          ewdp->D5phi *= 1 / pow(All.BoxSize, 6);
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 6
          ewdp->D6phi *= 1 / pow(All.BoxSize, 7);
#endif
#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 7
          ewdp->D7phi *= 1 / pow(All.BoxSize, 8);
#endif
        }

  mpi_printf("EWALD: Initialization of periodic boundaries finished.\n");

  ewald_is_initialized = 1;

  if(Shmem.Island_NTask != Shmem.World_NTask)
    {
      // We actually have multiple shared memory nodes in which we set aside one MPI rank for shared memory communication.
      // In this case, move the ewaldtable to the communication rank in order to consume this memory only once on the node

      if(Shmem.Island_ThisTask == 0)
        {
          size_t tab_len = sizeof(ewald_data) * (ENX + 1) * (ENY + 1) * (ENZ + 1);

          MPI_Send(&tab_len, sizeof(tab_len), MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_EWALD_ALLOC, MPI_COMM_WORLD);
          MPI_Send(Ewd, tab_len, MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_DMOM, MPI_COMM_WORLD);
        }

      Mem.myfree(Ewd);

      ptrdiff_t off;
      MPI_Bcast(&off, sizeof(ptrdiff_t), MPI_BYTE, Shmem.Island_NTask - 1, Shmem.SharedMemComm);

      Ewd = (ewald_data *)((char *)Shmem.SharedMemBaseAddr[Shmem.Island_NTask - 1] + off);
    }

#ifdef EWALD_TEST
  test_interpolation_accuracy();
#endif
}

void ewald::ewald_gridlookup(const MyIntPosType *p_intpos, const MyIntPosType *target_intpos, enum interpolate_options flag,
                             ewald_data &fper)
{
  // we determine the closest available point in our Ewald look-up table

  static MyIntPosType const halflenX   = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - (EWLEVEL + 1) - MAX_LONG_X_BITS);
  static MyIntPosType const intlenX    = halflenX << 1;
  static MyIntPosType const ewaldmaskX = ~(intlenX - 1);

  static MyIntPosType const halflenY   = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - (EWLEVEL + 1) - MAX_LONG_Y_BITS);
  static MyIntPosType const intlenY    = halflenY << 1;
  static MyIntPosType const ewaldmaskY = ~(intlenY - 1);

  static MyIntPosType const halflenZ   = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1) - (EWLEVEL + 1) - MAX_LONG_Z_BITS);
  static MyIntPosType const intlenZ    = halflenZ << 1;
  static MyIntPosType const ewaldmaskZ = ~(intlenZ - 1);

  MyIntPosType temppos[3] = {p_intpos[0] - target_intpos[0], p_intpos[1] - target_intpos[1], p_intpos[2] - target_intpos[2]};

  constrain_intpos(temppos);

  MyIntPosType gridpos[3];
  gridpos[0] = (temppos[0] + halflenX) & ewaldmaskX;
  gridpos[1] = (temppos[1] + halflenY) & ewaldmaskY;
  gridpos[2] = (temppos[2] + halflenZ) & ewaldmaskZ;

  vector<double> off;
  nearest_image_intpos_to_pos(temppos, gridpos, off.da);

  int i = (gridpos[0] >> (BITS_FOR_POSITIONS - (EWLEVEL + 1) - MAX_LONG_X_BITS));
  int j = (gridpos[1] >> (BITS_FOR_POSITIONS - (EWLEVEL + 1) - MAX_LONG_Y_BITS));
  int k = (gridpos[2] >> (BITS_FOR_POSITIONS - (EWLEVEL + 1) - MAX_LONG_Z_BITS));

  int signx = 1, signy = 1, signz = 1;

#if defined(GRAVITY_TALLBOX) && (GRAVITY_TALLBOX == 0)
  if(p_intpos[0] < target_intpos[0])
    {
      signx = -1;
      i     = ENX - i;
    }
#else
  if(i > ENX)
    {
      i = 2 * ENX - i;
      signx = -1;
    }
  else if(i == ENX && gridpos[0] < temppos[0])
    signx = -1;
#endif

#if defined(GRAVITY_TALLBOX) && (GRAVITY_TALLBOX == 1)
  if(p_intpos[1] < target_intpos[1])
    {
      signx = -1;
      j     = ENY - i;
    }
#else
  if(j > ENY)
    {
      j = 2 * ENY - j;
      signy = -1;
    }
  else if(j == ENY && gridpos[1] < temppos[1])
    signy = -1;
#endif

#if defined(GRAVITY_TALLBOX) && (GRAVITY_TALLBOX == 2)
  if(p_intpos[2] < target_intpos[2])
    {
      signz = -1;
      k     = ENZ - k;
    }
#else
  if(k > ENZ)
    {
      k = 2 * ENZ - k;
      signz = -1;
    }
  else if(k == ENZ && gridpos[2] < temppos[2])
    signz = -1;
#endif

  fper = Ewd[ewd_offset(i, j, k)];

  /* change signs as needed */

  fper.D1phi[0] *= signx;
  fper.D1phi[1] *= signy;
  fper.D1phi[2] *= signz;

  fper.D2phi[qXY] *= signx * signy;
  fper.D2phi[qXZ] *= signx * signz;
  fper.D2phi[qYZ] *= signy * signz;

  fper.D3phi[dXXX] *= signx;
  fper.D3phi[dXXY] *= signy;
  fper.D3phi[dXXZ] *= signz;
  fper.D3phi[dXYY] *= signx;
  fper.D3phi[dXYZ] *= signx * signy * signz;
  fper.D3phi[dXZZ] *= signx;
  fper.D3phi[dYYY] *= signy;
  fper.D3phi[dYYZ] *= signz;
  fper.D3phi[dYZZ] *= signy;
  fper.D3phi[dZZZ] *= signz;

#if EWALD_TAYLOR_ORDER == 3

  fper.D4phi[sXXXY] *= signx * signy;
  fper.D4phi[sXYYY] *= signx * signy;
  fper.D4phi[sXXXZ] *= signx * signz;
  fper.D4phi[sXZZZ] *= signx * signz;
  fper.D4phi[sYYYZ] *= signy * signz;
  fper.D4phi[sYZZZ] *= signy * signz;
  fper.D4phi[sXXYZ] *= signy * signz;
  fper.D4phi[sXYYZ] *= signx * signz;
  fper.D4phi[sXYZZ] *= signx * signy;

  // now Taylor corrections

  fper.D0phi += fper.D1phi * off + 0.5 * ((fper.D2phi * off) * off) + (1.0 / 6) * (((fper.D3phi * off) * off) * off);
  fper.D1phi += fper.D2phi * off + 0.5 * ((fper.D3phi * off) * off) + (1.0 / 6) * (((fper.D4phi * off) * off) * off);

  if(flag == POINTMASS)
    return;

#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 5
  fper.D5phi[rXXXXX] *= signx;
  fper.D5phi[rYYYYY] *= signy;
  fper.D5phi[rZZZZZ] *= signz;

  fper.D5phi[rXXXXY] *= signy;
  fper.D5phi[rXXXXZ] *= signz;
  fper.D5phi[rXYYYY] *= signx;
  fper.D5phi[rXZZZZ] *= signx;
  fper.D5phi[rYYYYZ] *= signz;
  fper.D5phi[rYZZZZ] *= signy;

  fper.D5phi[rXXXYY] *= signx;
  fper.D5phi[rXXXZZ] *= signx;
  fper.D5phi[rXXYYY] *= signy;
  fper.D5phi[rXXZZZ] *= signz;
  fper.D5phi[rYYYZZ] *= signy;
  fper.D5phi[rYYZZZ] *= signz;

  fper.D5phi[rXXYZZ] *= signy;
  fper.D5phi[rXXYYZ] *= signz;
  fper.D5phi[rXYYZZ] *= signx;

  fper.D5phi[rXXXYZ] *= signx * signy * signz;
  fper.D5phi[rXYYYZ] *= signx * signy * signz;
  fper.D5phi[rXYZZZ] *= signx * signy * signz;

  fper.D2phi += fper.D3phi * off + 0.5 * ((fper.D4phi * off) * off) + (1.0 / 6) * (((fper.D5phi * off) * off) * off);
#endif

#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 6
  fper.D6phi[pXXXXXY] *= signx * signy;
  fper.D6phi[pXXXXXZ] *= signx * signz;
  fper.D6phi[pXXXXYZ] *= signy * signz;
  fper.D6phi[pXXXYYY] *= signx * signy;
  fper.D6phi[pXXXYYZ] *= signx * signz;
  fper.D6phi[pXXXYZZ] *= signx * signy;
  fper.D6phi[pXXXZZZ] *= signx * signz;
  fper.D6phi[pXXYYYZ] *= signy * signz;
  fper.D6phi[pXXYZZZ] *= signy * signz;
  fper.D6phi[pXYYYYY] *= signx * signy;
  fper.D6phi[pXYYYYZ] *= signx * signz;
  fper.D6phi[pXYYYZZ] *= signx * signy;
  fper.D6phi[pXYYZZZ] *= signx * signz;
  fper.D6phi[pXYZZZZ] *= signx * signy;
  fper.D6phi[pXZZZZZ] *= signx * signz;
  fper.D6phi[pYYYYYZ] *= signy * signz;
  fper.D6phi[pYYYZZZ] *= signy * signz;
  fper.D6phi[pYZZZZZ] *= signy * signz;

  fper.D3phi += fper.D4phi * off + 0.5 * ((fper.D5phi * off) * off) + (1.0 / 6) * (((fper.D6phi * off) * off) * off);
#endif

#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 7
  fper.D7phi[tXXXXXXX] *= signx;
  fper.D7phi[tXXXXXXY] *= signy;
  fper.D7phi[tXXXXXXZ] *= signz;
  fper.D7phi[tXXXXXYY] *= signx;
  fper.D7phi[tXXXXXYZ] *= signx * signy * signz;
  fper.D7phi[tXXXXXZZ] *= signx;
  fper.D7phi[tXXXXYYY] *= signy;
  fper.D7phi[tXXXXYYZ] *= signz;
  fper.D7phi[tXXXXYZZ] *= signy;
  fper.D7phi[tXXXXZZZ] *= signz;
  fper.D7phi[tXXXYYYY] *= signx;
  fper.D7phi[tXXXYYYZ] *= signx * signy * signz;
  fper.D7phi[tXXXYYZZ] *= signx;
  fper.D7phi[tXXXYZZZ] *= signx * signy * signz;
  fper.D7phi[tXXXZZZZ] *= signx;
  fper.D7phi[tXXYYYYY] *= signy;
  fper.D7phi[tXXYYYYZ] *= signz;
  fper.D7phi[tXXYYYZZ] *= signy;
  fper.D7phi[tXXYYZZZ] *= signz;
  fper.D7phi[tXXYZZZZ] *= signy;
  fper.D7phi[tXXZZZZZ] *= signz;
  fper.D7phi[tXYYYYYY] *= signx;
  fper.D7phi[tXYYYYYZ] *= signx * signy * signz;
  fper.D7phi[tXYYYYZZ] *= signx;
  fper.D7phi[tXYYYZZZ] *= signx * signy * signz;
  fper.D7phi[tXYYZZZZ] *= signx;
  fper.D7phi[tXYZZZZZ] *= signx * signy * signz;
  fper.D7phi[tXZZZZZZ] *= signx;
  fper.D7phi[tYYYYYYY] *= signy;
  fper.D7phi[tYYYYYYZ] *= signz;
  fper.D7phi[tYYYYYZZ] *= signy;
  fper.D7phi[tYYYYZZZ] *= signz;
  fper.D7phi[tYYYZZZZ] *= signy;
  fper.D7phi[tYYZZZZZ] *= signz;
  fper.D7phi[tYZZZZZZ] *= signy;
  fper.D7phi[tZZZZZZZ] *= signz;

  fper.D4phi += fper.D5phi * off + 0.5 * ((fper.D6phi * off) * off) + (1.0 / 6) * (((fper.D7phi * off) * off) * off);
  fper.D5phi += fper.D6phi * off + 0.5 * ((fper.D7phi * off) * off);
#endif

#else

  // only second order Taylor expansion, i.e. EWALD_TAYLOR_ORDER==2

  // now Taylor corrections
  fper.D0phi += fper.D1phi * off + 0.5 * ((fper.D2phi * off) * off);
  fper.D1phi += fper.D2phi * off + 0.5 * ((fper.D3phi * off) * off);

  if(flag == POINTMASS)
    return;

#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 4
  fper.D4phi[sXXXY] *= signx * signy;
  fper.D4phi[sXYYY] *= signx * signy;
  fper.D4phi[sXXXZ] *= signx * signz;
  fper.D4phi[sXZZZ] *= signx * signz;
  fper.D4phi[sYYYZ] *= signy * signz;
  fper.D4phi[sYZZZ] *= signy * signz;
  fper.D4phi[sXXYZ] *= signy * signz;
  fper.D4phi[sXYYZ] *= signx * signz;
  fper.D4phi[sXYZZ] *= signx * signy;

  fper.D2phi += fper.D3phi * off + 0.5 * ((fper.D4phi * off) * off);
#endif

#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 5
  fper.D5phi[rXXXXX] *= signx;
  fper.D5phi[rYYYYY] *= signy;
  fper.D5phi[rZZZZZ] *= signz;

  fper.D5phi[rXXXXY] *= signy;
  fper.D5phi[rXXXXZ] *= signz;
  fper.D5phi[rXYYYY] *= signx;
  fper.D5phi[rXZZZZ] *= signx;
  fper.D5phi[rYYYYZ] *= signz;
  fper.D5phi[rYZZZZ] *= signy;

  fper.D5phi[rXXXYY] *= signx;
  fper.D5phi[rXXXZZ] *= signx;
  fper.D5phi[rXXYYY] *= signy;
  fper.D5phi[rXXZZZ] *= signz;
  fper.D5phi[rYYYZZ] *= signy;
  fper.D5phi[rYYZZZ] *= signz;

  fper.D5phi[rXXYZZ] *= signy;
  fper.D5phi[rXXYYZ] *= signz;
  fper.D5phi[rXYYZZ] *= signx;

  fper.D5phi[rXXXYZ] *= signx * signy * signz;
  fper.D5phi[rXYYYZ] *= signx * signy * signz;
  fper.D5phi[rXYZZZ] *= signx * signy * signz;

  fper.D3phi += fper.D4phi * off + 0.5 * ((fper.D5phi * off) * off);
#endif

#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 6
  fper.D6phi[pXXXXXY] *= signx * signy;
  fper.D6phi[pXXXXXZ] *= signx * signz;
  fper.D6phi[pXXXXYZ] *= signy * signz;
  fper.D6phi[pXXXYYY] *= signx * signy;
  fper.D6phi[pXXXYYZ] *= signx * signz;
  fper.D6phi[pXXXYZZ] *= signx * signy;
  fper.D6phi[pXXXZZZ] *= signx * signz;
  fper.D6phi[pXXYYYZ] *= signy * signz;
  fper.D6phi[pXXYZZZ] *= signy * signz;
  fper.D6phi[pXYYYYY] *= signx * signy;
  fper.D6phi[pXYYYYZ] *= signx * signz;
  fper.D6phi[pXYYYZZ] *= signx * signy;
  fper.D6phi[pXYYZZZ] *= signx * signz;
  fper.D6phi[pXYZZZZ] *= signx * signy;
  fper.D6phi[pXZZZZZ] *= signx * signz;
  fper.D6phi[pYYYYYZ] *= signy * signz;
  fper.D6phi[pYYYZZZ] *= signy * signz;
  fper.D6phi[pYZZZZZ] *= signy * signz;

  fper.D4phi += fper.D5phi * off + 0.5 * ((fper.D6phi * off) * off);
#endif

#if(HIGHEST_NEEDEDORDER_EWALD_DPHI + EWALD_TAYLOR_ORDER) >= 7
  fper.D7phi[tXXXXXXX] *= signx;
  fper.D7phi[tXXXXXXY] *= signy;
  fper.D7phi[tXXXXXXZ] *= signz;
  fper.D7phi[tXXXXXYY] *= signx;
  fper.D7phi[tXXXXXYZ] *= signx * signy * signz;
  fper.D7phi[tXXXXXZZ] *= signx;
  fper.D7phi[tXXXXYYY] *= signy;
  fper.D7phi[tXXXXYYZ] *= signz;
  fper.D7phi[tXXXXYZZ] *= signy;
  fper.D7phi[tXXXXZZZ] *= signz;
  fper.D7phi[tXXXYYYY] *= signx;
  fper.D7phi[tXXXYYYZ] *= signx * signy * signz;
  fper.D7phi[tXXXYYZZ] *= signx;
  fper.D7phi[tXXXYZZZ] *= signx * signy * signz;
  fper.D7phi[tXXXZZZZ] *= signx;
  fper.D7phi[tXXYYYYY] *= signy;
  fper.D7phi[tXXYYYYZ] *= signz;
  fper.D7phi[tXXYYYZZ] *= signy;
  fper.D7phi[tXXYYZZZ] *= signz;
  fper.D7phi[tXXYZZZZ] *= signy;
  fper.D7phi[tXXZZZZZ] *= signz;
  fper.D7phi[tXYYYYYY] *= signx;
  fper.D7phi[tXYYYYYZ] *= signx * signy * signz;
  fper.D7phi[tXYYYYZZ] *= signx;
  fper.D7phi[tXYYYZZZ] *= signx * signy * signz;
  fper.D7phi[tXYYZZZZ] *= signx;
  fper.D7phi[tXYZZZZZ] *= signx * signy * signz;
  fper.D7phi[tXZZZZZZ] *= signx;
  fper.D7phi[tYYYYYYY] *= signy;
  fper.D7phi[tYYYYYYZ] *= signz;
  fper.D7phi[tYYYYYZZ] *= signy;
  fper.D7phi[tYYYYZZZ] *= signz;
  fper.D7phi[tYYYZZZZ] *= signy;
  fper.D7phi[tYYZZZZZ] *= signz;
  fper.D7phi[tYZZZZZZ] *= signy;
  fper.D7phi[tZZZZZZZ] *= signz;

  fper.D5phi += fper.D6phi * off + 0.5 * ((fper.D7phi * off) * off);
#endif

#endif
}

/*! \brief This function computes the potential correction term by means of Ewald
 *  summation.
 *
 *  \param x, y, z contains the distance vector for which the correction
 *  term should be computed
 *  \return the correction term
 */
double ewald::ewald_D0(double x, double y, double z)
{
  static int printed = 0;

  double D0 = 0.0;

#ifndef GRAVITY_TALLBOX

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D0 table: qxmax=%d qymax=%d qzmax=%d   nxmax=%d nymax=%d nzmax=%d\n", qxmax, qymax, qzmax, nxmax, nymax,
                 nzmax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv = (r > 0) ? 1.0 / r : 0.0;

          double g0;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g0 = -erfc(alpha * r) * rinv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* for small r:
               *
               *   [1- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               */

              if((alpha * r) < 0.5)
                {
                  g0 = 2.0 * pow(alpha, 1) / sqrt(M_PI) *
                       (1.0 - pow(alpha * r, 2) / 3.0 + pow(alpha * r, 4) / 10.0 - pow(alpha * r, 6) / 42.0 +
                        pow(alpha * r, 8) / 216.0 - pow(alpha * r, 10) / 1320.0);
                }
              else
                {
                  g0 = erf(alpha * r) * rinv;
                }
            }

          D0 += g0;
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          if(nx != 0 || ny != 0 || nz != 0)
            {
              double kx    = (2.0 * M_PI * LONG_X) * nx;
              double ky    = (2.0 * M_PI * LONG_Y) * ny;
              double kz    = (2.0 * M_PI * LONG_Z) * nz;
              double k2    = kx * kx + ky * ky + kz * kz;
              double kdotx = (x * kx + y * ky + z * kz);

              D0 += -4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * cos(kdotx);
            }
        }

  D0 += M_PI * (LONG_X * LONG_Y * LONG_Z) / (alpha * alpha);

#else
  /* in the tallbox case, the third dimension, z, is assumed to be the non-periodic one */

  double leff = sqrt(BOXX * BOXY);
  double alpha = 2.0 / leff;

  int qxmax = (int)(8.0 / (BOXX * alpha) + 0.5);
  int qymax = (int)(8.0 / (BOXY * alpha) + 0.5);

  int nxmax = (int)(2.0 * alpha * BOXX + 0.5);
  int nymax = (int)(2.0 * alpha * BOXY + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D0 table: qxmax=%d qymax=%d   nxmax=%d nymax=%d\n", qxmax, qymax, nxmax, nymax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      {
        double dx = x - nx * BOXX;
        double dy = y - ny * BOXY;
        double r = sqrt(dx * dx + dy * dy + z * z);

        double rinv = (r > 0) ? 1.0 / r : 0.0;

        double g0;

        if(nx != 0 || ny != 0)
          {
            g0 = -erfc(alpha * r) * rinv;
          }
        else
          {
            /* we add the 1/r term here to the (0|0) entry */

            if((alpha * r) < 0.5)
              {
                g0 = 2.0 * pow(alpha, 1) / sqrt(M_PI) *
                     (1.0 - pow(alpha * r, 2) / 3.0 + pow(alpha * r, 4) / 10.0 - pow(alpha * r, 6) / 42.0 + pow(alpha * r, 8) / 216.0 -
                      pow(alpha * r, 10) / 1320.0);
              }
            else
              {
                g0 = erf(alpha * r) * rinv;
              }
          }

        D0 += g0;
      }

  double alpha2 = alpha * alpha;

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double kx = (2.0 * M_PI / BOXX) * nx;
            double ky = (2.0 * M_PI / BOXY) * ny;
            double k2 = kx * kx + ky * ky;
            double k = sqrt(k2);

            double ex = exp(-k * z);  // note: z positive here

            if(ex > 1.0e-60)  // to prevent divisions by zero due to underflows
              D0 += -M_PI / (BOXX * BOXY) * cos(kx * x + ky * y) / k * (specerf(z, k, alpha) + specerf(-z, k, alpha));
          }
      }

  D0 += 2.0 * alpha / sqrt(M_PI) + 2 * sqrt(M_PI) / (BOXX * BOXY) * (exp(-alpha2 * z * z) / alpha + sqrt(M_PI) * z * erf(alpha * z));

#endif

  return D0;
}

double ewald::specerf(double z, double k, double alpha) { return exp(k * z) * erfc(k / (2 * alpha) + alpha * z); }

double ewald::d_specerf(double z, double k, double alpha)
{
  return -2 * alpha / (sqrt(M_PI) * exp(pow(k / (2 * alpha), 2) + pow(alpha * z, 2))) + k * specerf(z, k, alpha);
}

double ewald::dd_specerf(double z, double k, double alpha)
{
  return +4 * pow(alpha, 3) * z / (sqrt(M_PI) * exp(pow(k / (2 * alpha), 2) + pow(alpha * z, 2))) + k * d_specerf(z, k, alpha);
}

double ewald::ddd_specerf(double z, double k, double alpha)
{
  return +4 * pow(alpha, 3) / (sqrt(M_PI) * exp(pow(k / (2 * alpha), 2) + pow(alpha * z, 2))) -
         8 * pow(alpha, 5) * z * z / (sqrt(M_PI) * exp(pow(k / (2 * alpha), 2) + pow(alpha * z, 2))) + k * dd_specerf(z, k, alpha);
}

/*! \brief This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 *
 *  \param x, y, z contains the distance vector for which the correction
 *  force should be computed
 *  \param force  array will containing the correction force
 */
vector<double> ewald::ewald_D1(double x, double y, double z)
{
  static int printed = 0;

  vector<double> D1 = 0.0;

#ifndef GRAVITY_TALLBOX

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D1 table: qxmax=%d qymax=%d qzmax=%d   nxmax=%d nymax=%d nzmax=%d\n", qxmax, qymax, qzmax, nxmax, nymax,
                 nzmax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv  = (r > 0) ? 1.0 / r : 0.0;
          double r2inv = rinv * rinv;
          double r3inv = r2inv * rinv;

          double g1;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g1 = (erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* Note, for small r:
               *
               *   [1/- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               *
               *   Hence for r = 0:
               *
               *   g0 =  2     * alpha   / sqrt(pi)
               *   g1 = -4/3   * alpha^3 / sqrt(pi)
               *   g2 =  8/5   * alpha^5 / sqrt(pi)
               *   g3 = -16/7  * alpha^7 / sqrt(pi)
               *   g4 =  32/9  * alpha^9 / sqrt(pi)
               *   g5 = -64/11 * alpha^11/ sqrt(pi)
               */

              if((alpha * r) < 0.5)
                {
                  g1 = 4.0 * pow(alpha, 3) / sqrt(M_PI) *
                       (-1.0 / 3.0 + pow(alpha * r, 2) / 5.0 - pow(alpha * r, 4) / 14.0 + pow(alpha * r, 6) / 54.0 -
                        pow(alpha * r, 8) / 264.0 + pow(alpha * r, 10) / 1560.0);
                }
              else
                {
                  g1 = (-erf(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;
                }
            }

          D1 += g1 * dxyz;
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          double kx = (2.0 * M_PI * LONG_X) * nx;
          double ky = (2.0 * M_PI * LONG_Y) * ny;
          double kz = (2.0 * M_PI * LONG_Z) * nz;
          double k2 = kx * kx + ky * ky + kz * kz;

          if(k2 > 0)
            {
              double kdotx = (x * kx + y * ky + z * kz);
              double val   = 4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * sin(kdotx);
              D1[0] += kx * val;
              D1[1] += ky * val;
              D1[2] += kz * val;
            }
        }
#else
  /* this is the case with periodicity only in two dimensions */

  double leff = sqrt(BOXX * BOXY);
  double alpha = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 / (BOXX * alpha) + 0.5);
  int qymax = (int)(8.0 / (BOXY * alpha) + 0.5);

  int nxmax = (int)(2.0 * alpha * BOXX + 0.5);
  int nymax = (int)(2.0 * alpha * BOXY + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D1 table: qxmax=%d qymax=%d    nxmax=%d nymax=%d\n", qxmax, qymax, nxmax, nymax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      {
        double dx = x - nx * BOXX;
        double dy = y - ny * BOXY;
        double dz = z;

        vector<double> dxyz(dx, dy, dz);

        double r2 = dx * dx + dy * dy + dz * dz;
        double r = sqrt(r2);

        double rinv = (r > 0) ? 1.0 / r : 0.0;
        double r2inv = rinv * rinv;
        double r3inv = r2inv * rinv;

        double g1;

        if(nx != 0 || ny != 0)
          {
            g1 = (erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;
          }
        else
          {
            /* we add the 1/r term here to the (0|0) entry */

            if((alpha * r) < 0.5)
              {
                g1 = 4.0 * pow(alpha, 3) / sqrt(M_PI) *
                     (-1.0 / 3.0 + pow(alpha * r, 2) / 5.0 - pow(alpha * r, 4) / 14.0 + pow(alpha * r, 6) / 54.0 -
                      pow(alpha * r, 8) / 264.0 + pow(alpha * r, 10) / 1560.0);
              }
            else
              {
                g1 = (-erf(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;
              }
          }

        D1 += g1 * dxyz;
      }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double kx = (2.0 * M_PI / BOXX) * nx;
            double ky = (2.0 * M_PI / BOXY) * ny;
            double k2 = kx * kx + ky * ky;
            double k = sqrt(k2);

            double ex = exp(-k * z);  // note: z positive here

            if(ex > 1.0e-60)  // to prevent divisions by zero due to underflows
              {
                double val = M_PI / (BOXX * BOXY) / k * (specerf(z, k, alpha) + specerf(-z, k, alpha));

                D1[0] += kx * val * sin(kx * x + ky * y);
                D1[1] += ky * val * sin(kx * x + ky * y);
                D1[2] += -M_PI / (BOXX * BOXY) * cos(kx * x + ky * y) / k * (d_specerf(z, k, alpha) - d_specerf(-z, k, alpha));
              }
          }
      }

  D1[2] += 2.0 * M_PI / (BOXX * BOXY) * erf(alpha * z);
#endif

  return D1;  // now in dimensionless form;
}

symtensor2<double> ewald::ewald_D2(double x, double y, double z)
{
  static int printed = 0;

  symtensor2<double> D2 = 0.0;

#ifndef GRAVITY_TALLBOX

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D2 table: qxmax=%d qymax=%d qzmax=%d   nxmax=%d nymax=%d nzmax=%d\n", qxmax, qymax, qzmax, nxmax, nymax,
                 nzmax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv  = (r > 0) ? 1.0 / r : 0.0;
          double r2inv = rinv * rinv;
          double r3inv = r2inv * rinv;
          double r5inv = r3inv * r2inv;

          double g1, g2;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g1 = (erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;

              g2 = -(3.0 * erfc(alpha * r) + (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* Note, for small r:
               *
               *   [1/- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               *
               *   Hence for r = 0:
               *
               *   g0 =  2     * alpha   / sqrt(pi)
               *   g1 = -4/3   * alpha^3 / sqrt(pi)
               *   g2 =  8/5   * alpha^5 / sqrt(pi)
               *   g3 = -16/7  * alpha^7 / sqrt(pi)
               *   g4 =  32/9  * alpha^9 / sqrt(pi)
               *   g5 = -64/11 * alpha^11/ sqrt(pi)
               */

              if((alpha * r) < 0.5)
                {
                  g1 = 4.0 * pow(alpha, 3) / sqrt(M_PI) *
                       (-1.0 / 3.0 + pow(alpha * r, 2) / 5.0 - pow(alpha * r, 4) / 14.0 + pow(alpha * r, 6) / 54.0 -
                        pow(alpha * r, 8) / 264.0 + pow(alpha * r, 10) / 1560.0);

                  g2 = 8.0 * pow(alpha, 5) / sqrt(M_PI) *
                       (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                        pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);
                }
              else
                {
                  g1 = (-erf(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;

                  g2 = (3.0 * erf(alpha * r) - (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;
                }
            }

          D2 += g2 * (dxyz % dxyz);
          D2[qXX] += g1;
          D2[qYY] += g1;
          D2[qZZ] += g1;
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          if(nx != 0 || ny != 0 || nz != 0)
            {
              double kx = (2.0 * M_PI * LONG_X) * nx;
              double ky = (2.0 * M_PI * LONG_Y) * ny;
              double kz = (2.0 * M_PI * LONG_Z) * nz;
              double k2 = kx * kx + ky * ky + kz * kz;

              double kdotx = (x * kx + y * ky + z * kz);
              double val   = 4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * cos(kdotx);

              vector<double> kxyz(kx, ky, kz);

              D2 += (val * kxyz) % kxyz;
            }
        }
#else
  /* this is the case with periodicity only in two dimensions */
  /* this is the case with periodicity only in two dimensions */

  double leff = sqrt(BOXX * BOXY);
  double alpha = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 / (BOXX * alpha) + 0.5);
  int qymax = (int)(8.0 / (BOXY * alpha) + 0.5);

  int nxmax = (int)(2.0 * alpha * BOXX + 0.5);
  int nymax = (int)(2.0 * alpha * BOXY + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D2 table: qxmax=%d qymax=%d    nxmax=%d nymax=%d\n", qxmax, qymax, nxmax, nymax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      {
        double dx = x - nx * BOXX;
        double dy = y - ny * BOXY;
        double dz = z;

        vector<double> dxyz(dx, dy, dz);

        double r2 = dx * dx + dy * dy + dz * dz;
        double r = sqrt(r2);

        double rinv = (r > 0) ? 1.0 / r : 0.0;
        double r2inv = rinv * rinv;
        double r3inv = r2inv * rinv;
        double r5inv = r3inv * r2inv;

        double g1, g2;

        if(nx != 0 || ny != 0)
          {
            g1 = (erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;

            g2 = -(3.0 * erfc(alpha * r) + (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;
          }
        else
          {
            /* we add the 1/r term here to the (0|0) entry */

            if((alpha * r) < 0.5)
              {
                g1 = 4.0 * pow(alpha, 3) / sqrt(M_PI) *
                     (-1.0 / 3.0 + pow(alpha * r, 2) / 5.0 - pow(alpha * r, 4) / 14.0 + pow(alpha * r, 6) / 54.0 -
                      pow(alpha * r, 8) / 264.0 + pow(alpha * r, 10) / 1560.0);

                g2 = 8.0 * pow(alpha, 5) / sqrt(M_PI) *
                     (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                      pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);
              }
            else
              {
                g1 = (-erf(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;

                g2 = (3.0 * erf(alpha * r) - (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;
              }
          }

        D2 += g2 * (dxyz % dxyz);
        D2[qXX] += g1;
        D2[qYY] += g1;
        D2[qZZ] += g1;
      }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double kx = (2.0 * M_PI / BOXX) * nx;
            double ky = (2.0 * M_PI / BOXY) * ny;
            double k2 = kx * kx + ky * ky;
            double k = sqrt(k2);

            double ex = exp(-k * z);  // note: z positive here

            if(ex > 1.0e-60)  // to prevent divisions by zero due to underflows
              {
                double val = M_PI / (BOXX * BOXY) / k * (specerf(z, k, alpha) + specerf(-z, k, alpha));
                double dzval = M_PI / (BOXX * BOXY) / k * (d_specerf(z, k, alpha) - d_specerf(-z, k, alpha));

                D2[qXX] += kx * kx * val * cos(kx * x + ky * y);
                D2[qXY] += kx * ky * val * cos(kx * x + ky * y);
                D2[qXZ] += kx * dzval * sin(kx * x + ky * y);
                D2[qYY] += ky * ky * val * cos(kx * x + ky * y);
                D2[qYZ] += ky * dzval * sin(kx * x + ky * y);
                D2[qZZ] += -M_PI / (BOXX * BOXY) * cos(kx * x + ky * y) / k * (dd_specerf(z, k, alpha) + dd_specerf(-z, k, alpha));
              }
          }
      }

  D2[qZZ] += 4.0 * alpha * sqrt(M_PI) / (BOXX * BOXY) * exp(-pow(alpha * z, 2));
#endif

  return D2;
}

symtensor3<double> ewald::ewald_D3(double x, double y, double z)
{
  static int printed = 0;

  symtensor3<double> D3 = 0.0;

#ifndef GRAVITY_TALLBOX

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D3 table: qxmax=%d qymax=%d qzmax=%d   nxmax=%d nymax=%d nzmax=%d\n", qxmax, qymax, qzmax, nxmax, nymax,
                 nzmax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv  = (r > 0) ? 1.0 / r : 0.0;
          double r2inv = rinv * rinv;
          double r3inv = r2inv * rinv;
          double r4inv = r2inv * r2inv;
          double r5inv = r2inv * r3inv;
          double r7inv = r3inv * r4inv;

          double g2, g3;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g2 = -(3.0 * erfc(alpha * r) + (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;

              g3 = (15.0 * erfc(alpha * r) +
                    (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r7inv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* Note, for small r:
               *
               *   [1/- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               *
               *   Hence for r = 0:
               *
               *   g0 =  2     * alpha   / sqrt(pi)
               *   g1 = -4/3   * alpha^3 / sqrt(pi)
               *   g2 =  8/5   * alpha^5 / sqrt(pi)
               *   g3 = -16/7  * alpha^7 / sqrt(pi)
               *   g4 =  32/9  * alpha^9 / sqrt(pi)
               *   g5 = -64/11 * alpha^11/ sqrt(pi)
               */
              if((alpha * r) < 0.5)
                {
                  g2 = 8.0 * pow(alpha, 5) / sqrt(M_PI) *
                       (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                        pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);

                  g3 = 16.0 * pow(alpha, 7) / sqrt(M_PI) *
                       (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                        pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);
                }
              else
                {
                  g2 = (3.0 * erf(alpha * r) - (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;

                  g3 = (-15.0 * erf(alpha * r) +
                        (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r7inv;
                }
            }

          symtensor2<double> aux2 = dxyz % dxyz;
          symtensor3<double> aux3;

          setup_D3(ADD, D3, dxyz, aux2, aux3, g2, g3);
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          if(nx != 0 || ny != 0 || nz != 0)
            {
              double kx = (2.0 * M_PI * LONG_X) * nx;
              double ky = (2.0 * M_PI * LONG_Y) * ny;
              double kz = (2.0 * M_PI * LONG_Z) * nz;
              double k2 = kx * kx + ky * ky + kz * kz;

              double kdotx = (x * kx + y * ky + z * kz);
              double val   = -4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * sin(kdotx);

              vector<double> kxyz(kx, ky, kz);

              D3 += (val * kxyz) % (kxyz % kxyz);
            }
        }
#else
  /* this is the case with periodicity only in two dimensions */
  /* this is the case with periodicity only in two dimensions */
  /* this is the case with periodicity only in two dimensions */

  double leff = sqrt(BOXX * BOXY);
  double alpha = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 / (BOXX * alpha) + 0.5);
  int qymax = (int)(8.0 / (BOXY * alpha) + 0.5);

  int nxmax = (int)(2.0 * alpha * BOXX + 0.5);
  int nymax = (int)(2.0 * alpha * BOXY + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D2 table: qxmax=%d qymax=%d    nxmax=%d nymax=%d\n", qxmax, qymax, nxmax, nymax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      {
        double dx = x - nx * BOXX;
        double dy = y - ny * BOXY;
        double dz = z;

        vector<double> dxyz(dx, dy, dz);

        double r2 = dx * dx + dy * dy + dz * dz;
        double r = sqrt(r2);

        double rinv = (r > 0) ? 1.0 / r : 0.0;
        double r2inv = rinv * rinv;
        double r3inv = r2inv * rinv;
        double r4inv = r2inv * r2inv;
        double r5inv = r2inv * r3inv;
        double r7inv = r3inv * r4inv;

        double g2, g3;

        if(nx != 0 || ny != 0)
          {
            g2 = -(3.0 * erfc(alpha * r) + (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;

            g3 = (15.0 * erfc(alpha * r) +
                  (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                 r7inv;
          }
        else
          {
            if((alpha * r) < 0.5)
              {
                g2 = 8.0 * pow(alpha, 5) / sqrt(M_PI) *
                     (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                      pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);

                g3 = 16.0 * pow(alpha, 7) / sqrt(M_PI) *
                     (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                      pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);
              }
            else
              {
                g2 = (3.0 * erf(alpha * r) - (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;

                g3 = (-15.0 * erf(alpha * r) +
                      (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                     r7inv;
              }
          }

        symtensor2<double> aux2 = dxyz % dxyz;
        symtensor3<double> aux3;

        setup_D3(ADD, D3, dxyz, aux2, aux3, g2, g3);
      }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      {
        if(nx != 0 || ny != 0)
          {
            double kx = (2.0 * M_PI / BOXX) * nx;
            double ky = (2.0 * M_PI / BOXY) * ny;
            double k2 = kx * kx + ky * ky;
            double k = sqrt(k2);

            double ex = exp(-k * z);  // note: z positive here

            if(ex > 1.0e-60)  // to prevent divisions by zero due to underflows
              {
                double val = M_PI / (BOXX * BOXY) / k * (specerf(z, k, alpha) + specerf(-z, k, alpha));
                double dzval = M_PI / (BOXX * BOXY) / k * (d_specerf(z, k, alpha) - d_specerf(-z, k, alpha));
                double dzdzval = M_PI / (BOXX * BOXY) / k * (dd_specerf(z, k, alpha) + dd_specerf(-z, k, alpha));

                D3[dXXX] += -kx * kx * kx * val * sin(kx * x + ky * y);
                D3[dXXY] += -kx * kx * ky * val * sin(kx * x + ky * y);
                D3[dXXZ] += kx * kx * dzval * cos(kx * x + ky * y);
                D3[dXYY] += -kx * ky * ky * val * sin(kx * x + ky * y);
                D3[dXYZ] += kx * ky * dzval * cos(kx * x + ky * y);
                D3[dXZZ] += kx * dzdzval * sin(kx * x + ky * y);
                D3[dYYY] += -ky * ky * ky * val * sin(kx * x + ky * y);
                D3[dYYZ] += ky * ky * dzval * cos(kx * x + ky * y);
                D3[dYZZ] += ky * dzdzval * sin(kx * x + ky * y);
                D3[dZZZ] += -M_PI / (BOXX * BOXY) * cos(kx * x + ky * y) / k * (ddd_specerf(z, k, alpha) - ddd_specerf(-z, k, alpha));
              }
          }
      }

  D3[dZZZ] += -8.0 * pow(alpha, 3) * z * sqrt(M_PI) / (BOXX * BOXY) * exp(-pow(alpha * z, 2));

#endif

  return D3;
}

symtensor4<double> ewald::ewald_D4(double x, double y, double z)
{
  static int printed = 0;

#ifdef GRAVITY_TALLBOX
  Terminate("GRAVITY_TALLBOX is not implemented for MULTIPOLE_ORDER >= 4");
#endif

  symtensor4<double> D4 = 0.0;

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D4 table: qxmax=%d qymax=%d qzmax=%d   nxmax=%d nymax=%d nzmax=%d\n", qxmax, qymax, qzmax, nxmax, nymax,
                 nzmax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv  = (r > 0) ? 1.0 / r : 0.0;
          double r2inv = rinv * rinv;
          double r3inv = r2inv * rinv;
          double r4inv = r2inv * r2inv;
          double r5inv = r2inv * r3inv;
          double r7inv = r3inv * r4inv;
          double r9inv = r4inv * r5inv;

          double g2, g3, g4;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g2 = -(3.0 * erfc(alpha * r) + (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;

              g3 = (15.0 * erfc(alpha * r) +
                    (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r7inv;

              g4 = -(105.0 * erfc(alpha * r) +
                     (210.0 * alpha * r + 140.0 * pow(alpha * r, 3) + 56.0 * pow(alpha * r, 5) + 16.0 * pow(alpha * r, 7)) /
                         sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r9inv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* Note, for small r:
               *
               *   [1/- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               *
               *   Hence for r = 0:
               *
               *   g0 =  2     * alpha   / sqrt(pi)
               *   g1 = -4/3   * alpha^3 / sqrt(pi)
               *   g2 =  8/5   * alpha^5 / sqrt(pi)
               *   g3 = -16/7  * alpha^7 / sqrt(pi)
               *   g4 =  32/9  * alpha^9 / sqrt(pi)
               *   g5 = -64/11 * alpha^11/ sqrt(pi)
               */

              if((alpha * r) < 0.5)
                {
                  g2 = 8.0 * pow(alpha, 5) / sqrt(M_PI) *
                       (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                        pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);

                  g3 = 16.0 * pow(alpha, 7) / sqrt(M_PI) *
                       (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                        pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);

                  g4 = 32.0 * pow(alpha, 9) / sqrt(M_PI) *
                       (1.0 / 9.0 - pow(alpha * r, 2) / 11.0 + pow(alpha * r, 4) / 26.0 - pow(alpha * r, 6) / 90.0 +
                        pow(alpha * r, 8) / 408.0 - pow(alpha * r, 10) / 2280.0);
                }
              else
                {
                  g2 = (3.0 * erf(alpha * r) - (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;

                  g3 = (-15.0 * erf(alpha * r) +
                        (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r7inv;

                  g4 = (105.0 * erf(alpha * r) -
                        (210.0 * alpha * r + 140.0 * pow(alpha * r, 3) + 56.0 * pow(alpha * r, 5) + 16.0 * pow(alpha * r, 7)) /
                            sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r9inv;
                }
            }

          symtensor2<double> aux2 = dxyz % dxyz;
          symtensor3<double> aux3 = dxyz % aux2;
          symtensor4<double> aux4;

          setup_D4(ADD, D4, dxyz, aux2, aux3, aux4, g2, g3, g4);
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          if(nx != 0 || ny != 0 || nz != 0)
            {
              double kx = (2.0 * M_PI * LONG_X) * nx;
              double ky = (2.0 * M_PI * LONG_Y) * ny;
              double kz = (2.0 * M_PI * LONG_Z) * nz;
              double k2 = kx * kx + ky * ky + kz * kz;

              double kdotx = (x * kx + y * ky + z * kz);
              double val   = -4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * cos(kdotx);

              vector<double> kxyz(kx, ky, kz);

              D4 += (val * kxyz) % ((kxyz % (kxyz % kxyz)));
            }
        }

  return D4;
}

symtensor5<double> ewald::ewald_D5(double x, double y, double z)
{
  static int printed = 0;

#ifdef GRAVITY_TALLBOX
  Terminate("GRAVITY_TALLBOX is not implemented for MULTIPOLE_ORDER >= 4");
#endif

  symtensor5<double> D5 = 0.0;

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D5 table: qxmax=%d qymax=%d qzmax=%d   nxmax=%d nymax=%d nzmax=%d\n", qxmax, qymax, qzmax, nxmax, nymax,
                 nzmax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv   = (r > 0) ? 1.0 / r : 0.0;
          double r2inv  = rinv * rinv;
          double r3inv  = r2inv * rinv;
          double r4inv  = r2inv * r2inv;
          double r5inv  = r2inv * r3inv;
          double r7inv  = r3inv * r4inv;
          double r9inv  = r4inv * r5inv;
          double r11inv = r4inv * r7inv;

          double g3, g4, g5;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g3 = (15.0 * erfc(alpha * r) +
                    (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r7inv;

              g4 = -(105.0 * erfc(alpha * r) +
                     (210.0 * alpha * r + 140.0 * pow(alpha * r, 3) + 56.0 * pow(alpha * r, 5) + 16.0 * pow(alpha * r, 7)) /
                         sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r9inv;

              g5 = (945.0 * erfc(alpha * r) + (1890.0 * alpha * r + 1260.0 * pow(alpha * r, 3) + 504.0 * pow(alpha * r, 5) +
                                               144.0 * pow(alpha * r, 7) + 32.0 * pow(alpha * r, 9)) /
                                                  sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r11inv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* Note, for small r:
               *
               *   [1/- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               *
               *   Hence for r = 0:
               *
               *   g0 =  2     * alpha   / sqrt(pi)
               *   g1 = -4/3   * alpha^3 / sqrt(pi)
               *   g2 =  8/5   * alpha^5 / sqrt(pi)
               *   g3 = -16/7  * alpha^7 / sqrt(pi)
               *   g4 =  32/9  * alpha^9 / sqrt(pi)
               *   g5 = -64/11 * alpha^11/ sqrt(pi)
               */

              if((alpha * r) < 0.5)
                {
                  g3 = 16.0 * pow(alpha, 7) / sqrt(M_PI) *
                       (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                        pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);

                  g4 = 32.0 * pow(alpha, 9) / sqrt(M_PI) *
                       (1.0 / 9.0 - pow(alpha * r, 2) / 11.0 + pow(alpha * r, 4) / 26.0 - pow(alpha * r, 6) / 90.0 +
                        pow(alpha * r, 8) / 408.0 - pow(alpha * r, 10) / 2280.0);

                  g5 = 64.0 * pow(alpha, 11) / sqrt(M_PI) *
                       (-1.0 / 11.0 + pow(alpha * r, 2) / 13.0 - pow(alpha * r, 4) / 30.0 + pow(alpha * r, 6) / 102.0 -
                        pow(alpha * r, 8) / 456.0 + pow(alpha * r, 10) / 2520.0);
                }
              else
                {
                  g3 = (-15.0 * erf(alpha * r) +
                        (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r7inv;

                  g4 = (105.0 * erf(alpha * r) -
                        (210.0 * alpha * r + 140.0 * pow(alpha * r, 3) + 56.0 * pow(alpha * r, 5) + 16.0 * pow(alpha * r, 7)) /
                            sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r9inv;

                  g5 = (-945.0 * erf(alpha * r) + (1890.0 * alpha * r + 1260.0 * pow(alpha * r, 3) + 504.0 * pow(alpha * r, 5) +
                                                   144.0 * pow(alpha * r, 7) + 32.0 * pow(alpha * r, 9)) /
                                                      sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r11inv;
                }
            }

          symtensor3<double> aux3 = dxyz % (dxyz % dxyz);
          symtensor4<double> aux4 = dxyz % aux3;
          symtensor5<double> aux5;

          setup_D5(ADD, D5, dxyz, aux3, aux4, aux5, g3, g4, g5);
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          double kx = (2.0 * M_PI * LONG_X) * nx;
          double ky = (2.0 * M_PI * LONG_Y) * ny;
          double kz = (2.0 * M_PI * LONG_Z) * nz;
          double k2 = kx * kx + ky * ky + kz * kz;

          if(k2 > 0)
            {
              double kdotx = (x * kx + y * ky + z * kz);
              double val   = 4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * sin(kdotx);

              vector<double> kxyz(kx, ky, kz);

              D5 += (val * kxyz) % (kxyz % ((kxyz % (kxyz % kxyz))));
            }
        }

  return D5;
}

symtensor6<double> ewald::ewald_D6(double x, double y, double z)
{
  static int printed = 0;

#ifdef GRAVITY_TALLBOX
  Terminate("GRAVITY_TALLBOX is not implemented for MULTIPOLE_ORDER >= 4");
#endif

  symtensor6<double> D6 = 0.0;

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D6 table: qxmax=%d qymax=%d qzmax=%d   nxmax=%d nymax=%d nzmax=%d\n", qxmax, qymax, qzmax, nxmax, nymax,
                 nzmax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv   = (r > 0) ? 1.0 / r : 0.0;
          double r2inv  = rinv * rinv;
          double r3inv  = r2inv * rinv;
          double r4inv  = r2inv * r2inv;
          double r5inv  = r2inv * r3inv;
          double r7inv  = r3inv * r4inv;
          double r9inv  = r4inv * r5inv;
          double r11inv = r4inv * r7inv;
          double r13inv = r4inv * r9inv;

          double g3, g4, g5, g6;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g3 = (15.0 * erfc(alpha * r) +
                    (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r7inv;

              g4 = -(105.0 * erfc(alpha * r) +
                     (210.0 * alpha * r + 140.0 * pow(alpha * r, 3) + 56.0 * pow(alpha * r, 5) + 16.0 * pow(alpha * r, 7)) /
                         sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r9inv;

              g5 = (945.0 * erfc(alpha * r) + (1890.0 * alpha * r + 1260.0 * pow(alpha * r, 3) + 504.0 * pow(alpha * r, 5) +
                                               144.0 * pow(alpha * r, 7) + 32.0 * pow(alpha * r, 9)) /
                                                  sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r11inv;

              g6 = (-10395.0 * erfc(alpha * r) -
                    2.0 *
                        (10395.0 * alpha * r + 6930.0 * pow(alpha * r, 3) + 2772.0 * pow(alpha * r, 5) + 792.0 * pow(alpha * r, 7) +
                         176.0 * pow(alpha * r, 9) + 32.0 * pow(alpha * r, 11)) /
                        sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r13inv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* Note, for small r:
               *
               *   [1/- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               *
               *   Hence for r = 0:
               *
               *   g0 =  2     * alpha   / sqrt(pi)
               *   g1 = -4/3   * alpha^3 / sqrt(pi)
               *   g2 =  8/5   * alpha^5 / sqrt(pi)
               *   g3 = -16/7  * alpha^7 / sqrt(pi)
               *   g4 =  32/9  * alpha^9 / sqrt(pi)
               *   g5 = -64/11 * alpha^11/ sqrt(pi)
               */

              if((alpha * r) < 0.5)
                {
                  g3 = 16.0 * pow(alpha, 7) / sqrt(M_PI) *
                       (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                        pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);

                  g4 = 32.0 * pow(alpha, 9) / sqrt(M_PI) *
                       (1.0 / 9.0 - pow(alpha * r, 2) / 11.0 + pow(alpha * r, 4) / 26.0 - pow(alpha * r, 6) / 90.0 +
                        pow(alpha * r, 8) / 408.0 - pow(alpha * r, 10) / 2280.0);

                  g5 = 64.0 * pow(alpha, 11) / sqrt(M_PI) *
                       (-1.0 / 11.0 + pow(alpha * r, 2) / 13.0 - pow(alpha * r, 4) / 30.0 + pow(alpha * r, 6) / 102.0 -
                        pow(alpha * r, 8) / 456.0 + pow(alpha * r, 10) / 2520.0);

                  g6 = 128.0 * pow(alpha, 13) / sqrt(M_PI) *
                       (1.0 / 13.0 - pow(alpha * r, 2) / 15.0 + pow(alpha * r, 4) / 34.0 - pow(alpha * r, 6) / 114.0 +
                        pow(alpha * r, 8) / 504.0 - pow(alpha * r, 10) / 2760.0);
                }
              else
                {
                  g3 = (-15.0 * erf(alpha * r) +
                        (30.0 * alpha * r + 20.0 * pow(alpha * r, 3) + 8.0 * pow(alpha * r, 5)) / sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r7inv;

                  g4 = (105.0 * erf(alpha * r) -
                        (210.0 * alpha * r + 140.0 * pow(alpha * r, 3) + 56.0 * pow(alpha * r, 5) + 16.0 * pow(alpha * r, 7)) /
                            sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r9inv;

                  g5 = (-945.0 * erf(alpha * r) + (1890.0 * alpha * r + 1260.0 * pow(alpha * r, 3) + 504.0 * pow(alpha * r, 5) +
                                                   144.0 * pow(alpha * r, 7) + 32.0 * pow(alpha * r, 9)) /
                                                      sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r11inv;

                  g6 = (10395.0 * erf(alpha * r) -
                        2.0 *
                            (10395.0 * alpha * r + 6930.0 * pow(alpha * r, 3) + 2772.0 * pow(alpha * r, 5) +
                             792.0 * pow(alpha * r, 7) + 176.0 * pow(alpha * r, 9) + 32.0 * pow(alpha * r, 11)) /
                            sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r13inv;
                }
            }

          setup_D6(ADD, D6, dxyz, g3, g4, g5, g6);
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          double kx = (2.0 * M_PI * LONG_X) * nx;
          double ky = (2.0 * M_PI * LONG_Y) * ny;
          double kz = (2.0 * M_PI * LONG_Z) * nz;
          double k2 = kx * kx + ky * ky + kz * kz;

          if(k2 > 0)
            {
              double kdotx = (x * kx + y * ky + z * kz);
              double val   = 4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * cos(kdotx);

              vector<double> kxyz(kx, ky, kz);

              D6 += (val * kxyz) % (kxyz % (kxyz % ((kxyz % (kxyz % kxyz)))));
            }
        }

  return D6;
}

symtensor7<double> ewald::ewald_D7(double x, double y, double z)
{
  static int printed = 0;

#ifdef GRAVITY_TALLBOX
  Terminate("GRAVITY_TALLBOX is not implemented for MULTIPOLE_ORDER >= 4");
#endif

  symtensor7<double> D7 = 0.0;

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  if(printed == 0)
    {
      mpi_printf("EWALD: D7 table: qxmax=%d qymax=%d qzmax=%d   nxmax=%d nymax=%d nzmax=%d\n", qxmax, qymax, qzmax, nxmax, nymax,
                 nzmax);
      printed = 1;
    }

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv   = (r > 0) ? 1.0 / r : 0.0;
          double r2inv  = rinv * rinv;
          double r3inv  = r2inv * rinv;
          double r4inv  = r2inv * r2inv;
          double r5inv  = r2inv * r3inv;
          double r7inv  = r3inv * r4inv;
          double r9inv  = r4inv * r5inv;
          double r11inv = r4inv * r7inv;
          double r13inv = r4inv * r9inv;
          double r15inv = r4inv * r11inv;

          double g4, g5, g6, g7;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g4 = -(105.0 * erfc(alpha * r) +
                     (210.0 * alpha * r + 140.0 * pow(alpha * r, 3) + 56.0 * pow(alpha * r, 5) + 16.0 * pow(alpha * r, 7)) /
                         sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r9inv;

              g5 = (945.0 * erfc(alpha * r) + (1890.0 * alpha * r + 1260.0 * pow(alpha * r, 3) + 504.0 * pow(alpha * r, 5) +
                                               144.0 * pow(alpha * r, 7) + 32.0 * pow(alpha * r, 9)) /
                                                  sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r11inv;

              g6 = (-10395.0 * erfc(alpha * r) -
                    2.0 *
                        (10395.0 * alpha * r + 6930.0 * pow(alpha * r, 3) + 2772.0 * pow(alpha * r, 5) + 792.0 * pow(alpha * r, 7) +
                         176.0 * pow(alpha * r, 9) + 32.0 * pow(alpha * r, 11)) /
                        sqrt(M_PI) * exp(-alpha2 * r2)) *
                   r13inv;

              g7 =
                  (135135.0 * erfc(alpha * r) + 2.0 *
                                                    (135135.0 * alpha * r + 90090.0 * pow(alpha * r, 3) + 36036.0 * pow(alpha * r, 5) +
                                                     10296.0 * pow(alpha * r, 7) + 2288.0 * pow(alpha * r, 9) +
                                                     416.0 * pow(alpha * r, 11) + 64.0 * pow(alpha * r, 13)) /
                                                    sqrt(M_PI) * exp(-alpha2 * r2)) *
                  r15inv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* Note, for small r:
               *
               *   [1/- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               *
               *   Hence for r = 0:
               *
               *   g0 =  2     * alpha   / sqrt(pi)
               *   g1 = -4/3   * alpha^3 / sqrt(pi)
               *   g2 =  8/5   * alpha^5 / sqrt(pi)
               *   g3 = -16/7  * alpha^7 / sqrt(pi)
               *   g4 =  32/9  * alpha^9 / sqrt(pi)
               *   g5 = -64/11 * alpha^11/ sqrt(pi)
               */

              if((alpha * r) < 0.5)
                {
                  g4 = 32.0 * pow(alpha, 9) / sqrt(M_PI) *
                       (1.0 / 9.0 - pow(alpha * r, 2) / 11.0 + pow(alpha * r, 4) / 26.0 - pow(alpha * r, 6) / 90.0 +
                        pow(alpha * r, 8) / 408.0 - pow(alpha * r, 10) / 2280.0);

                  g5 = 64.0 * pow(alpha, 11) / sqrt(M_PI) *
                       (-1.0 / 11.0 + pow(alpha * r, 2) / 13.0 - pow(alpha * r, 4) / 30.0 + pow(alpha * r, 6) / 102.0 -
                        pow(alpha * r, 8) / 456.0 + pow(alpha * r, 10) / 2520.0);

                  g6 = 128.0 * pow(alpha, 13) / sqrt(M_PI) *
                       (1.0 / 13.0 - pow(alpha * r, 2) / 15.0 + pow(alpha * r, 4) / 34.0 - pow(alpha * r, 6) / 114.0 +
                        pow(alpha * r, 8) / 504.0 - pow(alpha * r, 10) / 2760.0);

                  g7 = 256.0 * pow(alpha, 15) / sqrt(M_PI) *
                       (-1.0 / 15.0 + pow(alpha * r, 2) / 17.0 - pow(alpha * r, 4) / 38.0 + pow(alpha * r, 6) / 126.0 -
                        pow(alpha * r, 8) / 552.0 + pow(alpha * r, 10) / 3000.0);
                }
              else
                {
                  g4 = (105.0 * erf(alpha * r) -
                        (210.0 * alpha * r + 140.0 * pow(alpha * r, 3) + 56.0 * pow(alpha * r, 5) + 16.0 * pow(alpha * r, 7)) /
                            sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r9inv;

                  g5 = (-945.0 * erf(alpha * r) + (1890.0 * alpha * r + 1260.0 * pow(alpha * r, 3) + 504.0 * pow(alpha * r, 5) +
                                                   144.0 * pow(alpha * r, 7) + 32.0 * pow(alpha * r, 9)) /
                                                      sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r11inv;

                  g6 = (10395.0 * erf(alpha * r) -
                        2.0 *
                            (10395.0 * alpha * r + 6930.0 * pow(alpha * r, 3) + 2772.0 * pow(alpha * r, 5) +
                             792.0 * pow(alpha * r, 7) + 176.0 * pow(alpha * r, 9) + 32.0 * pow(alpha * r, 11)) /
                            sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r13inv;

                  g7 = (-135135.0 * erf(alpha * r) +
                        2.0 *
                            (135135.0 * alpha * r + 90090.0 * pow(alpha * r, 3) + 36036.0 * pow(alpha * r, 5) +
                             10296.0 * pow(alpha * r, 7) + 2288.0 * pow(alpha * r, 9) + 416.0 * pow(alpha * r, 11) +
                             64.0 * pow(alpha * r, 13)) /
                            sqrt(M_PI) * exp(-alpha2 * r2)) *
                       r15inv;
                }
            }

          setup_D7(ADD, D7, dxyz, g4, g5, g6, g7);
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          double kx = (2.0 * M_PI * LONG_X) * nx;
          double ky = (2.0 * M_PI * LONG_Y) * ny;
          double kz = (2.0 * M_PI * LONG_Z) * nz;
          double k2 = kx * kx + ky * ky + kz * kz;

          if(k2 > 0)
            {
              double kdotx = (x * kx + y * ky + z * kz);
              double val   = -4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * sin(kdotx);

              vector<double> kxyz(kx, ky, kz);

              D7 += (val * kxyz) % (kxyz % (kxyz % (kxyz % (kxyz % (kxyz % kxyz)))));
            }
        }

  return D7;
}
