/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file ewald_test.cc
 *
 *  \brief some unit test routines for the table look-up in the ewald correction
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
#include "../sort/cxxsort.h"
#include "../system/system.h"

#ifdef EWALD_TEST

/*!
 *  We here use trilinear interpolation
 *  to get it from the precomputed tables, which contain one octant
 *  around the target particle at the origin. The other octants are
 *  obtained from it by exploiting the symmetry properties.
 *
 *  \param dx x component of the distance between the two particles
 *  \param dx y component of the distance between the two particles
 *  \param dx z component of the distance between the two particles
 *  \param fper pointer to array containing the correction force
 *  \param fper[4] contains the correction potential
 */
void ewald::ewald_corr(double dx, double dy, double dz, enum interpolate_options flag, ewald_data &fper)
{
  if(!ewald_is_initialized)
    Terminate("How come that Ewald tables are not initialized?");

  int signx, signy, signz;

  if(dx < 0)
    {
      dx    = -dx;
      signx = -1;
    }
  else
    signx = +1;

  if(dy < 0)
    {
      dy    = -dy;
      signy = -1;
    }
  else
    signy = +1;

  if(dz < 0)
    {
      dz    = -dz;
      signz = -1;
    }
  else
    signz = +1;

  ewald_interpolate(dx, dy, dz, flag, fper);

  /* change signs as needed */

  fper.D1phi[0] *= signx;
  fper.D1phi[1] *= signy;
  fper.D1phi[2] *= signz;

  if(flag == POINTMASS)
    return;

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

  fper.D4phi[sXXXY] *= signx * signy;
  fper.D4phi[sXYYY] *= signx * signy;
  fper.D4phi[sXXXZ] *= signx * signz;
  fper.D4phi[sXZZZ] *= signx * signz;
  fper.D4phi[sYYYZ] *= signy * signz;
  fper.D4phi[sYZZZ] *= signy * signz;
  fper.D4phi[sXXYZ] *= signy * signz;
  fper.D4phi[sXYYZ] *= signx * signz;
  fper.D4phi[sXYZZ] *= signx * signy;

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
}

void ewald::ewald_corr_exact(double dx, double dy, double dz, enum interpolate_options flag, ewald_data &fper)
{
  double fac = 1.0 / All.BoxSize;
  double x   = dx * fac;
  double y   = dy * fac;
  double z   = dz * fac;

  fper.D0phi = pow(fac, 1) * ewald_D0(x, y, z);
  fper.D1phi = pow(fac, 2) * ewald_D1(x, y, z);

  if(flag == POINTMASS)
    return;

  fper.D2phi = pow(fac, 3) * ewald_D2(x, y, z);
  fper.D3phi = pow(fac, 4) * ewald_D3(x, y, z);
  fper.D4phi = pow(fac, 5) * ewald_D4(x, y, z);
  fper.D5phi = pow(fac, 6) * ewald_D5(x, y, z);
}

void ewald::ewald_interpolate(double dx, double dy, double dz, enum interpolate_options flag, ewald_data &fper)
{
  double u = dx * Ewd_fac_intp[0];
  int i    = (int)u;
  if(i >= ENX)
    i = ENX - 1;
  u -= i;

  double v = dy * Ewd_fac_intp[1];
  int j    = (int)v;
  if(j >= ENY)
    j = ENY - 1;
  v -= j;

  double w = dz * Ewd_fac_intp[2];
  int k    = (int)w;
  if(k >= ENZ)
    k = ENZ - 1;
  w -= k;

  double f1 = (1 - u) * (1 - v) * (1 - w);
  double f2 = (1 - u) * (1 - v) * (w);
  double f3 = (1 - u) * (v) * (1 - w);
  double f4 = (1 - u) * (v) * (w);
  double f5 = (u) * (1 - v) * (1 - w);
  double f6 = (u) * (1 - v) * (w);
  double f7 = (u) * (v) * (1 - w);
  double f8 = (u) * (v) * (w);

  ewald_data &C1 = Ewd[ewd_offset(i, j, k)];
  ewald_data &C2 = Ewd[ewd_offset(i, j, k + 1)];
  ewald_data &C3 = Ewd[ewd_offset(i, j + 1, k)];
  ewald_data &C4 = Ewd[ewd_offset(i, j + 1, k + 1)];
  ewald_data &C5 = Ewd[ewd_offset(i + 1, j, k)];
  ewald_data &C6 = Ewd[ewd_offset(i + 1, j, k + 1)];
  ewald_data &C7 = Ewd[ewd_offset(i + 1, j + 1, k)];
  ewald_data &C8 = Ewd[ewd_offset(i + 1, j + 1, k + 1)];

#ifdef EVALPOTENTIAL
  fper.D0phi =
      f1 * C1.D0phi + f2 * C2.D0phi + f3 * C3.D0phi + f4 * C4.D0phi + f5 * C5.D0phi + f6 * C6.D0phi + f7 * C7.D0phi + f8 * C8.D0phi;
#endif
  fper.D1phi =
      f1 * C1.D1phi + f2 * C2.D1phi + f3 * C3.D1phi + f4 * C4.D1phi + f5 * C5.D1phi + f6 * C6.D1phi + f7 * C7.D1phi + f8 * C8.D1phi;

  if(flag == POINTMASS)
    return;

  fper.D2phi =
      f1 * C1.D2phi + f2 * C2.D2phi + f3 * C3.D2phi + f4 * C4.D2phi + f5 * C5.D2phi + f6 * C6.D2phi + f7 * C7.D2phi + f8 * C8.D2phi;

  fper.D3phi =
      f1 * C1.D3phi + f2 * C2.D3phi + f3 * C3.D3phi + f4 * C4.D3phi + f5 * C5.D3phi + f6 * C6.D3phi + f7 * C7.D3phi + f8 * C8.D3phi;

  fper.D4phi =
      f1 * C1.D4phi + f2 * C2.D4phi + f3 * C3.D4phi + f4 * C4.D4phi + f5 * C5.D4phi + f6 * C6.D4phi + f7 * C7.D4phi + f8 * C8.D4phi;

  fper.D5phi =
      f1 * C1.D5phi + f2 * C2.D5phi + f3 * C3.D5phi + f4 * C4.D5phi + f5 * C5.D5phi + f6 * C6.D5phi + f7 * C7.D5phi + f8 * C8.D5phi;
}

void ewald::test_interpolation_accuracy(void)
{
  init_rng(42);

  printf("\n\n");

  double errsum0 = 0;
  double errmax0 = 0;
  double count0  = 0;
  for(int i = 0; i < 1000; i++)
    {
      double x = get_random_number() - 0.5;
      double y = get_random_number() - 0.5;
      double z = get_random_number() - 0.5;

      double r = sqrt(x * x + y * y + z * z);

      double D0phi_exact = (1.0 / All.BoxSize) * ewald_D0(x, y, z);

      ewald_data ew;
      ewald_corr(x * All.BoxSize, y * All.BoxSize, z * All.BoxSize, ewald::MULTIPOLES, ew);

      double err = fabs((D0phi_exact - ew.D0phi) / D0phi_exact);

      errsum0 += err;
      if(err > errmax0)
        errmax0 = err;
      count0++;

      if(err > 0.1)
        printf("%4d   r=%g D0_exact=%g  D0_interpol=%g  rel_error=%g\n", i, r, D0phi_exact, ew.D0phi, err);
    }

  double errsum1 = 0;
  double errmax1 = 0;
  double count1  = 0;
  for(int i = 0; i < 1000; i++)
    {
      double x = get_random_number() - 0.5;
      double y = get_random_number() - 0.5;
      double z = get_random_number() - 0.5;

      double r = sqrt(x * x + y * y + z * z);

      vector<double> D1phi_exact = pow(All.BoxSize, -2) * ewald_D1(x, y, z);

      ewald_data ew;
      ewald_corr(x * All.BoxSize, y * All.BoxSize, z * All.BoxSize, ewald::MULTIPOLES, ew);

      double norm = D1phi_exact.norm();

      for(int j = 0; j < 3; j++)
        {
          double err = fabs(D1phi_exact[j] - ew.D1phi[j]) / norm;
          errsum1 += err;
          if(err > errmax1)
            errmax1 = err;
          count1++;

          if(err > 0.1)
            {
              printf("%4d  r=%g bin=%d  D1_exact[%d]=%g  D1_interpol[%d]=%g  rel_error=%g\n", i, r,
                     (int)(r * All.BoxSize * Ewd_fac_intp[0]), j, D1phi_exact[j], j, ew.D1phi[j], err);
            }
        }
    }

  double errsum2 = 0;
  double errmax2 = 0;
  double count2  = 0;
  for(int i = 0; i < 1000; i++)
    {
      double x = get_random_number() - 0.5;
      double y = get_random_number() - 0.5;
      double z = get_random_number() - 0.5;

      double r = sqrt(x * x + y * y + z * z);

      symtensor2<double> D2phi_exact = pow(All.BoxSize, -3) * ewald_D2(x, y, z);

      ewald_data ew;
      ewald_corr(x * All.BoxSize, y * All.BoxSize, z * All.BoxSize, ewald::MULTIPOLES, ew);

      double norm = D2phi_exact.norm();

      for(int j = 0; j < 6; j++)
        {
          double err = fabs((D2phi_exact[j] - ew.D2phi[j]) / norm);
          errsum2 += err;
          if(err > errmax2)
            errmax2 = err;
          count2++;

          if(err > 0.1)
            printf("%4d  r=%g  D2_exact[%d]=%g  D2_interpol[%d]=%g  rel_error=%g\n", i, r, j, D2phi_exact[j], j, ew.D2phi[j], err);
        }
    }

  double errsum3 = 0;
  double errmax3 = 0;
  double count3  = 0;
  for(int i = 0; i < 1000; i++)
    {
      double x = get_random_number() - 0.5;
      double y = get_random_number() - 0.5;
      double z = get_random_number() - 0.5;

      double r = sqrt(x * x + y * y + z * z);

      symtensor3<double> D3phi_exact = pow(All.BoxSize, -4) * ewald_D3(x, y, z);

      ewald_data ew;
      ewald_corr(x * All.BoxSize, y * All.BoxSize, z * All.BoxSize, ewald::MULTIPOLES, ew);

      double norm = D3phi_exact.norm();

      for(int j = 0; j < 10; j++)
        {
          double err = fabs((D3phi_exact[j] - ew.D3phi[j]) / norm);
          errsum3 += err;
          if(err > errmax3)
            errmax3 = err;
          count3++;

          if(err > 0.1)
            printf("%4d   r=%g  D3_exact[%d]=%g  D3_interpol[%d]=%g  rel_error=%g\n", i, r, j, D3phi_exact[j], j, ew.D3phi[j], err);
        }
    }

  double errsum4 = 0;
  double errmax4 = 0;
  double count4  = 0;
  for(int i = 0; i < 1000; i++)
    {
      double x = get_random_number() - 0.5;
      double y = get_random_number() - 0.5;
      double z = get_random_number() - 0.5;

      double r = sqrt(x * x + y * y + z * z);

      symtensor4<double> D4phi_exact = pow(All.BoxSize, -5) * ewald_D4(x, y, z);

      ewald_data ew;
      ewald_corr(x * All.BoxSize, y * All.BoxSize, z * All.BoxSize, ewald::MULTIPOLES, ew);

      double norm = D4phi_exact.norm();

      for(int j = 0; j < 15; j++)
        {
          double err = fabs((D4phi_exact[j] - ew.D4phi[j]) / norm);
          errsum4 += err;
          if(err > errmax4)
            errmax4 = err;
          count4++;

          if(err > 0.1)
            printf("%4d  r=%g  D4_exact[%d]=%g  D4_interpol[%d]=%g  rel_error=%g\n", i, r, j, D4phi_exact[j], j, ew.D4phi[j], err);
        }
    }

  double errsum5 = 0;
  double errmax5 = 0;
  double count5  = 0;
  for(int i = 0; i < 1000; i++)
    {
      double x = get_random_number() - 0.5;
      double y = get_random_number() - 0.5;
      double z = get_random_number() - 0.5;

      symtensor5<double> D5phi_exact = pow(All.BoxSize, -6) * ewald_D5(x, y, z);

      ewald_data ew;
      ewald_corr(x * All.BoxSize, y * All.BoxSize, z * All.BoxSize, ewald::MULTIPOLES, ew);

      double norm = D5phi_exact.norm();

      for(int j = 0; j < 21; j++)
        {
          double err = fabs((D5phi_exact[j] - ew.D5phi[j]) / norm);
          errsum5 += err;
          if(err > errmax5)
            errmax5 = err;
          count5++;
        }
    }

  printf("\n\n");

  printf("D0:   max error = %g   mean error=%g\n", errmax0, errsum0 / count0);
  printf("D1:   max error = %g   mean error=%g\n", errmax1, errsum1 / count1);
  printf("D2:   max error = %g   mean error=%g\n", errmax2, errsum2 / count2);
  printf("D3:   max error = %g   mean error=%g\n", errmax3, errsum3 / count3);
  printf("D4:   max error = %g   mean error=%g\n", errmax4, errsum4 / count4);
  printf("D5:   max error = %g   mean error=%g\n", errmax5, errsum5 / count5);

  printf("\n\n");

  {
    double errsum0 = 0;
    double errmax0 = 0;
    double count0  = 0;
    for(int i = 0; i < 1000; i++)
      {
        double x = get_random_number() - 0.5;
        double y = get_random_number() - 0.5;
        double z = get_random_number() - 0.5;

        double D0phi_exact = (1.0 / All.BoxSize) * ewald_D0(x, y, z);

        ewald_data ew;
        double posd[3] = {x * All.BoxSize, y * All.BoxSize, z * All.BoxSize};
        MyIntPosType pos[3];
        pos_to_signedintpos(posd, (MySignedIntPosType *)pos);
        MyIntPosType ref[3] = {0, 0, 0};
        ewald_gridlookup(pos, ref, ewald::MULTIPOLES, ew);

        double err = fabs((D0phi_exact - ew.D0phi) / D0phi_exact);

        errsum0 += err;
        if(err > errmax0)
          errmax0 = err;
        count0++;
      }

    double errsum1 = 0;
    double errmax1 = 0;
    double count1  = 0;
    for(int i = 0; i < 1000; i++)
      {
        double x = get_random_number() - 0.5;
        double y = get_random_number() - 0.5;
        double z = get_random_number() - 0.5;

        vector<double> D1phi_exact = pow(All.BoxSize, -2) * ewald_D1(x, y, z);

        ewald_data ew;
        double posd[3] = {x * All.BoxSize, y * All.BoxSize, z * All.BoxSize};
        MyIntPosType pos[3];
        pos_to_signedintpos(posd, (MySignedIntPosType *)pos);
        MyIntPosType ref[3] = {0, 0, 0};
        ewald_gridlookup(pos, ref, ewald::MULTIPOLES, ew);

        double norm = D1phi_exact.norm();

        for(int j = 0; j < 3; j++)
          {
            double err = fabs(D1phi_exact[j] - ew.D1phi[j]) / norm;

            errsum1 += err;
            if(err > errmax1)
              errmax1 = err;
            count1++;
          }
      }

    double errsum2 = 0;
    double errmax2 = 0;
    double count2  = 0;
    for(int i = 0; i < 1000; i++)
      {
        double x = get_random_number() - 0.5;
        double y = get_random_number() - 0.5;
        double z = get_random_number() - 0.5;

        symtensor2<double> D2phi_exact = pow(All.BoxSize, -3) * ewald_D2(x, y, z);

        ewald_data ew;
        double posd[3] = {x * All.BoxSize, y * All.BoxSize, z * All.BoxSize};
        MyIntPosType pos[3];
        pos_to_signedintpos(posd, (MySignedIntPosType *)pos);
        MyIntPosType ref[3] = {0, 0, 0};
        ewald_gridlookup(pos, ref, ewald::MULTIPOLES, ew);

        double norm = D2phi_exact.norm();

        for(int j = 0; j < 6; j++)
          {
            double err = fabs((D2phi_exact[j] - ew.D2phi[j]) / norm);
            errsum2 += err;
            if(err > errmax2)
              errmax2 = err;
            count2++;
          }
      }

    double errsum3 = 0;
    double errmax3 = 0;
    double count3  = 0;
    for(int i = 0; i < 1000; i++)
      {
        double x = get_random_number() - 0.5;
        double y = get_random_number() - 0.5;
        double z = get_random_number() - 0.5;

        symtensor3<double> D3phi_exact = pow(All.BoxSize, -4) * ewald_D3(x, y, z);

        ewald_data ew;
        double posd[3] = {x * All.BoxSize, y * All.BoxSize, z * All.BoxSize};
        MyIntPosType pos[3];
        pos_to_signedintpos(posd, (MySignedIntPosType *)pos);
        MyIntPosType ref[3] = {0, 0, 0};
        ewald_gridlookup(pos, ref, ewald::MULTIPOLES, ew);

        double norm = D3phi_exact.norm();

        for(int j = 0; j < 10; j++)
          {
            double err = fabs((D3phi_exact[j] - ew.D3phi[j]) / norm);
            errsum3 += err;
            if(err > errmax3)
              errmax3 = err;
            count3++;
          }
      }

    double errsum4 = 0;
    double errmax4 = 0;
    double count4  = 0;
    for(int i = 0; i < 1000; i++)
      {
        double x = get_random_number() - 0.5;
        double y = get_random_number() - 0.5;
        double z = get_random_number() - 0.5;

        symtensor4<double> D4phi_exact = pow(All.BoxSize, -5) * ewald_D4(x, y, z);

        ewald_data ew;
        double posd[3] = {x * All.BoxSize, y * All.BoxSize, z * All.BoxSize};
        MyIntPosType pos[3];
        pos_to_signedintpos(posd, (MySignedIntPosType *)pos);
        MyIntPosType ref[3] = {0, 0, 0};
        ewald_gridlookup(pos, ref, ewald::MULTIPOLES, ew);

        double norm = D4phi_exact.norm();

        for(int j = 0; j < 15; j++)
          {
            double err = fabs((D4phi_exact[j] - ew.D4phi[j]) / norm);
            errsum4 += err;
            if(err > errmax4)
              errmax4 = err;
            count4++;
          }
      }

    double errsum5 = 0;
    double errmax5 = 0;
    double count5  = 0;
    for(int i = 0; i < 1000; i++)
      {
        double x = get_random_number() - 0.5;
        double y = get_random_number() - 0.5;
        double z = get_random_number() - 0.5;

        symtensor5<double> D5phi_exact = pow(All.BoxSize, -6) * ewald_D5(x, y, z);

        ewald_data ew;
        double posd[3] = {x * All.BoxSize, y * All.BoxSize, z * All.BoxSize};
        MyIntPosType pos[3];
        pos_to_signedintpos(posd, (MySignedIntPosType *)pos);
        MyIntPosType ref[3] = {0, 0, 0};
        ewald_gridlookup(pos, ref, ewald::MULTIPOLES, ew);

        double norm = D5phi_exact.norm();

        for(int j = 0; j < 21; j++)
          {
            double err = fabs((D5phi_exact[j] - ew.D5phi[j]) / norm);
            errsum5 += err;
            if(err > errmax5)
              errmax5 = err;
            count5++;
          }
      }

    printf("Grid look-up: \n\n");

    printf("D0:   max error = %g   mean error=%g\n", errmax0, errsum0 / count0);
    printf("D1:   max error = %g   mean error=%g\n", errmax1, errsum1 / count1);
    printf("D2:   max error = %g   mean error=%g\n", errmax2, errsum2 / count2);
    printf("D3:   max error = %g   mean error=%g\n", errmax3, errsum3 / count3);
    printf("D4:   max error = %g   mean error=%g\n", errmax4, errsum4 / count4);
    printf("D5:   max error = %g   mean error=%g\n", errmax5, errsum5 / count5);
  }

  Terminate("stop");
}

ewaldtensor6<double> ewald::ewald_P6(void)
{
#ifdef GRAVITY_TALLBOX
  Terminate("GRAVITY_TALLBOX is not implemented");
#endif

  ewaldtensor6<double> P6 = 0.0;

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = -nx * (1.0 / LONG_X);
          double dy = -ny * (1.0 / LONG_Y);
          double dz = -nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv  = (r > 0) ? 1.0 / r : 0.0;
          double r2inv = rinv * rinv;
          double r3inv = r2inv * rinv;
          double r4inv = r2inv * r2inv;
          double r7inv = r3inv * r4inv;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              // derivatives of f(r)  = -Erfc[alpha * r] /r

              double ar   = alpha * r;
              double ar2  = ar * ar;
              double ar3  = ar * ar * ar;
              double ar5  = ar3 * ar * ar;
              double ar7  = ar5 * ar * ar;
              double ar9  = ar7 * ar * ar;
              double ar11 = ar9 * ar * ar;
              double xir2 = pow(dx * rinv, 2);
              double xir4 = pow(dx * rinv, 4);
              double xir6 = pow(dx * rinv, 6);

              double yir2 = pow(dy * rinv, 2);
              double zir2 = pow(dz * rinv, 2);

              P6.XXXXXX += (450.0 * ar + 300.0 * ar3 + 120.0 * ar5 - xir2 * (9450.0 * ar + 6300.0 * ar3 + 2520.0 * ar5 + 720.0 * ar7) +
                            xir4 * (28350.0 * ar + 18900.0 * ar3 + 7560.0 * ar5 + 2160.0 * ar7 + 480.0 * ar9) -
                            xir6 * (20790.0 * ar + 13860.0 * ar3 + 5544.0 * ar5 + 1584.0 * ar7 + 352.0 * ar9 + 64.0 * ar11)) /
                               sqrt(M_PI) * exp(-ar2) * r7inv +
                           erfc(ar) * (225.0 - 4725.0 * xir2 + 14175.0 * xir4 - 10395.0 * xir6) * r7inv;

              P6.XXXXYY += (90.0 * ar + 60.0 * ar3 + 24.0 * ar5 - xir2 * (1260.0 * ar + 840.0 * ar3 + 336.0 * ar5 + 96.0 * ar7) +
                            xir4 * (1890.0 * ar + 1260.0 * ar3 + 504.0 * ar5 + 144.0 * ar7 + 32.0 * ar9) -
                            yir2 * (630.0 * ar + 420.0 * ar3 + 168.0 * ar5 + 48.0 * ar7) +
                            xir2 * yir2 * (11340.0 * ar + 7560.0 * ar3 + 3024.0 * ar5 + 864.0 * ar7 + 192.0 * ar9) -
                            xir4 * yir2 * (20790.0 * ar + 13860.0 * ar3 + 5544.0 * ar5 + 1584.0 * ar7 + 352.0 * ar9 + 64.0 * ar11)) /
                               sqrt(M_PI) * exp(-ar2) * r7inv +
                           erfc(ar) *
                               (45.0 - 630.0 * xir2 + 945.0 * xir4 - 315.0 * yir2 + 5670.0 * xir2 * yir2 - 10395.0 * xir4 * yir2) *
                               r7inv;

              P6.XXYYZZ +=
                  (30.0 * ar + 20.0 * ar3 + 8.0 * ar5 - (xir2 + yir2 + zir2) * (210.0 * ar + 140.0 * ar3 + 56.0 * ar5 + 16.0 * ar7) +
                   +(xir2 * yir2 + xir2 * zir2 + yir2 * zir2) * (1890.0 * ar + 1260.0 * ar3 + 504.0 * ar5 + 144.0 * ar7 + 32.0 * ar9) -
                   xir2 * yir2 * zir2 * (20790.0 * ar + 13860.0 * ar3 + 5544.0 * ar5 + 1584.0 * ar7 + 352.0 * ar9 + 64.0 * ar11)) /
                      sqrt(M_PI) * exp(-ar2) * r7inv +
                  erfc(ar) *
                      (15.0 - 105.0 * (xir2 + yir2 + zir2) + 945.0 * (xir2 * yir2 + xir2 * zir2 + yir2 * zir2) -
                       10395.0 * xir2 * yir2 * zir2) *
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
               */

              P6.XXXXXX += (-240.0) * pow(alpha, 7) / (7.0 * sqrt(M_PI));
              P6.XXXXYY += (-48.0) * pow(alpha, 7) / (7.0 * sqrt(M_PI));
              P6.XXYYZZ += 0;
            }
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

              double val = 4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2));

              P6.XXXXXX += val * pow(kx, 6);
              P6.XXXXYY += val * pow(kx, 4) * pow(ky, 2);
              P6.XXYYZZ += val * pow(kx, 2) * pow(ky, 2) * pow(kz, 2);
            }
        }

  return P6;
}

ewaldtensor8<double> ewald::ewald_P8(void)
{
#ifdef GRAVITY_TALLBOX
  Terminate("GRAVITY_TALLBOX is not implemented");
#endif

  ewaldtensor8<double> P8 = 0.0;

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = -nx * (1.0 / LONG_X);
          double dy = -ny * (1.0 / LONG_Y);
          double dz = -nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv  = (r > 0) ? 1.0 / r : 0.0;
          double r2inv = rinv * rinv;
          double r3inv = r2inv * rinv;
          double r4inv = r2inv * r2inv;
          double r5inv = r2inv * r3inv;
          double r9inv = r4inv * r5inv;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              // derivatives of f(r)  = -Erfc[alpha * r] /r

              double ar   = alpha * r;
              double ar2  = ar * ar;
              double ar3  = ar * ar * ar;
              double ar5  = ar3 * ar * ar;
              double ar7  = ar5 * ar * ar;
              double ar9  = ar7 * ar * ar;
              double ar11 = ar9 * ar * ar;
              double ar13 = ar11 * ar * ar;
              double ar15 = ar13 * ar * ar;
              double xir2 = pow(dx * rinv, 2);
              double xir4 = pow(dx * rinv, 4);
              double xir6 = pow(dx * rinv, 6);
              double xir8 = pow(dx * rinv, 8);

              double yir2 = pow(dy * rinv, 2);
              double yir4 = pow(dy * rinv, 4);
              double zir2 = pow(dz * rinv, 2);

              P8.XXXXXXXX +=
                  (-(22050.0 * ar + 14700.0 * ar3 + 5880.0 * ar5 + 1680.0 * ar7) +
                   xir2 * (793800.0 * ar + 529200.0 * ar3 + 211680.0 * ar5 + 60480.0 * ar7 + 13440.0 * ar9) -
                   xir4 * (4365900.0 * ar + 2910600.0 * ar3 + 1164240.0 * ar5 + 332640.0 * ar7 + 73920.0 * ar9 + 13440.0 * ar11) +
                   xir6 * (7567560.0 * ar + 5045040.0 * ar3 + 2018016.0 * ar5 + 576576.0 * ar7 + 128128.0 * ar9 + 23296.0 * ar11 +
                           3584 * ar13) -
                   xir8 * (4054050.0 * ar + 2702700.0 * ar3 + 1081080.0 * ar5 + 308880.0 * ar7 + 68640.0 * ar9 + 12480.0 * ar11 +
                           1920.0 * ar13 + 256.0 * ar15)) /
                      sqrt(M_PI) * exp(-ar2) * r9inv +
                  erfc(ar) * (-11025.0 + 396900.0 * xir2 - 2182950.0 * xir4 + 3783780.0 * xir6 - 2027025.0 * xir8) * r9inv;

              P8.XXXXXXYY =
                  (-(3150.0 * ar + 2100.0 * ar3 + 840.0 * ar5 + 240.0 * ar7) +
                   xir2 * (85050.0 * ar + 56700.0 * ar3 + 22680.0 * ar5 + 6480.0 * ar7 + 1440.0 * ar9) -
                   xir4 * (311850.0 * ar + 207900.0 * ar3 + 83160.0 * ar5 + 23760.0 * ar7 + 5280.0 * ar9 + 960.0 * ar11) +
                   xir6 *
                       (270270.0 * ar + 180180.0 * ar3 + 72072.0 * ar5 + 20592.0 * ar7 + 4576.0 * ar9 + 832.0 * ar11 + 128.0 * ar13) +
                   yir2 * (28350.0 * ar + 18900 * ar3 + 7560.0 * ar5 + 2160 * ar7 + 480 * ar9) -
                   xir2 * yir2 * (935550 * ar + 623700.0 * ar3 + 249480.0 * ar5 + 71280.0 * ar7 + 15840.0 * ar9 + 2880.0 * ar11) +
                   xir4 * yir2 *
                       (4054050.0 * ar + 2702700.0 * ar3 + 1081080.0 * ar5 + 308880.0 * ar7 + 68640.0 * ar9 + 12480.0 * ar11 +
                        1920.0 * ar13) -
                   xir6 * yir2 *
                       (4054050.0 * ar + 2702700.0 * ar3 + 1081080.0 * ar5 + 308880.0 * ar7 + 68640.0 * ar9 + 12480.0 * ar11 +
                        1920.0 * ar13 + 256.0 * ar15)) /
                      sqrt(M_PI) * exp(-ar2) * r9inv +
                  erfc(ar) *
                      (-1575.0 + 42525.0 * xir2 - 155925.0 * xir4 + 135135.0 * xir6 + 14175.0 * yir2 - 467775.0 * xir2 * yir2 +
                       2027025.0 * xir4 * yir2 - 2027025 * xir6 * yir2) *
                      r9inv;

              P8.XXXXYYYY =
                  (-(1890.0 * ar + 1260.0 * ar3 + 504.0 * ar5 + 144.0 * ar7) +
                   (xir2 + yir2) * (34020.0 * ar + 22680.0 * ar3 + 9072.0 * ar5 + 2592.0 * ar7 + 576.0 * ar9) -
                   (xir4 + yir4) * (62370.0 * ar + 41580.0 * ar3 + 16632.0 * ar5 + 4752.0 * ar7 + 1056.0 * ar9 + 192.0 * ar11) +
                   -xir2 * yir2 * (748440.0 * ar + 498960.0 * ar3 + 199584.0 * ar5 + 57024.0 * ar7 + 12672.0 * ar9 + 2304.0 * ar11) +

                   (xir4 * yir2 + xir2 * yir4) * (1621620.0 * ar + 1081080.0 * ar3 + 432432.0 * ar5 + 123552.0 * ar7 + 27456.0 * ar9 +
                                                  4992.0 * ar11 + 768.0 * ar13) -
                   xir4 * yir4 *
                       (4054050.0 * ar + 2702700.0 * ar3 + 1081080.0 * ar5 + 308880.0 * ar7 + 68640.0 * ar9 + 12480.0 * ar11 +
                        1920.0 * ar13 + 256.0 * ar15)) /
                      sqrt(M_PI) * exp(-ar2) * r9inv +
                  erfc(ar) *
                      (-945.0 + 17010.0 * (xir2 + yir2) - 31185.0 * (xir4 + yir4) - 374220.0 * xir2 * yir2 +
                       810810.0 * (xir4 * yir2 + xir2 * yir4) - 2027025.0 * xir4 * yir4) *
                      r9inv;

              P8.XXXXYYZZ = (-(630.0 * ar + 420.0 * ar3 + 168.0 * ar5 + 48.0 * ar7) +
                             (xir2) * (11340.0 * ar + 7560.0 * ar3 + 3024.0 * ar5 + 864.0 * ar7 + 192.0 * ar9) -
                             (xir4) * (20790.0 * ar + 13860.0 * ar3 + 5544.0 * ar5 + 1584.0 * ar7 + 352.0 * ar9 + 64.0 * ar11) +
                             +(yir2 + zir2) * (5670.0 * ar + 3780.0 * ar3 + 1512.0 * ar5 + 432.0 * ar7 + 96.0 * ar9) -
                             (xir2 * yir2 + xir2 * zir2) *
                                 (124740 * ar + 83160.0 * ar3 + 33264.0 * ar5 + 9504 * ar7 + 2112.0 * ar9 + 384.0 * ar11) +
                             (xir4 * yir2 + xir4 * zir2) * (270270.0 * ar + 180180.0 * ar3 + 72072 * ar5 + 20592.0 * ar7 +
                                                            4576.0 * ar9 + 832.0 * ar11 + 128.0 * ar13) -
                             yir2 * zir2 * (62370.0 * ar + 41580 * ar3 + 16632.0 * ar5 + 4752 * ar7 + 1056.0 * ar9 + 192.0 * ar11) +
                             (xir2 * yir2 * zir2) * (1621620.0 * ar + 1081080.0 * ar3 + 432432.0 * ar5 + 123552.0 * ar7 +
                                                     27456.0 * ar9 + 4992.0 * ar11 + 768.0 * ar13) -
                             xir4 * yir2 * zir2 *
                                 (4054050.0 * ar + 2702700.0 * ar3 + 1081080.0 * ar5 + 308880.0 * ar7 + 68640.0 * ar9 +
                                  12480.0 * ar11 + 1920.0 * ar13 + 256.0 * ar15)) /
                                sqrt(M_PI) * exp(-ar2) * r9inv +
                            erfc(ar) *
                                (-315.0 + 5670.0 * xir2 - 10395.0 * xir4 + 2835.0 * (yir2 + zir2) -
                                 62370.0 * (xir2 * yir2 + xir2 * zir2) + 135135.0 * xir4 * (yir2 + zir2) - 31185.0 * yir2 * zir2 +
                                 810810.0 * xir2 * yir2 * zir2 - 2027025 * xir4 * yir2 * zir2) *
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
               */

              P8.XXXXXXXX += 1120.0 * pow(alpha, 9) / (3.0 * sqrt(M_PI));
              P8.XXXXXXYY += 160.0 * pow(alpha, 9) / (3.0 * sqrt(M_PI));
              P8.XXXXYYYY += 0;
              P8.XXXXYYZZ += 32.0 * pow(alpha, 9) / (3.0 * sqrt(M_PI));
            }
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

              double val = -4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2));

              P8.XXXXXXXX += val * pow(kx, 8);
              P8.XXXXXXYY += val * pow(kx, 6) * pow(ky, 2);
              P8.XXXXYYYY += val * pow(kx, 4) * pow(ky, 4);
              P8.XXXXYYZZ += val * pow(kx, 4) * pow(ky, 2) * pow(kz, 2);
            }
        }

  return P8;
}

ewaldtensor10<double> ewald::ewald_P10(void)
{
#ifdef GRAVITY_TALLBOX
  Terminate("GRAVITY_TALLBOX is not implemented");
#endif

  ewaldtensor10<double> P10 = 0.0;

  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = -nx * (1.0 / LONG_X);
          double dy = -ny * (1.0 / LONG_Y);
          double dz = -nz * (1.0 / LONG_Z);

          vector<double> dxyz(dx, dy, dz);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv   = (r > 0) ? 1.0 / r : 0.0;
          double r2inv  = rinv * rinv;
          double r3inv  = r2inv * rinv;
          double r4inv  = r2inv * r2inv;
          double r7inv  = r3inv * r4inv;
          double r11inv = r4inv * r7inv;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              double ar   = alpha * r;
              double ar2  = ar * ar;
              double ar3  = ar * ar * ar;
              double ar5  = ar3 * ar * ar;
              double ar7  = ar5 * ar * ar;
              double ar9  = ar7 * ar * ar;
              double ar11 = ar9 * ar * ar;
              double ar13 = ar11 * ar * ar;
              double ar15 = ar13 * ar * ar;
              double ar17 = ar15 * ar * ar;
              double ar19 = ar17 * ar * ar;

              double xir2  = pow(dx * rinv, 2);
              double xir4  = pow(dx * rinv, 4);
              double xir6  = pow(dx * rinv, 6);
              double xir8  = pow(dx * rinv, 8);
              double xir10 = pow(dx * rinv, 10);

              double yir2 = pow(dy * rinv, 2);
              double yir4 = pow(dy * rinv, 4);
              double zir2 = pow(dz * rinv, 2);

              P10.XXXXXXXXXX +=
                  ((1786050.0 * ar + 1190700.0 * ar3 + 476280.0 * ar5 + 136080.0 * ar7 + 30240 * ar9) -
                   xir2 * (98232750.0 * ar + 65488500.0 * ar3 + 26195400.0 * ar5 + 7484400.0 * ar7 + 1663200.0 * ar9 + 302400 * ar11) +
                   xir4 * (851350500.0 * ar + 567567000.0 * ar3 + 227026800.0 * ar5 + 64864800.0 * ar7 + 14414400.0 * ar9 +
                           2620800.0 * ar11 + 403200 * ar13) -
                   xir6 * (2554051500.0 * ar + 1702701000.0 * ar3 + 681080400.0 * ar5 + 194594400.0 * ar7 + 43243200.0 * ar9 +
                           7862400.0 * ar11 + 1209600.0 * ar13 + 161280 * ar15) +
                   xir8 * (3101348250.0 * ar + 2067565500.0 * ar3 + 827026200.0 * ar5 + 236293200.0 * ar7 + 52509600.0 * ar9 +
                           9547200.0 * ar11 + 1468800.0 * ar13 + 195840.0 * ar15 + 23040 * ar17) -
                   xir10 * (1309458150 * ar + 872972100 * ar3 + 349188840 * ar5 + 99768240 * ar7 + 22170720 * ar9 + 4031040 * ar11 +
                            620160 * ar13 + 82688 * ar15 + 9728 * ar17 + 1024.0 * ar19)) /
                      sqrt(M_PI) * exp(-ar2) * r11inv +
                  erfc(ar) *
                      (893025.0 - 49116375.0 * xir2 + 425675250.0 * xir4 - 1277025750 * xir6 + 1550674125.0 * xir8 -
                       654729075.0 * xir10) *
                      r11inv;

              P10.XXXXXXXXYY +=
                  ((198450.0 * ar + 132300.0 * ar3 + 52920.0 * ar5 + 15120.0 * ar7 + 3360.0 * ar9) -
                   xir2 * (8731800.0 * ar + 5821200.0 * ar3 + 2328480.0 * ar5 + 665280.0 * ar7 + 147840.0 * ar9 + 26880.0 * ar11) +
                   xir4 * (56756700.0 * ar + 37837800.0 * ar3 + 15135120.0 * ar5 + 4324320.0 * ar7 + 960960.0 * ar9 + 174720.0 * ar11 +
                           26880.0 * ar13) -
                   xir6 * (113513400.0 * ar + 75675600.0 * ar3 + 30270240.0 * ar5 + 8648640.0 * ar7 + 1921920.0 * ar9 +
                           349440.0 * ar11 + 53760.0 * ar13 + 7168.0 * ar15) +
                   xir8 * (68918850.0 * ar + 45945900.0 * ar3 + 18378360.0 * ar5 + 5250960.0 * ar7 + 1166880.0 * ar9 +
                           212160.0 * ar11 + 32640.0 * ar13 + 4352.0 * ar15 + 512.0 * ar17) -
                   yir2 * (2182950.0 * ar + 1455300.0 * ar3 + 582120.0 * ar5 + 166320.0 * ar7 + 36960.0 * ar9 + 6720.0 * ar11) +
                   xir2 * yir2 *
                       (113513400.0 * ar + 75675600.0 * ar3 + 30270240.0 * ar5 + 8648640.0 * ar7 + 1921920.0 * ar9 + 349440.0 * ar11 +
                        53760.0 * ar13) -
                   xir4 * yir2 *
                       (851350500.0 * ar + 567567000.0 * ar3 + 227026800.0 * ar5 + 64864800.0 * ar7 + 14414400.0 * ar9 +
                        2620800.0 * ar11 + 403200.0 * ar13 + 53760.0 * ar15) +
                   xir6 * yir2 *
                       (1929727800.0 * ar + 1286485200.0 * ar3 + 514594080.0 * ar5 + 147026880.0 * ar7 + 32672640.0 * ar9 +
                        5940480.0 * ar11 + 913920.0 * ar13 + 121856.0 * ar15 + 14336.0 * ar17) -
                   xir8 * yir2 *
                       (1309458150.0 * ar + 872972100.0 * ar3 + 349188840.0 * ar5 + 99768240 * ar7 + 22170720 * ar9 + 4031040 * ar11 +
                        620160 * ar13 + 82688 * ar15 + 9728 * ar17 + 1024.0 * ar19)) /
                      sqrt(M_PI) * exp(-ar2) * r11inv +
                  erfc(ar) *
                      (99225.0 - 4365900.0 * xir2 + 28378350.0 * xir4 - 56756700 * xir6 + 34459425.0 * xir8 - 1091475.0 * yir2 +
                       56756700.0 * xir2 * yir2 - 425675250.0 * xir4 * yir2 + 964863900.0 * xir6 * yir2 - 654729075.0 * xir8 * yir2) *
                      r11inv;

              P10.XXXXXXYYYY +=
                  ((85050.0 * ar + 56700.0 * ar3 + 22680.0 * ar5 + 6480.0 * ar7 + 1440.0 * ar9) -
                   xir2 * (2806650.0 * ar + 1871100.0 * ar3 + 748440.0 * ar5 + 213840.0 * ar7 + 47520.0 * ar9 + 8640.0 * ar11) +
                   xir4 * (12162150.0 * ar + 8108100.0 * ar3 + 3243240.0 * ar5 + 926640.0 * ar7 + 205920.0 * ar9 + 37440.0 * ar11 +
                           5760.0 * ar13) -
                   xir6 * (12162150.0 * ar + 8108100.0 * ar3 + 3243240.0 * ar5 + 926640.0 * ar7 + 205920.0 * ar9 + 37440.0 * ar11 +
                           5760.0 * ar13 + 768.0 * ar15) -
                   yir2 * (1871100.0 * ar + 1247400.0 * ar3 + 498960.0 * ar5 + 142560.0 * ar7 + 31680.0 * ar9 + 5760.0 * ar11) +
                   xir2 * yir2 *
                       (72972900.0 * ar + 48648600 * ar3 + 19459440.0 * ar5 + 5559840.0 * ar7 + 1235520.0 * ar9 + 224640.0 * ar11 +
                        34560.0 * ar13) -
                   xir4 * yir2 *
                       (364864500.0 * ar + 243243000.0 * ar3 + 97297200.0 * ar5 + 27799200.0 * ar7 + 6177600.0 * ar9 +
                        1123200.0 * ar11 + 172800.0 * ar13 + 23040.0 * ar15) +
                   xir6 * yir2 *
                       (413513100.0 * ar + 275675400.0 * ar3 + 110270160.0 * ar5 + 31505760.0 * ar7 + 7001280.0 * ar9 +
                        1272960.0 * ar11 + 195840.0 * ar13 + 26112.0 * ar15 + 3072.0 * ar17) +
                   yir4 * (4054050.0 * ar + 2702700.0 * ar3 + 1081080.0 * ar5 + 308880.0 * ar7 + 68640.0 * ar9 + 12480.0 * ar11 +
                           1920 * ar13) -
                   xir2 * yir4 *
                       (182432250.0 * ar + 121621500.0 * ar3 + 48648600.0 * ar5 + 13899600 * ar7 + 3088800 * ar9 + 561600 * ar11 +
                        86400 * ar13 + 11520 * ar15) +
                   xir4 * yir4 *
                       (1033782750.0 * ar + 689188500.0 * ar3 + 275675400.0 * ar5 + 78764400 * ar7 + 17503200 * ar9 + 3182400 * ar11 +
                        489600 * ar13 + 65280 * ar15 + 7680 * ar17) -
                   xir6 * yir4 *
                       (1309458150.0 * ar + 872972100.0 * ar3 + 349188840.0 * ar5 + 99768240.0 * ar7 + 22170720.0 * ar9 +
                        4031040.0 * ar11 + 620160.0 * ar13 + 82688.0 * ar15 + 9728.0 * ar17 + 1024.0 * ar19)) /
                      sqrt(M_PI) * exp(-ar2) * r11inv +
                  erfc(ar) *
                      (42525.0 - 1403325 * xir2 + 6081075.0 * xir4 - 6081075 * xir6 - 935550.0 * yir2 + 36486450.0 * xir2 * yir2 -
                       182432250.0 * xir4 * yir2 + 206756550.0 * xir6 * yir2 + 2027025.0 * yir4 - 91216125.0 * xir2 * yir4 +
                       516891375.0 * xir4 * yir4 - 654729075.0 * xir6 * yir4) *
                      r11inv;

              P10.XXXXXXYYZZ +=
                  ((28350 * ar + 18900 * ar3 + 7560 * ar5 + 2160 * ar7 + 480 * ar9) -
                   xir2 * (935550 * ar + 623700 * ar3 + 249480 * ar5 + 71280 * ar7 + 15840 * ar9 + 2880 * ar11) +
                   xir4 * (4054050 * ar + 2702700 * ar3 + 1081080 * ar5 + 308880 * ar7 + 68640 * ar9 + 12480 * ar11 + 1920 * ar13) -
                   xir6 * (4054050 * ar + 2702700 * ar3 + 1081080 * ar5 + 308880 * ar7 + 68640 * ar9 + 12480 * ar11 + 1920 * ar13 +
                           256 * ar15) -
                   (yir2 + zir2) * (311850 * ar + 207900 * ar3 + 83160 * ar5 + 23760 * ar7 + 5280 * ar9 + 960 * ar11) +
                   (xir2 * yir2 + xir2 * zir2) *
                       (12162150 * ar + 8108100 * ar3 + 3243240 * ar5 + 926640 * ar7 + 205920 * ar9 + 37440 * ar11 + 5760 * ar13) -
                   (xir4 * yir2 + xir4 * zir2) * (60810750 * ar + 40540500 * ar3 + 16216200 * ar5 + 4633200 * ar7 + 1029600 * ar9 +
                                                  187200 * ar11 + 28800 * ar13 + 3840 * ar15) +
                   (xir6 * yir2 + xir6 * zir2) * (68918850 * ar + 45945900 * ar3 + 18378360 * ar5 + 5250960 * ar7 + 1166880 * ar9 +
                                                  212160 * ar11 + 32640 * ar13 + 4352 * ar15 + 512 * ar17) +
                   yir2 * zir2 *
                       (4054050 * ar + 2702700 * ar3 + 1081080.0 * ar5 + 308880.0 * ar7 + 68640.0 * ar9 + 12480.0 * ar11 +
                        1920 * ar13) -
                   xir2 * yir2 * zir2 *
                       (182432250.0 * ar + 121621500.0 * ar3 + 48648600.0 * ar5 + 13899600 * ar7 + 3088800 * ar9 + 561600 * ar11 +
                        86400 * ar13 + 11520 * ar15) +
                   xir4 * yir2 * zir2 *
                       (1033782750.0 * ar + 689188500.0 * ar3 + 275675400.0 * ar5 + 78764400 * ar7 + 17503200 * ar9 + 3182400 * ar11 +
                        489600 * ar13 + 65280 * ar15 + 7680 * ar17) -
                   xir6 * yir2 * zir2 *
                       (1309458150.0 * ar + 872972100.0 * ar3 + 349188840.0 * ar5 + 99768240.0 * ar7 + 22170720.0 * ar9 +
                        4031040.0 * ar11 + 620160.0 * ar13 + 82688.0 * ar15 + 9728.0 * ar17 + 1024.0 * ar19)) /
                      sqrt(M_PI) * exp(-ar2) * r11inv +
                  erfc(ar) *
                      (14175.0 - 467775 * xir2 + 2027025 * xir4 - 2027025 * xir6 - 155925.0 * (yir2 + zir2) +
                       6081075.0 * xir2 * (yir2 + zir2) - 30405375.0 * xir4 * (yir2 + zir2) + 34459425.0 * xir6 * (yir2 + zir2) +
                       2027025.0 * yir2 * zir2 - 91216125.0 * xir2 * yir2 * zir2 + 516891375.0 * xir4 * yir2 * zir2 -
                       654729075.0 * xir6 * yir2 * zir2) *
                      r11inv;

              P10.XXXXYYYYZZ +=
                  ((17010 * ar + 11340 * ar3 + 4536 * ar5 + 1296 * ar7 + 288 * ar9) -
                   (xir2 + yir2) * (374220 * ar + 249480 * ar3 + 99792 * ar5 + 28512 * ar7 + 6336 * ar9 + 1152 * ar11) +
                   (xir4 + yir4) * (810810 * ar + 540540 * ar3 + 216216 * ar5 + 61776 * ar7 + 13728 * ar9 + 2496 * ar11 + 384 * ar13) +
                   xir2 * yir2 *
                       (9729720 * ar + 6486480 * ar3 + 2594592 * ar5 + 741312 * ar7 + 164736 * ar9 + 29952 * ar11 + 4608 * ar13) -
                   (xir4 * yir2 + xir2 * yir4) * (24324300 * ar + 16216200 * ar3 + 6486480 * ar5 + 1853280 * ar7 + 411840 * ar9 +
                                                  74880 * ar11 + 11520 * ar13 + 1536 * ar15) +
                   xir4 * yir4 *
                       (68918850 * ar + 45945900 * ar3 + 18378360 * ar5 + 5250960 * ar7 + 1166880 * ar9 + 212160 * ar11 +
                        32640 * ar13 + 4352 * ar15 + 512.0 * ar17) -
                   zir2 * (187110 * ar + 124740 * ar3 + 49896 * ar5 + 14256 * ar7 + 3168 * ar9 + 576 * ar11) +
                   (xir2 + yir2) * zir2 *
                       (4864860 * ar + 3243240 * ar3 + 1297296 * ar5 + 370656 * ar7 + 82368 * ar9 + 14976 * ar11 + 2304 * ar13) -
                   (xir4 + yir4) * zir2 *
                       (12162150 * ar + 8108100 * ar3 + 3243240 * ar5 + 926640 * ar7 + 205920 * ar9 + 37440 * ar11 + 5760 * ar13 +
                        768 * ar15) -
                   xir2 * yir2 * zir2 *
                       (145945800 * ar + 97297200 * ar3 + 38918880 * ar5 + 11119680 * ar7 + 2471040 * ar9 + 449280 * ar11 +
                        69120 * ar13 + 9216 * ar15) +
                   (xir4 * yir2 + xir2 * yir4) * zir2 *
                       (413513100 * ar + 275675400 * ar3 + 110270160 * ar5 + 31505760 * ar7 + 7001280 * ar9 + 1272960 * ar11 +
                        195840 * ar13 + 26112 * ar15 + 3072 * ar17) -
                   xir4 * yir4 * zir2 *
                       (1309458150.0 * ar + 872972100.0 * ar3 + 349188840.0 * ar5 + 99768240.0 * ar7 + 22170720.0 * ar9 +
                        4031040.0 * ar11 + 620160.0 * ar13 + 82688.0 * ar15 + 9728.0 * ar17 + 1024.0 * ar19)) /
                      sqrt(M_PI) * exp(-ar2) * r11inv +
                  erfc(ar) *
                      (8505.0 - 187110 * (xir2 + yir2) + 405405 * (xir4 + yir4) + 4864860 * xir2 * yir2 -
                       12162150 * (xir4 * yir2 + xir2 * yir4) + 34459425 * xir4 * yir4 - 93555 * zir2 +
                       2432430 * (xir2 + yir2) * zir2 - 6081075 * (xir4 + yir4) * zir2 - 72972900 * xir2 * yir2 * zir2 +
                       206756550 * (xir4 * yir2 + xir2 * yir4) * zir2 - 654729075.0 * xir4 * yir4 * zir2) *
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
               */

              P10.XXXXXXXXXX += (-60480.0) * pow(alpha, 11) / (11.0 * sqrt(M_PI));
              P10.XXXXXXXXYY += (-6720.0) * pow(alpha, 11) / (11.0 * sqrt(M_PI));
              P10.XXXXXXYYYY += (-2880.0) * pow(alpha, 11) / (11.0 * sqrt(M_PI));
              P10.XXXXXXYYZZ += (-960.0) * pow(alpha, 11) / (11.0 * sqrt(M_PI));
              P10.XXXXYYYYZZ += (-576.0) * pow(alpha, 11) / (11.0 * sqrt(M_PI));
            }
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

              double val = 4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2));

              P10.XXXXXXXXXX += val * pow(kx, 10);
              P10.XXXXXXXXYY += val * pow(kx, 8) * pow(ky, 2);
              P10.XXXXXXYYYY += val * pow(kx, 6) * pow(ky, 4);
              P10.XXXXXXYYZZ += val * pow(kx, 6) * pow(ky, 2) * pow(kz, 2);
              P10.XXXXYYYYZZ += val * pow(kx, 4) * pow(ky, 4) * pow(kz, 2);
            }
        }

  return P10;
}

#endif
