/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file test_symtensors.cc
 *
 *  \brief some test routines for the symmetric tensor implementation
 */

#include "gadgetconfig.h"

#ifdef DEBUG_SYMTENSORS

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/symtensors.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"

static bool compare_list(const int &a, const int &b) { return a < b; }

static void symtensor_test_tensor4_contraction_with_tensor1(void)
{
  symtensor4<double> T4;
  vector<double> T1;

  int I4[3][3][3][3];
  int I3[3][3][3];
  int I1[3] = {0, 1, 2};

  int count = 0;

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
          {
            if(i <= j && j <= k && k <= l)
              {
                I4[i][j][k][l] = count++;
              }
            else
              {
                int li[] = {i, j, k, l};

                mycxxsort(li, li + 4, compare_list);

                I4[i][j][k][l] = I4[li[0]][li[1]][li[2]][li[3]];
              }

            // printf("%c%c%c%c:  %d\n", 'X' + i, 'X' + j, 'X' + k, 'X' + l, I4[i][j][k][l]);
          }

  count = 0;

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        {
          if(i <= j && j <= k)
            {
              I3[i][j][k] = count++;
            }
          else
            {
              int li[] = {i, j, k};

              mycxxsort(li, li + 3, compare_list);

              I3[i][j][k] = I3[li[0]][li[1]][li[2]];
            }

          // printf("%c%c:  %d\n", 'X' + i, 'X' + j, I2[i][j]);
        }

  for(int i = 0; i < 15; i++)
    T4[i] = get_random_number();

  for(int i = 0; i < 3; i++)
    T1[i] = get_random_number();

  /* now fill in the whole general matrix based on the symmetric one */
  double M4[3][3][3][3];
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
          M4[i][j][k][l] = T4[I4[i][j][k][l]];

  /* now fill in the whole general matrix based on the symmetric one */
  double M1[3];
  for(int i = 0; i < 3; i++)
    M1[i] = T1[I1[i]];

  /* result of full matrix reduction */
  double R3[3][3][3];
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        {
          R3[i][j][k] = 0;

          for(int l = 0; l < 3; l++)
            R3[i][j][k] += M4[i][j][k][l] * M1[l];
        }

  /* now let's compare the result */

  symtensor3<double> S3 = T4 * T1;  // reduction of symmetric 4-tensor with symmetric 2-tensor

  printf("Result comparison:\n");

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        {
          if(i <= j && j <= k)
            {
              printf("%c%c%c:  %g   %g\n", 'X' + i, 'X' + j, 'X' + k, S3[I3[i][j][k]], R3[i][j][k]);
            }
        }

  printf("\n");
}

static void symtensor_test_tensor4_contraction_with_tensor2(void)
{
  symtensor4<double> T4;
  symtensor2<double> T2;

  int I4[3][3][3][3];
  int I2[3][3];

  int count = 0;

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
          {
            if(i <= j && j <= k && k <= l)
              {
                I4[i][j][k][l] = count++;
              }
            else
              {
                int li[] = {i, j, k, l};

                mycxxsort(li, li + 4, compare_list);

                I4[i][j][k][l] = I4[li[0]][li[1]][li[2]][li[3]];
              }

            // printf("%c%c%c%c:  %d\n", 'X' + i, 'X' + j, 'X' + k, 'X' + l, I4[i][j][k][l]);
          }

  // printf("count=%d\n\n", count);

  count = 0;

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      {
        if(i <= j)
          {
            I2[i][j] = count++;
          }
        else
          {
            int li[] = {i, j};

            mycxxsort(li, li + 2, compare_list);

            I2[i][j] = I2[li[0]][li[1]];
          }

        // printf("%c%c:  %d\n", 'X' + i, 'X' + j, I2[i][j]);
      }

  // printf("count=%d\n\n", count);

  for(int i = 0; i < 15; i++)
    T4[i] = get_random_number();

  for(int i = 0; i < 6; i++)
    T2[i] = get_random_number();

  /* now fill in the whole general matrix based on the symmetric one */
  double M4[3][3][3][3];
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
          M4[i][j][k][l] = T4[I4[i][j][k][l]];

  /* now fill in the whole general matrix based on the symmetric one */
  double M2[3][3];
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      M2[i][j] = T2[I2[i][j]];

  /* result of full matrix reduction */
  double R2[3][3];
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      {
        R2[i][j] = 0;

        for(int k = 0; k < 3; k++)
          for(int l = 0; l < 3; l++)
            R2[i][j] += M4[i][j][k][l] * M2[k][l];
      }

  /* now let's compare the result */

  symtensor2<double> S2 = T4 * T2;  // reduction of symmetric 4-tensor with symmetric 2-tensor

  printf("Result comparison:\n");

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      {
        if(i <= j)
          {
            printf("%c%c:  %g   %g\n", 'X' + i, 'X' + j, S2[I2[i][j]], R2[i][j]);
          }
      }

  printf("\n");
}

static void symtensor_test_tensor4_contraction_with_tensor3(void)
{
  symtensor4<double> T4;
  symtensor3<double> T3;

  int I4[3][3][3][3];
  int I3[3][3][3];

  int count = 0;

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
          {
            if(i <= j && j <= k && k <= l)
              {
                I4[i][j][k][l] = count++;
              }
            else
              {
                int li[] = {i, j, k, l};

                mycxxsort(li, li + 4, compare_list);

                I4[i][j][k][l] = I4[li[0]][li[1]][li[2]][li[3]];
              }

            // printf("%c%c%c%c:  %d\n", 'X' + i, 'X' + j, 'X' + k, 'X' + l, I4[i][j][k][l]);
          }

  // printf("count=%d\n\n", count);

  count = 0;

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        {
          if(i <= j && j <= k)
            {
              I3[i][j][k] = count++;
            }
          else
            {
              int li[] = {i, j, k};

              mycxxsort(li, li + 3, compare_list);

              I3[i][j][k] = I3[li[0]][li[1]][li[2]];
            }

          // printf("%c%c:  %d\n", 'X' + i, 'X' + j, I2[i][j]);
        }

  // printf("count=%d\n\n", count);

  for(int i = 0; i < 15; i++)
    T4[i] = get_random_number();

  for(int i = 0; i < 10; i++)
    T3[i] = get_random_number();

  /* now fill in the whole general matrix based on the symmetric one */
  double M4[3][3][3][3];
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
          M4[i][j][k][l] = T4[I4[i][j][k][l]];

  /* now fill in the whole general matrix based on the symmetric one */
  double M3[3][3][3];
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
        M3[i][j][k] = T3[I3[i][j][k]];

  /* result of full matrix reduction */
  double R1[3];
  for(int i = 0; i < 3; i++)
    {
      R1[i] = 0;

      for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
          for(int l = 0; l < 3; l++)
            R1[i] += M4[i][j][k][l] * M3[j][k][l];
    }

  /* now let's compare the result */

  vector<double> S1 = T4 * T3;  // reduction of symmetric 4-tensor with symmetric 2-tensor

  printf("Result comparison:\n");

  for(int i = 0; i < 3; i++)
    printf("%c:  %g   %g\n", 'X' + i, S1[i], R1[i]);

  printf("\n");
}

void symtensor_test(void)
{
  symtensor_test_tensor4_contraction_with_tensor1();
  symtensor_test_tensor4_contraction_with_tensor2();
  symtensor_test_tensor4_contraction_with_tensor3();

  Terminate("Done with test");
}

#endif
