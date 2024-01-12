/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  system.h
 *
 *  \brief declares functions for various low level helper routines
 */

#ifndef SYSTEM_H
#define SYSTEM_H

#include "gadgetconfig.h"

#include <gsl/gsl_rng.h>
#include <stdio.h>

extern gsl_rng *random_generator; /**< the random number generator used */

void myflush(FILE *fstream);

void enable_core_dumps_and_fpu_exceptions(void);

void permutate_chunks_in_list(int ncount, int *list);

void subdivide_evenly(long long N, int pieces, int index_bin, long long *first, int *count);
void subdivide_evenly(int N, int pieces, int index, int *first, int *count);
void subdivide_evenly_get_bin(int N, int pieces, int index, int *bin);

void init_rng(int thistask);
double get_random_number(void);

int my_fls(unsigned int x);

template <typename T>
inline T square(T const value)
{
  return value * value;
}

#endif
