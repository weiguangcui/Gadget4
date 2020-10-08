/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file macros.h
 *
 *  \brief defines macros for run termination or for warnings
 */

#ifndef MACROS_H
#define MACROS_H

#ifdef MPI_HYPERCUBE_ALLGATHERV
#define MPI_Allgatherv MPI_hypercube_Allgatherv
#endif

#define Terminate(...)                                                                                                      \
  {                                                                                                                         \
    {                                                                                                                       \
      char termbuf1__[8000], termbuf2__[8000];                                                                              \
      int thistask;                                                                                                         \
      MPI_Comm_rank(MPI_COMM_WORLD, &thistask);                                                                             \
      sprintf(termbuf1__, "Code termination on task=%d, function %s(), file %s, line %d", thistask, __FUNCTION__, __FILE__, \
              __LINE__);                                                                                                    \
      sprintf(termbuf2__, __VA_ARGS__);                                                                                     \
      printf("%s: %s\n", termbuf1__, termbuf2__);                                                                           \
      fflush(stdout);                                                                                                       \
      MPI_Abort(MPI_COMM_WORLD, 1);                                                                                         \
    }                                                                                                                       \
    exit(0);                                                                                                                \
  }
#define warn(...)                                                                                                                \
  {                                                                                                                              \
    char termbuf1__[8000], termbuf2__[8000];                                                                                     \
    int thistask;                                                                                                                \
    MPI_Comm_rank(MPI_COMM_WORLD, &thistask);                                                                                    \
    sprintf(termbuf1__, "Code warning on task=%d, function %s(), file %s, line %d", thistask, __FUNCTION__, __FILE__, __LINE__); \
    sprintf(termbuf2__, __VA_ARGS__);                                                                                            \
    printf("%s: %s\n", termbuf1__, termbuf2__);                                                                                  \
    myflush(stdout);                                                                                                             \
    FILE *fd__ = fopen("WARNINGS", "w");                                                                                         \
    fclose(fd__);                                                                                                                \
  }

#endif
