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

#include "gadgetconfig.h"

#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "../system/system.h"

#define TERMINATE_STATUS EXIT_FAILURE
#define TERMINATE_MSG "TERMINATE: ******!!!!!******  Code termination on task=%d, function %s(), file %s, line %d: "
#define Terminate(...)                                               \
  do                                                                 \
    {                                                                \
      int thistask;                                                  \
      MPI_Comm_rank(MPI_COMM_WORLD, &thistask);                      \
      printf(TERMINATE_MSG, thistask, __func__, __FILE__, __LINE__); \
      printf(__VA_ARGS__);                                           \
      printf("\n");                                                  \
      fflush(stdout);                                                \
      MPI_Abort(MPI_COMM_WORLD, TERMINATE_STATUS);                   \
      exit(TERMINATE_STATUS);                                        \
    }                                                                \
  while(0)

#define WARNINGS_FILE_NAME "WARNINGS"
#define WARN_MSG "WARNING: Code warning on task=%d, function %s(), file %s, line %d: "
#define warn(...)                                                         \
  do                                                                      \
    {                                                                     \
      int thistask;                                                       \
      MPI_Comm_rank(MPI_COMM_WORLD, &thistask);                           \
      printf(WARN_MSG, thistask, __func__, __FILE__, __LINE__);           \
      printf(__VA_ARGS__);                                                \
      printf("\n");                                                       \
      myflush(stdout);                                                    \
      FILE *const warn_fd = fopen(WARNINGS_FILE_NAME, "a");               \
      fprintf(warn_fd, WARN_MSG, thistask, __func__, __FILE__, __LINE__); \
      fprintf(warn_fd, __VA_ARGS__);                                      \
      fprintf(warn_fd, "\n");                                             \
      fclose(warn_fd);                                                    \
    }                                                                     \
  while(0)

/* define an "assert" macro which outputs MPI rank (we do NOT want to call
 * MPI_Abort, because then the assertion failure isn't caught in the debugger) */
#define ASSERT_MSG "Assertion failure!\n\ttask=%d, function %s(), file %s, line %d:\n\t%s\n"
#ifdef NDEBUG
#define myassert(cond)
#else
#define myassert(cond)                                                       \
  do                                                                         \
    {                                                                        \
      if(!(cond))                                                            \
        {                                                                    \
          int thistask;                                                      \
          MPI_Comm_rank(MPI_COMM_WORLD, &thistask);                          \
          printf(ASSERT_MSG, thistask, __func__, __FILE__, __LINE__, #cond); \
          myflush(stdout);                                                   \
          assert(0);                                                         \
        }                                                                    \
    }                                                                        \
  while(0)
#endif

#endif
