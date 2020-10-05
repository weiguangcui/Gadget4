/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  mpi_types.cc
 *
 *  \brief implements some user defined MPI types for collectives
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../mpi_utils/mpi_utils.h"

static void min_MPI_MyIntPosType(void *in, void *inout, int *len, MPI_Datatype *type)
{
  /* we here trust that this is called only for the correct type  */

  MyIntPosType *invalues    = (MyIntPosType *)in;
  MyIntPosType *inoutvalues = (MyIntPosType *)inout;

  for(int i = 0; i < *len; i++)
    if(invalues[i] < inoutvalues[i])
      inoutvalues[i] = invalues[i];
}

static void max_MPI_MyIntPosType(void *in, void *inout, int *len, MPI_Datatype *type)
{
  /* we here trust that this is called only for the correct type  */

  MyIntPosType *invalues    = (MyIntPosType *)in;
  MyIntPosType *inoutvalues = (MyIntPosType *)inout;

  for(int i = 0; i < *len; i++)
    if(invalues[i] > inoutvalues[i])
      inoutvalues[i] = invalues[i];
}

static void min_MPI_MySignedIntPosType(void *in, void *inout, int *len, MPI_Datatype *type)
{
  /* we here trust that this is called only for the correct type  */

  MySignedIntPosType *invalues    = (MySignedIntPosType *)in;
  MySignedIntPosType *inoutvalues = (MySignedIntPosType *)inout;

  for(int i = 0; i < *len; i++)
    if(invalues[i] < inoutvalues[i])
      inoutvalues[i] = invalues[i];
}

static void max_MPI_MySignedIntPosType(void *in, void *inout, int *len, MPI_Datatype *type)
{
  /* we here trust that this is called only for the correct type  */

  MySignedIntPosType *invalues    = (MySignedIntPosType *)in;
  MySignedIntPosType *inoutvalues = (MySignedIntPosType *)inout;

  for(int i = 0; i < *len; i++)
    if(invalues[i] > inoutvalues[i])
      inoutvalues[i] = invalues[i];
}

void my_mpi_types_init(void)
{
  /* create our new data type */
  MPI_Type_contiguous(sizeof(MyIntPosType), MPI_BYTE, &MPI_MyIntPosType);
  MPI_Type_commit(&MPI_MyIntPosType);

  /* create our operators */
  MPI_Op_create(min_MPI_MyIntPosType, 1, &MPI_MIN_MyIntPosType);
  MPI_Op_create(max_MPI_MyIntPosType, 1, &MPI_MAX_MyIntPosType);
  MPI_Op_create(min_MPI_MySignedIntPosType, 1, &MPI_MIN_MySignedIntPosType);
  MPI_Op_create(max_MPI_MySignedIntPosType, 1, &MPI_MAX_MySignedIntPosType);
}
