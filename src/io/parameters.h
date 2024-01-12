/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file parameters.h
 *
 *  \brief declares a class for dealing with the parameter file
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "gadgetconfig.h"

#include "../data/dtypes.h"
#include "../mpi_utils/setcomm.h"

#define PARAM_DOUBLE 1
#define PARAM_STRING 2
#define PARAM_INT 3

#define PARAM_FIXED 0
#define PARAM_CHANGEABLE 1

#define MAXLEN_PARAM_TAG 50    /**< maximum length of the tag of a parameter in the parameter file */
#define MAXLEN_PARAM_VALUE 200 /**< maximum length of the value of a parameter in the parameter file */
#define MAX_PARAMETERS 300     /**< maximum number of parameters in the parameter file */

class parameters : public setcomm
{
 public:
  // constructors
  parameters() : setcomm("delayed init") {}
  parameters(MPI_Comm comm) : setcomm(comm) {}

  int read_parameter_file(const char *fname);

  void add_param(const char *name, void *buf, int type, int flag);

  void write_used_parameters(const char *dirname, const char *fname);

  int NParameters = 0;

  char ParametersTag[MAX_PARAMETERS][MAXLEN_PARAM_TAG];
  void *ParametersValue[MAX_PARAMETERS];
  char ParametersType[MAX_PARAMETERS];
  char ParametersChangeable[MAX_PARAMETERS];
  int ParameterSequence[MAX_PARAMETERS];
};

#endif /* PARAMETERS_H */
