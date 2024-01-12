/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file parameters.cc
 *
 *  \brief parses the parameter file
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../io/io.h"
#include "../io/parameters.h"
#include "../main/main.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../system/system.h"

void parameters::add_param(const char *name, void *buf, int type, int flag)
{
  if(NParameters > MAX_PARAMETERS)
    Terminate("exceeded MAX_PARAMETERS=%d", MAX_PARAMETERS);

  if(strlen(name) > MAXLEN_PARAM_TAG - 1)
    Terminate("parameter '%s' too long", name);

  strcpy(ParametersTag[NParameters], name);
  ParametersValue[NParameters]      = buf;
  ParametersType[NParameters]       = type;
  ParametersChangeable[NParameters] = flag;
  NParameters++;
}

/*! \brief This function parses the parameter file.
 *
 *  Each parameter is defined by a keyword (`tag'), and can be either
 *  of type douple, int, or character string. Three arrays containing the name,
 *  type and address of the parameter are filled first. The routine then parses
 *  the parameter file and fills the referenced variables. The routine makes sure that
 *  each parameter appears exactly once in the parameter file, otherwise
 *  error messages are produced that complain about the missing parameters.
 *  Basic checks are performed on the supplied parameters in the end.
 *
 *  \param fname The file name of the parameter file
 */
int parameters::read_parameter_file(const char *fname)
{
  FILE *fd, *fdout;
  int bufsize = MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200;
  char buf[bufsize];
  int param_handled[MAX_PARAMETERS];
  int errorFlag = 0;

  for(int i = 0; i < MAX_PARAMETERS; i++)
    {
      param_handled[i]     = 0;
      ParameterSequence[i] = -1;
    }

  if(sizeof(long long) != 8)
    Terminate("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");

  if(sizeof(int) != 4)
    Terminate("\nType `int' is not 32 bit on this platform. Stopping.\n\n");

  if(sizeof(float) != 4)
    Terminate("\nType `float' is not 32 bit on this platform. Stopping.\n\n");

  if(sizeof(double) != 8)
    Terminate("\nType `double' is not 64 bit on this platform. Stopping.\n\n");

  if(ThisTask == 0) /* read parameter file on process 0 */
    {
      if((fd = fopen(fname, "r")))
        {
          snprintf(buf, bufsize, "%s%s", fname, "-usedvalues");
          if(!(fdout = fopen(buf, "w")))
            {
              printf("error opening file '%s' \n", buf);
              errorFlag = 1;
            }
          else
            {
              printf("Obtaining parameters from file '%s':\n\n", fname);
              int cnt = 0;
              while(!feof(fd))
                {
                  int bufsize1 = MAXLEN_PARAM_TAG + 200;
                  int bufsize2 = MAXLEN_PARAM_VALUE + 200;
                  int bufsize3 = MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 400;
                  char buf1[bufsize1], buf2[bufsize2], buf3[bufsize3];

                  *buf = 0;
                  fgets(buf, MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200, fd);
                  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                    continue;

                  if(buf1[0] == '%' || buf1[0] == '#')
                    continue;

                  int j = -1;
                  for(int i = 0; i < NParameters; i++)
                    if(strcmp(buf1, ParametersTag[i]) == 0)
                      {
                        if(param_handled[i] == 0)
                          {
                            j                        = i;
                            param_handled[i]         = 1;
                            ParameterSequence[cnt++] = i;
                            break;
                          }
                        else
                          {
                            j = -2;
                            break;
                          }
                      }

                  if(j >= 0)
                    {
                      switch(ParametersType[j])
                        {
                          case PARAM_DOUBLE:
                            *((double *)ParametersValue[j]) = atof(buf2);
                            snprintf(buf3, bufsize3, "%%-%ds%%g\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, *((double *)ParametersValue[j]));
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, *((double *)ParametersValue[j]));
                            break;
                          case PARAM_STRING:
                            if(strcmp(buf2, "OUTPUT_DIR") == 0)
                              {
                                if(getenv("OUTPUT_DIR"))
                                  strcpy(buf2, getenv("OUTPUT_DIR"));
                                else
                                  Terminate("no environment variable OUTPUT_DIR found");
                              }
                            strcpy((char *)ParametersValue[j], buf2);
                            snprintf(buf3, bufsize3, "%%-%ds%%s\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, buf2);
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, buf2);
                            break;
                          case PARAM_INT:
                            *((int *)ParametersValue[j]) = atoi(buf2);
                            snprintf(buf3, bufsize3, "%%-%ds%%d\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, *((int *)ParametersValue[j]));
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, *((int *)ParametersValue[j]));
                            break;
                        }
                    }
                  else if(j == -2)
                    {
                      fprintf(stdout, "Error in file %s:   Tag '%s' multiply defined.\n", fname, buf1);
                      errorFlag = 1;
                    }
                  else
                    {
                      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed\n", fname, buf1);
                      errorFlag = 1;
                    }
                }
              fclose(fd);
              fclose(fdout);
              printf("\n");
            }
        }
      else
        {
          printf("Parameter file %s not found.\n", fname);
          errorFlag = 1;
        }

      for(int i = 0; i < NParameters; i++)
        {
          if(param_handled[i] != 1)
            {
              printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", ParametersTag[i], fname);
              errorFlag = 1;
            }
        }
    }

  return errorFlag;
}

void parameters::write_used_parameters(const char *dirname, const char *fname)
{
  if(ThisTask == 0)
    {
      mkdir(dirname, 02755);
      char buf[MAXLEN_PATH_EXTRA];
      snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", dirname, fname);
      FILE *fdout = fopen(buf, "w");
      if(!fdout)
        Terminate("Can't open file '%s'", buf);

      for(int i = 0; i < NParameters; i++)
        {
          int j = ParameterSequence[i];

          if(j >= 0)
            {
              int bufsize3 = MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 400;
              char buf3[bufsize3];

              switch(ParametersType[j])
                {
                  case PARAM_DOUBLE:
                    snprintf(buf3, bufsize3, "%%-%ds%%g\n", MAXLEN_PARAM_TAG);
                    fprintf(fdout, buf3, ParametersTag[j], *((double *)ParametersValue[j]));
                    break;
                  case PARAM_STRING:
                    snprintf(buf3, bufsize3, "%%-%ds%%s\n", MAXLEN_PARAM_TAG);
                    fprintf(fdout, buf3, ParametersTag[j], (char *)ParametersValue[j]);
                    break;
                  case PARAM_INT:
                    snprintf(buf3, bufsize3, "%%-%ds%%d\n", MAXLEN_PARAM_TAG);
                    fprintf(fdout, buf3, ParametersTag[j], *((int *)ParametersValue[j]));
                    break;
                }
            }
        }

      fclose(fdout);
    }
}
