/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file calc_checksum.cc
 *
 *  \brief auxiliary routines to compute MD5 checksums for blocks of memory
 */

#include "gadgetconfig.h"

#ifdef DEBUG_MD5

#include <mpi.h>
#include <stdio.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../debug_md5/Md5.h"
#include "../logs/logs.h"
#include "../main/simulation.h"

void logs::block_checksum(void *base, size_t bytes, int res[4])
{
  MD5_CTX sum;
  union
  {
    unsigned char digest[16];
    int val[4];
  } u;

  MD5Init(&sum);
  MD5UpdateLong(&sum, (unsigned char *)base, bytes);
  MD5Final(&sum);

  for(int i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  for(int i = 0; i < 4; i++)
    res[i] = u.val[i];
}

void logs::calc_memory_checksum(const char *msg, void *base, size_t bytes)
{
  MD5_CTX sum;
  union
  {
    unsigned char digest[16];
    int val[4];
  } u, uglob;

  MD5Init(&sum);
  MD5UpdateLong(&sum, (unsigned char *)base, bytes);
  MD5Final(&sum);

  for(int i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob.val, 4, MPI_INT, MPI_SUM, Communicator);

  if(ThisTask == 0)
    {
      fprintf(FdDebug, "\n");
      fprintf(FdDebug, "Step=%8d  P[]     %s:   ", All.NumCurrentTiStep, msg);
      for(int i = 0; i < 16; i++)
        fprintf(FdDebug, "%02x", uglob.digest[i]);
      fprintf(FdDebug, "\n");
    }
}

void logs::log_debug_md5(const char *msg)
{
  MD5_CTX sum;
  union
  {
    unsigned char digest[16];
    int val[4];
  } u, uglob_P, uglob_SphP;

  MD5Init(&sum);
  MD5UpdateLong(&sum, (unsigned char *)Sp->P, Sp->NumPart * sizeof(particle_data));
  MD5Final(&sum);

  for(int i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob_P.val, 4, MPI_INT, MPI_SUM, Communicator);

  MD5Init(&sum);
  MD5UpdateLong(&sum, (unsigned char *)Sp->SphP, Sp->NumGas * sizeof(sph_particle_data));
  MD5Final(&sum);

  for(int i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob_SphP.val, 4, MPI_INT, MPI_SUM, Communicator);

  if(ThisTask == 0)
    {
      fprintf(FdDebug, "\n");
      fprintf(FdDebug, "Step=%8d  P[]     %s:   ", All.NumCurrentTiStep, msg);
      for(int i = 0; i < 16; i++)
        fprintf(FdDebug, "%02x", uglob_P.digest[i]);
      fprintf(FdDebug, "\n");

      fprintf(FdDebug, "               SphP[]  %s:   ", msg);
      for(int i = 0; i < 16; i++)
        fprintf(FdDebug, "%02x", uglob_SphP.digest[i]);
      fprintf(FdDebug, "\n");

      fflush(FdDebug);
    }
}

#endif
