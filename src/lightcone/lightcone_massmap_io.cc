/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file lightcone_massmap_io.cc
 *
 * \brief routines for I/O of lightcone mass projections
 */

#include "gadgetconfig.h"

#if defined(LIGHTCONE) && defined(LIGHTCONE_MASSMAPS)

#include <gsl/gsl_rng.h>
#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../lightcone/lightcone.h"
#include "../lightcone/lightcone_massmap_io.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../system/system.h"

/*!
 * \brief Function for field registering.
 *
 * For init_field( arguments read the documentation of init_field(.
 * Don't forget to add the new IO_FLAG to io_private.h
 */

lightcone_massmap_io::lightcone_massmap_io(mmparticles *Mp_ptr, lightcone *LightCone_ptr, MPI_Comm comm, int format)
    : IO_Def(comm, format)
{
  Mp        = Mp_ptr;
  LightCone = LightCone_ptr;

  this->N_IO_Fields  = 0;
  this->N_DataGroups = 1;
  this->header_size  = sizeof(header);
  this->header_buf   = &header;
  this->type_of_file = FILE_IS_MASSMAP;
  snprintf(this->info, MAXLEN_PATH, "LIGHTCONE: writing mass map data");

  init_field("MAMP", "Mass", MEM_DOUBLE, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_MM, &LightCone->MassMap[0], NULL, MASSMAPS, 1, 0., -1.,
             0., 1., 0., All.UnitMass_in_g, true);
}

void lightcone_massmap_io::lightcone_massmap_save(int num)
{
  char buf[MAXLEN_PATH_EXTRA];

  selected_bnd = num;

  long long NumLP_tot = LightCone->Mp->NpixLoc;
  MPI_Allreduce(MPI_IN_PLACE, &NumLP_tot, 1, MPI_LONG_LONG, MPI_SUM, Communicator);
  mpi_printf("\nLIGHTCONE: writing lightcone massmap files #%d ... (Npix_tot = %lld,  ascale %g to %g)\n", num, NumLP_tot,
             LightCone->MassMapBoundariesAscale[num], LightCone->MassMapBoundariesAscale[num + 1]);

  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          char buf[MAXLEN_PATH_EXTRA];
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s/mapsdir_%03d", All.OutputDir, num);
          mkdir(buf, 02755);
        }
      MPI_Barrier(Communicator);
    }

  if(All.NumFilesPerSnapshot > 1)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/mapsdir_%03d/%s_%03d", All.OutputDir, num, "maps", num);
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/%s_%03d", All.OutputDir, "maps", num);

  write_multiple_files(buf, All.NumFilesPerSnapshot);

  mpi_printf("LIGHTCONE: done with writing mass map.\n");

  /* clear the massmap again */
  memset(LightCone->MassMap, 0, LightCone->Mp->NpixLoc * sizeof(double));
}

void lightcone_massmap_io::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type)
{
  /* determine global and local pixel numbers */

  n_type[0] = LightCone->Mp->NpixLoc;

  /* determine particle numbers of each type in file */
  if(ThisTask == writeTask)
    {
      ntot_type[0] = n_type[0];

      for(int task = writeTask + 1; task <= lastTask; task++)
        {
          long long nn;
          MPI_Recv(&nn, 1, MPI_LONG_LONG, task, TAG_LOCALN, Communicator, MPI_STATUS_IGNORE);
          ntot_type[0] += nn;
        }

      for(int task = writeTask + 1; task <= lastTask; task++)
        MPI_Send(&ntot_type[0], 1, MPI_LONG_LONG, task, TAG_N, Communicator);
    }
  else
    {
      MPI_Send(&n_type[0], 1, MPI_LONG_LONG, writeTask, TAG_LOCALN, Communicator);
      MPI_Recv(&ntot_type[0], 1, MPI_LONG_LONG, writeTask, TAG_N, Communicator, MPI_STATUS_IGNORE);
    }

  /* fill file header */
  header.nside      = All.LightConeMassMapsNside; /* healpix nside */
  header.npix_local = ntot_type[0];
  header.npix_total = LightCone->Mp->Npix;
  header.num_files  = All.NumFilesPerSnapshot;

  header.AscaleStart  = LightCone->MassMapBoundariesAscale[selected_bnd];
  header.AscaleEnd    = LightCone->MassMapBoundariesAscale[selected_bnd + 1];
  header.ComDistStart = LightCone->MassMapBoundariesComDist[selected_bnd];
  header.ComDistEnd   = LightCone->MassMapBoundariesComDist[selected_bnd + 1];
}

void lightcone_massmap_io::write_header_fields(hid_t handle)
{
  write_scalar_attribute(handle, "Nside", &header.nside, H5T_NATIVE_INT);

  write_scalar_attribute(handle, "NpixLocal", &header.npix_local, H5T_NATIVE_INT);

  write_scalar_attribute(handle, "NpixTotal", &header.npix_total, H5T_NATIVE_INT);

  write_scalar_attribute(handle, "NumFiles", &header.num_files, H5T_NATIVE_INT);

  write_scalar_attribute(handle, "AscaleStart", &header.AscaleStart, H5T_NATIVE_DOUBLE);

  write_scalar_attribute(handle, "AscaleEnd", &header.AscaleEnd, H5T_NATIVE_DOUBLE);

  write_scalar_attribute(handle, "ComDistStart", &header.ComDistStart, H5T_NATIVE_DOUBLE);

  write_scalar_attribute(handle, "ComDistEnd", &header.ComDistEnd, H5T_NATIVE_DOUBLE);
}

void lightcone_massmap_io::set_filenr_in_header(int numfiles) { header.num_files = numfiles; }

void lightcone_massmap_io::get_datagroup_name(int type, char *buf)
{
  if(type == 0)
    snprintf(buf, MAXLEN_PATH, "/Maps");
  else
    Terminate("should not get here");
}

void *lightcone_massmap_io::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_MM:
        return (void *)(LightCone->MassMap + index);

      default:
        Terminate("we don't expect to get here");
    }

  return NULL;
}

void lightcone_massmap_io::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *nloc_part,
                                            long long *npart, int *nstart)
{
  /* empty */
}

void lightcone_massmap_io::read_header_fields(const char *fname)
{ /* empty */
}

void lightcone_massmap_io::read_increase_numbers(int type, int n_for_this_task)
{ /* empty */
}

int lightcone_massmap_io::get_filenr_from_header(void)
{
  /* empty */
  return 0;
}

int lightcone_massmap_io::get_type_of_element(int index)
{
  /* empty */
  return 0;
}

void lightcone_massmap_io::set_type_of_element(int index, int type)
{ /* empty */
}

#endif
