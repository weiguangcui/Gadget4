/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file lightcone_particle_io.cc
 *
 * \brief routines for I/O of lightcone particles
 */

#include "gadgetconfig.h"

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES)

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
#include "../lightcone/lightcone_particle_io.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../system/system.h"

/*!
 * \brief Function for field registering.
 *
 * For init_field( arguments read the documentation of init_field(.
 * Don't forget to add the new IO_FLAG to io_private.h
 */
#ifdef MERGERTREE
lightcone_particle_io::lightcone_particle_io(lcparticles *Lp_ptr, lightcone *LightCone_ptr, mergertree *MergerTree_ptr, MPI_Comm comm,
                                             int format)
    : IO_Def(comm, format)
#else
lightcone_particle_io::lightcone_particle_io(lcparticles *Lp_ptr, lightcone *LightCone_ptr, MPI_Comm comm, int format)
    : IO_Def(comm, format)
#endif
{
  Lp        = Lp_ptr;
  LightCone = LightCone_ptr;
#ifdef MERGERTREE
  MergerTree = MergerTree_ptr;
#endif

  this->N_IO_Fields  = 0;
  this->N_DataGroups = NTYPES + 2;
  /* the +1 data group is a tree table used only for storing reordered lightcone data,
     the +2 data group is a coarse healpix table used only for storing normal lightcone data */

  this->header_size  = sizeof(header);
  this->header_buf   = &header;
  this->type_of_file = FILE_IS_LIGHTCONE;
  snprintf(this->info, MAXLEN_PATH, "LIGHTCONE: writing particle lightcone data");

  init_field("POS ", "Coordinates", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_LC, NULL, io_func_pos, ALL_TYPES, 1, 1.,
             -1., 1., 0., 0., All.UnitLength_in_cm, true);

#ifdef OUTPUT_VELOCITIES_IN_HALF_PRECISION
  init_field("VEL ", "Velocities", MEM_MY_FLOAT, FILE_HALF, READ_IF_PRESENT, 3, A_LC, &Lp->P[0].Vel[0], NULL, ALL_TYPES, 1, 0.5, 0.,
             0., 0., 1., All.UnitVelocity_in_cm_per_s);
#else
  init_field("VEL ", "Velocities", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_LC, &Lp->P[0].Vel[0], NULL, ALL_TYPES, 1, 0.5,
             0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
#endif

#ifdef LIGHTCONE_OUTPUT_ACCELERATIONS
#ifdef OUTPUT_ACCELERATIONS_IN_HALF_PRECISION
  All.accel_normalize_fac = 10.0 * All.Hubble * (100.0 * 1.0e5 / All.UnitVelocity_in_cm_per_s);
  init_field("ACCE", "Acceleration", MEM_MY_FLOAT, FILE_HALF, READ_IF_PRESENT, 3, A_LC, NULL, io_func_accel, ALL_TYPES, 1, -2.0, 1, -1,
             0, 2, All.accel_normalize_fac * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#else
  init_field("ACCE", "Acceleration", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_LC, &Lp->P[0].GravAccel[0], NULL, ALL_TYPES,
             1, -2.0, 1, -1, 0, 2, All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#endif
#endif

  init_field("ID  ", "ParticleIDs", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, READ_IF_PRESENT, 1, A_LC, NULL, io_func_id, ALL_TYPES, 0, 0, 0, 0,
             0, 0, 0, true);

#ifndef LEAN
  init_field("MASS", "Masses", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_LC, NULL, io_func_mass, MASS_BLOCK, 1, 0., -1.,
             0., 1., 0., All.UnitMass_in_g);
#endif

#ifndef LEAN
  init_field("ASCL", "Ascale", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_LC, &Lp->P[0].Ascale, NULL, ALL_TYPES, 0, 0, 0, 0,
             0, 0, 0);
#endif

#ifdef REARRANGE_OPTION
  init_field("MTRI", "TreeID", MEM_INT64, FILE_INT64, SKIP_ON_READ, 1, A_LC, &Lp->P[0].TreeID, NULL, ALL_TYPES, 0, 0, 0, 0, 0, 0, 0);
#endif

#if defined(SUBFIND) && defined(SUBFIND_STORE_LOCAL_DENSITY)
  init_field("SFDE", "SubfindDensity", MEM_MY_FLOAT, All.RestartFlag != RST_CREATEICS ? FILE_MY_IO_FLOAT : FILE_NONE, SKIP_ON_READ, 1,
             A_PS, &Lp->PS[0].SubfindDensity, 0, ALL_TYPES, /* subfind density */
             1, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);

  init_field("SFHS", "SubfindHsml", MEM_MY_FLOAT, All.RestartFlag != RST_CREATEICS ? FILE_MY_IO_FLOAT : FILE_NONE, SKIP_ON_READ, 1,
             A_PS, &Lp->PS[0].SubfindHsml, 0, ALL_TYPES, /* subfind hsml */
             1, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

  init_field("SFVD", "SubfindVelDisp", MEM_MY_FLOAT, All.RestartFlag != RST_CREATEICS ? FILE_MY_IO_FLOAT : FILE_NONE, SKIP_ON_READ, 1,
             A_PS, &Lp->PS[0].SubfindVelDisp, 0, ALL_TYPES, /* subfind velocity dispersion */
             1, 0., 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
#endif

#ifdef MERGERTREE
  init_field("MTRL", "ParticleCount", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_TT, &MergerTree->TreeTable[0].HaloCount, NULL, TREETABLE,
             0, 0, 0, 0, 0, 0, 0, true);
  init_field("MTRS", "ParticleFirst", MEM_INT64, FILE_INT64, SKIP_ON_READ, 1, A_TT, &MergerTree->TreeTable[0].FirstHalo, NULL,
             TREETABLE, 0, 0, 0, 0, 0, 0, 0, true);
  init_field("MTRI", "TreeID", MEM_INT64, FILE_INT64, SKIP_ON_READ, 1, A_TT, &MergerTree->TreeTable[0].TreeID, NULL, TREETABLE, 0, 0,
             0, 0, 0, 0, 0, true);
#endif

  init_field("HPHT", "ParticleCount", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_MM, &Lp->HealPixTab_PartCount[0], NULL, HEALPIXTAB, 0, 0,
             0, 0, 0, 0, 0, true);
}

void lightcone_particle_io::lightcone_read(int num, int conenr)
{
  double t0 = Logs.second();

  Lp->TotNumPart = 0;

  char fname[MAXLEN_PATH_EXTRA];

  if(All.NumFilesPerSnapshot > 1)
    snprintf(fname, MAXLEN_PATH_EXTRA, "%s/lightcone_%02d/conedir_%04d/%s_%04d", All.OutputDir, conenr, num, "conesnap", num);
  else
    snprintf(fname, MAXLEN_PATH_EXTRA, "%s/lightcone_%02d/%s_%04d", All.OutputDir, conenr, "conesnap", num);

  int num_files = find_files(fname, fname);

  reset_io_byte_count();

  for(int rep = 0; rep < 2; rep++)
    {
      Lp->NumPart = 0;

      read_files_driver(fname, rep, num_files);

      /* now do the memory allocation */
      if(rep == 0)
        {
          Lp->MaxPart = Lp->NumPart;
          Lp->allocate_memory(); /* allocate lightcone particle storage */
        }
    }

  MPI_Barrier(Communicator);

  long long byte_count = get_io_byte_count(), byte_count_all;
  sumup_longs(1, &byte_count, &byte_count_all, Communicator);

  double t1 = Logs.second();

  mpi_printf("LIGHTCONE-READ: reading done. Took %g sec, total size %g MB, corresponds to effective I/O rate of %g MB/sec\n",
             Logs.timediff(t0, t1), byte_count_all / (1024.0 * 1024.0), byte_count_all / (1024.0 * 1024.0) / Logs.timediff(t0, t1));

  mpi_printf("\nLIGHTCONE-READ: Total number of particles :  %lld\n\n", Lp->TotNumPart);
}

void lightcone_particle_io::lightcone_save(int num, int conenr, bool reordered_flag)
{
  char buf[MAXLEN_PATH_EXTRA];

  cone         = conenr; /* note: cone is here a variable of the class, NOT a local variable */
  reorder_flag = reordered_flag;

  /* determine local and global particle numbers for current lightcone */

  int n_type[NTYPES];

  /* determine global particle numbers for file header */
  for(int n = 0; n < NTYPES; n++)
    n_type[n] = 0;

  for(int n = 0; n < Lp->NumPart; n++)
    {
      if(reorder_flag || LightCone->lightcone_is_cone_member(n, cone))
        n_type[Lp->P[n].getType()]++;
    }

  sumup_large_ints(NTYPES, n_type, ntot_type_all, Communicator);

  if(!reorder_flag)
    {
      /* prepare healpix look-up table */
      Lp->Npix = nside2npix(LIGHTCONE_ORDER_NSIDE);

      subdivide_evenly(Lp->Npix, NTask, ThisTask, &Lp->FirstPix, &Lp->NpixLoc);

      Lp->HealPixTab_PartCount = (int *)Mem.mymalloc("HealPixTab_PartCount", Lp->NpixLoc * sizeof(int));

      int *tmp_PartCount = (int *)Mem.mymalloc_clear("tmp_PartCount", Lp->Npix * sizeof(int));
      for(int n = 0; n < Lp->NumPart; n++)
        {
          if(LightCone->lightcone_is_cone_member(n, cone))
            tmp_PartCount[Lp->P[n].ipnest]++;
        }

      MPI_Allreduce(MPI_IN_PLACE, tmp_PartCount, Lp->Npix, MPI_INT, MPI_SUM, Communicator);

      memcpy(Lp->HealPixTab_PartCount, tmp_PartCount + Lp->FirstPix, Lp->NpixLoc * sizeof(int));

      Mem.myfree(tmp_PartCount);
    }
  else
    Lp->Npix = Lp->NpixLoc = 0;

  mpi_printf("\nLIGHTCONE: writing cone=%d\n", conenr);

  char lname[MAXLEN_PATH];
  if(reordered_flag)
    snprintf(lname, MAXLEN_PATH, "lightcone_treeorder");
  else
    snprintf(lname, MAXLEN_PATH, "lightcone");

  if(ThisTask == 0)
    {
      char buf[MAXLEN_PATH_EXTRA];
      snprintf(buf, MAXLEN_PATH_EXTRA, "%s/%s_%02d", All.OutputDir, lname, cone);
      mkdir(buf, 02755);
    }
  MPI_Barrier(Communicator);

  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          char buf[MAXLEN_PATH_EXTRA];
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s/%s_%02d/conedir_%04d", All.OutputDir, lname, cone, num);
          mkdir(buf, 02755);
        }
      MPI_Barrier(Communicator);
    }

  if(All.NumFilesPerSnapshot > 1)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/%s_%02d/conedir_%04d/%s_%04d", All.OutputDir, lname, cone, num, "conesnap", num);
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/%s_%02d/%s_%04d", All.OutputDir, lname, cone, "conesnap", num);

  write_multiple_files(buf, All.NumFilesPerSnapshot);

  if(!reorder_flag)
    {
      Mem.myfree(Lp->HealPixTab_PartCount);
      Lp->HealPixTab_PartCount = NULL;
    }
}

void lightcone_particle_io::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type)
{
  /* determine global and local particle numbers */
  for(int n = 0; n < NTYPES + 2; n++)
    n_type[n] = 0;

  for(int n = 0; n < Lp->NumPart; n++)
    if(reorder_flag || LightCone->lightcone_is_cone_member(n, cone))
      if(Lp->P[n].getType() < NTYPES)
        n_type[Lp->P[n].getType()]++;

#ifdef MERGERTREE
  n_type[NTYPES + 0] = MergerTree->Ntrees;
#else
  n_type[NTYPES + 0] = 0;
#endif

  n_type[NTYPES + 1] = Lp->NpixLoc;

  /* determine particle numbers of each type in file */
  if(ThisTask == writeTask)
    {
      for(int n = 0; n < NTYPES + 2; n++)
        ntot_type[n] = n_type[n];

      for(int task = writeTask + 1; task <= lastTask; task++)
        {
          long long nn[NTYPES + 2];
          MPI_Recv(&nn[0], NTYPES + 2, MPI_LONG_LONG, task, TAG_LOCALN, Communicator, MPI_STATUS_IGNORE);
          for(int n = 0; n < NTYPES + 2; n++)
            ntot_type[n] += nn[n];
        }

      for(int task = writeTask + 1; task <= lastTask; task++)
        MPI_Send(&ntot_type[0], NTYPES + 2, MPI_LONG_LONG, task, TAG_N, Communicator);
    }
  else
    {
      MPI_Send(&n_type[0], NTYPES + 2, MPI_LONG_LONG, writeTask, TAG_LOCALN, Communicator);
      MPI_Recv(&ntot_type[0], NTYPES + 2, MPI_LONG_LONG, writeTask, TAG_N, Communicator, MPI_STATUS_IGNORE);
    }

  /* fill file header */

  for(int n = 0; n < NTYPES; n++)
    {
      header.npart[n]      = ntot_type[n];
      header.npartTotal[n] = ntot_type_all[n];
    }

  if(reorder_flag)
    {
      header.Ntrees = ntot_type[NTYPES];
#ifdef MERGERTREE
      header.TotNtrees = MergerTree->TotNtrees;
#else
      header.TotNtrees = 0;
#endif
      header.Npix    = 0;
      header.TotNpix = 0;
    }
  else
    {
      header.Ntrees    = 0;
      header.TotNtrees = 0;

      header.Npix    = ntot_type[NTYPES + 1];
      header.TotNpix = Lp->Npix;
    }

  header.num_files = All.NumFilesPerSnapshot;

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  int oindex = LightCone->Cones[cone].OriginIndex;
  for(int n = 0; n < 3; n++)
    header.Origin[n] = LightCone->ConeOrigins[oindex].PosOrigin[n];
#endif
}

void lightcone_particle_io::write_header_fields(hid_t handle)
{
  write_vector_attribute(handle, "NumPart_ThisFile", header.npart, H5T_NATIVE_UINT64, NTYPES);
  write_vector_attribute(handle, "NumPart_Total", header.npartTotal, H5T_NATIVE_UINT64, NTYPES);
  write_scalar_attribute(handle, "NumFiles", &header.num_files, H5T_NATIVE_INT);

  if(header.TotNtrees > 0)
    {
      write_scalar_attribute(handle, "Ntrees_ThisFile", &header.Ntrees, H5T_NATIVE_UINT64);
      write_scalar_attribute(handle, "Ntrees_Total", &header.TotNtrees, H5T_NATIVE_UINT64);
    }

  if(header.TotNpix > 0)
    {
      write_scalar_attribute(handle, "Npix_ThisFile", &header.Npix, H5T_NATIVE_UINT32);
      write_scalar_attribute(handle, "Npix_Total", &header.TotNpix, H5T_NATIVE_UINT32);
    }

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  write_vector_attribute(handle, "Origin", header.Origin, H5T_NATIVE_DOUBLE, 3);
#endif
}

void lightcone_particle_io::set_filenr_in_header(int numfiles) { header.num_files = numfiles; }

void lightcone_particle_io::get_datagroup_name(int type, char *buf)
{
  if(type < NTYPES)
    snprintf(buf, MAXLEN_PATH, "/PartType%d", type);
  else if(type == NTYPES)
    snprintf(buf, MAXLEN_PATH, "/TreeTable");
  else if(type == NTYPES + 1)
    snprintf(buf, MAXLEN_PATH, "/HealPixHashTable");
  else
    Terminate("wrong group");
}

void *lightcone_particle_io::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_LC:
        return (void *)(Lp->P + index);

      case A_PS:
        return (void *)(Lp->PS + index);

#ifdef MERGERTREE
      case A_TT:
        return (void *)(MergerTree->TreeTable + index);
#endif

      case A_MM:
        return (void *)(Lp->HealPixTab_PartCount + index);

      default:
        Terminate("we don't expect to get here");
    }

  return NULL;
}

void lightcone_particle_io::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *n_type,
                                             long long *ntot_type, int *nstart)
{
  n_type[NTYPES]        = 0;
  ntot_type[NTYPES]     = 0;
  n_type[NTYPES + 1]    = 0;
  ntot_type[NTYPES + 1] = 0;

  if(Lp->TotNumPart == 0)
    {
      if(header.num_files <= 1)
        for(int i = 0; i < NTYPES; i++)
          header.npartTotal[i] = header.npart[i];

      Lp->TotNumPart = 0;

      for(int i = 0; i < NTYPES; i++)
        Lp->TotNumPart += header.npartTotal[i];
    }

  if(ThisTask == readTask)
    {
      if(filenr == 0 && nstart == NULL)
        {
          for(int type = 0; type < NTYPES; type++)
            {
              mpi_printf("READ-LIGHTCONE: Type %d:         %8lld  (tot=%15lld)\n", type, (long long)header.npart[type],
                         (long long)header.npartTotal[type]);
            }
          mpi_printf("\n");
        }
    }

  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

  long long nall = 0;
  for(int type = 0; type < NTYPES; type++)
    {
      ntot_type[type] = header.npart[type];

      long long n_in_file = header.npart[type];
      int ntask           = lastTask - readTask + 1;
      int n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
        n_for_this_task++;

      n_type[type] = n_for_this_task;

      nall += n_for_this_task;
    }

  if(nstart)
    {
      memmove(&Lp->P[nall], &Lp->P[0], Lp->NumPart * sizeof(lightcone_particle_data));

      *nstart = 0;
    }
}

void lightcone_particle_io::read_header_fields(const char *fname)
{
  for(int i = 0; i < NTYPES; i++)
    {
      header.npart[i]      = 0;
      header.npartTotal[i] = 0;
    }

  header.Ntrees    = 0;
  header.TotNtrees = 0;

  int ntypes = NTYPES;

  hid_t hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t handle    = my_H5Gopen(hdf5_file, "/Header");

  /* check if the file in question actually has this number of types */
  hid_t hdf5_attribute = my_H5Aopen_name(handle, "NumPart_ThisFile");
  hid_t space          = H5Aget_space(hdf5_attribute);
  hsize_t dims, len;
  H5Sget_simple_extent_dims(space, &dims, &len);
  H5Sclose(space);
  if(len != (size_t)ntypes)
    Terminate("Length of NumPart_ThisFile attribute (%d) does not match NTYPES(ICS) (%d)", (int)len, ntypes);
  my_H5Aclose(hdf5_attribute, "NumPart_ThisFile");

  /* now read the header fields */
  read_vector_attribute(handle, "NumPart_ThisFile", header.npart, H5T_NATIVE_UINT64, ntypes);
  read_vector_attribute(handle, "NumPart_Total", header.npartTotal, H5T_NATIVE_UINT64, ntypes);
  read_scalar_attribute(handle, "NumFiles", &header.num_files, H5T_NATIVE_INT);

  my_H5Gclose(handle, "/Header");
  my_H5Fclose(hdf5_file, fname);
}

void lightcone_particle_io::read_increase_numbers(int type, int n_for_this_task) { Lp->NumPart += n_for_this_task; }

int lightcone_particle_io::get_filenr_from_header(void) { return header.num_files; }

int lightcone_particle_io::get_type_of_element(int index)
{
  if(index < 0 || index >= Lp->NumPart)
    Terminate("index = %d  Lp->NumPart=%d", index, Lp->NumPart);

  if(reorder_flag || LightCone->lightcone_is_cone_member(index, cone))
    return Lp->P[index].getType();
  else
    return -1; /* this will skip this particle */
}

void lightcone_particle_io::set_type_of_element(int index, int type)
{
  if(type < NTYPES)
    Lp->P[index].setType(type);
}

#endif
