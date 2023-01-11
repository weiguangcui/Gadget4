/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file fof_io.cc
 *
 *  \brief routines for implementing the I/O routines concerned with the group catalogues
 */

#include "gadgetconfig.h"

#ifdef FOF

#include <hdf5.h>
#include <mpi.h>
#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../fof/fof.h"
#include "../fof/fof_io.h"
#include "../gitversion/version.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

/*! \file fof_io.c
 *  \brief parallel I/O routines for the FOF and SUBFIND group finders
 */

template <typename partset>
fof_io<partset>::fof_io(fof<partset> *FoF_ptr, MPI_Comm comm, int format) : IO_Def(comm, format)
{
  FoF = FoF_ptr;

  this->N_IO_Fields  = 0;
  this->N_DataGroups = 3;
  this->header_size  = sizeof(catalogue_header);
  this->header_buf   = &catalogue_header;
  this->type_of_file = FILE_IS_GROUPCAT;
  snprintf(this->info, MAXLEN_PATH, "FOF/SUBFIND: writing group catalogue");

  init_field("FLEN", "GroupLen", mem_len_type, file_len_type, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].Len, NULL, GROUPS, 0, 0, 0, 0, 0,
             0, 0, true);

  init_field("FMAS", "GroupMass", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].Mass, NULL, GROUPS, 0, 0, 0,
             0, 0, 0, 0);

  init_field("FPOS", "GroupPos", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_G, &FoF->Group[0].Pos[0], NULL, GROUPS, 0, 0,
             0, 0, 0, 0, 0);

  init_field("FVEL", "GroupVel", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_G, &FoF->Group[0].Vel[0], NULL, GROUPS, 0, 0, 0,
             0, 0, 0, 0);

  init_field("FLTY", "GroupLenType", mem_len_type, file_len_type, READ_IF_PRESENT, NTYPES, A_G, &FoF->Group[0].LenType[0], NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0, true);

  init_field("FLOT", "GroupOffsetType", MEM_INT64, FILE_INT64, READ_IF_PRESENT, NTYPES, A_G, &FoF->Group[0].OffsetType[0], NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0, true);

  init_field("FMTY", "GroupMassType", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, READ_IF_PRESENT, NTYPES, A_G, &FoF->Group[0].MassType[0], NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("FMAA", "GroupAscale", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].Ascale, NULL, GROUPS, 0,
             0, 0, 0, 0, 0, 0);

#if defined(SUBFIND_ORPHAN_TREATMENT)
  init_field("FPEN", "GroupLenPrevMostBnd", MEM_INT, FILE_INT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].LenPrevMostBnd, NULL, GROUPS, 0,
             0, 0, 0, 0, 0, 0, true);
#endif

#ifdef STARFORMATION
  init_field("FSFR", "GroupSFR", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_G, &FoF->Group[0].Sfr, NULL, GROUPS, 0, 0, 0, 0, 0,
             0, 0);
#endif

#ifdef MERGERTREE
  init_field("GFLO", "GroupFileOffset", MEM_MY_FILEOFFSET, FILE_NONE, SKIP_ON_READ, 1, A_G, &FoF->Group[0].FileOffset, NULL, GROUPS, 0,
             0, 0, 0, 0, 0, 0, true);
#endif

#ifdef SUBFIND
  // ------------------ additional Group fields -------------------------------------

  init_field("FMM2", "Group_M_Mean200", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].M_Mean200, NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);
  init_field("FRM2", "Group_R_Mean200", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].R_Mean200, NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("FMC2", "Group_M_Crit200", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].M_Crit200, NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);
  init_field("FRC2", "Group_R_Crit200", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].R_Crit200, NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("FMC5", "Group_M_Crit500", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].M_Crit500, NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);
  init_field("FRC5", "Group_R_Crit500", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].R_Crit500, NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("FMT2", "Group_M_TopHat200", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].M_TopHat200, NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);
  init_field("FMR2", "Group_R_TopHat200", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].R_TopHat200, NULL,
             GROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("FNSH", "GroupNsubs", MEM_INT, FILE_INT, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].Nsubs, NULL, GROUPS, 0, 0, 0, 0, 0, 0, 0,
             true);

  init_field("FFSH", "GroupFirstSub", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_G, &FoF->Group[0].FirstSub, NULL, GROUPS, 0, 0, 0,
             0, 0, 0, 0, true);

  // ------------------- genuine Subhalo fields ------------------------------------

  init_field("SGNR", "SubhaloGroupNr", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].GroupNr, NULL, SUBGROUPS, 0, 0,
             0, 0, 0, 0, 0, true);

  init_field("SGRG", "SubhaloRankInGr", MEM_INT, FILE_INT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].SubRankInGr, NULL, SUBGROUPS, 0,
             0, 0, 0, 0, 0, 0, true);

  init_field("SLEN", "SubhaloLen", mem_len_type, file_len_type, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].Len, NULL, SUBGROUPS, 0, 0,
             0, 0, 0, 0, 0, true);

  init_field("SMAS", "SubhaloMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].Mass, NULL, SUBGROUPS, 0,
             0, 0, 0, 0, 0, 0);

  init_field("SPOS", "SubhaloPos", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_S, &FoF->Subhalo[0].Pos[0], NULL, SUBGROUPS,
             0, 0, 0, 0, 0, 0, 0);

  init_field("SVEL", "SubhaloVel", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_S, &FoF->Subhalo[0].Vel[0], NULL, SUBGROUPS,
             0, 0, 0, 0, 0, 0, 0);

  init_field("SLTY", "SubhaloLenType", mem_len_type, file_len_type, READ_IF_PRESENT, NTYPES, A_S, &FoF->Subhalo[0].LenType[0], NULL,
             SUBGROUPS, 0, 0, 0, 0, 0, 0, 0, true);

  init_field("SLOT", "SubhaloOffsetType", MEM_INT64, FILE_INT64, READ_IF_PRESENT, NTYPES, A_S, &FoF->Subhalo[0].OffsetType[0], NULL,
             SUBGROUPS, 0, 0, 0, 0, 0, 0, 0, true);

  init_field("SMTY", "SubhaloMassType", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, NTYPES, A_S, &FoF->Subhalo[0].MassType[0],
             NULL, SUBGROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("SCMP", "SubhaloCM", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_S, &FoF->Subhalo[0].CM[0], NULL, SUBGROUPS, 0,
             0, 0, 0, 0, 0, 0);

  init_field("SSPI", "SubhaloSpin", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_S, &FoF->Subhalo[0].Spin[0], NULL, SUBGROUPS,
             0, 0, 0, 0, 0, 0, 0);

  init_field("SSPI", "SubhaloVelDisp", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].SubVelDisp, NULL,
             SUBGROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("SVMX", "SubhaloVmax", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].SubVmax, NULL, SUBGROUPS,
             0, 0, 0, 0, 0, 0, 0);

  init_field("SVRX", "SubhaloVmaxRad", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].SubVmaxRad, NULL,
             SUBGROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("SHMR", "SubhaloHalfmassRad", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].SubHalfMassRad,
             NULL, SUBGROUPS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SHMT", "SubhaloHalfmassRadType", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, NTYPES, A_S,
             &FoF->Subhalo[0].SubHalfMassRadType[0], NULL, SUBGROUPS, 0, 0, 0, 0, 0, 0, 0);

  init_field("SIDM", "SubhaloIDMostbound", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].SubMostBoundID,
             NULL, SUBGROUPS, 0, 0, 0, 0, 0, 0, 0, true);

  init_field("SPRT", "SubhaloParentRank", MEM_INT, FILE_INT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].SubParentRank, NULL, SUBGROUPS,
             0, 0, 0, 0, 0, 0, 0, true);

#ifdef MERGERTREE
  init_field("SFLO", "SubhaloFileOffset", MEM_MY_FILEOFFSET, FILE_NONE, SKIP_ON_READ, 1, A_S, &FoF->Subhalo[0].FileOffset, NULL,
             SUBGROUPS, 0, 0, 0, 0, 0, 0, 0, true);
#endif

#if defined(SUBFIND_ORPHAN_TREATMENT)
  init_field("LPMO", "SubhaloLenPrevMostBnd", MEM_INT, FILE_INT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].SubhaloLenPrevMostBnd, NULL,
             SUBGROUPS, 0, 0, 0, 0, 0, 0, 0, true);
#endif

#ifdef STARFORMATION
  init_field("SSFR", "SubhaloSFR", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].Sfr, NULL, SUBGROUPS, 0,
             0, 0, 0, 0, 0, 0);
  init_field("SSFG", "SubhaloGasMassSFR", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_S, &FoF->Subhalo[0].GasMassSfr, NULL,
             SUBGROUPS, 0, 0, 0, 0, 0, 0, 0);
#endif

#endif
}

template <typename partset>
int fof_io<partset>::get_type_of_element(int index)
{
  return 0; /* not needed here */
}

template <typename partset>
void fof_io<partset>::set_type_of_element(int index, int type)
{
  /* empty */
}

/* prepare list of ids with assigned group numbers */

template <typename partset>
void fof_io<partset>::fof_subfind_save_groups(int num, const char *basename, const char *grpcat_dirbasename)
{
#ifdef DO_NOT_PRODUCE_BIG_OUTPUT
  mpi_printf("FOF/SUBFIND: We skip saving group catalogues.\n");
  return;
#endif

  char buf[MAXLEN_PATH_EXTRA];

  double t0 = Logs.second();
  reset_io_byte_count();

  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s/%s_%03d", All.OutputDir, grpcat_dirbasename, num);
          mkdir(buf, 02755);
        }
      MPI_Barrier(Communicator);
    }

  if(All.NumFilesPerSnapshot > 1)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/%s_%03d/%s_%03d", All.OutputDir, grpcat_dirbasename, num, basename, num);
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, basename, num);

  write_multiple_files(buf, All.NumFilesPerSnapshot);

  long long byte_count = get_io_byte_count(), byte_count_all;
  sumup_longs(1, &byte_count, &byte_count_all, Communicator);

  double t1 = Logs.second();

  mpi_printf("FOF/SUBFIND: Group catalogues saved. took = %g sec, total size %g MB, corresponds to effective I/O rate of %g MB/sec\n",
             Logs.timediff(t0, t1), byte_count_all / (1024.0 * 1024.0), byte_count_all / (1024.0 * 1024.0) / Logs.timediff(t0, t1));
}

template <typename partset>
void fof_io<partset>::fof_subfind_load_groups(int num)
{
  FoF->TotNgroups = 0;

  char fname[MAXLEN_PATH_EXTRA], fname_multiple[MAXLEN_PATH_EXTRA];

  snprintf(fname_multiple, MAXLEN_PATH_EXTRA, "%s/groups_%03d/%s_%03d", All.OutputDir, num, "fof_subhalo_tab", num);
  snprintf(fname, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, "fof_subhalo_tab", num);

  int num_files = find_files(fname, fname_multiple);

  if(num_files > 1)
    strcpy(fname, fname_multiple);

  /* we repeat reading the headers of the files two times. In the first iteration, only the
   * particle numbers ending up on each processor are assembled, followed by memory allocation.
   * In the second iteration, the data is actually read in.
   */

  for(int rep = 0; rep < 2; rep++)
    {
      FoF->Ngroups   = 0;
      FoF->Nsubhalos = 0;

      read_files_driver(fname, rep, num_files);

      /* now do the memory allocation */
      if(rep == 0)
        {
          FoF->Group = (typename fof<partset>::group_properties *)Mem.mymalloc_movable(
              &FoF->Group, "Group", FoF->Ngroups * sizeof(typename fof<partset>::group_properties));
#ifdef SUBFIND
          FoF->Subhalo = (typename fof<partset>::subhalo_properties *)Mem.mymalloc_movable(
              &FoF->Subhalo, "Subhalo", FoF->Nsubhalos * sizeof(typename fof<partset>::subhalo_properties));
#endif
        }
    }

  MPI_Barrier(Communicator);

  mpi_printf("\nFOF/SUBFIND: reading done.\n");
  mpi_printf("FOF/SUBFIND: Total number of groups read:    %lld\n", (long long int)FoF->TotNgroups);
  mpi_printf("FOF/SUBFIND: Total number of subhalos read: %lld\n\n", (long long int)FoF->TotNsubhalos);

  MPI_Allreduce(MPI_IN_PLACE, &LegacyFormat, 1, MPI_INT, MPI_MAX, Communicator);
}

template <typename partset>
void fof_io<partset>::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type)
{
  /* determine group/id numbers of each type in file */

  n_type[0] = FoF->Ngroups;
  n_type[1] = FoF->Nsubhalos;
  n_type[2] = FoF->Nids;

  if(ThisTask == writeTask)
    {
      for(int n = 0; n < 3; n++)
        ntot_type[n] = n_type[n];

      for(int task = writeTask + 1; task <= lastTask; task++)
        {
          long long nn[3];
          MPI_Recv(&nn[0], 3, MPI_LONG_LONG, task, TAG_LOCALN, Communicator, MPI_STATUS_IGNORE);
          for(int n = 0; n < 3; n++)
            ntot_type[n] += nn[n];
        }

      for(int task = writeTask + 1; task <= lastTask; task++)
        MPI_Send(&ntot_type[0], 3, MPI_LONG_LONG, task, TAG_N, Communicator);
    }
  else
    {
      MPI_Send(&n_type[0], 3, MPI_LONG_LONG, writeTask, TAG_LOCALN, Communicator);
      MPI_Recv(&ntot_type[0], 3, MPI_LONG_LONG, writeTask, TAG_N, Communicator, MPI_STATUS_IGNORE);
    }

  /* fill file header */

  catalogue_header.Ngroups   = ntot_type[0];
  catalogue_header.Nsubhalos = ntot_type[1];
  catalogue_header.Nids      = ntot_type[2];

  catalogue_header.TotNgroups   = FoF->TotNgroups;
  catalogue_header.TotNsubhalos = FoF->TotNsubhalos;
  catalogue_header.TotNids      = FoF->TotNids;

  catalogue_header.num_files = All.NumFilesPerSnapshot;

  catalogue_header.time = All.Time;
  if(All.ComovingIntegrationOn)
    catalogue_header.redshift = 1.0 / All.Time - 1;
  else
    catalogue_header.redshift = 0;
  catalogue_header.BoxSize = All.BoxSize;
}

template <typename partset>
void fof_io<partset>::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *n_type,
                                       long long *ntot_type, int *nstart)
{
  if(ThisTask == readTask)
    {
      if(filenr == 0 && nstart == NULL)
        {
          mpi_printf(
              "\nREAD-FOF: filenr=%d, '%s' contains:\n"
              "READ-FOF: Type 0 (groups):    %8lld\n"
              "READ-FOF: Type 1 (subhalos):  %8lld\n"
              "READ-FOF: Type 2 (ids):       %8lld\n",
              filenr, fname, catalogue_header.Ngroups, catalogue_header.Nsubhalos, catalogue_header.Nids);
        }
    }

  if(FoF->TotNgroups == 0)
    {
      FoF->TotNgroups   = catalogue_header.TotNgroups;
      FoF->TotNsubhalos = catalogue_header.TotNsubhalos;
      FoF->TotNids      = catalogue_header.TotNids;
    }

  FoF->Redshift = catalogue_header.redshift;
  FoF->Time     = catalogue_header.time;

  for(int k = 0; k < 3; k++)
    n_type[k] = ntot_type[k] = 0;

  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

  {
    ntot_type[0] = catalogue_header.Ngroups;

    long long n_in_file = catalogue_header.Ngroups;
    int ntask           = lastTask - readTask + 1;
    int n_for_this_task = n_in_file / ntask;
    if((ThisTask - readTask) < (n_in_file % ntask))
      n_for_this_task++;

    n_type[0] = n_for_this_task;

    if(nstart)
      memmove(&FoF->Group[n_for_this_task], &FoF->Group[0], FoF->Ngroups * sizeof(typename fof<partset>::group_properties));
  }

#ifdef SUBFIND
  {
    ntot_type[1] = catalogue_header.Nsubhalos;

    long long n_in_file = catalogue_header.Nsubhalos;
    int ntask           = lastTask - readTask + 1;
    int n_for_this_task = n_in_file / ntask;
    if((ThisTask - readTask) < (n_in_file % ntask))
      n_for_this_task++;

    n_type[1] = n_for_this_task;

    if(nstart)
      memmove(&FoF->Subhalo[n_for_this_task], &FoF->Subhalo[0], FoF->Nsubhalos * sizeof(typename fof<partset>::subhalo_properties));
  }
#endif

  if(nstart)
    *nstart = 0;
}

template <typename partset>
void fof_io<partset>::write_header_fields(hid_t handle)
{
  write_scalar_attribute(handle, "Ngroups_ThisFile", &catalogue_header.Ngroups, H5T_NATIVE_UINT64);
  write_scalar_attribute(handle, "Nsubhalos_ThisFile", &catalogue_header.Nsubhalos, H5T_NATIVE_UINT64);
  write_scalar_attribute(handle, "Nids_ThisFile", &catalogue_header.Nids, H5T_NATIVE_UINT64);

  write_scalar_attribute(handle, "Ngroups_Total", &catalogue_header.TotNgroups, H5T_NATIVE_UINT64);
  write_scalar_attribute(handle, "Nsubhalos_Total", &catalogue_header.TotNsubhalos, H5T_NATIVE_UINT64);
  write_scalar_attribute(handle, "Nids_Total", &catalogue_header.TotNids, H5T_NATIVE_UINT64);

  write_scalar_attribute(handle, "NumFiles", &catalogue_header.num_files, H5T_NATIVE_INT);

  write_scalar_attribute(handle, "Time", &catalogue_header.time, H5T_NATIVE_DOUBLE);
  write_scalar_attribute(handle, "Redshift", &catalogue_header.redshift, H5T_NATIVE_DOUBLE);
  write_scalar_attribute(handle, "BoxSize", &catalogue_header.BoxSize, H5T_NATIVE_DOUBLE);

  write_string_attribute(handle, "Git_commit", GIT_COMMIT);
  write_string_attribute(handle, "Git_date", GIT_DATE);
}

/*! \brief This function reads the snapshot header in case of hdf5 files (i.e. format 3)
 *
 * \param fname file name of the snapshot as given in the parameter file
 */
template <typename partset>
void fof_io<partset>::read_header_fields(const char *fname)
{
  memset(&catalogue_header, 0, sizeof(fof_subfind_header));

  hid_t hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t handle    = my_H5Gopen(hdf5_file, "/Header");

  /* now read the header fields */
  read_scalar_attribute(handle, "Ngroups_ThisFile", &catalogue_header.Ngroups, H5T_NATIVE_UINT64);

  if(read_scalar_attribute(handle, "Nsubhalos_ThisFile", "Nsubgroups_ThisFile", &catalogue_header.Nsubhalos, H5T_NATIVE_UINT64))
    LegacyFormat = 1;

  read_scalar_attribute(handle, "Nids_ThisFile", &catalogue_header.Nids, H5T_NATIVE_UINT64);

  read_scalar_attribute(handle, "Ngroups_Total", &catalogue_header.TotNgroups, H5T_NATIVE_UINT64);
  read_scalar_attribute(handle, "Nsubhalos_Total", "Nsubgroups_Total", &catalogue_header.TotNsubhalos, H5T_NATIVE_UINT64);
  read_scalar_attribute(handle, "Nids_Total", &catalogue_header.TotNids, H5T_NATIVE_UINT64);

  read_scalar_attribute(handle, "NumFiles", &catalogue_header.num_files, H5T_NATIVE_INT);

  read_scalar_attribute(handle, "Time", &catalogue_header.time, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "Redshift", &catalogue_header.redshift, H5T_NATIVE_DOUBLE);

  read_scalar_attribute(handle, "BoxSize", &catalogue_header.BoxSize, H5T_NATIVE_DOUBLE);

  my_H5Gclose(handle, "/Header");
  my_H5Fclose(hdf5_file, fname);
}

template <typename partset>
int fof_io<partset>::get_filenr_from_header(void)
{
  return catalogue_header.num_files;
}

template <typename partset>
void fof_io<partset>::set_filenr_in_header(int numfiles)
{
  catalogue_header.num_files = numfiles;
}

template <typename partset>
void fof_io<partset>::read_increase_numbers(int type, int n_for_this_task)
{
  switch(type)
    {
      case 0:
        FoF->Ngroups += n_for_this_task;
        break;
      case 1:
        FoF->Nsubhalos += n_for_this_task;
        break;
      case 2:
        FoF->Nids += n_for_this_task;
        break;
      default:
        Terminate("wrong group");
        break;
    }
}

template <typename partset>
void fof_io<partset>::get_datagroup_name(int type, char *buf)
{
  switch(type)
    {
      case 0:
        snprintf(buf, MAXLEN_PATH, "/Group");
        break;
      case 1:
        snprintf(buf, MAXLEN_PATH, "/Subhalo");
        break;
      case 2:
        snprintf(buf, MAXLEN_PATH, "/IDs");
        break;
      default:
        Terminate("wrong group: type=%d", type);
        break;
    }
}

template <typename partset>
void *fof_io<partset>::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_G:
        return (void *)(FoF->Group + index);

#ifdef SUBFIND
      case A_S:
        return (void *)(FoF->Subhalo + index);
#endif

      default:
        Terminate("we don't expect to get here");
    }

  return NULL;
}

template <typename partset>
void fof_io<partset>::fof_subfind_prepare_ID_list(void)
{
  double t0 = Logs.second();

  ID_list = (id_list *)Mem.mymalloc("ID_list", sizeof(id_list) * FoF->Nids);

  long long nids = 0;
  for(int i = 0; i < FoF->Tp->NumPart; i++)
    {
      if(FoF->Tp->PS[i].GroupNr.get() < HALONR_MAX)
        {
          if(nids >= FoF->Nids)
            Terminate("nids >= Nids");

          ID_list[nids].GroupNr = FoF->Tp->PS[i].GroupNr;
          ID_list[nids].Type    = FoF->Tp->P[i].getType();
          ID_list[nids].ID      = FoF->Tp->P[i].ID.get();
#ifdef SUBFIND
          ID_list[nids].SubRankInGr = FoF->Tp->PS[i].SubRankInGr;
          ID_list[nids].BindingEgy  = FoF->Tp->PS[i].v.DM_BindingEnergy;
#endif
          nids++;
        }
    }

  long long totNids;
  sumup_longs(1, &nids, &totNids, Communicator);
  if(totNids != FoF->TotNids)
    Terminate("Task=%d Nids=%lld totNids=%lld TotNids=%lld\n", ThisTask, FoF->Nids, totNids, FoF->TotNids);

  /* sort the particle IDs according to group-number, and optionally subhalo number and binding energy  */
  mycxxsort_parallel(ID_list, ID_list + FoF->Nids, fof_subfind_compare_ID_list, Communicator);

  double t1 = Logs.second();
  mpi_printf("FOF/SUBFIND: Particle/cell IDs in groups globally sorted. took = %g sec\n", Logs.timediff(t0, t1));
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof_io<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof_io<lcparticles>;
#endif

#endif
