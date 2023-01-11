/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  io_progenitors.cc
 *
 *  \brief routines for I/O of the progenitor links between two subhalo catalogues
 */

#include "gadgetconfig.h"

#ifdef MERGERTREE

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
#include "../fof/fof.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mergertree/io_progenitors.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

progenitors_io::progenitors_io(mergertree *MergerTree_ptr, MPI_Comm comm, int format) : IO_Def(comm, format)
{
  MergerTree = MergerTree_ptr;

  this->N_IO_Fields  = 0;
  this->N_DataGroups = 1;
  this->header_size  = sizeof(header);
  this->header_buf   = &header;
  this->type_of_file = FILE_IS_PROGCAT;
  snprintf(this->info, MAXLEN_PATH, "MERGERTREE: writing progenitor information");

  init_field("PSNR", "ProgSubhaloNr", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_PROG, NULL, io_func_progsubhalonr, CURRSUBS, 0, 0,
             0, 0, 0, 0, 0, true);

  init_field("FPNR", "FirstProgSubhaloNr", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_PROG, NULL, io_func_firstprogsubhalonr,
             CURRSUBS, 0, 0, 0, 0, 0, 0, 0, true);

  init_field("NDNR", "NextDescSubhaloNr", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_PROG, NULL, io_func_nextdescsubhalonr, CURRSUBS,
             0, 0, 0, 0, 0, 0, 0, true);

  init_field("SHNR", "SubhaloNr", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_PROG, &MergerTree->Progenitors[0].SubhaloNr, NULL,
             CURRSUBS, 0, 0, 0, 0, 0, 0, 0, true);  // this field can in principle be deleted

  init_field("PFLO", "ProgFileOffset", MEM_MY_FILEOFFSET, FILE_NONE, SKIP_ON_READ, 1, A_PROG, &MergerTree->Progenitors[0].FileOffset,
             NULL, CURRSUBS, 0, 0, 0, 0, 0, 0, 0, true);
}

void progenitors_io::mergertree_save_progenitors(int num)
{
  char buf[MAXLEN_PATH_EXTRA];

  /* write FirstProgenitor */
  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          snprintf(buf, MAXLEN_PATH, "%s/groups_%03d", All.OutputDir, num);
          mkdir(buf, 02755);
        }
      MPI_Barrier(Communicator);
    }

  if(All.NumFilesPerSnapshot > 1)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/groups_%03d/%s_%03d", All.OutputDir, num, "subhalo_prog", num);
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, "subhalo_prog", num);

  write_multiple_files(buf, All.NumFilesPerSnapshot);
}

void progenitors_io::mergertree_read_progenitors(int num)
{
  char fname[MAXLEN_PATH_EXTRA], fname_multiple[MAXLEN_PATH_EXTRA];

  snprintf(fname_multiple, MAXLEN_PATH_EXTRA, "%s/groups_%03d/%s_%03d", All.OutputDir, num, "subhalo_prog", num);
  snprintf(fname, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, "subhalo_prog", num);

  TotNsubhalos = 0;

  int num_files = find_files(fname, fname_multiple);

  if(num_files > 1)
    strcpy(fname, fname_multiple);

  /* we repeat reading the headers of the files two times. In the first iteration, only the
   * particle numbers ending up on each processor are assembled, followed by memory allocation.
   * In the second iteration, the data is actually read in.
   */
  for(int rep = 0; rep < 2; rep++)
    {
      Nsubhalos = 0;

      read_files_driver(fname, rep, num_files);

      /* now do the memory allocation */
      if(rep == 0)
        {
          MergerTree->Progenitors = (mergertree::prog_list *)Mem.mymalloc_movable(&MergerTree->Progenitors, "Progenitors",
                                                                                  Nsubhalos * sizeof(mergertree::prog_list));
        }
    }

  MPI_Barrier(Communicator);
}

void progenitors_io::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type)
{
  /* determine group/id numbers of each type in file */

  n_type[0] = MergerTree->CurrNsubhalos;

  if(ThisTask == writeTask)
    {
      for(int n = 0; n < 1; n++)
        ntot_type[n] = n_type[n];

      for(int task = writeTask + 1; task <= lastTask; task++)
        {
          long long nn[3];
          MPI_Recv(&nn[0], 1, MPI_LONG_LONG, task, TAG_LOCALN, Communicator, MPI_STATUS_IGNORE);
          for(int n = 0; n < 1; n++)
            ntot_type[n] += nn[n];
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

  header.Nsubhalos    = ntot_type[0];
  header.TotNsubhalos = MergerTree->CurrTotNsubhalos;
  header.num_files    = All.NumFilesPerSnapshot;
}

void progenitors_io::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *n_type,
                                      long long *ntot_type, int *nstart)
{
  if(ThisTask == readTask)
    {
      if(filenr == 0 && nstart == NULL)
        {
          mpi_printf("\nREAD-PROGENITORS: filenr=%d, '%s' contains  (subhalos):  %8d  out of a total of %lld\n", filenr, fname,
                     header.Nsubhalos, header.TotNsubhalos);
        }
    }

  if(TotNsubhalos == 0)
    TotNsubhalos = header.TotNsubhalos;

  for(int k = 0; k < 1; k++)
    n_type[k] = ntot_type[k] = 0;

  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

  ntot_type[0] = header.Nsubhalos;

  long long n_in_file = header.Nsubhalos;
  int ntask           = lastTask - readTask + 1;
  int n_for_this_task = n_in_file / ntask;
  if((ThisTask - readTask) < (n_in_file % ntask))
    n_for_this_task++;

  n_type[0] = n_for_this_task;

  if(nstart)
    {
      memmove(&MergerTree->Progenitors[n_for_this_task], &MergerTree->Progenitors[0], Nsubhalos * sizeof(mergertree::prog_list));
      *nstart = 0;
    }
}

void progenitors_io::write_header_fields(hid_t handle)
{
  write_scalar_attribute(handle, "Nsubhalos_ThisFile", &header.Nsubhalos, H5T_NATIVE_UINT64);

  write_scalar_attribute(handle, "Nsubhalos_Total", &header.TotNsubhalos, H5T_NATIVE_UINT64);

  write_scalar_attribute(handle, "NumFiles", &header.num_files, H5T_NATIVE_INT);
}

/*! \brief This function reads the snapshot header in case of hdf5 files (i.e. format 3)
 *
 * \param fname file name of the snapshot as given in the parameter file
 */
void progenitors_io::read_header_fields(const char *fname)
{
  memset(&header, 0, sizeof(io_header));

  hid_t hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t handle    = my_H5Gopen(hdf5_file, "/Header");

  /* now read the header fields */
  read_scalar_attribute(handle, "Nsubhalos_ThisFile", "Nsubgroups_ThisFile", &header.Nsubhalos, H5T_NATIVE_UINT64);

  read_scalar_attribute(handle, "Nsubhalos_Total", "Nsubgroups_Total", &header.TotNsubhalos, H5T_NATIVE_UINT64);

  read_scalar_attribute(handle, "NumFiles", &header.num_files, H5T_NATIVE_INT);

  my_H5Gclose(handle, "/Header");
  my_H5Fclose(hdf5_file, fname);
}

int progenitors_io::get_filenr_from_header(void) { return header.num_files; }

void progenitors_io::set_filenr_in_header(int numfiles) { header.num_files = numfiles; }

void progenitors_io::read_increase_numbers(int type, int n_for_this_task)
{
  switch(type)
    {
      case 0:
        Nsubhalos += n_for_this_task;
        break;
      default:
        Terminate("wrong group");
        break;
    }
}

void progenitors_io::get_datagroup_name(int type, char *buf)
{
  switch(type)
    {
      case 0:
        snprintf(buf, MAXLEN_PATH, "/Subhalo");
        break;
      default:
        Terminate("wrong group");
        break;
    }
}

int progenitors_io::get_type_of_element(int index)
{
  /* empty */
  return 0;
}

void progenitors_io::set_type_of_element(int index, int type)
{ /* empty */
}

void *progenitors_io::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_PROG:
        return (void *)(MergerTree->Progenitors + index);
      default:
        Terminate("strange, we don't expect to get here");
    }
}

#endif
