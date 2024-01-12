/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  io_readtrees_mbound.cc
 *
 *  \brief routines for I/O of most-bound particles belonging to merger trees
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
#include "../mergertree/io_readtrees_mbound.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

/*
 struct treehalo_ids_type
  {
    MyIDType SubMostBoundID;
    long long TreeID;
  };
  treehalo_ids_type *HaloIDdata;
*/

readtrees_mbound_io::readtrees_mbound_io(mergertree *MergerTree_ptr, MPI_Comm comm, int format) : IO_Def(comm, format)
{
  MergerTree = MergerTree_ptr;

  this->N_IO_Fields  = 0;
  this->N_DataGroups = 2;
  this->header_size  = sizeof(header);
  this->header_buf   = &header;
  this->type_of_file = FILE_IS_TREECAT;
  snprintf(this->info, MAXLEN_PATH, "MERGERTREE: reading/writing mergertrees");

  init_field("MTRL", "Length", MEM_INT, FILE_INT, READ_IF_PRESENT, 1, A_TT, &MergerTree->TreeTable[0].HaloCount, NULL, TREELENGTH, 0,
             0, 0, 0, 0, 0, 0);
  init_field("MTRS", "StartOffset", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_TT, &MergerTree->TreeTable[0].FirstHalo, NULL,
             TREELENGTH, 0, 0, 0, 0, 0, 0, 0);
  init_field("MTRI", "TreeID", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_TT, &MergerTree->TreeTable[0].TreeID, NULL, TREELENGTH, 0,
             0, 0, 0, 0, 0, 0);

  /* just read most-bound ID */

  init_field("TRID", "TreeID", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_TID, &MergerTree->HaloIDdata[0].TreeID, NULL, TREEHALOS, 0,
             0, 0, 0, 0, 0, 0);

  init_field("SIDM", "SubhaloIDMostbound", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, READ_IF_PRESENT, 1, A_TID,
             &MergerTree->HaloIDdata[0].SubMostBoundID, NULL, TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
}

void readtrees_mbound_io::read_trees_mostbound(void)
{
  double t0 = Logs.second();

  MergerTree->TotNtrees = 0;
  MergerTree->TotNhalos = 0;

  char fname[MAXLEN_PATH_EXTRA];

  if(All.NumFilesPerSnapshot > 1)
    snprintf(fname, MAXLEN_PATH_EXTRA, "%s/treedata/%s", All.OutputDir, "trees");
  else
    snprintf(fname, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "trees");

  int num_files = find_files(fname, fname);

  reset_io_byte_count();

  for(int rep = 0; rep < 2; rep++)
    {
      MergerTree->Ntrees = 0;
      MergerTree->Nhalos = 0;

      read_files_driver(fname, rep, num_files);

      /* now do the memory allocation */
      if(rep == 0)
        {
          MergerTree->TreeTable  = (halotrees_table *)Mem.mymalloc_movable(&MergerTree->TreeTable, "TreeTable",
                                                                           (MergerTree->Ntrees + 1) * sizeof(halotrees_table));
          MergerTree->HaloIDdata = (treehalo_ids_type *)Mem.mymalloc_movable(&MergerTree->HaloIDdata, "HaloIDdata",
                                                                             (MergerTree->Nhalos + 1) * sizeof(treehalo_ids_type));
        }
    }

  MPI_Barrier(Communicator);

  long long byte_count = get_io_byte_count(), byte_count_all;
  sumup_longs(1, &byte_count, &byte_count_all, Communicator);

  double t1 = Logs.second();

  mpi_printf("MERGERTREE-READ: reading done. Took %g sec, total size %g MB, corresponds to effective I/O rate of %g MB/sec\n",
             Logs.timediff(t0, t1), byte_count_all / (1024.0 * 1024.0), byte_count_all / (1024.0 * 1024.0) / Logs.timediff(t0, t1));

  mpi_printf("\nMERGERTREE-READ: Total number of trees=%ldd, total number of halos=%lld\n\n", MergerTree->TotNtrees,
             MergerTree->TotNhalos);

  MPI_Barrier(Communicator);

  for(int i = 0; i < MergerTree->Ntrees; i++)
    {
      if(MergerTree->HaloIDdata[i].TreeID > MergerTree->TotNtrees)
        Terminate("i=%d  MergerTree->Ntrees=%d  MergerTree->HaloIDdata[i].TreeID=%lld   MergerTree->TotNtrees=%lld ", i,
                  MergerTree->Ntrees, MergerTree->HaloIDdata[i].TreeID, MergerTree->TotNtrees);
    }
}

void readtrees_mbound_io::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type) {}

void readtrees_mbound_io::write_header_fields(hid_t handle) {}

int readtrees_mbound_io::get_filenr_from_header(void) { return header.num_files; }

void readtrees_mbound_io::set_filenr_in_header(int numfiles) { header.num_files = numfiles; }

void readtrees_mbound_io::read_increase_numbers(int type, int n_for_this_task)
{
  switch(type)
    {
      case 0:
        MergerTree->Ntrees += n_for_this_task;
        break;
      case 1:
        MergerTree->Nhalos += n_for_this_task;
        break;

      default:
        Terminate("wrong group");
        break;
    }
}

void readtrees_mbound_io::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *n_type,
                                           long long *ntot_type, int *nstart)
{
  if(ThisTask == readTask)
    {
      if(filenr == 0 && nstart == NULL)
        {
          mpi_printf(
              "\nREAD-TREES: filenr=%d, '%s' contains %lld trees out of a total of %lld, and %lld halos out of a total of %lld\n",
              filenr, fname, header.Ntrees, header.TotNtrees, header.Nhalos, header.TotNhalos);
        }
    }

  if(MergerTree->TotNtrees == 0)
    {
      MergerTree->TotNtrees = header.TotNtrees;
      MergerTree->TotNhalos = header.TotNhalos;
    }

  for(int k = 0; k < 2; k++)
    n_type[k] = ntot_type[k] = 0;

  {
    ntot_type[0]        = header.Ntrees;
    long long n_in_file = header.Ntrees;
    int ntask           = lastTask - readTask + 1;
    int n_for_this_task = n_in_file / ntask;
    if((ThisTask - readTask) < (n_in_file % ntask))
      n_for_this_task++;
    n_type[0] = n_for_this_task;

    if(nstart)
      memmove(&MergerTree->TreeTable[n_for_this_task], &MergerTree->TreeTable[0], MergerTree->Ntrees * sizeof(halotrees_table));
  }

  {
    ntot_type[1]        = header.Nhalos;
    long long n_in_file = header.Nhalos;
    int ntask           = lastTask - readTask + 1;
    int n_for_this_task = n_in_file / ntask;
    if((ThisTask - readTask) < (n_in_file % ntask))
      n_for_this_task++;
    n_type[1] = n_for_this_task;

    if(nstart)
      memmove(&MergerTree->HaloIDdata[n_for_this_task], &MergerTree->HaloIDdata[0], MergerTree->Nhalos * sizeof(treehalo_ids_type));
  }

  if(nstart)
    *nstart = 0;
}

void readtrees_mbound_io::read_header_fields(const char *fname)
{
  memset(&header, 0, sizeof(header));

  hid_t hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t handle    = my_H5Gopen(hdf5_file, "/Header");

  read_scalar_attribute(handle, "Ntrees_ThisFile", &header.Ntrees, H5T_NATIVE_UINT64);
  read_scalar_attribute(handle, "Ntrees_Total", &header.TotNtrees, H5T_NATIVE_UINT64);

  read_scalar_attribute(handle, "Nhalos_ThisFile", &header.Nhalos, H5T_NATIVE_UINT64);
  read_scalar_attribute(handle, "Nhalos_Total", &header.TotNhalos, H5T_NATIVE_UINT64);

  read_scalar_attribute(handle, "NumFiles", &header.num_files, H5T_NATIVE_INT);

  read_scalar_attribute(handle, "LastSnapShotNr", &header.lastsnapshotnr, H5T_NATIVE_INT);

  my_H5Gclose(handle, "/Header");
  my_H5Fclose(hdf5_file, fname);
}

void readtrees_mbound_io::get_datagroup_name(int type, char *buf)
{
  switch(type)
    {
      case 0:
        snprintf(buf, MAXLEN_PATH, "/TreeTable");
        break;
      case 1:
        snprintf(buf, MAXLEN_PATH, "/TreeHalos");
        break;
      default:
        Terminate("wrong group");
        break;
    }
}

int readtrees_mbound_io::get_type_of_element(int index)
{
  /* empty */
  return 0;
}

void readtrees_mbound_io::set_type_of_element(int index, int type)
{ /* empty */
}

void *readtrees_mbound_io::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_TT:
        return (void *)(MergerTree->TreeTable + index);
      case A_TID:
        return (void *)(MergerTree->HaloIDdata + index);
      default:
        Terminate("strange, we don't expect to get here");
    }
}

#endif
