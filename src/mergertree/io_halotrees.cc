/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  io_halotrees.cc
 *
 *  \brief routines to stores the constructed trees made of subhalos
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
#include "../mergertree/io_halotrees.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

halotrees_io::halotrees_io(mergertree *MergerTree_ptr, MPI_Comm comm, int format) : IO_Def(comm, format)
{
  MergerTree = MergerTree_ptr;

  this->N_IO_Fields  = 0;
  this->N_DataGroups = 3;
  this->header_size  = sizeof(header);
  this->header_buf   = &header;
  this->type_of_file = FILE_IS_TREECAT;
  snprintf(this->info, MAXLEN_PATH, "MERGERTREE: writing mergertrees");

  /* overview table for trees in the file */

  init_field("MTRL", "Length", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_TT, &MergerTree->TreeTable[0].HaloCount, NULL, TREELENGTH, 0, 0,
             0, 0, 0, 0, 0);
  init_field("MTRS", "StartOffset", MEM_INT64, FILE_INT64, SKIP_ON_READ, 1, A_TT, &MergerTree->TreeTable[0].FirstHalo, NULL,
             TREELENGTH, 0, 0, 0, 0, 0, 0, 0);
  init_field("MTRI", "TreeID", MEM_INT64, FILE_INT64, SKIP_ON_READ, 1, A_TT, &MergerTree->TreeTable[0].TreeID, NULL, TREELENGTH, 0, 0,
             0, 0, 0, 0, 0);

  /* link pointers of each subhalo in the trees */

  init_field("TMPR", "TreeMainProgenitor", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeMainProgenitor, NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0, true);
  init_field("TDES", "TreeDescendant", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeDescendant, NULL, TREEHALOS,
             0, 0, 0, 0, 0, 0, 0, true);
  init_field("TFPR", "TreeFirstProgenitor", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeFirstProgenitor, NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0, true);
  init_field("TNPR", "TreeNextProgenitor", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeNextProgenitor, NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0, true);
  init_field("TFHF", "TreeFirstHaloInFOFgroup", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeFirstHaloInFOFgroup,
             NULL, TREEHALOS, 0, 0, 0, 0, 0, 0, 0, true);
  init_field("TNHF", "TreeNextHaloInFOFgroup", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeNextHaloInFOFgroup,
             NULL, TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("TPRO", "TreeProgenitor", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeProgenitor, NULL, TREEHALOS,
             0, 0, 0, 0, 0, 0, 0, true);
  init_field("TFDE", "TreeFirstDescendant", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeFirstDescendant, NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0, true);
  init_field("TNDE", "TreeNextDescendant", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeNextDescendant, NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0, true);

  /* properties of each subhalo in the trees */

  init_field("SLEN", "SubhaloLen", mem_len_type, file_len_type, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].SubProp.Len, NULL,
             SUBGROUPS, 0, 0, 0, 0, 0, 0, 0, true);
  init_field("MASS", "SubhaloMass", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].SubProp.Mass, NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("FMC2", "Group_M_Crit200", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].SubProp.M_Crit200,
             NULL, TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SPOS", "SubhaloPos", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 3, A_H, &MergerTree->Halos[0].SubProp.Pos[0], NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SVEL", "SubhaloVel", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 3, A_H, &MergerTree->Halos[0].SubProp.Vel[0], NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SSPI", "SubhaloVelDisp", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].SubProp.SubVelDisp,
             NULL, TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SVMX", "SubhaloVmax", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].SubProp.SubVmax, NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SVRX", "SubhaloVmaxRad", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].SubProp.SubVmaxRad,
             NULL, TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SHMR", "SubhaloHalfmassRad", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_H,
             &MergerTree->Halos[0].SubProp.SubHalfMassRad, NULL, TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SSPI", "SubhaloSpin", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 3, A_H, &MergerTree->Halos[0].SubProp.Spin[0], NULL,
             TREEHALOS, 0, 0, 0, 0, 0, 0, 0);
  init_field("SIDM", "SubhaloIDMostbound", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, READ_IF_PRESENT, 1, A_H,
             &MergerTree->Halos[0].SubProp.SubMostBoundID, NULL, TREEHALOS, 0, 0, 0, 0, 0, 0, 0, true);

  /* where this subhalo came from */

  init_field("TSNP", "SnapNum", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].SnapNum, NULL, TREEHALOS, 0, 0, 0, 0, 0,
             0, 0, true);
  init_field("TSNR", "SubhaloNr", MEM_INT64, FILE_INT64, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].SubhaloNr, NULL, TREEHALOS, 0, 0,
             0, 0, 0, 0, 0, true);
  init_field("GRNR", "GroupNr", MEM_INT64, FILE_INT64, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].GroupNr, NULL, TREEHALOS, 0, 0, 0,
             0, 0, 0, 0, true);

  /* the fields TreeID and TreeIndex are in principle redundant but kept for convenience */

  init_field("TRID", "TreeID", MEM_INT64, FILE_INT64, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeID, NULL, TREEHALOS, 0, 0, 0, 0,
             0, 0, 0, true);
  init_field("TRIX", "TreeIndex", MEM_INT, FILE_INT, SKIP_ON_READ, 1, A_H, &MergerTree->Halos[0].TreeIndex, NULL, TREEHALOS, 0, 0, 0,
             0, 0, 0, 0, true);

  /**** output times */

  init_field("REDS", "Redshift", MEM_DOUBLE, FILE_DOUBLE, SKIP_ON_READ, 1, A_CT, &MergerTree->CatTimes[0].Redshift, NULL, TREETIMES, 0,
             0, 0, 0, 0, 0, 0);

  init_field("OUTT", "Time", MEM_DOUBLE, FILE_DOUBLE, SKIP_ON_READ, 1, A_CT, &MergerTree->CatTimes[0].Time, NULL, TREETIMES, 0, 0, 0,
             0, 0, 0, 0);
}

void halotrees_io::halotrees_save_trees(void)
{
  char buf[MAXLEN_PATH_EXTRA];

  /* write trees */
  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s/treedata", All.OutputDir);
          mkdir(buf, 02755);
        }
      MPI_Barrier(Communicator);
    }

  if(All.NumFilesPerSnapshot > 1)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/treedata/%s", All.OutputDir, "trees");
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "trees");

  write_multiple_files(buf, All.NumFilesPerSnapshot);
}

void halotrees_io::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type)
{
  /* determine group/id numbers of each type in file */

  n_type[0] = MergerTree->Ntrees;
  n_type[1] = MergerTree->Nhalos;
  if(ThisTask == writeTask)
    n_type[2] = MergerTree->LastSnapShotNr + 1;
  else
    n_type[2] = 0;

  if(ThisTask == writeTask)
    {
      for(int n = 0; n < N_DataGroups; n++)
        ntot_type[n] = n_type[n];

      for(int task = writeTask + 1; task <= lastTask; task++)
        {
          long long nn[N_DataGroups];
          MPI_Recv(&nn[0], N_DataGroups, MPI_LONG_LONG, task, TAG_LOCALN, Communicator, MPI_STATUS_IGNORE);
          for(int n = 0; n < N_DataGroups; n++)
            ntot_type[n] += nn[n];
        }

      for(int task = writeTask + 1; task <= lastTask; task++)
        MPI_Send(&ntot_type[0], N_DataGroups, MPI_LONG_LONG, task, TAG_N, Communicator);
    }
  else
    {
      MPI_Send(&n_type[0], N_DataGroups, MPI_LONG_LONG, writeTask, TAG_LOCALN, Communicator);
      MPI_Recv(&ntot_type[0], N_DataGroups, MPI_LONG_LONG, writeTask, TAG_N, Communicator, MPI_STATUS_IGNORE);
    }

  /* fill file header */
  header.Ntrees = ntot_type[0];
  header.Nhalos = ntot_type[1];

  header.TotNtrees = MergerTree->TotNtrees;
  header.TotNhalos = MergerTree->TotNhalos;

  header.lastsnapshotnr = MergerTree->LastSnapShotNr;

  header.num_files = All.NumFilesPerSnapshot;
}

void halotrees_io::write_header_fields(hid_t handle)
{
  write_scalar_attribute(handle, "Ntrees_ThisFile", &header.Ntrees, H5T_NATIVE_UINT64);
  write_scalar_attribute(handle, "Ntrees_Total", &header.TotNtrees, H5T_NATIVE_UINT64);

  write_scalar_attribute(handle, "Nhalos_ThisFile", &header.Nhalos, H5T_NATIVE_UINT64);
  write_scalar_attribute(handle, "Nhalos_Total", &header.TotNhalos, H5T_NATIVE_UINT64);

  write_scalar_attribute(handle, "NumFiles", &header.num_files, H5T_NATIVE_INT);

  write_scalar_attribute(handle, "LastSnapShotNr", &header.lastsnapshotnr, H5T_NATIVE_INT);
}

int halotrees_io::get_filenr_from_header(void) { return header.num_files; }

void halotrees_io::set_filenr_in_header(int numfiles) { header.num_files = numfiles; }

void halotrees_io::read_increase_numbers(int type, int n_for_this_task)
{
  switch(type)
    {
      case 0:
        MergerTree->Ntrees += n_for_this_task;
        break;
      case 1:
        MergerTree->Nhalos += n_for_this_task;
        break;
      case 2:
        MergerTree->CatTimes += n_for_this_task;
        break;

      default:
        Terminate("wrong group");
        break;
    }
}

void halotrees_io::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *n_type, long long *ntot_type,
                                    int *nstart)
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
      memmove(&MergerTree->Halos[n_for_this_task], &MergerTree->Halos[0], MergerTree->Nhalos * sizeof(treehalo_type));
  }

  if(nstart)
    *nstart = 0;
}

void halotrees_io::read_header_fields(const char *fname)
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

void halotrees_io::get_datagroup_name(int type, char *buf)
{
  switch(type)
    {
      case 0:
        snprintf(buf, MAXLEN_PATH, "/TreeTable");
        break;
      case 1:
        snprintf(buf, MAXLEN_PATH, "/TreeHalos");
        break;
      case 2:
        snprintf(buf, MAXLEN_PATH, "/TreeTimes");
        break;
      default:
        Terminate("wrong group");
        break;
    }
}

int halotrees_io::get_type_of_element(int index)
{
  /* empty */
  return 0;
}

void halotrees_io::set_type_of_element(int index, int type)
{ /* empty */
}

void *halotrees_io::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_TT:
        return (void *)(MergerTree->TreeTable + index);
      case A_H:
        return (void *)(MergerTree->Halos + index);
      case A_CT:
        return (void *)(MergerTree->CatTimes + index);
      default:
        Terminate("strange, we don't expect to get here");
    }
}

#endif
