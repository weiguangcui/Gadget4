/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  io_treelinks.cc
 *
 *  \brief routines for the I/O needed for the merger treelink files
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
#include "../mergertree/io_treelinks.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

treelinks_io::treelinks_io(mergertree *MergerTree_ptr, MPI_Comm comm, int format) : IO_Def(comm, format)
{
  MergerTree = MergerTree_ptr;

  this->N_IO_Fields  = 0;
  this->N_DataGroups = 1;
  this->header_size  = sizeof(header);
  this->header_buf   = &header;
  this->type_of_file = FILE_IS_TREELINK;
  snprintf(this->info, MAXLEN_PATH, "TREELINK: writing treelink information");

  init_field("TRNR", "TreeID", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_TL, &MergerTree->TreeLink[0].TreeID, NULL, TREELINK, 0, 0,
             0, 0, 0, 0, 0, true);

  init_field("TRIX", "TreeIndex", MEM_INT, FILE_INT, READ_IF_PRESENT, 1, A_TL, &MergerTree->TreeLink[0].TreeIndex, NULL, TREELINK, 0,
             0, 0, 0, 0, 0, 0, true);
}

void treelinks_io::treelinks_save(int num)
{
  char buf[MAXLEN_PATH_EXTRA];

  /* write treelink info */
  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s/groups_%03d", All.OutputDir, num);
          mkdir(buf, 02755);
        }
      MPI_Barrier(Communicator);
    }

  if(All.NumFilesPerSnapshot > 1)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/groups_%03d/%s_%03d", All.OutputDir, num, "subhalo_treelink", num);
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, "subhalo_treelink", num);

  write_multiple_files(buf, All.NumFilesPerSnapshot);
}

void treelinks_io::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type)
{
  /* determine group/id numbers of each type in file */

  n_type[0] = Nsubhalos;

  if(ThisTask == writeTask)
    {
      for(int n = 0; n < 1; n++)
        ntot_type[n] = n_type[n];

      for(int task = writeTask + 1; task <= lastTask; task++)
        {
          long long nn[1];
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
  header.TotNsubhalos = TotNsubhalos;
  header.num_files    = All.NumFilesPerSnapshot;
}

void treelinks_io::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *n_type, long long *ntot_type,
                                    int *nstart)
{
  /* empty */
}

void treelinks_io::write_header_fields(hid_t handle)
{
  write_scalar_attribute(handle, "Nsubhalos_ThisFile", &header.Nsubhalos, H5T_NATIVE_UINT64);

  write_scalar_attribute(handle, "Nsubhalos_Total", &header.TotNsubhalos, H5T_NATIVE_UINT64);

  write_scalar_attribute(handle, "NumFiles", &header.num_files, H5T_NATIVE_INT);
}

void treelinks_io::read_header_fields(const char *fname)
{ /* empty */
}

int treelinks_io::get_filenr_from_header(void) { return header.num_files; }

void treelinks_io::set_filenr_in_header(int numfiles) { header.num_files = numfiles; }

void treelinks_io::read_increase_numbers(int type, int n_for_this_task)
{ /* empty */
}

void treelinks_io::get_datagroup_name(int type, char *buf)
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

int treelinks_io::get_type_of_element(int index)
{
  /* empty */
  return 0;
}

void treelinks_io::set_type_of_element(int index, int type)
{ /* empty */
}

void *treelinks_io::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_TL:
        return (void *)(MergerTree->TreeLink + index);
      default:
        Terminate("strange, we don't expect to get here");
    }
}

#endif
