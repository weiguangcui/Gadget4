/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  io_readtrees_mbound.h
 *
 *  \brief definition of I/O class for most-bound particles belonging to merger trees
 */

#ifndef READTREES_MBOUND_IO_H
#define READTREES_MBOUND_IO_H

#include "gadgetconfig.h"

#ifdef MERGERTREE

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../fof/fof.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

class readtrees_mbound_io : public IO_Def
{
 private:
  mergertree *MergerTree;

 public:
  readtrees_mbound_io(mergertree *MergerTree_ptr, MPI_Comm comm, int format);

  void read_trees_mostbound(void);

  /* supplied virtual functions */
  void fill_file_header(int writeTask, int lastTask, long long *nloc_part, long long *npart);
  void read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *nloc_part, long long *npart,
                        int *nstart);
  void get_datagroup_name(int grnr, char *gname);
  void write_header_fields(hid_t);
  void read_header_fields(const char *fname);
  void read_increase_numbers(int type, int n_for_this_task);
  int get_filenr_from_header(void);
  void set_filenr_in_header(int);
  void *get_base_address_of_structure(enum arrays array, int index);
  int get_type_of_element(int index);
  void set_type_of_element(int index, int type);

  /** Header for the standard file format.
   */

  struct io_header
  {
    long long Nhalos;
    long long TotNhalos;

    long long Ntrees;
    long long TotNtrees;

    int num_files;
    int lastsnapshotnr;
  };
  io_header header;

  typedef mergertree::treehalo_ids_type treehalo_ids_type;
};

#endif

#endif /* HALOTREES_IO_H */
