/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  io_treelinks.h
 *
 *  \brief declaration of the class used for I/O of the merger treelink files
 */

#ifndef TREELINKS_IO_H
#define TREELINKS_IO_H

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

class treelinks_io : public IO_Def
{
 private:
  mergertree *MergerTree;

 public:
  treelinks_io(mergertree *MergerTree_ptr, MPI_Comm comm, int format);

  void treelinks_save(int num);

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
    long long Nsubhalos;
    long long TotNsubhalos;
    int num_files;
  };
  io_header header;

  int Nsubhalos;
  long long TotNsubhalos;

 private:
};

#endif

#endif /* TREELINKS_IO_H */
