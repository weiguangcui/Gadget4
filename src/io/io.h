/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file io.h
 *
 * \brief declarations of I/O enumerations and a base class for I/O in Gadget
 */

#ifndef IO_H
#define IO_H

#include "gadgetconfig.h"

#include <hdf5.h>
#ifdef LIGHTCONE_PARTICLES
#include <chealpix.h>
#endif

#include "../data/simparticles.h"
#include "../fof/fof.h"
#include "../io/io_streamcount.h"
#include "../mpi_utils/setcomm.h"

#define LABEL_LEN 4
#define DATASETNAME_LEN 256

enum arrays
{
  A_NONE,
  A_SPHP,
  A_P,
  A_PS,
  A_G,
  A_S,
  A_ID,
  A_DESC,
  A_PROG,
  A_MTRP,
  A_H,
  A_TT,
  A_CT,
  A_TL,
  A_LC,
  A_MM,
  A_IDS,
  A_TID,
};

enum file_contents
{
  FILE_IS_SNAPSHOT,
  FILE_IS_GROUPCAT,
  FILE_IS_DESCCAT,
  FILE_IS_PROGCAT,
  FILE_IS_TREECAT,
  FILE_IS_TREELINK,
  FILE_IS_LIGHTCONE,
  FILE_IS_MASSMAP,
  FILE_IS_GALSNAP
};

enum types_in_file
{
  FILE_NONE,
  FILE_INT,
  FILE_INT64,
  FILE_MY_IO_FLOAT,
  FILE_HALF,
  FILE_MY_ID_TYPE,
  FILE_MY_INTPOS_TYPE,
  FILE_DOUBLE,
  FILE_FLOAT,
};

enum types_in_memory
{
  MEM_INT,
  MEM_INT64,
  MEM_MY_ID_TYPE,
  MEM_MY_INTPOS_TYPE,
  MEM_FLOAT,
  MEM_DOUBLE,
  MEM_MY_FLOAT,
  MEM_MY_DOUBLE,
  MEM_MY_FILEOFFSET,
};

#ifdef FOF_ALLOW_HUGE_GROUPLENGTH
const types_in_memory mem_len_type = MEM_INT64;
const types_in_file file_len_type  = FILE_INT64;
#else
const types_in_memory mem_len_type = MEM_INT;
const types_in_file file_len_type  = FILE_INT;
#endif

enum e_typelist
{
  GAS_ONLY      = 1,
  STARS_ONLY    = 16,
  GAS_AND_STARS = 17,
  ALL_TYPES     = ((1 << NTYPES) - 1),
  MASS_BLOCK    = -1,
  AGE_BLOCK     = -2,
  Z_BLOCK       = -3,
  GROUPS        = -4,
  SUBGROUPS     = -5,
  ID_BLOCK      = -6,
  PREVSUBS      = -7,
  CURRSUBS      = -8,
  TREELENGTH    = -9,
  TREEHALOS     = -10,
  TREELINK      = -11,
  MASSMAPS      = -12,
  TREETABLE     = -13,
  TREETIMES     = -14,
  GALSNAPS      = -16,
  HEALPIXTAB    = -17,
};

enum read_flags
{
  READ_IF_PRESENT,
  SKIP_ON_READ
};

class IO_Def : public io_streamcount, public setcomm
{
 public:
  IO_Def(MPI_Comm comm, int format) : setcomm(comm)
  {
    determine_compute_nodes();
    file_format = format;
  }

  virtual ~IO_Def();

  int N_DataGroups  = 0;
  int N_IO_Fields   = 0;
  int Max_IO_Fields = 0;

  /* functions that need to be provided by the specific module */

  virtual void fill_file_header(int writeTask, int lastTask, long long *nloc_part, long long *npart) = 0;
  virtual void read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *nloc_part, long long *npart,
                                int *nstart)                                                         = 0;
  virtual void get_datagroup_name(int grnr, char *gname)                                             = 0;
  virtual void write_header_fields(hid_t)                                                            = 0;
  virtual void read_header_fields(const char *fname)                                                 = 0;
  virtual void read_increase_numbers(int type, int n_for_this_task)                                  = 0;
  virtual int get_filenr_from_header(void)                                                           = 0;
  virtual void set_filenr_in_header(int)                                                             = 0;
  virtual void *get_base_address_of_structure(enum arrays array, int index)                          = 0;
  virtual int get_type_of_element(int index)                                                         = 0;
  virtual void set_type_of_element(int index, int type)                                              = 0;

  void init_field(const char *label, const char *datasetname, enum types_in_memory type_in_memory,
                  enum types_in_file type_in_file_output, enum read_flags read_flag, int values_per_block, enum arrays array,
                  void *pointer_to_field, void (*io_func)(IO_Def *, int, int, void *, int), int typelist_bitmask, int hasunits,
                  double a, double h, double L, double M, double V, double c, bool compression_on = false);

  int find_files(const char *fname, const char *fname_multiple);
  void read_files_driver(const char *fname, int rep, int numfiles);
  void write_multiple_files(char *fname, int numfilesperdump, int append_flag = 0, int chunk_size = 0);
  void write_compile_time_options_in_hdf5(hid_t handle);
  void read_segment(const char *fname, int type, long long offset, long long count, int numfiles);
  void read_single_file_segment(const char *fname, int filenr, int type, long long offset, unsigned long long count,
                                long long storage_offset, int numfiles);

  void alloc_and_read_ntype_in_files(const char *fname, int num_files);

  size_t header_size;
  void *header_buf;

  long long *ntype_in_files;
  char info[MAXLEN_PATH];

#if defined(MERGERTREE)
  typedef fof<simparticles>::treehalo_t treehalo_type;
#endif

  /*
   * variables for new input/output functionality
   */

  enum file_contents type_of_file;

 private:
  struct IO_Field
  {
    enum types_in_memory type_in_memory;
    enum types_in_file type_in_file_output;
    enum read_flags read_flag;
    int values_per_block;
    int write_block;
    int read_block;
    char label[LABEL_LEN + 1];
    char datasetname[DATASETNAME_LEN + 1];
    void (*io_func)(IO_Def *, int, int, void *, int);
    int typelist;
    bool compression_on;
    enum arrays array;
    size_t offset;

    char hasunit;
    double a;
    double h;
    double L;
    double M;
    double V;
    double c;
  };
  IO_Field *IO_Fields;

  void write_file(char *fname, int writeTask, int lastTask, void *CommBuffer, int numfilesperdump, int chunksize);
  void read_file(const char *fname, int filenr, int readTask, int lastTask, void *CommBuffer);
  void append_file(char *fname, int writeTask, int lastTask, void *CommBuffer, int numfilesperdump, int chunksize);

  int get_values_per_blockelement(int blocknr);
  void get_dataset_name(int blocknr, char *buf);
  long long get_particles_in_block(int blocknr, long long *npartinfile, int *typelist);
  int get_bytes_per_memory_blockelement(int blocknr, int mode);
  hid_t get_hdf5_outputtype_of_block(int blocknr);
  hid_t get_hdf5_memorytype_of_block(int blocknr);
  void get_Tab_IO_Label(int blocknr, char *label);
  void type_cast_data(char *src, int src_bytes_per_element, char *target, int target_bytes_per_element, int len, int blocknr);
  void distribute_file(int nfiles, int *filenr, int *master, int *last);
  void share_particle_number_in_file(const char *fname, int filenr, int readTask, int lastTask);
  void find_block(char *label, FILE *fd);
  void fill_write_buffer(int blocknr, int *startindex, int pc, int type, void *CommBuffer);
  void empty_read_buffer(int blocknr, int offset, int pc, int type, long long nprevious, void *CommBuffer);

  void write_dataset_attributes(hid_t hdf5_dataset, int blocknr);
  void write_parameters_attributes_in_hdf5(hid_t handle);
  void rename_file_to_bak_if_it_exists(char *fname);

  void polling(int numfilesperdump);

  int files_started;
  int files_completed;
  int file_format;

  struct seq_data
  {
    int thistask;
    int rankinnode;
    int thisnode;
    bool operator<(const seq_data &b) const
    {
      if(rankinnode < b.rankinnode)
        return true;
      if(rankinnode > b.rankinnode)
        return false;
      if(thisnode < b.thisnode)
        return true;
      if(thisnode > b.thisnode)
        return false;
      return thistask < b.thistask;
    }
  };
  seq_data *seq;
};

#ifdef GADGET2_HEADER
#define NTYPES_HEADER 6
#else
#define NTYPES_HEADER NTYPES
#endif

#define FLAG_ZELDOVICH_ICS 1
#define FLAG_SECOND_ORDER_ICS 2

#endif
