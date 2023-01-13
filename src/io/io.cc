/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file io.cc
 *
 * \brief main routines for driving I/O in Gadget's three file formats for snapshots, group catalogues etc.
 */

#include "gadgetconfig.h"

#include <errno.h>
#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <algorithm>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../fof/fof.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../io/parameters.h"
#include "../lightcone/lightcone.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

#define HALF_ROUND_STYLE 1
#include "../half/half.hpp"
using half_float::half;

/* local functions */

IO_Def::~IO_Def()
{
  if(N_IO_Fields > 0)
    Mem.myfree_movable(IO_Fields);
}

/*!
 * \param field Specifies the field as an enumeration type iofields (io_private.h), e.g. IO_POS.
 * \param label The label of the dataset (4 characters, for obsolete old format=2)
 * \param datasetname The name of the hdf5 dataset (maximum 256 characters, for format=3)
 * \param type_in_memory The type of the field in the memory
 * \param type_in_file_output The output type in the file if it is written (upon read, the type found will be converted to the memory
 * type) \param read_flag This flags whether the field should be ignored upon read on (use SKIP_ON_READ, else READ_IF_PRESENT). \param
 * values_per_block The number of values per field, e.g. 1 for mass, 3 for velocities \param array The array in which the value is
 * stored \param pointer_to_field A Pointer to the field in one of the global arrays, e.g. &SphP[0].Density, or &P[0].Vel[0] \param
 * io_func Alternatively, if the value to output/input is not a simple field, you can define a function which handles i/o \param
 * typelist_bitmask Specifies for which particle type the field is present, e.g. 1+2+8 => field present for particle types 0,1,3 (or
 * use ALL_TYPES, GAS_ONLY,...)
 */
void IO_Def::init_field(const char *label, const char *datasetname, enum types_in_memory type_in_memory,
                        enum types_in_file type_in_file_output, enum read_flags read_flag, int values_per_block, enum arrays array,
                        void *pointer_to_field, void (*io_func)(IO_Def *, int, int, void *, int), int typelist_bitmask, int flagunit,
                        double a, double h, double L, double M, double V, double c, bool compression_on)
{
  const int alloc_step = 5;

  if(values_per_block < 1)  // if we have no values, we don't register the field
    return;

  if(Max_IO_Fields == 0)
    {
      IO_Fields     = (IO_Field *)Mem.mymalloc_movable(&IO_Fields, "IO_Fields", alloc_step * sizeof(IO_Field));
      Max_IO_Fields = alloc_step;
    }
  else if(Max_IO_Fields == N_IO_Fields)
    {
      Max_IO_Fields = ((Max_IO_Fields / alloc_step) + 1) * alloc_step;
      IO_Fields     = (IO_Field *)Mem.myrealloc_movable(IO_Fields, Max_IO_Fields * sizeof(IO_Field));
    }

  int n           = N_IO_Fields++;
  IO_Field *field = &IO_Fields[n];

  strncpy(field->label, label, LABEL_LEN);
  field->label[LABEL_LEN] = 0;

  strncpy(field->datasetname, datasetname, DATASETNAME_LEN);
  field->datasetname[DATASETNAME_LEN] = 0;

  field->type_in_memory      = type_in_memory;
  field->type_in_file_output = type_in_file_output;
  field->read_flag           = read_flag;
  field->values_per_block    = values_per_block;
  field->typelist            = typelist_bitmask;
#ifdef ALLOW_HDF5_COMPRESSION
  field->compression_on = compression_on;
#else
  field->compression_on = false;
#endif

  field->array   = array;
  field->io_func = io_func;

  if(array == A_NONE)
    {
      field->offset = 0;
    }
  else
    {
      field->offset = (size_t)pointer_to_field - (size_t)get_base_address_of_structure(array, 0);
    }

  field->hasunit = flagunit;
  field->a       = a;
  field->h       = h;
  field->L       = L;
  field->M       = M;
  field->V       = V;
  field->c       = c;
}

/*! \brief This function determines on how many files a given snapshot or group/desc catalogue is distributed.
 *
 *  \param fname file name of the snapshot as given in the parameter file
 */
int IO_Def::find_files(const char *fname, const char *fname_multiple)
{
  FILE *fd;
  char buf[MAXLEN_PATH_EXTRA], buf1[MAXLEN_PATH_EXTRA];
  int dummy, files_found = 0;

  if(file_format == FILEFORMAT_HDF5)
    {
      snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d.hdf5", fname_multiple, 0);
      snprintf(buf1, MAXLEN_PATH_EXTRA, "%s.hdf5", fname);
    }
  else
    {
      snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d", fname_multiple, 0);
      snprintf(buf1, MAXLEN_PATH_EXTRA, "%s", fname);
    }

  memset(header_buf, 0, header_size);

  if(ThisTask == 0)
    {
      if((fd = fopen(buf, "r")))
        {
          if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
            {
              if(file_format == FILEFORMAT_LEGACY2)
                {
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                }

              my_fread(&dummy, sizeof(dummy), 1, fd);
              my_fread(header_buf, header_size, 1, fd);
              my_fread(&dummy, sizeof(dummy), 1, fd);
            }
          fclose(fd);

          if(file_format == FILEFORMAT_HDF5)
            read_header_fields(buf);

          files_found = 1;
        }
    }

  MPI_Bcast(header_buf, header_size, MPI_BYTE, 0, Communicator);
  MPI_Bcast(&files_found, 1, MPI_INT, 0, Communicator);

  if(get_filenr_from_header() > 0)
    return get_filenr_from_header();

  if(ThisTask == 0)
    {
      if((fd = fopen(buf1, "r")))
        {
          if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
            {
              if(file_format == FILEFORMAT_LEGACY2)
                {
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                  my_fread(&dummy, sizeof(dummy), 1, fd);
                }

              my_fread(&dummy, sizeof(dummy), 1, fd);
              my_fread(header_buf, header_size, 1, fd);
              my_fread(&dummy, sizeof(dummy), 1, fd);
            }
          fclose(fd);

          if(file_format == FILEFORMAT_HDF5)
            read_header_fields(buf1);

          set_filenr_in_header(1);

          files_found = 1;
        }
    }

  MPI_Bcast(header_buf, header_size, MPI_BYTE, 0, Communicator);
  MPI_Bcast(&files_found, 1, MPI_INT, 0, Communicator);

  if(get_filenr_from_header() > 0)
    return get_filenr_from_header();

  if(files_found != 0)
    Terminate("Have found IC files, but number of files in header seems to be zero\n");

  Terminate("\nCan't find files, neither as '%s'\nnor as '%s'\n", buf, buf1);

  return 0;
}

void IO_Def::read_files_driver(const char *fname, int rep, int num_files)
{
  if(rep == 0)
    {
      ntype_in_files =
          (long long *)Mem.mymalloc_movable(&ntype_in_files, "ntype_in_files", num_files * N_DataGroups * sizeof(long long));
      memset(ntype_in_files, 0, num_files * N_DataGroups * sizeof(long long));
    }

  void *CommBuffer = Mem.mymalloc("CommBuffer", COMMBUFFERSIZE);

  int rest_files = num_files;

  while(rest_files > NTask)
    {
      char buf[MAXLEN_PATH_EXTRA];

      snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d", fname, ThisTask + (rest_files - NTask));
      if(file_format == FILEFORMAT_HDF5)
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));

      int ngroups = NTask / All.MaxFilesWithConcurrentIO;
      if((NTask % All.MaxFilesWithConcurrentIO))
        ngroups++;
      int groupMaster = (ThisTask / ngroups) * ngroups;

      for(int gr = 0; gr < ngroups; gr++)
        {
          if(ThisTask == (groupMaster + gr)) /* ok, it's this processor's turn */
            {
              if(rep == 0)
                share_particle_number_in_file(buf, ThisTask + (rest_files - NTask), ThisTask, ThisTask);
              else
                read_file(buf, ThisTask + (rest_files - NTask), ThisTask, ThisTask, CommBuffer);
            }
          MPI_Barrier(Communicator);
        }

      rest_files -= NTask;
    }

  if(rest_files > 0)
    {
      int masterTask, filenr, lastTask;

      distribute_file(rest_files, &filenr, &masterTask, &lastTask);

      char buf[MAXLEN_PATH_EXTRA];

      if(num_files > 1)
        {
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d", fname, filenr);
          if(file_format == FILEFORMAT_HDF5)
            snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d.hdf5", fname, filenr);
        }
      else
        {
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s", fname);
          if(file_format == FILEFORMAT_HDF5)
            snprintf(buf, MAXLEN_PATH_EXTRA, "%s.hdf5", fname);
        }

      int ngroups = rest_files / All.MaxFilesWithConcurrentIO;
      if((rest_files % All.MaxFilesWithConcurrentIO))
        ngroups++;

      for(int gr = 0; gr < ngroups; gr++)
        {
          if((filenr / All.MaxFilesWithConcurrentIO) == gr) /* ok, it's this processor's turn */
            {
              if(rep == 0)
                share_particle_number_in_file(buf, filenr, masterTask, lastTask);
              else
                read_file(buf, filenr, masterTask, lastTask, CommBuffer);
            }
          MPI_Barrier(Communicator);
        }
    }

  Mem.myfree(CommBuffer);

  if(rep == 0)
    MPI_Allreduce(MPI_IN_PLACE, ntype_in_files, num_files * N_DataGroups, MPI_LONG_LONG, MPI_MAX, Communicator);
  else
    {
      /* we are done */
      Mem.myfree_movable(ntype_in_files);
      ntype_in_files = NULL;
    }
}

/*! \brief This function distributes the particle numbers in the file fname
 *  to tasks 'readTask' to 'lastTask', and calculates the number of particles each task gets.
 *
 *  \param fname filename to be read
 *  \param readTask task responsible for reading the file fname
 *  \param lastTask last Task which gets data contained in the file
 *  \param readTypes readTypes is a bitfield that
 *  determines what particle types to read, only if the bit
 *  corresponding to a particle type is set, the corresponding data is
 *  loaded, otherwise its particle number is set to zero. (This is
 *  only implemented for HDF5 files.)
 */
void IO_Def::share_particle_number_in_file(const char *fname, int filenr, int readTask, int lastTask)
{
  long long n_type[N_DataGroups], npart[N_DataGroups];
  unsigned int blksize1, blksize2;

  if(ThisTask == readTask)
    {
      if(file_format == FILEFORMAT_HDF5)
        {
          read_header_fields(fname);
        }
      else if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
        {
          FILE *fd = 0;

          if(!(fd = fopen(fname, "r")))
            Terminate("can't open file `%s' for reading initial conditions.\n", fname);

          if(file_format == FILEFORMAT_LEGACY2)
            {
              char label[LABEL_LEN];
              int nextblock;
              my_fread(&blksize1, sizeof(int), 1, fd);
              my_fread(&label, sizeof(char), LABEL_LEN, fd);
              my_fread(&nextblock, sizeof(int), 1, fd);
              printf("%s Reading header => '%c%c%c%c' (%d byte)\n", info, label[0], label[1], label[2], label[3], nextblock);
              my_fread(&blksize2, sizeof(int), 1, fd);
            }

          my_fread(&blksize1, sizeof(int), 1, fd);
          my_fread(header_buf, header_size, 1, fd);
          my_fread(&blksize2, sizeof(int), 1, fd);

#ifdef GADGET2_HEADER
          if(blksize1 != 256 || blksize2 != 256)
            Terminate("incorrect GADGET2 header format, blksize1=%d blksize2=%d  header_size=%d\n", blksize1, blksize2,
                      (int)header_size);
#else
          if(blksize1 != blksize2)
            Terminate("incorrect header format, blksize1=%d blksize2=%d  header_size=%d \n%s \n", blksize1, blksize2, (int)header_size,
                      blksize1 == 256 ? "You may need to set GADGET2_HEADER" : "");
#endif
          fclose(fd);
        }
      else
        Terminate("illegal format");

      for(int task = readTask + 1; task <= lastTask; task++)
        MPI_Ssend(header_buf, header_size, MPI_BYTE, task, TAG_HEADER, Communicator);
    }
  else
    MPI_Recv(header_buf, header_size, MPI_BYTE, readTask, TAG_HEADER, Communicator, MPI_STATUS_IGNORE);

  read_file_header(fname, filenr, readTask, lastTask, n_type, npart, NULL);

  if(ThisTask == readTask)
    {
      mpi_printf("READIC: Reading file `%s' on task=%d and distribute it to %d to %d.\n", fname, ThisTask, readTask, lastTask);
      myflush(stdout);
    }

  for(int type = 0; type < N_DataGroups; type++)
    {
      ntype_in_files[filenr * N_DataGroups + type] = npart[type];

      long long n_in_file = npart[type];
      int ntask           = lastTask - readTask + 1;
      int n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
        n_for_this_task++;

      read_increase_numbers(type, n_for_this_task);
    }
}

/*! \brief This function fills the write buffer with particle data.
 *
 *  New output blocks can in principle be added here.
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param startindex pointer containing the offset in the write buffer
 *  \param pc number of particle to be put in the buffer
 *  \param type particle type
 */
void IO_Def::fill_write_buffer(int blocknr, int *startindex, int pc, int type, void *CommBuffer)
{
  if(blocknr < 0 || blocknr >= N_IO_Fields)
    Terminate("something is wrong here: blocknr=%d N_IO_Fields=%d", blocknr, N_IO_Fields);

  IO_Field *field = &IO_Fields[blocknr];

  int *intp        = (int *)CommBuffer;
  long long *longp = (long long *)CommBuffer;
  float *floatp    = (float *)CommBuffer;
  double *doublep  = (double *)CommBuffer;

  MyIDType *ip          = (MyIDType *)CommBuffer;
  MyFloat *fp           = (MyFloat *)CommBuffer;
  MyDouble *dp          = (MyDouble *)CommBuffer;
  MyIntPosType *intposp = (MyIntPosType *)CommBuffer;

  int pindex = *startindex;

  for(int n = 0; n < pc; pindex++)
    {
      if(type_of_file == FILE_IS_SNAPSHOT && type < NTYPES && get_type_of_element(pindex) != type)
        continue;
#ifdef LIGHTCONE
      if(type_of_file == FILE_IS_LIGHTCONE && type < NTYPES && get_type_of_element(pindex) != type)
        continue;
#endif
      if(field->io_func)
        {
          switch(field->type_in_memory)
            {
              case MEM_INT:
                field->io_func(this, pindex, field->values_per_block, intp, 0);
                intp += field->values_per_block;
                break;
              case MEM_INT64:
                field->io_func(this, pindex, field->values_per_block, longp, 0);
                longp += field->values_per_block;
                break;
              case MEM_MY_ID_TYPE:
                field->io_func(this, pindex, field->values_per_block, ip, 0);
                ip += field->values_per_block;
                break;
              case MEM_MY_INTPOS_TYPE:
                field->io_func(this, pindex, field->values_per_block, intposp, 0);
                intposp += field->values_per_block;
                break;
              case MEM_FLOAT:
                field->io_func(this, pindex, field->values_per_block, floatp, 0);
                floatp += field->values_per_block;
                break;
              case MEM_DOUBLE:
                field->io_func(this, pindex, field->values_per_block, doublep, 0);
                doublep += field->values_per_block;
                break;
              case MEM_MY_FLOAT:
                field->io_func(this, pindex, field->values_per_block, fp, 0);
                fp += field->values_per_block;
                break;
              case MEM_MY_DOUBLE:
                field->io_func(this, pindex, field->values_per_block, dp, 0);
                dp += field->values_per_block;
                break;
              default:
                Terminate("ERROR in fill_write_buffer: Type not found!\n");

                break;
            }
        }
      else
        {
          void *array_pos;

          switch(field->array)
            {
              case A_NONE:
                array_pos = NULL;
                break;

              default:
                array_pos = get_base_address_of_structure(field->array, pindex);
                break;
            }

          for(int k = 0; k < field->values_per_block; k++)
            {
              switch(field->type_in_memory)
                {
                  case MEM_INT:
                    *intp++ = *((int *)((size_t)array_pos + field->offset + k * sizeof(int)));
                    break;

                  case MEM_INT64:
                    *longp++ = *((long long *)((size_t)array_pos + field->offset + k * sizeof(long long)));
                    break;

                  case MEM_MY_ID_TYPE:
                    *ip++ = *((MyIDType *)((size_t)array_pos + field->offset + k * sizeof(MyIDType)));
                    break;

                  case MEM_MY_INTPOS_TYPE:
                    *intposp++ = *((MyIntPosType *)((size_t)array_pos + field->offset + k * sizeof(MyIntPosType)));
                    break;

                  case MEM_FLOAT:
                    *floatp++ = *((float *)((size_t)array_pos + field->offset + k * sizeof(float)));
                    break;

                  case MEM_DOUBLE:
                    *doublep++ = *((double *)((size_t)array_pos + field->offset + k * sizeof(double)));
                    break;

                  case MEM_MY_FLOAT:
                    *fp++ = *((MyFloat *)((size_t)array_pos + field->offset + k * sizeof(MyFloat)));
                    break;

                  case MEM_MY_DOUBLE:
                    *dp++ = *((MyDouble *)((size_t)array_pos + field->offset + k * sizeof(MyDouble)));
                    break;

                  default:
                    Terminate("ERROR in fill_write_buffer: Type not found!\n");
                    break;
                }
            }
        }

      n++;
    }

  *startindex = pindex;
}

/*! \brief This function reads out the io buffer that was filled with particle data.
 *
 * The data in the io buffer is put in the appropriate places of the particle structures.
 *
 * \param blocknr data block present in io buffer
 * \param offset particle corresponding to the first element in io buffer
 * \param pc number of elements in the io buffer
 * \param type If blocknr=IO_POS P[n].Type is set to type
 */
void IO_Def::empty_read_buffer(int blocknr, int offset, int pc, int type, long long nprevious, void *CommBuffer)
{
  IO_Field *field = &IO_Fields[blocknr];

  int *intp        = (int *)CommBuffer;
  long long *longp = (long long *)CommBuffer;
  float *floatp    = (float *)CommBuffer;
  double *doublep  = (double *)CommBuffer;

  MyIDType *ip          = (MyIDType *)CommBuffer;
  MyFloat *fp           = (MyFloat *)CommBuffer;
  MyDouble *dp          = (MyDouble *)CommBuffer;
  MyIntPosType *intposp = (MyIntPosType *)CommBuffer;

  if(field->read_flag != SKIP_ON_READ || field->type_in_memory == MEM_MY_FILEOFFSET)
    {
      for(int n = 0; n < pc; n++)
        {
          if(field->io_func)
            {
              switch(field->type_in_memory)
                {
                  case MEM_INT:
                    field->io_func(this, offset + n, field->values_per_block, intp, 1);
                    intp += field->values_per_block;
                    break;
                  case MEM_INT64:
                    field->io_func(this, offset + n, field->values_per_block, longp, 1);
                    longp += field->values_per_block;
                    break;
                  case MEM_MY_ID_TYPE:
                    field->io_func(this, offset + n, field->values_per_block, ip, 1);
                    ip += field->values_per_block;
                    break;
                  case MEM_MY_INTPOS_TYPE:
                    field->io_func(this, offset + n, field->values_per_block, intposp, 1);
                    intposp += field->values_per_block;
                    break;
                  case MEM_FLOAT:
                    field->io_func(this, offset + n, field->values_per_block, floatp, 1);
                    floatp += field->values_per_block;
                    break;
                  case MEM_DOUBLE:
                    field->io_func(this, offset + n, field->values_per_block, doublep, 1);
                    doublep += field->values_per_block;
                    break;
                  case MEM_MY_FLOAT:
                    field->io_func(this, offset + n, field->values_per_block, fp, 1);
                    fp += field->values_per_block;
                    break;
                  case MEM_MY_DOUBLE:
                    field->io_func(this, offset + n, field->values_per_block, dp, 1);
                    dp += field->values_per_block;
                    break;
                  case MEM_MY_FILEOFFSET:
                    Terminate("undefined");
                    break;
                }
            }
          else
            {
              void *array_pos;
              switch(field->array)
                {
                  case A_NONE:
                    array_pos = 0;
                    break;

                  default:
                    array_pos = get_base_address_of_structure(field->array, offset + n);
                    break;
                }

              for(int k = 0; k < field->values_per_block; k++)
                {
                  switch(field->type_in_memory)
                    {
                      case MEM_INT:
                        *((int *)((size_t)array_pos + field->offset + k * sizeof(int))) = *intp++;
                        break;
                      case MEM_INT64:
                        *((long long *)((size_t)array_pos + field->offset + k * sizeof(long long))) = *longp++;
                        break;
                      case MEM_MY_ID_TYPE:
                        *((MyIDType *)((size_t)array_pos + field->offset + k * sizeof(MyIDType))) = *ip++;
                        break;
                      case MEM_MY_INTPOS_TYPE:
                        *((MyIntPosType *)((size_t)array_pos + field->offset + k * sizeof(MyIntPosType))) = *intposp++;
                        break;
                      case MEM_FLOAT:
                        *((float *)((size_t)array_pos + field->offset + k * sizeof(float))) = *floatp++;
                        break;
                      case MEM_DOUBLE:
                        *((double *)((size_t)array_pos + field->offset + k * sizeof(double))) = *doublep++;
                        break;
                      case MEM_MY_FLOAT:
                        *((MyFloat *)((size_t)array_pos + field->offset + k * sizeof(MyFloat))) = *fp++;
                        break;
                      case MEM_MY_DOUBLE:
                        *((MyDouble *)((size_t)array_pos + field->offset + k * sizeof(MyDouble))) = *dp++;
                        break;
                      case MEM_MY_FILEOFFSET:
                        *((long long *)((size_t)array_pos + field->offset + k * sizeof(long long))) = nprevious++;
                        break;
                      default:
                        Terminate("ERROR: Type not found!\n");
                        break;
                    }
                }
            }
        }
    }
}

void IO_Def::polling(int numfilesperdump)
{
  if(ThisTask == 0)
    if(files_completed < numfilesperdump)
      {
        MPI_Status status;
        int flag;

        /* now check for a completion message  */
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_KEY, Communicator, &flag, &status);

        if(flag)
          {
            int source = status.MPI_SOURCE;

            int dummy;
            MPI_Recv(&dummy, 1, MPI_INT, source, TAG_KEY, Communicator, MPI_STATUS_IGNORE);
            files_completed++;

            if(files_started < numfilesperdump)
              {
                /* send start signal */
                MPI_Ssend(&ThisTask, 1, MPI_INT, seq[files_started++].thistask, TAG_N, Communicator);
              }
          }
      }
}

/* driver routine for outputting multiple files, scheduled in optimal order under the constraint not to write more than a certain
 * number of files simultanuously */
void IO_Def::write_multiple_files(char *fname, int numfilesperdump, int append_flag, int chunksize)
{
  if(ThisTask == 0)
    if(!(seq = (seq_data *)Mem.mymalloc("seq", NTask * sizeof(seq_data))))
      Terminate("can't allocate seq_data");

  void *CommBuffer = Mem.mymalloc("CommBuffer", COMMBUFFERSIZE);

  /* assign processors to output files */
  int filenr, masterTask, lastTask;
  distribute_file(numfilesperdump, &filenr, &masterTask, &lastTask);

  char buf[MAXLEN_PATH_EXTRA];
  if(numfilesperdump > 1)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d", fname, filenr);
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s", fname);

  seq_data seq_loc;
  seq_loc.thistask   = ThisTask;
  seq_loc.rankinnode = RankInThisNode;
  seq_loc.thisnode   = ThisNode;

  if(masterTask != ThisTask)
    seq_loc.thistask = -1;

  MPI_Gather(&seq_loc, sizeof(seq_data), MPI_BYTE, seq, sizeof(seq_data), MPI_BYTE, 0, Communicator);

  if(ThisTask == 0)
    {
      int count = NTask;
      for(int i = 0; i < count; i++)
        {
          if(seq[i].thistask < 0)
            {
              count--;
              seq[i] = seq[count];
              i--;
            }
        }
      if(count != numfilesperdump)
        Terminate("count=%d != numfilesperdump=%d", count, numfilesperdump);

      std::sort(seq, seq + numfilesperdump);

      files_started   = 0;
      files_completed = 0;

      for(int i = 1; i < std::min<int>(All.MaxFilesWithConcurrentIO, numfilesperdump); i++)
        {
          files_started++;
          MPI_Ssend(&ThisTask, 1, MPI_INT, seq[i].thistask, TAG_N, Communicator);
        }

      files_started++;
      if(append_flag)
        append_file(buf, masterTask, lastTask, CommBuffer, numfilesperdump, chunksize);
      else
        write_file(buf, masterTask, lastTask, CommBuffer, numfilesperdump, chunksize);
      files_completed++;

      if(files_started < numfilesperdump)
        {
          /* send start signal */
          MPI_Ssend(&ThisTask, 1, MPI_INT, seq[files_started++].thistask, TAG_N, Communicator);
        }

      while(files_completed < numfilesperdump)
        polling(numfilesperdump);
    }
  else if(masterTask == ThisTask)
    {
      /* wait for start signal */
      int dummy;
      MPI_Recv(&dummy, 1, MPI_INT, 0, TAG_N, Communicator, MPI_STATUS_IGNORE); /* wait until we are told to start */

      if(append_flag)
        append_file(buf, masterTask, lastTask, CommBuffer, numfilesperdump, chunksize);
      else
        write_file(buf, masterTask, lastTask, CommBuffer, numfilesperdump, chunksize);

      /* send back completion notice */
      MPI_Ssend(&ThisTask, 1, MPI_INT, 0, TAG_KEY, Communicator);
    }
  else
    {
      if(append_flag)
        append_file(buf, masterTask, lastTask, CommBuffer, numfilesperdump, chunksize);
      else
        write_file(buf, masterTask, lastTask, CommBuffer, numfilesperdump, chunksize);
    }

  Mem.myfree(CommBuffer);

  if(ThisTask == 0)
    Mem.myfree(seq);
}

/*! \brief Actually write the snapshot file to the disk
 *
 *  This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the SPH particles.  If
 *  cooling is enabled, mean molecular weight and neutral hydrogen abundance
 *  are written for the gas particles. This is followed by the SPH smoothing
 *  length and further blocks of information, depending on included physics
 *  and compile-time flags.
 *
 *  \param fname string containing the file name
 *  \param writeTask the b of the task in a writing group that which is responsible
 *         for the output operations
 *  \param lastTask the rank of the last task in a writing group
 *
 */
void IO_Def::write_file(char *fname, int writeTask, int lastTask, void *CommBuffer, int numfilesperdump, int chunksize)
{
  int typelist[N_DataGroups];
  long long n_type[N_DataGroups], npart[N_DataGroups], pcsum = 0;
  char label[LABEL_LEN + 1];
  unsigned int blksize, bytes_per_blockelement_in_file = 0;
  FILE *fd        = 0;
  hid_t hdf5_file = 0, hdf5_grp[N_DataGroups], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_dataspace_in_file = 0, hdf5_dataset = 0, hdf5_prop = 0;
  hsize_t dims[2], count[2], start[2];
  int rank             = 0;
  hid_t hdf5_paramsgrp = 0;
  hid_t hdf5_configgrp = 0;

#define SKIP                                 \
  {                                          \
    my_fwrite(&blksize, sizeof(int), 1, fd); \
  }

  fill_file_header(writeTask, lastTask, n_type, npart);

  /* open file and write header */
  if(ThisTask == writeTask)
    {
      if(file_format == FILEFORMAT_HDF5)
        {
          char buf[MAXLEN_PATH_EXTRA];
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s.hdf5", fname);
          mpi_printf("%s file: '%s' (file 1 of %d)\n", info, fname, numfilesperdump);

          rename_file_to_bak_if_it_exists(buf);

          hdf5_file = my_H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

          hdf5_headergrp = my_H5Gcreate(hdf5_file, "/Header", 0);
          write_header_fields(hdf5_headergrp);

          hdf5_paramsgrp = my_H5Gcreate(hdf5_file, "/Parameters", 0);
          write_parameters_attributes_in_hdf5(hdf5_paramsgrp);

          hdf5_configgrp = my_H5Gcreate(hdf5_file, "/Config", 0);
          write_compile_time_options_in_hdf5(hdf5_configgrp);

          for(int type = 0; type < N_DataGroups; type++)
            {
              if(npart[type] > 0)
                {
                  get_datagroup_name(type, buf);
                  hdf5_grp[type] = my_H5Gcreate(hdf5_file, buf, 0);
                }
            }
        }
      else
        {
          rename_file_to_bak_if_it_exists(fname);

          if(!(fd = fopen(fname, "w")))
            Terminate("can't open file `%s' for writing.\n", fname);

          mpi_printf("%s file: '%s' (file 1 of %d)\n", info, fname, numfilesperdump);

          if(file_format == FILEFORMAT_LEGACY2)
            {
              blksize = sizeof(int) + 4 * sizeof(char);
              SKIP;
              my_fwrite((const void *)"HEAD", sizeof(char), 4, fd);
              int nextblock = header_size + 2 * sizeof(int);
              my_fwrite(&nextblock, sizeof(int), 1, fd);
              SKIP;
            }

          blksize = header_size;
          SKIP;
          my_fwrite(header_buf, header_size, 1, fd);
          SKIP;
        }
    }

  for(int blocknr = 0; blocknr < N_IO_Fields; blocknr++)
    {
      polling(numfilesperdump);

      if(IO_Fields[blocknr].type_in_file_output != FILE_NONE)
        {
          unsigned int bytes_per_blockelement = get_bytes_per_memory_blockelement(blocknr, 0);
          int blockmaxlen                     = (int)(COMMBUFFERSIZE / bytes_per_blockelement);
          long long npart_in_block            = get_particles_in_block(blocknr, npart, typelist);
          hid_t hdf5_memory_datatype          = get_hdf5_memorytype_of_block(blocknr);
          char dname[MAXLEN_PATH];
          get_dataset_name(blocknr, dname);

          if(npart_in_block > 0)
            {
              mpi_printf("%s block %d (%s)...\n", info, blocknr, dname);

              if(ThisTask == writeTask)
                {
                  if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
                    {
                      bytes_per_blockelement_in_file =
                          IO_Fields[blocknr].values_per_block * H5Tget_size(get_hdf5_outputtype_of_block(blocknr));

                      if(file_format == FILEFORMAT_LEGACY2)
                        {
                          blksize = sizeof(int) + LABEL_LEN * sizeof(char);
                          SKIP;
                          get_Tab_IO_Label(blocknr, label);
                          my_fwrite(label, sizeof(char), LABEL_LEN, fd);
                          int nextblock = npart_in_block * bytes_per_blockelement_in_file + 2 * sizeof(int);
                          my_fwrite(&nextblock, sizeof(int), 1, fd);
                          SKIP;
                        }

                      blksize = npart_in_block * bytes_per_blockelement_in_file;
                      SKIP;
                    }
                }

              for(int type = 0; type < N_DataGroups; type++)
                {
                  if(typelist[type])
                    {
                      if(ThisTask == writeTask && file_format == FILEFORMAT_HDF5 && npart[type] > 0)
                        {
                          hid_t hdf5_file_datatype = get_hdf5_outputtype_of_block(blocknr);

                          dims[0] = npart[type];
                          dims[1] = get_values_per_blockelement(blocknr);
                          if(dims[1] == 1)
                            rank = 1;
                          else
                            rank = 2;

                          if(chunksize > 0)
                            {
                              hsize_t maxdims[2]     = {H5S_UNLIMITED, dims[1]};
                              hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, maxdims);

                              /* Modify dataset creation properties, i.e. enable chunking  */
                              hdf5_prop             = H5Pcreate(H5P_DATASET_CREATE);
                              hsize_t chunk_dims[2] = {0, dims[1]};
                              chunk_dims[0]         = chunksize;
                              H5Pset_chunk(hdf5_prop, rank, chunk_dims);

                              hdf5_dataset =
                                  my_H5Dcreate(hdf5_grp[type], dname, hdf5_file_datatype, hdf5_dataspace_in_file, hdf5_prop);
                            }
                          else if(IO_Fields[blocknr].compression_on)
                            {
                              hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);

                              /* Modify dataset creation properties, i.e. enable compression  */
                              hdf5_prop             = H5Pcreate(H5P_DATASET_CREATE);
                              hsize_t chunk_dims[2] = {COMPRESSION_CHUNKSIZE, dims[1]};
                              if(chunk_dims[0] > dims[0])
                                chunk_dims[0] = dims[0];
                              H5Pset_chunk(hdf5_prop, rank, chunk_dims); /* set chunk size */
                              H5Pset_shuffle(hdf5_prop);                 /* reshuffle bytes to get better compression ratio */
                              H5Pset_deflate(hdf5_prop, 9);              /* gzip compression level 9 */
                              if(H5Pall_filters_avail(hdf5_prop))
                                hdf5_dataset =
                                    my_H5Dcreate(hdf5_grp[type], dname, hdf5_file_datatype, hdf5_dataspace_in_file, hdf5_prop);
                              else
                                Terminate("HDF5: Compression not available!\n");
                            }
                          else
                            {
                              hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);
                              hdf5_dataset =
                                  my_H5Dcreate(hdf5_grp[type], dname, hdf5_file_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
                            }

                          write_dataset_attributes(hdf5_dataset, blocknr);

                          byte_count += dims[0] * dims[1] * my_H5Tget_size(hdf5_file_datatype); /* for I/O performance measurement */

                          pcsum = 0;
                        }

                      for(int task = writeTask, offset = 0; task <= lastTask; task++)
                        {
                          long long n_for_this_task;

                          if(task == ThisTask)
                            {
                              n_for_this_task = n_type[type];

                              for(int p = writeTask; p <= lastTask; p++)
                                if(p != ThisTask)
                                  MPI_Send(&n_for_this_task, sizeof(n_for_this_task), MPI_BYTE, p, TAG_NFORTHISTASK, Communicator);
                            }
                          else
                            MPI_Recv(&n_for_this_task, sizeof(n_for_this_task), MPI_BYTE, task, TAG_NFORTHISTASK, Communicator,
                                     MPI_STATUS_IGNORE);

                          while(n_for_this_task > 0)
                            {
                              long long pc = n_for_this_task;

                              if(pc > blockmaxlen)
                                pc = blockmaxlen;

                              if(ThisTask == task)
                                fill_write_buffer(blocknr, &offset, pc, type, CommBuffer);

                              if(ThisTask == writeTask && task != writeTask)
                                MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, Communicator,
                                         MPI_STATUS_IGNORE);

                              if(ThisTask != writeTask && task == ThisTask)
                                MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask, TAG_PDATA, Communicator);

                              if(ThisTask == writeTask)
                                {
                                  if(file_format == FILEFORMAT_HDF5)
                                    {
                                      start[0] = pcsum;
                                      start[1] = 0;

                                      count[0] = pc;
                                      count[1] = get_values_per_blockelement(blocknr);
                                      pcsum += pc;

                                      my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                      dims[0]               = pc;
                                      dims[1]               = get_values_per_blockelement(blocknr);
                                      hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                                      my_H5Dwrite(hdf5_dataset, hdf5_memory_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                                  H5P_DEFAULT, CommBuffer, dname);

                                      my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
                                    }
                                  else
                                    {
                                      if(bytes_per_blockelement_in_file != bytes_per_blockelement)
                                        {
                                          char *CommAuxBuffer =
                                              (char *)Mem.mymalloc("CommAuxBuffer", bytes_per_blockelement_in_file * pc);
                                          type_cast_data((char *)CommBuffer, bytes_per_blockelement, (char *)CommAuxBuffer,
                                                         bytes_per_blockelement_in_file, pc, blocknr);
                                          my_fwrite(CommAuxBuffer, bytes_per_blockelement_in_file, pc, fd);
                                          Mem.myfree(CommAuxBuffer);
                                        }
                                      else
                                        my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
                                    }
                                }

                              n_for_this_task -= pc;
                            }
                        }

                      if(ThisTask == writeTask && file_format == FILEFORMAT_HDF5 && npart[type] > 0)
                        {
                          if(file_format == FILEFORMAT_HDF5)
                            {
                              my_H5Dclose(hdf5_dataset, dname);
                              if(chunksize > 0 || IO_Fields[blocknr].compression_on)
                                my_H5Pclose(hdf5_prop);
                              my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                            }
                        }
                    }
                }

              if(ThisTask == writeTask)
                {
                  if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
                    SKIP;
                }
            }
        }
    }

  if(ThisTask == writeTask)
    {
      if(file_format == FILEFORMAT_HDF5)
        {
          char buf[MAXLEN_PATH];

          for(int type = N_DataGroups - 1; type >= 0; type--)
            if(npart[type] > 0)
              {
                get_datagroup_name(type, buf);
                my_H5Gclose(hdf5_grp[type], buf);
              }

          my_H5Gclose(hdf5_headergrp, "/Header");
          my_H5Gclose(hdf5_paramsgrp, "/Parameters");
          my_H5Gclose(hdf5_configgrp, "/Config");

          snprintf(buf, MAXLEN_PATH, "%s.hdf5", fname);
          my_H5Fclose(hdf5_file, buf);
        }
      else
        fclose(fd);
    }
}

void IO_Def::append_file(char *fname, int writeTask, int lastTask, void *CommBuffer, int numfilesperdump, int chunksize)
{
  int typelist[N_DataGroups];
  long long n_type[N_DataGroups], npart[N_DataGroups], n_previous[N_DataGroups], pcsum = 0;
  hid_t hdf5_file = 0, hdf5_grp[N_DataGroups], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  hsize_t dims[2], count[2], start[2];
  int rank = 0;

  if(file_format != FILEFORMAT_HDF5)
    Terminate("appending to files only works with HDF5 format\n");

  read_file_header(NULL, 0, writeTask, lastTask, n_type, npart, NULL);

  for(int n = 0; n < N_DataGroups; n++)
    n_previous[n] = n_type[n];

  fill_file_header(writeTask, lastTask, n_type, npart);

  /* open file and write header */
  if(ThisTask == writeTask)
    {
      char buf[MAXLEN_PATH_EXTRA];
      snprintf(buf, MAXLEN_PATH_EXTRA, "%s.hdf5", fname);

      hdf5_file = my_H5Fopen(buf, H5F_ACC_RDWR, H5P_DEFAULT);

      hdf5_headergrp = my_H5Gopen(hdf5_file, "/Header");
      write_header_fields(hdf5_headergrp);

      for(int type = 0; type < N_DataGroups; type++)
        {
          if(npart[type] > 0)
            {
              get_datagroup_name(type, buf);
              if(n_previous[type] == 0)
                hdf5_grp[type] = my_H5Gcreate(hdf5_file, buf, 0);
              else
                hdf5_grp[type] = my_H5Gopen(hdf5_file, buf);
            }
        }
    }

  for(int blocknr = 0; blocknr < N_IO_Fields; blocknr++)
    {
      polling(numfilesperdump);

      if(IO_Fields[blocknr].type_in_file_output != FILE_NONE)
        {
          unsigned int bytes_per_blockelement = get_bytes_per_memory_blockelement(blocknr, 0);
          int blockmaxlen                     = (int)(COMMBUFFERSIZE / bytes_per_blockelement);
          long long npart_in_block            = get_particles_in_block(blocknr, npart, typelist);
          hid_t hdf5_memory_datatype          = get_hdf5_memorytype_of_block(blocknr);
          char dname[MAXLEN_PATH];
          get_dataset_name(blocknr, dname);

          if(npart_in_block > 0)
            {
              mpi_printf("%s block %d (%s)...\n", info, blocknr, dname);

              for(int type = 0; type < N_DataGroups; type++)
                {
                  if(typelist[type])
                    {
                      if(ThisTask == writeTask && file_format == FILEFORMAT_HDF5 && npart[type] > 0)
                        {
                          hid_t hdf5_file_datatype = get_hdf5_outputtype_of_block(blocknr);

                          dims[0] = npart[type] + n_previous[type];
                          dims[1] = get_values_per_blockelement(blocknr);
                          if(dims[1] == 1)
                            rank = 1;
                          else
                            rank = 2;

                          hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);

                          if(n_previous[type] == 0)
                            {
                              if(chunksize > 0)
                                {
                                  hsize_t maxdims[2]     = {H5S_UNLIMITED, dims[1]};
                                  hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, maxdims);

                                  /* Modify dataset creation properties, i.e. enable chunking  */
                                  hid_t prop            = H5Pcreate(H5P_DATASET_CREATE);
                                  hsize_t chunk_dims[2] = {0, dims[1]};
                                  chunk_dims[0]         = chunksize;
                                  H5Pset_chunk(prop, rank, chunk_dims);

                                  hdf5_dataset = my_H5Dcreate(hdf5_grp[type], dname, hdf5_file_datatype, hdf5_dataspace_in_file, prop);
                                }
                              else
                                {
                                  hdf5_dataset =
                                      my_H5Dcreate(hdf5_grp[type], dname, hdf5_file_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
                                  write_dataset_attributes(hdf5_dataset, blocknr);
                                }
                            }
                          else
                            {
                              hdf5_dataset = my_H5Dopen_if_existing(hdf5_grp[type], dname);
                              my_H5Dset_extent(hdf5_dataset, dims);
                            }

                          byte_count += dims[0] * dims[1] * my_H5Tget_size(hdf5_file_datatype); /* for I/O performance measurement */

                          pcsum = 0;
                        }

                      for(int task = writeTask, offset = 0; task <= lastTask; task++)
                        {
                          long long n_for_this_task;

                          if(task == ThisTask)
                            {
                              n_for_this_task = n_type[type];

                              for(int p = writeTask; p <= lastTask; p++)
                                if(p != ThisTask)
                                  MPI_Send(&n_for_this_task, sizeof(n_for_this_task), MPI_BYTE, p, TAG_NFORTHISTASK, Communicator);
                            }
                          else
                            MPI_Recv(&n_for_this_task, sizeof(n_for_this_task), MPI_BYTE, task, TAG_NFORTHISTASK, Communicator,
                                     MPI_STATUS_IGNORE);

                          while(n_for_this_task > 0)
                            {
                              long long pc = n_for_this_task;

                              if(pc > blockmaxlen)
                                pc = blockmaxlen;

                              if(ThisTask == task)
                                fill_write_buffer(blocknr, &offset, pc, type, CommBuffer);

                              if(ThisTask == writeTask && task != writeTask)
                                MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, Communicator,
                                         MPI_STATUS_IGNORE);

                              if(ThisTask != writeTask && task == ThisTask)
                                MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask, TAG_PDATA, Communicator);

                              if(ThisTask == writeTask)
                                {
                                  start[0] = pcsum + n_previous[type];
                                  start[1] = 0;

                                  count[0] = pc;
                                  count[1] = get_values_per_blockelement(blocknr);
                                  pcsum += pc;

                                  my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                  dims[0]               = pc;
                                  dims[1]               = get_values_per_blockelement(blocknr);
                                  hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                                  my_H5Dwrite(hdf5_dataset, hdf5_memory_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                              H5P_DEFAULT, CommBuffer, dname);

                                  my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
                                }

                              n_for_this_task -= pc;
                            }
                        }

                      if(ThisTask == writeTask && npart[type] > 0)
                        {
                          my_H5Dclose(hdf5_dataset, dname);
                          my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                        }
                    }
                }
            }
        }
    }

  if(ThisTask == writeTask)
    {
      char buf[MAXLEN_PATH];

      for(int type = N_DataGroups - 1; type >= 0; type--)
        if(npart[type] > 0)
          {
            get_datagroup_name(type, buf);
            my_H5Gclose(hdf5_grp[type], buf);
          }

      my_H5Gclose(hdf5_headergrp, "/Header");

      snprintf(buf, MAXLEN_PATH, "%s.hdf5", fname);
      my_H5Fclose(hdf5_file, buf);
    }
}

void IO_Def::type_cast_data(char *src, int src_bytes_per_element, char *target, int target_bytes_per_element, int len, int blocknr)
{
  switch(IO_Fields[blocknr].type_in_memory)
    {
      case MEM_INT:
      case MEM_INT64:
      case MEM_MY_ID_TYPE:
      case MEM_MY_INTPOS_TYPE:
      case MEM_MY_FILEOFFSET:
        if(target_bytes_per_element > src_bytes_per_element)
          {
            if(target_bytes_per_element != 2 * src_bytes_per_element)
              Terminate("something is odd here: target_bytes_per_element=%d != 2 * src_bytes_per_element=%d", target_bytes_per_element,
                        src_bytes_per_element);

            int fac       = src_bytes_per_element / sizeof(int);
            int *sp       = (int *)src;
            long long *tp = (long long *)target;
            for(int i = 0; i < len * fac; i++)
              *tp++ = *sp++;
          }
        else
          {
            if(src_bytes_per_element != 2 * target_bytes_per_element)
              Terminate("something is odd here: src_bytes_per_element=%d != 2 * target_bytes_per_element=%d", src_bytes_per_element,
                        target_bytes_per_element);
            int fac       = src_bytes_per_element / sizeof(long long);
            long long *sp = (long long *)src;
            int *tp       = (int *)target;
            for(int i = 0; i < len * fac; i++)
              *tp++ = *sp++;
          }
        break;

      case MEM_FLOAT:
      case MEM_DOUBLE:
      case MEM_MY_FLOAT:
      case MEM_MY_DOUBLE:
        if((target_bytes_per_element % 8) == 0 &&
           (src_bytes_per_element % 4) == 0)  // target_bytes_per_element multiply of 8, src_bytes_per_element mulitple of 4
          {
            int fac    = src_bytes_per_element / sizeof(float);
            float *sp  = (float *)src;
            double *tp = (double *)target;
            for(int i = 0; i < len * fac; i++)
              *tp++ = *sp++;
          }
        else if((target_bytes_per_element % 4) == 0 && (src_bytes_per_element % 8) == 0)
          {
            int fac    = src_bytes_per_element / sizeof(double);
            double *sp = (double *)src;
            float *tp  = (float *)target;
            for(int i = 0; i < len * fac; i++)
              *tp++ = *sp++;
          }
        else if((target_bytes_per_element & 8) == 0 && (src_bytes_per_element % 2) == 0)
          {
            int fac    = src_bytes_per_element / sizeof(half);
            half *sp   = (half *)src;
            double *tp = (double *)target;
            for(int i = 0; i < len * fac; i++)
              *tp++ = *sp++;
          }
        else if((target_bytes_per_element % 2) == 0 && (src_bytes_per_element % 8) == 0)
          {
            int fac    = src_bytes_per_element / sizeof(double);
            double *sp = (double *)src;
            half *tp   = (half *)target;
            for(int i = 0; i < len * fac; i++)
              *tp++ = *sp++;
          }
        else
          {
            Terminate("Strange conversion requested: target_bytes_per_element=%d  src_bytes_per_element=%d\n",
                      target_bytes_per_element, src_bytes_per_element);
          }
        break;
    }
}

/*! \brief This function reads a file
 *
 *  This routine reads a single file. The data it contains is
 *  distributed to tasks 'readTask' to 'lastTask'.
 *
 *  \param fname filename to be read
 *  \param readTask task responsible for reading the file fname
 *  \param lastTask last Task which gets data contained in the file
 */
void IO_Def::read_file(const char *fname, int filenr, int readTask, int lastTask, void *CommBuffer)
{
  long long n_type[N_DataGroups], npart[N_DataGroups];
  int typelist[N_DataGroups];
  unsigned int blksize1, blksize2, bytes_per_blockelement_in_file = 0;
  hid_t hdf5_file = 0, hdf5_grp[N_DataGroups], hdf5_dataspace_in_file;
  hid_t hdf5_dataspace_in_memory, hdf5_dataset;
  FILE *fd = 0;

  /* open file and read header */

  if(ThisTask == readTask)
    {
      if(file_format == FILEFORMAT_HDF5)
        {
          read_header_fields(fname);
        }
      else if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
        {
          if(!(fd = fopen(fname, "r")))
            Terminate("can't open file `%s' for reading initial conditions.\n", fname);

          if(file_format == FILEFORMAT_LEGACY2)
            {
              int nextblock;
              char label[LABEL_LEN];
              my_fread(&blksize1, sizeof(int), 1, fd);
              my_fread(&label, sizeof(char), LABEL_LEN, fd);
              my_fread(&nextblock, sizeof(int), 1, fd);
              my_fread(&blksize2, sizeof(int), 1, fd);
            }

          my_fread(&blksize1, sizeof(int), 1, fd);
          my_fread(header_buf, header_size, 1, fd);
          my_fread(&blksize2, sizeof(int), 1, fd);
        }
      else
        Terminate("unknown ICFormat");

      for(int task = readTask + 1; task <= lastTask; task++)
        MPI_Ssend(header_buf, header_size, MPI_BYTE, task, TAG_HEADER, Communicator);
    }
  else
    MPI_Recv(header_buf, header_size, MPI_BYTE, readTask, TAG_HEADER, Communicator, MPI_STATUS_IGNORE);

  int nstart;
  read_file_header(fname, filenr, readTask, lastTask, n_type, npart, &nstart);

  if(ThisTask == readTask)
    {
      if(file_format == FILEFORMAT_HDF5)
        {
          hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

          for(int type = 0; type < N_DataGroups; type++)
            {
              if(npart[type] > 0)
                {
                  char buf[MAXLEN_PATH];
                  get_datagroup_name(type, buf);
                  hdf5_grp[type] = my_H5Gopen(hdf5_file, buf);
                }
            }
        }
    }

  for(int blocknr = 0; blocknr < N_IO_Fields; blocknr++)
    {
      if((IO_Fields[blocknr].read_flag != SKIP_ON_READ &&
          !(file_format == FILEFORMAT_LEGACY1 && All.RestartFlag == RST_BEGIN && type_of_file == FILE_IS_SNAPSHOT &&
            blocknr > 4) /* this second conditions allows short legacy ICs to be read in */
          ) ||
         IO_Fields[blocknr].type_in_memory == MEM_MY_FILEOFFSET)
        {
          unsigned int bytes_per_blockelement = get_bytes_per_memory_blockelement(blocknr, 1);
          int blockmaxlen                     = (int)(COMMBUFFERSIZE / bytes_per_blockelement);
          long long npart_in_block            = get_particles_in_block(blocknr, npart, &typelist[0]);
          hid_t hdf5_memory_datatype          = get_hdf5_memorytype_of_block(blocknr);
          char dname[MAXLEN_PATH];
          get_dataset_name(blocknr, dname);

          if(npart_in_block > 0)
            {
              if(filenr == 0)
                mpi_printf("READIC: reading block %d (%s)...\n", blocknr, dname);

              if(ThisTask == readTask && IO_Fields[blocknr].type_in_memory != MEM_MY_FILEOFFSET)
                {
                  if(file_format == FILEFORMAT_LEGACY2)
                    {
                      char expected_label[LABEL_LEN + 1];
                      get_Tab_IO_Label(blocknr, expected_label);
                      find_block(expected_label, fd);
                    }

                  if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
                    {
                      my_fread(&blksize1, sizeof(int), 1, fd);
                      bytes_per_blockelement_in_file = blksize1 / npart_in_block;
                    }
                }

              int offset = 0;

              for(int type = 0; type < N_DataGroups; type++)
                {
                  if(type_of_file != FILE_IS_SNAPSHOT)
                    offset = 0;

                  int n_in_file = npart[type];
                  int pcsum     = 0;

                  long long nprevious = 0;
                  for(int t = 0; t < type; t++)
                    for(int f = 0; f < get_filenr_from_header(); f++)
                      nprevious += ntype_in_files[f * N_DataGroups + t];

                  for(int nr = 0; nr < filenr; nr++)
                    nprevious += ntype_in_files[nr * N_DataGroups + type];

                  if(typelist[type] == 0)
                    {
                      int ntask           = lastTask - readTask + 1;
                      int n_for_this_task = n_in_file / ntask;
                      if((ThisTask - readTask) < (n_in_file % ntask))
                        n_for_this_task++;

                      offset += n_for_this_task;
                    }
                  else
                    {
                      for(int task = readTask; task <= lastTask; task++)
                        {
                          int ntask           = lastTask - readTask + 1;
                          int n_for_this_task = n_in_file / ntask;
                          if((task - readTask) < (n_in_file % ntask))
                            n_for_this_task++;

                          do
                            {
                              int pc = n_for_this_task;

                              if(pc > blockmaxlen)
                                pc = blockmaxlen;

                              if(ThisTask == readTask && IO_Fields[blocknr].type_in_memory != MEM_MY_FILEOFFSET)
                                {
                                  if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
                                    {
                                      if(bytes_per_blockelement_in_file != bytes_per_blockelement)
                                        {
                                          char *CommAuxBuffer =
                                              (char *)Mem.mymalloc("CommAuxBuffer", bytes_per_blockelement_in_file * pc);
                                          my_fread(CommAuxBuffer, bytes_per_blockelement_in_file, pc, fd);
                                          type_cast_data((char *)CommAuxBuffer, bytes_per_blockelement_in_file, (char *)CommBuffer,
                                                         bytes_per_blockelement, pc, blocknr);
                                          Mem.myfree(CommAuxBuffer);
                                        }
                                      else
                                        my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
                                    }

                                  if(file_format == FILEFORMAT_HDF5 && pc > 0)
                                    {
                                      hdf5_dataset = my_H5Dopen_if_existing(hdf5_grp[type], dname);
                                      int rank;
                                      hsize_t dims[2], count[2], start[2];

                                      dims[0] = npart[type];
                                      dims[1] = get_values_per_blockelement(blocknr);
                                      if(dims[1] == 1)
                                        rank = 1;
                                      else
                                        rank = 2;

                                      hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);

                                      dims[0]                  = pc;
                                      hdf5_dataspace_in_memory = my_H5Screate_simple(rank, dims, NULL);

                                      start[0] = pcsum;
                                      start[1] = 0;

                                      count[0] = pc;
                                      count[1] = get_values_per_blockelement(blocknr);
                                      pcsum += pc;

                                      my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                      // Test if dataset was present
                                      if(hdf5_dataset < 0)
                                        {
                                          // no, pad with zeros
                                          if((ThisTask == readTask) && (task == ThisTask))
                                            printf("\tDataset %s not present for particle type %d, using zero.\n", dname, type);
                                          memset(CommBuffer, 0, dims[0] * dims[1] * my_H5Tget_size(hdf5_memory_datatype));
                                        }
                                      else
                                        {
                                          hid_t hdf5_file_datatype = H5Dget_type(hdf5_dataset);
                                          byte_count += dims[0] * dims[1] *
                                                        my_H5Tget_size(hdf5_file_datatype); /* for I/O performance measurement */
                                          H5Tclose(hdf5_file_datatype);

                                          my_H5Dread(hdf5_dataset, hdf5_memory_datatype, hdf5_dataspace_in_memory,
                                                     hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer, dname);
                                          my_H5Dclose(hdf5_dataset, dname);
                                        }
                                      my_H5Sclose(hdf5_dataspace_in_memory, H5S_SIMPLE);
                                      my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                                    }
                                }

                              if(ThisTask == readTask && task != readTask && pc > 0)
                                MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, Communicator);

                              if(ThisTask != readTask && task == ThisTask && pc > 0)
                                MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask, TAG_PDATA, Communicator,
                                         MPI_STATUS_IGNORE);

                              if(ThisTask == task)
                                {
                                  if(blocknr == 0 && IO_Fields[blocknr].array == A_P)
                                    {
                                      for(int n = 0; n < pc; n++)
                                        set_type_of_element(nstart + offset + n, type); /* initialize type */
                                    }
#ifdef MERGERTREE
                                  if(blocknr == 0 && IO_Fields[blocknr].array == A_MTRP)
                                    {
                                      for(int n = 0; n < pc; n++)
                                        {
                                          set_type_of_element(nstart + offset + n, type); /* initialize type */
                                          //         MtrP[nstart + offset + n].Type = type;  /* initialize type */
                                        }
                                    }
#endif
                                  empty_read_buffer(blocknr, nstart + offset, pc, type, nprevious, CommBuffer);

                                  offset += pc;
                                }

                              n_for_this_task -= pc;
                              nprevious += pc;
                            }
                          while(n_for_this_task > 0);
                        }
                    }
                }
              if(ThisTask == readTask && IO_Fields[blocknr].type_in_memory != MEM_MY_FILEOFFSET)
                {
                  if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
                    {
                      my_fread(&blksize2, sizeof(int), 1, fd);
                      ;
                      if(blksize1 != blksize2)
                        {
                          char buf[MAXLEN_PATH_EXTRA];
                          snprintf(buf, MAXLEN_PATH,
                                   "incorrect block-sizes detected!\n Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask,
                                   blocknr, blksize1, blksize2);
                          if(blocknr == 2) /* block number 2 is always IDs */
                            strcat(buf, "Possible mismatch of 32bit and 64bit ID's in IC file and GADGET compilation !\n");
                          Terminate(buf);
                        }
                    }
                }
            }
        }
    }

  for(int type = 0; type < N_DataGroups; type++)
    {
      long long n_in_file = npart[type];
      int ntask           = lastTask - readTask + 1;
      int n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
        n_for_this_task++;

      read_increase_numbers(type, n_for_this_task);
    }

  if(ThisTask == readTask)
    {
      if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
        fclose(fd);

      if(file_format == FILEFORMAT_HDF5)
        {
          for(int type = N_DataGroups - 1; type >= 0; type--)
            if(npart[type] > 0)
              {
                char buf[MAXLEN_PATH];
                get_datagroup_name(type, buf);
                my_H5Gclose(hdf5_grp[type], buf);
              }
          my_H5Fclose(hdf5_file, fname);
        }
    }
}

/*! \brief This function assigns a certain number of tasks to each file.
 *
 *  These tasks are containing the content of that file after the ICs have been read
 *  The number of tasks per file is as homogeneous as possible.
 *  The number of files may at most be equal to the number of tasks.
 *
 *  \param nfiles Number of files of which the snapshot is distributed
 *  \param filenr contains the file number to which this task belongs
 *  \param master the number of the task responsible to read the file
 *  \param last number of the last task belonging to the same file as this task
 */
void IO_Def::distribute_file(int nfiles, int *filenr, int *master, int *last)
{
  int tasks_per_file = NTask / nfiles;
  int tasks_left     = NTask % nfiles;

  if(tasks_left == 0)
    {
      int group = ThisTask / tasks_per_file;
      *master   = group * tasks_per_file;
      *last     = (group + 1) * tasks_per_file - 1;
      *filenr   = group;
      return;
    }

  double tpf = ((double)NTask) / nfiles;

  *last = -1;
  for(int i = 0; i < nfiles; i++)
    {
      *master = *last + 1;
      *last   = (i + 1) * tpf;
      if(*last >= NTask)
        *last = *last - 1;
      if(*last < *master)
        Terminate("last < master");
      *filenr = i;

      if(i == nfiles - 1)
        *last = NTask - 1;

      if(ThisTask >= *master && ThisTask <= *last)
        return;
    }
}

/*! \brief This function tells the size in bytes of one data entry in each of the blocks
 *  defined for the output file.
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param mode used to distinguish whether the function is called in input
 *         mode (mode > 0) or in output mode (mode = 0). The size of one data
 *         entry may vary depending on the mode
 *  \return size of the data entry in bytes
 */
int IO_Def::get_bytes_per_memory_blockelement(int blocknr, int mode)
{
  if(blocknr < 0 || blocknr >= N_IO_Fields)
    Terminate("something is wrong here: blocknr=%d N_IO_Fields=%d", blocknr, N_IO_Fields);

  int bytes_per_blockelement = 0;

  IO_Field *field = &IO_Fields[blocknr];

  switch(field->type_in_memory)
    {
      case MEM_INT:
        bytes_per_blockelement = field->values_per_block * sizeof(int);
        break;
      case MEM_INT64:
        bytes_per_blockelement = field->values_per_block * sizeof(long long);
        break;
      case MEM_MY_ID_TYPE:
        bytes_per_blockelement = field->values_per_block * sizeof(MyIDType);
        break;
      case MEM_MY_INTPOS_TYPE:
        bytes_per_blockelement = field->values_per_block * sizeof(MyIntPosType);
        break;
      case MEM_FLOAT:
        bytes_per_blockelement = field->values_per_block * sizeof(float);
        break;
      case MEM_DOUBLE:
        bytes_per_blockelement = field->values_per_block * sizeof(double);
        break;
      case MEM_MY_FLOAT:
        bytes_per_blockelement = field->values_per_block * sizeof(MyFloat);
        break;
      case MEM_MY_DOUBLE:
        bytes_per_blockelement = field->values_per_block * sizeof(MyDouble);
        break;
      case MEM_MY_FILEOFFSET:
        bytes_per_blockelement = field->values_per_block * sizeof(long long);
        break;
    }

  return bytes_per_blockelement;
}

hid_t IO_Def::get_hdf5_outputtype_of_block(int blocknr)
{
  hid_t hdf5_datatype = 0;

  switch(IO_Fields[blocknr].type_in_file_output)
    {
      case FILE_INT:
        hdf5_datatype = H5T_NATIVE_INT;
        break;
      case FILE_INT64:
        hdf5_datatype = H5T_NATIVE_INT64;
        break;
      case FILE_MY_IO_FLOAT:
#ifdef OUTPUT_IN_DOUBLEPRECISION
        hdf5_datatype = H5T_NATIVE_DOUBLE;
#else
        hdf5_datatype = H5T_NATIVE_FLOAT;
#endif
        break;
      case FILE_MY_ID_TYPE:
#if defined(IDS_32BIT)
        hdf5_datatype = H5T_NATIVE_UINT32;
#elif defined(IDS_48BIT)
        hdf5_datatype = Int48_memtype;
#else
        hdf5_datatype = H5T_NATIVE_UINT64;
#endif
        break;
      case FILE_MY_INTPOS_TYPE:
#if defined(POSITIONS_IN_32BIT)
        hdf5_datatype = H5T_NATIVE_UINT32;
#elif defined(POSITIONS_IN_64BIT)
        hdf5_datatype = H5T_NATIVE_UINT64;
#else
        hdf5_datatype = Int128_memtype;
#endif
        break;
      case FILE_DOUBLE:
        hdf5_datatype = H5T_NATIVE_DOUBLE;
        break;
      case FILE_FLOAT:
        hdf5_datatype = H5T_NATIVE_FLOAT;
        break;
      case FILE_HALF:
        hdf5_datatype = Halfprec_memtype;
        break;
      case FILE_NONE:
        Terminate("undefined type");
        break;
    }

  return hdf5_datatype;
}

hid_t IO_Def::get_hdf5_memorytype_of_block(int blocknr)
{
  hid_t hdf5_datatype = 0;

  switch(IO_Fields[blocknr].type_in_memory)
    {
      case MEM_INT:
        hdf5_datatype = H5T_NATIVE_INT;
        break;
      case MEM_INT64:
        hdf5_datatype = H5T_NATIVE_INT64;
        break;
      case MEM_MY_ID_TYPE:
#ifdef IDS_32BIT
        hdf5_datatype = H5T_NATIVE_UINT32;
#else
        hdf5_datatype = H5T_NATIVE_UINT64;
#endif
        break;
      case MEM_MY_INTPOS_TYPE:
#ifdef POSITIONS_IN_32BIT
        hdf5_datatype = H5T_NATIVE_UINT32;
#elif defined(POSITIONS_IN_64BIT)
        hdf5_datatype = H5T_NATIVE_UINT64;
#else
        hdf5_datatype = Int128_memtype;
#endif
        break;
      case MEM_FLOAT:
        hdf5_datatype = H5T_NATIVE_FLOAT;
        break;
      case MEM_DOUBLE:
        hdf5_datatype = H5T_NATIVE_DOUBLE;
        break;
      case MEM_MY_FLOAT:
        hdf5_datatype = H5T_NATIVE_MYFLOAT;
        break;
      case MEM_MY_DOUBLE:
        hdf5_datatype = H5T_NATIVE_MYDOUBLE;
        break;
      case MEM_MY_FILEOFFSET:
        hdf5_datatype = H5T_NATIVE_INT64;
        break;
    }

  return hdf5_datatype;
}

/*! \brief This function determines the number of elements composing one data entry
 *  in each of the blocks defined for the output file.
 *
 *  Used only if output in HDF5 format is enabled
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \return number of elements of one data entry
 */
int IO_Def::get_values_per_blockelement(int blocknr)
{
  if(blocknr < 0 || blocknr >= N_IO_Fields)
    Terminate("something is wrong here: blocknr=%d N_IO_Fields=%d", blocknr, N_IO_Fields);

  return IO_Fields[blocknr].values_per_block;
}

/*! \brief Get particle number in an output block
 *
 *  This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param typelist array that contains the number of particles of each type in the block
 *  \return the total number of particles in the block
 */
long long IO_Def::get_particles_in_block(int blocknr, long long *npart_file, int *typelist)
{
  long long npart = 0;

  for(int i = 0; i < N_DataGroups; i++)
    typelist[i] = 0;

  switch(IO_Fields[blocknr].typelist)
    {
      case MASS_BLOCK:
        for(int i = 0; i < NTYPES; i++)
          {
            typelist[i] = (All.MassTable[i] == 0 && npart_file[i] > 0);
            npart += npart_file[i] * typelist[i];
          }
        return npart;
        break;

      case AGE_BLOCK:
        for(int i = 0; i < NTYPES; i++)
          {
            typelist[i] = (npart_file[i] > 0);

            if((file_format == FILEFORMAT_HDF5 && (i == 0 || i == 1 || i == 5)) || (file_format != FILEFORMAT_HDF5 && i != 4))
              typelist[i] = 0;

            npart += npart_file[i] * typelist[i];
          }
        return npart;
        break;

      case Z_BLOCK:
        for(int i = 0; i < NTYPES; i++)
          {
            typelist[i] = (npart_file[i] > 0);
            if((file_format == FILEFORMAT_HDF5 && (i == 1 || i == 5)) || (file_format != FILEFORMAT_HDF5 && (i != 0 && i != 4)))
              typelist[i] = 0;

            npart += npart_file[i] * typelist[i];
          }
        return npart;
        break;

      case GROUPS:
        npart       = npart_file[0];
        typelist[0] = 1;
        return npart;
        break;

      case SUBGROUPS:
        npart       = npart_file[1];
        typelist[1] = 1;
        return npart;
        break;

      case ID_BLOCK:
        npart       = npart_file[2];
        typelist[2] = 1;
        return npart;
        break;

      case CURRSUBS:
      case PREVSUBS:
      case TREELINK:
      case TREELENGTH:
      case MASSMAPS:
      case GALSNAPS:
        npart       = npart_file[0];
        typelist[0] = 1;
        return npart;
        break;

      case TREEHALOS:
        npart       = npart_file[1];
        typelist[1] = 1;
        return npart;
        break;

      case TREETIMES:
        npart       = npart_file[2];
        typelist[2] = 1;
        return npart;
        break;

      case TREETABLE:
        npart            = npart_file[NTYPES];
        typelist[NTYPES] = 1;
        return npart;
        break;

      case HEALPIXTAB:
        npart                = npart_file[NTYPES + 1];
        typelist[NTYPES + 1] = 1;
        return npart;
        break;

      default:
        for(int i = 0; i < N_DataGroups; i++)
          {
            if((IO_Fields[blocknr].typelist & (1 << i)) && npart_file[i] > 0)
              {
                typelist[i] = 1;
                npart += npart_file[i];
              }
            else
              typelist[i] = 0;
          }
        return npart;
        break;
    }

  return 0;
}

/*! \brief This function associates a short 4-character block name with each block number.
 *
 *   This is stored in front of each block for snapshot FileFormat=2.
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param label string containing the dataset name
 */
void IO_Def::get_Tab_IO_Label(int blocknr, char *label)
{
  if(blocknr < 0 || blocknr >= N_IO_Fields)
    Terminate("something is wrong here: blocknr=%d N_IO_Fields=%d", blocknr, N_IO_Fields);

  strcpy(label, IO_Fields[blocknr].label);
}

/*! \brief This function associates a dataset name with each block number.
 *
 *   This is needed to name the dataset if the output is written in HDF5 format
 *
 *  \param blocknr ID of the output block (i.e. position, velocities...)
 *  \param buf string containing the dataset name
 */
void IO_Def::get_dataset_name(int blocknr, char *buf)
{
  if(blocknr < 0 || blocknr >= N_IO_Fields)
    Terminate("something is wrong here: blocknr=%d N_IO_Fields=%d", blocknr, N_IO_Fields);

  strcpy(buf, IO_Fields[blocknr].datasetname);
}

void IO_Def::write_dataset_attributes(hid_t hdf5_dataset, int blocknr)
{
  if(blocknr < 0 || blocknr >= N_IO_Fields)
    Terminate("something is wrong here: blocknr=%d N_IO_Fields=%d", blocknr, N_IO_Fields);

  if(IO_Fields[blocknr].hasunit == 0)
    return;

  if(All.ComovingIntegrationOn)
    {
      write_scalar_attribute(hdf5_dataset, "a_scaling", &IO_Fields[blocknr].a, H5T_NATIVE_DOUBLE);
      write_scalar_attribute(hdf5_dataset, "h_scaling", &IO_Fields[blocknr].h, H5T_NATIVE_DOUBLE);
    }
  else
    {
      double zero = 0;
      write_scalar_attribute(hdf5_dataset, "a_scaling", &zero, H5T_NATIVE_DOUBLE);
      write_scalar_attribute(hdf5_dataset, "h_scaling", &zero, H5T_NATIVE_DOUBLE);
    }

  write_scalar_attribute(hdf5_dataset, "length_scaling", &IO_Fields[blocknr].L, H5T_NATIVE_DOUBLE);
  write_scalar_attribute(hdf5_dataset, "mass_scaling", &IO_Fields[blocknr].M, H5T_NATIVE_DOUBLE);
  write_scalar_attribute(hdf5_dataset, "velocity_scaling", &IO_Fields[blocknr].V, H5T_NATIVE_DOUBLE);

  write_scalar_attribute(hdf5_dataset, "to_cgs", &IO_Fields[blocknr].c, H5T_NATIVE_DOUBLE);
}

/*! \brief Write the parameters read from the parameter file into the HDF5 snapshot file
 *
 *  This function stores the parameter io_header as attributes belonging
 *  to the parameter group of the HDF5 file.
 *
 *  \param handle contains a reference to the parameter group
 */
void IO_Def::write_parameters_attributes_in_hdf5(hid_t handle)
{
  for(int i = 0; i < All.NParameters; i++)
    {
      switch(All.ParametersType[i])
        {
          case PARAM_DOUBLE:
            write_scalar_attribute(handle, All.ParametersTag[i], All.ParametersValue[i], H5T_NATIVE_DOUBLE);
            break;
          case PARAM_STRING:
            write_string_attribute(handle, All.ParametersTag[i], (const char *)All.ParametersValue[i]);
            break;
          case PARAM_INT:
            write_scalar_attribute(handle, All.ParametersTag[i], All.ParametersValue[i], H5T_NATIVE_INT);
            break;
        }
    }
}

/*---------------------- Routine find a block in a snapfile -------------------*/
void IO_Def::find_block(char *label, FILE *fd)
{
  unsigned int blocksize = 0, blksize;
  char blocklabel[5]     = {"    "};

#define FBSKIP                              \
  {                                         \
    my_fread(&blksize, sizeof(int), 1, fd); \
  }

  rewind(fd);

  while(!feof(fd) && blocksize == 0)
    {
      FBSKIP;

      if(blksize != 8)
        {
          Terminate("Incorrect Format (blksize=%u)!\n", blksize);
        }
      else
        {
          my_fread(blocklabel, LABEL_LEN * sizeof(char), 1, fd);
          my_fread(&blocksize, sizeof(int), 1, fd);

          FBSKIP;

          if(strncmp(label, blocklabel, LABEL_LEN) != 0)
            {
              fseek(fd, blocksize, 1);
              blocksize = 0;
            }
        }
    }
  if(feof(fd))
    Terminate("Block '%c%c%c%c' not found !\n", label[0], label[1], label[2], label[3]);
}

void IO_Def::read_segment(const char *fname, int type, long long offset, long long count, int num_files)
{
  long long nleft   = count;
  long long offleft = offset;
  long long nskip   = 0;

  for(int filenr = 0; filenr < num_files && nleft > 0; filenr++)
    {
      long long nloc = ntype_in_files[filenr * N_DataGroups + type];

      if(nloc > offleft)  // we may have something in this file
        {
          long long nread;

          if(nloc - offleft > nleft)  // there are more particles in the file then we need
            nread = nleft;
          else
            nread = nloc - offleft;

          /* now read partial list */
          read_single_file_segment(fname, filenr, type, offleft, nread, nskip, num_files);

          nleft -= nread;
          nskip += nread;
          offleft += nread;
        }

      offleft -= nloc;
    }

  if(nleft > 0)
    {
      for(int filenr = 0; filenr < num_files; filenr++)
        {
          long long nloc = ntype_in_files[filenr * N_DataGroups + type];
          printf("filenr=%d:  nloc=%lld\n", filenr, nloc);
        }
      Terminate("Not all desired entries read: nleft=%lld  type=%d\n", nleft, type);
    }
}

void IO_Def::read_single_file_segment(const char *basename, int filenr, int type, long long offset, unsigned long long count,
                                      long long storage_offset, int num_files)
{
  int bytes_per_blockelement_in_file = 0;
  hid_t hdf5_file = 0, hdf5_grp = 0, hdf5_dataspace_in_file;
  hid_t hdf5_dataspace_in_memory, hdf5_dataset;
  FILE *fd = 0;
  char fname[MAXLEN_PATH_EXTRA];

  if(num_files > 1)
    {
      if(file_format == FILEFORMAT_HDF5)
        snprintf(fname, MAXLEN_PATH_EXTRA, "%s.%d.hdf5", basename, filenr);
      else
        snprintf(fname, MAXLEN_PATH_EXTRA, "%s.%d", basename, filenr);
    }
  else
    {
      if(file_format == FILEFORMAT_HDF5)
        snprintf(fname, MAXLEN_PATH_EXTRA, "%s.hdf5", basename);
      else
        snprintf(fname, MAXLEN_PATH_EXTRA, "%s", basename);
    }

  /* open file  */
  if(file_format == FILEFORMAT_HDF5)
    {
      hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
      char buf[MAXLEN_PATH];
      get_datagroup_name(type, buf);
      hdf5_grp = my_H5Gopen(hdf5_file, buf);
    }
  else if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
    {
      if(!(fd = fopen(fname, "r")))
        Terminate("can't open file `%s' for reading initial conditions.\n", fname);

      unsigned int blksize1, blksize2;

      if(file_format == FILEFORMAT_LEGACY2)
        {
          int nextblock;
          char label[LABEL_LEN];
          my_fread(&blksize1, sizeof(int), 1, fd);
          my_fread(&label, sizeof(char), LABEL_LEN, fd);
          my_fread(&nextblock, sizeof(int), 1, fd);
          my_fread(&blksize2, sizeof(int), 1, fd);
        }

      my_fread(&blksize1, sizeof(int), 1, fd);
      my_fread(header_buf, header_size, 1, fd);
      my_fread(&blksize2, sizeof(int), 1, fd);
    }
  else
    Terminate("unknown ICFormat");

  long long npart[N_DataGroups];
  for(int i = 0; i < N_DataGroups; i++)
    npart[i] = ntype_in_files[filenr * N_DataGroups + i];

  for(int blocknr = 0; blocknr < N_IO_Fields; blocknr++)
    {
      if(IO_Fields[blocknr].type_in_memory != MEM_MY_FILEOFFSET)
        {
          unsigned int blksize1, blksize2;
          int typelist[N_DataGroups];
          int bytes_per_blockelement = get_bytes_per_memory_blockelement(blocknr, 1);
          long long npart_in_block   = get_particles_in_block(blocknr, npart, &typelist[0]);
          hid_t hdf5_memory_datatype = get_hdf5_memorytype_of_block(blocknr);
          char dname[MAXLEN_PATH];
          get_dataset_name(blocknr, dname);

          if(npart_in_block > 0 && typelist[type] > 0 && IO_Fields[blocknr].read_flag != SKIP_ON_READ)
            {
              if(file_format == FILEFORMAT_LEGACY2)
                {
                  char expected_label[LABEL_LEN + 1];
                  get_Tab_IO_Label(blocknr, expected_label);
                  find_block(expected_label, fd);
                }

              if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
                {
                  my_fread(&blksize1, sizeof(int), 1, fd);
                  bytes_per_blockelement_in_file = blksize1 / npart_in_block;
                }

              void *CommBuffer = (char *)Mem.mymalloc("CommBuffer", bytes_per_blockelement * count);

              if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
                {
                  fseek(fd, bytes_per_blockelement_in_file * offset, SEEK_CUR);

                  if(bytes_per_blockelement_in_file != bytes_per_blockelement)
                    {
                      char *CommAuxBuffer = (char *)Mem.mymalloc("CommAuxBuffer", bytes_per_blockelement_in_file * count);
                      my_fread(CommAuxBuffer, bytes_per_blockelement_in_file, count, fd);
                      type_cast_data((char *)CommAuxBuffer, bytes_per_blockelement_in_file, (char *)CommBuffer, bytes_per_blockelement,
                                     count, blocknr);
                      Mem.myfree(CommAuxBuffer);
                    }
                  else
                    my_fread(CommBuffer, bytes_per_blockelement, count, fd);

                  fseek(fd, bytes_per_blockelement_in_file * (npart_in_block - offset - count), SEEK_CUR);

                  my_fread(&blksize2, sizeof(int), 1, fd);
                  if(blksize1 != blksize2)
                    {
                      char buf[MAXLEN_PATH_EXTRA];
                      snprintf(buf, MAXLEN_PATH_EXTRA,
                               "incorrect block-sizes detected!\n Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask, blocknr,
                               blksize1, blksize2);
                      if(blocknr == 2) /* block number 2 is always IDs */
                        strcat(buf, "Possible mismatch of 32bit and 64bit ID's in IC file and GADGET compilation !\n");
                      Terminate(buf);
                    }
                }

              if(file_format == FILEFORMAT_HDF5)
                {
                  hdf5_dataset = my_H5Dopen_if_existing(hdf5_grp, dname);
                  int rank;
                  hsize_t dims[2], nelem[2], start[2];

                  dims[0] = npart[type];
                  dims[1] = get_values_per_blockelement(blocknr);
                  if(dims[1] == 1)
                    rank = 1;
                  else
                    rank = 2;

                  hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);

                  dims[0]                  = count;
                  hdf5_dataspace_in_memory = my_H5Screate_simple(rank, dims, NULL);

                  start[0] = offset;
                  start[1] = 0;

                  nelem[0] = count;
                  nelem[1] = get_values_per_blockelement(blocknr);

                  my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, nelem, NULL);

                  // Test if dataset was present
                  if(hdf5_dataset < 0)
                    {
                      // no, pad with zeros
                      printf("\tDataset %s not present for particle type %d, using zero.\n", dname, type);
                      memset(CommBuffer, 0, dims[0] * dims[1] * my_H5Tget_size(hdf5_memory_datatype));
                    }
                  else
                    {
                      hid_t hdf5_file_datatype = H5Dget_type(hdf5_dataset);
                      byte_count += dims[0] * dims[1] * my_H5Tget_size(hdf5_file_datatype); /* for I/O performance measurement */
                      H5Tclose(hdf5_file_datatype);

                      my_H5Dread(hdf5_dataset, hdf5_memory_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT,
                                 CommBuffer, dname);
                      my_H5Dclose(hdf5_dataset, dname);
                    }
                  my_H5Sclose(hdf5_dataspace_in_memory, H5S_SIMPLE);
                  my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                }

              empty_read_buffer(blocknr, storage_offset, count, type, 0, CommBuffer);

              Mem.myfree(CommBuffer);
            }
          else
            {
              if(file_format == FILEFORMAT_LEGACY1 && npart_in_block > 0)
                {
                  my_fread(&blksize1, sizeof(int), 1, fd);
                  bytes_per_blockelement_in_file = blksize1 / npart_in_block;
                  fseek(fd, bytes_per_blockelement_in_file * npart_in_block, SEEK_CUR);
                  my_fread(&blksize2, sizeof(int), 1, fd);
                  if(blksize1 != blksize2)
                    {
                      char buf[MAXLEN_PATH_EXTRA];
                      snprintf(buf, MAXLEN_PATH_EXTRA,
                               "incorrect block-sizes detected!\n Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask, blocknr,
                               blksize1, blksize2);
                      if(blocknr == 2) /* block number 2 is always IDs */
                        strcat(buf, "Possible mismatch of 32bit and 64bit ID's in IC file and GADGET compilation !\n");
                      Terminate(buf);
                    }
                }
            }
        }
    }

  if(file_format == FILEFORMAT_LEGACY1 || file_format == FILEFORMAT_LEGACY2)
    fclose(fd);

  if(file_format == FILEFORMAT_HDF5)
    {
      char buf[MAXLEN_PATH];
      get_datagroup_name(type, buf);
      my_H5Gclose(hdf5_grp, buf);
      my_H5Fclose(hdf5_file, fname);
    }

  read_increase_numbers(type, count);
}

void IO_Def::rename_file_to_bak_if_it_exists(char *fname)
{
  char fin[MAXLEN_PATH], buf[MAXLEN_PATH_EXTRA];
  
  strncpy(fin, fname, MAXLEN_PATH);
  fin[MAXLEN_PATH - 1] = 0;
  
  char *p = strrchr(fin, '/');
  if(p)
    {
      *p = 0;
      snprintf(buf, MAXLEN_PATH_EXTRA, "%s/bak-%s", fin, p + 1);
    }
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "bak-%s", fname);

  if(FILE *fcheck = fopen(fname, "r"))  // check if file already exists, if yes, try to rename the existing file
    {
      fclose(fcheck);

      if(!fopen(buf, "r"))  // only do the rename of the old file if the back-up file doesn't exist yet
        {
          mpi_printf("%s rename '%s' to '%s'\n", info, fname, buf);
          rename(fname, buf);
        }
    }
}

void IO_Def::alloc_and_read_ntype_in_files(const char *fname, int num_files)
{
  ntype_in_files = (long long *)Mem.mymalloc_movable(&ntype_in_files, "ntype_in_files", num_files * N_DataGroups * sizeof(long long));

  for(int filenr = 0; filenr < num_files; filenr++)
    {
      char buf[MAXLEN_PATH_EXTRA];

      if(num_files > 1)
        {
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d", fname, filenr);
          if(file_format == 3)
            snprintf(buf, MAXLEN_PATH_EXTRA, "%s.%d.hdf5", fname, filenr);
        }
      else
        {
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s", fname);
          if(file_format == 3)
            snprintf(buf, MAXLEN_PATH_EXTRA, "%s.hdf5", fname);
        }

      if(file_format == 3)
        {
          read_header_fields(buf);
        }
      else if(file_format == 1 || file_format == 2)
        {
          FILE *fd = 0;

          if(!(fd = fopen(buf, "r")))
            Terminate("can't open file `%s' for reading initial conditions.\n", fname);

          int blksize1, blksize2;

          if(file_format == 2)
            {
              char label[4];
              int nextblock;
              my_fread(&blksize1, sizeof(int), 1, fd);
              my_fread(&label, sizeof(char), 4, fd);
              my_fread(&nextblock, sizeof(int), 1, fd);
              printf("READIC: Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3], nextblock);
              my_fread(&blksize2, sizeof(int), 1, fd);
            }

          my_fread(&blksize1, sizeof(int), 1, fd);
          my_fread(header_buf, header_size, 1, fd);
          my_fread(&blksize2, sizeof(int), 1, fd);

          if(blksize1 != blksize2)
            Terminate("incorrect header format, blksize1=%d blksize2=%d  header_size=%d\n", blksize1, blksize2, (int)header_size);

          fclose(fd);
        }

      long long n_type[N_DataGroups], npart[N_DataGroups];

      read_file_header(fname, filenr, 0, 0, n_type, npart, NULL);

      for(int type = 0; type < N_DataGroups; type++)
        ntype_in_files[filenr * N_DataGroups + type] = npart[type];
    }
}
