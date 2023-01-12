/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_readid_io.cc
 *
 *  \brief routines to read particle IDs of group catalogues
 */

#include "gadgetconfig.h"

#ifdef SUBFIND_ORPHAN_TREATMENT

#include <errno.h>
#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../data/particle_data.h"
#include "../fof/fof.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mergertree/io_readsnap.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../subfind/subfind_readid_io.h"
#include "../system/system.h"

subreadid_io::subreadid_io(idstoredata *IdStore_ptr, MPI_Comm comm, int format) : IO_Def(comm, format)
{
  IdStore = IdStore_ptr;

  this->N_IO_Fields  = 0;
  this->N_DataGroups = NTYPES;
  this->header_size  = sizeof(header);
  this->header_buf   = &header;
  this->type_of_file = FILE_IS_SNAPSHOT;
  snprintf(this->info, MAXLEN_PATH, "MERGERTREE: reading snapshot IDs");

  init_field("ID  ", "ParticleIDs", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, READ_IF_PRESENT, 1, A_IDS, IdStore->ID, NULL, ALL_TYPES, 0, 0, 0,
             0, 0, 0, 0);
}

/*! \brief This function reads initial conditions that are in on of the default file formats
 * of Gadget.
 *
 * Snapshot files can be used as input files.  However, when a
 * snapshot file is used as input, not all the information in the header is
 * used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 * Alternatively, the code can be started with restartflag 2, then snapshots
 * from the code can be used as initial conditions-files without having to
 * change the parameter file.  For gas particles, only the internal energy is
 * read, the density and mean molecular weight will be recomputed by the code.
 * When InitGasTemp>0 is given, the gas temperature will be initialized to this
 * value assuming a mean molecular weight either corresponding to complete
 * neutrality, or full ionization.
 *
 * \param fname file name of the ICs
 * \param readTypes readTypes is a bitfield that
 * determines what particle types to read, only if the bit
 * corresponding to a particle type is set, the corresponding data is
 * loaded, otherwise its particle number is set to zero. (This is
 * only implemented for HDF5 files.)
 */
void subreadid_io::previously_bound_read_snap_ids(int num)
{
  if(All.ICFormat < 1 || All.ICFormat > 4)
    Terminate("ICFormat=%d not supported.\n", All.ICFormat);

  char fname[MAXLEN_PATH_EXTRA], fname_multiple[MAXLEN_PATH_EXTRA];
  snprintf(fname_multiple, MAXLEN_PATH_EXTRA, "%s/snapdir_%03d/%s-prevmostboundonly_%03d", All.OutputDir, num, All.SnapshotFileBase,
           num);
  snprintf(fname, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);

  TIMER_START(CPU_SNAPSHOT);

  int num_files = find_files(fname, fname_multiple);

  if(num_files > 1)
    strcpy(fname, fname_multiple);

  /* we repeat reading the headers of the files two times. In the first iteration, only the
   * particle numbers ending up on each processor are assembled, followed by memory allocation.
   * In the second iteration, the data is actually read in.
   */
  for(int rep = 0; rep < 2; rep++)
    {
      IdStore->NumPart = 0;

      read_files_driver(fname, rep, num_files);

      /* now do the memory allocation */
      if(rep == 0)
        {
          IdStore->ID = (MyIDType *)Mem.mymalloc_movable_clear(&IdStore->ID, "IdStore->ID", IdStore->NumPart * sizeof(MyIDType));
        }
    }

  MPI_Barrier(Communicator);

  mpi_printf("READSNAPID: reading done.\n");

  TIMER_STOP(CPU_SNAPSHOT);
}

void subreadid_io::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type)
{ /* empty */
}

void subreadid_io::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *n_type, long long *ntot_type,
                                    int *nstart)
{
  if(ThisTask == readTask)
    {
      if(filenr == 0 && nstart == NULL)
        {
          mpi_printf(
              "\nREADSNAPID: filenr=%d, '%s' contains:\n"
              "READSNAPID: Type 0 (gas):   %8lld  (tot=%15lld) masstab= %g\n",
              filenr, fname, (long long)header.npart[0], (long long)header.npartTotal[0], All.MassTable[0]);

          for(int type = 1; type < NTYPES; type++)
            {
              mpi_printf("READSNAPID: Type %d:         %8lld  (tot=%15lld) masstab= %g\n", type, (long long)header.npart[type],
                         (long long)header.npartTotal[type], All.MassTable[type]);
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
      memmove(&IdStore->ID[nall], &IdStore->ID[0], IdStore->NumPart * sizeof(MyIDType));

      *nstart = 0;
    }
}

void subreadid_io::write_header_fields(hid_t handle)
{ /* empty */
}

/*! \brief This function reads the snapshot header in case of hdf5 files (i.e. format 3)
 *
 * \param fname file name of the snapshot as given in the parameter file
 */
void subreadid_io::read_header_fields(const char *fname)
{
  for(int i = 0; i < NTYPES; i++)
    {
      header.npart[i]      = 0;
      header.npartTotal[i] = 0;
      header.mass[i]       = 0;
    }

  hsize_t ntypes = NTYPES;

  hid_t hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t handle    = my_H5Gopen(hdf5_file, "/Header");

  /* check if the file in question actually has this number of types */
  hid_t hdf5_attribute = my_H5Aopen_name(handle, "NumPart_ThisFile");
  hid_t space          = H5Aget_space(hdf5_attribute);
  hsize_t dims, len;
  H5Sget_simple_extent_dims(space, &dims, &len);
  H5Sclose(space);
  if(len != ntypes)
    Terminate("Length of NumPart_ThisFile attribute (%d) does not match NTYPES(ICS) (%d)", (int)len, (int)ntypes);
  my_H5Aclose(hdf5_attribute, "NumPart_ThisFile");

  /* now read the header fields */

#ifdef GADGET2_HEADER
  read_vector_attribute(handle, "NumPart_ThisFile", header.npart, H5T_NATIVE_UINT, ntypes);
  read_vector_attribute(handle, "NumPart_Total", header.npartTotal, H5T_NATIVE_UINT, ntypes);
  read_scalar_attribute(handle, "Omega0", &header.Omega0, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "OmegaLambda", &header.OmegaLambda, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "HubbleParam", &header.HubbleParam, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "Hubble", &header.Hubble, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "Flag_Sfr", &header.flag_sfr, H5T_NATIVE_INT);
  read_scalar_attribute(handle, "Flag_Cooling", &header.flag_cooling, H5T_NATIVE_INT);
  read_scalar_attribute(handle, "Flag_Feedback", &header.flag_feedback, H5T_NATIVE_INT);
  read_scalar_attribute(handle, "Flag_DoublePrecision", &header.flag_doubleprecision, H5T_NATIVE_INT);
#else
  read_vector_attribute(handle, "NumPart_ThisFile", header.npart, H5T_NATIVE_UINT64, ntypes);
  read_vector_attribute(handle, "NumPart_Total", header.npartTotal, H5T_NATIVE_UINT64, ntypes);
#endif

  read_scalar_attribute(handle, "BoxSize", &header.BoxSize, H5T_NATIVE_DOUBLE);
  read_vector_attribute(handle, "MassTable", header.mass, H5T_NATIVE_DOUBLE, ntypes);
  read_scalar_attribute(handle, "Time", &header.time, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "Redshift", &header.redshift, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "NumFilesPerSnapshot", &header.num_files, H5T_NATIVE_INT);

  my_H5Gclose(handle, "/Header");
  my_H5Fclose(hdf5_file, fname);
}

int subreadid_io::get_filenr_from_header(void) { return header.num_files; }

void subreadid_io::set_filenr_in_header(int numfiles) { header.num_files = numfiles; }

void subreadid_io::read_increase_numbers(int type, int n_for_this_task) { IdStore->NumPart += n_for_this_task; }

void subreadid_io::get_datagroup_name(int type, char *buf) { snprintf(buf, MAXLEN_PATH, "/PartType%d", type); }

int subreadid_io::get_type_of_element(int index) { return 0; /* empty */ }

void subreadid_io::set_type_of_element(int index, int type)
{ /* empty */
}

void *subreadid_io::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_IDS:
        return (void *)(IdStore->ID + index);
      default:
        Terminate("we don't expect to get here");
    }

  return NULL;
}
#endif
