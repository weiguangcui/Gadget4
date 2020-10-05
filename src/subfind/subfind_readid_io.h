/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_readid_io.h
 *
 *  \brief declaration of a class for reading particle IDs of group catalogues
 */

#ifndef SUBREADID_IO_H
#define SUBREADID_IO_H

#include "gadgetconfig.h"

#ifdef SUBFIND_ORPHAN_TREATMENT

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../data/particle_data.h"
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

class subreadid_io : public IO_Def
{
 private:
  idstoredata *IdStore;

 public:
  subreadid_io(idstoredata *IdStore_ptr, MPI_Comm comm, int format);

  void previously_bound_read_snap_ids(int num);

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
#ifdef GADGET2_HEADER
  struct io_header
  {
    int npart[NTYPES_HEADER];                       /**< number of particles of each type in this file */
    double mass[NTYPES_HEADER];                     /**< mass of particles of each type. If 0, then the masses are explicitly
                                                           stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                                    /**< time of snapshot file */
    double redshift;                                /**< redshift of snapshot file */
    int flag_sfr;                                   /**< flags whether the simulation was including star formation */
    int flag_feedback;                              /**< flags whether feedback was included (obsolete) */
    unsigned int npartTotal[NTYPES_HEADER];         /**< total number of particles of each type in this snapshot. This can be
                                               different from npart if one is dealing with a multi-file snapshot. */
    int flag_cooling;                               /**< flags whether cooling was included  */
    int num_files;                                  /**< number of files in multi-file snapshot */
    double BoxSize;                                 /**< box-size of simulation in case periodic boundaries were used */
    double Omega0;                                  /**< matter density in units of critical density */
    double OmegaLambda;                             /**< cosmological constant parameter */
    double HubbleParam;                             /**< Hubble parameter in units of 100 km/sec/Mpc */
    double Hubble;                                  /**< Hubble constant in internal units */
    unsigned int npartTotalHighWord[NTYPES_HEADER]; /**< High word of the total number of particles of each type */
    int flag_entropy_instead_u;                     /**< flags that IC-file contains entropy instead of u */
    int flag_doubleprecision;                       /**< flags that snapshot contains double-precision instead of single precision */
    int flag_ic_info;        /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                    or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                    For snapshots files, the value informs whether the simulation was evolved from
                                    Zeldoch or 2lpt ICs. Encoding is as follows:
                                      FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                      FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                */
    float lpt_scalingfactor; /*!< scaling factor for 2lpt initial conditions */
    char fill[48];           /**< fills to 256 Bytes, for compatability with Gadget2/3 */
  };
  io_header header; /**< holds header for snapshot files */
#else

  /* new simplified header */
  struct io_header
  {
    long long npart[NTYPES_HEADER]; /**< number of particles of each type in this file */
    long long npartTotal[NTYPES_HEADER];
    double mass[NTYPES_HEADER]; /**< mass of particles of each type. If 0, then the masses are explicitly
                                       stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                /**< time of snapshot file */
    double redshift;            /**< redshift of snapshot file */
    double BoxSize;             /**< box-size of simulation in case periodic boundaries were used */
    int num_files;              /**< number of files in multi-file snapshot */
  };
  io_header header; /**< holds header for snapshot files */

#endif
};

#endif

#endif /* READSNAP_IO_H */
