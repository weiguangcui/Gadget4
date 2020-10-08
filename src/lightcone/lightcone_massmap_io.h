/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file lightcone_massmap_io.h
 *
 * \brief declares class used for I/O of lightcone mass projections
 */

#ifndef LIGHTCONE_MASSMAP_IO_H
#define LIGHTCONE_MASSMAP_IO_H

#include "gadgetconfig.h"

#if defined(LIGHTCONE) && defined(LIGHTCONE_MASSMAPS)

#include "../data/intposconvert.h"
#include "../data/mmparticles.h"
#include "../io/io.h"
#include "../lightcone/lightcone.h"

class lightcone_massmap_io : public IO_Def
{
 public:
  lightcone_massmap_io(mmparticles *Mp_ptr, lightcone *LightCone_ptr, MPI_Comm comm, int format);

  void lightcone_massmap_save(int num);

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
    int nside; /* healpix nside */
    int npix_local;
    int npix_total;
    int num_files;
    double AscaleStart;
    double AscaleEnd;
    double ComDistStart;
    double ComDistEnd;
  };
  io_header header; /**< holds header for snapshot files */

 private:
  mmparticles *Mp;
  lightcone *LightCone;

  int selected_bnd;
};

#endif
#endif
