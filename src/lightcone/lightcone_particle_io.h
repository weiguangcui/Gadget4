/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file lightcone_particle_io.h
 *
 * \brief declares class for I/O of lightcone particles
 */

#ifndef LIGHTCONE_IO_H
#define LIGHTCONE_IO_H

#include "gadgetconfig.h"

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES)

#include "../data/intposconvert.h"
#include "../data/lcparticles.h"
#include "../io/io.h"
#include "../lightcone/lightcone.h"
#include "../mergertree/mergertree.h"

class lightcone_particle_io : public IO_Def
{
 public:
#ifdef MERGERTREE
  lightcone_particle_io(lcparticles *Lp_ptr, lightcone *LightCone_ptr, mergertree *MergerTree_ptr, MPI_Comm comm, int format);
#else
  lightcone_particle_io(lcparticles *Lp_ptr, lightcone *LightCone_ptr, MPI_Comm comm, int format);
#endif

  void lightcone_save(int num, int conenr, bool reordered_flag);
  void lightcone_read(int num, int conenr);

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
    long long npart[NTYPES]; /**< number of particles of each type in this file */
    long long npartTotal[NTYPES];

    long long Ntrees;
    long long TotNtrees;

    int Npix;
    int TotNpix;

    int num_files;

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
    double Origin[3];
#endif
  };
  io_header header; /**< holds header for snapshot files */

 private:
  lcparticles *Lp;
  lightcone *LightCone;
#ifdef MERGERTREE
  mergertree *MergerTree;
#endif

  int cone;
  bool reorder_flag;
  long long ntot_type_all[NTYPES];

  /*
   * special input/output functions for certain fields
   */

  static void io_func_pos(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    /* note: we know that components==3 here */
    lightcone_particle_io *thisobj = (lightcone_particle_io *)ptr;

    if(mode == 0)
      {
        MyDouble *out_buffer = (MyDouble *)buffer;

        double xyz[3];
        thisobj->Lp->signedintpos_to_pos((MySignedIntPosType *)thisobj->Lp->P[particle].IntPos, xyz);

        for(int k = 0; k < 3; k++)
          out_buffer[k] = xyz[k];
      }
    else
      {
        MyDouble *in_buffer = (MyDouble *)buffer;

        /* note: for non-periodic positions, the conversion to integer coordinates is undefined only after the initial read.
         * We therefore store the coordinates first in a temporary array */

        double xyz[3];

        for(int k = 0; k < 3; k++)
          xyz[k] = in_buffer[k];

        /* converts floating point representation to integers */
        thisobj->Lp->pos_to_signedintpos(xyz, (MySignedIntPosType *)thisobj->Lp->P[particle].IntPos);
      }
  }

#if defined(LIGHTCONE_OUTPUT_ACCELERATIONS) && defined(OUTPUT_ACCELERATIONS_IN_HALF_PRECISION)
  static void io_func_accel(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    lightcone_particle_io *thisobj = (lightcone_particle_io *)ptr;

    if(mode == 0)  // writing
      {
        MyFloat *out_buffer = (MyFloat *)buffer;
        for(int k = 0; k < 3; k++)
          out_buffer[k] = thisobj->Lp->P[particle].GravAccel[k] / All.accel_normalize_fac;
      }
    else  // reading
      {
        MyFloat *in_buffer = (MyFloat *)buffer;
        for(int k = 0; k < components; k++)
          thisobj->Lp->P[particle].GravAccel[k] = All.accel_normalize_fac * in_buffer[k];
      }
  }
#endif

  static void io_func_id(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    lightcone_particle_io *thisobj = (lightcone_particle_io *)ptr;

    if(mode == 0)
      {
        MyIDType *out_buffer = (MyIDType *)buffer;
        out_buffer[0]        = thisobj->Lp->P[particle].ID.get();
      }
    else
      {
        MyIDType *in_buffer = (MyIDType *)buffer;
        thisobj->Lp->P[particle].ID.set(in_buffer[0]);
      }
  }

  static void io_func_mass(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    lightcone_particle_io *thisobj = (lightcone_particle_io *)ptr;

    if(mode == 0)
      {
        MyDouble *out_buffer = (MyDouble *)buffer;
        out_buffer[0]        = thisobj->Lp->P[particle].getMass();
      }
    else
      {
        MyDouble *in_buffer = (MyDouble *)buffer;
        thisobj->Lp->P[particle].setMass(in_buffer[0]);
      }
  }
};

#endif

#endif /* LIGHTCONE_IO_H */
