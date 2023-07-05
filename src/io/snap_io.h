/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file snap_io.h
 *
 *  \brief declares class used for I/O of snapshot files
 */

#ifndef SNAP_READ_WRITE_H
#define SNAP_READ_WRITE_H

#include "gadgetconfig.h"

#include "../data/intposconvert.h"
#include "../data/simparticles.h"
#include "../io/io.h"
#include "../mergertree/mergertree.h"

class snap_io : public IO_Def
{
 public:
  void init_basic(simparticles *Sp_ptr);

  snap_io(simparticles *Sp_ptr, MPI_Comm comm, int format) : IO_Def(comm, format) { init_basic(Sp_ptr); }

#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
  void init_extra(simparticles *Sp_ptr, mergertree *MergerTree_ptr);
  snap_io(simparticles *Sp_ptr, mergertree *MergerTree_ptr, MPI_Comm comm, int format) : IO_Def(comm, format)
  {
    init_basic(Sp_ptr);
    init_extra(Sp_ptr, MergerTree_ptr);
  }
#endif

  void write_snapshot(int num, mysnaptype snap_type);
  void read_snapshot(int num, mysnaptype snap_type);

  int acquire_basic_treeinfo(int num, mysnaptype loc_snap_type);
  void free_basic_treeinfo(void);
  long long load_orphans(int num, long long treenr, int num_files);
  void free_orphans(void);

  void read_ic(const char *fname);

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
    int npart[NTYPES_HEADER];                      /**< number of particles of each type in this file */
    double mass[NTYPES_HEADER];                    /**< mass of particles of each type. If 0, then the masses are explicitly
                                                          stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                                   /**< time of snapshot file */
    double redshift;                               /**< redshift of snapshot file */
    int flag_sfr;                                  /**< flags whether the simulation was including star formation */
    int flag_feedback;                             /**< flags whether feedback was included (obsolete) */
    unsigned int npartTotalLowWord[NTYPES_HEADER]; /**< total number of particles of each type in this snapshot. This can be
                                       different from npart if one is dealing with a multi-file snapshot. */
    int flag_cooling;                              /**< flags whether cooling was included  */
    int num_files;                                 /**< number of files in multi-file snapshot */
    double BoxSize;                                /**< box-size of simulation in case periodic boundaries were used */
    double Omega0;                                 /**< matter density in units of critical density */
    double OmegaLambda;                            /**< cosmological constant parameter */
#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
    long long Ntrees;     // this replaces the storage space for HubbleParam
    long long TotNtrees;  // this replaces the storage space for Hubble
#else
    double HubbleParam; /**< little 'h' to scale units systems */
    double Hubble;      /**< Hubble constant in internal units */
#endif
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

    long long npartTotal[NTYPES_HEADER]; /**< fills to 256 Bytes, and for compatability with Gadget2/3 */
  };

#else
  /* new simplified header format */
  struct io_header
  {
    long long npart[NTYPES_HEADER];      /**< number of particles of each type in this file */
    long long npartTotal[NTYPES_HEADER]; /**< total number of particles of each type in this snapshot. This can be
                                           different from npart if one is dealing with a multi-file snapshot. */
    double mass[NTYPES_HEADER];          /**< mass of particles of each type. If 0, then the masses are explicitly
                                                stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                         /**< time of snapshot file */
    double redshift;                     /**< redshift of snapshot file */
    double BoxSize;                      /**< box-size of simulation in case periodic boundaries were used */
    int num_files;                       /**< number of files in multi-file snapshot */

    long long Ntrees;
    long long TotNtrees;
  };

#endif

  io_header header; /**< holds header for snapshot files */

 private:
  simparticles *Sp;
#ifdef MERGERTREE
  mergertree *MergerTree;
#endif

  mysnaptype snap_type;

  long long ntot_type_all[NTYPES]; /**< contains the global number of particles of each type in the snapshot file */

#ifndef OUTPUT_COORDINATES_AS_INTEGERS
  struct ptmp_data
  {
    MyDouble Pos[3];
  };
  ptmp_data *Ptmp;
#endif

  void snap_init_domain_mapping(void);
  void read_increase_particle_numbers(int type, int n_for_this_task);

  /*
   * special input/output functions for certain fields
   */
#ifndef OUTPUT_COORDINATES_AS_INTEGERS
  static void io_func_pos(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    /* note: we know that components==3 here */
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)
      {
        MyDouble *out_buffer = (MyDouble *)buffer;

        /* converts the integer coordinates to floating point */
        thisobj->Sp->intpos_to_pos(thisobj->Sp->P[particle].IntPos, out_buffer);
      }
    else
      {
        MyDouble *in_buffer = (MyDouble *)buffer;

        /* note: for non-periodic positions, the conversion to integer coordinates is undefined only after the initial read.
         * We therefore store the coordinates first in a temporary array */

        for(int k = 0; k < 3; k++)
          thisobj->Ptmp[particle].Pos[k] = in_buffer[k];

#if defined(NGENIC) && !defined(CREATE_GRID)  // This is meant to become active when a glass file is used for IC creation
        if(All.RestartFlag == RST_BEGIN || All.RestartFlag == RST_CREATEICS)
          {
            double fac = All.BoxSize / thisobj->header.BoxSize;
#ifdef TILING
            fac /= TILING;
#endif
            for(int k = 0; k < 3; k++)
              thisobj->Ptmp[particle].Pos[k] *= fac;  // scale the glass file to the right size
          }
#endif

#ifdef SQUASH_TEST
        thisobj->Ptmp[particle].Pos[1] *= 1.0 / 4;
        thisobj->Ptmp[particle].Pos[2] *= 1.0 / 16;
#endif
      }
  }
#endif

  static void io_func_intpos(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    /* note: we know that components==3 here */
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)
      {
        MyIntPosType *out_buffer = (MyIntPosType *)buffer;

        /* converts the integer coordinates to integer by subtracting a possible randomization shift */
        thisobj->Sp->intpos_to_intpos(thisobj->Sp->P[particle].IntPos, out_buffer);
      }
    else
      {
        MyIntPosType *in_buffer = (MyIntPosType *)buffer;

        for(int k = 0; k < 3; k++)
          thisobj->Sp->P[particle].IntPos[k] = in_buffer[k];
      }
  }

  static void io_func_vel(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)
      {
        MyFloat *out_buffer = (MyFloat *)buffer;
        for(int k = 0; k < 3; k++)
          {
            out_buffer[k] = thisobj->Sp->P[particle].Vel[k];

            /* we are using p = a^2 * xdot internally as velocity unit. Convert to legacy Gadget velocity units */
            if(All.RestartFlag != RST_CONVERTSNAP)
              out_buffer[k] *= sqrt(All.cf_a3inv);
          }
      }
    else
      {
        MyFloat *in_buffer = (MyFloat *)buffer;
        for(int k = 0; k < components; k++)
          {
            thisobj->Sp->P[particle].Vel[k] = in_buffer[k];
          }
      }
  }

  static void io_func_id(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)
      {
        MyIDType *out_buffer = (MyIDType *)buffer;
        out_buffer[0]        = thisobj->Sp->P[particle].ID.get();
      }
    else
      {
        MyIDType *in_buffer = (MyIDType *)buffer;
        thisobj->Sp->P[particle].ID.set(in_buffer[0]);
      }
  }

  static void io_func_mass(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)
      {
        MyDouble *out_buffer = (MyDouble *)buffer;
        out_buffer[0]        = thisobj->Sp->P[particle].getMass();
      }
    else
      {
        MyDouble *in_buffer = (MyDouble *)buffer;
        thisobj->Sp->P[particle].setMass(in_buffer[0]);
      }
  }

  static void io_func_u(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)
      {
        MyFloat *out_buffer = (MyFloat *)buffer;
        out_buffer[0]       = thisobj->Sp->get_utherm_from_entropy(particle);
      }
    else
      {
        MyFloat *in_buffer                  = (MyFloat *)buffer;
        thisobj->Sp->SphP[particle].Entropy = in_buffer[0];
      }
  }

#ifdef STARFORMATION
  static void io_func_sfr(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)
      {
        MyFloat *out_buffer = (MyFloat *)buffer;
        out_buffer[0]       = thisobj->Sp->SphP[particle].Sfr;
      }
    else
      {
        MyFloat *in_buffer              = (MyFloat *)buffer;
        thisobj->Sp->SphP[particle].Sfr = in_buffer[0];
      }
  }

  static void io_func_metallicity(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)
      {
        MyFloat *out_buffer = (MyFloat *)buffer;
        if(thisobj->Sp->P[particle].getType() == 0)
          {
            out_buffer[0] = thisobj->Sp->SphP[particle].Metallicity;
          }
        else
          {
            out_buffer[0] = thisobj->Sp->P[particle].Metallicity;
          }
      }
    else
      {
        MyFloat *in_buffer                   = (MyFloat *)buffer;
        thisobj->Sp->P[particle].Metallicity = in_buffer[0];
      }
  }
#endif

#ifdef OUTPUT_ACCELERATION
  static void io_func_accel(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj = (snap_io *)ptr;

    if(mode == 0)  // writing
      {
        MyFloat *out_buffer = (MyFloat *)buffer;
        if(All.RestartFlag != RST_CONVERTSNAP)
          for(int k = 0; k < 3; k++)
            out_buffer[k] = All.cf_a2inv * thisobj->Sp->P[particle].GravAccel[k];
        else
          for(int k = 0; k < 3; k++)
            out_buffer[k] = thisobj->Sp->P[particle].GravAccel[k];
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
        if(All.RestartFlag != RST_CONVERTSNAP)
          for(int k = 0; k < 3; k++)
            out_buffer[k] += All.cf_a2inv * thisobj->Sp->P[particle].GravPM[k];
        else
          for(int k = 0; k < 3; k++)
            out_buffer[k] += thisobj->Sp->P[particle].GravPM[k];
#endif

        for(int k = 0; k < 3; k++)
          out_buffer[k] /= All.accel_normalize_fac;
      }
    else  // reading
      {
        MyFloat *in_buffer = (MyFloat *)buffer;
        for(int k = 0; k < components; k++)
          thisobj->Sp->P[particle].GravAccel[k] = All.accel_normalize_fac * in_buffer[k];
      }
  }

#endif

#ifdef OUTPUT_PRESSURE
  static void io_func_pressure(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj    = (snap_io *)ptr;
    MyFloat *out_buffer = (MyFloat *)buffer;

    out_buffer[0] = thisobj->Sp->SphP[particle].get_pressure();
  }
#endif

#ifdef OUTPUT_TIMESTEP
  static void io_func_timestep(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj    = (snap_io *)ptr;
    MyFloat *out_buffer = (MyFloat *)buffer;

    out_buffer[0] = (thisobj->Sp->P[particle].TimeBinGrav ? (((integertime)1) << thisobj->Sp->P[particle].TimeBinGrav) : 0) *
                    All.Timebase_interval;
  }

  static void io_func_timestephydro(IO_Def *ptr, int particle, int components, void *buffer, int mode)
  {
    snap_io *thisobj    = (snap_io *)ptr;
    MyFloat *out_buffer = (MyFloat *)buffer;

    out_buffer[0] =
        (thisobj->Sp->P[particle].getTimeBinHydro() ? (((integertime)1) << thisobj->Sp->P[particle].getTimeBinHydro()) : 0) *
        All.Timebase_interval;
  }
#endif
};

#endif /* SNAP_READ_WRITE_H */
