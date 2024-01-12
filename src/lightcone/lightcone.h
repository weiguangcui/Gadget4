/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  lightcone.h
 *
 *  \brief declares a class for accumulating particles on the lightcone
 */

#ifndef LIGHTCONE_H
#define LIGHTCONE_H

#include "gadgetconfig.h"

#ifdef LIGHTCONE

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/lcparticles.h"
#include "../data/mmparticles.h"
#include "../data/simparticles.h"
#include "../data/symtensors.h"
#include "../mpi_utils/setcomm.h"
#include "../time_integration/driftfac.h"

#ifdef LIGHTCONE_MASSMAPS
#include <chealpix.h>
#endif

#ifndef LIGHTCONE_MAX_BOXREPLICAS
#define LIGHTCONE_MAX_BOXREPLICAS 1000
#endif

#ifndef LIGHTCONE_ORDER_NSIDE
#define LIGHTCONE_ORDER_NSIDE 256
#endif

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
#define LIGHTCONE_MAX_NUMBER_ORIGINS 32
#else
#define LIGHTCONE_MAX_NUMBER_ORIGINS 1
#endif

#define LC_TYPE_FULLSKY 0
#define LC_TYPE_OCTANT 1
#define LC_TYPE_PENCIL 2
#define LC_TYPE_DISK 3
#define LC_TYPE_SQUAREMAP 4

class lightcone : public parameters
{
 public:
  simparticles *Sp;

#if defined(LIGHTCONE_PARTICLES)
  lcparticles *Lp;
#endif
#if defined(LIGHTCONE_MASSMAPS)
  mmparticles *Mp;
#endif

 public:
#if defined(LIGHTCONE_PARTICLES) && defined(LIGHTCONE_MASSMAPS) /* both particles and massmaps */
  lightcone(MPI_Comm comm, simparticles *Sp_ptr, lcparticles *Lp_ptr, mmparticles *Mp_ptr) : parameters(comm)
  {
    Sp = Sp_ptr;
    Lp = Lp_ptr;
    Mp = Mp_ptr;

    Sp->LightCone = this;
  }
#else
#if defined(LIGHTCONE_PARTICLES)
  lightcone(MPI_Comm comm, simparticles *Sp_ptr, lcparticles *Lp_ptr) : parameters(comm) /* only particles */
  {
    Sp = Sp_ptr;
    Lp = Lp_ptr;

    Sp->LightCone = this;
  }
#else
#if defined(LIGHTCONE_MASSMAPS)
  lightcone(MPI_Comm comm, simparticles *Sp_ptr, mmparticles *Mp_ptr) : parameters(comm) /* only massmaps */
  {
    Sp = Sp_ptr;
    Mp = Mp_ptr;

    Sp->LightCone = this;
  }
#endif
#endif
#endif

  void lightcone_init_intposconverter(double linklength);

  int lightcone_test_for_particle_addition(particle_data *P, integertime time0, integertime time1, double dt_drift);

  void register_parameters(void);

  void makeimage(int argc, char **argv);

#ifdef LIGHTCONE_PARTICLES
  int Nlightcones;

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  int NlightconeOrigins;

  struct cone_origin
  {
    double PosOrigin[3];
  };
  cone_origin *ConeOrigins;
#endif

  struct cone_data
  {
    double Astart;
    double Aend;

    integertime Time_start;
    integertime Time_end;

    double ComDistStart;
    double ComDistEnd;

    int LcType;
    int OnlyMostBoundFlag;

    int OctantNr;

    vector<double> PencilDir;
    double PencilAngle;
    double PencilAngleRad;

    vector<double> DiskNormal;
    double DiskThickness;

    vector<double> SquareMapXdir;
    vector<double> SquareMapYdir;
    vector<double> SquareMapZdir;

    double SquareMapAngle;
    double SquareMapAngleRad;

    char Tag[100];

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
    int OriginIndex;
#endif
  };
  cone_data *Cones;

  double ConeGlobAstart;
  double ConeGlobAend;

  double ConeGlobTime_start;
  double ConeGlobTime_end;

  double ConeGlobComDistStart;
  double ConeGlobComDistEnd;

#endif

  struct boxlist
  {
    int i, j, k; /* displacement of principal box */
    double Rmin; /* minimum comoving distance of this box */
    double Rmax; /* minimum comoving distance of this box */
  };

  struct boxorigin
  {
    boxlist *BoxList;
    int NumBoxes;
  };

  boxorigin BoxOrigin[LIGHTCONE_MAX_NUMBER_ORIGINS];

  void lightcone_init_geometry(char *fname);
  void lightcone_add_position_particles(particle_data *P, double *pos, double ascale, int oindex);
  int lightcone_init_times(void);
  int lightcone_init_boxlist(void);
  bool lightcone_is_cone_member(int i, int cone);
  bool lightcone_is_cone_member_basic(double ascale, vector<double> &pos, bool previously, int cone);

  bool lightcone_box_at_corner_overlaps_at_least_with_one_cone(double *corner, double &rmin, double &rmax, int oindex);
  void lightcone_clear_boxlist(double ascale);

  static bool lightcone_compare_BoxList_Rmax(const boxlist &a, const boxlist &b)
  {
    return a.Rmax > b.Rmax; /* sort in descending order */
  }

#ifdef LIGHTCONE_MASSMAPS

  double *MassMap;

  int NumMassMapBoundaries;
  double *MassMapBoundariesAscale;
  integertime *MassMapBoundariesTime;
  double *MassMapBoundariesComDist;

  void lightcone_init_massmaps(void);
  void lightcone_massmap_binning(void);
  void lightcone_massmap_flush(int dump_allowed_flag);
  int lightcone_add_position_massmaps(particle_data *P, double *pos, double ascale);
  int lightcone_massmap_report_boundaries(void);

  static bool compare_doubles(const double &a, const double &b) { return a < b; }
#endif
};

#endif
#endif
