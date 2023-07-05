/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file lcparticles.h
 *
 *  \brief declares a class responsible for holding the (buffered) particles on the lightcone
 */

#ifndef LCPART_H
#define LCPART_H

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES)

#include "gadgetconfig.h"

#include <math.h>

#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/lightcone_particle_data.h"
#include "../data/macros.h"
#include "../data/mymalloc.h"
#include "../data/particle_data.h"
#include "../data/sph_particle_data.h"
#include "../mpi_utils/mpi_utils.h"
#include "../mpi_utils/setcomm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

class lcparticles : public intposconvert, public setcomm
{
 public:
  lcparticles(MPI_Comm comm) : setcomm(comm) {}

  int NumPart;    /**< number of particles on the LOCAL processor */
  int NumGas = 0; /* this is added here to simplify the template code */

  int MaxPart;        /**< This gives the maxmimum number of particles that can be stored on one processor. */
  int MaxPartSph = 0; /* this is added here to simplify the template code */

  long long TotNumPart;
  long long TotNumGas;

  typedef lightcone_particle_data pdata;

  /*! This structure holds all the information that is
   * stored for each particle of the simulation.
   */
  lightcone_particle_data *P; /*!< holds particle data on local processor */

  /* the following struture holds data that is stored for each SPH particle in addition to the collisionless
   * variables.
   */

  sph_particle_data *SphP = NULL; /* the current code does not yet support actual Sph particles on the lightcone */

  subfind_data *PS;

  int *HealPixTab_PartCount;
  int Npix;
  int FirstPix;
  int NpixLoc;

  /* This routine allocates memory for
   * particle storage, both the collisionless and the SPH particles.
   * The memory for the ordered binary tree of the timeline
   * is also allocated.
   */
  void allocate_memory(void)
  {
    if(MaxPart < BASENUMBER)
      MaxPart = BASENUMBER;

    P = (lightcone_particle_data *)Mem.mymalloc_movable_clear(&P, "P", MaxPart * sizeof(lightcone_particle_data));
  }

  void free_memory(void) { Mem.myfree(P); }

  void reallocate_memory_maxpart(int maxpartNew)
  {
    MaxPart = maxpartNew;

    if(MaxPart < BASENUMBER)
      MaxPart = BASENUMBER;

    P = (lightcone_particle_data *)Mem.myrealloc_movable(P, MaxPart * sizeof(lightcone_particle_data));

    /*
    if(NumPart > MaxPart)  // should be ok because this only happens when P has already been copied away
      Terminate("NumPart=%d > MaxPart=%d", NumPart, MaxPart);
      */
  }

  void reallocate_memory_maxpartsph(int maxpartNew)
  {
    // Don't need to do anything here.
  }

  bool TestIfAboveFillFactor(int SpMaxPart)
  {
    int max_in[2] = {NumPart, SpMaxPart}, max_out[2];
    MPI_Allreduce(max_in, max_out, 2, MPI_INT, MPI_MAX, Communicator);

    /* also recompute the total number of particles in buffer to have current value */
    sumup_large_ints(1, &NumPart, &TotNumPart, Communicator);

    if(max_out[0] > 0 && (All.Ti_Current >= TIMEBASE || max_out[0] >= LIGHTCONE_MAX_FILLFACTOR * max_out[1]))
      return true;
    else
      return false;
  }

  static bool compare_ID(const lightcone_particle_data &a, const lightcone_particle_data &b) { return a.ID.get() < b.ID.get(); }

  static bool compare_ipnest(const lightcone_particle_data &a, const lightcone_particle_data &b) { return a.ipnest < b.ipnest; }

#ifdef REARRANGE_OPTION
  static bool compare_TreeID_ID(const lightcone_particle_data &a, const lightcone_particle_data &b)
  {
    if(a.TreeID < b.TreeID)
      return true;

    if(a.TreeID > b.TreeID)
      return false;

    return a.ID.get() < b.ID.get();
  }
#endif

  void dump_particles(void) {}

  inline int drift_particle(lightcone_particle_data *P, sph_particle_data *SphP, integertime time1, bool ignore_light_cone = false)
  {
    return 0;
  }

  inline MyFloat get_Hsml(int i) { return 0; }

  inline MyFloat get_DtHsml(int i) { return 0; }

  inline MyFloat get_OldAcc(int i) { return 0; }

  inline MyFloat get_Csnd(int i) { return 0; }

  inline double get_utherm_from_entropy(int i) { return 0; }

  inline int getTimeBinSynchronized(int bin) { return 1; }

  void fill_active_gravity_list_with_all_particles(void) {}

#ifdef FOF
  MyIDStorage *MinID;
  int *Head, *Len, *Next, *Tail, *MinIDTask;
  MyFloat *fof_nearest_distance;
  MyFloat *fof_nearest_hsml;
#if defined(LIGHTCONE_PARTICLES_GROUPS)
  double *DistanceOrigin;
#endif

  struct bit_flags
  {
    unsigned char Nonlocal : 2, MinIDChanged : 2, Marked : 2, Changed : 2;
  } * Flags;

  double LinkL;

  void link_two_particles(int target, int j)
  {
    if(Head[target] != Head[j]) /* only if not yet linked */
      {
        int p, s;
        if(Len[Head[target]] > Len[Head[j]]) /* p group is longer */
          {
            p = target;
            s = j;
          }
        else
          {
            p = j;
            s = target;
          }
        Next[Tail[Head[p]]] = Head[s];

        Tail[Head[p]] = Tail[Head[s]];

        Len[Head[p]] += Len[Head[s]];

        if(MinID[Head[s]].get() < MinID[Head[p]].get())
          {
            MinID[Head[p]]     = MinID[Head[s]];
            MinIDTask[Head[p]] = MinIDTask[Head[s]];
          }

#if defined(LIGHTCONE_PARTICLES_GROUPS)
        if(DistanceOrigin[Head[p]] < DistanceOrigin[Head[s]])
          DistanceOrigin[Head[s]] = DistanceOrigin[Head[p]];
        else
          DistanceOrigin[Head[p]] = DistanceOrigin[Head[s]];
#endif

        int ss = Head[s];
        do
          Head[ss] = Head[p];
        while((ss = Next[ss]) >= 0);
      }
  }

#ifdef SUBFIND
  struct nearest_r2_data
  {
    double dist[2];
  } * R2Loc;

#endif
#endif
};

#endif

#endif
