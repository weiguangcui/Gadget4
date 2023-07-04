/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file simparticles.h
 *
 *  \brief class for organizing the storage of the actual simulation particles
 */

#ifndef SIMPART_H
#define SIMPART_H

#include "gadgetconfig.h"

#include <math.h>

#include "../data/allvars.h"
#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/macros.h"
#include "../data/mymalloc.h"
#include "../data/particle_data.h"
#include "../data/sph_particle_data.h"
#include "../main/main.h"
#include "../mpi_utils/mpi_utils.h"
#include "../mpi_utils/setcomm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

#ifdef LIGHTCONE
class lightcone;
#endif

class simparticles : public intposconvert, public setcomm
{
 public:
  simparticles(MPI_Comm comm) : setcomm(comm) {}

  int NumPart; /**< number of particles on the LOCAL processor */
  int NumGas;  /**< number of gas particles on the LOCAL processor  */

  int MaxPart;    /**< This gives the maxmimum number of particles that can be stored on one processor. */
  int MaxPartSph; /**< This gives the maxmimum number of SPH particles that can be stored on one processor. */

  long long TotNumPart; /**<  total particle numbers (global value) */
  long long TotNumGas;  /**<  total gas particle number (global value) */

  typedef particle_data pdata;

  /*! This structure holds all the information that is
   * stored for each particle of the simulation.
   */
  particle_data *P; /*!< holds particle data on local processor */

  /* the following struture holds data that is stored for each SPH particle in addition to the collisionless
   * variables.
   */
  sph_particle_data *SphP; /*!< holds SPH particle data on local processor */

  unsigned short int MarkerValue;

  subfind_data *PS;

  inline void copy_particle(particle_data *Ptarget, particle_data *Psource)
  {
    // we do this ugly trick here because the atomic_flag in particle_data has an implicitly deleted copy operator...
    // but we know what we are doing, and this is the easiest way at the moment to work around this in our case unnecessary protection
    memcpy(static_cast<void *>(Ptarget), static_cast<void *>(Psource), sizeof(particle_data));
  }

  static bool inline compare_IDs(const MyIDType &a, const MyIDType &b) { return a < b; }

#if defined(LIGHTCONE_PARTICLES_GROUPS) && defined(FOF)
  double *DistanceOrigin;
#endif

#ifdef SUBFIND_ORPHAN_TREATMENT
  idstoredata IdStore;
  static inline bool compare_SpP_ID(const particle_data &a, const particle_data &b) { return a.ID.get() < b.ID.get(); }
#endif

#ifdef LIGHTCONE
  lightcone *LightCone;
#endif

#ifdef FOF
  MyIDStorage *MinID;
  int *Len;  // this is here enough in 32bit because only the group segments on the local processor are treated
  int *Head, *Next, *Tail, *MinIDTask;
  MyFloat *fof_nearest_distance;
  MyFloat *fof_nearest_hsml;

  struct bit_flags
  {
    unsigned char Nonlocal : 2, MinIDChanged : 2, Marked : 2;
  } * Flags;

  double LinkL;

  inline void link_two_particles(int target, int j)
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

#ifdef PMGRID
  double Asmth[2], Rcut[2];
#endif

#if defined(PMGRID) && (!defined(PERIODIC) || defined(PLACEHIGHRESREGION))
  double TotalMeshSize[2]; /* this is in integer space but should be double here to protect against overflows */
  MySignedIntPosType Corner[2][3];
  MySignedIntPosType Xmintot[2][3], Xmaxtot[2][3];
  MyIntPosType MeshSize[2][3];
  MyIntPosType Left[2][3];
  MyIntPosType OldMeshSize[2];
  MyIntPosType ReferenceIntPos[2][3];
#endif

#if defined(RANDOMIZE_DOMAINCENTER_TYPES) || defined(PLACEHIGHRESREGION)
  MyIntPosType PlacingMask;
  MyIntPosType PlacingBlocksize;
#endif

#ifdef PLACEHIGHRESREGION
  inline int check_high_res_overlap(MyIntPosType *center, MyIntPosType halflen)
  {
    MyIntPosType intleft[3] = {center[0] - halflen - ReferenceIntPos[HIGH_MESH][0],
                               center[1] - halflen - ReferenceIntPos[HIGH_MESH][1],
                               center[2] - halflen - ReferenceIntPos[HIGH_MESH][2]};

    MyIntPosType intright[3] = {center[0] + halflen - ReferenceIntPos[HIGH_MESH][0],
                                center[1] + halflen - ReferenceIntPos[HIGH_MESH][1],
                                center[2] + halflen - ReferenceIntPos[HIGH_MESH][2]};

    MySignedIntPosType *left  = (MySignedIntPosType *)intleft;
    MySignedIntPosType *right = (MySignedIntPosType *)intright;

    if(right[0] <= Xmintot[HIGH_MESH][0] || left[0] >= Xmaxtot[HIGH_MESH][0] || right[1] <= Xmintot[HIGH_MESH][1] ||
       left[1] >= Xmaxtot[HIGH_MESH][1] || right[2] <= Xmintot[HIGH_MESH][2] || left[2] >= Xmaxtot[HIGH_MESH][2])
      return FLAG_OUTSIDE;
    else if(right[0] <= Xmaxtot[HIGH_MESH][0] && left[0] >= Xmintot[HIGH_MESH][0] && right[1] <= Xmaxtot[HIGH_MESH][1] &&
            left[1] >= Xmintot[HIGH_MESH][1] && right[2] <= Xmaxtot[HIGH_MESH][2] && left[2] >= Xmintot[HIGH_MESH][2])
      return FLAG_INSIDE;
    else
      return FLAG_BOUNDARYOVERLAP;
  }

  inline int check_high_res_point_location(MyIntPosType *intpos)
  {
    MyIntPosType relpos[3] = {intpos[0] - ReferenceIntPos[HIGH_MESH][0], intpos[1] - ReferenceIntPos[HIGH_MESH][1],
                              intpos[2] - ReferenceIntPos[HIGH_MESH][2]};

    MySignedIntPosType *pos = (MySignedIntPosType *)relpos;

    if(pos[0] < Xmintot[HIGH_MESH][0] || pos[0] >= Xmaxtot[HIGH_MESH][0] || pos[1] < Xmintot[HIGH_MESH][1] ||
       pos[1] >= Xmaxtot[HIGH_MESH][1] || pos[2] < Xmintot[HIGH_MESH][2] || pos[2] >= Xmaxtot[HIGH_MESH][2])
      return FLAG_OUTSIDE;
    else
      return FLAG_INSIDE;
  }

#endif

  int TimeBinSynchronized[TIMEBINS];
  TimeBinData TimeBinsHydro;
  TimeBinData TimeBinsGravity;

  int nsource;
  int *indexlist;

#ifdef STARFORMATION
  double TimeBinSfr[TIMEBINS];
#endif

  inline int getTimeBinSynchronized(int bin) { return TimeBinSynchronized[bin]; }

#ifdef REARRANGE_OPTION
  static bool compare_TreeID_ID(const particle_data &a, const particle_data &b)
  {
    if(a.TreeID < b.TreeID)
      return true;

    if(a.TreeID > b.TreeID)
      return false;

    return a.ID.get() < b.ID.get();
  }

  static bool compare_ID(const particle_data &a, const particle_data &b) { return a.ID.get() < b.ID.get(); }
#endif

  inline MyFloat get_DtHsml(int i) { return SphP[i].DtHsml; }

  inline MyFloat get_Csnd(int i) { return SphP[i].Csnd; }

  inline MyFloat get_OldAcc(int i) { return P[i].OldAcc; }

  /* sets the internal energy per unit mass of particle i  from its entropy */
  inline double get_utherm_from_entropy(int i)
  {
#ifdef ISOTHERM_EQS
    return SphP[i].Entropy;
#else
    double fact_entropy_to_u = pow(SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
    return SphP[i].Entropy * fact_entropy_to_u;
#endif
  }

  /* sets the entropy of particle i from its internal energy per unit mass */
  inline void set_entropy_from_utherm(double utherm, int i)
  {
    double fact_u_to_entropy = GAMMA_MINUS1 / pow(SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1);
    SphP[i].Entropy          = utherm * fact_u_to_entropy;
    SphP[i].EntropyPred      = SphP[i].Entropy;

#ifdef PRESSURE_ENTROPY_SPH
    SphP[i].EntropyToInvGammaPred = pow(SphP[i].EntropyPred, 1.0 / GAMMA);
#endif
  }

  void fill_active_gravity_list_with_all_particles(void)
  {
    TimeBinsGravity.NActiveParticles = 0;

    for(int i = 0; i < NumPart; i++)
      TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles++] = i;
  }

  /* This routine allocates memory for
   * particle storage, both the collisionless and the SPH particles.
   * The memory for the ordered binary tree of the timeline
   * is also allocated.
   */
  void allocate_memory(void)
  {
    /* Note: P and SphP are initialized to zero */
    P    = (particle_data *)Mem.mymalloc_movable_clear(&P, "P", MaxPart * sizeof(particle_data));
    SphP = (sph_particle_data *)Mem.mymalloc_movable_clear(&SphP, "SphP", MaxPartSph * sizeof(sph_particle_data));

    TimeBinsHydro.timebins_allocate();
    TimeBinsGravity.timebins_allocate();
  }

  void free_memory(void)
  {
    TimeBinsGravity.timebins_free();
    TimeBinsHydro.timebins_free();

    Mem.myfree(SphP);
    Mem.myfree(P);
  }

  void reallocate_memory_maxpart(int maxpartNew)
  {
    mpi_printf("ALLOCATE: Changing to MaxPart = %d\n", maxpartNew);

    P = (particle_data *)Mem.myrealloc_movable(P, maxpartNew * sizeof(particle_data));
    if(maxpartNew > MaxPart)
      memset(((char *)P) + MaxPart * sizeof(particle_data), 0, (maxpartNew - MaxPart) * sizeof(particle_data));
    MaxPart = maxpartNew;

    TimeBinsGravity.timebins_reallocate();
  }

  void reallocate_memory_maxpartsph(int maxpartsphNew)
  {
    mpi_printf("ALLOCATE: Changing to MaxPartSph = %d\n", maxpartsphNew);

    SphP = (sph_particle_data *)Mem.myrealloc_movable(SphP, maxpartsphNew * sizeof(sph_particle_data));
    if(maxpartsphNew > MaxPartSph)
      memset(((char *)SphP) + MaxPartSph * sizeof(sph_particle_data), 0, (maxpartsphNew - MaxPartSph) * sizeof(sph_particle_data));
    MaxPartSph = maxpartsphNew;

    TimeBinsHydro.timebins_reallocate();
  }

  /*! This function dumps some of the basic particle data to a file. In case
   *  the tree construction fails, this is called just before the run
   *  terminates with an error message. Examination of the generated file may
   *  then give clues to what caused the problem.
   */
  void dump_particles(void)
  {
    FILE *fd;
    char buffer[MAXLEN_PATH];
    snprintf(buffer, MAXLEN_PATH, "particles_%d.dat", ThisTask);
    if((fd = fopen(buffer, "w")))
      {
        fwrite(&NumPart, 1, sizeof(int), fd);
        for(int i = 0; i < NumPart; i++)
          fwrite(&P[i].IntPos[0], 3, sizeof(MyIntPosType), fd);
        for(int i = 0; i < NumPart; i++)
          fwrite(&P[i].Vel[0], 3, sizeof(MyFloat), fd);
        for(int i = 0; i < NumPart; i++)
          fwrite(&P[i].ID, 1, sizeof(MyIDStorage), fd);
        fclose(fd);
      }
  }

  /** \brief Print information relative to a particle / cell to standard output.
   *
   *  \param i particle / cell index
   */
  void print_particle_info(int i)
  {
    MyReal pos[3];
    intpos_to_pos(P[i].IntPos, pos); /* converts the integer coordinates to floating point */

    printf("Task=%d, ID=%llu, Type=%d, TimeBinGrav=%d, TimeBinHydro=%d, Mass=%g, pos=%g|%g|%g, vel=%g|%g|%g, OldAcc=%g\n", ThisTask,
           (unsigned long long)P[i].ID.get(), P[i].getType(), P[i].TimeBinGrav, P[i].getTimeBinHydro(), P[i].getMass(), pos[0], pos[1],
           pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].OldAcc);
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
    printf("GravAccel=%g|%g|%g, GravPM=%g|%g|%g, Soft=%g, SoftClass=%d\n", P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2],
           P[i].GravPM[0], P[i].GravPM[1], P[i].GravPM[2], All.ForceSoftening[P[i].getSofteningClass()], P[i].getSofteningClass());
#else
#ifndef LEAN
    printf("GravAccel=%g|%g|%g, Soft=%g, SoftType=%d\n", P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2],
           All.ForceSoftening[P[i].getSofteningClass()], P[i].getSofteningClass());
#endif
#endif

    if(P[i].getType() == 0)
      {
        printf("rho=%g, hsml=%g, entr=%g, csnd=%g\n", SphP[i].Density, SphP[i].Hsml, SphP[i].Entropy, SphP[i].get_sound_speed());
        printf("ID=%llu SphP[p].CurrentMaxTiStep=%g\n", (unsigned long long)P[i].ID.get(), SphP[i].CurrentMaxTiStep);
      }

    myflush(stdout);
  }

  /** \brief Print information relative to a particle / cell to standard output given its ID.
   *  *
   *   *  \param ID particle / cell ID
   *    */
  void print_particle_info_from_ID(MyIDType ID)
  {
    for(int i = 0; i < NumPart; i++)
      if(P[i].ID.get() == ID)
        print_particle_info(i);
  }

 public:
  inline int get_active_index(int idx)
  {
#ifdef HIERARCHICAL_GRAVITY
    return TimeBinsGravity.ActiveParticleList[idx];
#else
    return idx;
#endif
  }

  void reconstruct_timebins(void);
  integertime find_next_sync_point(void);
  void mark_active_timebins(void);
  void drift_all_particles(void);
  int drift_particle(particle_data *P, sph_particle_data *SphP, integertime time1, bool ignore_light_cone = false);
  void make_list_of_active_particles(void);
  integertime get_timestep_grav(int p);
  integertime get_timestep_hydro(int p);

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  integertime get_timestep_pm(void);
#endif

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
  void find_long_range_step_constraint(void);
#endif

  void timebins_get_bin_and_do_validity_checks(integertime ti_step, int *bin_new, int bin_old);

  void assign_hydro_timesteps(void);
  void timebin_cleanup_list_of_active_particles(void);

  int test_if_grav_timestep_is_too_large(int p, int bin);
  int get_timestep_bin(integertime ti_step);

#ifdef ADAPTIVE_HYDRO_SOFTENING
  int get_softeningtype_for_hydro_particle(int i)
  {
    double soft = All.GasSoftFactor * SphP[i].Hsml;

    if(soft <= All.ForceSoftening[NSOFTCLASSES])
      return NSOFTCLASSES;

    int k = 0.5 + log(soft / All.ForceSoftening[NSOFTCLASSES]) / log(All.AdaptiveHydroSofteningSpacing);
    if(k >= NSOFTCLASSES_HYDRO)
      k = NSOFTCLASSES_HYDRO - 1;

    return NSOFTCLASSES + k;
  }
#endif

#ifdef INDIVIDUAL_GRAVITY_SOFTENING

 public:
#if(INDIVIDUAL_GRAVITY_SOFTENING) & 2
#error "INDIVIDUAL_GRAVITY_SOFTENING may not include particle type 1 which is used as a reference point"
#endif

#if((INDIVIDUAL_GRAVITY_SOFTENING)&1) && defined(ADAPTIVE_HYDRO_SOFTENING)
#error "INDIVIDUAL_GRAVITY_SOFTENING may not include particle type 0 when ADAPTIVE_HYDRO_SOFTENING is used"
#endif

  int get_softening_type_from_mass(double mass)
  {
    int min_type   = -1;
    double eps     = get_desired_softening_from_mass(mass);
    double min_dln = MAX_FLOAT_NUMBER;

    for(int i = 0; i < NSOFTCLASSES; i++)
      {
        if(All.ForceSoftening[i] > 0)
          {
            double dln = fabs(log(eps) - log(All.ForceSoftening[i]));

            if(dln < min_dln)
              {
                min_dln  = dln;
                min_type = i;
              }
          }
      }

    if(min_type < 0)
      Terminate("min_type < 0");

    return min_type;
  }

  /*! \brief Returns the desired softening length depending on the particle mass with type 1 as a reference point
   *
   * \param mass particle mass
   * \return softening length for a particle of mass #mass
   */
  double get_desired_softening_from_mass(double mass)
  {
    return All.ForceSoftening[All.SofteningClassOfPartType[1]] * pow(mass / All.AvgType1Mass, 1.0 / 3);
  }

  /*! \brief Initializes the mass dependent softening calculation for Type 1 particles
   *
   * The average mass of Type 1 particles is calculated.
   */
  void init_individual_softenings(void)
  {
    int ndm     = 0;
    double mass = 0, masstot, massmin = MAX_DOUBLE_NUMBER, massmax = 0;
    long long ndmtot;

    for(int i = 0; i < NumPart; i++)
      if(P[i].getType() == 1)
        {
          ndm++;
          mass += P[i].getMass();

          if(massmin > P[i].getMass())
            massmin = P[i].getMass();

          if(massmax < P[i].getMass())
            massmax = P[i].getMass();
        }

    sumup_large_ints(1, &ndm, &ndmtot, Communicator);
    MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, Communicator);

    MPI_Allreduce(MPI_IN_PLACE, &massmin, 1, MPI_DOUBLE, MPI_MIN, Communicator);
    MPI_Allreduce(MPI_IN_PLACE, &massmax, 1, MPI_DOUBLE, MPI_MAX, Communicator);

    All.AvgType1Mass = masstot / ndmtot;

    mpi_printf("INIT: AvgType1Mass = %g   (min=%g max=%g) Ndm1tot=%lld\n", All.AvgType1Mass, massmin, massmax, ndmtot);

    if(massmax > 1.00001 * massmin)
      Terminate("Strange: Should use constant mass type-1 particles if INDIVIDUAL_GRAVITY_SOFTENING is used\n");

    if(All.ComovingIntegrationOn)
      {
        double rhomean_dm = (All.Omega0 - All.OmegaBaryon) * (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

        mpi_printf("INIT: For this AvgType1Mass, the mean particle spacing is %g and the assigned softening is %g\n",
                   pow(All.AvgType1Mass / rhomean_dm, 1.0 / 3), All.SofteningTable[All.SofteningClassOfPartType[1]]);
      }

    for(int i = 0; i < NumPart; i++)
      if(((1 << P[i].getType()) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[i].setSofteningClass(get_softening_type_from_mass(P[i].getMass()));
  }
#endif
};

#endif
