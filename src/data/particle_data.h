/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file particle_data.h
 *
 *  \brief declares a structure that holds the data stored for a single particle
 */

#ifndef PARTDATA_H
#define PARTDATA_H

#include "gadgetconfig.h"

#include <atomic>
#include <climits>

#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/idstorage.h"
#include "../data/intposconvert.h"
#include "../data/macros.h"
#include "../data/mymalloc.h"
#include "../data/symtensors.h"
#include "../mpi_utils/setcomm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/** This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data
{
  // Note that because the atomic_flag and atomic data-types in particle_data have an implicitly deleted copy operator,
  // we are using them encapsulated in a structure defined in dtypes.h that implements these copy and assignment operators explicitly.
  // This is fine because the code logic guarantees that there is no concurrent access when such copies of particle_data happen.

  particle_data() {}

  MyIntPosType IntPos[3];    /**< particle position at its current time, stored as an integer type */
  MyFloat Vel[3];            /**< particle velocity at its current time */
  vector<MyFloat> GravAccel; /**< particle acceleration due to gravity */
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
  MyFloat GravPM[3]; /**< particle acceleration due to long-range PM gravity force */
#endif

  copyable_atomic<integertime> Ti_Current; /**< current time on integer timeline */
  float OldAcc;                            /**< magnitude of old gravitational force. Used in relative opening criterion */
  int GravCost;                            /**< weight factors used for balancing the work-load */

#ifndef LEAN
 private:
  MyDouble Mass; /**< particle mass */
 public:
#endif

  MyIDStorage ID;           // 6-byte
  signed char TimeBinGrav;  // 1-byte
#ifndef LEAN
  signed char TimeBinHydro;
#endif
#if defined(MERGERTREE) && defined(SUBFIND)
  compactrank_t PrevRankInSubhalo;  // 1-byte
  MyHaloNrType PrevSubhaloNr;       // 6-byte
  approxlen PrevSizeOfSubhalo;      // 2-byte
#endif

#ifndef LEAN
 private:
  unsigned char Type; /**< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
 public:
#endif

#ifndef LEAN
  copyable_atomic_flag access;
#endif

#ifdef REARRANGE_OPTION
  unsigned long long TreeID;
#endif

#if NSOFTCLASSES > 1
 private:
  unsigned char
      SofteningClass : 7; /* we use only 7 bits here so that we can stuff 1 bit for ActiveFlag into it in the Tree_Points structure */
 public:
#endif

#if defined(PMGRID) && defined(PLACEHIGHRESREGION)
  unsigned char InsideOutsideFlag : 1;
#endif

#ifdef FORCETEST
  MyFloat GravAccelDirect[3]; /*!< particle acceleration calculated by direct summation */
  MyFloat PotentialDirect;
  MyFloat DistToID1;
#ifdef PMGRID
  MyFloat GravAccelShortRange[3];
  MyFloat PotentialShortRange;
#ifdef PLACEHIGHRESREGION
  MyFloat GravAccelVeryShortRange[3];
  MyFloat PotentialVeryShortRange;
  MyFloat PotentialHPM;
  MyFloat GravAccelHPM[3];
#endif
#endif
#ifdef FORCETEST_FIXEDPARTICLESET
  bool SelectedFlag;
#endif
#endif

#if defined(EVALPOTENTIAL) || defined(OUTPUT_POTENTIAL)
  MyFloat Potential; /**< gravitational potential */
#if defined(PMGRID)
  MyFloat PM_Potential;
#endif
#ifdef EXTERNALGRAVITY
  MyFloat ExtPotential;
#endif
#endif

#ifdef STARFORMATION
  MyFloat StellarAge;  /**< formation time of star particle */
  MyFloat Metallicity; /**< metallicity of gas or star particle */
#endif

  inline unsigned char getType(void)
  {
#ifdef LEAN
    return 1;
#else
    return Type;
#endif
  }

  inline unsigned char getTimeBinHydro(void)
  {
#ifndef LEAN
    return TimeBinHydro;
#else
    return 0;
#endif
  }

  inline void setTimeBinHydro(unsigned char bin)
  {
#ifndef LEAN
    TimeBinHydro = bin;
#endif
  }

  inline void setType(unsigned char type)
  {
#ifndef LEAN
    Type = type;
#endif
  }

  inline float getOldAcc(void) { return OldAcc; }

  inline int getGravCost(void) { return GravCost; }

  inline MyDouble getMass(void)
  {
#ifdef LEAN
    return All.PartMass;
#else
    return Mass;
#endif
  }

  inline void setMass(MyDouble mass)
  {
#ifndef LEAN
    Mass = mass;
#endif
  }

  inline integertime get_Ti_Current(void) { return Ti_Current; }

  inline signed char getTimeBinGrav(void) { return TimeBinGrav; }

  inline unsigned char getSofteningClass(void)
  {
#if NSOFTCLASSES > 1
    return SofteningClass;
#else
    return 0;
#endif
  }

  inline void setSofteningClass(unsigned char softclass)
  {
#if NSOFTCLASSES > 1
    SofteningClass = softclass;
#endif
  }

  inline double getAscale(void) { return All.Time; }

#if defined(LIGHTCONE_PARTICLES_GROUPS)
  inline void setFlagSaveDistance(void) {}
  inline void clearFlagSaveDistance(void) {}

  inline bool getFlagSaveDistance(void) { return true; }
#endif
};

struct subfind_data
{
  MyHaloNrType GroupNr;
#if defined(MERGERTREE)
  MyHaloNrType SubhaloNr;
  approxlen SizeOfSubhalo;
  compactrank_t RankInSubhalo;
#endif
  char DomainFlag;

  int OriginIndex, OriginTask;
  int TargetIndex, TargetTask;

#ifdef SUBFIND
  int SubRankInGr;

#ifndef SUBFIND_HBT
  struct nearest_ngb_data
  {
    location index[2];
    int count;
  };

  nearest_ngb_data nearest;

  int submark;
  int InvIndex;
#endif

#ifndef LEAN
  int Type;
  MyFloat Utherm;
#endif

#ifdef SUBFIND_STORE_LOCAL_DENSITY
  MyFloat SubfindHsml;     // search radius used for SUBFIND dark matter neighborhood
  MyFloat SubfindDensity;  // total matter density
  MyFloat SubfindVelDisp;  // 3D dark matter velocity dispersion
#endif

  union
  {
    struct
    {
      int originindex, origintask;

      union
      {
        MyFloat DM_Density;
        MyFloat DM_Potential;
      } u;

    } s;

    peanokey Key;
  } u;

  union
  {
    MyFloat DM_Hsml;
    MyFloat DM_BindingEnergy;
  } v;
#else
  /* this are fields defined when we have FOF without SUBFIND */
#ifndef LEAN
  int Type;
#endif
  union
  {
    peanokey Key;
  } u;
#endif
};

#ifdef SUBFIND_ORPHAN_TREATMENT
struct idstoredata
{
  int NumPart;
  MyIDType* ID;
};

#endif

#endif
