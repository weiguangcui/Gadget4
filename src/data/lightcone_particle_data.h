/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file lightcone_particle_data.h
 *
 *  \brief declares a structure for holding particle on the lightcone
 */

#ifndef LCPARTDATA_H
#define LCPARTDATA_H

#include "gadgetconfig.h"

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES)

#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/idstorage.h"
#include "../data/macros.h"
#include "../data/symtensors.h"

struct lightcone_particle_data
{
  MyIntPosType IntPos[3];
  MyFloat Vel[3];
#ifdef LIGHTCONE_OUTPUT_ACCELERATIONS
  vector<MyFloat> GravAccel;
#endif
  MyFloat Ascale;
  MyIDStorage ID;

  long ipnest;

#ifdef REARRANGE_OPTION
  unsigned long long TreeID;
#endif

#ifndef LEAN
 private:
  MyDouble Mass; /**< particle mass */
 public:
#endif

#ifndef LEAN
 private:
  unsigned char Type; /**< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
 public:
#endif

#if NSOFTCLASSES > 1
 private:
  unsigned char
      SofteningClass : 7; /* we use only 7 bits here so that we can stuff 1 bit for ActiveFlag into it in the Tree_Points structure */
 public:
#endif

#if defined(MERGERTREE) && defined(SUBFIND)
  // compactrank_t PrevRankInSubhalo;  // 1-byte
  MyHaloNrType PrevSubhaloNr;   // 6-byte
  approxlen PrevSizeOfSubhalo;  // 2-byte
#endif

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  unsigned char OriginIndex;
#endif

#ifdef LIGHTCONE_IMAGE_COMP_HSML_VELDISP
  int NumNgb;
  MyFloat Hsml;
  MyFloat VelDisp;
  MyFloat Density;
  MyFloat Vx;
  MyFloat Vy;
  MyFloat Vz;
#endif

#if defined(LIGHTCONE_PARTICLES_GROUPS) && defined(FOF)
 private:
  bool FlagSaveDistance;

 public:
  inline void setFlagSaveDistance(void) { FlagSaveDistance = true; }
  inline void clearFlagSaveDistance(void) { FlagSaveDistance = false; }
  inline bool getFlagSaveDistance(void) { return FlagSaveDistance; }

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

#endif

  inline unsigned char getType(void)
  {
#ifdef LEAN
    return 1;
#else
    return Type;
#endif
  }

  inline double getAscale(void) { return Ascale; }

  inline void setType(unsigned char type)
  {
#ifndef LEAN
    Type = type;
#endif
  }

  inline MyDouble getMass(void)
  {
#ifdef LEAN
    return All.PartMass;
#else
    return Mass;
#endif
  }

  inline integertime get_Ti_Current(void) { return 0; }

  inline void setMass(MyDouble mass)
  {
#ifndef LEAN
    Mass = mass;
#endif
  }

  inline float getOldAcc(void) { return 0; }

  inline signed char getTimeBinGrav(void) { return 0; }
  inline signed char getTimeBinHydro(void) { return 0; }
  inline int getGravCost(void) { return 0; }

#if defined(MERGERTREE) && defined(SUBFIND)
  inline void setPrevSubhaloNr(int nr) {}
  inline void setPrevRankInSubhalo(int nr) {}
  inline long long getPrevSubhaloNr(void) { return 0; }
  inline int getPrevRankInSubhalo(void) { return 0; }
#endif
};

#endif
#endif
