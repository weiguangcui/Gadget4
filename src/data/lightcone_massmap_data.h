/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file lightcone_massmap_data.h
 *
 *  \brief defines a structure to hold data for lightcone particles projected on healpix map
 */

#ifndef MMPARTDATA_H
#define MMPARTDATA_H

#if defined(LIGHTCONE) && defined(LIGHTCONE_MASSMAPS)

#include "gadgetconfig.h"

#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/macros.h"

struct lightcone_massmap_data
{
  MyFloat Ascale;
  int PixIndex; /* global healpix index */
  int Task;

#ifndef LEAN
 private:
  MyDouble Mass;

 public:
#endif

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
};

#endif
#endif
