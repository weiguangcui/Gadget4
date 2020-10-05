/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file grav_forcetest.h
 *
 *  \brief declares a class needed for the force test calculations
 */

#ifndef GRAV_FORCETEST_H
#define GRAV_FORCETEST_H

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/simparticles.h"
#include "../domain/domain.h"
#include "../gravtree/gravtree.h"

#define TESTGRID 384

#ifdef LONG_X_BITS
#if TESTGRID != ((TESTGRID / LONG_X) * LONG_X)
#error "TESTGRID must be a multiple of the stretch factor in the x-direction"
#endif
#endif

#ifdef LONG_Y_BITS
#if TESTGRID != ((TESTGRID / LONG_Y) * LONG_Y)
#error "TESTGRID must be a multiple of the stretch factor in the y-direction"
#endif
#endif

#ifdef LONG_Z_BITS
#if TESTGRID != ((TESTGRID / LONG_Z) * LONG_Z)
#error "TESTGRID must be a multiple of the stretch factor in the x-direction"
#endif
#endif

#define TESTGRIDX ((TESTGRID / LONG_X))
#define TESTGRIDY ((TESTGRID / LONG_Y))
#define TESTGRIDZ ((TESTGRID / LONG_Z))

#define TESTGRIDZ2 (TESTGRIDZ / 2 + 1)

#define ASMTH_DIRECT 30.0

class gravtest
{
 private:
  simparticles *Sp;
  gravtree<simparticles> *GravTree;
  domain<simparticles> *D;

 public:
  gravtest(simparticles *Sp_ptr, gravtree<simparticles> *GravTree_ptr, domain<simparticles> *Domain_ptr)
  {
    Sp       = Sp_ptr;
    GravTree = GravTree_ptr;
    D        = Domain_ptr;
  }

  void gravity_forcetest(int timebin);
};

#endif /* GRAV_FORCETEST_H */
