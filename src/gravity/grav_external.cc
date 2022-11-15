/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file gravity_external.cc
 *
 * \brief can add an optional external gravity field
 */

#include "gadgetconfig.h"

#ifdef EXTERNALGRAVITY

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fmm/fmm.h"
#include "../gravtree/gravtree.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

void sim::gravity_external(void)
{
#ifdef PERIODIC
  // pick middle of (stretched) box, but could also choose other point
  vector<double> pos_center{0.5 * All.BoxSize / LONG_X, 0.5 * All.BoxSize / LONG_Y, 0.5 * All.BoxSize / LONG_Z};
#else
  // here pick origin
  vector<double> pos_center{0, 0, 0};
#endif
  MyIntPosType intpos_center[3];
  Sp.pos_to_intpos(pos_center.da, intpos_center);

  for(int i = 0; i < Sp.TimeBinsGravity.NActiveParticles; i++)
    {
      int target = Sp.TimeBinsGravity.ActiveParticleList[i];

#if defined(EVALPOTENTIAL) || defined(OUTPUT_POTENTIAL)
      Sp.P[target].ExtPotential = 0;
#endif

#ifdef EXTERNALGRAVITY_STATICHQ
      {
        vector<double> pos;
        Sp.nearest_image_intpos_to_pos(Sp.P[target].IntPos, intpos_center, pos.da);

        double r = sqrt(pos.r2());

        double m = All.Mass_StaticHQHalo * pow(r / (r + All.A_StaticHQHalo), 2);

        if(r > 0)
          Sp.P[target].GravAccel += (-All.G * m / (r * r * r)) * pos;

#if defined(EVALPOTENTIAL) || defined(OUTPUT_POTENTIAL)
        Sp.P[target].ExtPotential += (-All.G * All.Mass_StaticHQHalo / (r + All.A_StaticHQHalo));
#endif
      }
#endif
    }
}

#endif
