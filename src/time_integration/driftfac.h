/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  driftfac.h
 *
 *  \brief declares a class for supporting cosmological drift/kick factor calculations
 */

#ifndef DRIFTFAC_H
#define DRIFTFAC_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/dtypes.h"
#include "../main/main.h"
#include "gadgetconfig.h"

class driftfac
{
 public:
  void init_drift_table(void);
  double get_drift_factor(integertime time0, integertime time1);
  double get_gravkick_factor(integertime time0, integertime time1);
  double get_hydrokick_factor(integertime time0, integertime time1);
  double get_comoving_distance(integertime time0);
  double get_comoving_distance_for_scalefactor(double ascale);
  double get_scalefactor_for_comoving_distance(double dist);
  integertime get_gravkick_factor_inverse(double fac);

  static double hubble_function(double a)
  {
    double hubble_a = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a) + All.OmegaLambda;

    hubble_a = All.Hubble * sqrt(hubble_a);

    return hubble_a;
  }

 private:
#define DRIFT_TABLE_LENGTH 1000

  /** table for the cosmological drift factors */
  double DriftTable[DRIFT_TABLE_LENGTH];

  /** table for the cosmological kick factor for gravitational forces */
  double GravKickTable[DRIFT_TABLE_LENGTH];

  /** table for the cosmological kick factor for hydrodynmical forces */
  double HydroKickTable[DRIFT_TABLE_LENGTH];

  double logTimeBegin;
  double logTimeMax;

  static double drift_integ(double a, void *param)
  {
    double h = hubble_function(a);

    return 1 / (h * a * a * a);
  }

  static double gravkick_integ(double a, void *param)
  {
    double h = hubble_function(a);

    return 1 / (h * a * a);
  }

  static double hydrokick_integ(double a, void *param)
  {
    double h = hubble_function(a);

    return 1 / (h * pow(a, 3 * GAMMA_MINUS1) * a);
  }
};

extern driftfac Driftfac;

#endif
