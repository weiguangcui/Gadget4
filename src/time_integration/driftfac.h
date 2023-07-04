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

#include "gadgetconfig.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/dtypes.h"
#include "../main/main.h"

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

  static double E_of_a(double a) { return sqrt(All.Omega0 / pow(a, 3) + All.OmegaCurvature / pow(a, 2) + All.OmegaLambda); }

  static double E_of_a_diff(double a) { return (-3 * All.Omega0 / pow(a, 4) + -2 * All.OmegaCurvature / pow(a, 3)) / (2 * E_of_a(a)); }

  double linear_growth_factor(double astart, double aend);

  static double get_Omega0_a(double a) { return All.Omega0 / pow(a, 3) * pow(All.Hubble / hubble_function(a), 2); }
  static double get_OmegaLambda_a(double a) { return All.OmegaLambda * pow(All.Hubble / hubble_function(a), 2); }
  static double get_OmegaCurvature_a(double a) { return All.OmegaCurvature / pow(a, 2) * pow(All.Hubble / hubble_function(a), 2); }
  static double get_OmegaMatter_a(double a) { return get_Omega0_a(a); }

  static double hubble_function(double a) { return All.Hubble * E_of_a(a); }

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

  double linear_growth_simple(double a);
  double linear_growth_ode(double a);

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

  static double growth_simple_int(const double a, void *param)
  {
    if(a == 0)
      return 0;
    else
      return 1.0 / pow(a * hubble_function(a) / All.Hubble, 3);
  }

  static int growth_ode_int(double a, const double y[], double dyda[], void *params)
  {
    dyda[0] = y[1];
    dyda[1] = -(3.0 / a + E_of_a_diff(a) / E_of_a(a)) * y[1] + 1.5 * get_OmegaMatter_a(a) / pow(a, 2) * y[0];
    return GSL_SUCCESS;
  }
};

extern driftfac Driftfac;

#endif
