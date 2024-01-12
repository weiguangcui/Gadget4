/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  driftfac.cc
 *
 *  \brief tabulates cosmological drift/kick factors for fast look-up
 */

#include "gadgetconfig.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../time_integration/driftfac.h"

void driftfac::init_drift_table(void)
{
#define WORKSIZE 100000

  gsl_function F;
  gsl_integration_workspace *workspace;

  logTimeBegin = log(All.TimeBegin);
  logTimeMax   = log(All.TimeMax);

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  for(int i = 0; i < DRIFT_TABLE_LENGTH; i++)
    {
      double result, abserr;

      F.function = &driftfac::drift_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
                          1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      DriftTable[i] = result;

      F.function = &gravkick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
                          1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      GravKickTable[i] = result;

      F.function = &hydrokick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), 0,
                          1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      HydroKickTable[i] = result;
    }

  gsl_integration_workspace_free(workspace);
}

/*! This function integrates the cosmological prefactor for a drift
 *   step between time0 and time1. The value returned is
 *  \f[ \int_{a_0}^{a_1} \frac{{\rm d}a}{a^3 * H(a)}
 *  \f]
 *
 *  A lookup-table is used for reasons of speed.
 */
double driftfac::get_drift_factor(integertime time0, integertime time1)
{
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  double a1 = logTimeBegin + time0 * All.Timebase_interval;
  double a2 = logTimeBegin + time1 * All.Timebase_interval;

  double u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  int i1    = (int)u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  double df1;
  if(i1 <= 1)
    df1 = u1 * DriftTable[0];
  else
    df1 = DriftTable[i1 - 1] + (DriftTable[i1] - DriftTable[i1 - 1]) * (u1 - i1);

  double u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  int i2    = (int)u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  double df2;
  if(i2 <= 1)
    df2 = u2 * DriftTable[0];
  else
    df2 = DriftTable[i2 - 1] + (DriftTable[i2] - DriftTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}

double driftfac::get_gravkick_factor(integertime time0, integertime time1)
{
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  double a1 = logTimeBegin + time0 * All.Timebase_interval;
  double a2 = logTimeBegin + time1 * All.Timebase_interval;

  double u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  int i1    = (int)u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  double df1;
  if(i1 <= 1)
    df1 = u1 * GravKickTable[0];
  else
    df1 = GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1);

  double u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  int i2    = (int)u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  double df2;
  if(i2 <= 1)
    df2 = u2 * GravKickTable[0];
  else
    df2 = GravKickTable[i2 - 1] + (GravKickTable[i2] - GravKickTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}

double driftfac::get_hydrokick_factor(integertime time0, integertime time1)
{
  static integertime last_time0 = -1, last_time1 = -1;
  static double last_value;

  if(time0 == last_time0 && time1 == last_time1)
    return last_value;

  /* note: will only be called for cosmological integration */

  double a1 = logTimeBegin + time0 * All.Timebase_interval;
  double a2 = logTimeBegin + time1 * All.Timebase_interval;

  double u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  int i1    = (int)u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  double df1;
  if(i1 <= 1)
    df1 = u1 * HydroKickTable[0];
  else
    df1 = HydroKickTable[i1 - 1] + (HydroKickTable[i1] - HydroKickTable[i1 - 1]) * (u1 - i1);

  double u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  int i2    = (int)u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  double df2;
  if(i2 <= 1)
    df2 = u2 * HydroKickTable[0];
  else
    df2 = HydroKickTable[i2 - 1] + (HydroKickTable[i2] - HydroKickTable[i2 - 1]) * (u2 - i2);

  last_time0 = time0;
  last_time1 = time1;

  return last_value = (df2 - df1);
}

double driftfac::get_comoving_distance(integertime time0)
{
  /* we just need to multiply this with the speed of light to get the correct cosmological factor */
  double fac = get_gravkick_factor(time0, TIMEBASE);

  return fac * (CLIGHT / All.UnitVelocity_in_cm_per_s);
}

double driftfac::get_comoving_distance_for_scalefactor(double ascale)
{
  integertime time0 = log(ascale / All.TimeBegin) / All.Timebase_interval;

  /* we just need to multiply this with the speed of light to get the correct cosmological factor */
  double fac = get_gravkick_factor(time0, TIMEBASE);

  return fac * (CLIGHT / All.UnitVelocity_in_cm_per_s);
}

double driftfac::get_scalefactor_for_comoving_distance(double dist)
{
  double fac = dist / (CLIGHT / All.UnitVelocity_in_cm_per_s);

  integertime time0 = get_gravkick_factor_inverse(fac);

  double ascale = All.TimeBegin * exp(time0 * All.Timebase_interval);

  return ascale;
}

integertime driftfac::get_gravkick_factor_inverse(double fac)
{
  integertime time1 = TIMEBASE;

  double a2 = logTimeBegin + time1 * All.Timebase_interval;

  double u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  int i2    = (int)u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  double df2;
  if(i2 <= 1)
    df2 = u2 * GravKickTable[0];
  else
    df2 = GravKickTable[i2 - 1] + (GravKickTable[i2] - GravKickTable[i2 - 1]) * (u2 - i2);

  double df1 = df2 - fac;

  int i0 = 0;
  int i1 = DRIFT_TABLE_LENGTH - 1;

  if(df1 < 0 || df1 > GravKickTable[i1])
    Terminate("out of range:  df1=%g  GravKickTable[i0]=%g  GravKickTable[i1]=%g\n", df1, GravKickTable[i0], GravKickTable[i1]);

  double u1;

  if(df1 <= GravKickTable[0])
    u1 = df1 / GravKickTable[0];
  else
    {
      while(i1 - i0 > 1)
        {
          int im = (i0 + i1) / 2;
          if(df1 < GravKickTable[im])
            i1 = im;
          else
            i0 = im;
        }

      u1 = (df1 - GravKickTable[i0]) / (GravKickTable[i1] - GravKickTable[i0]) + i1;
    }

  double a1 = u1 * (logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH + logTimeBegin;

  integertime time0 = (a1 - logTimeBegin) / All.Timebase_interval;

  return time0;
}

double driftfac::linear_growth_factor(const double astart, const double aend)
{
  // if we have a simple cosmology without radiation, with matter, a cosmological constant, and flat space,
  // we can apply the simple integral solution:

  ///  return linear_growth_simple(aend) / linear_growth_simple(astart);

  // but in the more general case, we need to integrate the ODE directly to find the correct
  // growth factor.

  return linear_growth_ode(aend) / linear_growth_ode(astart);
}

double driftfac::linear_growth_ode(const double a)
{
  gsl_odeiv2_system sys = {growth_ode_int, NULL, 2, NULL};
  gsl_odeiv2_driver *d  = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

  double amin = 0.00001;  // start time

  double epsilon = 3 + 2 * amin * E_of_a_diff(amin) / E_of_a(amin);
  double w       = 1 - get_OmegaMatter_a(amin);
  double n       = 0.25 * (-1 - epsilon + sqrt(pow(1 + epsilon, 2) + 24 * (1 - w)));

  // ICs
  double D[2]     = {1.0, n / amin};
  double acurrent = amin;

  int status = gsl_odeiv2_driver_apply(d, &acurrent, a, D);

  if(status != GSL_SUCCESS)
    Terminate("error, ODE return value=%d\n", status);

  gsl_odeiv2_driver_free(d);

  return D[0];
}

double driftfac::linear_growth_simple(const double a)
{
  /* this simple integral solution is only correct for LCDM variants, i.e.
   * no radition, flat space, and simple cosmological constant
   */
  const int worksize = WORKSIZE;

  double result, abserr;
  gsl_function F;

  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(worksize);
  F.function                           = &growth_simple_int;

  gsl_integration_qag(&F, 0, a, 0, 1.0e-8, worksize, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  const double hubble_a = hubble_function(a) / All.Hubble;

  return hubble_a * result;
}
