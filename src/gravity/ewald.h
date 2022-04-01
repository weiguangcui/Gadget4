/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file ewald.h
 *
 *  \brief definition of class that implements the functionality needed for Ewald corrections
 */

#ifndef EWALD_H
#define EWALD_H

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../data/symtensors.h"
#include "../gravity/ewaldtensors.h"
#include "../gravtree/gravtree.h"
#include "../io/io.h"
#include "../system/system.h"

/*!
 *
 *  This file contains the definitions of the Ewald correction lookup table code.
 */

#if defined(GRAVITY_TALLBOX) && defined(FMM)
#error "FFM and GRAVITY_TALLBOX cannot be used together"
#endif

#define EWLEVEL 6

#define EN (1 << EWLEVEL)  //!<  Base dimension of cubical Ewald lookup table, for one octant

#define ENX (DBX * EN)
#define ENY (DBY * EN)
#define ENZ (DBZ * EN)

/*!  \brief Ewald correction functionality.
 *
 * This class collects all the functionality provided for Ewald table lookups.
 */
class ewald : public intposconvert, public io_streamcount, public setcomm
{
 public:
  ewald(void) : setcomm("delayed init") {}

  void ewald_init(void);

  struct ewald_header
  {
    int resx, resy, resz, varsize, ewaldtype;
  };

  enum interpolate_options
  {
    POINTMASS,
    MULTIPOLES,
    LONGRANGE_FORCETEST
  };

  void ewald_corr(double dx, double dy, double dz, enum interpolate_options, ewald_data &fper);
  void ewald_corr_exact(double dx, double dy, double dz, enum interpolate_options flag, ewald_data &fper);
  void ewald_gridlookup(const MyIntPosType *p_intpos, const MyIntPosType *target_intpos, enum interpolate_options flag,
                        ewald_data &fper);

  double ewald_D0(double x, double y, double z);
  vector<double> ewald_D1(double x, double y, double z);
  symtensor2<double> ewald_D2(double x, double y, double z);
  symtensor3<double> ewald_D3(double x, double y, double z);
  symtensor4<double> ewald_D4(double x, double y, double z);
  symtensor5<double> ewald_D5(double x, double y, double z);
  symtensor6<double> ewald_D6(double x, double y, double z);
  symtensor7<double> ewald_D7(double x, double y, double z);

  ewaldtensor6<double> ewald_P6(void);
  ewaldtensor8<double> ewald_P8(void);
  ewaldtensor10<double> ewald_P10(void);

 private:
  ewald_data *Ewd;  // points to an [ENX + 1][ENY + 1][ENZ + 1] array

  inline int ewd_offset(int i, int j, int k) { return (i * (ENY + 1) + j) * (ENZ + 1) + k; }
  inline double specerf(double z, double k, double alpha);
  inline double d_specerf(double z, double k, double alpha);
  inline double dd_specerf(double z, double k, double alpha);
  inline double ddd_specerf(double z, double k, double alpha);

  double Ewd_fac_intp[3];

  /*
   *   in D0phi we store the correction potential:
   *
   *      phi = 1/x + pi/alpha^2 - sum_q (erfc(alpha |x-q|)/|x-q|)  - 4pi/V sum_k exp(-k^2/(4alpha^2))/k^2 cos(k*x)
   *
   *   in D1phi we store the correction force (first derivative of correction potential)
   *
   *      dphi/dx_i
   *
   *   in D2Phi we store the correction tensor (second derivatives of correction potential)
   *
   *      d2phi/(dx_i dx_j)
   *
   *   in D3phi we store the third order correction tensor (third derivatives of correction potential)
   *
   *      d3phi/(dx_i dx_j dx_k)
   */

  int ewald_is_initialized = 0;

  void test_interpolation_accuracy(void);

  void ewald_interpolate(double dx, double dy, double dz, enum interpolate_options, ewald_data &fper);

 public:
#if defined(EVALPOTENTIAL) && defined(FMM) && defined(PERIODIC) && !defined(PMGRID)
  double ewald_gridlookup_origin_D0(void) { return Ewd->D0phi; }
#endif
};

extern ewald Ewald;

#endif
