/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file second_order_ics.cc
 *
 *  \brief produce actual ICs from special 2nd order LPT ICs created by Adrian Jenkin's code
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../gravtree/gravtree.h"
#include "../io/io.h"
#include "../main/simulation.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"

#ifdef SECOND_ORDER_LPT_ICS

double sim::F1_Omega(double a)
{
  double omega_a = All.Omega0 / (All.Omega0 + a * (1 - All.Omega0 - All.OmegaLambda) + a * a * a * All.OmegaLambda);

  return pow(omega_a, 5.0 / 9);
}

double sim::F2_Omega(double a)
{
  double omega_a = All.Omega0 / (All.Omega0 + a * (1 - All.Omega0 - All.OmegaLambda) + a * a * a * All.OmegaLambda);

  return 2.0 * pow(omega_a, 6.0 / 11);
}

void sim::second_order_ic_correction(void)
{
  if(executed)
    return;

  mpi_printf("SECOND_ORDER_LPT_ICS: Now producing ICs based on 2nd-order LPT\n");

  particle_data *P = Sp.P;

  /* first, back up masses and set special 2nd LPT masses (which have been read in in the field OldAcc) */

  double *mass_bak = (double *)Mem.mymalloc("mass_bak", Sp.NumPart * sizeof(double));

  for(int i = 0; i < Sp.NumPart; i++)
    {
      mass_bak[i] = P[i].getMass();
      P[i].setMass(P[i].OldAcc);
    }

    /* now do the gravity computation */

#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
  gravity_long_range_force();
#endif
  gravity(0);

  if(All.TypeOfOpeningCriterion == 1)
    gravity(0);

  /* correct the ICs accordingly */

  double a = All.TimeBegin;

  double hubble = Driftfac.hubble_function(a);

  double fac1 = 1.0 / (a * a * hubble * F1_Omega(a)); /* this factor converts the Zeldovich peculiar velocity
                                                 (expressed as  a^2*dX/dt  to comoving displacement */
  double fac2 = All.LptScalingfactor;

  double fac3 = fac2 * a * a * hubble * F2_Omega(a);

  mpi_printf("SECOND_ORDER_LPT_ICS: fac1=%g  fac2=%g  fac3=%g\n", fac1, fac2, fac3);

  for(int i = 0; i < Sp.NumPart; i++)
    {
      double posdiff[3];
      for(int j = 0; j < 3; j++)
        posdiff[j] = fac1 * P[i].Vel[j]; /* Zeldovich displacement */

      MyIntPosType delta[3];
      Sp.pos_to_signedintpos(posdiff, (MySignedIntPosType *)delta);

      for(int j = 0; j < 3; j++)
        P[i].IntPos[j] += delta[j];

      double acc[3];
#if defined(PMGRID) && defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
      for(int j = 0; j < 3; j++)
        acc[j] = P[i].GravPM[j] + P[i].GravAccel[j];
#else
      for(int j = 0; j < 3; j++)
        acc[j] = P[i].GravAccel[j];
#endif

      for(int j = 0; j < 3; j++)
        posdiff[j] = fac2 * acc[j]; /* second order lpt displacement */

      Sp.pos_to_signedintpos(posdiff, (MySignedIntPosType *)delta);

      for(int j = 0; j < 3; j++)
        P[i].IntPos[j] += delta[j];

      for(int j = 0; j < 3; j++)
        P[i].Vel[j] += fac3 * acc[j]; /* second order lpt velocity correction */
    }

  /* now restore the correct masses */

  for(int i = 0; i < Sp.NumPart; i++)
    {
      P[i].setMass(mass_bak[i]);
      P[i].OldAcc = 0;
    }

  Mem.myfree(mass_bak);

  mpi_printf("SECOND_ORDER_LPT_ICS: finished,\n");

  executed = 1;

  if(All.TypeOfOpeningCriterion == 1)
    All.RelOpeningCriterionInUse = 0;

  /* put in an extra domain decomposition because particle positions have been shifted */

  NgbTree.treefree();
  Domain.domain_free();
  Domain.domain_decomposition(STANDARD);
  NgbTree.treeallocate(Sp.NumGas, &Sp, &Domain);
  NgbTree.treebuild(Sp.NumGas, NULL);
}

#endif
