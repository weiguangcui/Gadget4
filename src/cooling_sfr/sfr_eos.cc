/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file sfr_eos.cc
 *
 *  \brief Star formation rate routines for the effective multi-phase model
 */

#include "gadgetconfig.h"

#ifdef STARFORMATION

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../logs/logs.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/** \brief Main driver for star formation and gas cooling.
 *
 *  This function loops over all the active gas cells. If a given cell
 *  meets the criteria for star formation to be active the multi-phase
 *  model is activated, the properties of the cell are updated according to
 *  the latter and the star formation rate computed. In the other case, the
 *  standard isochoric cooling is applied to the gas cell by calling the function
 *  cool_sph_particle() and the star formation rate is set to 0.
 */
void coolsfr::cooling_and_starformation(simparticles *Sp)
{
  TIMER_START(CPU_COOLING_SFR);

  /* clear the SFR stored in the active timebins */
  for(int bin = 0; bin < TIMEBINS; bin++)
    if(Sp->TimeBinSynchronized[bin])
      Sp->TimeBinSfr[bin] = 0;

  All.set_cosmo_factors_for_current_time();

  gas_state gs    = GasState;
  do_cool_data dc = DoCoolData;

  for(int i = 0; i < Sp->TimeBinsHydro.NActiveParticles; i++)
    {
      int target = Sp->TimeBinsHydro.ActiveParticleList[i];
      if(Sp->P[target].getType() == 0)
        {
          if(Sp->P[target].getMass() == 0 && Sp->P[target].ID.get() == 0)
            continue; /* skip cells that have been swallowed or eliminated */

          double dens = Sp->SphP[target].Density;

          double dt =
              (Sp->P[target].getTimeBinHydro() ? (((integertime)1) << Sp->P[target].getTimeBinHydro()) : 0) * All.Timebase_interval;
          /*  the actual time-step */

          double dtime = All.cf_atime * dt / All.cf_atime_hubble_a;

          /* check whether conditions for star formation are fulfilled.
           *
           * f=1  normal cooling
           * f=0  star formation
           */
          int flag = 1; /* default is normal cooling */

          if(dens * All.cf_a3inv >= All.PhysDensThresh)
            flag = 0;

          if(All.ComovingIntegrationOn)
            if(dens < All.OverDensThresh)
              flag = 1;

          if(flag == 1) /* normal implicit isochoric cooling */
            {
              Sp->SphP[target].Sfr = 0;
              cool_sph_particle(Sp, target, &gs, &dc);
            }

          if(flag == 0) /* active star formation */
            {
              double tsfr = sqrt(All.PhysDensThresh / (dens * All.cf_a3inv)) * All.MaxSfrTimescale;

              double factorEVP = pow(dens * All.cf_a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

              double egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

              double ne = Sp->SphP[target].Ne;

              double tcool        = GetCoolingTime(egyhot, dens * All.cf_a3inv, &ne, &gs, &dc);
              Sp->SphP[target].Ne = ne;

              double y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

              double x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

              double egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

              double cloudmass = x * Sp->P[target].getMass();
              double utherm    = Sp->get_utherm_from_entropy(target);

              if(tsfr < dtime)
                tsfr = dtime;

              if(dt > 0)
                {
                  if(Sp->P[target].getTimeBinHydro()) /* upon start-up, we need to protect against dt==0 */
                    {
                      double trelax     = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
                      double egycurrent = utherm;

                      double unew;
                      if(egycurrent > egyeff)
                        {
                          double dtcool = dtime;
                          ne            = Sp->SphP[target].Ne;

                          unew = DoCooling(egycurrent, dens * All.cf_a3inv, dtcool, &ne, &gs, &dc);

                          if(unew < egyeff)
                            {
                              unew = egyeff;
                            }
                        }
                      else
                        unew = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));

                      double du = unew - utherm;
                      if(unew < All.MinEgySpec)
                        du = All.MinEgySpec - utherm;

                      utherm += du;
                      Sp->set_entropy_from_utherm(utherm, target);
                      Sp->SphP[target].DtEntropy = 0.0;

#ifdef OUTPUT_COOLHEAT
                      if(dtime > 0)
                        Sp->SphP[target].CoolHeat = du * Sp->P[target].getMass() / dtime;
#endif
                      Sp->SphP[target].set_thermodynamic_variables();
                    }
                }

              if(utherm > 1.01 * egyeff)
                Sp->SphP[target].Sfr = 0;
              else
                {
                  /* note that we convert the star formation rate to solar masses per year */
                  Sp->SphP[target].Sfr =
                      (1 - All.FactorSN) * cloudmass / tsfr * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
                }

              Sp->TimeBinSfr[Sp->P[target].getTimeBinHydro()] += Sp->SphP[target].Sfr;
            }
        }
    } /* end of main loop over active particles */

  TIMER_STOP(CPU_COOLING_SFR);
}

/** \brief Initialize the parameters of effective multi-phase model.
 *
 *   In particular this function computes the value of PhysDensThresh, that is
 *   the physical density threshold above which star formation is active, if its
 *   value was set to 0 in the parameter file.
 */
void coolsfr::init_clouds(void)
{
  gas_state gs    = GasState;
  do_cool_data dc = DoCoolData;

  if(All.PhysDensThresh == 0)
    {
      double A0 = All.FactorEVP;

      double egyhot = All.EgySpecSN / A0;

      double meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* note: assuming FULL ionization */

      double u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

      double dens;
      if(All.ComovingIntegrationOn)
        dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
        dens = 1.0e6 * (1.0e-29 / All.UnitDensity_in_cgs);

      if(All.ComovingIntegrationOn)
        {
          All.Time = 1.0; /* to be guaranteed to get z=0 rate */
          All.set_cosmo_factors_for_current_time();
          IonizeParams();
        }

      double ne = 1.0;
      SetZeroIonization();

      double tcool = GetCoolingTime(egyhot, dens, &ne, &gs, &dc);

      double coolrate = egyhot / tcool / dens;

      double x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh =
          x / pow(1 - x, 2) * (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold) / (All.MaxSfrTimescale * coolrate);

      mpi_printf(
          "\nA0= %g  \nComputed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\nEXPECTED FRACTION OF COLD GAS AT THRESHOLD = "
          "%g\n\ntcool=%g dens=%g egyhot=%g\n",
          A0, All.PhysDensThresh, All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs), x, tcool, dens,
          egyhot);

      dens = All.PhysDensThresh * 10;

      double neff;
      do
        {
          double tsfr      = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
          double factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
          egyhot           = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

          ne = 0.5;

          tcool = GetCoolingTime(egyhot, dens, &ne, &gs, &dc);

          double y      = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
          x             = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
          double egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

          double peff = GAMMA_MINUS1 * dens * egyeff;

          double fac = 1 / (log(dens * 1.025) - log(dens));
          dens *= 1.025;

          neff = -log(peff) * fac;

          tsfr      = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
          factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
          egyhot    = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

          ne = 0.5;

          tcool = GetCoolingTime(egyhot, dens, &ne, &gs, &dc);

          y      = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
          x      = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
          egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

          peff = GAMMA_MINUS1 * dens * egyeff;

          neff += log(peff) * fac;
        }
      while(neff > 4.0 / 3);

      double thresholdStarburst = dens;

      mpi_printf("Run-away sets in for dens=%g\nDynamic range for quiescent star formation= %g\n", thresholdStarburst,
                 thresholdStarburst / All.PhysDensThresh);

      integrate_sfr();

      if(ThisTask == 0)
        {
          double sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

          printf("Isotherm sheet central density: %g   z0=%g\n", M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
                 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
          myflush(stdout);
        }

      if(All.ComovingIntegrationOn)
        {
          All.Time = All.TimeBegin;
          All.set_cosmo_factors_for_current_time();
          IonizeParams();
        }
    }
}

/** \brief Compute the effective equation of state for the gas and
 *         the integrated SFR per unit area.
 *
 *  This function computes the effective equation of state for the gas and
 *  the integrated SFR per unit area. It saves the results into two files:
 *  eos.txt for the equation of state and sfrrate.txt for the integrated SFR.
 *  In the latter case, the SFR is determined by integrating along the vertical
 *  direction the gas density of an infinite self-gravitating isothermal sheet.
 *  The integrated gas density is saved as well, so effectively sfrrate.txt
 *  contains the Kennicutt-Schmidt law of the star formation model.
 */
void coolsfr::integrate_sfr(void)
{
  double meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* note: assuming FULL ionization */
  double u4         = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
  gas_state gs    = GasState;
  do_cool_data dc = DoCoolData;

  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0; /* to be guaranteed to get z=0 rate */
      All.set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  FILE *fd = (WriteMiscFiles && (ThisTask == 0)) ? fopen("eos.txt", "w") : NULL;

  for(double rho = All.PhysDensThresh; rho <= 1000 * All.PhysDensThresh; rho *= 1.1)
    {
      double tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

      double factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

      double egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

      double ne = 1.0;

      double tcool = GetCoolingTime(egyhot, rho, &ne, &gs, &dc);

      double y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
      double x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      double egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      double P0 = GAMMA_MINUS1 * rho * egyeff;

      if(WriteMiscFiles && (ThisTask == 0))
        {
          fprintf(fd, "%g %g\n", rho, P0);
        }
    }

  if(WriteMiscFiles && (ThisTask == 0))
    {
      fclose(fd);
      fd = fopen("sfrrate.txt", "w");
    }

  for(double rho0 = All.PhysDensThresh; rho0 <= 10000 * All.PhysDensThresh; rho0 *= 1.02)
    {
      double rho = rho0;
      double q   = 0;
      double dz  = 0.001;

      double sigma = 0, sigmasfr = 0, sigma_u4 = 0, x = 0;

      while(rho > 0.0001 * rho0)
        {
          double tsfr, P0, gam;
          if(rho > All.PhysDensThresh)
            {
              tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

              double factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

              double egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

              double ne = 1.0;

              double tcool = GetCoolingTime(egyhot, rho, &ne, &gs, &dc);

              double y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
              x        = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

              double egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

              P0        = GAMMA_MINUS1 * rho * egyeff;
              double P1 = P0;

              double rho2       = 1.1 * rho;
              double tsfr2      = sqrt(All.PhysDensThresh / rho2) * All.MaxSfrTimescale;
              double factorEVP2 = pow(rho2 / All.PhysDensThresh, -0.8) * All.FactorEVP;
              double egyhot2    = All.EgySpecSN / (1 + factorEVP2) + All.EgySpecCold;

              double tcool2  = GetCoolingTime(egyhot2, rho2, &ne, &gs, &dc);
              double y2      = tsfr2 / tcool2 * egyhot2 / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
              double x2      = 1 + 1 / (2 * y2) - sqrt(1 / y2 + 1 / (4 * y2 * y2));
              double egyeff2 = egyhot2 * (1 - x2) + All.EgySpecCold * x2;
              double P2      = GAMMA_MINUS1 * rho2 * egyeff2;

              gam = log(P2 / P1) / log(rho2 / rho);
            }
          else
            {
              tsfr = 0;

              P0  = GAMMA_MINUS1 * rho * u4;
              gam = 1.0;

              sigma_u4 += rho * dz;
            }

          double drho = q;
          double dq   = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P0) * rho * rho * rho;

          sigma += rho * dz;
          if(tsfr > 0)
            {
              sigmasfr += (1 - All.FactorSN) * rho * x / tsfr * dz;
            }

          rho += drho * dz;
          q += dq * dz;
        }

      sigma *= 2; /* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;

      sigma *= All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);
      sigmasfr *= All.HubbleParam * All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * (SEC_PER_YEAR / All.UnitTime_in_s);
      sigma_u4 *= All.HubbleParam * (All.UnitMass_in_g / SOLAR_MASS) * PARSEC * PARSEC / (All.UnitLength_in_cm * All.UnitLength_in_cm);

      if(WriteMiscFiles && (ThisTask == 0))
        {
          fprintf(fd, "%g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4);
        }
    }

  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      All.set_cosmo_factors_for_current_time();
      IonizeParams();
    }

  if(WriteMiscFiles && (ThisTask == 0))
    fclose(fd);
}

/** \brief Set the appropriate units for the parameters of the multi-phase model.
 */
void coolsfr::set_units_sfr(void)
{
  All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  double meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC); /* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
}
#endif
