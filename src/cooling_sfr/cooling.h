/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file cooling.h
 *
 *  \brief defines a class for dealing with cooling and star formation
 */

#ifndef COOLING_H
#define COOLING_H

#include "gadgetconfig.h"

#ifdef COOLING

#include "../data/simparticles.h"
#include "../mpi_utils/setcomm.h"

class coolsfr : public setcomm
{
 public:
  coolsfr(MPI_Comm comm) : setcomm(comm) {}

  double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer);

  void InitCool(void);
  void IonizeParams(void);

  void cooling_only(simparticles *Sp);

#ifdef STARFORMATION
  void sfr_create_star_particles(simparticles *Sp);

  void set_units_sfr(void);

  void cooling_and_starformation(simparticles *Sp);

  void init_clouds(void);
#endif

 private:
#define NCOOLTAB 2000

  /* data for gas state */
  struct gas_state
  {
    double ne, necgs, nHcgs;
    double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
    double gJH0ne, gJHe0ne, gJHepne;
    double nH0, nHp, nHep, nHe0, nHepp;
    double XH, yhelium;
    double mhboltz;
    double ethmin; /* minimum internal energy for neutral gas */
  };

  /* tabulated rates */
  struct rate_table
  {
    double BetaH0, BetaHep, Betaff;
    double AlphaHp, AlphaHep, Alphad, AlphaHepp;
    double GammaeH0, GammaeHe0, GammaeHep;
  };

  /* photo-ionization/heating rate table */
  struct photo_table
  {
    float variable;       /* logz for UVB */
    float gH0, gHe, gHep; /* photo-ionization rates */
    float eH0, eHe, eHep; /* photo-heating rates */
  };

  /* current interpolated photo-ionization/heating rates */
  struct photo_current
  {
    char J_UV;
    double gJH0, gJHep, gJHe0, epsH0, epsHep, epsHe0;
  };

  /* cooling data */
  struct do_cool_data
  {
    double u_old_input, rho_input, dt_input, ne_guess_input;
  };

  gas_state GasState;      /**< gas state */
  do_cool_data DoCoolData; /**< cooling data */

  rate_table *RateT;      /**< tabulated rates */
  photo_table *PhotoTUVB; /**< photo-ionization/heating rate table for UV background*/
  photo_current pc;       /**< current interpolated photo rates */

  double Tmin = 1.0; /**< min temperature in log10 */
  double Tmax = 9.0; /**< max temperature in log10 */
  double deltaT;     /**< log10 of temperature spacing in the interpolation tables */
  int NheattabUVB;   /**< length of UVB photo table */

#ifdef COOLING
  double DoCooling(double u_old, double rho, double dt, double *ne_guess, gas_state *gs, do_cool_data *DoCool);
  double GetCoolingTime(double u_old, double rho, double *ne_guess, gas_state *gs, do_cool_data *DoCool);
  void cool_sph_particle(simparticles *Sp, int i, gas_state *gs, do_cool_data *DoCool);

  void SetZeroIonization(void);
#endif

  void integrate_sfr(void);

  double CoolingRate(double logT, double rho, double *nelec, gas_state *gs, const do_cool_data *DoCool);
  double CoolingRateFromU(double u, double rho, double *ne_guess, gas_state *gs, const do_cool_data *DoCool);
  void find_abundances_and_rates(double logT, double rho, double *ne_guess, gas_state *gs, const do_cool_data *DoCool);
  void IonizeParamsUVB(void);
  void ReadIonizeParams(char *fname);

  double convert_u_to_temp(double u, double rho, double *ne_guess, gas_state *gs, const do_cool_data *DoCool);

  void MakeRateTable(void);

#ifdef STARFORMATION
  const int WriteMiscFiles = 1;

  void make_star(simparticles *Sp, int i, double prob, MyDouble mass_of_star, double *sum_mass_stars);
  void spawn_star_from_sph_particle(simparticles *Sp, int igas, double birthtime, int istar, MyDouble mass_of_star);
  void convert_sph_particle_into_star(simparticles *Sp, int i, double birthtime);

  int stars_spawned;           /**< local number of star particles spawned in the time step */
  int tot_stars_spawned;       /**< global number of star paricles spawned in the time step */
  int stars_converted;         /**< local number of gas cells converted into stars in the time step */
  int tot_stars_converted;     /**< global number of gas cells converted into stars in the time step */
  int altogether_spawned;      /**< local number of star+wind particles spawned in the time step */
  int tot_altogether_spawned;  /**< global number of star+wind particles spawned in the time step */
  double cum_mass_stars = 0.0; /**< cumulative mass of stars created in the time step (global value) */
#endif
};

#endif
#endif
