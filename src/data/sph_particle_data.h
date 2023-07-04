/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file sph_particle_data.h
 *
 *  \brief defines the structure holding the extra hydrodynamic data for a single SPH particle
 */

#ifndef SPHPARTDATA_H
#define SPHPARTDATA_H

#include "gadgetconfig.h"

#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/macros.h"
#include "../data/mymalloc.h"
#include "../mpi_utils/setcomm.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/* in this structure, all SPH variables are put that are needed for passive
 * particles in the hydro force calculation. Only this part will be sent
 * to other nodes if needed
 */
struct sph_particle_data_hydrocore
{
  MyFloat Hsml;               /*!< current smoothing length */
  MyFloat DhsmlDensityFactor; /*!< correction factor needed in entropy formulation of SPH */
  MyFloat VelPred[3];         /*!< predicted SPH particle velocity at the current time, needed if particle is inactive */

  MyFloat DivVel;  /*!< local velocity divergence */
  MyFloat CurlVel; /*!< local velocity curl */
  MyFloat Csnd;

  MyFloat Density;  /*!< current baryonic mass density of particle */
  MyFloat Pressure; /*!< current pressure */

#ifdef TIMEDEP_ART_VISC
  MyFloat Alpha; /*!< time-dependend viscosity parameter */
#endif
#ifdef PRESSURE_ENTROPY_SPH
  MyFloat EntropyToInvGammaPred;     /*!< current entropy function A to the power 1 / gamma */
  MyFloat DhsmlDerivedDensityFactor; /*!< additional correction factor needed for pressure formulation of SPH */
  MyFloat PressureSphDensity;        /* current density derived from the pressure estimate */
#endif
};

/** Holds data that is stored for each sph particle in addition to
    the collisionless variables.
 */
struct sph_particle_data : public sph_particle_data_hydrocore
{
  MyFloat Entropy;     /*!< value of the entropic function */
  MyFloat EntropyPred; /*!< predicted entropy at current time, needed if the particle is inactive */

  MyFloat HydroAccel[3]; /*!< acceleration due to hydrodynamical forces */
#ifdef HIERARCHICAL_GRAVITY
  MyFloat FullGravAccel[3]; /*!< most recent full calculation of gravitational acceleration, used to advanced VelPred */
#endif
  MyFloat DtEntropy; /*!< rate of change of entropy */
  MyFloat DtDensity; /*!< rate of change of density, needed to predict densities for passive particles */
  MyFloat DtHsml;    /*!< rate of change of smoothing length, needed to predict hsml for passive particles */

  MyFloat NumNgb; /*!< effective number of neighbours used in density estimation loop (note: this could be changed to a temporary
                     variable in density) */

  MyFloat Rot[3]; /*!< local velocity curl */

  MyFloat MaxSignalVel; /*!< maximum signal velocity */
  MyFloat CurrentMaxTiStep;

#ifdef PRESSURE_ENTROPY_SPH
  MyFloat
      DtPressureSphDensity; /*!< rate of change of the pressure derived density, needed to predict densities for passive particles */
#endif

#ifdef TIMEDEP_ART_VISC
  MyFloat DivVelOld; /* local velocity gradient from the previous time step */
  MyFloat decayVel;  /* decay velocity for the viscosity parameter */
#endif

#ifdef IMPROVED_VELOCITY_GRADIENTS
  struct
  {
    MyFloat dx_dx;
    MyFloat dx_dy;
    MyFloat dx_dz;
    MyFloat dy_dy;
    MyFloat dy_dz;
    MyFloat dz_dz;
  } dpos; /* contains the matrix elements needed for the improved gradient estimate */

  MyFloat dvel[NUMDIMS][NUMDIMS]; /* contains the velocity gradients */
#endif

#ifdef STARFORMATION
  MyFloat Metallicity;
  MyFloat MassMetallicity;
#endif

#ifdef COOLING
  MyFloat Ne; /*!< free electron fraction, expressed as local electron number density normalized to the hydrogen number density. Gives
                 indirectly mean molecular weight. */
#endif

#ifdef OUTPUT_COOLHEAT
  MyFloat CoolHeat;
#endif

#ifdef STARFORMATION
  MyFloat Sfr;
#endif

  inline MyFloat get_sound_speed(void)
  {
    MyFloat csnd;

    if(Density > 0)
      csnd = sqrt(static_cast<MyReal>(GAMMA) * Pressure / Density);
    else
      csnd = 0;

    return csnd;
  }

  /* compute the pressure of particle i */
  inline MyFloat get_pressure(void)
  {
#ifndef PRESSURE_ENTROPY_SPH
    return EntropyPred * pow(Density, (MyFloat)GAMMA);
#else
    return pow(EntropyToInvGammaPred * PressureSphDensity, GAMMA);
#endif
  }

  inline void set_thermodynamic_variables(void)
  {
    Pressure = get_pressure();

    if(Pressure < 0)
      Terminate("Pressure=%g  rho=%g  entr=%g entrpred=%g\n", Pressure, Density, Entropy, EntropyPred);

    Csnd = get_sound_speed();
  }

  inline MyFloat get_Hsml() { return Hsml; }

#ifdef IMPROVED_VELOCITY_GRADIENTS
  void set_velocity_gradients(void);
#endif

#ifdef TIMEDEP_ART_VISC
  void set_viscosity_coefficient(double dt);
#endif
};

#endif
