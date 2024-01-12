/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/
/*! \file allvars.h
 *  \brief declares a structure for global parameters and variables.
 *
 *  This file declares a structure for holding global parameters and variables and objects.
 *  Further variables should be added here. The actual instance of this object is provided
 *  by the file 'allvars.cc'.
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include "gadgetconfig.h"

#include <math.h>

#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/macros.h"
#include "../io/parameters.h"

/** Data which is the SAME for all tasks (mostly code parameters read
 * from the parameter file).  Holding this data in a structure is
 * convenient for writing/reading the restart file, and it allows the
 * introduction of new global variables in a simple way. The only
 * thing to do is to introduce them into this structure.
 */
struct global_data_all_processes : public parameters
{
#if defined(COOLING)
  char TreecoolFile[255];
#endif

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  double AvgType1Mass;
#endif

  double TopNodeFactor;

  int ICFormat; /**< selects different versions of IC file-format */

  int SnapFormat; /**< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /**< number of files in multi-file snapshot dumps */
  int MaxFilesWithConcurrentIO; /**< maximum number of files that may be written simultaneously when
                                      writing/reading restart-files, or when writing snapshot files */

  double TreeAllocFactor; /**< Each processor allocates a number of nodes which is TreeAllocFactor times
                               the maximum(!) number of particles.  Note: A typical local tree for N
                               particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor; /**< Each processor allocates a number of nodes which is TreeAllocFactor times
                                  the maximum(!) number of particles.  Note: A typical local tree for N
                                  particles needs usually about ~0.65*N nodes. */

  double NgbTreeAllocFactor; /**< Each processor allocates a number of nodes for the neighbor search which is NgbTreeAllocFactor times
                                   the maximum(!) number of gas particles.  Note: A typical local tree for N
                                   particles needs usually about ~0.65*N nodes. */

  double ForeignNodeAllocFactor;

  double ForeignPointAllocFactor;

  enum restart_options RestartFlag; /**< taken from command line used to start code. 0 is normal start-up from
                                         initial conditions, 1 is resuming a run from a set of restart files, while 2
                                         marks a restart from a snapshot file. */

  int MaxMemSize; /**< size of maximum memory consumption in MB */

  /* some SPH parameters */

  int DesNumNgb; /**< Desired number of SPH neighbours */

#ifdef LEAN
  MyDouble PartMass;
#endif

#ifdef TIMEDEP_ART_VISC
  double AlphaMin; /*!< Minimum of allowed viscosity parameter */
#endif

#ifdef SUBFIND
  int DesLinkNgb;
#endif

  long long GlobalNSynchronizedHydro;
  long long GlobalNSynchronizedGravity;

  double MaxNumNgbDeviation; /**< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst; /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;      /**< may be used to set the temperature in the IC's */
  double InitGasU;         /**< the same, but converted to thermal energy per unit mass */
  double MinEgySpec;       /**< the minimum allowed temperature expressed as energy per unit mass */

  /* some force counters  */

  long long TotNumOfForces;     /**< counts total number of force computations  */
  long long TotNumDirectForces; /**< counts total number of direct force computations  */
  long long TotNumDensity;      /**< counts total number of SPH density calculations  */
  long long TotNumHydro;        /**< counts total number of SPH hydro force calculations  */

  /* various cosmological factors that are only a function of the current scale factor, and in non-comoving runs are set to 1 */
  double cf_atime, cf_atime2, cf_ainv, cf_a2inv, cf_a3inv, cf_afac1, cf_afac2, cf_afac3, cf_hubble_a, cf_atime_hubble_a,
      cf_atime2_hubble_a, cf_redshift;
  /* Hubble rate at the current time, valid both for comoving and non-comoving intergation */
  double cf_H;

  double accel_normalize_fac; /* used in I/O to normalize accelerations if reduced precision storage is used */

  /* system of units  */
  double UnitTime_in_s;            /**< factor to convert internal time unit to seconds/h */
  double UnitMass_in_g;            /**< factor to convert internal mass unit to grams/h */
  double UnitVelocity_in_cm_per_s; /**< factor to convert internal velocity unit to cm/sec */
  double UnitLength_in_cm;         /**< factor to convert internal length unit to cm/h */
  double UnitPressure_in_cgs;      /**< factor to convert internal pressure unit to cgs units (little 'h' still
                                   around!) */
  double UnitDensity_in_cgs;       /**< factor to convert internal length unit to g/cm^3*h^2 */
  double UnitCoolingRate_in_cgs;   /**< factor to convert internal cooling rate to cgs units */
  double UnitEnergy_in_cgs;        /**< factor to convert internal energy to cgs units */
  double UnitTime_in_Megayears;    /**< factor to convert internal time to megayears/h */
  double UnitTime_in_years;        /**< factor to convert internal time to years/h */
  double GravityConstantInternal;  /**< If set to zero in the parameterfile, the internal value of the
                                   gravitational constant is set to the Newtonian value based on the system of
                                   units specified. Otherwise the value provided is taken as internal gravity
                                   constant G. */
  double G;                        /**< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;         /**< Hubble-constant in internal units */
  double Omega0;         /**< matter density in units of the critical density (at z=0) */
  double OmegaLambda;    /**< vaccum energy density relative to crictical density (at z=0) */
  double OmegaBaryon;    /**< baryon density in units of the critical density (at z=0) */
  double OmegaCurvature; /**< curvature relative to crictical density (at z=0) */
  double HubbleParam; /**< little `h', i.e. can be used to scale unit system to absorb uncertain value of Hubble constant.  Only needed
                       * to get absolute physical values for cooling physics
                       */

  double BoxSize; /**< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;  /**< flags that comoving integration is enabled */
  int TypeOfOpeningCriterion; /**< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
                                   criterion */
  int OutputListOn;           /**< flags that output times are listed in a specified file */

  int LowestActiveTimeBin;
  int HighestActiveTimeBin;
  int LowestOccupiedTimeBin;
  int HighestOccupiedTimeBin;
  int LowestOccupiedGravTimeBin;
  int HighestOccupiedGravTimeBin;
  int HighestSynchronizedTimeBin;
  int SmallestTimeBinWithDomainDecomposition;
  double ActivePartFracForNewDomainDecomp;
  /* parameters determining output frequency */

#if defined(TREEPM_NOTIMESPLIT) || defined(PLACEHIGHRESREGION)
  double ActivePartFracForPMinsteadOfEwald;
#endif

  int SnapshotFileCount;        /**< number of snapshot that is written next */
  double TimeBetSnapshot;       /**< simulation time interval between snapshot files */
  double TimeOfFirstSnapshot;   /**< simulation time of first snapshot files */
  double CpuTimeBetRestartFile; /**< cpu-time between regularly generated restart files */
  double TimeLastRestartFile;   /**< cpu-time when last restart-file was written */
  double TimeBetStatistics;     /**< simulation time interval between computations of energy statistics */
  double TimeLastStatistics;    /**< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;         /**< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time;      /**< current time of the simulation */
  double TimeBegin; /**< time of initial conditions of the simulation */
  double TimeStep;  /**< difference between current times of previous and current timestep */
  double TimeMax;   /**< marks the point of time until the simulation is to be evolved */
  double TimeOld;   /**< time of previous synchronization point, needed only for logging purposes */

  /* variables for organizing discrete timeline */

  double Timebase_interval;  /**< factor to convert from floating point time interval to integer timeline */
  integertime Ti_Current;    /**< current time on integer timeline */
  integertime Ti_nextoutput; /**< next output time on integer timeline */
  integertime Ti_lastoutput;

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  integertime PM_Ti_endstep, PM_Ti_begstep;
#endif

#if defined(EVALPOTENTIAL) && defined(PMGRID) && defined(PERIODIC)
  double TotalMass;
#endif

  inline double get_absolutetime_from_integertime(integertime ti)
  {
    if(ComovingIntegrationOn)
      return TimeBegin * exp(ti * Timebase_interval);
    else
      return TimeBegin + ti * Timebase_interval;
  }

  char DumpFlag_nextoutput;

  integertime Ti_begstep[TIMEBINS]; /**< marks start of current step of each timebin on integer timeline */

#ifdef FORCE_EQUAL_TIMESTEPS
  integertime GlobalTimeStep;
#endif

  /* variables that keep track of CPU consumption */

  double TimeLimitCPU;

  double CPUForLastPMExecution;

  /* tree code opening criterion */

  double ErrTolTheta;    /**< BH tree opening angle */
  double ErrTolThetaMax; /**< maximum BH tree opening angle when relative criterion is in use */
  double ErrTolForceAcc; /**< parameter for relative opening criterion in tree walk */

  char RelOpeningCriterionInUse; /**< flags that we now an old acceleration and the relative opening criterion can be used */

  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy; /**< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                                 timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep; /**< minimum allowed timestep. Normally, the simulation terminates if the
                               timestep determined by the timestep criteria falls below this limit. */
  double MaxSizeTimestep; /**< maximum allowed timestep */

  double CourantFac; /**< SPH-Courant factor */

  int CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];

  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */

  int SofteningClassOfPartType[NTYPES];

  double SofteningComoving[NSOFTCLASSES]; /**< comoving gravitational softening lengths for each softeniung type */
  double SofteningMaxPhys[NSOFTCLASSES];  /**< maximum physical gravitational softening lengths for each softening type */

  double SofteningTable[NSOFTCLASSES +
                        NSOFTCLASSES_HYDRO]; /**< current (comoving) gravitational softening lengths for each softening type */
  double ForceSoftening[NSOFTCLASSES + NSOFTCLASSES_HYDRO +
                        2]; /**< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */

#ifdef ADAPTIVE_HYDRO_SOFTENING
  double MinimumComovingHydroSoftening;
  double AdaptiveHydroSofteningSpacing;
  double GasSoftFactor;
#endif

  /** If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[NTYPES];

  /* some filenames */
  char InitCondFile[MAXLEN_PATH], OutputDir[MAXLEN_PATH], SnapshotFileBase[MAXLEN_PATH], OutputListFilename[MAXLEN_PATH];

  /** table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];
  char OutputListFlag[MAXLEN_OUTPUTLIST];
  int OutputListLength; /**< number of times stored in table of desired output times */

#ifdef SECOND_ORDER_LPT_ICS
  double LptScalingfactor;
#endif

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  double DtDisplacement;
#endif

#ifdef LIGHTCONE

#ifdef LIGHTCONE_PARTICLES
  char LightConeDefinitionFile[MAXLEN_PATH];
  int LightconeFileCount;

#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  char LightConeOriginsFile[MAXLEN_PATH];
#endif
#endif

#ifdef LIGHTCONE_MASSMAPS
  int LightConeMassMapsNside;
  double LightConeMassMapThickness;
  double LightConeMassMapMaxRedshift;
  int CurrentMassMapBoundary;
#endif

  int LightConeImageConeNr;
  double LightConeImageLengthX;
  double LightConeImageLengthY;
  double LightConeImageCornerX;
  double LightConeImageCornerY;
  int LightConeImagePixelsX;
  int LightConeImagePixelsY;
  char LightConeImagePicName[MAXLEN_PATH];
  int LightConeImageFirstConeDir;
  int LightConeImageLastConeDir;
#endif

#ifdef STARFORMATION /* star formation and feedback sector */
  double CritOverDensity;
  double CritPhysDensity;
  double OverDensThresh;
  double PhysDensThresh;
  double EgySpecSN;
  double EgySpecCold;
  double FactorEVP;
  double TempSupernova;
  double TempClouds;
  double MaxSfrTimescale;
  double FactorSN;
  MyIDType MaxID;
#endif

#ifdef REDUCE_FLUSH
  double FlushCpuTimeDiff;
  double FlushLast;
#endif

#ifdef NGENIC
  int NSample;
  int SphereMode;
  int PowerSpectrumType;
  int ReNormalizeInputSpectrum;
  double PrimordialIndex;
  double ShapeGamma;
  double Sigma8;
  char PowerSpectrumFile[MAXLEN_PATH];
  double InputSpectrum_UnitLength_in_cm;
  int NgenicSeed;
#endif

#ifdef CREATE_GRID
  int GridSize;
#endif

#ifdef EXTERNALGRAVITY_STATICHQ
  double A_StaticHQHalo;
  double Mass_StaticHQHalo;
#endif

  void set_cosmo_factors_for_current_time(void);
  void register_parameters(void);
  void read_outputlist(char *fname);
  void some_parameter_checks(void);

  inline char *get_data_ptr(void) { return (char *)this + sizeof(parameters); }

  inline size_t get_data_size(void) { return sizeof(global_data_all_processes) - sizeof(parameters); }
};

extern global_data_all_processes All;

#endif
