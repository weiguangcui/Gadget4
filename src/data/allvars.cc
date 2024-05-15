/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file allvars.cc
 *
 *  \brief instance and code for an object dealing with global parameters and variables
 */

// clang-format off
#include "gadgetconfig.h"
// clang-format on

#include "../data/allvars.h"

#include "../data/constants.h"
#include "../data/dtypes.h"
#include "../data/macros.h"
#include "../time_integration/driftfac.h"

void global_data_all_processes::set_cosmo_factors_for_current_time(void)
{
  if(ComovingIntegrationOn)
    {
      cf_atime    = Time;
      cf_atime2   = Time * Time;
      cf_ainv     = 1 / Time;
      cf_a2inv    = 1 / (Time * Time);
      cf_a3inv    = 1 / (Time * Time * Time);
      cf_afac1    = pow(Time, 3 * GAMMA_MINUS1);
      cf_afac2    = 1 / pow(Time, 3 * GAMMA - 2);
      cf_afac3    = pow(Time, 3 * (1 - GAMMA) / 2.0);
      cf_hubble_a = cf_H = Driftfac.hubble_function(Time);
      cf_atime_hubble_a  = Time * cf_hubble_a;
      cf_atime2_hubble_a = Time * Time * cf_hubble_a;
      cf_redshift        = 1 / Time - 1;
    }
  else
    {
      cf_atime           = 1;
      cf_atime2          = 1;
      cf_ainv            = 1;
      cf_a2inv           = 1;
      cf_a3inv           = 1;
      cf_afac1           = 1;
      cf_afac2           = 1;
      cf_afac3           = 1;
      cf_hubble_a        = 1;
      cf_H               = 0;
      cf_atime_hubble_a  = 1;
      cf_atime2_hubble_a = 1;
      cf_redshift        = 0;
    }
}

void global_data_all_processes::register_parameters(void)
{
  add_param("InitCondFile", InitCondFile, PARAM_STRING, PARAM_FIXED);

  add_param("OutputDir", OutputDir, PARAM_STRING, PARAM_CHANGEABLE);
  add_param("SnapshotFileBase", SnapshotFileBase, PARAM_STRING, PARAM_CHANGEABLE);
  add_param("OutputListFilename", OutputListFilename, PARAM_STRING, PARAM_CHANGEABLE);
  add_param("OutputListOn", &OutputListOn, PARAM_INT, PARAM_CHANGEABLE);

  add_param("Omega0", &Omega0, PARAM_DOUBLE, PARAM_FIXED);
  add_param("OmegaBaryon", &OmegaBaryon, PARAM_DOUBLE, PARAM_FIXED);
  add_param("OmegaLambda", &OmegaLambda, PARAM_DOUBLE, PARAM_FIXED);
  add_param("Hubble", &Hubble, PARAM_DOUBLE, PARAM_FIXED);
  add_param("HubbleParam", &HubbleParam, PARAM_DOUBLE, PARAM_FIXED);
  add_param("BoxSize", &BoxSize, PARAM_DOUBLE, PARAM_FIXED);

  add_param("MaxMemSize", &MaxMemSize, PARAM_INT, PARAM_CHANGEABLE);
  add_param("TimeOfFirstSnapshot", &TimeOfFirstSnapshot, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("CpuTimeBetRestartFile", &CpuTimeBetRestartFile, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("TimeBetStatistics", &TimeBetStatistics, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("TimeBegin", &TimeBegin, PARAM_DOUBLE, PARAM_FIXED);
  add_param("TimeMax", &TimeMax, PARAM_DOUBLE, PARAM_FIXED); /* can be changed nevertheless through special function */

  add_param("TimeBetSnapshot", &TimeBetSnapshot, PARAM_DOUBLE, PARAM_CHANGEABLE);

  add_param("UnitVelocity_in_cm_per_s", &UnitVelocity_in_cm_per_s, PARAM_DOUBLE, PARAM_FIXED);
  add_param("UnitLength_in_cm", &UnitLength_in_cm, PARAM_DOUBLE, PARAM_FIXED);
  add_param("UnitMass_in_g", &UnitMass_in_g, PARAM_DOUBLE, PARAM_FIXED);
  add_param("GravityConstantInternal", &GravityConstantInternal, PARAM_DOUBLE, PARAM_FIXED);

  add_param("ErrTolIntAccuracy", &ErrTolIntAccuracy, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("ErrTolTheta", &ErrTolTheta, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("ErrTolThetaMax", &ErrTolThetaMax, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("ErrTolForceAcc", &ErrTolForceAcc, PARAM_DOUBLE, PARAM_CHANGEABLE);

  add_param("MaxSizeTimestep", &MaxSizeTimestep, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("MinSizeTimestep", &MinSizeTimestep, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("ArtBulkViscConst", &ArtBulkViscConst, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("CourantFac", &CourantFac, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("DesNumNgb", &DesNumNgb, PARAM_INT, PARAM_CHANGEABLE);
  add_param("TopNodeFactor", &TopNodeFactor, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("ActivePartFracForNewDomainDecomp", &ActivePartFracForNewDomainDecomp, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("MaxNumNgbDeviation", &MaxNumNgbDeviation, PARAM_DOUBLE, PARAM_CHANGEABLE);

  add_param("ComovingIntegrationOn", &ComovingIntegrationOn, PARAM_INT, PARAM_FIXED);

  add_param("ICFormat", &ICFormat, PARAM_INT, PARAM_CHANGEABLE);
  add_param("SnapFormat", &SnapFormat, PARAM_INT, PARAM_CHANGEABLE);

  add_param("NumFilesPerSnapshot", &NumFilesPerSnapshot, PARAM_INT, PARAM_CHANGEABLE);
  add_param("MaxFilesWithConcurrentIO", &MaxFilesWithConcurrentIO, PARAM_INT, PARAM_CHANGEABLE);

  add_param("TypeOfOpeningCriterion", &TypeOfOpeningCriterion, PARAM_INT, PARAM_FIXED);

  add_param("TimeLimitCPU", &TimeLimitCPU, PARAM_DOUBLE, PARAM_CHANGEABLE);

  add_param("InitGasTemp", &InitGasTemp, PARAM_DOUBLE, PARAM_FIXED);
  add_param("MinEgySpec", &MinEgySpec, PARAM_DOUBLE, PARAM_CHANGEABLE);

  for(int i = 0; i < NSOFTCLASSES; i++)
    {
      char buf_l[MAXLEN_PARAM_TAG];
      snprintf(buf_l, MAXLEN_PARAM_TAG, "SofteningComovingClass%d", i);
      add_param(buf_l, &SofteningComoving[i], PARAM_DOUBLE, PARAM_FIXED);
    }

  for(int i = 0; i < NSOFTCLASSES; i++)
    {
      char buf_l[MAXLEN_PARAM_TAG];
      snprintf(buf_l, MAXLEN_PARAM_TAG, "SofteningMaxPhysClass%d", i);
      add_param(buf_l, &SofteningMaxPhys[i], PARAM_DOUBLE, PARAM_FIXED);
    }

  for(int i = 0; i < NTYPES; i++)
    {
      char buf_l[MAXLEN_PARAM_TAG];
      snprintf(buf_l, MAXLEN_PARAM_TAG, "SofteningClassOfPartType%d", i);
      add_param(buf_l, &SofteningClassOfPartType[i], PARAM_INT, PARAM_FIXED);
    }

#if defined(TREEPM_NOTIMESPLIT) || defined(PLACEHIGHRESREGION)
  add_param("ActivePartFracForPMinsteadOfEwald", &ActivePartFracForPMinsteadOfEwald, PARAM_DOUBLE, PARAM_CHANGEABLE);
#endif

#ifdef ADAPTIVE_HYDRO_SOFTENING
  add_param("MinimumComovingHydroSoftening", &MinimumComovingHydroSoftening, PARAM_DOUBLE, PARAM_FIXED);
  add_param("AdaptiveHydroSofteningSpacing", &AdaptiveHydroSofteningSpacing, PARAM_DOUBLE, PARAM_FIXED);
  add_param("GasSoftFactor", &GasSoftFactor, PARAM_DOUBLE, PARAM_FIXED);
#endif

#ifdef TIMEDEP_ART_VISC
  add_param("ViscosityAlphaMin", &AlphaMin, PARAM_DOUBLE, PARAM_CHANGEABLE);
#endif

#ifdef SUBFIND
  add_param("DesLinkNgb", &DesLinkNgb, PARAM_INT, PARAM_CHANGEABLE);
#endif

#ifdef LIGHTCONE_PARTICLES
  add_param("LightConeDefinitionFile", LightConeDefinitionFile, PARAM_STRING, PARAM_CHANGEABLE);
#ifdef LIGHTCONE_MULTIPLE_ORIGINS
  add_param("LightConeOriginsFile", LightConeOriginsFile, PARAM_STRING, PARAM_CHANGEABLE);
#endif
#endif

#ifdef LIGHTCONE_MASSMAPS
  add_param("LightConeMassMapsNside", &LightConeMassMapsNside, PARAM_INT, PARAM_FIXED);
  add_param("LightConeMassMapThickness", &LightConeMassMapThickness, PARAM_DOUBLE, PARAM_CHANGEABLE);
  add_param("LightConeMassMapMaxRedshift", &LightConeMassMapMaxRedshift, PARAM_DOUBLE, PARAM_CHANGEABLE);
#endif

#ifdef REDUCE_FLUSH
  add_param("FlushCpuTimeDiff", &FlushCpuTimeDiff, PARAM_DOUBLE, PARAM_CHANGEABLE);
#endif

#ifdef COOLING
  add_param("TreecoolFile", TreecoolFile, PARAM_STRING, PARAM_CHANGEABLE);
#endif

#ifdef STARFORMATION
  add_param("CritOverDensity", &CritOverDensity, PARAM_DOUBLE, PARAM_FIXED);
  add_param("CritPhysDensity", &CritPhysDensity, PARAM_DOUBLE, PARAM_FIXED);
  add_param("FactorSN", &FactorSN, PARAM_DOUBLE, PARAM_FIXED);
  add_param("FactorEVP", &FactorEVP, PARAM_DOUBLE, PARAM_FIXED);
  add_param("TempSupernova", &TempSupernova, PARAM_DOUBLE, PARAM_FIXED);
  add_param("TempClouds", &TempClouds, PARAM_DOUBLE, PARAM_FIXED);
  add_param("MaxSfrTimescale", &MaxSfrTimescale, PARAM_DOUBLE, PARAM_FIXED);
#endif

#ifdef NGENIC
  add_param("NSample", &NSample, PARAM_INT, PARAM_FIXED);
  add_param("SphereMode", &SphereMode, PARAM_INT, PARAM_FIXED);
  add_param("PowerSpectrumType", &PowerSpectrumType, PARAM_INT, PARAM_FIXED);
  add_param("ReNormalizeInputSpectrum", &ReNormalizeInputSpectrum, PARAM_INT, PARAM_FIXED);
  add_param("PrimordialIndex", &PrimordialIndex, PARAM_DOUBLE, PARAM_FIXED);
  add_param("ShapeGamma", &ShapeGamma, PARAM_DOUBLE, PARAM_FIXED);
  add_param("Sigma8", &Sigma8, PARAM_DOUBLE, PARAM_FIXED);
  add_param("PowerSpectrumFile", PowerSpectrumFile, PARAM_STRING, PARAM_FIXED);
  add_param("InputSpectrum_UnitLength_in_cm", &InputSpectrum_UnitLength_in_cm, PARAM_DOUBLE, PARAM_FIXED);
  add_param("Seed", &NgenicSeed, PARAM_INT, PARAM_FIXED);
#endif

#ifdef CREATE_GRID
  add_param("GridSize", &GridSize, PARAM_INT, PARAM_FIXED);
#endif

#ifdef EXTERNALGRAVITY_STATICHQ
  add_param("A_StaticHQHalo", &A_StaticHQHalo, PARAM_DOUBLE, PARAM_FIXED);
  add_param("Mass_StaticHQHalo", &Mass_StaticHQHalo, PARAM_DOUBLE, PARAM_FIXED);
#endif
}

/*! \brief This function reads a table with a list of desired output times.
 *
 *  The table does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 *
 *  \param fname The file name of the outputlist
 */
void global_data_all_processes::read_outputlist(char *fname)
{
  if(ThisTask == 0)
    {
      FILE *fd;

      if(!(fd = fopen(fname, "r")))
        {
          Terminate("can't read output list in file '%s'\n", fname);
        }

      OutputListLength = 0;

      while(1)
        {
          char buf[512];
          if(fgets(buf, 500, fd) != buf)
            break;

          int flag;
          int count = sscanf(buf, " %lg %d ", &OutputListTimes[OutputListLength], &flag);

          if(count == 1)
            flag = 1;

          if(count == 1 || count == 2)
            {
              if(OutputListLength >= MAXLEN_OUTPUTLIST)
                Terminate("\ntoo many entries in output-list. You should increase MAXLEN_OUTPUTLIST=%d.\n", MAXLEN_OUTPUTLIST);

              OutputListFlag[OutputListLength] = flag;
              OutputListLength++;
            }
        }

      fclose(fd);

      mpi_printf("\nfound %d times in output-list.\n", OutputListLength);

#ifndef OUTPUT_NON_SYNCHRONIZED_ALLOWED
      // now test that the output spacing is not smaller than the maximum timestep
      if(OutputListLength > 1)
        {
          double dtmin = 0;
          for(int i = 0; i < OutputListLength - 1; i++)
            {
              double dt;
              if(All.ComovingIntegrationOn)
                dt = log(OutputListTimes[i + 1]) - log(OutputListTimes[i]);
              else
                dt = OutputListTimes[i + 1] - OutputListTimes[i];

              if(i == 0)
                dtmin = dt;
              else
                dtmin = std::min<double>(dt, dtmin);
            }
          if(dtmin < All.MaxSizeTimestep)
            Terminate(
                "\nMinimal spacing in desired output list is dlna/dt = %g, which is smaller than MaxSizeTimestep=%g - this could lead "
                "to missing outputs, which you probably don't want. Better to correct this.\n",
                dtmin, All.MaxSizeTimestep);
        }
#endif
    }

  /* tell all other processes */
  MPI_Bcast(get_data_ptr(), get_data_size(), MPI_BYTE, 0, Communicator);
}

void global_data_all_processes::some_parameter_checks(void)
{
  if(MaxFilesWithConcurrentIO > NTask)
    {
      mpi_printf("NOTICE: MaxFilesWithConcurrentIO has been reduced to the number of processors\n");
      MaxFilesWithConcurrentIO = NTask;
    }

  if(MaxFilesWithConcurrentIO == 0)
    {
      mpi_printf("NOTICE: MaxFilesWithConcurrentIO has been set to be equal to the number of processors\n");
      MaxFilesWithConcurrentIO = NTask;
    }

  if(SnapFormat < 1 || SnapFormat > 3)
    Terminate("Unsupported File-Format: SnapFormat = %d\n", SnapFormat);

  if(NTask < NumFilesPerSnapshot)
    {
      mpi_printf("WARNING: Number of processors less than 'NumFilesPerSnapshot=%d' - reducing this to NumFilesPerSnapshot=%d\n",
                 NumFilesPerSnapshot, NTask);
      NumFilesPerSnapshot = NTask;
    }

  for(int i = 0; i < NTYPES; i++)
    {
      if(SofteningClassOfPartType[i] >= NSOFTCLASSES || SofteningClassOfPartType[i] < 0)
        Terminate("SofteningClassOfPartType%d  invalid (NSOFTCLASSES=%d)\n", i, NSOFTCLASSES);
    }
}
