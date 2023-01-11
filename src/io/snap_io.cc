/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file snap_io.cc
 *
 *  \brief routines for I/O of snapshot files
 */

#include "gadgetconfig.h"

#include <errno.h>
#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../fof/fof.h"
#include "../gitversion/version.h"
#include "../gravtree/gravtree.h"
#include "../io/hdf5_util.h"
#include "../io/io.h"
#include "../io/snap_io.h"
#include "../lightcone/lightcone.h"
#include "../logs/logs.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/peano.h"
#include "../src/pm/pm.h"
#include "../system/system.h"

/*!
 * \brief Function for field registering.
 *
 * For init_field arguments read the documentation of init_field.
 * Don't forget to add the new IO_FLAG to io_private.h
 */
void snap_io::init_basic(simparticles *Sp_ptr)
{
  Sp = Sp_ptr;

  this->N_IO_Fields  = 0;
  this->N_DataGroups = NTYPES + 1;  // the last data group is a tree table used only for storing/reading tree-reordered particle data
  this->header_size  = sizeof(header);
  this->header_buf   = &header;
  this->type_of_file = FILE_IS_SNAPSHOT;
  snprintf(this->info, MAXLEN_PATH, "SNAPSHOT: writing snapshot");

#ifdef OUTPUT_COORDINATES_AS_INTEGERS
  init_field("IPOS", "IntCoordinates", MEM_MY_INTPOS_TYPE, FILE_MY_INTPOS_TYPE, READ_IF_PRESENT, 3, A_P, NULL, io_func_intpos,
             ALL_TYPES, 1, 1., -1., 1., 0., 0., All.UnitLength_in_cm * Sp->FacIntToCoord, true);
#else
  init_field("POS ", "Coordinates", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_P, NULL, io_func_pos, ALL_TYPES, 1, 1., -1.,
             1., 0., 0., All.UnitLength_in_cm, true);
#endif

#ifdef OUTPUT_VELOCITIES_IN_HALF_PRECISION
  init_field("VEL ", "Velocities", MEM_MY_FLOAT, FILE_HALF, READ_IF_PRESENT, 3, A_NONE, NULL, io_func_vel, ALL_TYPES, 1, 0.5, 0., 0.,
             0., 1., All.UnitVelocity_in_cm_per_s);
#else
  init_field("VEL ", "Velocities", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_NONE, NULL, io_func_vel, ALL_TYPES, 1, 0.5,
             0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
#endif

  init_field("ID  ", "ParticleIDs", MEM_MY_ID_TYPE, FILE_MY_ID_TYPE, READ_IF_PRESENT, 1, A_P, NULL, io_func_id, ALL_TYPES, 0, 0, 0, 0,
             0, 0, 0, true);

#ifndef LEAN
  init_field("MASS", "Masses", MEM_MY_DOUBLE, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_P, NULL, io_func_mass, MASS_BLOCK, 1, 0., -1.,
             0., 1., 0., All.UnitMass_in_g);
#endif

#ifdef SECOND_ORDER_LPT_ICS
  init_field("LPTM", "SecondOrderICMasses", MEM_FLOAT, FILE_NONE, READ_IF_PRESENT, 1, A_P, &Sp->P[0].OldAcc, NULL, ALL_TYPES, 0, 0, 0,
             0, 0, 0, 0);
#endif

  init_field("U   ", "InternalEnergy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_NONE, NULL, io_func_u, GAS_ONLY, 1, 0.,
             0., 0., 0., 2., All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);

  init_field("RHO ", "Density", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].Density, NULL, GAS_ONLY, 1,
             -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);

  init_field("HSML", "SmoothingLength", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].Hsml, NULL, GAS_ONLY,
             1, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

#ifdef OUTPUT_ACCELERATION
#ifdef OUTPUT_ACCELERATIONS_IN_HALF_PRECISION
  All.accel_normalize_fac = 10.0 * All.Hubble * (100.0 * 1.0e5 / All.UnitVelocity_in_cm_per_s);

  init_field("ACCE", "Acceleration", MEM_MY_FLOAT, FILE_HALF, SKIP_ON_READ, 3, A_NONE, 0, io_func_accel, ALL_TYPES, 1, -2.0, 1, -1, 0,
             2, All.accel_normalize_fac * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#else
  All.accel_normalize_fac = 1.0;

  init_field("ACCE", "Acceleration", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 3, A_NONE, 0, io_func_accel, ALL_TYPES, 1, -2.0, 1,
             -1, 0, 2, All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm);
#endif

  /* hydro acceleration */
  init_field("HACC", "HydroAcceleration", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 3, A_SPHP, &Sp->SphP[0].HydroAccel, 0,
             GAS_ONLY, 0, 0, 0, 0, 0, 0, 0);
#endif

#ifdef STARFORMATION

  init_field("SFR ", "StarFormationRate", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, All.RestartFlag == RST_FOF ? READ_IF_PRESENT : SKIP_ON_READ,
             1, A_NONE, 0, io_func_sfr, GAS_ONLY, 1, 0., 0., -1., 1., 1., SOLAR_MASS / SEC_PER_YEAR);

  init_field("AGE ", "StellarFormationTime", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_P, &Sp->P[0].StellarAge, NULL,
             AGE_BLOCK, /* stellar formation time */
             0, 0, 0, 0, 0, 0, 0);

  init_field("Z   ", "Metallicity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_NONE, 0, io_func_metallicity,
             Z_BLOCK, /* gas and star metallicity */
             0, 0, 0, 0, 0, 0, 0);
#endif

#if defined(PRESSURE_ENTROPY_SPH) && defined(OUTPUT_PRESSURE_SPH_DENSITY)
  init_field("PRHO", "PressureSphDensity", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].PressureSphDensity,
             NULL, GAS_ONLY, /* Pressure density */
             1, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);
#endif

#ifdef OUTPUT_PRESSURE
  init_field("PRES", "Pressure", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_NONE, 0, io_func_pressure,
             GAS_ONLY, /* particle pressure */
             1, -3 * GAMMA, 2., -3., 1., 2., All.UnitDensity_in_cgs * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
#endif

#if defined(TIMEDEP_ART_VISC) && defined(OUTPUT_VISCOSITY_PARAMETER)
  init_field("ALP ", "ArtificialViscosityParameter", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].Alpha,
             NULL, GAS_ONLY, 1, -3., 2., -3., 1., 0., 1);
#endif

#ifdef OUTPUT_ENTROPY
  init_field("ENTR", "Entropy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].Entropy, 0,
             GAS_ONLY, /* particle entropy */
             0, 0, 0, 0, 0, 0, 0);
#endif

#ifdef COOLING
  init_field("NE  ", "ElectronAbundance", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].Ne, 0,
             GAS_ONLY, /* electron abundance */
             0, 0, 0, 0, 0, 0, 0);
#endif

#ifdef OUTPUT_POTENTIAL
  init_field("POT ", "Potential", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_P, &Sp->P[0].Potential, 0,
             ALL_TYPES, /* potential */
             1, -1., 0., 0., 0., 2., All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
#endif

#ifdef OUTPUT_CHANGEOFENTROPY
  init_field("ENDT", "RateOfChangeOfEntropy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].DtEntropy, 0,
             GAS_ONLY, /* particle entropy change */
             0, 0, 0, 0, 0, 0, 0);
#endif

#ifdef OUTPUT_TIMESTEP
  init_field("TSTP", "TimeStep", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_NONE, 0, io_func_timestep,
             ALL_TYPES, /* time step */
             0, 0, 0, 0, 0, 0, 0);

  init_field("TSTH", "TimeStepHydro", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, SKIP_ON_READ, 1, A_NONE, 0, io_func_timestephydro,
             GAS_ONLY, /* hydro time step */
             0, 0, 0, 0, 0, 0, 0);
#endif

#ifdef OUTPUT_DIVVEL
  init_field("DIVV", "VelocityDivergence", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].DivVel, 0,
             GAS_ONLY, /* hydro velocity divergence */
             0, 0, 0, 0, 0, 0, 0);
#endif

#ifdef OUTPUT_CURLVEL
  init_field("ROTV", "VelocityCurl", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].CurlVel, 0,
             GAS_ONLY, /* absolute value of rot v */
             0, 0, 0, 0, 0, 0, 0);
#endif

#ifdef OUTPUT_VELOCITY_GRADIENT
#ifndef IMPROVED_VELOCITY_GRADIENTS
#error "The option OUTPUT_VELOCITY_GRADIENT requires IMPROVED_VELOCITY_GRADIENTS"
#endif
  init_field("GRAV", "VelocityGradient", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 9, A_SPHP, &Sp->SphP[0].dvel[0][0], 0,
             GAS_ONLY, 0, 0, 0, 0, 0, 0, 0);
#endif

#ifdef OUTPUT_COOLHEAT
  init_field("COHE", "CoolingHeatingEnergy", MEM_MY_FLOAT, FILE_MY_IO_FLOAT, READ_IF_PRESENT, 1, A_SPHP, &Sp->SphP[0].CoolHeat, 0,
             GAS_ONLY, 0, 0, 0, 0, 0, 0, 0);
#endif

#if defined(SUBFIND) && defined(SUBFIND_STORE_LOCAL_DENSITY)

  init_field("SFDE", "SubfindDensity", MEM_MY_FLOAT, All.RestartFlag != RST_CREATEICS ? FILE_MY_IO_FLOAT : FILE_NONE, SKIP_ON_READ, 1,
             A_PS, &Sp->PS[0].SubfindDensity, 0, ALL_TYPES, /* subfind density */
             1, -3., 2., -3., 1., 0., All.UnitDensity_in_cgs);

  init_field("SFHS", "SubfindHsml", MEM_MY_FLOAT, All.RestartFlag != RST_CREATEICS ? FILE_MY_IO_FLOAT : FILE_NONE, SKIP_ON_READ, 1,
             A_PS, &Sp->PS[0].SubfindHsml, 0, ALL_TYPES, /* subfind hsml */
             1, 1., -1., 1., 0., 0., All.UnitLength_in_cm);

  init_field("SFVD", "SubfindVelDisp", MEM_MY_FLOAT, All.RestartFlag != RST_CREATEICS ? FILE_MY_IO_FLOAT : FILE_NONE, SKIP_ON_READ, 1,
             A_PS, &Sp->PS[0].SubfindVelDisp, 0, ALL_TYPES, /* subfind velocity dispersion */
             1, 0., 0., 0., 0., 1., All.UnitVelocity_in_cm_per_s);
#endif

#if defined(GADGET2_HEADER) && defined(REARRANGE_OPTION) && defined(MERGERTREE)
  Terminate("GADGET2_HEADER does not work together with REARRANGE_OPTION\n");
#endif
}

#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
void snap_io::init_extra(simparticles *Sp_ptr, mergertree *MergerTree_ptr)
{
  Sp         = Sp_ptr;
  MergerTree = MergerTree_ptr;

  init_field("MTRI", "TreeID", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_P, &Sp->P[0].TreeID, NULL, ALL_TYPES, 0, 0, 0, 0, 0, 0, 0);

  init_field("MTRL", "ParticleCount", MEM_INT, FILE_INT, READ_IF_PRESENT, 1, A_TT, &MergerTree->TreeTable[0].HaloCount, NULL,
             TREETABLE, 0, 0, 0, 0, 0, 0, 0);
  init_field("MTRS", "ParticleFirst", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_TT, &MergerTree->TreeTable[0].FirstHalo, NULL,
             TREETABLE, 0, 0, 0, 0, 0, 0, 0);
  init_field("MTRI", "TreeID", MEM_INT64, FILE_INT64, READ_IF_PRESENT, 1, A_TT, &MergerTree->TreeTable[0].TreeID, NULL, TREETABLE, 0,
             0, 0, 0, 0, 0, 0);
}
#endif

void snap_io::read_snapshot(int num, mysnaptype loc_snap_type)
{
  snap_type = loc_snap_type;

  char buf[MAXLEN_PATH_EXTRA];

  if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT)
    {
      if(All.NumFilesPerSnapshot > 1)
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s/snapdir_%03d/%s-prevmostboundonly_%03d", All.OutputDir, num, All.SnapshotFileBase, num);
      else
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s-prevmostboundonly_%03d", All.OutputDir, All.SnapshotFileBase, num);
    }
  else
    {
      if(All.NumFilesPerSnapshot > 1)
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s/snapdir_%03d/%s_%03d", All.OutputDir, num, All.SnapshotFileBase, num);
      else
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);
    }

  read_ic(buf);
}

/*! \brief This function reads initial conditions that are in on of the default file formats
 * of Gadget.
 *
 * Snapshot files can be used as input files.  However, when a
 * snapshot file is used as input, not all the information in the header is
 * used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 * Alternatively, the code can be started with All.RestartFlag 2, then snapshots
 * from the code can be used as initial conditions-files without having to
 * change the parameter file.  For gas particles, only the internal energy is
 * read, the density and mean molecular weight will be recomputed by the code.
 * When InitGasTemp>0 is given, the gas temperature will be initialized to this
 * value assuming a mean molecular weight either corresponding to complete
 * neutrality, or full ionization.
 *
 * \param fname file name of the ICs
 * \param readTypes readTypes is a bitfield that
 * determines what particle types to read, only if the bit
 * corresponding to a particle type is set, the corresponding data is
 * loaded, otherwise its particle number is set to zero. (This is
 * only implemented for HDF5 files.)
 */
void snap_io::read_ic(const char *fname)
{
  if(All.ICFormat < 1 || All.ICFormat > 4)
    Terminate("ICFormat = %d not supported.\n", All.ICFormat);

  TIMER_START(CPU_SNAPSHOT);

  double t0 = Logs.second();

  Sp->TotNumPart = 0;

  int num_files = find_files(fname, fname);

  reset_io_byte_count();

  /* we repeat reading the headers of the files two times. In the first iteration, only the
   * particle numbers ending up on each processor are assembled, followed by memory allocation.
   * In the second iteration, the data is actually read in.
   */
  for(int rep = 0; rep < 2; rep++)
    {
      Sp->NumPart = 0;
      Sp->NumGas  = 0;

      read_files_driver(fname, rep, num_files);

      /* now do the memory allocation */
      if(rep == 0)
        {
          int max_load, max_sphload;
          MPI_Allreduce(&Sp->NumPart, &max_load, 1, MPI_INT, MPI_MAX, Communicator);
          MPI_Allreduce(&Sp->NumGas, &max_sphload, 1, MPI_INT, MPI_MAX, Communicator);

#ifdef TILING
          max_load *= TILING * TILING * TILING;
          max_sphload *= TILING * TILING * TILING;
#endif

          Sp->MaxPart    = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
          Sp->MaxPartSph = max_sphload / (1.0 - 2 * ALLOC_TOLERANCE);

          Sp->allocate_memory(); /* allocate particle storage */

#ifndef OUTPUT_COORDINATES_AS_INTEGERS
          Ptmp = (ptmp_data *)Mem.mymalloc_movable(&Ptmp, "Ptmp", Sp->NumPart * sizeof(ptmp_data));
#endif
        }
    }

  MPI_Barrier(Communicator);

  long long byte_count = get_io_byte_count(), byte_count_all;
  sumup_longs(1, &byte_count, &byte_count_all, Communicator);

  double t1 = Logs.second();

  mpi_printf("READIC: reading done. Took %g sec, total size %g MB, corresponds to effective I/O rate of %g MB/sec\n",
             Logs.timediff(t0, t1), byte_count_all / (1024.0 * 1024.0), byte_count_all / (1024.0 * 1024.0) / Logs.timediff(t0, t1));

  mpi_printf("\nREADIC: Total number of particles :  %lld\n\n", Sp->TotNumPart);

  snap_init_domain_mapping();

#ifndef OUTPUT_COORDINATES_AS_INTEGERS
  Mem.myfree(Ptmp);
#endif

#ifdef TILING

  MyIntPosType halfboxlen = ((MyIntPosType)1) << ((BITS_FOR_POSITIONS - 1));

  MyIntPosType len = (halfboxlen / TILING) + (halfboxlen / TILING);

  Sp->TotNumGas *= TILING * TILING * TILING;
  Sp->TotNumPart *= TILING * TILING * TILING;

  int add_numgas = Sp->NumGas * (TILING * TILING * TILING - 1);
  /* create a gap behind the existing gas particles where we will insert the new gas particles */
  memmove(static_cast<void *>(Sp->P + Sp->NumGas + add_numgas), static_cast<void *>(Sp->P + Sp->NumGas),
          (Sp->NumPart - Sp->NumGas) * sizeof(simparticles::pdata));

  int off = 0;

  for(int i = 0; i < TILING; i++)
    for(int j = 0; j < TILING; j++)
      for(int k = 0; k < TILING; k++)
        if(i != 0 || j != 0 || k != 0)
          for(int n = 0; n < Sp->NumGas; n++)
            {
              Sp->P[Sp->NumGas + off] = Sp->P[n];

              Sp->P[Sp->NumGas + off].IntPos[0] += i * len;
              Sp->P[Sp->NumGas + off].IntPos[1] += j * len;
              Sp->P[Sp->NumGas + off].IntPos[2] += k * len;

              off++;
            }

  if(off != add_numgas)
    Terminate("this should not happen");

  Sp->NumGas += add_numgas;
  Sp->NumPart += add_numgas;

  int add_numpart = (Sp->NumPart - Sp->NumGas) * (TILING * TILING * TILING - 1);

  off = 0;

  for(int i = 0; i < TILING; i++)
    for(int j = 0; j < TILING; j++)
      for(int k = 0; k < TILING; k++)
        if(i != 0 || j != 0 || k != 0)
          for(int n = Sp->NumGas; n < Sp->NumPart; n++)
            {
              Sp->P[Sp->NumPart + off] = Sp->P[n];

              Sp->P[Sp->NumPart + off].IntPos[0] += i * len;
              Sp->P[Sp->NumPart + off].IntPos[1] += j * len;
              Sp->P[Sp->NumPart + off].IntPos[2] += k * len;

              off++;
            }

  if(off != add_numpart)
    Terminate("this should not happen");

  Sp->NumPart += add_numpart;
#endif

#ifdef GADGET2_HEADER
#ifndef INITIAL_CONDITIONS_CONTAIN_ENTROPY
  if(header.flag_entropy_instead_u)
    Terminate("Initial condition file contains entropy, but INITIAL_CONDITIONS_CONTAIN_ENTROPY is not set\n");
#else
  if(!header.flag_entropy_instead_u)
    Terminate("Initial condition file contains uthermal, but INITIAL_CONDITIONS_CONTAIN_ENTROPY is set\n");
#endif
#endif

  TIMER_STOP(CPU_SNAPSHOT);
}

/*! This routine initializes the global domain mapping between integer coordinates and real space coordinates.
 *  If periodic is on, the extent is the box size. Otherwise we look at the maximum extent of the particles.
 */
void snap_io::snap_init_domain_mapping(void)
{
#ifdef PERIODIC

  Sp->RegionLen     = All.BoxSize;
  Sp->FacCoordToInt = pow(2.0, BITS_FOR_POSITIONS) / Sp->RegionLen;
  Sp->FacIntToCoord = Sp->RegionLen / pow(2.0, BITS_FOR_POSITIONS);

#if defined(NGENIC) && !defined(CREATE_GRID)
  if(All.RestartFlag == RST_BEGIN || All.RestartFlag == RST_CREATEICS)
    {
      // Make sure that the velocities are zero when a glass file is fed to IC creation
      mpi_printf("READIC: Setting velocities in glass file to zero.\n");
      for(int i = 0; i < Sp->NumPart; i++)
        for(int k = 0; k < 3; k++)
          Sp->P[i].Vel[k] = 0;
    }
#endif

#else

  double posmin[3], posmax[3];
  for(int k = 0; k < 3; k++)
    {
      posmin[k] = MAX_REAL_NUMBER;
      posmax[k] = -MAX_REAL_NUMBER;
    }

  for(int i = 0; i < Sp->NumPart; i++)
    for(int k = 0; k < 3; k++)
      {
        if(Ptmp[i].Pos[k] < posmin[k])
          posmin[k] = Ptmp[i].Pos[k];

        if(Ptmp[i].Pos[k] > posmax[k])
          posmax[k] = Ptmp[i].Pos[k];
      }

  double xyz[6] = {posmin[0], posmin[1], posmin[2], -posmax[0], -posmax[1], -posmax[2]};
  double xyz_glob[6];

  MPI_Allreduce(xyz, xyz_glob, 6, MPI_DOUBLE, MPI_MIN, Communicator);

  mpi_printf("READIC: Region covered with particles: (%g %g %g) -> (%g %g %g)\n", xyz_glob[0], xyz_glob[1], xyz_glob[2], -xyz_glob[3],
             -xyz_glob[4], -xyz_glob[5]);

  Sp->RegionLen = 0;
  for(int j = 0; j < 3; j++)
    if(-xyz_glob[j + 3] - xyz_glob[j] > Sp->RegionLen)
      Sp->RegionLen = -xyz_glob[j + 3] - xyz_glob[j];

  Sp->RegionLen *= 4.0;

  mpi_printf("READIC: Initial root node size: %g\n", Sp->RegionLen);

  for(int j = 0; j < 3; j++)
    {
      Sp->RegionCenter[j] = 0.5 * (xyz_glob[j] - xyz_glob[j + 3]);
      Sp->RegionCorner[j] = Sp->RegionCenter[j] - 0.5 * Sp->RegionLen;
    }

  Sp->FacCoordToInt = pow(2.0, BITS_FOR_POSITIONS) / Sp->RegionLen;
  Sp->FacIntToCoord = Sp->RegionLen / pow(2.0, BITS_FOR_POSITIONS);

#endif

#ifndef OUTPUT_COORDINATES_AS_INTEGERS
  for(int i = 0; i < Sp->NumPart; i++)
    Sp->pos_to_intpos(Ptmp[i].Pos, Sp->P[i].IntPos); /* converts floating point representation to integers */
#endif
}

int snap_io::get_type_of_element(int index)
{
  if(snap_type == NORMAL_SNAPSHOT)
    {
      return Sp->P[index].getType();
    }
  else if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT)
    {
      if(Sp->P[index].ID.is_previously_most_bound())
        return Sp->P[index].getType();
      else
        return -1;  // marks that particle has type that is not written out
    }
  else if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT_REORDERED)
    {
      return Sp->P[index].getType();
    }
  else
    Terminate("can't be");
}

void snap_io::set_type_of_element(int index, int type)
{
  if(type < NTYPES)
    Sp->P[index].setType(type);
}

void *snap_io::get_base_address_of_structure(enum arrays array, int index)
{
  switch(array)
    {
      case A_SPHP:
        return (void *)(Sp->SphP + index);
      case A_P:
        return (void *)(Sp->P + index);
      case A_PS:
        return (void *)(Sp->PS + index);
#ifdef LIGHTCONE_PARTICLES
      case A_LC:
        Terminate("a");  // return (void *) (Lp->P + index);
#endif
#ifdef LIGHTCONE_MASSMAPS
      case A_MM:
        Terminate("b");  // return (void *) (Lp->P + index);
#endif
#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
      case A_TT:
        return (void *)(MergerTree->TreeTable + index);
#endif
      default:
        Terminate("we don't expect to get here");
    }

  return NULL;
}

/*! \brief Save snapshot to disk
 *
 * This function writes a snapshot of the particle distribution to one or
 * several files. If NumFilesPerSnapshot>1, the snapshot is distributed
 * into several files, which are written simultaneously. Each file contains
 * data from a group of processors of size roughly NTask/NumFilesPerSnapshot.
 * \param num the snapshot number
 */
void snap_io::write_snapshot(int num, mysnaptype loc_snap_type)
{
#ifdef DO_NOT_PRODUCE_BIG_OUTPUT
  if(snap_type != MOST_BOUND_PARTICLE_SNAPHOT)
    {
      mpi_printf("\nSNAPSHOT: We skip writing snapshot file #%d @ time %g\n", num, All.Time);
      return;
    }
#endif

  snap_type = loc_snap_type;

  TIMER_START(CPU_SNAPSHOT);

  mpi_printf("\nSNAPSHOT: writing snapshot file #%d @ time %g ... \n", num, All.Time);

  double t0 = Logs.second();
  reset_io_byte_count();

  All.set_cosmo_factors_for_current_time();

  int n_type[NTYPES];

  /* determine global and local particle numbers */
  for(int n = 0; n < NTYPES; n++)
    n_type[n] = 0;

  for(int n = 0; n < Sp->NumPart; n++)
    {
      int type = get_type_of_element(n);
      if(type >= 0)
        n_type[type]++;
    }

  sumup_large_ints(NTYPES, n_type, ntot_type_all, Communicator);

  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          char buf[MAXLEN_PATH_EXTRA];
          snprintf(buf, MAXLEN_PATH_EXTRA, "%s/snapdir_%03d", All.OutputDir, num);
          mkdir(buf, 02755);
        }
      MPI_Barrier(Communicator);
    }

  char buf[MAXLEN_PATH_EXTRA];
  if(All.NumFilesPerSnapshot > 1)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s/snapdir_%03d/%s_%03d", All.OutputDir, num, All.SnapshotFileBase, num);
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);

  if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT)
    {
      if(All.NumFilesPerSnapshot > 1)
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s/snapdir_%03d/%s-prevmostboundonly_%03d", All.OutputDir, num, All.SnapshotFileBase, num);
      else
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s-prevmostboundonly_%03d", All.OutputDir, All.SnapshotFileBase, num);
    }
  else if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT_REORDERED)
    {
      if(All.NumFilesPerSnapshot > 1)
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s/snapdir_%03d/%s-prevmostboundonly-treeorder_%03d", All.OutputDir, num,
                 All.SnapshotFileBase, num);
      else
        snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s-prevmostboundonly-treeorder_%03d", All.OutputDir, All.SnapshotFileBase, num);
    }

  /* now write the files */
  write_multiple_files(buf, All.NumFilesPerSnapshot);

  long long byte_count = get_io_byte_count(), byte_count_all;
  sumup_longs(1, &byte_count, &byte_count_all, Communicator);

  double t1 = Logs.second();

  mpi_printf("SNAPSHOT: done with writing snapshot.  Took %g sec, total size %g MB, corresponds to effective I/O rate of %g MB/sec\n",
             Logs.timediff(t0, t1), byte_count_all / (1024.0 * 1024.0), byte_count_all / (1024.0 * 1024.0) / Logs.timediff(t0, t1));

  All.Ti_lastoutput = All.Ti_Current;

  TIMER_STOP(CPU_SNAPSHOT);
}

void snap_io::fill_file_header(int writeTask, int lastTask, long long *n_type, long long *ntot_type)
{
  /* determine global and local particle numbers */
  for(int n = 0; n < NTYPES + 1; n++)
    n_type[n] = 0;

  for(int n = 0; n < Sp->NumPart; n++)
    {
      int type = get_type_of_element(n);
      if(type >= 0)
        n_type[type]++;
    }

#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
  if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT_REORDERED)
    n_type[NTYPES] = MergerTree->Ntrees;
#endif

  /* determine particle numbers of each type in file */
  if(ThisTask == writeTask)
    {
      for(int n = 0; n < NTYPES + 1; n++)
        ntot_type[n] = n_type[n];

      for(int task = writeTask + 1; task <= lastTask; task++)
        {
          long long nn[NTYPES + 1];
          MPI_Recv(&nn[0], NTYPES + 1, MPI_LONG_LONG, task, TAG_LOCALN, Communicator, MPI_STATUS_IGNORE);
          for(int n = 0; n < NTYPES + 1; n++)
            ntot_type[n] += nn[n];
        }

      for(int task = writeTask + 1; task <= lastTask; task++)
        MPI_Send(&ntot_type[0], NTYPES + 1, MPI_LONG_LONG, task, TAG_N, Communicator);
    }
  else
    {
      MPI_Send(&n_type[0], NTYPES + 1, MPI_LONG_LONG, writeTask, TAG_LOCALN, Communicator);
      MPI_Recv(&ntot_type[0], NTYPES + 1, MPI_LONG_LONG, writeTask, TAG_N, Communicator, MPI_STATUS_IGNORE);
    }

  /* fill file header */

  for(int n = 0; n < NTYPES; n++)
    {
      header.npart[n]      = ntot_type[n];
      header.npartTotal[n] = ntot_type_all[n];
    }

#if defined(MERGERTREE) && !defined(GADGET2_HEADER)
  if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT_REORDERED)
    {
      header.Ntrees = ntot_type[NTYPES];
#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
      header.TotNtrees = MergerTree->TotNtrees;
#endif
    }
  else
    {
      header.Ntrees    = 0;
      header.TotNtrees = 0;
    }
#endif

  for(int n = 0; n < NTYPES; n++)
    header.mass[n] = All.MassTable[n];

  header.time = All.Time;

  if(All.ComovingIntegrationOn)
    header.redshift = 1.0 / All.Time - 1;
  else
    header.redshift = 0;

  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize   = All.BoxSize;

#ifdef GADGET2_HEADER
  for(int n = 0; n < NTYPES; n++)
    header.npartTotalLowWord[n] = ntot_type_all[n];

  header.flag_sfr      = 0;
  header.flag_feedback = 0;
  header.flag_cooling  = 0;

#ifdef COOLING
  header.flag_cooling = 1;
#endif

#ifdef STARFORMATION
  header.flag_sfr      = 1;
  header.flag_feedback = 1;
#endif
  header.Omega0      = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;

#if !(defined(REARRANGE_OPTION) && defined(MERGERTREE))
  header.HubbleParam = All.HubbleParam;
  header.Hubble      = All.Hubble;
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
  header.flag_doubleprecision = 1;
#else
  header.flag_doubleprecision = 0;
#endif
#endif
}

void snap_io::read_file_header(const char *fname, int filenr, int readTask, int lastTask, long long *n_type, long long *ntot_type,
                               int *nstart)
{
  n_type[NTYPES]    = 0;
  ntot_type[NTYPES] = 0;

#ifdef GADGET2_HEADER
  for(int i = 0; i < NTYPES_HEADER; i++)
    if(header.npartTotalLowWord[i] > 0)
      header.npartTotal[i] = header.npartTotalLowWord[i];  // + (((long long)header.npartTotalHighWord[i]) << 32);
#endif

  if(Sp->TotNumPart == 0)
    {
      if(header.num_files <= 1)
        for(int i = 0; i < NTYPES; i++)
          header.npartTotal[i] = header.npart[i];

      Sp->TotNumGas  = header.npartTotal[0];
      Sp->TotNumPart = 0;

      for(int i = 0; i < NTYPES; i++)
        Sp->TotNumPart += header.npartTotal[i];

#ifdef GADGET2_HEADER
#if defined(SECOND_ORDER_LPT_ICS)
      if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
        {
          mpi_printf("READIC:  Second order ICs detected. Will complete them first before starting run.\n");
          All.LptScalingfactor = header.lpt_scalingfactor;
        }
      else
        Terminate("READIC:  No second order ICs detected even though you activated SECOND_ORDER_LPT_ICS.\n");
#else
      if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
        Terminate("Detected second order ICs but SECOND_ORDER_LPT_ICS is not enabled.\n");
#endif
#endif

#ifdef GENERATE_GAS_IN_ICS
      if(All.RestartFlag == RST_BEGIN)
        {
#ifdef SPLIT_PARTICLE_TYPE
          for(int i = 0; i < NTYPES; i++)
            if((1 << i) & (SPLIT_PARTICLE_TYPE))
              {
                Sp->TotNumGas += header.npartTotal[i];
                Sp->TotNumPart += header.npartTotal[i];
              }
#else
          Sp->TotNumGas += header.npartTotal[1];
          Sp->TotNumPart += header.npartTotal[1];
#endif
        }
#endif

      for(int i = 0; i < NTYPES; i++)
        All.MassTable[i] = header.mass[i];

      if(All.RestartFlag == RST_BEGIN || All.RestartFlag == RST_RESUME || All.RestartFlag == RST_CREATEICS)
        All.Time = All.TimeBegin;
      else
        All.Time = All.TimeBegin = header.time;

      All.set_cosmo_factors_for_current_time();
    }

  if(ThisTask == readTask)
    {
      if(filenr == 0 && nstart == NULL)
        {
          mpi_printf(
              "\nREADIC: filenr=%d, '%s' contains:\n"
              "READIC: Type 0 (gas):   %8lld  (tot=%15lld) masstab= %g\n",
              filenr, fname, (long long)header.npart[0], (long long)header.npartTotal[0], All.MassTable[0]);

          for(int type = 1; type < NTYPES; type++)
            {
              mpi_printf("READIC: Type %d:         %8lld  (tot=%15lld) masstab= %g\n", type, (long long)header.npart[type],
                         (long long)header.npartTotal[type], All.MassTable[type]);
            }
          mpi_printf("\n");
        }
    }

  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

  long long nall = 0;
  for(int type = 0; type < NTYPES; type++)
    {
      ntot_type[type] = header.npart[type];

      long long n_in_file = header.npart[type];
      int ntask           = lastTask - readTask + 1;
      int n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
        n_for_this_task++;

      n_type[type] = n_for_this_task;

      nall += n_for_this_task;
    }

#if defined(MERGERTREE) && !defined(GADGET2_HEADER)
  if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT_REORDERED)
    ntot_type[NTYPES] = header.Ntrees;
#endif

  if(nstart)
    {
      memmove(static_cast<void *>(&Sp->P[Sp->NumGas + nall]), static_cast<void *>(&Sp->P[Sp->NumGas]),
              (Sp->NumPart - Sp->NumGas) * sizeof(particle_data));
#ifndef OUTPUT_COORDINATES_AS_INTEGERS
      memmove(&Ptmp[Sp->NumGas + nall], &Ptmp[Sp->NumGas], (Sp->NumPart - Sp->NumGas) * sizeof(ptmp_data));
#endif
      *nstart = Sp->NumGas;
    }
}

/*! \brief Write the fields contained in the header group of the HDF5 snapshot file
 *
 *  This function stores the fields of the structure io_header as attributes belonging
 *  to the header group of the HDF5 file.
 *
 *  \param handle contains a reference to the header group
 */
void snap_io::write_header_fields(hid_t handle)
{
#ifdef GADGET2_HEADER
  write_vector_attribute(handle, "NumPart_ThisFile", header.npart, H5T_NATIVE_UINT, NTYPES);
#else
  write_vector_attribute(handle, "NumPart_ThisFile", header.npart, H5T_NATIVE_UINT64, NTYPES);
#endif
  write_vector_attribute(handle, "NumPart_Total", header.npartTotal, H5T_NATIVE_UINT64, NTYPES);
  write_vector_attribute(handle, "MassTable", header.mass, H5T_NATIVE_DOUBLE, NTYPES);
  write_scalar_attribute(handle, "Time", &header.time, H5T_NATIVE_DOUBLE);
  write_scalar_attribute(handle, "Redshift", &header.redshift, H5T_NATIVE_DOUBLE);
  write_scalar_attribute(handle, "BoxSize", &header.BoxSize, H5T_NATIVE_DOUBLE);
  write_scalar_attribute(handle, "NumFilesPerSnapshot", &header.num_files, H5T_NATIVE_INT);
  write_string_attribute(handle, "Git_commit", GIT_COMMIT);
  write_string_attribute(handle, "Git_date", GIT_DATE);

#if defined(REARRANGE_OPTION) && defined(MERGERTREE)
  if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT_REORDERED)
    {
      write_scalar_attribute(handle, "Ntrees_ThisFile", &header.Ntrees, H5T_NATIVE_UINT64);
      write_scalar_attribute(handle, "Ntrees_Total", &header.TotNtrees, H5T_NATIVE_UINT64);
    }
#endif
}

/*! \brief This function reads the snapshot header in case of hdf5 files (i.e. format 3)
 *
 * \param fname file name of the snapshot as given in the parameter file
 */
void snap_io::read_header_fields(const char *fname)
{
  for(int i = 0; i < NTYPES; i++)
    {
      header.npart[i]      = 0;
      header.npartTotal[i] = 0;
      header.mass[i]       = 0;
    }

  int ntypes = NTYPES;

  hid_t hdf5_file = my_H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t handle    = my_H5Gopen(hdf5_file, "/Header");

  /* check if the file in question actually has this number of types */
  hid_t hdf5_attribute = my_H5Aopen_name(handle, "NumPart_ThisFile");
  hid_t space          = H5Aget_space(hdf5_attribute);
  hsize_t dims, len;
  H5Sget_simple_extent_dims(space, &dims, &len);
  H5Sclose(space);
  if(len != (size_t)ntypes)
    Terminate("Length of NumPart_ThisFile attribute (%d) does not match NTYPES(ICS) (%d)", (int)len, ntypes);
  my_H5Aclose(hdf5_attribute, "NumPart_ThisFile");

  /* now read the header fields */

#ifdef GADGET2_HEADER
  read_vector_attribute(handle, "NumPart_ThisFile", header.npart, H5T_NATIVE_UINT, ntypes);
#else
  read_vector_attribute(handle, "NumPart_ThisFile", header.npart, H5T_NATIVE_UINT64, ntypes);
#endif

  read_vector_attribute(handle, "NumPart_Total", header.npartTotal, H5T_NATIVE_UINT64, ntypes);

#if defined(MERGERTREE) && !defined(GADGET2_HEADER)
  if(snap_type == MOST_BOUND_PARTICLE_SNAPHOT_REORDERED)
    {
      read_scalar_attribute(handle, "Ntrees_ThisFile", &header.Ntrees, H5T_NATIVE_UINT64);
      read_scalar_attribute(handle, "Ntrees_Total", &header.TotNtrees, H5T_NATIVE_UINT64);
    }
#endif

  read_vector_attribute(handle, "MassTable", header.mass, H5T_NATIVE_DOUBLE, ntypes);
  read_scalar_attribute(handle, "Time", &header.time, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "Redshift", &header.redshift, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "BoxSize", &header.BoxSize, H5T_NATIVE_DOUBLE);
  read_scalar_attribute(handle, "NumFilesPerSnapshot", &header.num_files, H5T_NATIVE_INT);

#if defined(GADGET2_HEADER) && defined(SECOND_ORDER_LPT_ICS)
  read_scalar_attribute(handle, "LptScalingfactor", &header.lpt_scalingfactor, H5T_NATIVE_FLOAT);
#endif

  my_H5Gclose(handle, "/Header");
  my_H5Fclose(hdf5_file, fname);
}

int snap_io::get_filenr_from_header(void) { return header.num_files; }

void snap_io::set_filenr_in_header(int numfiles) { header.num_files = numfiles; }

void snap_io::read_increase_numbers(int type, int n_for_this_task)
{
  if(type < NTYPES)
    Sp->NumPart += n_for_this_task;

  if(type == 0)
    Sp->NumGas += n_for_this_task;
}

void snap_io::get_datagroup_name(int type, char *buf)
{
  if(type < NTYPES)
    snprintf(buf, MAXLEN_PATH, "/PartType%d", type);
  else if(type == NTYPES)
    snprintf(buf, MAXLEN_PATH, "/TreeTable");
  else
    Terminate("wrong group");
}
