/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file init.cc
 *
 *  \brief code for initialization of a simulation from initial conditions
 */

// clang-format off
#include "gadgetconfig.h"
// clang-format on

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../cooling_sfr/cooling.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fof/fof.h"
#include "../gravtree/gravtree.h"
#include "../io/io.h"
#include "../io/snap_io.h"
#include "../logs/timer.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../ngenic/ngenic.h"
#include "../pm/pm.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind_readid_io.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

using namespace std;

/*! \brief Prepares the loaded initial conditions for the run
 *
 *  It is only called if RestartFlag != RST_RESUME. Various counters and variables are initialized.
 *  Entries of the particle data structures not read from initial conditions are
 *  initialized or converted and a initial domain decomposition is performed.
 *  If SPH particles are present, the initial SPH smoothing lengths are determined.
 */
void sim::init(int RestartSnapNum)
{
#ifdef NGENIC
  if(All.RestartFlag == RST_CREATEICS || All.RestartFlag == RST_BEGIN)
    {
      Ngenic.ngenic_displace_particles();

      if(All.RestartFlag == RST_CREATEICS)
        {
          double fac = 1 / sqrt(All.cf_a3inv);
          for(int i = 0; i < Sp.NumPart; i++)
            for(int k = 0; k < 3; k++)
              Sp.P[i].Vel[k] *= fac;

          strcat(All.SnapshotFileBase, "_ics");
          mpi_printf("Start writing file %s\nRestartSnapNum %d\n", All.SnapshotFileBase, 0);
          snap_io Snap(&Sp, Communicator, All.SnapFormat); /* get an I/O object */
          Snap.write_snapshot(0, NORMAL_SNAPSHOT);
          endrun();
        }
    }
#else
  if(All.RestartFlag == RST_CREATEICS)
    {
      Terminate("Compile with option NGENIC to create cosmological initial conditions");
    }
#endif

#ifdef LIGHTCONE_PARTICLES
  Lp.MaxPart = LIGHTCONE_ALLOC_FAC * Sp.MaxPart;
  Lp.NumPart = 0;
  Lp.allocate_memory();
#endif

#ifdef LIGHTCONE_MASSMAPS
  LightCone.Mp->Npix = nside2npix(All.LightConeMassMapsNside);

  subdivide_evenly(LightCone.Mp->Npix, NTask, ThisTask, &LightCone.Mp->FirstPix, &LightCone.Mp->NpixLoc);

  Mp.MaxPart = LIGHTCONE_MASSMAP_ALLOC_FAC * (Sp.TotNumPart / NTask);
  Mp.NumPart = 0;
  Mp.allocate_memory();
  LightCone.MassMap = (double *)Mem.mymalloc_movable_clear(&LightCone.MassMap, "MassMap", LightCone.Mp->NpixLoc * sizeof(double));
#endif

  /* this makes sure that masses are initialized in the case that the mass-block
     is empty for this particle type */

  for(int i = 0; i < Sp.NumPart; i++)
    if(All.MassTable[Sp.P[i].getType()] != 0)
      {
#ifndef LEAN
        Sp.P[i].setMass(All.MassTable[Sp.P[i].getType()]);
#else
        All.PartMass = All.MassTable[Sp.P[i].getType()];
#endif
      }

#if NSOFTCLASSES > 1
  for(int i = 0; i < Sp.NumPart; i++)
    Sp.P[i].setSofteningClass(All.SofteningClassOfPartType[Sp.P[i].getType()]);
#endif

#ifdef GENERATE_GAS_IN_ICS
  if(All.RestartFlag == RST_BEGIN)
    {
      /* determine maximum ID */
      MyIDType maxid = 0;
      for(int i = 0; i < Sp.NumPart; i++)
        if(Sp.P[i].ID.get() > maxid)
          maxid = Sp.P[i].ID.get();

      MyIDType *tmp = (MyIDType *)Mem.mymalloc("tmp", NTask * sizeof(MyIDType));
      MPI_Allgather(&maxid, sizeof(MyIDType), MPI_BYTE, tmp, sizeof(MyIDType), MPI_BYTE, Communicator);

      for(int i = 0; i < NTask; i++)
        if(tmp[i] > maxid)
          maxid = tmp[i];

      Mem.myfree(tmp);

      int count = 0;
      for(int i = 0; i < Sp.NumPart; i++)
#ifdef SPLIT_PARTICLE_TYPE
        if((1 << Sp.P[i].getType()) & (SPLIT_PARTICLE_TYPE))
#else
        if(Sp.P[i].getType() == 1)
#endif
          count++;

      int *numpart_list = (int *)Mem.mymalloc("numpart_list", NTask * sizeof(int));
      MPI_Allgather(&count, 1, MPI_INT, numpart_list, 1, MPI_INT, Communicator);

      maxid++;

      for(int i = 0; i < ThisTask; i++)
        maxid += numpart_list[i];

      Mem.myfree(numpart_list);

      Domain.domain_resize_storage(count + Sp.NumPart, count, 0);

      memmove(Sp.P + count, Sp.P, sizeof(particle_data) * Sp.NumPart);

      Sp.NumPart += count;
      Sp.NumGas += count;

      if(Sp.NumGas > Sp.MaxPartSph)
        Terminate("Task=%d ends up getting more SPH particles (%d) than allowed (%d)\n", ThisTask, Sp.NumGas, Sp.MaxPartSph);

      if(Sp.NumPart > Sp.MaxPart)
        Terminate("Task=%d ends up getting more particles (%d) than allowed (%d)\n", ThisTask, Sp.NumPart, Sp.MaxPart);

      double fac = All.OmegaBaryon / All.Omega0;
      double rho = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      int j = 0;
      for(int i = count; i < Sp.NumPart; i++)
#ifdef SPLIT_PARTICLE_TYPE
        if((1 << Sp.P[i].getType()) & (SPLIT_PARTICLE_TYPE))
#else
        if(Sp.P[i].getType() == 1)
#endif
          {
            double d = pow(Sp.P[i].getMass() / rho, 1.0 / 3);
            double a = 0.5 * All.OmegaBaryon / All.Omega0 * d;
            double b = 0.5 * (All.Omega0 - All.OmegaBaryon) / All.Omega0 * d;

            MyIntPosType delta_a[3];
            double aa[3] = {a, a, a};
            Sp.pos_to_signedintpos(aa, (MySignedIntPosType *)delta_a);

            MyIntPosType delta_b[3];
            double bb[3] = {b, b, b};
            Sp.pos_to_signedintpos(bb, (MySignedIntPosType *)delta_b);

            Sp.P[j] = Sp.P[i];

            Sp.P[j].setMass(Sp.P[j].getMass() * fac);
            Sp.P[i].setMass(Sp.P[i].getMass() * (1 - fac));

#ifdef SECOND_ORDER_LPT_ICS
            Sp.P[j].OldAcc *= fac;
            Sp.P[i].OldAcc *= (1 - fac);
#endif

            Sp.P[j].setType(0);
#if NSOFTCLASSES > 1
            Sp.P[j].setSofteningClass(All.SofteningClassOfPartType[0]);
#endif
            Sp.P[j].ID.set(maxid++);
            Sp.P[i].IntPos[0] += delta_a[0];
            Sp.P[i].IntPos[1] += delta_a[1];
            Sp.P[i].IntPos[2] += delta_a[2];
            Sp.P[j].IntPos[0] -= delta_b[0];
            Sp.P[j].IntPos[1] -= delta_b[1];
            Sp.P[j].IntPos[2] -= delta_b[2];

            j++;
          }

      All.MassTable[0] = 0;

#ifdef SPLIT_PARTICLE_TYPE
      for(int i = 1; i < NTYPES; i++)
        if((1 << i) & (SPLIT_PARTICLE_TYPE))
          All.MassTable[i] *= (1 - fac);
#else
      All.MassTable[1] *= (1 - fac);
#endif

      mpi_printf("\nGENERATE_GAS_IN_ICS: Generated gas particles from DM particle distribution.  TotNumGas=%lld\n\n", Sp.TotNumGas);
    }
#endif

#ifdef STARFORMATION
  if(All.RestartFlag == RST_BEGIN)
    {
      if(All.MassTable[STAR_TYPE] == 0 && All.MassTable[0] > 0)
        {
          All.MassTable[0] = 0;
        }
    }
#endif

  double u_init = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.InitGasTemp;
  u_init *= All.UnitMass_in_g / All.UnitEnergy_in_cgs; /* unit conversion */

  double molecular_weight;
  if(All.InitGasTemp > 1.0e4) /* assuming FULL ionization */
    molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  else /* assuming NEUTRAL GAS */
    molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);

  u_init /= molecular_weight;

  All.InitGasU = u_init;

  if(All.RestartFlag == RST_BEGIN)
    {
      if(All.InitGasTemp > 0)
        {
          for(int i = 0; i < Sp.NumGas; i++)
            {
              if(ThisTask == 0 && i == 0 && Sp.SphP[i].Entropy == 0)
                mpi_printf("READIC: Initializing u from InitGasTemp !\n");

              if(Sp.SphP[i].Entropy == 0)
                Sp.SphP[i].Entropy = All.InitGasU;
              /* Note: the coversion to entropy will be done in the function init(),
                 after the densities have been computed */
            }
        }
    }

  for(int i = 0; i < Sp.NumGas; i++)
    Sp.SphP[i].Entropy = std::max<double>(All.MinEgySpec, Sp.SphP[i].Entropy);

#ifdef COOLING
  CoolSfr.IonizeParams();
#endif

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current        = 0;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current        = 0;
    }

  All.set_cosmo_factors_for_current_time();

  All.NumCurrentTiStep  = 0; /* setup some counters */
  All.SnapshotFileCount = 0;

  if(All.RestartFlag == RST_STARTFROMSNAP)
    {
      if(RestartSnapNum < 0)
        All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;
      else
        All.SnapshotFileCount = RestartSnapNum + 1;
    }

  All.TotNumOfForces     = 0;
  All.TotNumDirectForces = 0;
  All.TotNumDensity      = 0;
  All.TotNumHydro        = 0;

  All.TopNodeAllocFactor = 0.08;
  All.TreeAllocFactor    = 0.3;
  All.NgbTreeAllocFactor = 0.7;

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

#if defined(EVALPOTENTIAL) && defined(PMGRID) && defined(PERIODIC)
  double mass_sum = 0;

  for(int i = 0; i < Sp.NumPart; i++)
    mass_sum += Sp.P[i].getMass();

  MPI_Allreduce(&mass_sum, &All.TotalMass, 1, MPI_DOUBLE, MPI_SUM, Communicator);
#endif

  if(All.ComovingIntegrationOn) /*  change to new velocity variable */
    {
      double afac = sqrt(All.Time) * All.Time;

      for(int i = 0; i < Sp.NumPart; i++)
        {
          for(int j = 0; j < 3; j++)
            Sp.P[i].Vel[j] *= afac; /* for dm/gas particles, p = a^2 xdot */
        }
    }

  for(int i = 0; i < TIMEBINS; i++)
    All.Ti_begstep[i] = 0;

#if defined(PMGRID) && !defined(TREEPM_NOTIMESPLIT)
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

  for(int i = 0; i < Sp.NumPart; i++) /*  start-up initialization with non-zero values where required */
    {
#ifndef LEAN
      Sp.P[i].access.clear();
#endif
#ifdef MERGERTREE
      Sp.P[i].PrevSubhaloNr.set(HALONR_MAX);
#endif
    }

  for(int i = 0; i < TIMEBINS; i++)
    Sp.TimeBinSynchronized[i] = 1;

  Sp.reconstruct_timebins();

  for(int i = 0; i < Sp.NumGas; i++) /* initialize sph_properties with non-zero values where required */
    {
      Sp.SphP[i].EntropyPred = Sp.SphP[i].Entropy;

      for(int j = 0; j < 3; j++)
        Sp.SphP[i].VelPred[j] = Sp.P[i].Vel[j];

      if(All.RestartFlag == RST_BEGIN)
        {
#ifdef COOLING
          Sp.SphP[i].Ne = 1.0;
#endif
        }

#ifdef PRESSURE_ENTROPY_SPH
      Sp.SphP[i].EntropyToInvGammaPred = pow(Sp.SphP[i].EntropyPred, 1.0 / GAMMA);
#endif

#ifdef TIMEDEP_ART_VISC
#ifdef HIGH_ART_VISC_START
      Sp.SphP[i].Alpha = All.ArtBulkViscConst;
#else
      Sp.SphP[i].Alpha = All.AlphaMin;
#endif
#endif
    }

#ifdef RECREATE_UNIQUE_IDS
  recreate_unique_ids();
#endif

  test_id_uniqueness();

  Domain.domain_decomposition(STANDARD); /* do initial domain decomposition (gives equal numbers of particles) */

  GravTree.set_softenings();

#ifdef ADAPTIVE_HYDRO_SOFTENING
  mpi_printf("INIT: Adaptive hydro softening, minimum gravitational softening for SPH particles: %g\n",
             All.MinimumComovingHydroSoftening);
  mpi_printf("INIT: Adaptive hydro softening, maximum gravitational softening for SPH particles: %g\n",
             All.MinimumComovingHydroSoftening * pow(All.AdaptiveHydroSofteningSpacing, NSOFTCLASSES_HYDRO - 1));
  mpi_printf("INIT: Adaptive hydro softening, number of softening values: %d\n", NSOFTCLASSES_HYDRO);
#endif

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  Sp.init_individual_softenings();
#endif

  if(All.RestartFlag == RST_FOF)
    {
#ifdef FOF

#if defined(SUBFIND) && defined(MERGERTREE)
      // we are reading the previous subhalo catalogue, if available, to assign the previous subhalo length to particles
      MergerTree.get_previous_size_of_subhalo_for_each_particle(RestartSnapNum - 1);
#endif

      Sp.PS = (subfind_data *)Mem.mymalloc_movable(&Sp.PS, "PS", Sp.MaxPart * sizeof(subfind_data));
      memset(Sp.PS, 0, Sp.MaxPart * sizeof(subfind_data));

      for(int i = 0; i < Sp.NumGas; i++)
        {
          if(ThisTask == 0 && i == 0)
            printf("INIT: Converting u -> entropy   All.cf_a3inv=%g\n", All.cf_a3inv);

          Sp.SphP[i].Entropy = GAMMA_MINUS1 * Sp.SphP[i].Entropy / pow(Sp.SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1);
        }

      /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
      for(int i = 0; i < Sp.NumPart; i++)
        {
          Sp.PS[i].OriginTask  = ThisTask;
          Sp.PS[i].OriginIndex = i;
        }
      fof<simparticles> FoF{Communicator, &Sp, &Domain};
      FoF.fof_fof(RestartSnapNum, "fof", "groups", 0);

      {
        All.DumpFlag_nextoutput = 1;
        snap_io Snap(&Sp, Communicator, All.SnapFormat); /* get an I/O object */
        Snap.write_snapshot(RestartSnapNum, NORMAL_SNAPSHOT);
      }

#ifdef SUBFIND_ORPHAN_TREATMENT
      {
        /* now read the IDs of the most bound particles of a previously existing special dump, with the idea being to
         * be able also on postprocessing to construct a cumulative list of all particles that used to be a most particles
         * in any previous dump. This can be accomplished by computing the group catalogues consecutively from 000 to the last
         * snapshoyt number
         */

        if(RestartSnapNum > 0)
          {
            subreadid_io SnapIDread(&Sp.IdStore, Communicator, All.SnapFormat);
            SnapIDread.previously_bound_read_snap_ids(RestartSnapNum - 1);

            FoF.subfind_match_ids_of_previously_most_bound_ids(&Sp);
          }

        snap_io Snap(&Sp, Communicator, All.SnapFormat);
        Snap.write_snapshot(RestartSnapNum, MOST_BOUND_PARTICLE_SNAPHOT); /* write special snapshot file */
      }
#endif

#endif
      endrun();
    }

  /* build neighbor tree */
  NgbTree.treeallocate(Sp.NumGas, &Sp, &Domain);
  NgbTree.treebuild(Sp.NumGas, NULL);

  if(All.RestartFlag == RST_POWERSPEC)
    {
#if defined(PMGRID) && defined(PERIODIC)

      PM.calculate_power_spectra(RestartSnapNum);
#else
      mpi_printf("\nThis option (Power Spectrum) only works for PERIODIC and enabled PMGRID.\n\n");
#endif
      endrun();
    }

  All.Ti_Current = 0;

  setup_smoothinglengths();

  /* at this point, the entropy variable actually contains the
   * internal energy, read in from the initial conditions file.
   * Once the density has been computed, we can convert to entropy.
   */
#ifdef PRESSURE_ENTROPY_SPH
#ifndef INITIAL_CONDITIONS_CONTAIN_ENTROPY
  NgbTree.setup_entropy_to_invgamma();
#endif
#endif

  for(int i = 0; i < Sp.NumGas; i++)
    {
#ifndef INITIAL_CONDITIONS_CONTAIN_ENTROPY

      if(ThisTask == 0 && i == 0)
        printf("INIT: Converting u -> entropy\n");

#if !defined(PRESSURE_ENTROPY_SPH) && !defined(ISOTHERM_EQS)
      Sp.SphP[i].Entropy = GAMMA_MINUS1 * Sp.SphP[i].Entropy / pow(Sp.SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1);
#endif
      Sp.SphP[i].EntropyPred = Sp.SphP[i].Entropy;

#endif
      /* The predicted entropy values have been already set for all SPH formulation */
      /* so it should be ok computing pressure and csound now */
      Sp.SphP[i].set_thermodynamic_variables();
    }

  if(All.ComovingIntegrationOn)
    {
#ifdef PERIODIC
      if(All.RestartFlag == RST_BEGIN || All.RestartFlag == RST_RESUME || All.RestartFlag == RST_STARTFROMSNAP ||
         All.RestartFlag == RST_CREATEICS)
        {
          /* can't do this check when not all particles are loaded */
          check_omega();
        }
      else
        {
          mpi_printf("Skipping Omega check since not all particles are loaded\n");
        }
#endif
    }

#ifdef STARFORMATION
  /* initialize absolute masses in materials */
  for(int i = 0; i < Sp.NumGas; i++)
    {
      Sp.SphP[i].Metallicity = Sp.P[i].Metallicity;  // set above

      Sp.SphP[i].MassMetallicity = Sp.SphP[i].Metallicity * Sp.P[i].getMass();
    }
#endif

    // tree_based_timesteps_set_initialmaxtistep();

#ifdef DEBUG_MD5
  Logs.log_debug_md5("AFTER-INIT");
#endif

  return;
}

#ifdef PERIODIC
/*! \brief This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.
 *
 * If discrepant, the run is terminated.
 */
void sim::check_omega(void)
{
  double mass = 0;

  for(int i = 0; i < Sp.NumPart; i++)
    mass += Sp.P[i].getMass();

  double masstot;
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, Communicator);

  double omega = masstot * (LONG_Z * LONG_Y * LONG_Z) / (All.BoxSize * All.BoxSize * All.BoxSize) /
                 (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));
  if(fabs(omega - All.Omega0) > 1.0e-2)
    {
      mpi_printf(
          "\n\nI've found something odd!\nThe mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the "
          "parameterfile.\n\nI better stop.\n",
          omega, All.Omega0);
      endrun();
    }
}
#endif

/*! \brief This function is used to find an initial smoothing length for each SPH
 *  particle.
 *
 *  It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void sim::setup_smoothinglengths(void)
{
  Sp.TimeBinsGravity.NActiveParticles = 0;

  for(int i = 0; i < Sp.NumGas; i++)
    Sp.TimeBinsGravity.ActiveParticleList[Sp.TimeBinsGravity.NActiveParticles++] = i;

  sumup_large_ints(1, &Sp.TimeBinsGravity.NActiveParticles, &Sp.TimeBinsGravity.GlobalNActiveParticles, Communicator);

  if(Sp.TimeBinsGravity.GlobalNActiveParticles > 0)
    {
      mpi_printf("INIT: Setup smoothing lengths.\n");

      GravTree.treeallocate(Sp.NumPart, &Sp, &Domain);
      GravTree.treebuild(Sp.TimeBinsGravity.NActiveParticles, Sp.TimeBinsGravity.ActiveParticleList);

      for(int i = 0; i < Sp.NumGas; i++)
        {
          int no = GravTree.Father[i];

          while(10 * All.DesNumNgb * Sp.P[i].getMass() > GravTree.get_nodep(no)->mass)
            {
              int p = GravTree.get_nodep(no)->father;

              if(p < 0)
                break;

              no = p;
            }

          double len;
          if(GravTree.get_nodep(no)->level > 0)
            len = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - GravTree.get_nodep(no)->level)) * Sp.FacIntToCoord;
          else
            len = Sp.RegionLen;

          Sp.SphP[i].Hsml = pow(3.0 / (4 * M_PI) * All.DesNumNgb * Sp.P[i].getMass() / GravTree.get_nodep(no)->mass, 1.0 / 3) * len;
        }

      Sp.TimeBinsHydro.NActiveParticles = 0;
      for(int i = 0; i < Sp.NumGas; i++)
        Sp.TimeBinsGravity.ActiveParticleList[Sp.TimeBinsHydro.NActiveParticles++] = i;

      NgbTree.density(Sp.TimeBinsGravity.ActiveParticleList, Sp.TimeBinsHydro.NActiveParticles);

#ifdef PRESSURE_ENTROPY_SPH
      for(int i = 0; i < Sp.NumGas; i++)
        Sp.SphP[i].PressureSphDensity = Sp.SphP[i].Density;
#endif

      GravTree.treefree();
    }
}

void sim::recreate_unique_ids(void)
{
  mpi_printf("INIT: Setting new unique IDs.\n");

  int *numpart_list = (int *)Mem.mymalloc("numpart_list", NTask * sizeof(int));

  MPI_Allgather(&Sp.NumPart, 1, MPI_INT, numpart_list, 1, MPI_INT, Communicator);

  MyIDType id = 1;

  for(int i = 0; i < ThisTask; i++)
    id += numpart_list[i];

  for(int i = 0; i < Sp.NumPart; i++)
    Sp.P[i].ID.set(id++);

  Mem.myfree(numpart_list);
}

/*! \brief This function checks for unique particle  IDs
 *
 *  The particle IDs are copied to an array and then sorted among all tasks.
 *  This array is then checked for duplicates. In that case the code terminates.
 */
void sim::test_id_uniqueness(void)
{
  mpi_printf("INIT: Testing ID uniqueness...\n");

  double t0 = Logs.second();

  MyIDType *ids       = (MyIDType *)Mem.mymalloc("ids", (Sp.NumPart + 1) * sizeof(MyIDType));
  MyIDType *ids_first = (MyIDType *)Mem.mymalloc("ids_first", NTask * sizeof(MyIDType));
  int *num_list       = (int *)Mem.mymalloc("num_list", NTask * sizeof(int));

  for(int i = 0; i < Sp.NumPart; i++)
    ids[i] = Sp.P[i].ID.get();

  mycxxsort_parallel(ids, ids + Sp.NumPart, Sp.compare_IDs, Communicator);

  for(int i = 1; i < Sp.NumPart; i++)
    {
      if(ids[i] == ids[i - 1])
        Terminate("non-unique ID=%lld found on task=%d (i=%d Sp.NumPart=%d type=%d)\n", (long long)ids[i], ThisTask, i, Sp.NumPart,
                  Sp.P[i].getType());
    }

  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, Communicator);
  MPI_Allgather(&Sp.NumPart, 1, MPI_INT, num_list, 1, MPI_INT, Communicator);

  int next_non_empty_task = ThisTask + 1;

  while(next_non_empty_task < NTask)
    if(num_list[next_non_empty_task] == 0)
      next_non_empty_task++;
    else
      break;

  if(Sp.NumPart > 0 && next_non_empty_task < NTask)
    {
      if(ids[Sp.NumPart - 1] == ids_first[next_non_empty_task])
        Terminate("non-unique ID=%lld found on task=%d\n", (long long)ids[Sp.NumPart - 1], ThisTask);
    }

  Mem.myfree(num_list);
  Mem.myfree(ids_first);
  Mem.myfree(ids);

  double t1 = Logs.second();

  mpi_printf("INIT: success.  took=%g sec\n\n", Logs.timediff(t0, t1));
}
