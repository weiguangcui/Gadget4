/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pm_periodic.h
 *
 *  \brief declaration of a class used for periodic PM-force calculations
 */

#ifndef PM_PERIODIC_H
#define PM_PERIODIC_H

#include "gadgetconfig.h"

#include <gsl/gsl_integration.h>

#if defined(PMGRID) || defined(NGENIC)

#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <math.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../data/simparticles.h"
#include "../domain/domain.h"
#include "../logs/timer.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm_mpi_fft.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

class pm_periodic : public pm_mpi_fft
{
 public:
  pm_periodic(MPI_Comm comm) : setcomm(comm), pm_mpi_fft(comm) {}

#if defined(PMGRID) && defined(PERIODIC)

 private:
#ifdef LONG_X_BITS
#if PMGRID != ((PMGRID / LONG_X) * LONG_X)
#error "PMGRID must be a multiple of the stretch factor in the x-direction"
#endif
#endif

#ifdef LONG_Y_BITS
#if PMGRID != ((PMGRID / LONG_Y) * LONG_Y)
#error "PMGRID must be a multiple of the stretch factor in the y-direction"
#endif
#endif

#ifdef LONG_Z_BITS
#if PMGRID != ((PMGRID / LONG_Z) * LONG_Z)
#error "PMGRID must be a multiple of the stretch factor in the x-direction"
#endif
#endif

#define GRIDX ((PMGRID / LONG_X) * DBX + DBX_EXTRA)
#define GRIDY ((PMGRID / LONG_Y) * DBY + DBY_EXTRA)
#define GRIDZ ((PMGRID / LONG_Z) * DBZ + DBZ_EXTRA)

#define INTCELL ((~((MyIntPosType)0)) / PMGRID + 1)

#if(GRIDX > 1024) || (GRIDY > 1024) || (GRIDZ > 1024)
  typedef long long large_array_offset; /* use a larger data type in this case so that we can always address all cells of the 3D grid
                                           with a single index */
#else
  typedef int large_array_offset;
#endif

#ifdef NUMPART_PER_TASK_LARGE
  typedef long long large_numpart_type; /* if there is a risk that the local particle number times 8 overflows a 32-bit integer, this
                                           data type should be used */
#else
  typedef int large_numpart_type;
#endif

  /* variables for power spectrum estimation */
#ifndef BINS_PS
#define BINS_PS 4000 /* number of bins for power spectrum computation */
#endif
#ifndef POWERSPEC_FOLDFAC
#define POWERSPEC_FOLDFAC 16 /* folding factor to obtain an estimate of the power spectrum on very small scales */
#endif

  char power_spec_fname[MAXLEN_PATH_EXTRA];

  int NSource;

  void pmforce_measure_powerspec(int flag, int *typeflag);
  void pmforce_do_powerspec(int *typeflag);
#if defined(GRAVITY_TALLBOX)
  void pmforce_setup_tallbox_kernel(void);
  double pmperiodic_tallbox_long_range_potential(double x, double y, double z);
#endif

  fft_plan myplan; /*!< In this structure, various bookkeeping variables for the distributed FFTs are stored */

  /*! \var maxfftsize
   *  \brief maximum size of the local fft grid among all tasks
   */
  long long maxfftsize;

  /*! \var rhogrid
   *  \brief This array hold the local part of the density field and
   *  after the FFTs the local part of the potential
   *
   *  \var forcegrid
   *  \brief This array will contain the force field
   *
   *  \var workspace
   *  \brief Workspace array used during the FFTs
   */
  fft_real *rhogrid, *forcegrid, *workspace;

  /*! \brief Array containing the FFT of #rhogrid
   *
   *  This pointer points to the same array as #rhogrid,
   *  because in-place FFTs are used.
   */
  fft_complex *fft_of_rhogrid;

#if defined(GRAVITY_TALLBOX)
  fft_real *kernel; /*!< If the tallbox option is used, the code will construct and store the k-space Greens function by FFTing it from
                       real space */
  fft_complex *fft_of_kernel;
#endif

#ifdef PM_ZOOM_OPTIMIZED

  /*! \brief This structure links the particles to the mesh cells, to which they contribute their mass
   *
   * Each particle will have eight items of this structure in the #part array.
   * For each of the eight mesh cells the CIC assignment will contribute,
   * one item of this struct exists.
   */

 public:
  struct part_slab_data
  {
    large_array_offset globalindex; /*!< index in the global density mesh */
    large_numpart_type partindex;   /*!< contains the local particle index shifted by 2^3, the first three bits encode to which part of
                                       the CIC assignment this item belongs to */
    large_array_offset localindex;  /*!< index to a local copy of the corresponding mesh cell of the global density array (used during
                                       local mass and force assignment) */
  };
  part_slab_data *part; /*!< array of part_slab_data linking the local particles to their mesh cells */

  /* realize the comparison function as a functor, so that it can have an internal state (here the data array for which we sort indices
   */
  struct pm_periodic_sortindex_comparator
  {
   private:
    const part_slab_data *data;

   public:
    pm_periodic_sortindex_comparator(const part_slab_data *data_) : data(data_) {}

    bool operator()(const large_numpart_type &a, const large_numpart_type &b) const
    {
      return data[a].globalindex < data[b].globalindex;
    }
  };

 private:
  size_t *localfield_sendcount, *localfield_first, *localfield_offset, *localfield_recvcount;
  large_array_offset *localfield_globalindex, *import_globalindex;
  fft_real *localfield_data, *import_data;

  void pmforce_zoom_optimized_prepare_density(int mode, int *typelist);
  void pmforce_zoom_optimized_readout_forces_or_potential(fft_real *grid, int dim);

#else

  struct partbuf
  {
#ifndef LEAN
    MyFloat Mass;
#endif
    MyIntPosType IntPos[3];
  };
  partbuf *partin, *partout;

  size_t nimport, nexport;

  size_t *Sndpm_count, *Sndpm_offset;
  size_t *Rcvpm_count, *Rcvpm_offset;

  void pmforce_uniform_optimized_prepare_density(int mode, int *typelist);

  void pmforce_uniform_optimized_readout_forces_or_potential_xy(fft_real *grid, int dim);

  void pmforce_uniform_optimized_readout_forces_or_potential_xz(fft_real *grid, int dim);
  void pmforce_uniform_optimized_readout_forces_or_potential_zy(fft_real *grid, int dim);
#endif

 public:
  simparticles *Sp;

  void pm_init_periodic(simparticles *Sp_ptr);
  void pmforce_periodic(int mode, int *typelist);

  void calculate_power_spectra(int num);

  static double growthfactor_integrand(double a, void *param)
  {
    return pow(a / (All.Omega0 + (1 - All.Omega0 - All.OmegaLambda) * a + All.OmegaLambda * a * a * a), 1.5);
  }

  double linear_growth_factor(double astart, double aend) { return linear_growth(aend) / linear_growth(astart); }

  double linear_growth(double a)
  {
    double hubble_a = sqrt(All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a) + All.OmegaLambda);

    const int worksize = 100000;

    double result, abserr;
    gsl_function F;

    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(worksize);
    F.function                           = &growthfactor_integrand;

    gsl_integration_qag(&F, 0, a, 0, 1.0e-8, worksize, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return hubble_a * result;
  }
#endif
};

#endif

#endif
