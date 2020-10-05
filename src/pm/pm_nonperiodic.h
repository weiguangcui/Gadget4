/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  pm_nonperiodic.h
 *
 *  \brief declaration of a class used for non-periodic long-range PM force calculation
 */

#ifndef PM_NONPERIODIC_H
#define PM_NONPERIODIC_H

#include "gadgetconfig.h"

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

class pm_nonperiodic : public pm_mpi_fft
{
 public:
  pm_nonperiodic(MPI_Comm comm) : setcomm(comm), pm_mpi_fft(comm) {}

#if defined(PMGRID) && (!defined(PERIODIC) || defined(PLACEHIGHRESREGION))

 private:
#if defined(LONG_X_BITS) || defined(LONG_Y_BITS) || defined(LONG_Z_BITS)
#error "LONG_X/Y/Z_BITS not supported for the non-periodic FFT gravity code"
#endif

#if defined(GRAVITY_TALLBOX)
#error "GRAVITY_TALLBOX not supported for the non-periodic FFT gravity code"
#endif

#if(HRPMGRID > 1024)
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

  /* short-cut macros for accessing different 3D arrays */

  int NSource;

  fft_plan myplan; /*!< In this structure, various bookkeeping variables for the distributed FFTs are stored */

  /*! \var maxfftsize
   *  \brief maximum size of the local fft grid among all tasks
   */
  size_t maxfftsize;

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

  fft_real *kernel[2];
  fft_complex *fft_of_kernel[2];

 public:
  simparticles *Sp;

  void pm_init_nonperiodic(simparticles *Sp_ptr);
  void pm_init_regionsize(void);
  int pmforce_nonperiodic(int grnr);
  void pm_setup_nonperiodic_kernel(void);
  void pmforce_nonperiodic_zoom_optimized_prepare_density(int grnr);

 private:
#ifdef PM_ZOOM_OPTIMIZED

  void pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(int grnr, int dim);

  /*! \brief This structure links the particles to the mesh cells, to which they contribute their mass
   *
   * Each particle will have eight items of this structure in the #part array.
   * For each of the eight mesh cells the CIC assignment will contribute,
   * one item of this struct exists.
   */
  struct part_slab_data
  {
    large_array_offset globalindex; /*!< index in the global density mesh */
    large_numpart_type partindex;   /*!< contains the local particle index shifted by 2^3, the first three bits encode to which part of
                                       the CIC assignment this item belongs to */
    large_array_offset localindex;  /*!< index to a local copy of the corresponding mesh cell of the global density array (used during
                                       local mass and force assignment) */
  };
  part_slab_data *part; /*!< array of part_slab_data linking the local particles to their mesh cells */

  size_t *localfield_sendcount, *localfield_first, *localfield_offset, *localfield_recvcount;
  large_array_offset *localfield_globalindex, *import_globalindex;
  fft_real *localfield_data, *import_data;
  large_numpart_type num_on_grid;

  /* realize the comparison function as a functor, so that it can have an internal state (here the data array for which we sort indices
   */
  struct pm_nonperiodic_sortindex_comparator
  {
   private:
    const part_slab_data *data;

   public:
    pm_nonperiodic_sortindex_comparator(const part_slab_data *data_) : data(data_) {}

    bool operator()(const large_numpart_type &a, const large_numpart_type &b) const
    {
      return data[a].globalindex < data[b].globalindex;
    }
  };

#else

  /*
   *  Here come the routines for a different communication algorithm that is better suited for a homogenuously loaded boxes.
   */
  struct partbuf
  {
    MyIntPosType IntPos[3];
    MyFloat Mass;
  };
  partbuf *partin, *partout;

  size_t nimport, nexport;

  size_t *Sndpm_count, *Sndpm_offset;
  size_t *Rcvpm_count, *Rcvpm_offset;

  void pmforce_nonperiodic_uniform_optimized_prepare_density(int grnr);
  void pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(int grnr, int dim);

#endif

#endif
};

#endif
