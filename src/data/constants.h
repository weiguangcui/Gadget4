/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file constants.h
 *
 *  \brief declares global constants and checks plausibility of configuration
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define GADGET_VERSION "4.0" /* code version string */

#define FILEFORMAT_LEGACY1 1
#define FILEFORMAT_LEGACY2 2
#define FILEFORMAT_HDF5 3

#define MODE_LOCAL_PARTICLES 0
#define MODE_IMPORTED_PARTICLES 1
#define MODE_DEFAULT 2
#define MODE_LOCAL_NO_EXPORT 3

#define FIRST_HALF_STEP 0
#define SECOND_HALF_STEP 1

#define FLAG_OUTSIDE 0
#define FLAG_INSIDE 1
#define FLAG_BOUNDARYOVERLAP 2

#define LOW_MESH 0  /* low-res  mesh selector */
#define HIGH_MESH 1 /* high-res mesh selector */

#define MAX_THREADS 128

#ifndef DIRECT_SUMMATION_THRESHOLD
#define DIRECT_SUMMATION_THRESHOLD 500
#endif

#define NUMBER_OF_MEASUREMENTS_TO_RECORD 6

#define MAX_FIRST_ELEMENTS_CONSIDERED \
  5 /* This sets the number of lowest loaded tasks to be considered for assignment of next domain patch */

#define COMMBUFFERSIZE (32 * 1024LL * 1024LL)

#ifndef MPI_MESSAGE_SIZELIMIT_IN_MB
#define MPI_MESSAGE_SIZELIMIT_IN_MB 200
#endif

#define MPI_MESSAGE_SIZELIMIT_IN_BYTES ((MPI_MESSAGE_SIZELIMIT_IN_MB)*1024LL * 1024LL)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TO_MBYTE_FAC (1.0 / (1024.0 * 1024.0))

#ifndef LIGHTCONE_ALLOC_FAC
#define LIGHTCONE_ALLOC_FAC 0.1
#endif

#ifndef LIGHTCONE_MASSMAP_ALLOC_FAC
#define LIGHTCONE_MASSMAP_ALLOC_FAC 1.0
#endif

#ifndef LIGHTCONE_MAX_FILLFACTOR
#define LIGHTCONE_MAX_FILLFACTOR 0.9
#endif

#ifndef ALLOC_TOLERANCE
#define ALLOC_TOLERANCE 0.2
#endif

#define ALLOC_STARBH_ROOM 0.02

#define MAX_FLOAT_NUMBER 1e37
#define MIN_FLOAT_NUMBER 1e-37
#define MAX_DOUBLE_NUMBER 1e306
#define MIN_DOUBLE_NUMBER 1e-306
#define SMALLNUM 1e-60

#ifdef DOUBLEPRECISION
#if(DOUBLEPRECISION == 2)
#define MAX_REAL_NUMBER MAX_FLOAT_NUMBER
#define MIN_REAL_NUMBER MIN_FLOAT_NUMBER
#else
#define MAX_REAL_NUMBER MAX_DOUBLE_NUMBER
#define MIN_REAL_NUMBER MIN_DOUBLE_NUMBER
#endif
#else
#define MAX_REAL_NUMBER MAX_FLOAT_NUMBER
#define MIN_REAL_NUMBER MIN_FLOAT_NUMBER
#endif

#if !defined(GAMMA) && !defined(ISOTHERM_EQS)
#define GAMMA (5.0 / 3) /**< adiabatic index of simulated gas */
#endif

#ifdef ISOTHERM_EQS
#if defined(GAMMA)
#error "ISOTHERM_EQS overwrites your definition of GAMMA"
#endif
#undef GAMMA
#define GAMMA 1.0
#endif

#define GAMMA_MINUS1 (GAMMA - 1)

#define HYDROGEN_MASSFRAC 0.76 /**< mass fraction of hydrogen, relevant only for radiative cooling */

#define METAL_YIELD 0.02 /**< effective metal yield for star formation */

/* ... often used physical constants (cgs units; NIST 2010) */

#define GRAVITY 6.6738e-8
#define SOLAR_MASS 1.989e33
#define BOLTZMANN 1.38065e-16
#define CLIGHT 2.99792458e10

#define PARSEC 3.085678e18
#define PROTONMASS 1.67262178e-24
#define HUBBLE 3.2407789e-18 /* in h/sec */

#define SEC_PER_MEGAYEAR 3.15576e13
#define SEC_PER_YEAR 3.15576e7

#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif

#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
#endif

#ifndef FOF_LINKLENGTH
#define FOF_LINKLENGTH 0.2
#endif

#ifndef FOF_GROUP_MIN_LEN
#define FOF_GROUP_MIN_LEN 32
#endif

#if defined(PMGRID) && !defined(HRPMGRID)
#define HRPMGRID PMGRID
#endif

#if !defined(RANDOMIZE_DOMAINCENTER_TYPES) && defined(PLACEHIGHRESREGION)
#define RANDOMIZE_DOMAINCENTER_TYPES PLACEHIGHRESREGION
#endif

#if defined(SUBFIND) && !defined(SELFGRAVITY)
#error "Running SUBFIND without SELFGRAVITY enabled does not make sense."
#endif

#if defined(SUBFIND) && !defined(FOF)
#error "Running SUBFIND without FOF is not possible."
#endif

#if defined(TILING) && !defined(RECREATE_UNIQUE_IDS)
#error "Running with TILING requires RECREATE_UNIQUE_IDS"
#endif

#if defined(MULTIPOLE_ORDER)
#if(MULTIPOLE_ORDER < 1) || (MULTIPOLE_ORDER > 5)
#error "MULTIPOLE_ORDER must be either 1, 2, 3, 4, or 5"
#endif
#else
#define MULTIPOLE_ORDER 1
#endif

#if defined(FORCETEST) && !defined(EVALPOTENTIAL)
#error "Running with FORCETEST requires EVALPOTENTIAL."
#endif

#if defined(DEBUG_ENABLE_FPU_EXCEPTIONS) && !defined(__linux__)
#warning "DEBUG_ENABLE_FPU_EXCEPTIONS only works under Linux."
#undef DEBUG_ENABLE_FPU_EXCEPTIONS
#endif

#if defined(HOST_MEMORY_REPORTING) && !defined(__linux__)
#warning "HOST_MEMORY_REPORTING only works under Linux."
#undef HOST_MEMORY_REPORTING
#endif

#if !defined(HOST_MEMORY_REPORTING) && defined(__linux__)
#define HOST_MEMORY_REPORTING  // let's switch it always on under Linux
#endif

#if defined(STARFORMATION) && !defined(COOLING)
#error "STARFORMATION requires COOLING"
#endif

#if defined(FORCE_EQUAL_TIMESTEPS) && defined(HIERARCHICAL_GRAVITY)
#error "FORCE_EQUAL_TIMESTEPS cannot be used together with HIERARCHICAL_GRAVITY"
#endif

#if defined(EXTRAPOTTERM) && !defined(EVALPOTENTIAL)
#error "EXTRAPOTTERM makes only sense for EVALPOTENTIAL"
#endif

#if defined(MERGERTREE) && !defined(SUBFIND)
#error "MERGERTREE requires SUBFIND."
#endif

#if defined(POWERSPEC_ON_OUTPUT) && !(defined(PERIODIC) && defined(PMGRID))
#error "The option POWERSPEC_ON_OUTPUT requires PMGRID and PERIODIC."
#endif

#if defined(CREATE_GRID) && !defined(NGENIC)
#error "CREATE_GRID only makes sense with NGENIC"
#endif

#if defined(OUTPUT_COORDINATES_AS_INTEGERS) && !defined(PERIODIC)
#error "The OUTPUT_COORDINATES_AS_INTEGERS option is only allowed when PERIODIC is on"
#endif

#if defined(ALLOW_DIRECT_SUMMATION) && !defined(HIERARCHICAL_GRAVITY)
#error "The option ALLOW_DIRECT_SUMMATION is only availble when HIERARCHICAL_GRAVITY is used"
#endif

#if defined(PMGRID) && !defined(PERIODIC) && !defined(TREEPM_NOTIMESPLIT)
#error "If PMGRID is used without PERIODIC, TREEPM_NOTIMESPLIT needs to be activated"
#endif

#if defined(PMGRID) && defined(HIERARCHICAL_GRAVITY) && !defined(TREEPM_NOTIMESPLIT)
#error "If PMGRID is used together with HIERARCHICAL_GRAVITY, you also need to use TREEPM_NOTIMESPLIT"
#endif

#if defined(PLACEHIGHRESREGION) && !defined(RANDOMIZE_DOMAINCENTER)
#error "PLACEHIGHRESREGION requires RANDOMIZE_DOMAINCENTER."
#endif

#if defined(PLACEHIGHRESREGION) && !defined(PMGRID)
#error "PLACEHIGHRESREGION requires PMGRID."
#endif

#if defined(TREEPM_NOTIMESPLIT) && !defined(PMGRID)
#error "The option TREEPM_NOTIMESPLIT requires PMGRID."
#endif

#if defined(HRPMGRID) && !defined(PMGRID)
#error "It doesn't make sense to set HRPMGRID without having PMGRID."
#endif

#if defined(POSITIONS_IN_64BIT) && defined(POSITIONS_IN_128BIT)
#error "The options POSITIONS_IN_64BIT and POSITIONS_IN_128BIT should not be activated together."
#endif

#if !defined(EVALPOTENTIAL) && defined(FORCETEST)
#error "When you enable FORCETEST you should also switch on EVALPOTENTIAL"
#endif

#if(defined(LONG_X_BITS) || defined(LONG_Y_BITS) || defined(LONG_Z_BITS)) && !defined(PERIODIC)
#error "LONG_X/Y/Z_BITS requires the PERIODIC option"
#endif

#if defined(LIGHTCONE_PARTICLES) && !defined(LIGHTCONE)
#error "The option LIGHTCONE_PARTICLES requires LIGHTCONE"
#endif

#if defined(LIGHTCONE_MASSMAPS) && !defined(LIGHTCONE)
#error "The option LIGHTCONE_MASSMAPS requires LIGHTCONE"
#endif

#if defined(LIGHTCONE) && !defined(LIGHTCONE_PARTICLES) && !defined(LIGHTCONE_MASSMAPS)
#error "The option LIGHTCONE requires selection of at least one of LIGHTCONE_PARTICLES or LIGHTCONE_MASSMAPS"
#endif

#if defined(SUBFIND_HBT) && !defined(MERGERTREE)
#error "The option SUBFIND_HBT requires MERGERTREE"
#endif

#if defined(LIGHTCONE_PARTICLES) && defined(LIGHTCONE_OUTPUT_ACCELERATIONS) && defined(PMGRID) && \
    (!defined(TREEPM_NOTIMESPLIT) || defined(HIERARCHICAL_GRAVITY))
#error "LIGHTCONE_OUTPUT_ACCELERATIONS only works with PMGRID if TREEPM_NOTIMESPLIT is used and HIERARCHICAL_GRAVITY is not used"
#endif

#if defined(EXTERNALGRAVITY_STATICHQ) && !defined(EXTERNALGRAVITY)
#error "EXTERNALGRAVITY_STATICHQ only works when EXTERNALGRAVITY is activated"
#endif

#if defined(LIGHTCONE_MULTIPLE_ORIGINS) && defined(LIGHTCONE_PARTICLES_GROUPS)
#error "Presently, the option LIGHTCONE_MULTIPLE_ORIGINS cannot be used yet together with LIGHTCONE_PARTICLES_GROUPS"
#endif

#if defined(LIGHTCONE_MULTIPLE_ORIGINS) && defined(LIGHTCONE_MASSMAPS)
#error "Presently, the option LIGHTCONE_MULTIPLE_ORIGINS cannot be used yet together with LIGHTCONE_MASSMAPS"
#endif

#ifndef ASMTH
/** ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25
#endif
#ifndef RCUT
/** RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT 7.0
#endif

#ifndef MAXLEN_OUTPUTLIST
#define MAXLEN_OUTPUTLIST 1100 /**< maxmimum number of entries in output list */
#endif

#define MAXLEN_PATH 512        /**< maximum length of various filenames (full path) */
#define MAXLEN_PATH_EXTRA 2048 /**< maximum length of various filenames, plus extra space */

#define BASENUMBER 100

#define MAXITER 10000

#ifndef NTYPES
#define NTYPES 6
#endif

#ifndef NSOFTCLASSES
#define NSOFTCLASSES NTYPES
#endif

#ifdef ADAPTIVE_HYDRO_SOFTENING
#ifndef NSOFTCLASSES_HYDRO
#define NSOFTCLASSES_HYDRO 64
#endif
#else
#undef NSOFTCLASSES_HYDRO
#define NSOFTCLASSES_HYDRO 0
#endif

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
typedef long long integertime;
#define TIMEBINS 60
#define TIMEBASE                                                                                           \
  (((long long)1) << TIMEBINS) /* The simulated timespan is mapped onto the integer interval [0,TIMESPAN], \
                                *  where TIMESPAN needs to be a power of 2. */
#else
typedef int integertime;
#define TIMEBINS 29
#define TIMEBASE (1 << TIMEBINS)
#endif

#if(NSOFTCLASSES + NSOFTCLASSES_HYDRO) >= 128
#error "(NSOFTCLASSES + NSOFTCLASSES_HYDRO) must be smaller than 128"
#endif

#if NSOFTCLASSES < 1
#error "NSOFTCLASSES must be at least 1"
#endif

#ifdef GADGET2_HEADER
#if NTYPES > 6
#error "NTYPES may not be larger than 6 if GADGET2_HEADER is set"
#endif
#endif

#ifndef STAR_TYPE
#define STAR_TYPE 4
#endif

#if defined(STARFORMATION) && (STAR_TYPE >= NTYPES)
#error "STAR_TYPE must be an available type according to the set NTYPES"
#endif

#ifdef ONEDIMS
#define NUMDIMS 1
#define KERNEL_COEFF_1 (4.0 / 3)
#define KERNEL_COEFF_2 (8.0)
#define KERNEL_COEFF_3 (24.0)
#define KERNEL_COEFF_4 (16.0)
#define KERNEL_COEFF_5 (8.0 / 3)
#define KERNEL_COEFF_6 (-8.0)
#define NORM_COEFF 2.0
#else
#ifndef TWODIMS
#define NUMDIMS 3                     /**< For 3D-normalized kernel */
#define KERNEL_COEFF_1 2.546479089470 /**< Coefficients for SPH spline kernel and its derivative */
#define KERNEL_COEFF_2 15.278874536822
#define KERNEL_COEFF_3 45.836623610466
#define KERNEL_COEFF_4 30.557749073644
#define KERNEL_COEFF_5 5.092958178941
#define KERNEL_COEFF_6 (-15.278874536822)
#define NORM_COEFF 4.188790204786 /**< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#define NUMDIMS 2                                 /**< For 2D-normalized kernel */
#define KERNEL_COEFF_1 (5.0 / 7 * 2.546479089470) /**< Coefficients for SPH spline kernel and its derivative */
#define KERNEL_COEFF_2 (5.0 / 7 * 15.278874536822)
#define KERNEL_COEFF_3 (5.0 / 7 * 45.836623610466)
#define KERNEL_COEFF_4 (5.0 / 7 * 30.557749073644)
#define KERNEL_COEFF_5 (5.0 / 7 * 5.092958178941)
#define KERNEL_COEFF_6 (5.0 / 7 * (-15.278874536822))
#define NORM_COEFF M_PI /**< Coefficient for kernel normalization. */
#endif
#endif /* ONEDIMS */

#define SOFTFAC1 (32.0 / 3) /**< Coefficients for gravitational softening */
#define SOFTFAC2 32.0
#define SOFTFAC3 (-38.4)
#define SOFTFAC4 (-2.8)
#define SOFTFAC5 (16.0 / 3)
#define SOFTFAC6 6.4
#define SOFTFAC7 (-9.6)
#define SOFTFAC8 (64.0 / 3)
#define SOFTFAC9 (-48.0)
#define SOFTFAC10 38.4
#define SOFTFAC11 (-32.0 / 3)
#define SOFTFAC12 (-1.0 / 15)
#define SOFTFAC13 (-3.2)
#define SOFTFAC14 (1.0 / 15)
#define SOFTFAC15 (-16.0)
#define SOFTFAC16 9.6
#define SOFTFAC17 (-64.0 / 30)
#define SOFTFAC18 128.0
#define SOFTFAC19 (-115.2)
#define SOFTFAC20 (64.0 / 3)
#define SOFTFAC21 (-96.0)
#define SOFTFAC22 115.2
#define SOFTFAC23 (-128.0 / 3)
#define SOFTFAC24 (4.0 / 30)

#define SOFTFAC30 (32.0 / 3)
#define SOFTFAC31 (-576.0 / 5)
#define SOFTFAC32 (128.0)
#define SOFTFAC33 (-1152.0 / 5)
#define SOFTFAC34 (384.0)
#define SOFTFAC35 (2.0 * 384.0)

#define SOFTFAC40 (64.0 / 3)
#define SOFTFAC41 (2.0 / 15)
#define SOFTFAC42 (-96.0)
#define SOFTFAC43 (576.0 / 5)
#define SOFTFAC44 (-128.0 / 3)
#define SOFTFAC45 (-96.0)
#define SOFTFAC46 (-2.0 / 5)
#define SOFTFAC47 (1152.0 / 5)
#define SOFTFAC48 (-128.0)
#define SOFTFAC49 (8.0 / 5)
#define SOFTFAC50 (-256.0)
#define SOFTFAC51 (-8.0)

#define SQRT_PI 1.772453850906           /* sqrt(M_PI) */
#define FACT1 0.366025403785             /* FACT1 = 0.5 * (sqrt(3)-1) */
#define FACTSQRT3HALF 0.866025403785     /* sqrt(3)/2 */
#define FACTSQRT3 (2.0 * 0.866025403785) /* sqrt(3) */
#endif
