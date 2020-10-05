
##################################################
#  Enable/Disable compile-time options as needed #
##################################################



#----------------------------------------Parallelization options

IMPOSE_PINNING
#IMPOSE_PINNING_OVERRIDE_MODE
#EXPLICIT_VECTORIZATION         # This uses AVX at selected places through the vectorclass C++ library
#PRESERVE_SHMEM_BINARY_INVARIANCE
#SIMPLE_DOMAIN_AGGREGATION

#--------------------------------------- Basic operation mode of code
PERIODIC
#TWODIMS
#ONEDIMS
#LONG_X_BITS=2
#LONG_Y_BITS=2
#LONG_Z_BITS=1

#NTYPES=6                       # Number of particle types. Note that this may only be changed from the default value of 6 if
                                # HDF5 snapshot files are used.

#GADGET2_HEADER                 # allows reading in of Gadget2/3 snapshots by using this code's header format for snaphot file formats 1 and 2
#SECOND_ORDER_LPT_ICS           # treats second order LPT ICs generated with Adrian Jenkin's code
#LEAN                           # selects a special 'lean' mode of code operation, which is optimized for extreme memory saving

#--------------------------------------- Gravity calculation
SELFGRAVITY                     # switch on for self-gravity
HIERARCHICAL_GRAVITY
#FMM                            # enables FMM tree algorithm
ALLOW_DIRECT_SUMMATION          # allows calculation of small number of particles with direct summation
#EXTERNALGRAVITY                # switch on for external potential
#EVALPOTENTIAL                  # computes gravitational potential
#EXTRAPOTTERM                   # if this is set, the extra multipole term needed in principle for the potential is computed even though it does not enter the force
#EXTRA_HIGH_EWALD_ACCURACY      # if this is activate, a third order Taylor expansion is used for the Ewald corrections 
 
#MULTIPOLE_ORDER=2              # enables in tree and/or FMM for a value of 2 monopoles+dipoles (this is default), 3 gives quadrupoles, 4 octupoles, and 5 hexadecupoles
#RANDOMIZE_DOMAINCENTER

#--------------------------------------- TreePM Options
PMGRID=512
#ASMTH=1.25
#RCUT=6.0
#NTAB=128                       # size of short-range look-up table
#TREEPM_NOTIMESPLIT             # if this is activated, long-range and short-range gravity are time-integrated on a single timestep
#PLACEHIGHRESREGION=2
#HRPMGRID=512                   # High-res PM grid size (optional, default is HRPMGRID=PMGRID)

#FFT_COLUMN_BASED
#PM_ZOOM_OPTIMIZED
#GRAVITY_TALLBOX=2              # this can be used to set-up gravity with two-dimensional boxes
#TREE_NUM_BEFORE_NODESPLIT=10   # Optional number that tells how many particles are allowed in a tree node before it is split

#--------------------------------------- Treatment of gravitational softening
#INDIVIDUAL_GRAVITY_SOFTENING=4+8+16+32
#NSOFTCLASSES=4                   # Number of different softening values. Normally equal to number of particle types, but can be chosen differently if DECOUPLE_TYPE_AND_SOFTTYPE is set.
#ADAPTIVE_HYDRO_SOFTENING



#--------------------------------------- SPH treatmeant and formulation
#PRESSURE_ENTROPY_SPH           # enables the Hopkins 2013 pressure-entropy formulation
#OUTPUT_PRESSURE_SPH_DENSITY    # output also density computed in the pressure-entropy forumulation
#GAMMA=1.4
#ISOTHERM_EQS
#REUSE_HYDRO_ACCELERATIONS_FROM_PREVIOUS_STEP
IMPROVED_VELOCITY_GRADIENTS # use higher-order gradients of the velocity according to Hu et. al (2014)
VISCOSITY_LIMITER_FOR_LARGE_TIMESTEPS #limits the acceleration due to the viscosity  

#--------------------------------------- SPH kernels
CUBIC_SPLINE_KERNEL           # uses the cubic spline kernel
#WENDLAND_C2_KERNEL           # the Wendland C2 kernel from Dehnen & Aly 2012
#WENDLAND_C4_KERNEL           # the Wendland C4 kernel from Dehnen & Aly 2012
#WENDLAND_C6_KERNEL           # the Wendland C6 kernel from Dehnen & Aly 2012
WENDLAND_BIAS_CORRECTION      # reduces self-contribution for Wendland kernels according to Dehnen & Aly 2012



#--------------------------------------- SPH viscosity options
#TIMEDEP_ART_VISC               # Enables time dependend viscosity
#HIGH_ART_VISC_START            # Start with high rather than low viscosity
#NO_SHEAR_VISCOSITY_LIMITER     # Turns of the shear viscosity supression
OUTPUT_VISCOSITY_PARAMETER


#--------------------------------------- Extra physics
#COOLING
#STARFORMATION


#--------------------------------------- Time integration options
#FORCE_EQUAL_TIMESTEPS          # this chooses a variable but global timestep


#---------------------------------------- Single/double precision and data types
POSITIONS_IN_32BIT              # if set, use 32-integers for positions
POSITIONS_IN_64BIT              # if set, use 64-integers for positions
POSITIONS_IN_128BIT             # if set, use 128-integers for positions
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
#OUTPUT_IN_DOUBLEPRECISION      # output files will be written in double precision
#ENLARGE_DYNAMIC_RANGE_IN_TIME  # This extends the dynamic range of the integer timeline from 32 to 64 bit
#IDS_32BIT                      # Selects 32-bit IDs for internal storage (default)
#IDS_48BIT                      # Selects 48-bit IDs for internal storage
#IDS_64BIT                      # Selects 64-bit IDs for internal storage
#USE_SINGLEPRECISION_INTERNALLY  # reduces default double precision for most internal computations to single precision

#---------------------------------------- Output/Input options
#REDUCE_FLUSH
#OUTPUT_VELOCITY_GRADIENT
#OUTPUT_PRESSURE
#OUTPUT_ENTROPY
#OUTPUT_CHANGEOFENTROPY
#OUTPUT_POTENTIAL
#OUTPUT_ACCELERATION
#OUTPUT_TIMESTEP
#OUTPUT_DIVVEL                  # output  velocity divergence
#OUTPUT_CURLVEL                 # output  velocity curl
#OUTPUT_COOLHEAT                # output actual energy loss/gain in cooling/heating routine
#POWERSPEC_ON_OUTPUT
#OUTPUT_NON_SYNCHRONIZED_ALLOWED
#OUTPUT_VELOCITIES_IN_HALF_PRECISION
#OUTPUT_ACCELERATIONS_IN_HALF_PRECISION
#OUTPUT_COORDINATES_AS_INTEGERS
#ALLOW_HDF5_COMPRESSION

#---------------------------------------- On the fly FOF groupfinder
#FOF                            # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2       # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32  # 2^type for the types linked to nearest primaries
##FOF_GROUP_MIN_LEN=32          # default is 32
##FOF_LINKLENGTH=0.16           # Linkinglength for FoF (default=0.2)
#FOF_ALLOW_HUGE_GROUPLENGTH     # if this is set, groups and subhalos may have more than 2 billion particles in length

#---------------------------------------- Subfind
SUBFIND                        # enables substructure finder
#SUBFIND_STORE_LOCAL_DENSITY    # will calculate local densities and velocity dispersions for *all* particles
                                # selected by FOF primary and seconday link types (not only for particles in gropus), 
                                # and store them in snapshots
#SUBFIND_ORPHAN_TREATMENT       # creates special snapshots with formerly most bound particles
#SUBFIND_HBT                    # use previous subhalo catalogue instead of density excursions to define subhalo candidates

#---------------------------------------- Merger tree code
MERGERTREE                      # enables on-the-fly calculation of descendants, and merger tree construction from group catalogues

#---------------------------------------- On-the-fly lightcone creation
#LIGHTCONE
#LIGHTCONE_PARTICLES
#LIGHTCONE_MASSMAPS
#LIGHTCONE_PARTICLES_GROUPS
#LIGHTCONE_OUTPUT_ACCELERATIONS
#LIGHTCONE_IMAGE_COMP_HSML_VELDISP

#REARRANGE_OPTION

#--------------------------------------- IC generation
#NGENIC=256                     # generate cosmological ICs, set NGENIC to the FFT grid size used for IC generation
#CREATE_GRID                    # start with a regular cartesian DM particle grid, instead of reading a glass file (for NGENIC)
#GENERATE_GAS_IN_ICS            # add SPH particles to dark matter only ICs
#SPLIT_PARTICLE_TYPE=4+8        # particle types to be split if GENERATE_GAS_IN_ICS is activated
#NGENIC_FIX_MODE_AMPLITUDES     # when activated, this leaves the mode amplitudes at sqrt(P(k)), instead of sampling from a Rayleigh distribution
#NGENIC_MIRROR_PHASES           # if this is activated, all phases are turned by 180 degrees
#NGENIC_2LPT                    # applies 2LPT instead of just Zeldovich
#NGENIC_TEST                    # can be used to create ICs, measure the power spectrum, and then stop

#---------------------------------------- MPI related settings
#ISEND_IRECV_IN_DOMAIN
#USE_MPIALLTOALLV_IN_DOMAINDECOMP
#MPI_HYPERCUBE_ALLGATHERV       # some MPI-libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPI_MESSAGE_SIZELIMIT_IN_MB=200
#NUMPART_PER_TASK_LARGE         # Set this if the number of particle per task is very large (so that more than 2GB data is comprised just by the particle data)

#MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY=64
#NUMBER_OF_MPI_LISTENERS_PER_NODE=1

#---------------------------------------- Testing and Debugging options
DEBUG                           # enables core-dumps
#DEBUG_ENABLE_FPU_EXCEPTIONS    # tries to enable FPU exceptions
#DEBUG_SYMTENSORS               # carries out a few unit tests related to the tensor algebra routines of the code
HOST_MEMORY_REPORTING           # reports after start-up the available system memory by analyzing /proc/meminfo
#FORCETEST=0.001                # calculates for given fraction of particles direct summation forces to check accuracy of tree force
#FORCETEST_FIXEDPARTICLESET     # if this set, always the same particle receive a force accuracy check during a run
#FORCETEST_TESTFORCELAW=1       # Special option for measuring the effective force law, can be set to 1 or 2 for TreePM+HighRes
#VTUNE_INSTRUMENT
#DEBUG_MD5                      # make MD5 check sum routines available for debugging
#SQUASH_TEST
#DOMAIN_SPECIAL_CHECK
#RECREATE_UNIQUE_IDS            # Before carrying out the ID uniqueness test, this sets new IDs. Can be used if IC files contain broken IDs.
#NO_STOP_BELOW_MINTIMESTEP      # do not stop when the code wants to go below minimum timestep
#TILING=2                       # duplicated an IC in each dimension
#DO_NOT_PRODUCE_BIG_OUTPUT      # for scaling tests, one may disable the writing of restart files and snapshot dumps
#MEASURE_TOTAL_MOMENTUM
#ENABLE_HEALTHTEST
#EWALD_TEST                     # a development test for the Ewald tables
#STOP_AFTER_STEP=10             # ends a simulation after the specified timestep (to simplify scalability tests)

#TREE_NO_SAFETY_BOX             # when set, this disables the geometric 'near node' protection

#---------------------------------------- Postprocessing options

#LGALAXIES                      # Semi-analytic L-Galaxies code

