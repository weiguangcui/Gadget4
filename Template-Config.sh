
#########################################
#  Pick compile-time options as needed  #
#  from this template.                  #
#########################################


#--------------------------------------- Basic operation mode of code

#PERIODIC                                     # enables periodic boundary condistions
#NTYPES=6                                     # number of particle types 
#RANDOMIZE_DOMAINCENTER                       # shifts the particle distribution randomly each step to reduce correlations of force errors in time
#RANDOMIZE_DOMAINCENTER_TYPES                 # can be used to set a zoom region via a particle type which is then never placed across large node boundaries
#LEAN                                         # selects a special 'lean' mode of code operation, which is optimized for aggressive memory saving
#LONG_X_BITS=2                                # can be used to reduce periodic box-dimension in x-direction relative to nominal box size
#LONG_Y_BITS=2                                # can be used to reduce periodic box-dimension in y-direction relative to nominal box size
#LONG_Z_BITS=1                                # can be used to reduce periodic box-dimension in z-direction relative to nominal box size
#GRAVITY_TALLBOX=2                            # this can be used to set-up gravity with two periodicity only in two of three dimensions
#TWODIMS                                      # restricts SPH formulation to two dimensions (for tests)
#ONEDIMS                                      # restricts SPH formulation to one dimension (for tests)
#GADGET2_HEADER                               # allows reading of snapshots with GADGET-2/3 header format 
#SECOND_ORDER_LPT_ICS                         # treats second order LPT ICs generated with Adrian Jenkin's code


#--------------------------------------- Gravity calculation

SELFGRAVITY                                   # switch to enable self-gravity of particles (typically always on)
#FMM                                          # enables Fast Multipole Method instead of one-sided tree algorithm
#MULTIPOLE_ORDER=2                            # sets the multipole order of Tree or FMM computations
#EVALPOTENTIAL                                # computes gravitational potential besides force
#EXTRAPOTTERM                                 # this computes an extra multipole term for the potential which is not needed for the forces
#EXTRA_HIGH_EWALD_ACCURACY                    # this uses third-order instead of second-order Taylor expansion to interpolate Ewald corrections from table 
#ALLOW_DIRECT_SUMMATION                       # allows calculation of direct summation gravity force if only a tiny number of particles as active 
#EXTERNALGRAVITY                              # switches on inclusion of external gravitational potential
#EXTERNALGRAVITY_STATICHQ                     # example for a simple external potential due to a Hernquist halo

#--------------------------------------- TreePM Options

#PMGRID=512                                   # basic mesh size for TreePM calculations
#ASMTH=1.25                                   # sets the smoothing of the PM force and thus the split-scale between Tree/FMM force and PM force
#RCUT=6.0                                     # cut-off radius beyond which Tree or FMM evaluations are stopped in TreePM / FMM-PM
#NTAB=128                                     # size of short-range look-up table
#PLACEHIGHRESREGION=2                         # bitmask selecting the particle types that are used to define the location of a secondary PM mesh
#HRPMGRID=512                                 # dimension of high-res PM grid (optional, default is HRPMGRID=PMGRID)
#FFT_COLUMN_BASED                             # uses a column-based FFT algorithm instead of the default slab-based one
#PM_ZOOM_OPTIMIZED                            # selects a communication strategy in the PM code that is better balanced for zoom simulations
#TREE_NUM_BEFORE_NODESPLIT=4                  # number of particles are are at most allowed in a tree node before it is split (can be 1)


#--------------------------------------- Time integration options

#TREEPM_NOTIMESPLIT                           # if this is activated, long-range and short-range gravity are time-integrated on a common timestep
#HIERARCHICAL_GRAVITY                         # enables hierarchical time integration of the gravity 
#FORCE_EQUAL_TIMESTEPS                        # this chooses a global timestep for all particles


#--------------------------------------- Treatment of gravitational softening

#INDIVIDUAL_GRAVITY_SOFTENING=4+8+16+32       # bitmasks which selects the particle type(s) which pick their softening class based on particle mass
#NSOFTCLASSES=4                               # number of different softening classes
#ADAPTIVE_HYDRO_SOFTENING                     # makes SPH gas particles pick an adaptive gravitational softening proportional to their SPH smoothing lengths


#--------------------------------------- SPH treatmeant and formulation

#REUSE_HYDRO_ACCELERATIONS_FROM_PREVIOUS_STEP # does not recompute the pressure forces after application of source functions 
#VISCOSITY_LIMITER_FOR_LARGE_TIMESTEPS        # limits the acceleration due to the viscosity  
#PRESSURE_ENTROPY_SPH                         # enables the Hopkins (2013) pressure-entropy formulation, other density-entropy is used
#GAMMA=1.4                                    # sets the adiabatic index
#ISOTHERM_EQS                                 # selects an isothermal equation of state
#IMPROVED_VELOCITY_GRADIENTS                  # use higher-order gradients of the velocities according to Hu et. al (2014)


#--------------------------------------- SPH kernels

#CUBIC_SPLINE_KERNEL                          # uses the cubic spline kernel (default)
#WENDLAND_C2_KERNEL                           # use the Wendland C2 kernel from Dehnen & Aly 2012
#WENDLAND_C4_KERNEL                           # use the Wendland C4 kernel from Dehnen & Aly 2012
#WENDLAND_C6_KERNEL                           # use the Wendland C6 kernel from Dehnen & Aly 2012
#WENDLAND_BIAS_CORRECTION                     # reduces self-contribution for Wendland kernels according to Dehnen & Aly 2012


#--------------------------------------- SPH viscosity options

#TIMEDEP_ART_VISC                             # enables time dependend viscosity
#HIGH_ART_VISC_START                          # start with high rather than low viscosity
#NO_SHEAR_VISCOSITY_LIMITER                   # turns of the shear viscosity supression


#--------------------------------------- Extra physics

#COOLING                                      # Enables radiative atomic cooling by hydrogen and helium
#STARFORMATION                                # Enables star formation with the Springel & Hernquist (2003) model


#---------------------------------------- Single/double precision and data types

#POSITIONS_IN_32BIT                           # if set, use 32-integers for positions  (default for single precision)
#POSITIONS_IN_64BIT                           # if set, use 64-integers for positions  (default for double precision)
#POSITIONS_IN_128BIT                          # if set, use 128-integers for positions (may only work on some systems)
DOUBLEPRECISION=1                             # if activated and set to 1, use double precision internally, for 2 use mixed precision, otherwise single precision
#DOUBLEPRECISION_FFTW                         # if set, carries out FFTs in double precision, otherwise single precision
#USE_SINGLEPRECISION_INTERNALLY               # reduces default double precision for most internal computations to single precision
#ENLARGE_DYNAMIC_RANGE_IN_TIME                # eextend the dynamic range of the integer timeline from 32 to 64 bits
#IDS_32BIT                                    # selects 32-bit IDs for internal storage (default)
#IDS_48BIT                                    # selects 48-bit IDs for internal storage
#IDS_64BIT                                    # selects 64-bit IDs for internal storage
#OUTPUT_IN_DOUBLEPRECISION                    # output files will be written in double precision


#---------------------------------------- Output/Input options

INITIAL_CONDITIONS_CONTAIN_ENTROPY
#OUTPUT_VELOCITY_GRADIENT                     # output velocity gradients
#OUTPUT_PRESSURE                              # output gas pressure   
#OUTPUT_ENTROPY                               # output gas entropy
#OUTPUT_CHANGEOFENTROPY                       # output rate of change of entropy
#OUTPUT_POTENTIAL                             # output gravitational potential
#OUTPUT_ACCELERATION                          # output total acceleration
#OUTPUT_TIMESTEP                              # output timesteps
#OUTPUT_DIVVEL                                # output velocity divergence
#OUTPUT_CURLVEL                               # output velocity curl
#OUTPUT_COOLHEAT                              # output actual energy loss/gain in cooling/heating routine
#OUTPUT_PRESSURE_SPH_DENSITY                  # output also density computed in the pressure-entropy forumulation
#OUTPUT_VISCOSITY_PARAMETER                   # output current viscosity parameter in SPH particle output (only for time-dependent viscosity)
#OUTPUT_NON_SYNCHRONIZED_ALLOWED              # allow snapshot creation also at steps that are not fully synchronized
#OUTPUT_VELOCITIES_IN_HALF_PRECISION          # special option to store velocities in reduced precision
#OUTPUT_ACCELERATIONS_IN_HALF_PRECISION       # special option to store acclerations in reduced precision
#OUTPUT_COORDINATES_AS_INTEGERS               # special option to store coordinates as integers that are used internally
#POWERSPEC_ON_OUTPUT                          # computes a matter power spectrum when the code writes a snapshot output
#ALLOW_HDF5_COMPRESSION                       # applies HDF5 loss-less compression to selected output fields
#REDUCE_FLUSH                                 # do not flush the I/O streams of the log-files every system step


#---------------------------------------- On the fly FOF groupfinder

#FOF                                          # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2                     # bitmask, 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32             # bitmask, 2^type for the type(s) linked to nearest primary particle
#FOF_GROUP_MIN_LEN=32                         # minimum group length (default is 32)
#FOF_LINKLENGTH=0.2                           # linking length for FoF (default=0.2)
#FOF_ALLOW_HUGE_GROUPLENGTH                   # set if there are groups and subhalos with more than 2 billion particles in length


#---------------------------------------- Subfind

#SUBFIND                                      # enables substructure finder
#SUBFIND_HBT                                  # use previous subhalo catalogue instead of density excursions to define subhalo candidates
#SUBFIND_STORE_LOCAL_DENSITY                  # calculates local densities and velocity dispersions for all particles and stores them in snapshots
#SUBFIND_ORPHAN_TREATMENT                     # creates special snapshots with formerly most bound particles


#---------------------------------------- Merger tree code

#MERGERTREE                                   # enables on-the-fly calculation of descendants, and merger tree construction from group catalogues


#---------------------------------------- On-the-fly lightcone creation

#LIGHTCONE                                    # master option for lightcone output creation
#LIGHTCONE_PARTICLES                          # produces particle lightcones
#LIGHTCONE_MASSMAPS                           # produces mass shells on the lightcone
#LIGHTCONE_PARTICLES_GROUPS                   # computes groups for particles buffered on the lightcone
#LIGHTCONE_PARTICLES_SKIP_SAVING              # prevents that particle data is saved along with the found groups on the lightcone
#LIGHTCONE_OUTPUT_ACCELERATIONS               # stores accelerations for particles on lightcone
#LIGHTCONE_IMAGE_COMP_HSML_VELDISP            # option for computing densities and smoothing length for lightcones in postprocessing
#LIGHTCONE_MULTIPLE_ORIGINS                   # switch this on if you want to be able to define lightcone origins different from (0,0,0)
#REARRANGE_OPTION                             # special option to reorder lightcone data in mergertree order


#--------------------------------------- IC creation
 
#NGENIC=256                                   # generate cosmological ICs, set NGENIC to the FFT grid size used for IC generation
#NGENIC_2LPT                                  # applies 2LPT instead of just Zeldovich approximation
#CREATE_GRID                                  # start with a regular Cartesian DM particle grid, instead of reading a glass file (for NGENIC)
#GENERATE_GAS_IN_ICS                          # add SPH particles to created or read dark matter only ICs
#SPLIT_PARTICLE_TYPE=4+8                      # specifies particle types to be split if GENERATE_GAS_IN_ICS is activated
#NGENIC_FIX_MODE_AMPLITUDES                   # when activated, this leaves the mode amplitudes at sqrt(P(k)), instead of sampling from a Rayleigh distribution
#NGENIC_MIRROR_PHASES                         # if this is activated, all phases are turned by 180 degrees
#NGENIC_TEST                                  # can be used to create ICs, measure the power spectrum, and then stop


#----------------------------------------Parallelization options

#IMPOSE_PINNING                               # enables pinning of MPI processes to CPU cores
#IMPOSE_PINNING_OVERRIDE_MODE                 # tries to do the pinning even if a prior pinning is detected
#PRESERVE_SHMEM_BINARY_INVARIANCE             # preserve binary invariance of results despite machine weather, at the price of more tree walk overhead 
#EXPLICIT_VECTORIZATION                       # use AVX at selected places in SPH kernels through the vectorclass C++ library
#SIMPLE_DOMAIN_AGGREGATION                    # this is an experimental modification of the domain decomposition algorithm (can either help or harm performance)


#---------------------------------------- MPI related settings

#NUMBER_OF_MPI_LISTENERS_PER_NODE=1           # set such that the number of MPI-ranks per node and listener is maller than MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY 
#MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY=64    # default is 64, but can also be set to 32
#NUMPART_PER_TASK_LARGE                       # set this if the number of particles per task is so large that more than 2 GB are comprised just by particle data
#USE_MPIALLTOALLV_IN_DOMAINDECOMP             # replaces hypercube communication in domain particle exchance with a single MPI_Allgatherv (can be less stable)
#MPI_HYPERCUBE_ALLGATHERV                     # if your MPI-library uses too much internal storage for MPI_Allgatherv, this uses a hypercube as a work-around
#MPI_MESSAGE_SIZELIMIT_IN_MB=200              # limit the message size of very large MPI transfers
#MPI_HYPERCUBE_ALLTOALL                       # use a robust hyercube for MPI_Alltoall instead the native algorithm if the MPI library
#ISEND_IRECV_IN_DOMAIN                        # uses asynchronous communication instead of synchronous communication in hypercube pattern (can be less stable)
#ALLOCATE_SHARED_MEMORY_VIA_POSIX             # if this is set, do use POSIX directly to allocated shared memory instead of MPI-3 calls
#OLDSTYLE_SHARED_MEMORY_ALLOCATION            # disables new memory allocation mechanism via memfd_create()


#---------------------------------------- Testing and Debugging options

#HOST_MEMORY_REPORTING                        # reports after start-up the available system memory by analyzing /proc/meminfo (active by default on Linux)
#ENABLE_HEALTHTEST                            # can be used to enable some auto-checking of the machine you're running on
#DEBUG                                        # enables core-dumps in case MPI_Init() has disabled them
#DEBUG_ENABLE_FPU_EXCEPTIONS                  # tries to enable FPU exceptions so that a bad floating point operation triggers a core dump
#DEBUG_SYMTENSORS                             # carries out a few unit tests related to the tensor algebra routines of the code
#FORCETEST=0.001                              # calculates for given fraction of particles direct summation forces to check accuracy of tree force
#FORCETEST_FIXEDPARTICLESET                   # if this set, always the same particle receive a force accuracy check during a run
#FORCETEST_TESTFORCELAW=1                     # special option for measuring the effective force law, can be set to 1 or 2 for TreePM+HighRes
#VTUNE_INSTRUMENT                             # outputs auxiliary info to simplify Vtune performance analysis
#DEBUG_MD5                                    # make MD5 check sum routines available for debugging and code development
#SQUASH_TEST                                  # special check-option for code development
#DOMAIN_SPECIAL_CHECK                         # special check-option for code development
#RECREATE_UNIQUE_IDS                          # before carrying out the ID uniqueness test, this sets new IDs (e.g. to fix broken ICs)
#NO_STOP_BELOW_MINTIMESTEP                    # do not stop when the code wants to go below a prescribed minimum timestep
#TILING=2                                     # duplicate an IC in each dimension for scaling tests
#DO_NOT_PRODUCE_BIG_OUTPUT                    # for scaling tests, disable the writing of restart files and snapshot dumps
#MEASURE_TOTAL_MOMENTUM                       # special test option for code development
#EWALD_TEST                                   # a development test for the Ewald tables
#STOP_AFTER_STEP=10                           # ends a simulation after the specified timestep (to simplify scalability tests)
#TREE_NO_SAFETY_BOX                           # when set, this disables the geometric 'near node' protection

