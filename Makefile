#/*******************************************************************************
# * \copyright   This file is part of the GADGET4 N-body/SPH code developed
# * \copyright   by Volker Springel. Copyright (C) 2014, 2015 by Volker Springel
# * \copyright   (volker.springel@h-its.org) and all contributing authors.
# *******************************************************************************/
#
# You might be looking for the compile-time Makefile options of the
# code if you are familiar with Gadget2...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the
#  file "Template-Config.sh" only. Usually, they should be added
#  there in the disabled/default version.
#
#  "Config.sh" should not be checked in to the repository.
#
#  Note: It is possible to override the default name of the
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make CONFIG=MyNewConf.sh EXEC=Gadget4_new
#
#-----------------------------------------------------------------
#
# You might also be looking for the target system SYSTYPE option
#
# It has also moved to a separate file.
#
# To build the code, do the following:
#
# (A) set the SYSTYPE variable in your .bashrc (or similar file):
#
#        e.g. export SYSTYPE=Magny
# or
#
# (B) set SYSTYPE in Makefile.systype
#     This file has priority over your shell variable:
#
#    (1) Copy the file "Template-Makefile.systype"  to  "Makefile.systype"
#
#        cp Template-Makefile.systype Makefile.systype
#
#    (2) Uncomment your system in  "Makefile.systype".
#
# For the chosen system type, an if-clause should be defined below,
# loading short definitions of library path names and/or compiler
# names and options from the buildsystem/ directory. A new system
# type should also be added to Template-Makefile.systype
#


ifdef DIR
EXEC = $(DIR)/Gadget4
CONFIG = $(DIR)/Config.sh
BUILD_DIR = $(DIR)/build
else
EXEC   = Gadget4
CONFIG   = Config.sh
BUILD_DIR = build
endif


SRC_DIR = src

###################
#determine SYSTYPE#
###################
ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif




$(info Build configuration:)
$(info SYSTYPE: $(SYSTYPE))
$(info CONFIG: $(CONFIG))
$(info EXEC: $(EXEC))
$(info )


PYTHON   = python

RESULT     := $(shell CONFIG=$(CONFIG) PYTHON=$(PYTHON) BUILD_DIR=$(BUILD_DIR) SRC_DIR=$(SRC_DIR) CURDIR=$(CURDIR) make -f buildsystem/Makefile.config)
$(info $(RESULT))
CONFIGVARS := $(shell cat $(BUILD_DIR)/gadgetconfig.h)

RESULT     := $(shell SRC_DIR=$(SRC_DIR) BUILD_DIR=$(BUILD_DIR) ./buildsystem/git_version.sh)

##########################
#define available Systems#
##########################
ifeq ($(SYSTYPE),"Generic-gcc")
include buildsystem/Makefile.gen.libs
include buildsystem/Makefile.comp.gcc
endif
ifeq ($(SYSTYPE),"Generic-intel")
include buildsystem/Makefile.comp.gcc-paranoia
include buildsystem/Makefile.gen.libs
endif

ifeq ($(SYSTYPE),"SuperMUC-NG")
include buildsystem/Makefile.comp.supermuc-ng
include buildsystem/Makefile.path.supermuc-ng
endif

ifeq ($(SYSTYPE),"SuperMUC-NG-OpenMPI")
include buildsystem/Makefile.comp.supermuc-ng-openmpi
include buildsystem/Makefile.path.supermuc-ng
endif

ifeq ($(SYSTYPE),"SuperMUC-NG-GCC")
include buildsystem/Makefile.comp.supermuc-ng-gcc
include buildsystem/Makefile.path.supermuc-ng-gcc
endif

ifeq ($(SYSTYPE), "Generic-gcc-single")
include buildsystem/Makefile.comp.gcc
include buildsystem/Makefile.gen.libs
endif

ifeq ($(SYSTYPE), "Generic-intel-single")
include buildsystem/Makefile.comp.gcc-paranoia
include buildsystem/Makefile.gen.libs
endif

ifeq ($(SYSTYPE),"Darwin")
include buildsystem/Makefile.comp.gcc
include buildsystem/Makefile.path.macports
endif

ifeq ($(SYSTYPE),"Magny")
include buildsystem/Makefile.comp.gcc
include buildsystem/Makefile.path.magny
endif

ifeq ($(SYSTYPE),"Freya")
include buildsystem/Makefile.comp.freya
include buildsystem/Makefile.path.freya
endif

#module load  gcc/7.2    gsl/2.2  hdf5-serial/gcc/1.8.18   fftw/gcc/3.3.6   
ifeq ($(SYSTYPE),"FreyaOpenMPI")
include buildsystem/Makefile.comp.freyaopenmpi
include buildsystem/Makefile.path.freya
endif


ifeq ($(SYSTYPE),"Cobra")
include buildsystem/Makefile.comp.cobra
include buildsystem/Makefile.path.cobra
endif

ifeq ($(SYSTYPE),"RavenOpenMPI")
include buildsystem/Makefile.comp.ravenopenmpi
include buildsystem/Makefile.path.cobra
endif

ifeq ($(SYSTYPE),"CobraOpenMPI")
include buildsystem/Makefile.comp.cobraopenmpi
include buildsystem/Makefile.path.cobra
endif

ifeq ($(SYSTYPE),"Haswell")
include buildsystem/Makefile.comp.gcc
include buildsystem/Makefile.path.haswell
endif

ifeq ($(SYSTYPE),"gcc-paranoia")
include buildsystem/Makefile.comp.gcc-paranoia
include buildsystem/Makefile.path.mpa_desktop
endif


ifeq ($(SYSTYPE),"libs")
include buildsystem/makefile.comp.gcc
include buildsystem/Makefile.path.libs
endif

ifeq ($(SYSTYPE),"hydra")
include buildsystem/Makefile.comp.gcc
include buildsystem/Makefile.path.hydra
endif

ifeq ($(SYSTYPE),"bwforcluster")
include buildsystem/Makefile.comp.gcc
include buildsystem/Makefile.path.bwforcluster
endif


ifndef LINKER
LINKER = $(CPP)
endif

##########################################
#determine the needed object/header files#
##########################################
SUBDIRS += .


SUBDIRS += main
OBJS    += main/begrun.o main/init.o main/main.o main/run.o
INCL    += main/main.h main/simulation.h


SUBDIRS += data
OBJS    += data/mymalloc.o data/allvars.o data/test_symtensors.o
INCL    += data/allvars.h data/dtypes.h data/mymalloc.h data/idstorage.h data/symtensors.h \
           data/intposconvert.h data/constants.h data/simparticles.h \
           data/macros.h data/particle_data.h data/sph_particle_data.h \
           data/lightcone_particle_data.h data/lightcone_massmap_data.h data/lcparticles.h data/mmparticles.h


SUBDIRS += domain
OBJS    += domain/domain.o domain/domain_balance.o domain/domain_box.o \
           domain/domain_exchange.o domain/domain_toplevel.o
INCL    += domain/domain.h


SUBDIRS += io
OBJS    += io/hdf5_util.o io/snap_io.o io/parameters.o \
           io/restart.o io/io.o io/test_io_bandwidth.o
INCL    += io/io.h io/hdf5_util.h io/snap_io.h io/parameters.h \
	         io/restart.h io/io_streamcount.h io/test_io_bandwidth.h


SUBDIRS += logs
OBJS    += logs/logs.o
INCL    += logs/logs.h logs/timer.h


SUBDIRS += gitversion
INCL    += gitversion/version.h


SUBDIRS += mpi_utils
OBJS    += mpi_utils/hypercube_allgatherv.o mpi_utils/mpi_types.o mpi_utils/mpi_vars.o mpi_utils/sums_and_minmax.o \
           mpi_utils/sizelimited_sendrecv.o mpi_utils/myalltoall.o mpi_utils/allreduce_sparse_double_sum.o mpi_utils/healthtest.o \
           mpi_utils/allreduce_debugcheck.o mpi_utils/shared_mem_handler.o 
INCL    += mpi_utils/mpi_utils.h mpi_utils/generic_comm.h mpi_utils/shared_mem_handler.h


SUBDIRS += pm
OBJS    += pm/pm_nonperiodic.o pm/pm_periodic.o \
           pm/pm_mpi_fft.o
INCL    += pm/pm.h pm/pm_mpi_fft.h pm/pm_periodic.h pm/pm_nonperiodic.h


SUBDIRS += sort
OBJS    += sort/peano.o
INCL    += sort/peano.h sort/cxxsort.h sort/parallel_sort.h


SUBDIRS += sph
OBJS    += sph/density.o sph/hydra.o sph/init_entropy.o sph/artificial_viscosity.o
INCL    += sph/kernel.h sph/sph.h


SUBDIRS += system
OBJS    += system/pinning.o system/system.o
INCL    += system/system.h system/pinning.h


SUBDIRS += time_integration
OBJS    += time_integration/driftfac.o time_integration/kicks.o \
           time_integration/predict.o time_integration/timestep.o \
           time_integration/timestep_treebased.o
INCL    += time_integration/timestep.h time_integration/driftfac.h


SUBDIRS += gravity
OBJS    += gravity/gravity.o gravity/ewald.o gravity/ewald_test.o \
           gravity/grav_forcetest.o gravity/grav_external.o \
           gravity/grav_direct.o gravity/second_order_ics.o
INCL    += gravity/ewald.h gravity/ewaldtensors.h gravity/grav_forcetest.h


SUBDIRS += tree
OBJS    += tree/tree.o
INCL    += tree/tree.h


SUBDIRS += gravtree
OBJS    += gravtree/gravtree_build.o gravtree/gravtree.o gravtree/gwalk.o
INCL    += gravtree/gravtree.h  gravtree/gwalk.h  


SUBDIRS += ngbtree
OBJS    += ngbtree/ngbtree_build.o 
INCL    += ngbtree/ngbtree.h


ifeq (EXPLICIT_VECTORIZATION,$(findstring EXPLICIT_VECTORIZATION,$(CONFIGVARS)))
SUBDIRS += vectorclass
OBJS    += vectorclass/instrset_detect.o
INCL    +=
endif


ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
OBJS    += cooling_sfr/cooling.o cooling_sfr/sfr_eos.o cooling_sfr/starformation.o
INCL    += cooling_sfr/cooling.h
SUBDIRS += cooling_sfr
endif


ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
OBJS    += fof/fof.o fof/fof_findgroups.o fof/fof_nearest.o fof/fof_io.o fof/foftree_build.o
INCL    += fof/fof.h  fof/fof_io.h  fof/foftree.h
SUBDIRS += fof
endif


ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
OBJS	+= subfind/subfind.o subfind/subfind_treepotential.o \
           subfind/subfind_processing.o subfind/subfind_density.o subfind/subfind_distribute.o subfind/subfind_findlinkngb.o \
           subfind/subfind_nearesttwo.o subfind/subfind_properties.o subfind/subfind_unbind.o subfind/subfind_history.o \
           subfind/subfind_so.o subfind/subfind_readid_io.o subfind/subfind_orphanids.o subfind/subfind_excursionset.o 
INCL	+= subfind/subfind.h subfind/subfind_readid_io.h
SUBDIRS += subfind
endif


ifeq (FMM,$(findstring FMM,$(CONFIGVARS)))
SUBDIRS += fmm
OBJS    += fmm/fmm.o
INCL    += fmm/fmm.h
endif


ifeq (MERGERTREE,$(findstring MERGERTREE,$(CONFIGVARS)))
OBJS    += mergertree/descendant.o mergertree/io_descendant.o mergertree/io_progenitors.o \
           mergertree/postproc_descendants.o mergertree/io_readsnap.o mergertree/halotrees.o \
           mergertree/io_halotrees.o mergertree/io_treelinks.o mergertree/io_readtrees_mbound.o \
           mergertree/rearrange.o
INCL    += mergertree/mergertree.h mergertree/io_descendant.h mergertree/io_progenitors.h \
           mergertree/io_readsnap.h mergertree/io_treelinks.h mergertree/io_halotrees.h \
           mergertree/io_readtrees_mbound.h
SUBDIRS += mergertree
endif


ifeq (LIGHTCONE,$(findstring LIGHTCONE,$(CONFIGVARS)))
SUBDIRS += lightcone
OBJS    += lightcone/lightcone.o lightcone/lightcone_particle_io.o lightcone/lightcone_massmap_io.o      
INCL    += lightcone/lightcone.h lightcone/lightcone_particle_io.h lightcone/lightcone_massmap_io.h 
endif


ifeq (LIGHTCONE_MASSMAPS,$(findstring LIGHTCONE_MASSMAPS,$(CONFIGVARS)))
MAPS_LIBS += -lchealpix -lcfitsio #-lcurl
endif


ifeq (LIGHTCONE_PARTICLES,$(findstring LIGHTCONE_PARTICLES,$(CONFIGVARS)))
MAPS_LIBS += -lchealpix -lcfitsio #-lcurl
endif


ifeq (NGENIC,$(findstring NGENIC,$(CONFIGVARS)))
OBJS    += ngenic/ngenic.o ngenic/power.o ngenic/grid.o
INCL    += ngenic/ngenic.h
SUBDIRS += ngenic
endif


ifeq (DEBUG_MD5,$(findstring DEBUG_MD5,$(CONFIGVARS)))
SUBDIRS += debug_md5
OBJS    += debug_md5/calc_checksum.o debug_md5/Md5.o
INCL    += debug_md5/Md5.h
endif


################################
#determine the needed libraries#
################################

# we only need fftw if PMGRID is turned on
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
FFTW_LIBS += -lfftw3
else
FFTW_LIBS += -lfftw3f
endif
else
ifeq (NGENIC, $(findstring NGENIC, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
FFTW_LIBS += -lfftw3
else
FFTW_LIBS += -lfftw3f
endif
else

endif
endif

ifeq (FORCETEST, $(findstring FORCETEST, $(CONFIGVARS)))
FFTW_LIBS += -lfftw3
endif


HWLOC_LIBS  += -lhwloc

ifneq (IMPOSE_PINNING,$(findstring IMPOSE_PINNING,$(CONFIGVARS)))
HWLOC_INCL =
HWLOC_LIBS =
endif

ifneq (VTUNE_INSTRUMENT,$(findstring VTUNE_INSTRUMENT,$(CONFIGVARS)))
VTUNE_INCL =
VTUNE_LIBS =
endif

GSL_LIBS   += -lgsl -lgslcblas
HDF5_LIBS  += -lhdf5 -lz
MATH_LIBS  = -lm

ifneq ($(SYSTYPE),"Darwin")
ifeq (ALLOCATE_SHARED_MEMORY_VIA_POSIX,$(findstring ALLOCATE_SHARED_MEMORY_VIA_POSIX,$(CONFIGVARS)))
SHMEM_LIBS  = -lrt
endif
endif

MAKEFILES = $(MAKEFILE_LIST) buildsystem/Makefile.config

##########################
#combine compiler options#
##########################

CFLAGS = $(OPTIMIZE) $(OPT) $(HDF5_INCL) $(GSL_INCL) $(FFTW_INCL) $(HWLOC_INCL) $(VTUNE_INCL) $(MAPS_INCL) -I$(BUILD_DIR) -I$(SRC_DIR)

LIBS = $(MATH_LIBS) $(HDF5_LIBS) $(GSL_LIBS) $(FFTW_LIBS) $(HWLOC_LIBS) $(VTUNE_LIBS) $(TEST_LIBS) $(MAPS_LIBS) $(SHMEM_LIBS)


SUBDIRS := $(addprefix $(BUILD_DIR)/,$(SUBDIRS))
OBJS := $(addprefix $(BUILD_DIR)/,$(OBJS)) $(BUILD_DIR)/compile_time_info.o $(BUILD_DIR)/compile_time_info_hdf5.o $(BUILD_DIR)/version.o
INCL := $(addprefix $(SRC_DIR)/,$(INCL)) $(BUILD_DIR)/gadgetconfig.h


TO_CHECK := $(addsuffix .check, $(OBJS) $(patsubst $(SRC_DIR)%, $(BUILD_DIR)%, $(INCL)) )
TO_CHECK +=  $(BUILD_DIR)/Makefile.check
CONFIG_CHECK = $(BUILD_DIR)/$(notdir $(CONFIG)).check
DOCS_CHECK = $(BUILD_DIR)/config.check  $(BUILD_DIR)/param.check


################
#create subdirs#
################
RESULT := $(shell mkdir -p $(SUBDIRS)  )

###########################################
#create info file for command line options#
###########################################
RESULT := $(shell echo 'static const char *compiler_flags="$(CPP) $(CFLAGS)";' > $(BUILD_DIR)/compiler-command-line-args.h )

#############
#build rules#
#############

all: check_docs check build

build: $(EXEC)

$(EXEC): $(OBJS)
	$(LINKER) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -f $(OBJS) $(EXEC)
	rm -f $(BUILD_DIR)/compile_time_info.cc $(BUILD_DIR)/compile_time_info_hdf5.cc $(BUILD_DIR)/gadgetconfig.h
	rm -f $(TO_CHECK) $(CONFIG_CHECK)
	rm -f $(BUILD_DIR)/version.cc

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCL) $(MAKEFILES)
	$(CPP) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc $(INCL) $(MAKEFILES)
	$(CPP) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_time_info.o: $(BUILD_DIR)/compile_time_info.cc $(MAKEFILES)
	$(CPP) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_time_info_hdf5.o: $(BUILD_DIR)/compile_time_info_hdf5.cc $(MAKEFILES)
	$(CPP) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/version.o: $(BUILD_DIR)/version.cc $(MAKEFILES)
	$(CPP) $(CFLAGS) -c $< -o $@

check: $(CONFIG_CHECK)

check_docs: $(DOCS_CHECK)

$(CONFIG_CHECK): $(TO_CHECK) $(CONFIG) buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 2 $(CONFIG) $(CONFIG_CHECK) defines_extra $(TO_CHECK)

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.cpp Template-Config.sh defines_extra buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.cc Template-Config.sh defines_extra buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.h.check: $(SRC_DIR)/%.h Template-Config.sh defines_extra buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.o.check: $(BUILD_DIR)/%.cc Template-Config.sh defines_extra buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.h.check: $(BUILD_DIR)/%.h Template-Config.sh defines_extra buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/Makefile.check: Makefile Template-Config.sh defines_extra buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 3 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/config.check: documentation/04_config-options.md Template-Config.sh buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 4 Template-Config.sh $@  $<

$(BUILD_DIR)/param.check: documentation/05_parameterfile.md $(SRC_DIR)/io/parameters.cc buildsystem/check.py
	@$(PYTHON) buildsystem/check.py 5 $(SRC_DIR)/data/allvars.cc $@  $<

.PHONY = all check build clean
