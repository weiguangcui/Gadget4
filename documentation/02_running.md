
Compilation and basic usage                            {#usage}
===========================

[TOC]

If you are already familiar with older versions of GADGET or with
AREPO, you will feel quickly at home in compiling and working with
GADGET-4, since the code follows very similar usage concepts. However,
GADGET-4 is now a C++ code, and has a more elaborate build concept
that in particular attempts to prevent basic forms of code rot. This
makes the build process slightly more advanced, as described below.



Compilation requirements                              {#compreq}
========================


GADGET-4 needs the following non-standard libraries for compilation:

1. **mpi**: The `Message Passing Interface' (version 3.0 or higher).
   Many vendor supplied versions exist, in addition to excellent open
   source implementations, e.g. MPICH <https://www.mpich.org>, or
   OpenMPI <http://www.open-mpi.org>.

2. **gsl**: The GNU scientific library. This open-source package can
   be obtained at <http://www.gnu.org/software/gsl>. GADGET-4 only
   needs this library for a few very simple cosmological integrations
   at start-up.

3. **fftw3**: The `Fastest Fourier Transform in the West'. This
   open-source package can be obtained at <http://www.fftw.org>. Note
   that the MPI-capable version of FFTW-3 is not explicitly required
   (it makes no difference whether it is available or not as GADGET-4
   implements its own communication routines when MPI is used). FFTW
   is only needed for simulations that use the TreePM algorithm, or if
   power spectra are estimated, or cosmological ICs are created.

4. **hdf5**: The `Hierarchical Data Format' (available at
   <http://hdf.ncsa.uiuc.edu/HDF5>. This library is needed when one
   wants to read or write snapshot files in HDF5 format.  It is highly
   recommended to use HDF5 instead of the plain binary output format
   that can be used as an alternative and still available for
   historical reasons.

5. **hwloc** : The `hardware locality library' is useful for allowing
   the code to probe the processor topology it is running on and
   enable a pinning to individual cores. This library is optional and
   only required if the IMPOSE_PINNING flag is set. Note that many MPI
   libraries nowadays in any cases enable pinning by default.

6. **vectorclass** : This is a C++ library that the code utilizes when
   explicit vectorization via the AVX instruction set is enabled in
   the SPH compute kernels. This is then implemented with the
   vectorclass library by Agner Fog. For simplicity, the source code
   of this library (which takes the form of C++ header files) is
   included in the GADGET-4 distribution in directory
   src/vectorclass. Note that potential updates and bug fixed in this
   library require a corresponding update of this subdirectory.

Compilation of GADGET-4 needs a working C++ compiler, supporting at
least the C++11 standard. For GNU gcc, this means version 4.x or
later.  The code also makes use of GNU-Make and Python as part of its
build process.  Those should hence be available as well.



Building the code                                          {#building}
=================


After obtaining GADGET-4 from the code repository, you will find a
bunch of subdirectories as well as few further files for the build
system in the top-level directory. The most important subdirectory is
`src/` , which contains the actual source code files, distributed into
subdirectories according to the functionality of each file.

The GADGET-4 code is controlled and configured by two different files,
one containing compile-time options, and one listing runtime
parameters. The default name for the file with compile-time options is
`Config.sh` , and the meaning of the different options specified there
is explained in a special section of this manual. The run-time options
of GADGET-4 are controlled by a parameter file, and are described also
in a separate section of this documentation.

One possible way to create the file `Config.sh` needed for compilation
is to make a copy of `Template-Config.sh` and then modify it as needed
by commenting in or commenting out the desired options. Given the
length of this template file, better overview is provided by
assembling only the enabled options in a file (which can also be
gleaned from a previous GADGET-4 run). Another option is to use one of
the Config.sh files that come with the example problems, and modify it
if needed. The code uses the GNU make utility for controlling the
build process, which is specified by the file `Makefile`.  Typing
`make` will attempt to build the executable `Gadget4` (hint: using
`make -j X`, where X is some number of threads, the build process can
be carried much faster in parallel on multi-core machines). It is also
possible to override the default names for the configuration file and
the executable by specifying them as parameters for this command, for
example as `make CONFIG=MyNewConf.sh EXEC=Gadget4new`. Also, one may
pass a directory to make, for example in the form `make
DIR=mysim/newmodels/runA`, and then the configuration file contained
in this directory will be used, the build process will be carried out
there, and the executable appears there as well. This can be very
useful to build and run different code configurations from a common
code base without risking to accidentally mix up the executables. Each
simulation should be organized into a different subdirectory in this
case, with configuration file and executable being placed into the
same subdirectory.

Often, one needs to specify (possibly non-standard) locations of the
libraries required by GADGET-4, or the name of the compiler and the
settings of compiler flags to be used in order to build the code. To
make the selection of the target system and changes between different
compiler setups relatively simple, these settings have been largely
decoupled from the `Makefile` , i.e. the `Makefile` itself should
usually not be changed in any significant way. Only if you want to add
additional source files to the code or introduce a new target computer
system a small change is needed there.

This is facilitated through the concept of a system type (which is the
target computer combined with compiler settings), which is selected
through the environment variable `SYSTYPE`, or by the file
`Makefile.systype`. The simplest way to create the file
`Makefile.systype` is to copy it from `Template-Makefile.systype`, and
then comment in/out the desired type. For the selected symbolic name
of the target system type, there should then be a short section in the
`Makefile` that includes usually two makefile stubs from the folder
`buildsystem/` which specify the paths to the different libraries if
they are not in standard locations, and determine the compiler name(s)
and compiler options that should be used. Once this is all set-up
properly, you can then switch to a different compiler setting by
toggling one line in `Makefile.systype`. Also, if you set the
environment variable `SYSTYPE` in your login-script (.profile,
.bashrc, etc.)  on your target computer (recommended for convenience),
you will normally not have to deal with the file `Makefile.systype` at
all and can compile right away.

To summarize, if you want to set-up GADGET-4 for compilation on a new,
not yet defined computer system, you have to go through the following
steps:

- Make sure that the needed libraries are available, either by
  compiling them yourself (they can all be installed in user space if
  needed) or by loading the right modules on systems that support the
  module-environment.

- Add a new symbolic name for your computer system to
  `Template-Makefile.systype` (this file is meant to keep a record of
  all defined machine types), and then select this name through the
  `SYSTYPE` environment variable or through the file
  `Makefile.systype`.

- In the file `Makefile`, create an if-clause in the "define available
  Systems" section, in which you include two short files from the
  `buildsystem` directory that define path names to the libraries that
  you need (i.e. if not in default locations), and that specify the
  compiler and its flags that you want to select. Follow the provided
  examples to create your own versions of these files, as needed. The
  reason why this information is normally not directly included in an
  if-clause in the `Makefile` (which would be possible too) is to
  avoid repetition of identical information in the Makefile (for
  example, because the same settings for the `gcc` compiler may be
  used on many different target systems). Instead, the same small
  sub-file from `buildsystem/` can be included for many different
  target systems.


We note that many compile-time options are introduced with `#ifdef` /
`#endif` statements in the code in order to encapsulate the
corresponding code extensions and allow them to be fully disabled
(i.e.~to be **completely** eliminated from the compiled code). This is
done for performance and memory reasons, which for us take precedence
over the detrimental impact on readability/clarity of the code brought
about by such compiler pragmas. Another, more problematic side effect
of using such compile-time symbols is however that typos in specifying
any of them may go easily undetected. To provide protection against
this, the code automatically runs a python-based check script over the
source code as part of the build process that will complain if any
undefined or misspelled compile time symbols are listed in
`Config.sh`, or, conversely, if (new) symbols in the source code are
present that are neither defined nor briefly explained in
`Template-Config.sh`.  This checking mechanism can be disabled if
desired by using `make build` instead of `make`, but we strongly
advise against this. Also, it is possible to deliberately exempt a
symbol from the checking procedure. This is needed, for example, for
header-file guards, and those symbols should be added to the file
`defines_extra`.

As a further extension of these checks, we have also added
functionality that demands that all compile-time and run-time options
are actually documented. To this end, if you add a new compile time
option yourself, this needs to be documented in
`documentation/04_config-options.md`, and if you add a new run-time
parameter, it needs to be documented in
`documentation/05_parameterfile.md`, otherwise the code will refuse to
compile and complain accordingly. Similarly, documented options that
do no longer exist in the code will lead to error messages when
compiling. This checking of the documentation can be disabled by using
`make build` instead of `make`, but as pointed out above, this
constitutes bad practice and should not be done.

Finally, note that the GADGET-4 code is now written in the C++
language, and the source files can only be compiled with a C++
compiler. While many useful and advanced features of the C++ language
are used (like templating and operator overloading), because GADGET-4
evolved from an older code base in C, it effectively represents in
many places a mixture of C++ and C-styles of coding.


Starting the code                                            {#starting}
=================

To start a simulation, invoke the executable with a command like

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpirun -np 32  ./Gadget4  param.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will start the simulation using 32 MPI ranks, and with simulation
parameters as specified in the parameter file (see below) of name
param.txt.  Note that on some systems, the command to launch parallel
execution may have a different name or syntax (e.g mpiexec, srun, poe,
etc). Consult the man-pages or the local support if in doubt.

The code does not need to be recompiled for a different number of
processors, or for a different problem size. It is also possible to
run GADGET-4 with a single processor only. In this case, the leading
`mpirun -np 1` can normally be omitted, and GADGET-4 will behave like
a serial code. It is still necessary to have MPI installed, though,
because even then the code will make some calls to the MPI library
(but none that actually do non-trivial communications). There is no
restriction for the processor number to be a power of two, even though
these partition sizes are slightly preferred, because they allow the
most efficient communication schemes. However, in regimes where the
code works well in terms of scalability, the communication times
should is subdominant anyhow, so this is not an important
consideration in practice.

When more than one MPI-rank is used, the code will use a hybrid
communication scheme in which data stored by different MPI processes
on the same compute node are accessed directly via shared-memory. The
code automatically detects groups of MPI ranks running on the same
node. If more than one node is in use, at least one MPI process on
each node is set aside for asynchronously serving incoming
communication requests from other nodes (if only a single
shared-memory is used, this is not done). This means that multi-node
jobs must have a minimum of two MPI ranks on each node. Today's
machines offer typically many more cores per node than 2, and their
full power is made available to GADGET-4 if one MPI rank is placed on
every core.

At least one MPI communication rank needs to be set aside for every 64
ranks on a shared memory node in multi-node jobs. If the number of MPI
ranks per shared memory node is larger than 64, one therefore needs to
specify the `NUMBER_OF_MPI_LISTENERS_PER_NODE=X` option, with X larger
than 1 (which is the default).

While GADGET-4 is running, it will print out many log-messages that
inform about the present steps taken by the code. When you start a
simulation interactively, these log-messages will appear on the
screen, but you can also redirect them to a file. For production runs
on a cluster controlled by a queuing system, you will usually have to
put the above start-up command of GADGET-4 into a batch script-file
that is submitted to the queuing system. In this case, the standard
output of GADGET-4 is usually automatically piped into a file.

For example, assuming that the batch system SLURM is in use, a
batch-script similar to

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/bash 
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vspringel@mpa-garching.mpg.de
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --job-name SB128

echo
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo

mpiexec -np  $SLURM_NPROCS  ./Gadget4 param.txt 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

could be used to start a job with 2 nodes with 40 MPI ranks each (i.e.
using 80 cores in total). Here a semi-automatic setting of the
appropriate value of the total number of cores through SLURM
environment variables has been used. Note that in this example,
GADGET-4 will only use 78 MPI processes for the computational work,
because on each node one process will be set aside by the code purely
for communication purposes.


Interrupting a run                                     {#interrupting}
==================


Often, a single submission to a queue system will not provide enough
CPU-time to finish the computation of a large production
simulation. Therefore, a running simulation can be interrupted after
every timestep, and resumed at the very same place later on. If the
CPU-time limit is specified correctly in the parameter file, the code
will interrupt itself automatically before the CPU-time limit is
reached, and write a set of restart-files. Actually, each processor
writes its own restart file. These restart-files can be used to resume
the simulation conveniently (see below) at the point where it was
interrupted.

Sometimes you might want to interrupt the code manually with the
possibility to continue it later on without losing any of the
calculations. This can be achieved by generating a file named `stop`
in the output-directory of a simulation, e.g.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo > /u/vrs/myoutput/stop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code will then write a restart-file and terminate itself after the
current timestep has been completed, and the file `stop` will be
erased automatically. If a simulation reaches the specified CPU-time
limit, it creates a file `cont` in the output directory (besides the
restart files). This can be used in a job batch script to detect that
a run has ended gracefully and is hence ready for continuation, at
which point the script may submit another job for continuation.
Modern queuing systems also allow job dependencies, which can be used
to submit a chain of jobs that are executed sequentially, but only of
the previous job finishes with an OK status. In case a run crashes for
some reason, the file `cont` is absent and hence a restart should not
be initiated until the cause for the problem has been investigated and
removed. Note that restart files are also generated when the last
timestep of a simulation has been completed and the final time has
been reached. They can be used if one later wants to extend the
simulated timespan beyond the original final time.

One can also instruct the code to write restart files in regular
intervals (see parameterfile description) without stopping the
simulation. This is meant to protect against system failures or code
problems in long-running simulations. If such a thing happens, one can
then resume the run from the most recent set of restart files, which
limits the amount of lost time. Furthermore, whenever a new set of
restart files is written, the old set is first renamed into a set of
backup restart files. Hence, if a system failure occurs while the new
restart files are written, the old set of files still exits and is
valid, and can hence be used for resuming the run. To this end, one
needs to rename the backup versions of the restart files and give them
the standard name again, followed by resuming the run with the restart
option (see below).


Restarting a run                                           {#restarting}
================


Restarting from restart-files
-----------------------------

To resume a run from restart-files, start the code with an optional
flag 1 after the name of the parameterfile, in the form

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpirun -np 32  ./Gadget4  param.txt  1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


This will continue the simulation with the set of restart files in the
output-directory, with the latter being specified in the
parameterfile.  Restarting in this fashion is transparent to the code,
i.e. the simulation behaves after the restart exactly as if it had not
been interrupted to begin with. Strictly speaking this is only true if
the `PRESERVE_SHMEM_BINARY_INVARIANCE` option is activated. This
prevents that sums of partial forces may be computed for certain
particles in different order when a run is repeated due to varying
machine weather. While mathematically equivalent, this will introduce
differences in floating point round-off that can be quickly amplified
in regimes of non-linear evolution.

When the code is started with the restart-flag, the parameterfile is
parsed, but only some of the parameters are allowed to be changed,
while any changes in the others will be ignored. Which parameters
these are is explained in the section about the parameterfile, and if
such a new value for a parameter is accepted this is also reflected in
the log-files of the run.

It is **important** to not forget the 1 if you want to resume a run --
otherwise the simulation will restart from scratch, i.e. by reading in
the initial conditions again! Also note that upon restarting from
restart-files, the number of MPI ranks used for a simulation cannot be
changed; if this is attempted an error message is issued. Note that
restart files can not necessarily be transferred from one computer
system to another one, or reused when the compiler is changed, because
the layout and padding of structures stored in the binary restart
files may then be different. Hence, if you want to continue a
simulation with a different number of MPI ranks, on another computer
architecture, or with a different compiler, restarting from a snapshot
file is the method of choice. This will be discussed next.


Restarting from snapshot-files
------------------------------

There are two possibilities to restart a simulation from a previous
snapshot file. In the first possibility, one simply adopts the
snapshot file as new initial conditions. Note that this option
requires changes in the parameterfile: You need to specify the
snapshot-file as initial-conditions-file, you also need to set
`TimeBegin` (see below) to the correct time corresponding to the
snapshot-file. In addition, you should change the base filename for
the snapshot files, since the counter for the outputs will start at 0
again, thereby possibly overwriting outputs you might already have
obtained. Once this is done, you can continue/restart the run from the
snapshot-file (without the optional 1).

An alternative to this somewhat contrived procedure is to restart the
code with a flag equal to `2` after the name of the parameterfile,
i.e. just like above for the restart from restart-files, but with the
`1` replaced by `2` , and with an additional further parameter that
specifies the number of the snapshot in the output file that you want
to start from. The parameterfile normally needs not to be changed in
this case, in particular, the code takes the new starting time from
the snapshot file that is read. In addition, this also will cause any
further snapshots that are generated to have higher sequence numbers
than this starting number, i.e. the code will number all further
snapshots starting with one plus the number of the snapshot that is
used as input file. For example, the command for restarting from
snapshot 7 would look as follows:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpirun -np 32  ./Gadget4  param.txt  2 7
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Note that the restart-from-snapshot-files option allows a change of
the number of processors used for the simulation. This is not possible
if you restart from restart-files. However, restarting from
restart-files is always the preferred method to resume a simulation,
because it is faster (various start-up calculations do not have to be
done), and it minimises perturbations in the time integration scheme
of a running simulation.


Starting postprocessing                              {#postproc}
=======================

It is also possible to apply certain postprocessing algorithms built
into GADGET-4 to a snapshot file residing in the output
directory. This works similarly to the restart from a snapshot file
option in that one specifies both a start-up option and a snapshot
number (as well as potentially further parameters). For example, to
calculate a FOF/Subfind group catalogue for a snapshot file with the
number 7, one could start the code as

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpirun -np 32  ./Gadget4  param.txt  3  7
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


because the `3` selects the group finding algorithm. If such a
postprocessing option is selected, the code will not carry out any
time evolution. Instead, the given snapshot is just read in, some
postprocessing is done, the output is written to the output directory,
and then the code ends. One can obtain a short list of the available
postprocessing options when the code is started without any
paramater. Other postprocessing options for example include the
calculation of a matter power spectrum, the conversion of a snapshot
file from one of the supported file formats to another, or the
creation of cosmological initial conditions (the latter is then
actually more of a preprocessing and not a postprocessing option).

The available postprocessing options are as follows:

RestartFlag  | Action
------------ | -----------------------------------------------------------------
      0      | Read initial conditions and start simulation
      1      | Read restart files and resume simulation
      2      | Restart from specified snapshot dump and resume simulation
      3      | Run FOF and optionally SUBFIND
      4      | Calculate a matter power spectrum for specified snapshot number
      5      | Convert snapshot file to different format
      6      | Create cosmological initial conditions
      7      | Create descendant/progenitor information when merger tree is done in postprocessing
      8      | Arrange halos in merger trees, using all group catalogues up to given snapshot number
      9      | Carry out an I/O bandwidth test to determine best setting for the number of concurrent reads/writes
      10     | Rearrange particle-lightcone data in merger tree order
      11     | Rearrange most-bound snapshot data in merger tree order

