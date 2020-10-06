
Guide to code changes
=====================

[TOC]

In the following, we give a set of assorted hints and recommendations
about modifying and/or extending the code. We also comment about some
differences with respect to GADGET-2/3 in terms of code usage and its
architecture.

Coding style guide                                     {#style}
==================

We strongly recommend to follow the general coding practices (even if
you don't like them) and architecture of GADGET-4 in modifying or
extending the code. Only then different extensions by different people
have a chance to happily coexist with each other. This concerns mostly
the following points:

Code extensions
---------------

- Non-standard code extensions should always be written such that they
  can be switched off if not needed, and have no side effects on
  existing code when this is done. Normally this means that they have
  to be enclosed in conditional compilation precompiler statements
  (#ifdef), especially if variables in the global structures of the
  code need to be allocated for the extension. However, if the
  extension's execution can also be controlled at run-time by simple
  variables, then consider introducing a parameterfile variable to
  control the extension. In general, the number of symbols (and this
  additional `Config.sh` options) to control conditional compilation
  should be kept to a minimum.

- Do not place any substantial piece of code belonging to your
  extension into existing functions of the code. Write your own
  classes or functions for a code extension, and only place a function
  call (if needed bracketed by an #ifdef) into the appropriate place
  of the primary code. Also, place your extension functions into
  separate source files.


General code-style principles
-----------------------------

- Code formatting: Be consistent with the code formatting of the main
  code, which is more or less GNU-style, and is obtained by running
  the code formatting tool "clang-format" of the clang compiler
  suite. The specific options used by GADGET-4 are defined in the
  hidden file ".clang-format" in the code's main directory. Running
  there something like

~~~~~~~~~~~~~~~~~
  clang-format-mp-6.0 -style=file -i src/*/*
~~~~~~~~~~~~~~~~~

  should then make sure that all your source file(s) have a consistent
  indention and code formatting.

- Name functions all in lower case as a "command" that is descriptive of what
  the function does. Different words should normally be separated by
  underscores, e.g. calculate_normal_vector_for_triangle(...) For all functions,
  input arguments should come first, output arguments last.

- Global variables (whose use should be kept to an absolute minimum,
  and the sins of GADGET-4 in this regard should not be repeated)
  start with an upper case character. Global variable names are nouns,
  with words separated by mixed lower/upper case characters, and
  contain no underscores. Use abbreviations that are clear,
  e.g. `NumForceCalculations`. If you have to use a global variable in
  your own code class or module, try to restrict its scope to the
  classes or files of your module through appropriate declarations.

- Local variables start with lowercase, and should have descriptive
  names, too, except for simple loop iterators and the like. Try to
  narrow the scope of a variable as much as possible, for example by
  declaring them inside a block where they are needed. Generally
  define them as close as possible to where they are needed for the
  first time. Declaration and initialization should be combined in one
  command where possible, e.g.
~~~~~~~~~~~~~
int n = get_particle_count();
~~~~~~~~~~~~~
  instead of
~~~~~~~~~~~~~
int n;
  ...
n = get_particle_count();
~~~~~~~~~~~~~

- Avoid repetition of code, i.e. do not use cut & paste to implement
  the same or similar functionality multiple times. This is always an
  indication that a new function is in order here. Break up long
  functions into smaller more manageable pieces.

- Preprocessor macros that have arguments should be avoided whenever
  possible, and really should only be used for a compelling
  reason. Note that in C++ you can use inlined and type-safe functions
  instead. In any case, preprocessor macros should be fully
  capitalized.

- Magic numbers and non-trivial numerical constants should be avoided
  in the actual code and instead be replaced by a symbolic constant
  declared within a header file, with a name all in
  uppercase. Consider enums instead of the use of numerical constants
  to distinguish between different cases.

- All warnings emitted by the compiler upon compilation with "-Wall"
  should be addressed. Unless there are extremely good reasons, the
  code should compile without any warnings left.

- Include consistent commenting in your code. The meaning of all
  global variables should be commented where they are introduced,
  ideally with doxygen syntax, e.g:

~~~~~~~~~~~~~
       int MyGlobalCount;   /*!< counts the number of timesteps */
~~~~~~~~~~~~~

- Functions should be preceded by a brief explanation of what the
  function does, including instruction, if any, about how the function
  may be used, e.g.:
~~~~~~~~~~~~~
      /*!
       * Insert the point P[i] into the partice list. Start
       * the search at the current particle t.
       */
       int insert_point(int i, int t)
       {
         ...
       }
~~~~~~~~~~~~~
	   
- You do not need to state here *how* the function achieves what it
  does; this can be stated if appropriate in comments in the function
  body. There, avoid superfluous comments that just reflect what's
  obvious from the code anyway, like 
~~~~~~~~~~~~~
       do_domain_decomposition();   /* call domain decomposition */
~~~~~~~~~~~~~
	   
- Instead, focus on comments that help one to quickly understand/check
  what the code tries to do at an algorithmic level. If complicated
  formulae are implemented, try to include in a comment a reference to
  the equation that is implemented.



Adding a config-option                              {#addconfigs}
======================

To add a new symbolic configuration option to the code, you need to
add it to the file `Template-Config.sh` besides just using it in the
code. In addition, you must add a short explanation of the option to
the file `documentation/04_config-options.md` following the syntax
used for the other configuration option. Only then the checking
scripts of the source code will be happy and accept the new option.
Note that for this the option really needs to be used somewhere in the
source code that is described as belonging to GADGET-4 by the
`Makefile`.


Adding parameters                                   {#addparams}
=================

To add a new parameter to the parameterfile of GADGET-4, the following steps
are necessary:

- Add the new variable to the structure `global_data_all_processes` in
  the file `data/allvars.h`.  The variable has to be of type double,
  int, or a character string.

- In the file `io/parameters.c` , add a function call to
  add_param(...) for the new parameter in the function
  `register_parameters()` . Simply follow the examples for the other
  parameters. If your parameter is optional and should only be active
  when a specific option is activated, bracket the add_param() call
  with a corresponding #ifdef/#endif

- Add a short explanation of the new parameter in the file
  `documentation/05_parameterfile.md` This is needed otherwise the
  code checking scripts will not accept the new parameter.

Once these steps are followed, the new parameter will be included in
the parsing of the parameterfile, and will also be stored in restart
files as well as in HDF5 output files.


Adding a source file                                 {#addsources}
====================

To add a new source file, create it in a subdirectory of `src/` , not
in `src/` itself. Preferably you create your own subdirectory for
this, say `src/mymodule` Then in the `Makefile` of the code, add
appropriate statements that list the file(s) and header file(s)
belonging to your new source file. For example, something like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBDIRS += mymodule
OBJS    += mymodule/mysource1.o  mymodule/mysource2.o
INCL    += mymodule/mymodule.h
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your new source file should also come with a header file, at the very
least containing a function prototype for one of your new functions
that is intended to be called from the main code. This header file can
then be included in the main code at an appropriate place to
facilitate a call to your new code. The header file should be
protected by a header file guard, i.e. you should add something like

    #ifndef MYMODULE_H
    #define MYMODULE_H

at the beginning of the header file, and

    #endif

at the end. To tell the code checks scripts that the constant
`MYMODULE_H` is not an undocumented code configuration parameter, you
also need to add it to the file `defines_extra` in the main directory
of the code. Importantly, you should also add your new source file to
version control system. Ideally, you make a git-branch of GADGET-4 and
then implement your extension/modification in this branch. Try to keep
the branch current with the main line of the code as long as it is
still kept alive as a separate branch. If possible, try to integrate
your feature into a main line of development for GADGET-4, for example
through a pull request, and if successful close your branch again. To
update the html-documentation of the code with your new configuration
and parameterfile options, and to also include your source in the
cross-referenced source code documentation, run the command `doxygen`
in the main directory of the code. (If doxygen should be unavailable
on your local machine, you can always easily install it.)

 

Migration to GADGET4                                 {#migration}
====================

Main changes in the code
------------------------

There have been many major changes in the internal workings between
previous versions of GADGET and GADGET-4, affecting nearly all parts
of the code. In fact, basically all parts of the code have been
completely rewritten since GADGET-3 (sometimes several times), or were
in any case substantially revised. These modifications were done to
improve the accuracy and the performance of the code, to reduce its
memory consumption, and to make it applicable to a wider range of
simulation types. Detailed explanations of the reasoning behind each
of these modifications is beyond the scope of this short set of
notes. Some of this information can be however found in the code paper
on GADGET-4.

Some of the most important changes in code design include:

- A new hierarchical timestepping option for gravity

- An optional FMM solver for gravity

- Availability of hybrid parallelization with shared memory use

- Support for very large simulation sizes

- Stretched boxes also for gravity, as well as mixed boundary
  conditions for gravity

- Inclusion of FOF and SUBFIND in the public code version, a new
  parallelization scheme for SUBFIND, and the addition of SUBFIND_HBT

- Inclusion of a power spectrum estimator

- An optional column-based FFT can be used

- Inclusion of an updated version of the N-GenIC initial conditions
  generator that also supports 2LPT

- Explicit vectorization support in some parts of the code

- New domain decomposition algorithm

- Support of two formulations of SPH

- On-the-fly merger tree construction

- Continuous light-cone output

- Integration of a simple cooling and star formation formulation in
  the public version of the code


Migrating old set-ups to GADGET-4
---------------------------------

If you have used earlier version of GADGET before, you will be
familiar with all the practicalities of running GADGET already. This
is because nothing much has really changed here in GADGET-4. In
particular, the mechanisms for starting, interrupting, and resuming a
simulation are still the same.

However, as described in this guide, there have been a number of
changes in the parameterfile that controls each simulation, and there
are more compilation options that are now more conveniently controlled
through the `Config.sh` file.  If you want to reuse an old
parameterfile with GADGET-4, you therefore have to make a few
modifications of it. If you simply start the new code with the old
parameterfile, you will get error messages that complain about
obsolete or missing parameters. Delete or add these parameters as
needed. Also, carefully review which makefile options are still
appropriate for a given simulation.  While many of the basic compile
time options available in previous versions of GADGET are still
around, others have become obsolete. The code will complain about
options that are set in `Config.sh` but don't exist any more, so
you should be automatically informed about this when you try to use an
old setup.  Importantly, however, the structure and format of the
snapshot files produced by the code has hardly changed (and the system
of units has not changed at all), meaning that you should still be
able to use old analysis software from previous versions of GADGET
essentially unchanged.


Notes on memory consumption of GADGET-4
=======================================


GADGET-4 uses an internal memory manager where a large block of memory
is allocated once at the beginning of the simulation, and the code
then satisfies all its internal allocations out of this block. The
size of this block corresponds to the maximum allowed memory usage per
MPI process, and it is set by the parameter file `MaxMemSize`.

The code will periodically output information to the file `memory.txt`
about all allocated memory blocks and their sizes. This happens
whenever a new high watermark in the memory allocation is reached. It
is best to check this table to get an accurate assessment of the
amount of memory the code at least needs to run successfully for a
certain particle number. For communication and particle exchanges, the
code will use buffers that automatically adjust their size to the
amount of still available memory. It is hence advantageous to set
`MaxMemSize` as large as possible for the amount of physical memory
available on the compute nodes. An overcommitment of the physical
memory of a node is prevented by the code (it then exists with an
error message), but this check only works if `HOST_MEMORY_REPORTING`
is enabled. At the beginning of the stdout-file, you can also find the
sizes of the most important code structures for the chosen
configuration. In particular this tells about how many bytes per
collisionless particle are needed, and how many bytes are in addition
needed to store one SPH particle.



Notes on various limits in GADGET-4
===================================

- Maximum number of particles per MPI-rank is restricted to 2^31 ~ 2
  Billion

- Maximum number of groups/subhalos per MPI-rank is restricted to 2^31
  ~ 2 Billion

- The total number of particles, and the total number of
  groups/subhalos, can however be much larger than this if a
  sufficiently large number of MPI ranks is used.
  
- Maximum particle number in a single group or subhalo is for the
  default setting 2^31 ~ 2 billion, but this can be enlarged if needed
  by setting the option `FOF_ALLOW_HUGE_GROUPLENGTH`
  
- The maximum number of subhalos per group is restricted nevertheless
  to 2^31 ~ 2 billion.
 
- If the legacy output formats 1/2 are used, each block in the output
  is restricted to 2 Gbyte in size or less. This also means that the
  maximum number of particles per single output file, and the maximum
  number of groups/subhalos per single group catalogue file can be at
  most of order 2 Billion, normally substantially less than
  that. While this restriction can be circumvented by splitting output
  over enough files, it provides one further motivation to use the
  HDF5 format!
  
- Maximum number of lightcones: 256
