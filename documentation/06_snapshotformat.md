
Snapshot file format
====================

[TOC]

The primary result of a simulation with GADGET-4 are snapshots, which
are simply dumps of the state of the system at certain times. GADGET-4
supports parallel output by distributing a snapshot into several
files, each written by a group of processors. This procedure allows an
easier handling of very large simulations; instead of having to deal
with one file of size of a dozens of GB, say, it is much easier to
have several files with a smaller size of a few hundred MB to a couple
of GB instead. Also, the time spent for I/O in the simulation code can
be reduced if several files are written in parallel.

Each particle dump consists of `k` files, where `k` is given by
`NumFilesPerSnapshot`. For `k > 1`, the filenames are of the form
`snapshot_XXX.Y`, where `XXX` stands for the number of the dump, `Y`
for the number of the file within the dump. Say we work on dump 7 with
`k=16` files, then the filenames are `snapshot_007.0` to
`snapshot_007.15`, and if HDF5 is in use (recommended!), they will be
`snapshot_007.0.hdf5` to `snapshot_007.15.hdf5`. They will actually be
stored in a directory called `snapdir_007`, created in the output
directory. For `k=1`, the filenames will just have the form
`snapshot_XXX` and no separate subdirectory is created for the
snapshot. The base name "snapshot" can be changed by setting
`SnapshotFileBase` to another value.

Each of the individual files of a given set of snapshot files contains
a variable number of particles, but the files all have the same basic
format (this is the case for all three fileformats supported by
GADGET), and all of them are in binary. A binary representation of the
particle data is the preferred choice, because it allows much faster
I/O than ASCII files, and in addition, the resulting files are much
smaller while still providing loss-less storage of the data.


Legacy Format 1                                      {#format1}
===============

In the original default file format of GADGET (selected with
`SnapFormat=1`), the data is organised in blocks, each containing a
certain information about the particles. For example, there is a block
for the coordinates, and one for the temperatures, etc. The sequence
of blocks in snapshots files of GADGET-4 (which for the most part is
compatible with older versions of GADGET) is given in the table
below. Not all blocks are necessarily present in all simulations, for
example, the blocks describing gas properties such as internal energy
or density will only be included in hydrodynamic simulations, and some
output blocks are an optional `Config.sh` option. The presence of the
mass block depends on whether or not particle masses are defined to be
constant for certain particle types by means of the `Massarr` table in
the file header. If a non-zero mass is defined there for a certain
particle type, particles of this type will not be listed with
individual masses in the mass block. If such fixed particle masses are
defined for all types that are present in the snapshot file, the mass
block will be completely absent.



Nr  | Format2-ID | HDF5 Identifier  | Block contents
----|------------|------------------|------------------------------------------------------
1   | HEAD       | Header           | File header
2   | POS        | Coordinates      | Particle positions
3   | VEL        | Velocities       | Particle velocities
4   | ID         | ParticleIDs      | Particle IDs
5   | MASS       | Masses           | Masses (only for particle types with variable masses)
6   | U          | InternalEnergy   | Thermal energy per unit mass (only SPH particles)
7   | RHO        | Density          | Density of SPH particles
8   | HSML       | SmoothingLength  | SPH smoothing length h
9   | POT        | Potential        | Gravitational potential of particles
10  | ACCE       | Acceleration     | Acceleration of particles
11  | ENDT       | RateOfChangeOfEn | Rate of change of entropic function of SPH particles
12  | TSTP       | TimeStep         | Timestep of particles



Within each block, the particles are ordered according to their
particle type, i.e. gas particles will come first (type 0), then
type-1 particles, followed by type-2 particles, and so on. However, it
is important to realize that the detailed sequence of particles within
the blocks may change from snapshot to snapshot. Also, a given
particle may not always be stored in the snapshot file with the same
sub-number among the files belonging to one set of snapshot files.
This is because particles may move around from one processor to
another during the course of a parallel simulation. In order to trace
a particle between different outputs, one therefore has to resort to
the particle IDs, which are intended to be used to label particles
uniquely. (In fact, the uniqueness of the IDs in the initial
conditions is checked upon start up of the code.)


The first block (number 1) has a special role, it is a header which
contains global information about the particle set, for example the
number of particles of each type, the number of files used for this
snapshot set, etc.  The fields of the file header for formats 1 and 2
in GADGET-4 is given in the table below:


Header Field | Type   | HDF5 name           | Comment
-------------|--------|---------------------|-------------------
Npart[6]     | uint   | NumPart_ThisFile    | The number of particles of each type in the present file.
Nall[6]      | uint64 | NumPart_Total       | Total number of particles of each type in the simulation.
Massarr[6]   | double | MassTable           | The mass of each particle type. If set to 0 for a type which is present, individual particle masses are stored for this type.
Time         | double | Time                | Time of output, or expansion factor for cosmological simulations.
Redshift     | double | Redshift            | z=1/a-1 (only set for cosmological integrations)
BoxSize      | double | BoxSize             | Gives the box size if periodic boundary conditions are used.
NumFiles     | int    | NumFilesPerSnapshot | Number of files in each snapshot.


In GADGET-1/2/3, the outline of the fields in the header has been
slightly different, in particular the number of particle types was
always fixed to 6 (now this is given by `NTYPES`), and the total
particle number was stored as a 32-bit integer, such that simulations
exceeding particle numbers of a couple of billion needed as special
(and ugly) extension of the header, where the high-order bits in the
particle numbers where stored in a separate entry. In addition,
cosmological parameters like `Omega0` and various flags informing
about enabled/disabled code features were redundantly stored there as
well. Finally, the total header length was filled to a total length of
exactly 256 bytes, with a view to reserve the extra space for future
extensions. However, as the block-size guards (see below) anyhow store
the length of the header in the file format, this has been a
superfluous restriction. In GADGET-4, this is hence lifted, and also
other clean-ups of the header structure are implemented (such as going
to 64-bit integers for the particle numbers). While this makes the new
file format incompatible with older versions of GADGET, backwards
compatibility can be enforced by setting the `GADGET2_HEADER` switch.

To allow an easy access of the data also in Fortran (this should not
be misunderstood as an encouragement to use Fortran), the blocks are
stored using the "unformatted binary" convention of most Fortran
implementations. In it, the data of each read or write statement is
bracketed by block-size fields, which give the length of the data
block in bytes. These block fields (which are two 4-byte integers, one
stored before the data block and one after) can be useful also in
C-code to check the consistency of the file structure, to make sure
that one reads at the right place, and to skip individual blocks
quickly by using the length information of the block size fields to
fast forward in the file without actually having to read the data. The
latter makes it possible to efficiently read only certain blocks, for
example, just temperatures and densities but no coordinates and
velocities. Note however that this block size structure imposes the
restriction that individual blocks may not be larger than 4 GB in file
formats 1 and 2.

Assuming that variables have been allocated/declared appropriately, a
possible read-statement of some of these blocks in Fortran could then
for example take the form:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read (1) npart, nall, massarr, a, redshift
read (1) pos
read (1)
read (1) id
read (1) masses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, the block containing the velocities, and several
fields in the header, would be skipped. Further read-statements may
follow to read additional blocks if present. In the file
`read_snapshot.f` included in the code distribution, you can find a
simple example of a more complete read-in routine.  There you will
also find a few lines that automatically generate the filenames, and
further hints how to read in only certain parts of the data.

When you use C, the block-size fields need to be read or skipped
explicitly, which is quite simple to do. This is demonstrated in the
file `read_snapshot.c`.  Note that you need to generate the
appropriate block field when you write snapshots in C, for example
when generating initial conditions for GADGET.

GADGET-4 will check explicitly that the block-size fields contain the
correct values and refuse to continue if not. It will also use the
block-size information to convert between single and double precision
(or vice versa), and between 32-bit and 64-bit, where appropriate,
i.e. when reading in files the code will autodetect if single or
double precision is used and do conversions if needed.


Legacy Format 2                                       {#format2}
===============

Whether or not a certain block is present in GADGET's snapshot file
format depends on the type of simulation, and the makefile
options. This makes it difficult to write reasonably general I/O
routines for analysis software when file formats 1 and 2 are used,
capable of dealing with output produced for a variety of makefile
settings. For example, if one wants to analyse the (optional) rate of
entropy production of SPH particles, then the location of the block in
the file changes if, for example, the output of the gravitational
potential is enabled or not. This problem becomes particularly acute
if more complicated baryonic physics modules or optional output fields
are added to GADGET runs.

To remedy this problem, a variant of the default fileformat was
implemented in GADGET (based on a suggestion by Klaus Dolag), which
can be selected by setting `SnapFormat=2`. Its only difference is that
each block is preceded by a small additional block which contains an
identifier in the form of a 4-character string (filled with spaces
where appropriate). This identifier is listed in the "Format2-ID" in
the above block table. Using it, one can write more flexible I/O
routines and analysis software that quickly fast-forwards to blocks of
interest in a simple way, without having to know where it is expected
in the snapshot file.


HDF5 file format                                    {#format3}
================

A yet more general solution is provided by HDF5, which is nowadays the
*strongly* recommended format for the code. If the hierarchical data
format is selected as format for snapshot files or initial conditions
files, GADGET-4 accesses files with low level HDF5
routines. Advantages of HDF5 lie in its portability (e.g. automatic
endianness conversion) and generality (which however results in
somewhat more complicated read and write statements), and in the
availability of a number of tools to manipulate or display HDF5
files. Also, attributes such as units or conversion factors can be
stored along-side data-sets. Together with HDF5 utility commands like
h5ls, h5dump, etc., this makes the file format self-documenting to a
significant degree. A wealth of further information about HDF5 can be
found on the website of the HDF5 library.

In the HDF5 file format of GADGET, the blocks are stored as datasets
with names as given in the column "HDF5 Identifiers" in the table
above. To make accessing different particle types easier, for each
particle type presented in the snapshot file, a separate data group is
defined, named `ParticleType0`, `ParticleType1`, and so on, and 
the blocks of the standard file format are appearing as datasets
within these groups. These data groups can be thought of as
subdirectories in the HDF5 file, with each being devoted to one
particular particle type. In addition, in HDF5 each of these datasets
is also equipped with attribute values that store conversion factors
to cgs units, as well as factors that inform whether one has to
multiply with a certain power of the scale factor and/or of the
dimensionless Hubble parameter `h` to get to truly physical
quantities.

Finally, the file header is also moved to a special `Header` group in
the HDF5 output. This group does not contain datasets, but only
attributes that hold the values for all the fields defined for file
format 1 and 2 of GADGET-4. The names of these attributes are listed
in the column "HDF5 name" in the table above.

To graphically explore the contents of HDF5 files, the program
HDFView, available for free from the HDF-Group, is one good
possibility. It allows an easy exploration of the structure and the
contents of the various HDF5 outputs produced by GADGET4. It also
readily reveals the type of data set fields, the values of attributes,
etc., and hence greatly facilitates the writing of corresponding
analysis scripts.

Alternatively, one can also explore the contents of HDF5 files at the
command line, using simple tools that are part of the HDF5
library. For example, to list contents of the header of a snapshot
file on the command line, use the following command:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h5dump -g Header  snapshot_007.3.hdf5
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To list, for example, which blocks are available for type-0 particles,
you can use the following command:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h5ls  snapshot_007.3.hdf5/ParticleType0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Just using `h5ls snapshot_007.3.hdf5` would tell you which particle types are
available as separate data groups.

Finally, GADGET-4 also stores all the parameters and config file
options with which it was run to produce the given snapshot in all its
HDF5 files. This is realized through Parameter attributes stored in
two special data groups called `Parameters` and `Config`.

To examine the parameter file settings that were used, you can hence
issue the command:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h5dump -g Parameters  snapshot_007.3.hdf5
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

And for retrieving the `Config.sh` options, you can use:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h5dump -g Config  snapshot_007.3.hdf5
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Here are a couple of additions noted on the format, units, and
variable types of the individual blocks in the snapshot files:

- Particle positions: A floating point 3-vector for each particle,
  giving the comoving coordinates in internal length units
  (corresponds to kpc/h if the above choice for the system of units is
  adopted). N=sum(header.Npart) is the total number of particles in
  the file. Note: For file formats 1 and 2, the particles are ordered
  in the file according to their type, so it is guaranteed that
  particles of type 0 come before those of type 1, and so on. However,
  within a type, the sequence of particles can change from dump to
  dump. For file format 3, different particle types appear in separate
  particle groups.

- Particle velocities: Particle velocities `u` are in internal
  velocity units (corresponds to km/sec if the default choice for the
  system of units is adopted). For cosmological simulations, peculiar
  velocities `v` are obtained from `u` by multiplying `u` with
  sqrt(a), i.e. v = u * sqrt(a).

- Particle identifiers: These are unsigned 32-bit or 64-bit integers
  (if `IDS_64BIT` is set) and intended to provide a unique
  identification of particles.  The order of particles may change
  between different output dumps, but the IDs can be easily used to
  bring the particles back into the original order for every dump.

- Variable particle masses: Single or double precision floats, with a
  total length Nm. Only stored for those particle types that have
  variable particle masses (indicated by zero entries in the
  corresponding `massarr` entry in the header). Nm is thus the sum of
  those `npart` entries that have vanishing `massarr` . If `Nm=0`,
  this block is not present at all.

- Internal energy: Single/double precision float values for
  Ngas. Internal energy per unit mass for the `Ngas=npart(0)` gas
  particles in the file. The block is only present for `Ngas>0`. Units
  are again in internal code units, i.e. for the standard system of
  units, `u` is given in `(km/sec)^2` .

- Density: The comoving density of SPH particles. Units are again in
  internal code units, i.e. for the above system of units, `rho` is
  given in `10^10 Msun/h / (kpc/h)^3`.

- Smoothing length: Smoothing length of the SPH particles. Given in
  comoving coordinates in internal length units.

- Gravitational potential: This block will only be present if it is
  explicitly enabled in the makefile.

- Accelerations: Likewise, only present if it is explicitly enabled in
  the makefile.

- Rate of entropy production: The rate of change of the entropic
  function of each gas particles. For adiabatic simulations, the
  entropy can only increase due to the implemented artificial
  viscosity. This block will only be present if it is explicitly
  enabled in the makefile.

- Timesteps of particles: Individual particle timesteps of
  particles. For cosmological simulations, the values stored here are
  d ln(a), i.e. intervals of the natural logarithm of the expansion
  factor. This block will only be present if it is explicitly enabled
  in the makefile.



Format of initial conditions                      {#formatICs}
============================


The possible file formats for initial conditions are the same as those
for snapshot files, and are selected with the `ICFormat`
parameter. However, only the blocks up to and including the gas
temperature (if gas particles are included) need to be present; gas
densities, SPH smoothing lengths and all further blocks need not be
provided and are ignored.

In preparing initial conditions for simulations with gas particles,
the temperature block can be filled with zero values, in which case
the initial gas temperature is set in the parameterfile with the
`InitGasTemp` parameter.  However, even when this is done, the
temperature block must still be present.  Note that the field `Time`
in the header will be ignored when GADGET-4 is reading an initial
conditions file. Instead, you have to set the time of the start of the
simulation with the `TimeBegin` parameter in the parameterfile.

