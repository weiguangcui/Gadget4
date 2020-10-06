
Example Setups
==============

[TOC]

To get an idea of some of the different possible usages of GADGET-4,
and as a starting point for your own simulations and numerical
experiments, we include a small set of examples with the code
distribution. These can be found in the folder `examples`. Each of the
examples has its own subdirectory containing the example's
configuration and parameterfile, and sometimes auxiliary files for the
specific runs.

To compile the code for one of the examples, you can either copy its
configuration file `Config.sh` to the main code directory and run
`make` there, or you execute `make` by passing it the `DIR` variable
with a value that points to the subdirectory of the example. The
executable for the example is then created right in the subdirectory
of the corresponding problem, allowing you to run it right there. This
is the recommended approach.

For convenience, the corresponding make-commands are contained in the
file `make-examples.sh`, with one line with a suitable call of `make`
for each of the examples. Copying the corresponding line and pasting
it as a command into the main directory of the code will then build
the executable for the example. Executing the full shell script will
compile all the examples at once (provided there are no compilation
problems, of course).

Below, we describe each of the examples and also give a few
suggestions for setup-variants that could be of interest. Note that
for some of the examples, you need to obtain initial conditions that
are available as part of a separate
[IC-package](https://wwwmpa.mpa-garching.mpg.de/gadget4/example_ics.tar)
(183 MB) on the GADGET-4 web-site. This also contains the ICs for the
examples adopted from the GADGET-2 code distribution (these ICs are
identical to those distributed with GADGET-2).


Cosmological DM-only simulation with IC creation         {#dm20box}
================================================

The setup in `DM-L50-N128` simulates a small box of comoving
side-length 50 Mpc/h using 128^3 dark matter particles. The initial
conditions are created on the fly upon start-up of the code, using
second order Lagrangian perturbation theory with a starting redshift
of z=63. The `LEAN` option and 32-bit arithmetic are enabled to
minimize memory consumption of the code.

Gravity is computed with the TreePM algorithm at expansion order p=3.
Three output times are defined, for which FOF group finding is
enabled, and power spectra are computed as well for the snapshots that
are produced. Also, the code is asked to compute a power spectrum for
each output.


Aquarius Milky-Way Zoom                                  {#aquarius}
=======================

The directory `DM-Zoom-Aq-C-5` contains a setup for a cosmological
zoom-simulation that follows the formation of a Milky Way-sized dark
matter halo. To carry out this example, you need access to the
`Aq-C-5-dm` initial conditions files, which stem from the
[Aquarius Project](http://adsabs.harvard.edu/abs/2008MNRAS.391.1685S)
and were also used as part of the
[Aquila Project](https://ui.adsabs.harvard.edu/abs/2012MNRAS.423.1726S). 
These IC files are available as part of the 
[IC-package](https://wwwmpa.mpa-garching.mpg.de/gadget4/example_ics.tar)
on the GADGET-4 web-site.

The example uses particle type 1 for the high-resolution dark matter
particles, while types 2 to 4 are employed for more massive boundary
particles of increasingly higher mass. For computing the gravitational
forces, the TreePM algorithm is used.

The code is instructed to create 16 dumps, with output times that are
specified through an input file. Before each snapshot is written, the
FOF group finding algorithm is run on the high-resolution particles,
followed by SUBFIND, and the snapshot dumps will be stored according
to group order.

Interesting variants of this setup:

- The default setup of this example does not enable a high resolution
  mesh, which could however be used in principle.

- If the resolution is as high and spartially concentrated as here,
  `HIERARCHICAL_GRAVITY` can also be of interest for these types of
  zoom simulations.

- If a merger tree is desired, simulations of this kind could be run
  with the `MERGERTREE` option.


Colliding galaxies with star formation                    {#dm20box}
======================================

This simulation with setup in the folder `CollidingGalaxiesSFR`
considers the collision of two compound galaxies made up of a dark
matter halo, a stellar disk and bulge, and cold gas in the disk that
undergoes star formation.

Radiative cooling due to helium and hydrogen is included. Star
formation and feedback is modelled with a simple subgrid
treatment. The simulation corresponds closely to the model of a
galaxy collisions considered in the code paper.


Santa Barbara cluster                                 {#santabarbara}
=====================

The Santa Barbara cluster is a hydrodynamical simulation of the
formation of a galaxy cluster, which was introduced originally in the
code comparison paper by
[Frenk et al. (1999)](https://ui.adsabs.harvard.edu/abs/1999ApJ...525..554F).

In this example, we consider the problem at 2 x 64^3 resolution, using
the pressure-based formulation of SPH (for the sake of change, and not
because we think this is necessarily to be recommended for this
problem).


Old examples from GADGET-2 {#gadget2examples}
==========================

Galaxy collision
----------------

This purely collisionless simulation in `G2-galaxy` runs two disk
galaxies into each other, leading to a merger between the galaxies.
Each galaxy consists of a stellar disk, and a massive and extended
dark matter halo. This example uses plain Newtonian physics, with
20000 disk and 40000 halo particles in total.

The setup provided corresponds to a traditional tree algorithm with
ordinary timestepping. Alternatives to this now possible with GADGET-4
include the use the FMM algorithm for gravity, and higher order
multipole expansion. Note that this example from GADGET-2 is almost
trivially small, and some of these alternatives will only show their
real advantages in much larger simulations of higher resolution.

To get a first idea whether the example has worked you may check for
energy conservation by analysing the log-file `energy.txt`. A simple
example for doing this is provided in the form of the IDL script file
`plot_energy.pro`.


Adiabatic collapse of a gas sphere
----------------------------------

This simulation in `G2-gassphere` considers the gravitational collapse
of a self-gravitating sphere of gas which initially has a 1/r density
profile and a very low temperature. The gas falls under its own weight
to the centre, where it bounces back and a strong shock wave that
moves outwards develops. This common test problem of SPH codes has
first been described by Gus Evrard.

The simulation uses Newtonian physics in a natural system of units
(G=1). The setup corresponds to vanilla density-based SPH with the
entropy formulation introduced by
[Springel & Hernquist (2002)](http://adsabs.harvard.edu/abs/2002MNRAS.333..649S).
You can use the IDL-script `plot_energy_gassphere.pro` to display the
evolution of thermal, kinetic and potential energy for the collapsing
gassphere. Note that this is a really tiny simulation of just 1472
particles.



Cosmological formation of a cluster of galaxies
-----------------------------------------------

This problem in `G2-cluster` uses collisionless dynamics in an
expanding universe. It is a small cluster simulation that has been
set-up (a long time ago) with Bepi Tormen's initial conditions
generator ZIC using vacuum boundaries and a multi-mass technique. The
simulation has a total of 276498 particles. In a central
high-resolution zone there are 140005 particles, surrounded by a
boundary region with two layers of different softening, the inner one
containing 39616 particles, and the outer one 96877 particles.

Note that while this simulation is a cosmological simulation in
comoving coordinates, it is unusual in that it doesn't use periodic
boundary conditions but rather follows a sphere of matter around the
origin with average density equal to the mean density. This technique
is not very commonly used any more.



Large-scale structure formation including gas
---------------------------------------------

This problem in `G2-lcdm-gas` consists of 32^3 dark matter, and 32^3
gas particles, following structure formation in a periodic box of 50
Mpc/h on a side in a LCDM universe.  Only adiabatic gas physics is
included, and the minimum temperature of the gas is set to 1000
K. This simple example uses grid initial conditions, where gas
particles are put at the centres of the grid outlined by the dark
matter particles. The simulation starts at z=10, and the code will
produce snapshot files at redshifts 5, 3, 2, 1, and 0.  As in the
other old GADGET-2 examples, the SPH setup is density-based SPH based
on the entropy formation.
