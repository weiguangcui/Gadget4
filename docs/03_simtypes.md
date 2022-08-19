
Types of simulations
====================

[TOC]

There are a number of qualitatively different simulation set-ups that
can be run with GADGET-4, which differ mostly in the type of boundary
conditions employed, in the physics that is included, whether or not
cosmological integrations with comoving coordinates are used, and in
the selection of numerical algorithms. A schematic overview of a
subset of these simulation types, concentrating on how gravity
calculations are done, is given in the following table.


<table>
<caption id="multi_row">Overview of simulation types</caption>
<tr><th>Type of Simulation  <th>Computational Method     <th>Remarks

<tr><td> <img src="../../documentation/img/type1.png" width="100">  Newtonian space
    <td> Gravity: Tree, SPH (optional), vacuum boundary conditions
	<td> OmegaLambda should be set to zero

<tr><td> <img src="../../documentation/img/type2.png" width="100">  Periodic long box 
    <td> No gravity, only SPH, periodic boundary conditions
	<td> SELFGRAVITY needs to be deactivated,  LONG_X/Y/Z may be
 set to scale the dimensions of the box 

<tr><td> <img src="../../documentation/img/type3.png" width="100">  Cosmological, physical coordinates
    <td> Gravity: Tree, SPH, vacuum boundaries
	<td> ComovingIntegrationOn set to zero

<tr><td> <img src="../../documentation/img/type4.png" width="100">  Cosmological, comoving coordinates
    <td> Gravity: Tree, SPH, vacuum boundaries
	<td> ComovingIntegrationOn set to one

<tr><td> <img src="../../documentation/img/type5.png" width="100">  Cosmological, comoving periodic box
    <td> Gravity: Tree with Ewald-correction, SPH, periodic boundaries
	<td> PERIODIC needs to be set

<tr><td> <img src="../../documentation/img/type6.png" width="100">  Cosmological, comoving coordinates, TreePM
    <td> Gravity: Tree with long range PM, SPH, vacuum boundaries
	<td> PMGRID needs to be set

<tr><td> <img src="../../documentation/img/type7.png" width="100">  Cosmological, comoving periodic box, TreePM
    <td> Gravity: Tree with long range PM, SPH, periodic boundaries
	<td> PERIODIC and PMGRID need to be set

<tr><td> <img src="../../documentation/img/type8.png" width="100">  Cosmological, comoving coordinates, TreePM, Zoom
    <td> Gravity: Tree with long-range and intermediate-range PM, SPH, vacuum boundaries
	<td> PMGRID and PLACEHIGHRESREGION need to be set

<tr><td> <img src="../../documentation/img/type9.png" width="100">  Cosmological, periodic comoving box, TreePM, Zoom
    <td> Gravity: Tree with long-range and intermediate-range PM, SPH, periodic boundaries
	<td> PERIODIC, PMGRID and PLACEHIGHRESREGION need to be set

<tr><td> <img src="../../documentation/img/type10.png" width="100"> Newtonian space, TreePM
    <td> Gravity: Tree with long-range PM, SPH, vacuum boundaries
	<td> PMGRID needs to be set

</table>



Cosmological Simulations                     {#cosmosims}
========================

Cosmological integrations with comoving coordinates are selected via
the parameter `ComovingIntegrationOn` in the parameterfile. Otherwise
the simulations always assume ordinary Newtonian space. The use of
self-gravity needs to be explicitly enabled through the `SELFGRAVITY`
compile-time switch, otherwise only external gravitational fields
and/or hydrodynamical processes are computed.  Periodic boundary
conditions and the various variants of the Tree, FMM, TreePM, and
FMM-PM algorithms require compile-time switches to be set
appropriately in the configuration file.

In particular, the TreePM algorithm is switched on by passing the
desired mesh-size at compile time via the Config.sh file to the
code. The relevant parameter is `PMGRID`, see below. Using an explicit
force split, the long-range force is then (normally) computed with
Fourier techniques, while the short-range force is calculated with the
tree. Because the tree needs only be walked locally, a speed-up can
arise, particularly for near to homogeneous particle distributions,
but not only restricted to them. Both periodic and non-periodic
boundary conditions are implemented for the TreePM and FMM-PM
approaches. In the non-periodic case, the code will internally compute
FFTs of size `HRPMGRID`. In order to accommodate the required
zero-padding, only half that size is actually used to cover the
high-res region. If `HRPMGRID` is not specified, a default value equal
to `PMGRID` is used.  For zoom-simulations, the code automatically
places the refined mesh-layer on the high-resolution region. This is
controlled with the `PLACEHIGHRESREGION` option.


Newtonian space                              {#newtonsims}
===============


Periodic boxes need not necessarily be cubical in GADGET-4. They can
be, if desired, stretched independently in any of the coordinate
directions through `LONG_X`, `LONG_Y` and `LONG_Z` options. This works
both for SPH simulations and also for simulations including
self-gravity. In the latter case, one may also choose one direction to
be non-periodic in self-gravity, realizing mixed boundary conditions
that allow one to simulate infinitely extended two-dimensional
systems.


Stretched boxes                               {#stretched}
===============

This can now also be done with gravity and periodic boundaries in just
two dimensions.


Two-dimensional simulations                   {#twodims}
===========================

It is also possible to run SPH simulations in two dimensions only,
which is primarily provided for test-purposes. Note that the code is
not really optimised for doing this; three coordinates are still
stored for all particles in this case, and the computations are
formally carried out as in the 3D case, except that all particles lie
in one coordinate plane, i.e. either with equal x, y, or
z-coordinates.
