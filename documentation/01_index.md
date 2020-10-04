
Introduction to GADGET-4                            {#mainpage}
========================

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ___    __    ____    ___  ____  ____       __
 / __)  /__\  (  _ \  / __)( ___)(_  _)___  /. |
( (_-. /(__)\  )(_) )( (_-. )__)   )( (___)(_  _)
 \___/(__)(__)(____/  \___/(____) (__)       (_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


GADGET-4 is a massively parallel code for N-body/hydrodynamical
cosmological simulations. It is a flexible code that can be applied to
a variety of different types of simulations, offering a number of
sophisticated simulation algorithms.  An account of the numerical
algorithms employed by the code is given in the original code paper,
subsequent publications, and this documentation.

GADGET-4 was written mainly by
[Volker Springel](mailto:vspringel@mpa-garching.mpg.de), with
important contributions and suggestions being made by numerous people,
including [Ruediger Pakmor](mailto:rpakmor@mpa-garching.mpg.de),
[Oliver Zier](mailto:ozier@mpa-garching.mpg.de), and
[Martin Reinecke](mailto:martin@mpa-garching.mpg.de).

Except in very simple test problems, one or more specialized code
options will typically be used. These are broadly termed modules and
are activated through compile-time flags. Each module may have several
optional compile-time options and/or parameter file values.


Overview and history                                    {#overview}
====================

In its current implementation, the simulation code GADGET-4 (**GA**
laxies with **D** ark matter and **G** as int **E** rac **T** - this
peculiar acronym hints at the code's origin as a tool for studying
galaxy collisions) supports collisionless simulations and smoothed
particle hydrodynamics on massively parallel computers.  All
communication between concurrent execution processes is done either
explicitly by means of the message passing interface (MPI), or
implicitly through shared-memory accesses on processes on multi-core
nodes. The code is mostly written in ISO C++ (assuming the C++11
standard), and should run on all parallel platforms that support at
least MPI-3.  So far, the compatibility of the code with current
Linux/UNIX-based platforms has been confirmed on a large number of
systems.

The code can be used for plain Newtonian dynamics, or for cosmological
integrations in arbitrary cosmologies, both with or without periodic
boundary conditions. Stretched periodic boxes, and special cases such
as simulations with two periodic dimensions and one non-periodic
dimension are supported as well.  The modeling of hydrodynamics is
optional. The code is adaptive both in space and in time, and its
Lagrangian character makes it particularly suitable for simulations of
cosmic structure formation. Several post-processing options such as
group- and substructure finding, or power spectrum estimation are built
in and can be carried out on the fly or applied to existing
snapshots. Through a built-in cosmological initial conditions
generator, it is also particularly easy to carry out cosmological
simulations. In addition, merger trees can be determined directly by
the code.

The main reference for numerical and algorithmic aspects of the code
is the paper "Simulating cosmic structure formation with the GADGET-4
code" (Springel et al., 2020, MNRAS, submitted), and references
therein. Further information on the previous public versions of GADGET
can be found in "The cosmological simulation code GADGET-2" (Springel,
2005, MNRAS, 364, 1105), and in "GADGET: A code for collisionless and
gas-dynamical cosmological simulations" (Springel, Yoshida & White,
2001, New Astronomy, 6, 51). It is recommended to read these papers
before attempting to use the code. This documentation provides
additional technical information about the code, hopefully in a
sufficiently self-contained fashion to allow anyone interested to
learn using the code for cosmological N-body/SPH simulations on
parallel machines.

Most core algorithms in GADGET-4 have been written by Volker Springel
and constitute evolved and improved versions of earlier
implementations in GADGET-2 and GADGET-3. Substantial contributions to
the code have also been made by all the authors of the GADGET-4 code
paper. Note that the code is made publicly available under the GNU
general public license. This implies that you may freely copy,
distribute, or modify the sources, but the copyright for the original
code remains with the authors. If you find the code useful for your
scientific work, we kindly ask you to include a reference to the code
paper on GADGET-4 in all studies that use simulations carried out with
the code.


Disclaimer                                             {#disclaimer}
==========

It is important to note that the performance and accuracy of the code
is a sensitive function of some of the code parameters. We also stress
that GADGET-4 comes without any warranty, and without any guarantee
that it produces correct results. If in doubt about something, reading
(and potentially improving) the source code is always the best
strategy to understand what is going on!

**Please also note the following:**

The numerical parameter values used in the examples contained in the
code distribution do not represent a specific recommendation by the
authors! In particular, we do not endorse these parameter settings in
any way as standard values, nor do we claim that they will provide
fully converged results for the example problems, or for other initial
conditions. We think that it is extremely difficult to make general
statements about what accuracy is sufficient for certain scientific
goals, especially when one desires to achieve it with the smallest
possible computational effort. For this reason we refrain from making
such recommendations. We encourage every simulator to find out for
herself/himself what integration settings are needed to achieve
sufficient accuracy for the system under study. We strongly recommend
to make convergence and resolution studies to establish the range of
validity and the uncertainty of any numerical result obtained with
GADGET-4.

