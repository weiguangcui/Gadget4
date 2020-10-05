#!/bin/bash


# Dark matter cosmological simulation with 128^3 resolution in a 50 Mpc/h box for
# which ICs are created on the fly, and the LEAN option is used.
make -j 8 DIR=examples/DM-L50-N128


# DM-only zoom simulation of the formation of a Milky Way-sized halo.
make -j 8 DIR=examples/DM-Zoom-Aq-C-5


# Simple galaxy collision with gas and ongoing star formation.
make -j 8 DIR=examples/CollidingGalaxiesSFR


# A simulation of the Santa Barbara cluster at 2 x 64^3 resolution using pressure-SPH.
make -j 8 DIR=examples/SantaBarbara-PSPH-64



# The following examples are realizations of the (small) test problems contained
# in the GADGET-2 distribution from Springel (2005).

# collisionless galaxy collision
make -j 8 DIR=examples/G2-galaxy

# cosmological simulation of a galaxy cluster
make -j 8 DIR=examples/G2-cluster 

# gravitational collapse of a cold gas cloud ("Evrard test")
make -j 8 DIR=examples/G2-gassphere 

# low resolution cosmological simulation with gas
make -j 8 DIR=examples/G2-lcdm-gas 
