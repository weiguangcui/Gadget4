/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file domain_box.cc
 *
 *  \brief routines for finding domain extension, random shifting, and periodic wrapping if needed
 */

#include "gadgetconfig.h"

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../pm/pm.h"
#include "../system/system.h"

using namespace std;

/*! \file domain_box.cc
 *  \brief finds extension of particle set and/or wraps them back into the fundamental periodic box
 */

/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
template <>
void domain<simparticles>::do_box_wrapping(void)
{
  if(Mode != STANDARD)
    return;

#ifdef RANDOMIZE_DOMAINCENTER
  /* remove previous shift vector */
  for(int i = 0; i < Tp->NumPart; i++)
    {
      for(int j = 0; j < 3; j++)
        Tp->P[i].IntPos[j] -= Tp->CurrentShiftVector[j];

      Tp->constrain_intpos(Tp->P[i].IntPos);
    }

  for(int j = 0; j < 3; j++)
    Tp->CurrentShiftVector[j] = 0;
#endif

#ifndef PERIODIC
  /* If we don't use periodic boundaries, check whether we lie outside the central 3/8 of the region chosen for the root node.
   * This is based on the notion of using 1/4 for the initial region, leaving 1/4 as a safe buffer. Once half of this buffer
   * is used up, we trigger an adjustment. */

  MyIntPosType leftbound  = 5 * (((MyIntPosType)1) << (BITS_FOR_POSITIONS - 4));  /* 5/16 of full box length */
  MyIntPosType rightbound = 11 * (((MyIntPosType)1) << (BITS_FOR_POSITIONS - 4)); /* 11/16 of full box length */

  int flag = 0;
  int iter = 0;

  do
    {
      flag = 0;
      for(int i = 0; i < Tp->NumPart; i++)
        for(int j = 0; j < 3; j++)
          {
            if(Tp->P[i].IntPos[j] < leftbound)
              flag = 1;

            if(Tp->P[i].IntPos[j] > rightbound)
              flag = 1;
          }

      MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, Communicator);

      if(flag)
        {
          domain_printf("DOMAIN: Simulation region has enlarged, need to adjusted mapping.\n");

          for(int i = 0; i < Tp->NumPart; i++)
            for(int j = 0; j < 3; j++)
              Tp->P[i].IntPos[j] = (Tp->P[i].IntPos[j] >> 1) + (((MyIntPosType)1) << (BITS_FOR_POSITIONS - 2));

          Tp->FacCoordToInt *= 0.5;
          Tp->FacIntToCoord *= 2.0;
          Tp->RegionLen *= 2.0;

          for(int j = 0; j < 3; j++)
            Tp->RegionCorner[j] = Tp->RegionCenter[j] - 0.5 * Tp->RegionLen;

#if defined(PMGRID) && (!defined(PERIODIC) || defined(PLACEHIGHRESREGION))
          /* make sure that we compute new kernels the next time we execute a non-periodic pm calculation */
          Tp->OldMeshSize[0] = 0;
          Tp->OldMeshSize[1] = 0;
#endif

          iter++;
        }

      if(iter > 5)
        Terminate("too many iterations");
    }
  while(flag);

#endif

#ifdef RANDOMIZE_DOMAINCENTER
  /* determine new shift vector */

#if defined(RANDOMIZE_DOMAINCENTER_TYPES) || defined(PLACEHIGHRESREGION)
  domain_find_type_extension();
#endif

  if(ThisTask == 0)
    {
#if defined(RANDOMIZE_DOMAINCENTER_TYPES) || defined(PLACEHIGHRESREGION)
      int count = 0;
#endif

      for(int j = 0; j < 3; j++)
        {
          Tp->CurrentShiftVector[j] = get_random_number() * pow(2.0, 32);
#ifdef POSITIONS_IN_64BIT
          Tp->CurrentShiftVector[j] <<= 32;
          Tp->CurrentShiftVector[j] += get_random_number() * pow(2.0, 32);
#endif
#ifdef POSITIONS_IN_128BIT
          Tp->CurrentShiftVector[j] <<= 32;
          Tp->CurrentShiftVector[j] += get_random_number() * pow(2.0, 32);
          Tp->CurrentShiftVector[j] <<= 32;
          Tp->CurrentShiftVector[j] += get_random_number() * pow(2.0, 32);
          Tp->CurrentShiftVector[j] <<= 32;
          Tp->CurrentShiftVector[j] += get_random_number() * pow(2.0, 32);
#endif

#if defined(RANDOMIZE_DOMAINCENTER_TYPES) || defined(PLACEHIGHRESREGION)
          MyIntPosType boxoff   = (Tp->CurrentShiftVector[j] & Tp->PlacingMask);
          MyIntPosType inboxoff = (Tp->CurrentShiftVector[j] - boxoff) % (Tp->PlacingBlocksize - domainInnersize);

          MyIntPosType off = domainXmintot[j] + domainReferenceIntPos[j] + Tp->CurrentShiftVector[j];
          Tp->CurrentShiftVector[j] -= off;  // now we have the high-res region aligned with the left box size
          Tp->CurrentShiftVector[j] += boxoff + inboxoff;

          if(domain_type_extension_overlap(j))
            {
              domain_printf("(Tp->PlacingBlocksize - domainInnersize)=%g  %g\n",
                            (Tp->PlacingBlocksize - domainInnersize) * Tp->FacIntToCoord, inboxoff * Tp->FacIntToCoord);
              domain_printf("DOMAIN: Need to draw shift vector again for j=%d\n", j);
              Terminate("we should not get here anymore\n");
              j--;  // causes a repeat of the loop for the same index
              count++;
              if(count > 1000)
                Terminate("too many repeats");
              continue;
            }
#endif
        }
    }

#ifdef GRAVITY_TALLBOX
  Tp->CurrentShiftVector[GRAVITY_TALLBOX] = 0;
#endif

  MPI_Bcast(Tp->CurrentShiftVector, 3 * sizeof(MyIntPosType), MPI_BYTE, 0, Communicator);

  for(int i = 0; i < Tp->NumPart; i++)
    {
      for(int j = 0; j < 3; j++)
        Tp->P[i].IntPos[j] += Tp->CurrentShiftVector[j];

      Tp->constrain_intpos(Tp->P[i].IntPos);
    }

  domain_printf("DOMAIN: New shift vector determined (%g %g %g)\n",
                ((MySignedIntPosType)Tp->CurrentShiftVector[0]) * Tp->FacIntToCoord,
                ((MySignedIntPosType)Tp->CurrentShiftVector[1]) * Tp->FacIntToCoord,
                ((MySignedIntPosType)Tp->CurrentShiftVector[2]) * Tp->FacIntToCoord);
#endif
}

#if defined(RANDOMIZE_DOMAINCENTER_TYPES) || defined(PLACEHIGHRESREGION)

template <typename partset>
int domain<partset>::domain_type_extension_overlap(int j)
{
  MyIntPosType *xmintot = (MyIntPosType *)domainXmintot;
  MyIntPosType *xmaxtot = (MyIntPosType *)domainXmaxtot;

  MyIntPosType xmin = xmintot[j] + domainReferenceIntPos[j] + Tp->CurrentShiftVector[j];
  MyIntPosType xmax = xmaxtot[j] + domainReferenceIntPos[j] + Tp->CurrentShiftVector[j];

  if((xmin & Tp->PlacingMask) != (xmax & Tp->PlacingMask))
    return 1;
  else
    return 0;
}

template <typename partset>
void domain<partset>::domain_find_type_extension(void)
{
  /* first, find a reference coordinate by selecting an arbitrary particle in the respective region. For definiteness, we choose the
   * first particle */

  int have_high_mesh = NTask; /* default is we don't have a particle */

  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(((1 << Tp->P[i].getType()) & (RANDOMIZE_DOMAINCENTER_TYPES)))
        {
          for(int j = 0; j < 3; j++)
            domainReferenceIntPos[j] = Tp->P[i].IntPos[j];

          have_high_mesh = ThisTask;
          break;
        }
    }

  int have_global[2] = {have_high_mesh, ThisTask};

  MPI_Allreduce(MPI_IN_PLACE, have_global, 1, MPI_2INT, MPI_MINLOC, Communicator);

  if(have_global[0] >= NTask)
    Terminate("have_global[0]=%d  >= NTask=%d: Don't we have any particle?  Note: RANDOMIZE_DOMAINCENTER_TYPES=%d is a bitmask",
              have_global[0], NTask, RANDOMIZE_DOMAINCENTER_TYPES);

  MPI_Bcast(domainReferenceIntPos, 3 * sizeof(MyIntPosType), MPI_BYTE, have_global[1], Communicator);

  /* find enclosing rectangle */

  MySignedIntPosType xmin[3], xmax[3];

  for(int j = 0; j < 3; j++)
    {
      xmin[j] = 0;
      xmax[j] = 0;
    }

  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(((1 << Tp->P[i].getType()) & (RANDOMIZE_DOMAINCENTER_TYPES)))
        {
          MyIntPosType diff[3] = {Tp->P[i].IntPos[0] - domainReferenceIntPos[0], Tp->P[i].IntPos[1] - domainReferenceIntPos[1],
                                  Tp->P[i].IntPos[2] - domainReferenceIntPos[2]};

          MySignedIntPosType *delta = (MySignedIntPosType *)diff;

          for(int j = 0; j < 3; j++)
            {
              if(delta[j] > xmax[j])
                xmax[j] = delta[j];
              if(delta[j] < xmin[j])
                xmin[j] = delta[j];
            }
        }
    }

  MPI_Allreduce(xmin, domainXmintot, 3, MPI_MyIntPosType, MPI_MIN_MySignedIntPosType, Communicator);
  MPI_Allreduce(xmax, domainXmaxtot, 3, MPI_MyIntPosType, MPI_MAX_MySignedIntPosType, Communicator);

  for(int j = 0; j < 3; j++)
    domainXmaxtot[j] += 1; /* so that all particles fulfill   xmin <= pos < xmax instead of xmin <= pos <= xmax*/

  domainInnersize = domainXmaxtot[0] - domainXmintot[0];

  if((MyIntPosType)(domainXmaxtot[1] - domainXmintot[1]) > domainInnersize)
    domainInnersize = domainXmaxtot[1] - domainXmintot[1];

  if((MyIntPosType)(domainXmaxtot[2] - domainXmintot[2]) > domainInnersize)
    domainInnersize = domainXmaxtot[2] - domainXmintot[2];

  domain_printf("DOMAIN: Shrink-wrap region size for RANDOMIZE_DOMAINCENTER_TYPES is %g\n", domainInnersize * Tp->FacIntToCoord);

  if(domainInnersize * Tp->FacIntToCoord >= 0.125 * All.BoxSize)
    Terminate("inappropriately big region selection for RANDOMIZE_DOMAINCENTER_TYPES");

  /* increase the region by at least 1/8 of its size to still allow some randomness in placing the particles within the high-res node
   */
  MyIntPosType ref_size = domainInnersize + (domainInnersize >> 3);

  Tp->PlacingBlocksize = 1;
  Tp->PlacingMask      = ~((MyIntPosType)0);

  for(int i = 0; i < BITS_FOR_POSITIONS; i++)
    {
      if(Tp->PlacingBlocksize >= ref_size)
        break;

      Tp->PlacingBlocksize <<= 1;
      Tp->PlacingMask <<= 1;
    }

  domain_printf("DOMAIN: We enlarge this to %g    (%g times smaller than boxsize)\n", Tp->PlacingBlocksize * Tp->FacIntToCoord,
                All.BoxSize / (Tp->PlacingBlocksize * Tp->FacIntToCoord));
}
#endif

#ifdef LIGHTCONE_PARTICLES
template <>
void domain<lcparticles>::do_box_wrapping(void)
{
}
#endif

#include "../data/simparticles.h"
template class domain<simparticles>;

#ifdef LIGHTCONE_PARTICLES
#include "../data/lcparticles.h"
template class domain<lcparticles>;
#endif
