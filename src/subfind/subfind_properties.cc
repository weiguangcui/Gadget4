/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_properties.cc
 *
 *  \brief determination of subhalo properties
 */

#include "gadgetconfig.h"

#ifdef SUBFIND

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
#include "../fof/fof.h"
#include "../gravtree/gravtree.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/cxxsort.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"
#include "../time_integration/driftfac.h"

template <>
void fof<simparticles>::subfind_get_factors(double &fac_vel_to_phys, double &fac_hubbleflow, double &fac_comov_to_phys)
{
  if(All.ComovingIntegrationOn)
    {
      fac_vel_to_phys   = 1.0 / All.Time;  // converts to physical velocity
      fac_hubbleflow    = Driftfac.hubble_function(All.Time);
      fac_comov_to_phys = All.Time;  // converts comoving distance to physical distance
    }
  else
    {  // non-comoving integration
      fac_vel_to_phys   = 1;
      fac_hubbleflow    = 0;
      fac_comov_to_phys = 1;
    }
}

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
template <>
void fof<lcparticles>::subfind_get_factors(double &fac_vel_to_phys, double &fac_hubbleflow, double &fac_comov_to_phys)
{
  fac_vel_to_phys   = 1.0;  // the velocities on the lightcone are already peculiar velocities
  fac_hubbleflow    = Driftfac.hubble_function(Ascale);
  fac_comov_to_phys = Ascale;
}
#endif

template <typename partset>
int fof<partset>::subfind_determine_sub_halo_properties(int *d, int num, subhalo_properties *subhalo, MPI_Comm Communicator)
{
  /* get the local communication context */
  int commNTask, commThisTask;
  MPI_Comm_size(Communicator, &commNTask);
  MPI_Comm_rank(Communicator, &commThisTask);

  double fac_vel_to_phys, fac_hubbleflow, fac_comov_to_phys;
  subfind_get_factors(fac_vel_to_phys, fac_hubbleflow, fac_comov_to_phys);

  typename partset::pdata *P = Tp->P;
  subfind_data *PS           = Tp->PS;

  double mass_tab[NTYPES], halfmassradtype[NTYPES];
  long long len_type[NTYPES];
  for(int j = 0; j < NTYPES; j++)
    {
      mass_tab[j]        = 0;
      len_type[j]        = 0;
      halfmassradtype[j] = 0;
    }

#ifdef STARFORMATION
  double sfr = 0;
#endif

  int minindex  = -1;
  double minpot = MAX_DOUBLE_NUMBER;
  for(int i = 0; i < num; i++)
    {
      int p = d[i];
      if(PS[p].u.s.u.DM_Potential < minpot || minindex == -1)
        {
          minpot   = PS[p].u.s.u.DM_Potential;
          minindex = p;
        }
      len_type[P[p].getType()]++;

#ifdef STARFORMATION
      if(P[p].getType() == 0)
        sfr += Tp->SphP[p].Sfr;
#endif
    }

  MPI_Allreduce(MPI_IN_PLACE, len_type, NTYPES, MPI_LONG_LONG, MPI_SUM, Communicator);

  double *minpotlist = (double *)Mem.mymalloc("minpotlist", commNTask * sizeof(double));
  MPI_Allgather(&minpot, 1, MPI_DOUBLE, minpotlist, 1, MPI_DOUBLE, Communicator);

  int mincpu = -1;
  minpot     = MAX_DOUBLE_NUMBER;
  for(int i = 0; i < commNTask; i++)
    if(minpotlist[i] < minpot)
      {
        mincpu = i;
        minpot = minpotlist[mincpu];
      }

  Mem.myfree(minpotlist);

  if(mincpu < 0)
    Terminate("mincpu < 0");

  MyIntPosType intpos[3];

  if(commThisTask == mincpu)
    for(int j = 0; j < 3; j++)
      intpos[j] = P[minindex].IntPos[j];

  MPI_Bcast(intpos, 3 * sizeof(MyIntPosType), MPI_BYTE, mincpu, Communicator);

#ifdef STARFORMATION
  MPI_Allreduce(MPI_IN_PLACE, &sfr, 1, MPI_DOUBLE, MPI_SUM, Communicator);
#endif

  /* pos[] now holds the position of minimum potential */
  /* we'll take it that as the center */

  /* determine the particle ID with the smallest binding energy */
  minindex = -1;
  minpot   = MAX_DOUBLE_NUMBER;
  for(int i = 0; i < num; i++)
    {
      int p = d[i];
      if(PS[p].v.DM_BindingEnergy < minpot || minindex == -1)
        {
          minpot   = PS[p].v.DM_BindingEnergy;
          minindex = p;
        }
    }

  MyIDType mostboundid;

  minpotlist = (double *)Mem.mymalloc("minpotlist", commNTask * sizeof(double));
  MPI_Allgather(&minpot, 1, MPI_DOUBLE, minpotlist, 1, MPI_DOUBLE, Communicator);

  mincpu = -1;
  minpot = MAX_DOUBLE_NUMBER;
  for(int i = 0; i < commNTask; i++)
    if(minpotlist[i] < minpot)
      {
        mincpu = i;
        minpot = minpotlist[mincpu];
      }

  Mem.myfree(minpotlist);

  if(mincpu < 0)
    Terminate("mincpu < 0");

  int marked = 0;

  if(commThisTask == mincpu)
    {
      mostboundid = P[minindex].ID.get();
#if defined(SUBFIND_ORPHAN_TREATMENT)
      if(!P[minindex].ID.is_previously_most_bound())
        {
          P[minindex].ID.mark_as_formerly_most_bound();
          marked++;
        }
#endif
    }

  MPI_Bcast(&mostboundid, sizeof(mostboundid), MPI_BYTE, mincpu, Communicator);

  /* let's get bulk velocity and the center-of-mass */
  /* here we still take all particles */

  double s[3] = {0, 0, 0}, v[3] = {0, 0, 0};
  double mass = 0;

#if defined(SUBFIND_ORPHAN_TREATMENT)
  int totprevmostboundlen = 0;
#endif

  for(int i = 0; i < num; i++)
    {
      int p = d[i];

      double dxyz[3];
      Tp->nearest_image_intpos_to_pos(P[p].IntPos, intpos, dxyz);

      s[0] += Tp->P[p].getMass() * dxyz[0];
      s[1] += Tp->P[p].getMass() * dxyz[1];
      s[2] += Tp->P[p].getMass() * dxyz[2];

      for(int j = 0; j < 3; j++)
        v[j] += Tp->P[p].getMass() * P[p].Vel[j];

      mass += Tp->P[p].getMass();

      int ptype = P[p].getType();
      mass_tab[ptype] += Tp->P[p].getMass();
#if defined(SUBFIND_ORPHAN_TREATMENT)
      if(P[p].ID.is_previously_most_bound())
        totprevmostboundlen++;
#endif
    }

  double stot[3], vtot[3], masstot, mass_tabtot[NTYPES];

  MPI_Allreduce(s, stot, 3, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(v, vtot, 3, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(mass_tab, mass_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, Communicator);

#if defined(SUBFIND_ORPHAN_TREATMENT)
  MPI_Allreduce(MPI_IN_PLACE, &totprevmostboundlen, 1, MPI_INT, MPI_SUM, Communicator);
#endif

  mass = masstot;
  for(int j = 0; j < 3; j++)
    {
      s[j] = stot[j];
      v[j] = vtot[j];
    }

  for(int j = 0; j < NTYPES; j++)
    mass_tab[j] = mass_tabtot[j];

  double vel[3];
  for(int j = 0; j < 3; j++)
    {
      s[j] /= mass; /* center of mass */
      v[j] /= mass;
      vel[j] = fac_vel_to_phys * v[j];
    }

  MySignedIntPosType off[3];
  Tp->pos_to_signedintpos(s, off);

  MyIntPosType int_cm[3];

  int_cm[0] = off[0] + intpos[0];
  int_cm[1] = off[1] + intpos[1];
  int_cm[2] = off[2] + intpos[2];

  double cm[3];
  Tp->intpos_to_pos(int_cm, cm);

  double disp = 0, lx = 0, ly = 0, lz = 0;

  sort_r2list *rr_list = (sort_r2list *)Mem.mymalloc("rr_list", sizeof(sort_r2list) * (num + 1));

  for(int i = 0; i < num; i++)
    {
      int p = d[i];

      double dx[3], dv[3];

      /* relative to center of mass */

      Tp->nearest_image_intpos_to_pos(P[p].IntPos, int_cm, dx);

      dx[0] *= fac_comov_to_phys;
      dx[1] *= fac_comov_to_phys;
      dx[2] *= fac_comov_to_phys;

      for(int j = 0; j < 3; j++)
        {
          dv[j] = fac_vel_to_phys * (P[p].Vel[j] - v[j]);
          dv[j] += fac_hubbleflow * dx[j];

          disp += Tp->P[p].getMass() * dv[j] * dv[j];
        }

      lx += Tp->P[p].getMass() * (dx[1] * dv[2] - dx[2] * dv[1]);
      ly += Tp->P[p].getMass() * (dx[2] * dv[0] - dx[0] * dv[2]);
      lz += Tp->P[p].getMass() * (dx[0] * dv[1] - dx[1] * dv[0]);

      /* for rotation curve computation, take minimum of potential as center */

      double dxyz[3];
      Tp->nearest_image_intpos_to_pos(P[p].IntPos, intpos, dxyz);

      for(int j = 0; j < 3; j++)
        dxyz[j] *= fac_comov_to_phys;

      double r2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];

      double r = sqrt(r2);

      rr_list[i].mass = Tp->P[p].getMass();
      rr_list[i].r    = r;
    }

  double spintot[3], disploc = disp;
  double spinloc[3] = {lx, ly, lz};
  MPI_Allreduce(spinloc, spintot, 3, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(&disploc, &disp, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  lx = spintot[0];
  ly = spintot[1];
  lz = spintot[2];

  double spin[3] = {lx / mass, ly / mass, lz / mass};

  double veldisp = sqrt(disp / (3 * mass)); /* convert to 1d velocity dispersion */

  mycxxsort_parallel(rr_list, rr_list + num, subfind_compare_dist_rotcurve, Communicator);

  /* calculate cumulative mass */
  for(int i = 1; i < num; i++)
    rr_list[i].mass = rr_list[i - 1].mass + rr_list[i].mass;

  double halfmassrad = 0;
  double max = 0, maxrad = 0;

  double mass_part = 0;
  if(num)
    mass_part = rr_list[num - 1].mass;
  double *masslist = (double *)Mem.mymalloc("masslist", commNTask * sizeof(double));
  MPI_Allgather(&mass_part, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, Communicator);

  double massbefore = 0;
  for(int i = 0; i < commThisTask; i++)
    massbefore += masslist[i];

  for(int i = 0; i < num; i++)
    rr_list[i].mass += massbefore;

  Mem.myfree(masslist);

  /* now calculate rotation curve maximum and half mass radius */

  double halfmassrad_loc  = 0;
  sort_r2list *rr_lowlist = (sort_r2list *)Mem.mymalloc("rr_lowlist", commNTask * sizeof(sort_r2list));
  sort_r2list low_element;
  if(num > 0)
    low_element = rr_list[0];
  else
    {
      low_element.mass = 0;
      low_element.r    = 0;
    }
  MPI_Allgather(&low_element, sizeof(sort_r2list), MPI_BYTE, rr_lowlist, sizeof(sort_r2list), MPI_BYTE, Communicator);

  rr_list[num].mass = 0;
  rr_list[num].r    = 0;

  for(int j = commThisTask + 1; j < commNTask; j++)
    if(rr_lowlist[j].mass > 0)
      {
        rr_list[num] = rr_lowlist[j];
        break;
      }

  Mem.myfree(rr_lowlist);

  int *numlist = (int *)Mem.mymalloc("numlist", commNTask * sizeof(int));
  MPI_Allgather(&num, 1, MPI_INT, numlist, 1, MPI_INT, Communicator);

  int nbefore = 0;
  for(int i = 0; i < commThisTask; i++)
    nbefore += numlist[i];

  for(int i = num - 1; i >= 0; i--)
    {
      if((i + nbefore) > 5 && rr_list[i].mass > max * rr_list[i].r)
        {
          max    = rr_list[i].mass / rr_list[i].r;
          maxrad = rr_list[i].r;
        }

      if(rr_list[i].mass < 0.5 * mass && rr_list[i + 1].mass >= 0.5 * mass)
        halfmassrad_loc = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
    }

  Mem.myfree(numlist);

  MPI_Allreduce(&halfmassrad_loc, &halfmassrad, 1, MPI_DOUBLE, MPI_MAX, Communicator);
  double *maxlist    = (double *)Mem.mymalloc("maxlist", commNTask * sizeof(double));
  double *maxradlist = (double *)Mem.mymalloc("maxradlist", commNTask * sizeof(double));
  MPI_Allgather(&max, 1, MPI_DOUBLE, maxlist, 1, MPI_DOUBLE, Communicator);
  MPI_Allgather(&maxrad, 1, MPI_DOUBLE, maxradlist, 1, MPI_DOUBLE, Communicator);

  max = maxrad = 0;
  for(int i = 0; i < commNTask; i++)
    {
      if(maxlist[i] > max)
        {
          max    = maxlist[i];
          maxrad = maxradlist[i];
        }
    }
  Mem.myfree(maxradlist);
  Mem.myfree(maxlist);

  double vmax    = sqrt(All.G * max);
  double vmaxrad = maxrad;

  Mem.myfree(rr_list);

  /* half mass radii for different types */

  int len_type_loc[NTYPES];
  for(int j = 0; j < NTYPES; j++)
    len_type_loc[j] = 0;

  for(int i = 0; i < num; i++)
    {
      int p     = d[i];
      int ptype = P[p].getType();
      len_type_loc[ptype]++;
    }

  for(int type = 0; type < NTYPES; type++)
    {
      sort_r2list *rr_list = (sort_r2list *)Mem.mymalloc("rr_list", sizeof(sort_r2list) * (len_type_loc[type] + 1));
      int itmp             = 0;
      for(int i = 0; i < num; i++)
        {
          int p = d[i];

          int ptype = P[p].getType();

          if(ptype == type)
            {
              double dxyz[3];
              Tp->nearest_image_intpos_to_pos(P[p].IntPos, intpos, dxyz);

              for(int j = 0; j < 3; j++)
                dxyz[j] *= fac_comov_to_phys;

              double r2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];
              double r  = sqrt(r2);

              rr_list[itmp].mass = Tp->P[p].getMass();
              rr_list[itmp].r    = r;
              itmp++;
            }
        }

      if(itmp != len_type_loc[type])
        Terminate("should not occur: %d %d", itmp, len_type_loc[type]);

      mycxxsort_parallel(rr_list, rr_list + len_type_loc[type], subfind_compare_dist_rotcurve, Communicator);

      /* calculate cumulative mass */
      for(int i = 1; i < len_type_loc[type]; i++)
        rr_list[i].mass = rr_list[i - 1].mass + rr_list[i].mass;

      double mass_part = 0;
      if(len_type_loc[type])
        mass_part = rr_list[len_type_loc[type] - 1].mass;
      double *masslist = (double *)Mem.mymalloc("masslist", commNTask * sizeof(double));
      MPI_Allgather(&mass_part, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, Communicator);

      double massbefore = 0;
      for(int i = 0; i < commThisTask; i++)
        massbefore += masslist[i];

      for(int i = 0; i < len_type_loc[type]; i++)
        rr_list[i].mass += massbefore;

      Mem.myfree(masslist);

      /* now calculate half mass radii */
      double halfmassrad_loc  = 0;
      sort_r2list *rr_lowlist = (sort_r2list *)Mem.mymalloc("rr_lowlist", commNTask * sizeof(sort_r2list));
      sort_r2list low_element;
      if(len_type_loc[type] > 0)
        low_element = rr_list[0];
      else
        {
          low_element.mass = 0;
          low_element.r    = 0;
        }

      MPI_Allgather(&low_element, sizeof(sort_r2list), MPI_BYTE, rr_lowlist, sizeof(sort_r2list), MPI_BYTE, Communicator);

      rr_list[len_type_loc[type]].mass = 0;
      rr_list[len_type_loc[type]].r    = 0;
      for(int j = commThisTask + 1; j < commNTask; j++)
        if(rr_lowlist[j].mass > 0)
          {
            rr_list[len_type_loc[type]] = rr_lowlist[j];
            break;
          }

      Mem.myfree(rr_lowlist);

      for(int i = len_type_loc[type] - 1; i >= 0; i--)
        {
          if(rr_list[i].mass < 0.5 * mass_tab[type] && rr_list[i + 1].mass >= 0.5 * mass_tab[type])
            halfmassrad_loc = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
        }

      MPI_Allreduce(&halfmassrad_loc, &halfmassradtype[type], 1, MPI_DOUBLE, MPI_MAX, Communicator);

      Mem.myfree(rr_list);
    }

    /* properties of star forming gas */
#ifdef STARFORMATION
  double gasMassSfr = 0;
  for(int i = 0; i < num; i++)
    {
      int p = d[i];

      if(P[p].getType() == 0)
        if(Tp->SphP[p].Sfr > 0)
          gasMassSfr += P[p].getMass();
    }
#endif

#ifdef STARFORMATION
  double gasMassSfrtot;
  MPI_Allreduce(&gasMassSfr, &gasMassSfrtot, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  gasMassSfr = gasMassSfrtot;
#endif

  long long totlen;
  sumup_large_ints(1, &num, &totlen, Communicator);

  /* now store the calculated properties in the subhalo structure */
  if(commThisTask == 0)
    {
      subhalo->Len  = totlen;
      subhalo->Mass = mass;

      for(int j = 0; j < NTYPES; j++)
        {
          subhalo->MassType[j]           = mass_tab[j];
          subhalo->LenType[j]            = len_type[j];
          subhalo->SubHalfMassRadType[j] = halfmassradtype[j];
        }

      double pos[3];
      fof_get_halo_position(intpos, pos);

      for(int j = 0; j < 3; j++)
        {
          subhalo->IntPos[j] = intpos[j];
          subhalo->Pos[j]    = pos[j];
          subhalo->Vel[j]    = vel[j];
          subhalo->CM[j]     = cm[j];
          subhalo->Spin[j]   = spin[j];
        }

      subhalo->SubMostBoundID = mostboundid;
      subhalo->SubVelDisp     = veldisp;
      subhalo->SubVmax        = vmax;
      subhalo->SubVmaxRad     = vmaxrad;
      subhalo->SubHalfMassRad = halfmassrad;

#if defined(SUBFIND_ORPHAN_TREATMENT)
      subhalo->SubhaloLenPrevMostBnd = totprevmostboundlen;
#endif
#ifdef STARFORMATION
      subhalo->Sfr        = sfr;
      subhalo->GasMassSfr = gasMassSfr;
#endif
    }

  return marked;
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
