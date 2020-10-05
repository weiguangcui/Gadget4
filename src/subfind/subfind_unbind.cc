/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind_unbind.cc
 *
 *  \brief carries out the gravitational unbinding of a subhalo candidate
 */

#include "gadgetconfig.h"

#ifdef SUBFIND

#include <gsl/gsl_math.h>
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
#include "../sort/peano.h"
#include "../subfind/subfind.h"
#include "../system/system.h"

#define MAX_UNBOUND_FRAC_BEFORE_BULK_VELOCITY_UPDATE 0.02
#define MAX_UNBOUND_FRAC_BEFORE_POTENTIAL_UPDATE 0.20

/* this function takes a list of particles given via their indices in d[] and subjects them
 * to a gravitational unbinding procedure. The number of bound particles is returned,
 * and the array d[] is updated accordingly.
 */
template <typename partset>
int fof<partset>::subfind_unbind(domain<partset> *D, MPI_Comm Communicator, int *d, int num)
{
  double fac_vel_to_phys, fac_hubbleflow, fac_comov_to_phys;
  subfind_get_factors(fac_vel_to_phys, fac_hubbleflow, fac_comov_to_phys);

  /* get the local communication context */
  int commNTask, commThisTask;
  MPI_Comm_size(Communicator, &commNTask);
  MPI_Comm_rank(Communicator, &commThisTask);

  typename partset::pdata *P = Tp->P;
  subfind_data *PS           = Tp->PS;

  /* we will start out by recomputing the potential for all particles based on all particles */
  int phaseflag = RECOMPUTE_ALL;

  int iter = 0;

  long long totnum = num;
  MPI_Allreduce(MPI_IN_PLACE, &totnum, 1, MPI_LONG_LONG, MPI_SUM, Communicator);

  int num_removed = 0;
  long long totremoved;

  int *dremoved  = (int *)Mem.mymalloc("dremoved", num * sizeof(int));
  double *potold = (double *)Mem.mymalloc("potold", num * sizeof(double));

  do
    {
      FoFGravTree.treeallocate(Tp->NumPart, Tp, D);

      if(phaseflag == RECOMPUTE_ALL)
        {
          FoFGravTree.treebuild(num, d);
        }
      else
        {
          FoFGravTree.treebuild(num_removed, dremoved);

          for(int i = 0; i < num; i++)
            potold[i] = PS[d[i]].u.s.u.DM_Potential;
        }

      /* let's compute the potential energy */
      subfind_potential_compute(D, num, d);

      FoFGravTree.treefree();

      if(phaseflag == RECOMPUTE_ALL)
        {
          /* subtract self-potential and convert to physical potential */
          for(int i = 0; i < num; i++)
            {
              int softtype = P[d[i]].getSofteningClass();

              PS[d[i]].u.s.u.DM_Potential += Tp->P[d[i]].getMass() / (All.ForceSoftening[softtype] / 2.8);
              PS[d[i]].u.s.u.DM_Potential *= All.G / fac_comov_to_phys;
            }
        }
      else
        {
          /* do not correct for self-potential and instead use calculated potential as correction to previous potential */
          for(int i = 0; i < num; i++)
            {
              PS[d[i]].u.s.u.DM_Potential *= All.G / fac_comov_to_phys;
              PS[d[i]].u.s.u.DM_Potential = potold[i] - PS[d[i]].u.s.u.DM_Potential;
            }
        }

      /* At this point, we have in num/d[] the list of still considered particles, and the potential is current */

      /* we will not unbind all particles with positive energy, until none are left. Because the kinetic energy depends on
       * the velocity frame, we recompute the bulk velocity and the center of mass after 2% of the particles have been removed
       */

      /* Determine in intpos[] the potential minimum among the particles,
       * which we take that as the center of the halo
       */
      int minindex = -1;
      struct
      {
        double pot;
        int rank;
      } local = {MAX_DOUBLE_NUMBER, commThisTask}, global;

      for(int i = 0; i < num; i++)
        if(PS[d[i]].u.s.u.DM_Potential < local.pot)
          {
            local.pot = PS[d[i]].u.s.u.DM_Potential;
            minindex  = d[i];
          }

      MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);

      MyIntPosType intpos[3]; /* potential minimum */

      if(commThisTask == global.rank)
        for(int j = 0; j < 3; j++)
          intpos[j] = P[minindex].IntPos[j];

      MPI_Bcast(intpos, 3 * sizeof(MyIntPosType), MPI_BYTE, global.rank, Communicator);

      /* we start with zero removed particles */
      num_removed = 0;

      long long totunbound;
      do
        {
          /* let's get bulk velocity and the center-of-mass */
          double massloc = 0, sloc[3] = {0, 0, 0}, vloc[3] = {0, 0, 0};

          for(int i = 0; i < num; i++)
            {
              int part_index = d[i];

              double dxyz[3];
              Tp->nearest_image_intpos_to_pos(P[part_index].IntPos, intpos, dxyz);

              for(int j = 0; j < 3; j++)
                sloc[j] += Tp->P[part_index].getMass() * dxyz[j];

              for(int j = 0; j < 3; j++)
                vloc[j] += Tp->P[part_index].getMass() * P[part_index].Vel[j];

              massloc += Tp->P[part_index].getMass();
            }

          double s[3], v[3], mass;
          MPI_Allreduce(sloc, s, 3, MPI_DOUBLE, MPI_SUM, Communicator);
          MPI_Allreduce(vloc, v, 3, MPI_DOUBLE, MPI_SUM, Communicator);
          MPI_Allreduce(&massloc, &mass, 1, MPI_DOUBLE, MPI_SUM, Communicator);

          for(int j = 0; j < 3; j++)
            {
              s[j] /= mass; /* center of mass offset relative to minimum potential */
              v[j] /= mass; /* center of mass velocity */
            }

          MySignedIntPosType off[3];
          Tp->pos_to_signedintpos(s, off);

          /* get integer version of absolute center of mass position */
          MyIntPosType int_cm[3];
          int_cm[0] = off[0] + intpos[0];
          int_cm[1] = off[1] + intpos[1];
          int_cm[2] = off[2] + intpos[2];

          double *bnd_energy = (double *)Mem.mymalloc("bnd_energy", num * sizeof(double));

          /* calculate the binding energies */
          for(int i = 0; i < num; i++)
            {
              int part_index = d[i];

              /* distance to center of mass */
              double dx[3];
              Tp->nearest_image_intpos_to_pos(P[part_index].IntPos, int_cm, dx);

              /* get physical velocity relative to center of mass */
              double dv[3];
              for(int j = 0; j < 3; j++)
                {
                  dv[j] = fac_vel_to_phys * (P[part_index].Vel[j] - v[j]);
                  dv[j] += fac_hubbleflow * fac_comov_to_phys * dx[j];
                }

              PS[part_index].v.DM_BindingEnergy =
                  PS[part_index].u.s.u.DM_Potential + 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
#ifndef LEAN
              if(P[part_index].getType() == 0)
                PS[part_index].v.DM_BindingEnergy += PS[part_index].Utherm;
#endif
              bnd_energy[i] = PS[part_index].v.DM_BindingEnergy;
            }

          /* sort by binding energy, highest energies (num_unbound / most weakly bound) first */
          mycxxsort_parallel(bnd_energy, bnd_energy + num, subfind_compare_binding_energy, Communicator);

          int *npart = (int *)Mem.mymalloc("npart", commNTask * sizeof(int));
          MPI_Allgather(&num, 1, MPI_INT, npart, 1, MPI_INT, Communicator);

          /* (global) index of limiting energy value for least tightly bound fraction  - those we may
           * remove at most in one iteration */
          long long j = std::max<long long>(5, (long long)(MAX_UNBOUND_FRAC_BEFORE_BULK_VELOCITY_UPDATE * totnum));

          /* find the processor where this lies */
          int task = 0;
          while(j >= npart[task])
            {
              j -= npart[task];
              task++;
            }

          double energy_limit = MAX_DOUBLE_NUMBER;

          if(commThisTask == task)
            energy_limit = bnd_energy[j];

          MPI_Allreduce(MPI_IN_PLACE, &energy_limit, 1, MPI_DOUBLE, MPI_MIN, Communicator);

          /* now unbind particles */
          int num_unbound = 0;

          for(int i = 0; i < num; i++)
            {
              int p = d[i];

              if(PS[p].v.DM_BindingEnergy > 0 && PS[p].v.DM_BindingEnergy > energy_limit)
                {
                  num_unbound++;

                  dremoved[num_removed++] = d[i];

                  d[i] = d[num - 1];
                  num--;
                  i--;
                }
            }

          Mem.myfree(npart);
          Mem.myfree(bnd_energy);

          totunbound = num_unbound;
          totremoved = num_removed;
          totnum     = num;

          MPI_Allreduce(MPI_IN_PLACE, &totunbound, 1, MPI_LONG_LONG, MPI_SUM, Communicator);
          MPI_Allreduce(MPI_IN_PLACE, &totremoved, 1, MPI_LONG_LONG, MPI_SUM, Communicator);
          MPI_Allreduce(MPI_IN_PLACE, &totnum, 1, MPI_LONG_LONG, MPI_SUM, Communicator);
        }
      while(totunbound > 0 && totnum >= All.DesLinkNgb && totremoved < MAX_UNBOUND_FRAC_BEFORE_POTENTIAL_UPDATE * totnum);

      iter++;

      if(iter > MAX_ITER_UNBIND)
        Terminate("too many iterations");

      if(phaseflag == RECOMPUTE_ALL)
        {
          if(totremoved > 0)
            phaseflag = UPDATE_ALL;
        }
      else
        {
          if(totremoved == 0)
            {
              phaseflag  = RECOMPUTE_ALL; /* this will make us repeat everything once more for all particles */
              totremoved = 1;             /* to make the code check once more all particles */
            }
        }
    }
  while(totremoved > 0 && totnum >= All.DesLinkNgb);

  Mem.myfree(potold);
  Mem.myfree(dremoved);

  return num;
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
