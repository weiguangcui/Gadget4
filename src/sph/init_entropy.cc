/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  init_entropy.cc
 *
 *  \brief initialization code for the entropy variable of the particles
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../io/io.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../sort/peano.h"
#include "../sph/kernel.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

#ifdef PRESSURE_ENTROPY_SPH

#define MAX_ITER_ENTROPY 100
#define ENTROPY_TOLERANCE 1.0e-5

/*! \file init_entropy.c
 *  \brief SPH entropy computation from internal energies for pressure-entropy formulation of SPH
 *
 */
void sph::init_entropy(void)
{
  TIMER_STORE;
  TIMER_START(CPU_DENSITY);

  D->mpi_printf("SPH-INIT-ENTROPY: Begin entropy calculation.  (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());

  /* Create list of targets. We do this here to simplify the treatment later on */
  int *targetList = (int *)Mem.mymalloc("targetlist", Tp->NumGas * sizeof(int));

  int ndensities = 0;

  for(int i = 0; i < Tp->TimeBinsHydro.NActiveParticles; i++)
    {
      int target = Tp->TimeBinsHydro.ActiveParticleList[i];
      if(target < 0 || target >= Tp->NumGas)
        Terminate("target=%d i=%d\n", target, i);
      targetList[ndensities++] = target;
    }

  int iter = 0;

  // let's grab at most half the still available memory for imported points and nodes
  int nspace = (0.33 * Mem.FreeBytes) / (sizeof(ngbnode) + 8 * sizeof(foreign_sphpoint_data));

  MaxForeignNodes  = nspace;
  MaxForeignPoints = 8 * nspace;
  NumForeignNodes  = 0;
  NumForeignPoints = 0;

  sum_NumForeignNodes  = 0;
  sum_NumForeignPoints = 0;

  /* the following two arrays hold imported tree nodes and imported points to augment the local tree */
  Foreign_Nodes  = (ngbnode *)Mem.mymalloc_movable(&Foreign_Nodes, "Foreign_Nodes", MaxForeignNodes * sizeof(ngbnode));
  Foreign_Points = (foreign_sphpoint_data *)Mem.mymalloc_movable(&Foreign_Points, "Foreign_Points",
                                                                 MaxForeignPoints * sizeof(foreign_sphpoint_data));



  double tstart = Logs.second();

  int global_left_particles = 0;

  MPI_Allreduce(&ndensities, &global_left_particles, 1, MPI_INT, MPI_SUM, D->Communicator);

  do
    {
      double t0 = Logs.second();

 /*  Since EntropyToInvGammaPred of remote particles can change, we have to import the particles in every iteration */

      MaxForeignNodes  = nspace;
      MaxForeignPoints = 8 * nspace;
      NumForeignNodes  = 0;
      NumForeignPoints = 0;

      sum_NumForeignNodes  = 0;
      sum_NumForeignPoints = 0;

      tree_initialize_leaf_node_access_info();

      max_ncycles = 0;

      prepare_shared_memory_access();

      /* now do the primary work with this call */
      densities_determine(ndensities, targetList);

      MPI_Allreduce(MPI_IN_PLACE, &max_ncycles, 1, MPI_INT, MPI_MAX, D->Communicator);

      cleanup_shared_memory_access();

      /* do final operations on results */

      double entropy_old;

      int npleft = 0;

      for(int i = 0; i < ndensities; i++)
        {
          int target = targetList[i];
          if(target >= 0)
            {
              if(Tp->P[target].getType() != 0)
                Terminate("P[target].getType() != 0");

              sph_particle_data *SphP = Tp->SphP;

              if(SphP[target].EntropyToInvGammaPred > 0 && SphP[target].Density > 0)
                {
                  entropy_old = SphP[target].Entropy;
                  SphP[target].PressureSphDensity /= SphP[target].EntropyToInvGammaPred;
                  SphP[target].Entropy =
                      GAMMA_MINUS1 * SphP[target].EntropyPred / pow(SphP[target].PressureSphDensity * All.cf_a3inv, GAMMA_MINUS1);
                  SphP[target].EntropyToInvGammaPred = pow(SphP[target].Entropy, 1.0 / GAMMA);
                }
              else
                {
                  entropy_old                        = SphP[target].Entropy;
                  SphP[target].PressureSphDensity    = 0;
                  SphP[target].Entropy               = 0;
                  SphP[target].EntropyToInvGammaPred = 0;
                }
              /* entropy has not converged yet */
              if(fabs(entropy_old - SphP[target].Entropy) > ENTROPY_TOLERANCE * entropy_old)
                targetList[npleft++] = target;
            }
        }

      ndensities = npleft;

      MPI_Allreduce(&ndensities, &global_left_particles, 1, MPI_INT, MPI_SUM, D->Communicator);

      double t1 = Logs.second();

      if(npleft > 0)
        {
          iter++;

          D->mpi_printf("SPH-INIT-ENTROPY: ngb iteration %4d: took %8.3f  , need to repeat for %012lld local particles.\n", iter,
                        Logs.timediff(t0, t1), npleft);

          if(iter > MAXITER)
            Terminate("failed to converge in neighbour iteration in density()\n");
        }
      else
        D->mpi_printf("SPH-INIT-ENTROPY: ngb iteration %4d: took %8.3f\n", ++iter, Logs.timediff(t0, t1));
    }
  while(global_left_particles > 0);


  TIMER_STOP(CPU_DENSITY);


  /* free temporary buffers */

  Mem.myfree(Foreign_Points);
  Mem.myfree(Foreign_Nodes);

  Mem.myfree(targetList);

  double tb = Logs.second();

  D->mpi_printf("SPH-INIT-ENTROPY: entropy calculation is done. took: %8.3f\n", Logs.timediff(tstart, tb));
}

/*! \brief This function is used to find the initial entropy^invgamma for each SPH
 *  particle in the pressure-entropy formulation of SPH in case the ICs
 *  file contains internal energies.
 */
void sph::setup_entropy_to_invgamma(void)
{
  All.set_cosmo_factors_for_current_time();

  /* Initialize entropy and entropy^invgamma with a fist guess coming from standard SPH density estimate. */
  /* EntropyPred is untouched since it contains the internal energies needed for the iterative process; it */
  /* will be set in init.c to the correct value.                                                          */
  for(int i = 0; i < Tp->NumGas; i++)
    {
      Tp->SphP[i].Entropy = GAMMA_MINUS1 * Tp->SphP[i].EntropyPred / pow(Tp->SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1);
      Tp->SphP[i].EntropyToInvGammaPred = pow(Tp->SphP[i].Entropy, 1.0 / GAMMA);
    }

  init_entropy();
}
#endif
