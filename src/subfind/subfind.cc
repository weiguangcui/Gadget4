/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  subfind.cc
 *
 *  \brief main routines for carrying out the SUBFIND or SUBFIND_HBT algorithms on a set of FOF groups
 */

#include "gadgetconfig.h"

#ifdef SUBFIND

#include <mpi.h>
#include <unistd.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fof/fof.h"
#include "../fof/fof_io.h"
#include "../gravtree/gravtree.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mergertree/mergertree.h"
#include "../mpi_utils/mpi_utils.h"
#include "../sort/cxxsort.h"
#include "../sort/parallel_sort.h"
#include "../subfind/subfind.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

/* main routine of SUBFIND algorithm, for output number 'num'.
 */
template <typename partset>
void fof<partset>::subfind_find_subhalos(int num, const char *basename, const char *grpcat_dirbasename)
{
  TIMER_START(CPU_SUBFIND);

  double tstart = Logs.second();

  mpi_printf("\nSUBFIND: We now execute a parallel version of SUBFIND.\n");

#ifdef MERGERTREE
  double lensum = 0, partcount = 0;
  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(Tp->P[i].PrevSizeOfSubhalo.get() > 0)
        partcount += 1;
      lensum += Tp->P[i].PrevSizeOfSubhalo.get();
      Tp->PS[i].SizeOfSubhalo.set(0);  // initialize
    }
  MPI_Allreduce(MPI_IN_PLACE, &partcount, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  MPI_Allreduce(MPI_IN_PLACE, &lensum, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  mpi_printf("SUBFIND: Previous subhalo catalogue had approximately a size %g, and the summed squared subhalo size was %g\n",
             partcount, lensum);
#endif

#ifndef SUBFIND_HBT

  /* let's determine the total matter density around each particle/cell. For each point,
   * we determine a smoothing radius based on a certain number of dark matter particles
   * from the primary link type.
   */
  FoFGravTree.treeallocate(Tp->NumPart, Tp, FoFDomain);

  TIMER_STOP(CPU_SUBFIND);
  FoFGravTree.treebuild(Tp->NumPart, NULL);
  TIMER_START(CPU_SUBFIND);

  subfind_density_hsml_guess();

  /* find smoothing lengths for primary and secondary points in groups (counting primary particles only for setting the smoothing
   * length) and then estimate total matter density around them
   */
  double cputime = subfind_density();

  mpi_printf("SUBFIND: iteration to correct primary neighbor count and density estimate took %g sec\n", cputime);

  FoFGravTree.treefree();
#endif

  All.set_cosmo_factors_for_current_time();

  /* fill internal energy into auxiliary PS[] array */
#ifndef LEAN
  for(int i = 0; i < Tp->NumPart; i++)
    if(Tp->P[i].getType() == 0)
      Tp->PS[i].Utherm = Tp->get_utherm_from_entropy(i);
    else
      Tp->PS[i].Utherm = 0;
#endif

  double GroupSizeThresholdFactor = 0.55;
  int ncount = 0, nprocs = 0;
  long long seriallen = 0;
  long long sum_seriallen;

  /* Let's set a fiducial size limit for the maximum group size before we select the collective subfind algorithm */
  do
    {
      GroupSizeThresholdFactor += 0.05;
      ncount = nprocs = seriallen = 0;

      MaxSerialGroupLen = (int)(GroupSizeThresholdFactor * Tp->TotNumPart / NTask);

      /* Count, how many groups are above this limit, and how many processors we need for them */
      for(int i = 0; i < Ngroups; i++)
        if(Group[i].Len > MaxSerialGroupLen)
          {
            ncount++;
            nprocs += ((Group[i].Len - 1) / MaxSerialGroupLen) + 1;
          }
        else
          seriallen += Group[i].Len;

      MPI_Allreduce(&ncount, &Ncollective, 1, MPI_INT, MPI_SUM, Communicator);
      MPI_Allreduce(&nprocs, &NprocsCollective, 1, MPI_INT, MPI_SUM, Communicator);
      sumup_longs(1, &seriallen, &sum_seriallen, Communicator);
    }
  while(NprocsCollective + ((sum_seriallen > 0) ? 1 : 0) > NTask);

  mpi_printf("SUBFIND: Number of FOF halos treated with collective SubFind algorithm = %d\n", Ncollective);
  mpi_printf("SUBFIND: Number of processors used in different partitions for the collective SubFind code = %d\n", NprocsCollective);
  mpi_printf("SUBFIND: (The adopted size-limit for the collective algorithm was %d particles, for threshold size factor %g)\n",
             MaxSerialGroupLen, GroupSizeThresholdFactor);
  mpi_printf("SUBFIND: The other %lld FOF halos are treated in parallel with serial code\n", TotNgroups - Ncollective);

  /* set up a global table that informs about the processor assignment of the groups that are treated collectively */
  ProcAssign = (proc_assign_data *)Mem.mymalloc_movable(&ProcAssign, "ProcAssign", Ncollective * sizeof(proc_assign_data));
  proc_assign_data *locProcAssign = (proc_assign_data *)Mem.mymalloc("locProcAssign", ncount * sizeof(proc_assign_data));

  ncount = 0;
  for(int i = 0; i < Ngroups; i++)
    if(Group[i].Len > MaxSerialGroupLen)
      {
        locProcAssign[ncount].GroupNr = Group[i].GroupNr;
        locProcAssign[ncount].Len     = Group[i].Len;
        ncount++;
      }

  /* gather the information on the collective groups accross all CPUs */
  int *recvcounts = (int *)Mem.mymalloc("recvcounts", sizeof(int) * NTask);
  int *bytecounts = (int *)Mem.mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)Mem.mymalloc("byteoffset", sizeof(int) * NTask);

  MPI_Allgather(&ncount, 1, MPI_INT, recvcounts, 1, MPI_INT, Communicator);

  for(int task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(proc_assign_data);

  byteoffset[0] = 0;
  for(int task = 1; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  myMPI_Allgatherv(locProcAssign, bytecounts[ThisTask], MPI_BYTE, ProcAssign, bytecounts, byteoffset, MPI_BYTE, Communicator);

  Mem.myfree(byteoffset);
  Mem.myfree(bytecounts);
  Mem.myfree(recvcounts);
  Mem.myfree(locProcAssign);

  /* make sure, the table is sorted in ascending group-number order */
  mycxxsort(ProcAssign, ProcAssign + Ncollective, subfind_compare_procassign_GroupNr);

  /* assign the processor sets for the collective groups and set disjoint color-flag to later split the processors into different
   * communicators */
  nprocs         = 0;
  CommSplitColor = ThisTask; /* by default, this places every processor in his own processor group */

  /* now we assign the same color for groups of CPUs that each do one halo collectively */
  for(int i = 0; i < Ncollective; i++)
    {
      ProcAssign[i].FirstTask = nprocs;
      ProcAssign[i].NTask     = ((ProcAssign[i].Len - 1) / MaxSerialGroupLen) + 1;
      nprocs += ProcAssign[i].NTask;

      if(ThisTask >= ProcAssign[i].FirstTask && ThisTask < (ProcAssign[i].FirstTask + ProcAssign[i].NTask))
        CommSplitColor = i;
    }

  /* Now assign a target task for each local group. For collective groups, the target task is the master in the CPU set, whereas
   * the serial ones are distributed in a round-robin fashion to the remaining CPUs.
   */
  for(int i = 0; i < Ngroups; i++)
    {
      if(Group[i].Len > MaxSerialGroupLen) /* we have a collective group */
        {
          if(Group[i].GroupNr >= Ncollective || Group[i].GroupNr < 0)
            Terminate("odd");
          Group[i].TargetTask = ProcAssign[Group[i].GroupNr].FirstTask;
        }
      else
        Group[i].TargetTask = ((Group[i].GroupNr - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective;
    }

  /* distribute the groups */
  subfind_distribute_groups();

  /* sort the local groups by group number */
  mycxxsort(Group, Group + Ngroups, fof_compare_Group_GroupNr);

  /* assign target CPUs for the particles in groups */
  /* the particles not in groups will be distributed such that a close to uniform particle load results */
  double t0           = Logs.second();
  int *count_loc_task = (int *)Mem.mymalloc_clear("count_loc_task", NTask * sizeof(int));
  int *count_task     = (int *)Mem.mymalloc("count_task", NTask * sizeof(int));
  int *count_free     = (int *)Mem.mymalloc("count_free", NTask * sizeof(int));
  int count_loc_free  = 0;

  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(Tp->PS[i].GroupNr.get() < HALONR_MAX) /* particle is in a group */
        {
          if(Tp->PS[i].GroupNr.get() < (MyIDType)Ncollective) /* we are in a collectively treated group */
            Tp->PS[i].TargetTask = ProcAssign[Tp->PS[i].GroupNr.get()].FirstTask + (i % ProcAssign[Tp->PS[i].GroupNr.get()].NTask);
          else /* otherwise the whole group will be treated by a single CPU */
            Tp->PS[i].TargetTask = ((Tp->PS[i].GroupNr.get() - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective;

          count_loc_task[Tp->PS[i].TargetTask]++;
        }
      else
        {
          Tp->PS[i].TargetTask = ThisTask;  // default is that we stay
          count_loc_free++;
        }

      Tp->PS[i].TargetIndex = 0; /* unimportant here */
    }

  /* get a list of how many unassigned and thus in principle movable particles every processor has */
  MPI_Allgather(&count_loc_free, 1, MPI_INT, count_free, 1, MPI_INT, Communicator);

  /* determine how many particles from groups will end up on each processor */
  MPI_Allreduce(count_loc_task, count_task, NTask, MPI_INT, MPI_SUM, Communicator);

  /* get the total particle count */
  long long sum = 0;
  for(int i = 0; i < NTask; i++)
    sum += count_task[i] + count_free[i];

  /* find out the average maximum load, and what the local head-room to it would be */
  int maxload = (sum + NTask - 1) / NTask;
  for(int i = 0; i < NTask; i++)
    {
      count_task[i] = maxload - count_task[i]; /* this is the amount that can fit on this task */
      if(count_task[i] < 0)
        count_task[i] = 0;
    }

  /* let's see how many particles can stay without degrading the load balance */
  for(int i = 0; i < NTask; i++)
    {
      if(count_free[i] <= count_task[i])
        {
          count_task[i] -= count_free[i];
          count_free[i] = 0;
        }
      else
        {
          count_free[i] -= count_task[i];
          count_task[i] = 0;
        }
    }

  /* now we determine the optimum target task for the subset of local particles that should be moved to another rank for better load
   * balance */
  int current_task = 0;
  for(int i = 0; i < ThisTask; i++)
    {
      while(count_free[i] > 0 && current_task < NTask)
        {
          if(count_free[i] < count_task[current_task])
            {
              count_task[current_task] -= count_free[i];
              count_free[i] = 0;
            }
          else
            {
              count_free[i] -= count_task[current_task];
              count_task[current_task] = 0;
              current_task++;
            }
        }
    }

  /* move a subset of particles such that close to uniform load is achieved */
  for(int i = 0; i < Tp->NumPart && count_free[ThisTask] > 0; i++)
    {
      if(Tp->PS[i].GroupNr.get() == HALONR_MAX)
        {
          /* particle not in a group. Could in principle stay but we move it such that a good load balance is obtained. */
          while(count_task[current_task] == 0 && current_task < NTask - 1)
            current_task++;

          Tp->PS[i].TargetTask = current_task;
          count_task[current_task]--;
          count_free[ThisTask]--;
        }
    }

  Mem.myfree(count_free);
  Mem.myfree(count_task);
  Mem.myfree(count_loc_task);

  /* report current balance */
  double balance = subfind_get_particle_balance();
  mpi_printf("SUBFIND: particle balance=%g\n", balance);

  /* distribute particles such that groups are completely on the CPU(s) that do the corresponding group(s) */
  FoFDomain->particle_exchange_based_on_PS(Communicator);
  double t1 = Logs.second();
  mpi_printf("SUBFIND: subfind_exchange() took %g sec\n", Logs.timediff(t0, t1));

  /* report balance that has been achieved */
  balance = subfind_get_particle_balance();
  mpi_printf("SUBFIND: particle balance for processing=%g\n", balance);

  /* lets estimate the maximum number of substructures we may find and need to store on the local CPU */
  if(ThisTask < NprocsCollective)
    {
      MaxNsubhalos = (ProcAssign[CommSplitColor].Len / ProcAssign[CommSplitColor].NTask) / All.DesLinkNgb;
    }
  else
    {
      long long nlocid = 0;
      for(int i = 0; i < Ngroups; i++)
        nlocid += Group[i].Len;

      MaxNsubhalos = nlocid / All.DesLinkNgb; /* should be a quite conservative upper limit */
    }

  Mem.myfree(ProcAssign);

  // some log-variables
  count_decisions           = 0;
  count_different_decisions = 0;

  /* allocate storage space for locally found subhalos */
  Nsubhalos = 0;
  Subhalo   = (subhalo_properties *)Mem.mymalloc_movable(&Subhalo, "Subhalo", MaxNsubhalos * sizeof(subhalo_properties));

  /* we can now split the communicator to give each collectively treated group its own processor set */
  MPI_Comm_split(Communicator, CommSplitColor, ThisTask, &SubComm);
  MPI_Comm_size(SubComm, &SubNTask);
  MPI_Comm_rank(SubComm, &SubThisTask);

  double ti0 = Logs.second();
  {
    domain<partset> SubDomain(SubComm, Tp);

    if(SubDomain.NumNodes != 0)
      Terminate("SubDomain.NumNodes=%d\n", SubDomain.NumNodes);

    if(CommSplitColor < Ncollective)
      {
        /* for the collective algorithm, we now do a further reshuffling to distribute the particles
         * in density-order to the processors. This is here only a performance improvement later on
         * for out distributed link-lists, reducing the number of MPI messages that need to be sent.
         * */

        /* allocated storage for an auxiliary array needed for sorting */
        sort_as_data *as = (sort_as_data *)Mem.mymalloc_movable(&as, "as", Tp->NumPart * sizeof(sort_as_data));

        for(int i = 0; i < Tp->NumPart; i++)
          {
            as[i].density = Tp->PS[i].u.s.u.DM_Density;
            as[i].origin  = (((long long)SubThisTask) << 32) + i;
          }

        /* sort according to density  */
        mycxxsort_parallel(as, as + Tp->NumPart, subfind_compare_as_density, SubComm);

        for(int i = 0; i < Tp->NumPart; i++)
          as[i].targettask = SubThisTask;

        /* revert to original order */
        mycxxsort_parallel(as, as + Tp->NumPart, subfind_compare_as_origin, SubComm);

        for(int i = 0; i < Tp->NumPart; i++)
          Tp->PS[i].TargetTask = as[i].targettask;

        Mem.myfree(as);

        SubDomain.particle_exchange_based_on_PS(SubComm);

        if(SubDomain.NumNodes != 0)
          Terminate("SubDomain.NumNodes=%d\n", SubDomain.NumNodes);

        /* now call the routine that does the actual processing of the groups */
        subfind_processing(&SubDomain, COLL_SUBFIND); /* we are one of the CPUs that does a collective group */
      }
    else
      subfind_processing(&SubDomain, SERIAL_SUBFIND); /* we have several groups in full to be done by the local CPU */
  }
  MPI_Barrier(Communicator);
  double ti1 = Logs.second();
  mpi_printf("SUBFIND: Processing overall took  (total time=%g sec)\n", Logs.timediff(ti0, ti1));

  /* free the communicator again */
  MPI_Comm_free(&SubComm);

#ifndef SUBFIND_HBT
  /* report some statistics */
  MPI_Allreduce(MPI_IN_PLACE, &count_decisions, 1, MPI_LONG_LONG, MPI_SUM, Communicator);
  MPI_Allreduce(MPI_IN_PLACE, &count_different_decisions, 1, MPI_LONG_LONG, MPI_SUM, Communicator);
  mpi_printf("SUBFIND: %lld out of %lld decisions, fraction %g, where influenced by previous subhalo size\n",
             count_different_decisions, count_decisions, ((double)count_different_decisions) / count_decisions);

#endif

  /* reestablish consistent global values for Sp.MaxPart/MaxPartSph in case they have diverged
   * in the subcommunicators
   */
  int max_load, max_loadsph;
  MPI_Allreduce(&Tp->MaxPart, &max_load, 1, MPI_INT, MPI_MAX, Communicator);
  MPI_Allreduce(&Tp->MaxPartSph, &max_loadsph, 1, MPI_INT, MPI_MAX, Communicator);

  /* do resize */
  Tp->reallocate_memory_maxpart(max_load);
  Tp->PS = (subfind_data *)Mem.myrealloc_movable(Tp->PS, Tp->MaxPart * sizeof(subfind_data));

  Tp->reallocate_memory_maxpartsph(max_loadsph);

  /* distribute particles back to original CPU */
  t0 = Logs.second();
  for(int i = 0; i < Tp->NumPart; i++)
    {
      Tp->PS[i].TargetTask  = Tp->PS[i].OriginTask;
      Tp->PS[i].TargetIndex = Tp->PS[i].OriginIndex;
    }

  FoFDomain->particle_exchange_based_on_PS(Communicator);
  t1 = Logs.second();
  if(ThisTask == 0)
    printf("SUBFIND: subfind_exchange() (for return to original CPU)  took %g sec\n", Logs.timediff(t0, t1));

  /* Now do he spherical overdensity calculations around group centers.
   * For this, we need a search tree with all particles again.
   */

  FoFGravTree.treeallocate(Tp->NumPart, Tp, FoFDomain);

  TIMER_STOP(CPU_SUBFIND);
  FoFGravTree.treebuild(Tp->NumPart, NULL);
  TIMER_START(CPU_SUBFIND);

  /* compute spherical overdensities for FOF groups */
  double cputimeso = subfind_overdensity();
  mpi_printf("SUBFIND: determining spherical overdensity masses took %g sec\n", cputimeso);

  FoFGravTree.treefree();

  /* sort the groups according to group/subhalo-number */
  t0 = Logs.second();
  mycxxsort_parallel(Group, Group + Ngroups, fof_compare_Group_GroupNr, Communicator);
  mycxxsort_parallel(Subhalo, Subhalo + Nsubhalos, subfind_compare_Subhalo_GroupNr_SubRankInGr, Communicator);
  t1 = Logs.second();
  mpi_printf("SUBFIND: assembled and ordered groups and subhalos (took %g sec)\n", Logs.timediff(t0, t1));

  sumup_large_ints(1, &Nsubhalos, &TotNsubhalos, Communicator);

  /* determine largest subhalo and total particle/cell count in substructures */
  long long lenmax = 0, glob_lenmax;
  long long totlen = 0;
  long long totsublength;
  for(int i = 0; i < Nsubhalos; i++)
    {
      totlen += Subhalo[i].Len;

      if(Subhalo[i].Len > lenmax)
        lenmax = Subhalo[i].Len;
    }
  sumup_longs(1, &totlen, &totsublength, Communicator);
  MPI_Reduce(&lenmax, &glob_lenmax, 1, MPI_LONG_LONG, MPI_MAX, 0, Communicator);

  /* set binding energy of unbound particles to zero (was overwritten with Hsml before) */
  for(int i = 0; i < Tp->NumPart; i++)
    if(Tp->PS[i].SubRankInGr == INT_MAX)
      Tp->PS[i].v.DM_BindingEnergy = 0;

  TIMER_START(CPU_SNAPSHOT);

  /* now do the output of the subhalo catalogue */
  subfind_save_final(num, basename, grpcat_dirbasename);

  TIMER_STOP(CPU_SNAPSHOT);

  double tend = Logs.second();

  if(ThisTask == 0)
    {
      printf("SUBFIND: Finished with SUBFIND.  (total time=%g sec)\n", Logs.timediff(tstart, tend));
      printf("SUBFIND: Total number of subhalos with at least %d particles: %lld\n", All.DesLinkNgb, TotNsubhalos);
      if(TotNsubhalos > 0)
        {
          printf("SUBFIND: Largest subhalo has %lld particles/cells.\n", glob_lenmax);
          printf("SUBFIND: Total number of particles/cells in subhalos: %lld\n", totsublength);
        }
    }

  Mem.myfree_movable(Subhalo);

  TIMER_STOP(CPU_SUBFIND);
}

template <typename partset>
void fof<partset>::subfind_save_final(int num, const char *basename, const char *grpcat_dirbasename)
{
  double t0 = Logs.second();

  long long totsubs = 0;

  /* fill in the FirstSub-values */
  for(int i = 0; i < Ngroups; i++)
    {
      if(i > 0)
        Group[i].FirstSub = Group[i - 1].FirstSub + Group[i - 1].Nsubs;
      else
        Group[i].FirstSub = 0;
      totsubs += Group[i].Nsubs;
    }

  long long *Send_count  = (long long *)Mem.mymalloc("Send_count", sizeof(long long) * this->NTask);
  long long *Send_offset = (long long *)Mem.mymalloc("Send_offset", sizeof(long long) * this->NTask);

  MPI_Allgather(&totsubs, 1, MPI_LONG_LONG, Send_count, 1, MPI_LONG_LONG, Communicator);

  Send_offset[0] = 0;

  for(int i = 1; i < NTask; i++)
    Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];

  for(int i = 0; i < Ngroups; i++)
    {
      if(Group[i].Nsubs > 0)
        Group[i].FirstSub += Send_offset[ThisTask];
      else
        Group[i].FirstSub = -1;
    }

  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  subfind_assign_subhalo_offsettype();

  fof_io<partset> FoF_IO{this, this->Communicator, All.SnapFormat};
  FoF_IO.fof_subfind_save_groups(num, basename, grpcat_dirbasename);

  double t1 = Logs.second();

  mpi_printf("SUBFIND: Subgroup catalogues saved. took = %g sec\n", Logs.timediff(t0, t1));
}

template <typename partset>
void fof<partset>::subfind_assign_subhalo_offsettype(void)
{
  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * this->NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * this->NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * this->NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * this->NTask);

  if(Nsubhalos > 0)
    {
      Subhalo[0].SubRankInGr = 0;
      for(int j = 0; j < NTYPES; j++)
        Subhalo[0].OffsetType[j] = 0;
    }

  for(int i = 1; i < Nsubhalos; i++)
    if(Subhalo[i].GroupNr != Subhalo[i - 1].GroupNr)
      {
        Subhalo[i].SubRankInGr = 0;
        for(int j = 0; j < NTYPES; j++)
          Subhalo[i].OffsetType[j] = 0;
      }
    else
      {
        Subhalo[i].SubRankInGr = Subhalo[i - 1].SubRankInGr + 1;
        for(int j = 0; j < NTYPES; j++)
          Subhalo[i].OffsetType[j] = Subhalo[i - 1].OffsetType[j] + Subhalo[i - 1].LenType[j];
      }

  struct subnr_info
  {
    int grnr;
    int cnt;
    long long stype_cum[NTYPES];
  };
  subnr_info subnr_data;

  if(Nsubhalos > 0)
    {
      subnr_data.grnr = Subhalo[Nsubhalos - 1].GroupNr;
      subnr_data.cnt  = Subhalo[Nsubhalos - 1].SubRankInGr + 1;
      for(int j = 0; j < NTYPES; j++)
        subnr_data.stype_cum[j] = Subhalo[Nsubhalos - 1].OffsetType[j] + Subhalo[Nsubhalos - 1].LenType[j];
    }
  else
    subnr_data.grnr = INT_MAX;

  subnr_info *info_all = (subnr_info *)Mem.mymalloc("info_all", NTask * sizeof(subnr_info));
  MPI_Allgather(&subnr_data, sizeof(subnr_info), MPI_BYTE, info_all, sizeof(subnr_info), MPI_BYTE, Communicator);

  if(Nsubhalos > 0)
    {
      int cnt = 0;
      long long stype_cum[NTYPES];
      for(int j = 0; j < NTYPES; j++)
        stype_cum[j] = 0;

      for(int i = ThisTask - 1; i >= 0; i--)
        if(info_all[i].grnr == Subhalo[0].GroupNr)
          {
            cnt += info_all[i].cnt;
            for(int j = 0; j < NTYPES; j++)
              stype_cum[j] += info_all[i].stype_cum[j];
          }

      for(int i = 0; i < Nsubhalos; i++)
        if(Subhalo[i].GroupNr == Subhalo[0].GroupNr)
          {
            Subhalo[i].SubRankInGr += cnt;
            for(int j = 0; j < NTYPES; j++)
              Subhalo[i].OffsetType[j] += stype_cum[j];
          }
        else
          break;
    }

  Mem.myfree(info_all);

  /* now need to send the subhalos to the processor that holds the parent group to inquire about the group
   * offset and then add this in.
   */

  /* tell everybody how many groups each processor holds */
  int *gcount = (int *)Mem.mymalloc("gcount", NTask * sizeof(int));
  MPI_Allgather(&Ngroups, 1, MPI_INT, gcount, 1, MPI_INT, Communicator);

  int nexport = 0, nimport = 0;

  struct group_info
  {
    int grindex;
    int subindex;
    long long OffsetType[NTYPES];
  };
  group_info *export_group_data = NULL, *import_group_data = NULL;

  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target                   = 0;
      long long first_gr_in_target = 0; /* this the first particle (of this type) on the current target processor */

      for(int i = 0; i < Nsubhalos; i++)
        {
          while(Subhalo[i].GroupNr >= first_gr_in_target + gcount[target])
            {
              if(target >= NTask)
                Terminate("target=%d  i=%d Nsubhalos=%d  Subhalo[i],GroupNr=%lld\n", target, i, Nsubhalos, Subhalo[i].GroupNr);

              first_gr_in_target += gcount[target];
              target++;
            }

          if(mode == 0)
            Send_count[target]++;
          else
            {
              int off                         = Send_offset[target] + Send_count[target]++;
              export_group_data[off].grindex  = Subhalo[i].GroupNr - first_gr_in_target;
              export_group_data[off].subindex = i;
            }
        }

      if(mode == 0)
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
          Recv_offset[0] = Send_offset[0] = 0;
          for(int j = 0; j < NTask; j++)
            {
              nimport += Recv_count[j];
              nexport += Send_count[j];
              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          export_group_data = (group_info *)Mem.mymalloc("export_group_data", nexport * sizeof(group_info));
          import_group_data = (group_info *)Mem.mymalloc("import_group_data", nimport * sizeof(group_info));
        }
    }

  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_group_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(group_info), MPI_BYTE, recvTask,
                       TAG_DENS_B, &import_group_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(group_info), MPI_BYTE,
                       recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
    }

  /* now let's go through the imported groups and assign the offsets */
  for(int i = 0; i < nimport; i++)
    for(int j = 0; j < NTYPES; j++)
      import_group_data[i].OffsetType[j] = Group[import_group_data[i].grindex].OffsetType[j];

  /* send stuff back */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&import_group_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(group_info), MPI_BYTE, recvTask,
                       TAG_DENS_B, &export_group_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(group_info), MPI_BYTE,
                       recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
    }

  /* add it in to subhalo offsets */
  for(int i = 0; i < nexport; i++)
    for(int j = 0; j < NTYPES; j++)
      Subhalo[export_group_data[i].subindex].OffsetType[j] += export_group_data[i].OffsetType[j];

  Mem.myfree(import_group_data);
  Mem.myfree(export_group_data);
  Mem.myfree(gcount);

  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);
}

template <typename partset>
void fof<partset>::subfind_redetermine_groupnr(void)
{
  /* tell everybody how many groups we have locally */
  int *ngroups_all = (int *)Mem.mymalloc("ngroups_all", NTask * sizeof(int));
  MPI_Allgather(&Ngroups, 1, MPI_INT, ngroups_all, 1, MPI_INT, Communicator);

  int nsubs_local = 0;
  for(int i = 0; i < Ngroups; i++)
    nsubs_local += Group[i].Nsubs;

  /* accumulate the corresponding subhalo numbers */
  int *nsubs_all = (int *)Mem.mymalloc("nsubs_all", NTask * sizeof(int));
  MPI_Allgather(&nsubs_local, 1, MPI_INT, nsubs_all, 1, MPI_INT, Communicator);

  /* Finally, also tell everybody how many subhalos we have locally */
  int *nsubhalos_all = (int *)Mem.mymalloc("nsubhalos_all", NTask * sizeof(int));
  MPI_Allgather(&Nsubhalos, 1, MPI_INT, nsubhalos_all, 1, MPI_INT, Communicator);

  long long subhalo_previous = 0;
  for(int i = 0; i < ThisTask; i++)
    subhalo_previous += nsubhalos_all[i];

  int nexport = 0, nimport = 0;

  struct group_info
  {
    long long subhalonr;
    long long groupnr;
    int subindex;
  };
  group_info *export_group_data = NULL, *import_group_data = NULL;

  int *Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * this->NTask);
  int *Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * this->NTask);
  int *Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * this->NTask);
  int *Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * this->NTask);

  for(int mode = 0; mode < 2; mode++)
    {
      for(int i = 0; i < NTask; i++)
        Send_count[i] = 0;

      int target               = 0;
      long long nsubs_previous = 0;

      for(int i = 0; i < Nsubhalos; i++)
        {
          long long subhalonr = subhalo_previous + i;

          while(subhalonr >= nsubs_previous + nsubs_all[target])
            {
              if(target >= NTask)
                Terminate("target=%d\n", target);

              nsubs_previous += nsubs_all[target];
              target++;
            }

          if(mode == 0)
            Send_count[target]++;
          else
            {
              int off                          = Send_offset[target] + Send_count[target]++;
              export_group_data[off].subhalonr = subhalonr;
              export_group_data[off].subindex  = i;
            }
        }

      if(mode == 0)
        {
          myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);
          Recv_offset[0] = Send_offset[0] = 0;
          for(int j = 0; j < NTask; j++)
            {
              nimport += Recv_count[j];
              nexport += Send_count[j];
              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          export_group_data = (group_info *)Mem.mymalloc("export_group_data", nexport * sizeof(group_info));
          import_group_data = (group_info *)Mem.mymalloc("import_group_data", nimport * sizeof(group_info));
        }
    }

  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&export_group_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(group_info), MPI_BYTE, recvTask,
                       TAG_DENS_B, &import_group_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(group_info), MPI_BYTE,
                       recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
    }

  /* now let's go through the imported subhalos and assign the group numbers */

  long long nsubs = 0;
  for(int i = 0; i < ThisTask; i++)
    nsubs += nsubs_all[i];

  long long group_previous = 0;
  for(int i = 0; i < ThisTask; i++)
    group_previous += ngroups_all[i];

  int gr = 0;
  for(int i = 0; i < nimport; i++)
    {
      while(import_group_data[i].subhalonr >= nsubs + Group[gr].Nsubs)
        {
          nsubs += Group[gr].Nsubs;
          gr++;

          if(gr >= Ngroups)
            Terminate("i=%d|%d gr=%d  >= Ngroups=%d\n", i, nimport, gr, Ngroups);
        }

      import_group_data[i].groupnr = group_previous + gr;
    }

  /* send stuff back */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          myMPI_Sendrecv(&import_group_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(group_info), MPI_BYTE, recvTask,
                       TAG_DENS_B, &export_group_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(group_info), MPI_BYTE,
                       recvTask, TAG_DENS_B, Communicator, MPI_STATUS_IGNORE);
    }

  /* assign the groupnr */
  for(int i = 0; i < nexport; i++)
    {
      /*
      if(Subhalo[export_group_data[i].subindex].GroupNr != export_group_data[i].groupnr)
        Terminate("Subhalo[export_group_data[i].subindex].GroupNr=%lld  export_group_data[i].groupnr=%lld\n",
                  Subhalo[export_group_data[i].subindex].GroupNr, export_group_data[i].groupnr);
      */
      Subhalo[export_group_data[i].subindex].GroupNr = export_group_data[i].groupnr;
    }

  Mem.myfree(import_group_data);
  Mem.myfree(export_group_data);
  Mem.myfree(Recv_offset);
  Mem.myfree(Recv_count);
  Mem.myfree(Send_offset);
  Mem.myfree(Send_count);

  Mem.myfree(nsubhalos_all);
  Mem.myfree(nsubs_all);
  Mem.myfree(ngroups_all);
}

template <typename partset>
double fof<partset>::subfind_get_particle_balance(void)
{
  int maxpart;
  long long sum;
  MPI_Allreduce(&Tp->NumPart, &maxpart, 1, MPI_INT, MPI_MAX, Communicator);
  sumup_large_ints(1, &Tp->NumPart, &sum, Communicator);
  return maxpart / (((double)sum) / NTask);
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
