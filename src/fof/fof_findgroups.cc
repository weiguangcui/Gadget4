/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file fof_findgroups.cc
 *
 *  \brief routines for identifying particle groups via friends-of-friends (FOF) linking
 */

#include "gadgetconfig.h"

#ifdef FOF

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../fof/fof.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/generic_comm.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
struct foffind_in : data_in_generic
{
#if defined(LIGHTCONE_PARTICLES_GROUPS)
  double DistanceOrigin;
#endif
  MyIntPosType IntPos[3];
  MyIDStorage MinID;
  int MinIDTask;
};

/* local data structure that holds results acquired on remote processors */
struct foffind_out
{
  char link_count_flag;
};

/* routine that fills the relevant particle/cell data into the input structure defined above */

template <typename T_tree, typename T_domain, typename T_partset>
class foffind_comm : public generic_comm<foffind_in, foffind_out, T_tree, T_domain, T_partset>
{
 public:
  typedef generic_comm<foffind_in, foffind_out, T_tree, T_domain, T_partset> gcomm;
  using gcomm::D;
  using gcomm::Thread;
  using gcomm::Tp;  // This makes sure that we can access Tp from the base class without having to use "this->Tp"
  using gcomm::Tree;

  /* need to call the base class constructor explicitly */
  foffind_comm(T_domain *dptr, T_tree *tptr, T_partset *pptr) : gcomm(dptr, tptr, pptr) {}

  void particle2in(foffind_in *in, int i) override
  {
    for(int k = 0; k < 3; k++)
      in->IntPos[k] = Tp->P[i].IntPos[k];

    in->MinID     = Tp->MinID[Tp->Head[i]];
    in->MinIDTask = Tp->MinIDTask[Tp->Head[i]];

#if defined(LIGHTCONE_PARTICLES_GROUPS)
    in->DistanceOrigin = Tp->DistanceOrigin[Tp->Head[i]];
#endif
  }

  void out2particle(foffind_out *out, int i, int mode) override
  {
    if(mode == MODE_LOCAL_PARTICLES) /* initial store */
      {
        /* nothing to be done here */
      }
    else /* combine */
      {
        if(out->link_count_flag)
          Tp->Flags[i].Marked = 1;
      }
  }

  int evaluate(int target, int mode, int thread_id, int action, foffind_in *in, int numnodes, node_info *firstnode,
               foffind_out &out) override
  {
    memset(&out, 0, sizeof(out));

#if defined(LIGHTCONE_PARTICLES_GROUPS)
    double target_DistanceOrigin = in->DistanceOrigin;
#else
    double target_DistanceOrigin = 0;
#endif

    int numngb = Tree->treefind_fof_primary(in->IntPos, Tp->LinkL, target, mode, &Thread, numnodes, firstnode, Thread.Ngblist,
                                            in->MinID, in->MinIDTask, target_DistanceOrigin);

    if(mode == MODE_IMPORTED_PARTICLES)
      {
        if(numngb > 0)
          out.link_count_flag = 1;
        else
          out.link_count_flag = 0;
      }

    return numngb;
  }
};

template <typename partset>
double fof<partset>::fof_find_groups(void)
{
  double tstart = Logs.second();

  mpi_printf("FOF: Start linking particles (presently allocated=%g MB)\n", Mem.getAllocatedBytesInMB());

  Tp->Flags = (typename partset::bit_flags *)Mem.mymalloc_clear("Flags", Tp->NumPart * sizeof(typename partset::bit_flags));

  FoFNgbTree.FullyLinkedNodePIndex = (int *)Mem.mymalloc("FullyLinkedNodePIndex", FoFNgbTree.NumNodes * sizeof(int));
  FoFNgbTree.FullyLinkedNodePIndex -= FoFNgbTree.MaxPart;

  for(int i = 0; i < FoFNgbTree.NumNodes; i++)
    {
      int no                               = i + FoFNgbTree.MaxPart;
      FoFNgbTree.FullyLinkedNodePIndex[no] = -1;
    }

  /* let's link all those primary particles which are in small enough nodes, only process local non-toplevel nodes */
  double taa = Logs.second();
  for(int i = 0; i < FoFNgbTree.NumNodes; i++)
    {
      int no = i + FoFNgbTree.MaxPart;

      if(FoFNgbTree.FullyLinkedNodePIndex[no] < 0)
        {
          if(FoFNgbTree.get_nodep(no)->level > 0)
            {
              double len = (((MyIntPosType)1) << (BITS_FOR_POSITIONS - FoFNgbTree.get_nodep(no)->level)) * Tp->FacIntToCoord;
              if(FACTSQRT3 * len < Tp->LinkL)  // all particles in this node can be linked
                {
                  int q = FoFNgbTree.treefind_fof_return_a_particle_in_cell_recursive(no);

                  if(q >= 0)
                    FoFNgbTree.fof_link_particles_in_cell_recursive(no, q);
                }
            }
        }
    }
  double tbb = Logs.second();
  mpi_printf("FOF: linking of small cells took %g sec\n", Logs.timediff(taa, tbb));

  /* first, link only among local particles */
  int *targetlist = (int *)Mem.mymalloc("targetlist", Tp->NumPart * sizeof(int));

  int npart = 0;
  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(is_type_primary_link_type(Tp->P[i].getType()))
        targetlist[npart++] = i;
    }

  /* create an object for handling the communication */
  foffind_comm<foftree<partset>, domain<partset>, partset> commpattern(FoFDomain, &FoFNgbTree, Tp);

  TIMER_STORE; /* active timer should be CPU_FOF */

  double t0 = Logs.second();

  commpattern.execute(npart, targetlist, MODE_LOCAL_NO_EXPORT, logs::CPU_FOFWALK, logs::CPU_FOFWALK, logs::CPU_FOFIMBAL);

  double t1 = Logs.second();

  int marked = 0;
  for(int i = 0; i < Tp->NumPart; i++)
    {
      if(is_type_primary_link_type(Tp->P[i].getType()))
        {
          if(Tp->Flags[i].Nonlocal)
            targetlist[marked++] = i;
        }
    }

  double dt = TIMER_DIFF(CPU_FOFWALK), dtmax, dtsum;
  MPI_Allreduce(&dt, &dtmax, 1, MPI_DOUBLE, MPI_MAX, Communicator);
  MPI_Allreduce(&dt, &dtsum, 1, MPI_DOUBLE, MPI_SUM, Communicator);

  long long totmarked, totnpart;
  sumup_large_ints(1, &marked, &totmarked, Communicator);
  sumup_large_ints(1, &npart, &totnpart, Communicator);
  mpi_printf(
      "FOF: local links done (took %g sec, avg-work=%g, imbalance=%g).\nFOF: Marked=%lld out of the %lld primaries "
      "which "
      "are linked\n",
      Logs.timediff(t0, t1), dtsum / NTask, dtmax / (dtsum / NTask), totmarked, totnpart);

  npart = marked;

  mpi_printf("FOF: begin linking across processors (presently allocated=%g MB) \n", Mem.getAllocatedBytesInMB());

  for(int i = 0; i < Tp->NumPart; i++)
    Tp->Flags[i].Marked = 1;

  long long link_across_tot;

  do
    {
      double t0 = Logs.second();

      npart = 0;
      for(int i = 0; i < Tp->NumPart; i++)
        {
          if(is_type_primary_link_type(Tp->P[i].getType()))
            {
              if(Tp->Flags[i].Nonlocal && Tp->Flags[i].Marked)
                targetlist[npart++] = i;
            }

          Tp->Flags[i].MinIDChanged = 0;
          Tp->Flags[i].Marked       = 0;
        }

      commpattern.execute(npart, targetlist, MODE_DEFAULT, logs::CPU_FOFWALK, logs::CPU_FOFWALK, logs::CPU_FOFIMBAL);

      int link_across = commpattern.Thread.Interactions;

      sumup_large_ints(1, &link_across, &link_across_tot, Communicator);

      long long ntot;
      sumup_large_ints(1, &npart, &ntot, Communicator);

      double t1 = Logs.second();

      mpi_printf("FOF: have done %15lld cross links (processed %14lld, took %g sec)\n", link_across_tot, ntot, Logs.timediff(t0, t1));

      /* let's check out which particles have changed their MinID */

      for(int i = 0; i < Tp->NumPart; i++)
        {
          if(Tp->Flags[i].Nonlocal)
            {
              if(Tp->Flags[Tp->Head[i]].MinIDChanged)
                Tp->Flags[i].Marked = 1;
            }
        }
    }
  while(link_across_tot > 0);

  Mem.myfree(targetlist);
  Mem.myfree(FoFNgbTree.FullyLinkedNodePIndex + FoFNgbTree.MaxPart);
  Mem.myfree(Tp->Flags);

  mpi_printf("FOF: Local groups found.\n");

  double tend = Logs.second();
  return Logs.timediff(tstart, tend);
}

/*! This function returns neighbors with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
template <typename partset>
int foftree<partset>::treefind_fof_primary(MyIntPosType *searchcenter, MyNgbTreeFloat hsml, int target, int mode, thread_data *thread,
                                           int numnodes, node_info *firstnode, int *ngblist, MyIDStorage target_MinID,
                                           int target_MinIDTask, double target_DistanceOrigin)
{
  if(mode == MODE_LOCAL_NO_EXPORT)
    Tp->Flags[target].Nonlocal = 0;

  MyNgbTreeFloat hsml2 = hsml * hsml;

  MyIntPosType search_min[3], search_range[3];
  MyIntPosType inthsml = hsml * Tp->FacCoordToInt;

  for(int i = 0; i < 3; i++)
    {
      search_min[i]   = searchcenter[i] - inthsml;
      search_range[i] = inthsml + inthsml;
    }

  int numngb = 0;

  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == MODE_LOCAL_PARTICLES || mode == MODE_LOCAL_NO_EXPORT)
        {
          no = MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k].Node;
          no = get_nodep(no)->nextnode; /* open it */
        }

      int shmrank = TreeSharedMem_ThisTask;

      while(no >= 0)
        {
          if(no < MaxPart) /* single particle */
            {
              if(shmrank != TreeSharedMem_ThisTask)
                Terminate("unexpected because in the present algorithm we are only allowed walk local branches");

              int p  = no;
              auto P = get_Pp(no, shmrank);

              no = get_nextnodep(shmrank)[no]; /* note: here shmrank cannot change */

              if(mode == MODE_LOCAL_PARTICLES) /* because we have already linked those in previous phase with MODE_LOCAL_NO_EXPORT */
                continue;

              double dx = ((MySignedIntPosType)(P->IntPos[0] - searchcenter[0])) * Tp->FacIntToCoord;
              double dd = dx * dx;
              if(dd > hsml2)
                continue;

              double dy = ((MySignedIntPosType)(P->IntPos[1] - searchcenter[1])) * Tp->FacIntToCoord;
              dd += dy * dy;
              if(dd > hsml2)
                continue;

              double dz = ((MySignedIntPosType)(P->IntPos[2] - searchcenter[2])) * Tp->FacIntToCoord;
              dd += dz * dz;
              if(dd > hsml2)
                continue;

              if(mode == MODE_IMPORTED_PARTICLES)
                {
#if defined(LIGHTCONE_PARTICLES_GROUPS)
                  if(Tp->DistanceOrigin[Tp->Head[p]] > target_DistanceOrigin)
                    {
                      Tp->DistanceOrigin[Tp->Head[p]]     = target_DistanceOrigin;
                      Tp->Flags[Tp->Head[p]].MinIDChanged = 1;
                      numngb++;
                    }
#endif
                  if(Tp->MinID[Tp->Head[p]].get() > target_MinID.get())
                    {
                      Tp->MinID[Tp->Head[p]]              = target_MinID;
                      Tp->MinIDTask[Tp->Head[p]]          = target_MinIDTask;
                      Tp->Flags[Tp->Head[p]].MinIDChanged = 1;
                      numngb++;
                    }
                }
              else if(mode == MODE_LOCAL_NO_EXPORT)
                {
                  if(Tp->Head[target] != Tp->Head[p])
                    {
                      int no_target = Father[target];
                      int no_p      = Father[p];

                      Tp->link_two_particles(target, p);

                      if(no_target >= 0)
                        {
                          if(FullyLinkedNodePIndex[no_target] >= 0)
                            {
                              if(Tp->Head[FullyLinkedNodePIndex[no_target]] != Tp->Head[target])
                                Terminate("how come that the node is fully linked, but not to us");
                            }
                          else
                            {
                              // this parent node was not yet fully linked to the final group: Is this the case now?
                              // If needed, check also all parent nodes
                              while(treefind_fof_check_single_node_for_full_linking(no_target))
                                no_target = get_nodep(no_target)->father;
                            }
                        }

                      if(no_p >= 0)
                        {
                          if(FullyLinkedNodePIndex[no_p] >= 0)
                            {
                              if(Tp->Head[FullyLinkedNodePIndex[no_p]] != Tp->Head[p])
                                Terminate("how come that the node is fully linked, but not to us");
                            }
                          else
                            {
                              // this parent node was not yet fully linked to the final group: Is this the case now?
                              // If needed, check also all parent nodes
                              while(treefind_fof_check_single_node_for_full_linking(no_p))
                                no_p = get_nodep(no_p)->father;
                            }
                        }
                    }
                }
              else
                Terminate("strange mode");
            }
          else if(no < MaxPart + MaxNodes) /* internal node */
            {
              fofnode *current = get_nodep(no, shmrank);

              if(current->level <= LEVEL_ALWAYS_OPEN)
                {
                  /* we always open the root node (its full node length couldn't be stored in the integer type */
                  no      = current->nextnode; /* no change in shmrank expected here */
                  shmrank = current->nextnode_shmrank;
                  continue;
                }
              /* check whether the node lies outside our search range */

              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;

                  if(FullyLinkedNodePIndex[no] >= 0)
                    {
                      int head = Tp->Head[FullyLinkedNodePIndex[no]];

                      if(head >= 0)
                        if(Tp->MinID[head].get() <= target_MinID.get())
                          {
#if defined(LIGHTCONE_PARTICLES_GROUPS)
                            if(Tp->DistanceOrigin[head] <= target_DistanceOrigin)
#endif
                              {
                                no      = current->sibling; /* the node can be discarded */
                                shmrank = current->sibling_shmrank;
                                continue;
                              }
                          }
                    }
                }
              else if(mode == MODE_LOCAL_PARTICLES)
                {
                  int p = current->nextnode;

                  /* in case the next node after opening is not a top-level node, we have either reached a leaf node or are in a local
                   * branch we need to do nothing if we would end up on different shared memory thread */
                  if(p < MaxPart || (p >= FirstNonTopLevelNode && p < MaxPart + MaxNodes))
                    {
                      if(current->nextnode_shmrank != TreeSharedMem_ThisTask)
                        {
                          int task = D->ThisTask + current->nextnode_shmrank - TreeSharedMem_ThisTask;

                          if(target >= 0) /* export */
                            tree_export_node_threads_by_task_and_node(task, no, target, thread);

                          no      = current->sibling; /* in case the node can be discarded */
                          shmrank = current->sibling_shmrank;
                          continue;
                        }
                    }

                  if(no >= FirstNonTopLevelNode)
                    {
                      /* we have a node with only local particles, hence we can skip it for mode == 0 */
                      no      = current->sibling; /* in case the node can be discarded */
                      shmrank = current->sibling_shmrank;
                      continue;
                    }
                }
              else if(mode == MODE_LOCAL_NO_EXPORT)
                {
                  int p = current->nextnode;

                  /* in case the next node after opening is not a top-level node, we have either reached a leaf node or are in a local
                   * branch we need to do nothing if we would end up on different shared memory thread */
                  if(p < MaxPart || (p >= FirstNonTopLevelNode && p < MaxPart + MaxNodes))
                    {
                      if(current->nextnode_shmrank != TreeSharedMem_ThisTask)
                        {
                          no      = current->sibling; /* in case the node can be discarded */
                          shmrank = current->sibling_shmrank;

                          MyIntPosType left[3], right[3];

                          left[0]  = current->range_min[0] - search_min[0];
                          right[0] = current->range_max[0] - search_min[0];

                          /* check whether we can stop walking along this branch */
                          if(left[0] > search_range[0] && right[0] > left[0])
                            continue;

                          left[1]  = current->range_min[1] - search_min[1];
                          right[1] = current->range_max[1] - search_min[1];

                          /* check whether we can stop walking along this branch */
                          if(left[1] > search_range[1] && right[1] > left[1])
                            continue;

                          left[2]  = current->range_min[2] - search_min[2];
                          right[2] = current->range_max[2] - search_min[2];

                          /* check whether we can stop walking along this branch */
                          if(left[2] > search_range[2] && right[2] > left[2])
                            continue;

                          Tp->Flags[target].Nonlocal = 1;
                          continue;
                        }
                    }

                  if(FullyLinkedNodePIndex[no] >= 0)
                    if(Tp->Head[target] == Tp->Head[FullyLinkedNodePIndex[no]])  // all particles in the node are linked to us anyhow
                      {
                        no      = current->sibling; /* in case the node can be discarded */
                        shmrank = current->sibling_shmrank;
                        continue;
                      }
                }

              MyIntPosType left[3], right[3];

              left[0]  = current->range_min[0] - search_min[0];
              right[0] = current->range_max[0] - search_min[0];

              /* check whether we can stop walking along this branch */
              if(left[0] > search_range[0] && right[0] > left[0])
                {
                  no      = current->sibling; /* in case the node can be discarded */
                  shmrank = current->sibling_shmrank;
                  continue;
                }

              left[1]  = current->range_min[1] - search_min[1];
              right[1] = current->range_max[1] - search_min[1];

              /* check whether we can stop walking along this branch */
              if(left[1] > search_range[1] && right[1] > left[1])
                {
                  no      = current->sibling; /* in case the node can be discarded */
                  shmrank = current->sibling_shmrank;
                  continue;
                }

              left[2]  = current->range_min[2] - search_min[2];
              right[2] = current->range_max[2] - search_min[2];

              /* check whether we can stop walking along this branch */
              if(left[2] > search_range[2] && right[2] > left[2])
                {
                  no      = current->sibling; /* in case the node can be discarded */
                  shmrank = current->sibling_shmrank;
                  continue;
                }

              no      = current->nextnode;         /* ok, we need to open the node */
              shmrank = current->nextnode_shmrank; /* ok, we need to open the node */
            }
          else /* pseudo particle */
            {
              if(mode == MODE_LOCAL_PARTICLES)
                {
                  if(target >= 0) /* if no target is given, export will not occur */
                    tree_export_node_threads(no, target, thread);
                }
              else if(mode == MODE_LOCAL_NO_EXPORT)
                {
                  Tp->Flags[target].Nonlocal = 1;
                }

              no = get_nextnodep(shmrank)[no - MaxNodes];
              /* note: here shmrank does not need to change */

              continue;
            }
        }
    }

  return numngb;
}

template <typename partset>
void foftree<partset>::fof_link_particles_in_cell_recursive(int no, int q)
{
  if(no >= MaxPart && no < MaxPart + MaxNodes) /* internal node */
    {
      FullyLinkedNodePIndex[no] = q;

      int p = get_nodep(no)->nextnode;

      /* in case the next node after opening is not a top-level node, we have either reached a leaf node or are in a local
       * branch. We need to do nothing if we would end up on different shared memory thread */
      if(p < MaxPart || (p >= FirstNonTopLevelNode && p < MaxPart + MaxNodes))
        {
          if(get_nodep(no)->nextnode_shmrank != TreeSharedMem_ThisTask)
            return;
        }

      while(p != get_nodep(no)->sibling)
        {
          if(p < MaxPart) /* a particle */
            {
              if(p != q)
                {
                  Tp->link_two_particles(p, q);  // link them if not already linked
                }

              p = Nextnode[p];
            }
          else if(p < MaxPart + MaxNodes) /* an internal node  */
            {
              fof_link_particles_in_cell_recursive(p, q);

              p = get_nodep(p)->sibling;
            }
          else /* a pseudo particle */
            p = Nextnode[p - MaxNodes];
        }
    }
}

template <typename partset>
int foftree<partset>::treefind_fof_return_a_particle_in_cell_recursive(int no)
{
  if(no >= MaxPart && no < MaxPart + MaxNodes) /* internal node */
    {
      int p = get_nodep(no)->nextnode;

      /* in case the next node after opening is not a top-level node, we have either reached a leaf node or are in a local
       * branch. We need to do nothing if we would end up on different shared memory thread */
      if(p < MaxPart || (p >= FirstNonTopLevelNode && p < MaxPart + MaxNodes))
        {
          if(get_nodep(no)->nextnode_shmrank != TreeSharedMem_ThisTask)
            return -1;
        }

      while(p != get_nodep(no)->sibling)
        {
          if(p < MaxPart) /* a particle */
            {
              return p;

              p = Nextnode[p];
            }
          else if(p < MaxPart + MaxNodes) /* an internal node  */
            {
              int ret = treefind_fof_return_a_particle_in_cell_recursive(p);

              if(ret >= 0)
                return ret;

              p = get_nodep(p)->sibling;
            }
          else /* a pseudo particle */
            p = Nextnode[p - MaxNodes];
        }
    }

  return -1;
}

template <typename partset>
int foftree<partset>::treefind_fof_check_single_node_for_full_linking(int no)
{
  if(no >= MaxPart && no < MaxPart + MaxNodes) /* internal node */
    {
      if(FullyLinkedNodePIndex[no] >= 0)  // already linked
        return 0;

      int head = -1; /* no particle yet */

      int p = get_nodep(no)->nextnode;

      /* in case the next node after opening is not a top-level node, we have either reached a leaf node or are in a local
       * branch. We need to do nothing if we would end up on different shared memory thread */
      if(p < MaxPart || (p >= FirstNonTopLevelNode && p < MaxPart + MaxNodes))
        {
          if(get_nodep(no)->nextnode_shmrank != TreeSharedMem_ThisTask)
            return 0;
        }

      while(p != get_nodep(no)->sibling)
        {
          if(p < MaxPart) /* a particle */
            {
              if(head == -1)
                head = Tp->Head[p];
              else if(head >= 0)
                {
                  if(head != Tp->Head[p])
                    {
                      head = -2;
                      break;
                    }
                }

              p = Nextnode[p];
            }
          else if(p < MaxPart + MaxNodes) /* an internal node  */
            {
              if(FullyLinkedNodePIndex[p] >= 0)
                {
                  if(head == -1)
                    head = Tp->Head[FullyLinkedNodePIndex[p]];
                  else if(head >= 0)
                    {
                      if(head != Tp->Head[FullyLinkedNodePIndex[p]])
                        {
                          head = -2;
                          break;
                        }
                    }
                }
              else
                {
                  if(treefind_fof_return_a_particle_in_cell_recursive(no) >= 0)
                    {
                      head = -2;
                      break;
                    }
                }

              p = get_nodep(p)->sibling;
            }
          else /* a pseudo particle */
            p = Nextnode[p - MaxNodes];
        }

      if(head >= 0)
        {
          FullyLinkedNodePIndex[no] = treefind_fof_return_a_particle_in_cell_recursive(no);

          if(Tp->Head[FullyLinkedNodePIndex[no]] != head)
            Terminate("no=%d  FullyLinkedNodePIndex[no]=%d   Tp->Head[FullyLinkedNodePIndex[no]]=%d   head=%d \n", no,
                      FullyLinkedNodePIndex[no], Tp->Head[FullyLinkedNodePIndex[no]], head);

          return 1;
        }
    }

  return 0;
}

/* now make sure that the following classes are really instantiated, otherwise we may get a linking problem */
#include "../data/simparticles.h"
template class fof<simparticles>;

#if defined(LIGHTCONE) && defined(LIGHTCONE_PARTICLES_GROUPS)
#include "../data/lcparticles.h"
template class fof<lcparticles>;
#endif

#endif
