/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  generic_comm.h
 *
 *  \brief defines a class with a generic communication pattern for parallel tree walks
 */

#ifndef GENERIC_COMM_H
#define GENERIC_COMM_H

#include "gadgetconfig.h"

#include "../domain/domain.h"
#include "../logs/logs.h"
#include "../mpi_utils/mpi_utils.h"

#define EXTRA_SPACE 16384

struct data_in_generic
{
  int Firstnode;
};

template <typename T_in, typename T_out, typename T_tree, typename T_domain, typename T_partset>
class generic_comm
{
 public:
  virtual void particle2in(T_in *in, int i)              = 0;  // Pure virtual function, *must* be overridden
  virtual void out2particle(T_out *out, int i, int mode) = 0;  // Pure virtual function, *must* be overridden
  virtual int evaluate(int target, int mode, int thread_id, int action, T_in *in, int numnodes, node_info *firstnode,
                       T_out &out)                       = 0;  // Pure virtual function, *must* be overridden

  T_domain *D;
  T_tree *Tree;
  T_partset *Tp;

  thread_data Thread;

  generic_comm(T_domain *dptr, T_tree *tptr, T_partset *pptr)
  {
    D    = dptr;
    Tree = tptr;
    Tp   = pptr;

    set_MaxNexport(); /* initialisation */
  }

  /* some public diagnostic info */

  long long SumNexport;

 private:
  enum logs::timers cpu_primary   = logs::CPU_NONE;
  enum logs::timers cpu_secondary = logs::CPU_NONE;
  enum logs::timers cpu_imbalance = logs::CPU_NONE;

  int Nactive;
  int *Targets;
  int actioncode;

  T_in *DataIn, *DataGet;
  T_out *DataResult, *DataOut;

  struct send_recv_counts
  {
    int Count;
    int CountNodes;
  };
  send_recv_counts *Send, *Recv;

  /** Array of NTask size of the offset into the send array where the
      objects to be sent to the specified task starts. */
  int *Send_offset,
      /** Array of NTask size of the number of objects to send to the
          tasks. */
      *Send_count,
      /** Array of NTask size of the number of objects to receive from the
          tasks. */
      *Recv_count,
      /** Array of NTask size of the offset into the receive array where the
          objects from the specified task starts. */
      *Recv_offset;

  int *Send_offset_nodes, *Send_count_nodes, *Recv_count_nodes, *Recv_offset_nodes;

  /** Buffer of size NTask used for flagging whether a particle needs to
    be exported to the other tasks. */
  int *Exportflag;
  /** Buffer of size NTask used for counting how many nodes are to be
      exported to the other tasks? */
  int *Exportnodecount;
  /** Buffer of size NTask used for holding the index into the
      DataIndexTable. */
  int *Exportindex;

  size_t ExportSpace;
  size_t MinSpace;
  int NextParticle;
  int Nexport, Nimport;
  int NexportNodes, NimportNodes;

  node_info *NodeInfoIn;
  node_info *NodeInfoGet;

  char callorigin[1000];

  void primary_loop(int mode)
  {
    if(cpu_primary != logs::CPU_NONE)
      Logs.timer_start(cpu_primary);

    int idx;

    int j;

    for(j = 0; j < D->NTask; j++)
      Thread.Exportflag[j] = -1;

    while(1)
      {
        if(Thread.ExportSpace < MinSpace)
          break;

        idx = NextParticle++;

        if(idx >= Nactive)
          break;

        int i = Targets[idx];
        if(i < 0)
          continue;

        T_in local;
        T_out out;
        particle2in(&local, i);
        local.Firstnode = 0;

        int num = evaluate(i, mode, 0, actioncode, &local, 1, NULL, out);

        out2particle(&out, i, mode);

        Thread.Interactions += num;
      }

    if(cpu_primary != logs::CPU_NONE)
      Logs.timer_stop(cpu_primary);
  }

  void secondary_loop(void)
  {
    if(cpu_secondary != logs::CPU_NONE)
      Logs.timer_start(cpu_secondary);

    /* now do the particles that were sent to us */
    int i, cnt = 0;

    {
      while(1)
        {
          i = cnt++;

          if(i >= Nimport)
            break;

          int numnodes;
          node_info *firstnode;
          generic_get_numnodes(i, &numnodes, &firstnode);

          T_in *in = &DataGet[i];
          T_out out;

          int num = evaluate(i, MODE_IMPORTED_PARTICLES, 0, actioncode, in, numnodes, firstnode, out);

          DataResult[i] = out;

          Thread.Interactions += num;
        }
    }

    if(cpu_secondary != logs::CPU_NONE)
      Logs.timer_stop(cpu_secondary);
  }

  /* this function determines how much buffer space we may use based on the memory that is locally still free,
   * and it computes how much memory may at most be needed to process a single particle. We will only continue with a particle
   * if this can still be safely processed.
   */

  /* this function does the memory allocation at the beginning of a loop over the remaining local particles.
   * The fields PartList[] and NodeList[] share the buffer space of size "ExportSpace" (in bytes).
   * Here PartList will be filled in from the beginning, while NodeList will be filled in from the end.
   * Since we do not know a priory the relative share of these two fields, we can make optimum use of
   * the available space in this way.
   */

  void set_MaxNexport(const char *file = __FILE__, int line = __LINE__)
  {
    ExportSpace = 0.5 * (Mem.FreeBytes); /* we just grab at most half of the still available memory here */

    if(ExportSpace <= Tp->NumPart * sizeof(int))
      {
        Mem.dump_memory_table();
        Terminate("It seems we have too little space left for properly sized ExportSpace... (%lld %lld)   Need more memory.\n",
                  (long long)ExportSpace, (long long)Tp->NumPart * sizeof(int));
      }

    ExportSpace -= Tp->NumPart * sizeof(int); /* to account for the neighbor list buffer that the process allocated */

    /* make the size a multiple both of data_partlist and data_nodelist */
    ExportSpace /= (sizeof(data_partlist) * sizeof(data_nodelist));
    ExportSpace *= (sizeof(data_partlist) * sizeof(data_nodelist));

    MinSpace = (D->NTask - 1) * (sizeof(data_partlist) + sizeof(T_in) + sizeof(T_out)) +
               D->NTopleaves * (sizeof(data_nodelist) + sizeof(int));

    snprintf(callorigin, MAXLEN_PATH_EXTRA, "%s|%d|", file, line);

    if(ExportSpace < MinSpace)
      {
        Mem.dump_memory_table();
        Terminate(
            "Bummer. Can't even safely process a single particle for the available memory. FreeBytes=%lld  ExportSpace=%lld  "
            "MinSpace=%lld  D->NTask=%d  NTopleaves=%d",
            (long long)Mem.FreeBytes, (long long)ExportSpace, (long long)MinSpace, D->NTask, D->NTopleaves);
      }
  }

  /* this function does the memory allocation at the beginning of a loop over the remaining local particles.
   * The fields PartList[] and NodeList[] share the buffer space of size "ExportSpace" (in bytes).
   * Here PartList will be filled in from the beginning, while NodeList will be filled in from the end.
   * Since we do not know a priory the relative share of these two fields, we can make optimum use of
   * the available space in this way.
   */
  void generic_alloc_partlist_nodelist_ngblist_threadbufs(void)
  {
    Thread.Nexport      = 0;
    Thread.NexportNodes = 0;
    Thread.ExportSpace  = ExportSpace;
    Thread.InitialSpace = ExportSpace;
    Thread.ItemSize     = (sizeof(data_partlist) + sizeof(T_in) + sizeof(T_out));

    Thread.PartList = (data_partlist *)Mem.mymalloc_movable_g(&Thread.PartList, "PartList", ExportSpace);
    /* note: the NodeList array will be attached to the end of this buffer, growing backwards */

    Thread.Ngblist     = (int *)Mem.mymalloc_movable_g(&Thread.Ngblist, "Ngblist", Tp->NumPart * sizeof(int));
    Thread.Shmranklist = (int *)Mem.mymalloc_movable_g(&Thread.Shmranklist, "Shmranklist", Tp->NumPart * sizeof(int));
    Thread.Exportflag  = Exportflag;
  }

  void generic_free_partlist_nodelist_ngblist_threadbufs(void)
  {
    Mem.myfree(Thread.Shmranklist);
    Mem.myfree(Thread.Ngblist);
    Mem.myfree(Thread.PartList);
    Thread.Shmranklist = NULL;
    Thread.Ngblist     = NULL;
    Thread.PartList    = NULL;
  }

  void generic_prepare_export_counts(void)
  {
    for(int j = 0; j < D->NTask; j++)
      {
        Send[j].Count      = 0;
        Send[j].CountNodes = 0;
      }

    Nexport      = 0;
    NexportNodes = 0;

    for(int j = 0; j < Thread.Nexport; j++)
      Send[Thread.PartList[j].Task].Count++;

    data_nodelist *nodelist = (data_nodelist *)(((char *)Thread.PartList) + Thread.InitialSpace);

    for(int j = 0; j < Thread.NexportNodes; j++)
      Send[nodelist[-1 - j].Task].CountNodes++;

    Nexport += Thread.Nexport;
    NexportNodes += Thread.NexportNodes;

    SumNexport += Nexport;
  }

  /* establish the Recv counts from the Send counts (effectively a big transpose)
   */
  void generic_prepare_import_counts(void)
  {
    /* our standard approach for this is to use an all-to-all communication. For very large processor counts,
     * this in principle becomes inefficient since mostly zeros need to be communicated.
     * we have also two option experimental communication routines that use a sparse=communication pattern instead.
     */
    /* the default */
    myMPI_Alltoall(Send, sizeof(send_recv_counts), MPI_BYTE, Recv, sizeof(send_recv_counts), MPI_BYTE, D->Communicator);
  }

  /* initialize offset tables that we need for the communication
   */
  void generic_prepare_export_offsets(void)
  {
    Send_offset[0]       = 0;
    Send_offset_nodes[0] = 0;

    for(int j = 1; j < D->NTask; j++)
      {
        Send_offset[j]       = Send_offset[j - 1] + Send[j - 1].Count;
        Send_offset_nodes[j] = Send_offset_nodes[j - 1] + Send[j - 1].CountNodes;
      }
  }

  /* organize the particle and node data for export in contiguous memory regions
   */
  void generic_prepare_particle_data_for_export(void)
  {
    int *rel_node_index = (int *)Mem.mymalloc_g("rel_node_index", D->NTask * sizeof(int));

    for(int j = 0; j < D->NTask; j++)
      {
        Send[j].Count      = 0;
        Send[j].CountNodes = 0;
        rel_node_index[j]  = 0;
      }

    data_nodelist *nodelist = (data_nodelist *)(((char *)Thread.PartList) + Thread.InitialSpace);

    for(int j = 0, jj = 0; j < Thread.Nexport; j++)
      {
        int task = Thread.PartList[j].Task;
        int off  = Send_offset[task] + Send[task].Count++;

        int target = Thread.PartList[j].Index;

        particle2in(&DataIn[off], target);
        DataIn[off].Firstnode = rel_node_index[task];

        if(j < Thread.Nexport - 1)
          if(Thread.PartList[j].Index == Thread.PartList[j + 1].Index)
            continue;

        while(jj < Thread.NexportNodes && Thread.PartList[j].Index == nodelist[-1 - jj].Index)
          {
            int task = nodelist[-1 - jj].Task;
            int off  = Send_offset_nodes[task] + Send[task].CountNodes++;

            NodeInfoIn[off] = nodelist[-1 - jj].NodeInfo;

            rel_node_index[task]++;
            jj++;
          }
      }

    Mem.myfree(rel_node_index);
  }

  /* driver routine to process the results that we obtained for a particle from a remote processor
   * by working on it with the supplied out2particle() routine
   */
  void generic_add_results_to_local(void)
  {
    for(int j = 0; j < D->NTask; j++)
      Send[j].Count = 0;

    for(int j = 0; j < Thread.Nexport; j++)
      {
        int task = Thread.PartList[j].Task;
        int off  = Send_offset[task] + Send[task].Count++;

        int target = Thread.PartList[j].Index;

        out2particle(&DataOut[off], target, MODE_IMPORTED_PARTICLES);
      }
  }

#ifndef OMIT_GENERIC_GET_NUMNODES

  /* this function is called in the actual tree walk routine to find out how the number and
   * starting index of the section in the node-list that needs to be processed for the imported particle
   */
  void generic_get_numnodes(int target, int *numnodes, node_info **firstnode)
  {
    if(target == Nimport - 1)
      *numnodes = NimportNodes - DataGet[target].Firstnode;
    else
      *numnodes = DataGet[target + 1].Firstnode - DataGet[target].Firstnode;

    if(*numnodes > Tree->NumNodes)
      {
        Terminate(
            "odd: target=%d  Nimport=%d  NimportNodes=%d numnodes=%d Tree->NumNodes=%d Tree->MaxNodes=%d "
            "DataGet[target].Firstnode=%d\n",
            target, Nimport, NimportNodes, *numnodes, Tree->NumNodes, Tree->MaxNodes, DataGet[target].Firstnode);
      }

    *firstnode = &NodeInfoGet[DataGet[target].Firstnode];
  }

#endif

  /* calculate how many space we need to allocate to safely process a certain number of
   * nodes and particles that are imported.
   */
  size_t generic_calc_import_storage(int nimport, int nimportnodes)
  {
    size_t needed = nimport * sizeof(T_in) + nimportnodes * sizeof(int) + nimport * sizeof(T_out);

    /* add some extra space to not go to the last byte */
    needed += EXTRA_SPACE;

    return needed;
  }

  /* this routine carries out the communication step in several phases if needed
   */
  void generic_multiple_phases(void)
  {
    int ncycles;

    for(int ngrpstart = 1; ngrpstart < (1 << D->PTask); ngrpstart += ncycles)
      {
        /* now decide how many cycles we can process in this iteration */
        ncycles = (1 << D->PTask) - ngrpstart;

        do
          {
            Nimport      = 0;
            NimportNodes = 0;

            for(int ngrp = ngrpstart; ngrp < ngrpstart + ncycles; ngrp++)
              {
                int recvTask = D->ThisTask ^ ngrp;

                if(recvTask < D->NTask)
                  {
                    if(Recv[recvTask].Count > 0)
                      {
                        Nimport += Recv[recvTask].Count;
                        NimportNodes += Recv[recvTask].CountNodes;
                      }
                  }
              }

            int flag = 0, flagall;

            if(generic_calc_import_storage(Nimport, NimportNodes) > Mem.FreeBytes)
              flag = 1;

            MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, D->Communicator);

            if(flagall)
              ncycles /= 2;
            else
              break;
          }
        while(ncycles > 0);

        if(ncycles == 0)
          Terminate(
              "Seems like we can't even do one cycle: ncycles=%d  ngrpstart=%d  Nimport=%d  NimportNodes=%d  FreeBytes=%lld  needed "
              "storage=%lld",
              ncycles, ngrpstart, Nimport, NimportNodes, (long long)Mem.FreeBytes,
              (long long)generic_calc_import_storage(Nimport, NimportNodes));

        if(ngrpstart == 1 && ncycles != ((1 << D->PTask) - ngrpstart) && D->ThisTask == 0)
          warn("need multiple import/export phases to avoid memory overflow");

        /* now allocated the import and results buffers */

        DataGet     = (T_in *)Mem.mymalloc_movable_g(&DataGet, "DataGet", Nimport * sizeof(T_in));
        NodeInfoGet = (node_info *)Mem.mymalloc_movable_g(&NodeInfoGet, "NodeInfoGet", NimportNodes * sizeof(node_info));
        DataResult  = (T_out *)Mem.mymalloc_movable_g(&DataResult, "DataResult", Nimport * sizeof(T_out));

        Nimport      = 0;
        NimportNodes = 0;

        /* exchange particle data */
        for(int ngrp = ngrpstart; ngrp < ngrpstart + ncycles; ngrp++)
          {
            int recvTask = D->ThisTask ^ ngrp;

            if(recvTask < D->NTask)
              {
                if(Send[recvTask].Count > 0 || Recv[recvTask].Count > 0)
                  {
                    size_t len = sizeof(T_in);

                    /* get the particles */
                    myMPI_Sendrecv(&DataIn[Send_offset[recvTask]], Send[recvTask].Count * len, MPI_BYTE, recvTask, TAG_HYDRO_A,
                                   &DataGet[Nimport], Recv[recvTask].Count * len, MPI_BYTE, recvTask, TAG_HYDRO_A, D->Communicator,
                                   MPI_STATUS_IGNORE);

                    /* get the node info */
                    myMPI_Sendrecv(&NodeInfoIn[Send_offset_nodes[recvTask]], Send[recvTask].CountNodes * sizeof(node_info), MPI_BYTE,
                                   recvTask, TAG_GRAV_B, &NodeInfoGet[NimportNodes], Recv[recvTask].CountNodes * sizeof(node_info),
                                   MPI_BYTE, recvTask, TAG_GRAV_B, D->Communicator, MPI_STATUS_IGNORE);

                    for(int k = 0; k < Recv[recvTask].Count; k++)
                      DataGet[Nimport + k].Firstnode += NimportNodes;

                    Nimport += Recv[recvTask].Count;
                    NimportNodes += Recv[recvTask].CountNodes;
                  }
              }
          }

        /* now do the actual work for the imported points */
        secondary_loop();

        /* send the results */
        Nimport      = 0;
        NimportNodes = 0;

        for(int ngrp = ngrpstart; ngrp < ngrpstart + ncycles; ngrp++)
          {
            int recvTask = D->ThisTask ^ ngrp;
            if(recvTask < D->NTask)
              {
                if(Send[recvTask].Count > 0 || Recv[recvTask].Count > 0)
                  {
                    size_t len = sizeof(T_out);

                    /* exchange the results */
                    myMPI_Sendrecv(&DataResult[Nimport], Recv[recvTask].Count * len, MPI_BYTE, recvTask, TAG_HYDRO_B,
                                   &DataOut[Send_offset[recvTask]], Send[recvTask].Count * len, MPI_BYTE, recvTask, TAG_HYDRO_B,
                                   D->Communicator, MPI_STATUS_IGNORE);

                    Nimport += Recv[recvTask].Count;
                    NimportNodes += Recv[recvTask].CountNodes;
                  }
              }
          }

        Mem.myfree(DataResult);
        Mem.myfree(NodeInfoGet);
        Mem.myfree(DataGet);
      }
  }

  /* this function deals with the communication step, and then processes the imported particles, and finally computes the results back.
   * if there is not enough memory available to hold all the data sent to us from other processors, we process the incoming data in
   * multiple stages, which will always be possible.
   */
  void generic_exchange(void)
  {
    /* set up Sendcount table */
    generic_prepare_export_counts();

    /* do the all-to-all exchange so that we have the Recvcount table as well */
    generic_prepare_import_counts();

    /* prepare offsets in export tables */
    generic_prepare_export_offsets();

    /* allocate particle data buffers */
    DataIn     = (T_in *)Mem.mymalloc_movable_g(&DataIn, "DataIn", Nexport * sizeof(T_in));
    NodeInfoIn = (node_info *)Mem.mymalloc_movable_g(&NodeInfoIn, "NodeInfoIn", NexportNodes * sizeof(node_info));
    DataOut    = (T_out *)Mem.mymalloc_movable_g(&DataOut, "DataOut", Nexport * sizeof(T_out));

    /* prepare particle data for export */
    generic_prepare_particle_data_for_export();

    /* export particles and process them, if needed in several installments */
    generic_multiple_phases();

    /* add the results to the local particles */
    generic_add_results_to_local();

    Mem.myfree(DataOut);
    Mem.myfree(NodeInfoIn);
    Mem.myfree(DataIn);
  }

  void generic_allocate_comm_tables(void)
  {
    ptrdiff_t off;

    off = (char *)Tree->Nodes - Mem.Base;
    MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreeNodes_offsets, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreeSharedMemComm);

    off = (char *)Tree->Points - Mem.Base;
    MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreePoints_offsets, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreeSharedMemComm);

    off = (char *)Tree->Nextnode - Mem.Base;
    MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreeNextnode_offsets, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreeSharedMemComm);

    off = (char *)Tp->P - Mem.Base;
    MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreeP_offsets, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreeSharedMemComm);

    off = (char *)Tp->PS - Mem.Base;
    MPI_Allgather(&off, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreePS_offsets, sizeof(ptrdiff_t), MPI_BYTE, Tree->TreeSharedMemComm);

    Exportflag      = (int *)Mem.mymalloc("Exportflag", ((((D->NTask - 1) / 16) + 1) * 16) * sizeof(int));
    Exportindex     = (int *)Mem.mymalloc("Exportindex", D->NTask * sizeof(int));
    Exportnodecount = (int *)Mem.mymalloc("Exportnodecount", D->NTask * sizeof(int));

    Send = (send_recv_counts *)Mem.mymalloc("Send", sizeof(send_recv_counts) * D->NTask);
    Recv = (send_recv_counts *)Mem.mymalloc("Recv", sizeof(send_recv_counts) * D->NTask);

    Send_count  = (int *)Mem.mymalloc("Send_count", sizeof(int) * D->NTask);
    Send_offset = (int *)Mem.mymalloc("Send_offset", sizeof(int) * D->NTask);
    Recv_count  = (int *)Mem.mymalloc("Recv_count", sizeof(int) * D->NTask);
    Recv_offset = (int *)Mem.mymalloc("Recv_offset", sizeof(int) * D->NTask);

    Send_count_nodes  = (int *)Mem.mymalloc("Send_count_nodes", sizeof(int) * D->NTask);
    Send_offset_nodes = (int *)Mem.mymalloc("Send_offset_nodes", sizeof(int) * D->NTask);
    Recv_count_nodes  = (int *)Mem.mymalloc("Recv_count_nodes", sizeof(int) * D->NTask);
    Recv_offset_nodes = (int *)Mem.mymalloc("Recv_offset_nodes", sizeof(int) * D->NTask);
  }

  void generic_free_comm_tables(void)
  {
    Mem.myfree(Recv_offset_nodes);
    Mem.myfree(Recv_count_nodes);
    Mem.myfree(Send_offset_nodes);
    Mem.myfree(Send_count_nodes);
    Mem.myfree(Recv_offset);
    Mem.myfree(Recv_count);
    Mem.myfree(Send_offset);
    Mem.myfree(Send_count);
    Mem.myfree(Recv);
    Mem.myfree(Send);
    Mem.myfree(Exportnodecount);
    Mem.myfree(Exportindex);
    Mem.myfree(Exportflag);
  }

 public:
  /* Implements a repeated loop over the local particles in the list, processing them with the local kernel function,
   * until we're done or the export buffer is full. Then we exchange the data, and process the imported ones with the provided kernel.
   * We repeat if neeed until all processors are done.
   */
  int execute(int nactive, int *targetlist, int action)
  {
    generic_allocate_comm_tables();

    Nactive    = nactive;
    Targets    = targetlist;
    actioncode = action;

    Thread.Interactions = 0;

    int ndone_flag, ndone, iter = 0;

    SumNexport = 0; /* can be queried as a book-keeping variable */

    NextParticle = 0; /* first particle index for this task */

    if(cpu_imbalance != logs::CPU_NONE)
      Logs.timer_start(cpu_imbalance);

    do
      {
        iter++;

        /* allocate buffers to arrange communication */
        generic_alloc_partlist_nodelist_ngblist_threadbufs();

        /* do local particles */
        if(action == MODE_DEFAULT)
          {
            primary_loop(MODE_LOCAL_PARTICLES);
          }
        else if(action == MODE_LOCAL_NO_EXPORT)
          {
            primary_loop(MODE_LOCAL_NO_EXPORT);
          }
        else
          {
            Terminate("unknown action code (%d) for primary loop", action);
          }

        /* do all necessary bookkeeping, data exchange, and processing of imported particles */
        generic_exchange();

        /* free the rest of the buffers */
        generic_free_partlist_nodelist_ngblist_threadbufs();

        /* check whether we are done */
        if(NextParticle >= nactive)
          ndone_flag = 1;
        else
          ndone_flag = 0;

        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, D->Communicator);
      }
    while(ndone < D->NTask);

    generic_free_comm_tables();

    if(cpu_imbalance != logs::CPU_NONE)
      Logs.timer_stop(cpu_imbalance);

    return iter;
  }

  int execute(int nactive, int *targetlist, int action, enum logs::timers a, enum logs::timers b, enum logs::timers c)
  {
    cpu_primary   = a;
    cpu_secondary = b;
    cpu_imbalance = c;

    return execute(nactive, targetlist, action);
  }
};

#endif
