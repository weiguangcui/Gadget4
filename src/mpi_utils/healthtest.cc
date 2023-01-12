/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  healthtest.cc
 *
 *  \brief routines for testing whether all compute nodes yield the full CPU and communication performance
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
#include "../logs/logs.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../system/system.h"

#define TEST_PACKET_SIZE_IN_MB 5
#define WORK_LOOP_COUNTER 50000000
#define WORK_NUMBER_OF_IPROBE_TESTS 1000000

#ifndef MAX_VARIATION_TOLERANCE
#define MAX_VARIATION_TOLERANCE 0.5
#endif

void sim::healthtest(void)
{
  mpi_printf("\n");

  measure_cpu_performance(Communicator);

  // Let's take a look at the communication speed in a global all-to-all data exchange realized through pairwise exchanges along a
  // hypercube
  if(NTask > 1)
    measure_hyper_cube_speed("Full hypercube:", Communicator);

  // Let's take a look at inter-node communication speed
  if(NumNodes > 1)
    {
      int CommSplitColor;

      if(RankInThisNode == 0)
        CommSplitColor = 0;
      else
        CommSplitColor = 1;

      MPI_Comm comm;
      MPI_Comm_split(Communicator, CommSplitColor, ThisTask, &comm);

      if(RankInThisNode == 0)
        measure_hyper_cube_speed("Internode cube:", comm);

      MPI_Comm_free(&comm);
    }

  // Now look at intra-node communication speed
  if(NumNodes < NTask)
    {
      int CommSplitColor = ThisNode;
      MPI_Comm comm;
      MPI_Comm_split(Communicator, CommSplitColor, ThisTask, &comm);

      measure_hyper_cube_speed("Intranode cube, 1st node:", comm);

      MPI_Comm_free(&comm);
    }

  measure_iprobe_performance("Iprobe for any message:");

  mpi_printf("\n");
}

double sim::measure_cpu_performance(MPI_Comm Communicator)
{
  int loc_ntask, loc_thistask, loc_ptask;

  double ta = Logs.second();

  MPI_Comm_rank(Communicator, &loc_thistask);
  MPI_Comm_size(Communicator, &loc_ntask);

  for(loc_ptask = 0; loc_ntask > (1 << loc_ptask); loc_ptask++)
    ;

  double sum = 0;

  MPI_Barrier(Communicator);

  double t0 = Logs.second();

  // do some computationally intense (but useless) work for a while
  for(int i = 0; i < WORK_LOOP_COUNTER; i++)
    sum += sin((i + 0.1) / WORK_LOOP_COUNTER) / (2.0 + cos(i - 0.1) / WORK_LOOP_COUNTER);

  double t1 = Logs.second();

  double tperf = Logs.timediff(t0, t1), tperfsum;

  MPI_Allreduce(&tperf, &tperfsum, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  double tavg = tperfsum / loc_ntask;

  struct
  {
    double t;
    int rank;
  } local = {tperf, ThisTask}, localnode = {tperf, ThisNode}, min_time, max_time, min_timenode, max_timenode;

  MPI_Allreduce(&local, &min_time, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&local, &max_time, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  MPI_Allreduce(&localnode, &min_timenode, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&localnode, &max_timenode, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  double variation = (max_time.t - min_time.t) / tavg;

  double tb = Logs.second();

  mpi_printf(
      "HEALTHTEST: %25s  %8.3f sec            %7.3f%%  variation   | Best=%g on Task=%d/Node=%d, Worst=%g on Task=%d/Node=%d, test "
      "took %g sec\n",
      "CPU performance:", tavg, 100.0 * variation, min_time.t, min_time.rank, min_timenode.rank, max_time.t, max_time.rank,
      max_timenode.rank, Logs.timediff(ta, tb));

  if(variation >= MAX_VARIATION_TOLERANCE)
    {
      char name_maxnode[MPI_MAX_PROCESSOR_NAME];
      int len;
      if(ThisTask == max_time.rank)
        MPI_Get_processor_name(name_maxnode, &len);

      MPI_Bcast(name_maxnode, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, max_time.rank, Communicator);

      char buf[1000 + MPI_MAX_PROCESSOR_NAME];
      snprintf(buf, 1000 + MPI_MAX_PROCESSOR_NAME, "processes_%s.txt", name_maxnode);

      mpi_printf("HEALTHTEST: We are dumping a process list to the file '%s'\n", buf);

      if(ThisTask == max_time.rank)
        {
          char cmd[10000 + MPI_MAX_PROCESSOR_NAME];
          snprintf(cmd, 10000 + MPI_MAX_PROCESSOR_NAME, "ps -ef >& %s", buf);
          system(cmd);
        }

      MPI_Barrier(Communicator);

      // only issue a warning for now instead of terminating the code
      warn(
          "\n\nHEALTHTEST: We issue a warning because the performance variation=%g of the CPUs lies above the prescribed tolerance "
          "MAX_VARIATION_TOLERANCE=%g, possibly indicating a machine problem. (sum=%g)\n",
          variation, MAX_VARIATION_TOLERANCE, sum);
    }

  return sum;
}

double sim::measure_hyper_cube_speed(const char *tag, MPI_Comm Communicator)
{
  double ta = Logs.second();

  int loc_ntask, loc_thistask, loc_ptask;

  MPI_Comm_rank(Communicator, &loc_thistask);
  MPI_Comm_size(Communicator, &loc_ntask);

  for(loc_ptask = 0; loc_ntask > (1 << loc_ptask); loc_ptask++)
    ;

  int bytecount = (TEST_PACKET_SIZE_IN_MB * 1024L * 1024L) / loc_ntask;

  double tall = 0;
  int count   = 0;

  char *sendbuf = (char *)Mem.mymalloc_clear("send", bytecount * sizeof(char));
  char *recvbuf = (char *)Mem.mymalloc_clear("recv", bytecount * sizeof(char));

  /* exchange the test data */
  for(int ngrp = 1; ngrp < (1 << loc_ptask); ngrp++)
    {
      int recvTask = loc_thistask ^ ngrp;

      MPI_Barrier(Communicator);

      if(recvTask < loc_ntask)
        {
          double t0 = Logs.second();
          myMPI_Sendrecv(sendbuf, bytecount, MPI_BYTE, recvTask, TAG_DENS_A, recvbuf, bytecount, MPI_BYTE, recvTask, TAG_DENS_A,
                         Communicator, MPI_STATUS_IGNORE);
          double t1 = Logs.second();

          tall += Logs.timediff(t0, t1);
          count++;
        }
    }

  Mem.myfree(recvbuf);
  Mem.myfree(sendbuf);

  double tperf = 0.5 * tall / count, tperfsum;

  MPI_Allreduce(&tperf, &tperfsum, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  double tavg = tperfsum / loc_ntask;

  struct
  {
    double t;
    int rank;
  } local = {tperf, ThisTask}, localnode = {tperf, ThisNode}, min_time, max_time, min_timenode, max_timenode;

  MPI_Allreduce(&local, &min_time, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&local, &max_time, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  MPI_Allreduce(&localnode, &min_timenode, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&localnode, &max_timenode, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  double tb = Logs.second();

  double variation = (bytecount / min_time.t - bytecount / max_time.t) / (bytecount / tavg);

  mpi_printf(
      "HEALTHTEST: %25s  %8.1f MB/s per pair  %7.3f%%  variation   | Best=%g on Task=%d/Node=%d, Worst=%g on Task=%d/Node=%d, test "
      "took %g sec\n",
      tag, bytecount / tavg * TO_MBYTE_FAC, 100.0 * variation, bytecount / min_time.t * TO_MBYTE_FAC, min_time.rank, min_timenode.rank,
      bytecount / max_time.t * TO_MBYTE_FAC, max_time.rank, max_timenode.rank, Logs.timediff(ta, tb));

  if(variation > MAX_VARIATION_TOLERANCE && ThisTask == 0)
    warn(
        "\nThe performance variation=%g of the communication speed lies above the prescribed tolerance MAX_VARIATION_TOLERANCE=%g, "
        "possibly indicating a machine problem.\n",
        variation, MAX_VARIATION_TOLERANCE);

  return tavg;
}

void sim::measure_iprobe_performance(const char *tag)
{
  double ta = Logs.second();

  for(int i = 0; i < WORK_NUMBER_OF_IPROBE_TESTS; i++)
    {
      int flag;
      MPI_Status status;

      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, Communicator, &flag, &status);
    }

  double tb = Logs.second();

  double tperf = Logs.timediff(ta, tb) / WORK_NUMBER_OF_IPROBE_TESTS;

  struct
  {
    double t;
    int rank;
  } local = {tperf, ThisTask}, min_time, max_time;

  MPI_Allreduce(&local, &min_time, 1, MPI_DOUBLE_INT, MPI_MINLOC, Communicator);
  MPI_Allreduce(&local, &max_time, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  double tperfsum;
  MPI_Allreduce(&tperf, &tperfsum, 1, MPI_DOUBLE, MPI_SUM, Communicator);
  double tavg = tperfsum / NTask;

  char name_minnode[MPI_MAX_PROCESSOR_NAME];
  char name_maxnode[MPI_MAX_PROCESSOR_NAME];

  int len;
  if(ThisTask == min_time.rank)
    MPI_Get_processor_name(name_minnode, &len);
  if(ThisTask == max_time.rank)
    MPI_Get_processor_name(name_maxnode, &len);

  MPI_Bcast(name_minnode, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, min_time.rank, Communicator);
  MPI_Bcast(name_maxnode, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, max_time.rank, Communicator);

  double variation = (max_time.t - min_time.t) / tavg;

  mpi_printf(
      "HEALTHTEST: %25s  %g s per MPI_Ip%7.3f%%  variation   | Best=%g on Task=%d/Node=%s, Worst=%g on Task=%d/Node=%s, test took %g "
      "sec\n",
      tag, tavg, 100.0 * variation, min_time.t, min_time.rank, name_minnode, max_time.t, max_time.rank, name_maxnode,
      Logs.timediff(ta, tb));
}
