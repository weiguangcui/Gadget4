/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file test_io_bandwidth.cc
 *
 * \brief test routines to identify optimum setting of MaxFilesWithConcurrentIO on given machine
 */

#include "gadgetconfig.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../io/io.h"
#include "../io/test_io_bandwidth.h"
#include "../lightcone/lightcone.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngbtree/ngbtree.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

void test_io_bandwidth::measure_io_bandwidth(void)
{
  /* create directory for test data */
  if(ThisTask == 0)
    {
      char buf[MAXLEN_PATH_EXTRA];
      snprintf(buf, MAXLEN_PATH_EXTRA, "%s/testdata", All.OutputDir);
      mkdir(buf, 02755);
    }
  MPI_Barrier(Communicator);

  All.MaxFilesWithConcurrentIO = NTask;

  while(All.MaxFilesWithConcurrentIO > 1)
    {
      write_test_data();

      All.MaxFilesWithConcurrentIO /= 2;
    }

  mpi_printf("\n\nTEST: Completed.\n");

  fflush(stdout);
}

void test_io_bandwidth::write_test_data(void)
{
  double t0 = Logs.second();
  reset_io_byte_count();

  mpi_printf("TEST: Writing test data...\n");

  /* now work the I/O of the files, controlled by scheduler to achieve optimum I/O bandwidth under the constraint of a maximum number
   * for the concurrent file access */
  work_files(MODUS_WRITE);

  long long byte_count = get_io_byte_count(), byte_count_all;
  sumup_longs(1, &byte_count, &byte_count_all, Communicator);

  double t1 = Logs.second();

  mpi_printf(
      "TEST: done.  MaxFilesWithConcurrentIO=%6d   load/save took %g sec, total size %g MB, corresponds to effective I/O rate of %g "
      "MB/sec\n",
      All.MaxFilesWithConcurrentIO, Logs.timediff(t0, t1), byte_count_all / (1024.0 * 1024.0),
      byte_count_all / (1024.0 * 1024.0) / Logs.timediff(t0, t1));

  /* now delete test data */
  char buf[MAXLEN_PATH_EXTRA];
  snprintf(buf, MAXLEN_PATH_EXTRA, "%s/testdata/%s.%d", All.OutputDir, "testdata", ThisTask);
  unlink(buf);
  MPI_Barrier(Communicator);
}

void test_io_bandwidth::polling(int modus)
{
  if(ThisTask == 0)
    if(files_completed < NTask)
      {
        MPI_Status status;
        int flag;

        /* now check for a completion message  */
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_KEY, Communicator, &flag, &status);

        if(flag)
          {
            int source = status.MPI_SOURCE;

            int dummy;
            MPI_Recv(&dummy, 1, MPI_INT, source, TAG_KEY, Communicator, MPI_STATUS_IGNORE);
            files_completed++;

            if(files_started < NTask)
              {
                /* send start signal */
                MPI_Ssend(&ThisTask, 1, MPI_INT, seq[files_started++].thistask, TAG_N, Communicator);
              }
          }
      }
}

void test_io_bandwidth::work_files(int modus)
{
  if(ThisTask == 0)
    if(!(seq = (seq_data *)malloc(NTask * sizeof(seq_data))))
      Terminate("can't allocate seq_data");

  seq_data seq_loc;
  seq_loc.thistask   = ThisTask;
  seq_loc.rankinnode = RankInThisNode;
  seq_loc.thisnode   = ThisNode;

  MPI_Gather(&seq_loc, sizeof(seq_data), MPI_BYTE, seq, sizeof(seq_data), MPI_BYTE, 0, Communicator);

  if(ThisTask == 0)
    {
      std::sort(seq, seq + NTask);

      files_started   = 0;
      files_completed = 0;

      for(int i = 1; i < All.MaxFilesWithConcurrentIO; i++)
        {
          files_started++;
          MPI_Ssend(&ThisTask, 1, MPI_INT, seq[i].thistask, TAG_N, Communicator);
        }

      files_started++;
      contents_restart_file(modus);
      files_completed++;

      if(files_started < NTask)
        {
          /* send start signal */
          MPI_Ssend(&ThisTask, 1, MPI_INT, seq[files_started++].thistask, TAG_N, Communicator);
        }

      while(files_completed < NTask)
        polling(modus);

      free(seq);
    }
  else
    {
      /* wait for start signal */
      int dummy;
      MPI_Recv(&dummy, 1, MPI_INT, 0, TAG_N, Communicator, MPI_STATUS_IGNORE); /* wait until we are told to start */

      contents_restart_file(modus);

      /* send back completion notice */
      MPI_Ssend(&ThisTask, 1, MPI_INT, 0, TAG_KEY, Communicator);
    }
}

void test_io_bandwidth::contents_restart_file(int modus)
{
  char buf[MAXLEN_PATH_EXTRA];
  snprintf(buf, MAXLEN_PATH_EXTRA, "%s/testdata/%s.%d", All.OutputDir, "testdata", ThisTask);

  if(modus == MODUS_READ)
    {
      if(!(fd = fopen(buf, "r")))
        {
          Terminate("TEST: File '%s' not found.\n", buf);
        }
    }
  else if(modus == MODUS_WRITE)
    {
      if(!(fd = fopen(buf, "w")))
        {
          Terminate("TEST: File '%s' cannot be opened.\n", buf);
        }
    }
  else
    Terminate("unknown modus\n");

  size_t len = BUF_IN_MB * 1024LL * 1024LL;

  char *p = (char *)Mem.mymalloc("p", len);

  byten(p, len, modus);

  Mem.myfree(p);

  fclose(fd);
}

void test_io_bandwidth::byten(void *x, size_t n, int modus)
{
  char *p = (char *)x;

  while(n > BLKSIZE)
    {
      byten_doit(p, BLKSIZE, modus);
      p += BLKSIZE;
      n -= BLKSIZE;
      polling(modus);
    }

  if(n > 0)
    byten_doit(p, n, modus);
}

/*! \brief reads/writes n bytes to a restart file
 */
void test_io_bandwidth::byten_doit(void *x, size_t n, int modus)
{
  if(modus == MODUS_READ)
    my_fread(x, n, 1, fd);
  else
    my_fwrite(x, n, 1, fd);
}
