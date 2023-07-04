/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file test_io_bandwidth.h
 *
 * \brief declares a class for I/O performance measurements
 */

#ifndef TEST_IO_BANDWIDTH_H
#define TEST_IO_BANDWIDTH_H

#include "gadgetconfig.h"

#include <mpi.h>

#define MODUS_WRITE 0
#define MODUS_READ 1

#define BUF_IN_MB 10

#define BLKSIZE (1024 * 1024)

#include "../io/io_streamcount.h"

class test_io_bandwidth : public io_streamcount, public virtual setcomm
{
 public:
  test_io_bandwidth(MPI_Comm comm) : setcomm(comm) {}

  void measure_io_bandwidth(void);

 private:
  FILE *fd;

  struct seq_data
  {
    int thistask;
    int rankinnode;
    int thisnode;

    bool operator<(const seq_data &other) const
    {
      if(rankinnode < other.rankinnode)
        return true;
      if(rankinnode > other.rankinnode)
        return false;
      if(thisnode < other.thisnode)
        return true;
      if(thisnode > other.thisnode)
        return false;
      return thistask < other.thistask;
    }
  };
  seq_data *seq;

  int files_started;
  int files_completed;

  void work_files(int modus);
  void contents_restart_file(int modus);
  void polling(int modus);
  void write_test_data(void);
  void byten(void *x, size_t n, int modus);
  void byten_doit(void *x, size_t n, int modus);
};

#endif
