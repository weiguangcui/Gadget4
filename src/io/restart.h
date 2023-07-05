/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file restart.h
 *
 * \brief declares class used for reading/writing of restart files
 */

#ifndef RESTART_H
#define RESTART_H

#define MODUS_WRITE 0
#define MODUS_READ 1

#define BLKSIZE (1024 * 1024)

#include "gadgetconfig.h"

#include "../io/io_streamcount.h"
#include "../main/simulation.h"

class restart : public io_streamcount, public setcomm
{
 public:
  restart(MPI_Comm comm) : setcomm(comm) /* constructor */ { determine_compute_nodes(); }

  void load(sim *Sim_ptr);
  void write(sim *Sim_ptr);

 private:
  sim *Sim;

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
  int files_concurrent;
  int files_groups;

  void do_restart(int modus);

  void readjust_timebase(double TimeMax_old, double TimeMax_new);
  void work_files(int modus);
  void contents_restart_file(int modus);
  void backup_restartfiles(int task);
  void polling(int modus);
  void in(int *x, int modus);
  void byten(void *x, size_t n, int modus);
  void byten_doit(void *x, size_t n, int modus);
};

#endif /* RESTART_H */
