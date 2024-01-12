/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  system.cc
 *
 *  \brief various low level helper routines
 */

#include "gadgetconfig.h"

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/resource.h>
#include <sys/statvfs.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../io/io.h"
#include "../logs/logs.h"
#include "../main/main.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../system/system.h"

gsl_rng *random_generator; /*!< the random number generator used */

void init_rng(int thistask)
{
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42 + thistask); /* start-up seed, choose different on each rank */
}

double get_random_number(void) { return gsl_rng_uniform(random_generator); }

void subdivide_evenly(long long N, int pieces, int index_bin, long long *first, int *count)
{
  if(N / pieces > INT_MAX)
    Terminate("overflow");

  int nbase      = N / pieces;
  int additional = N % pieces;
  *first         = index_bin * ((long long)nbase) + ((index_bin < additional) ? index_bin : additional);
  *count         = nbase + (index_bin < additional);
}

void subdivide_evenly(int N, int pieces, int index_bin, int *first, int *count)
{
  int nbase      = N / pieces;
  int additional = N % pieces;
  *first         = index_bin * nbase + ((index_bin < additional) ? index_bin : additional);
  *count         = nbase + (index_bin < additional);
}

void subdivide_evenly_get_bin(int N, int pieces, int index, int *bin)
{
  int nbase      = N / pieces;
  int additional = N % pieces;

  if(index < additional * (nbase + 1))
    *bin = index / (nbase + 1);
  else
    *bin = (index - additional) / nbase;
}

void permutate_chunks_in_list(int ncount, int *list)
{
#define WALK_N_PIECES 32 /*!< Number of sets, the chunks are divided into */
#define WALK_N_SIZE 500  /*!< Number of particles per chunk */

  int nchunk     = 1;      /*!< Number of chunk sets used */
  int nchunksize = ncount; /*!< Size of each chunk */

  if(ncount > WALK_N_PIECES * WALK_N_SIZE)
    {
      nchunk     = WALK_N_PIECES;
      nchunksize = WALK_N_SIZE;
    }

  int currentchunk = 0; /*!< Chunk set currently processed */

  int *chunked_TargetList = (int *)Mem.mymalloc("chunked_TargetList", ncount * sizeof(int));
  for(int n = 0, nextparticle = 0; n < ncount; n++)
    {
      int i = nextparticle;

      chunked_TargetList[n] = list[i];
      if(i < ncount)
        {
          nextparticle++;

          if((nextparticle % nchunksize) == 0)
            nextparticle += (nchunk - 1) * nchunksize;

          if(nextparticle >= ncount)
            {
              currentchunk++;
              if(currentchunk < nchunk)
                nextparticle = currentchunk * nchunksize;
            }
        }
    }

  for(int n = 0; n < ncount; n++)
    list[n] = chunked_TargetList[n];

  Mem.myfree(chunked_TargetList);
}

void myflush(FILE *fstream)
{
#ifdef REDUCE_FLUSH
  /* do nothing */
  (void)fstream;
#else
  fflush(fstream);
#endif
}

#ifdef DEBUG
#include <fenv.h>
void enable_core_dumps_and_fpu_exceptions(void)
{
#if defined(DEBUG_ENABLE_FPU_EXCEPTIONS) && defined(__linux__)
  /* enable floating point exceptions */

  extern int feenableexcept(int __excepts);
  feenableexcept(FE_DIVBYZERO | FE_INVALID);

  /* Note: FPU exceptions appear not to work properly
   * when the Intel C-Compiler for Linux is used
   */
#endif

  /* set core-dump size to infinity */
  rlimit rlim;
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  /* MPICH catches the signales SIGSEGV, SIGBUS, and SIGFPE....
   * The following statements reset things to the default handlers,
   * which will generate a core file.
   */
  signal(SIGSEGV, SIG_DFL);
  signal(SIGBUS, SIG_DFL);
  signal(SIGFPE, SIG_DFL);
  signal(SIGINT, SIG_DFL);
}
#endif

long long sim::report_comittable_memory(long long *MemTotal, long long *Committed_AS, long long *SwapTotal, long long *SwapFree)
{
  FILE *fd;
  char buf[1024];

  if((fd = fopen("/proc/meminfo", "r")))
    {
      while(1)
        {
          if(fgets(buf, 500, fd) != buf)
            break;

          if(bcmp(buf, "MemTotal", 8) == 0)
            {
              *MemTotal = atoll(buf + 10);
            }
          if(strncmp(buf, "Committed_AS", 12) == 0)
            {
              *Committed_AS = atoll(buf + 14);
            }
          if(strncmp(buf, "SwapTotal", 9) == 0)
            {
              *SwapTotal = atoll(buf + 11);
            }
          if(strncmp(buf, "SwapFree", 8) == 0)
            {
              *SwapFree = atoll(buf + 10);
            }
        }
      fclose(fd);
    }

  return (*MemTotal - *Committed_AS);
}

// returns kByte
long long sim::report_free_size_in_tmpfs(void)
{
  struct statvfs buf;

  if(statvfs("/dev/shm", &buf) == 0)
    {
      return (1024.0 * (uint64_t)buf.f_bsize * ((uint64_t)(buf.f_bavail))) * TO_MBYTE_FAC;
    }
  else
    return 0;
}

void sim::mpi_report_comittable_memory(void)
{
  long long *sizelist, maxsize[7], minsize[7];
  double avgsize[7];
  int i, imem, mintask[7], maxtask[7];
  long long Mem[7];
  char label[MAXLEN_PATH];

  Mem[0] = report_comittable_memory(&Mem[1], &Mem[2], &Mem[3], &Mem[4]);
  Mem[5] = Mem[1] - Mem[0];

  Mem[6] = report_free_size_in_tmpfs();

  MemoryOnNode       = Mem[1];
  SharedMemoryOnNode = Mem[6];

  for(imem = 0; imem < 7; imem++)
    {
      sizelist = (long long *)malloc(NTask * sizeof(long long));
      MPI_Allgather(&Mem[imem], sizeof(long long), MPI_BYTE, sizelist, sizeof(long long), MPI_BYTE, Communicator);

      for(i = 1, mintask[imem] = 0, maxtask[imem] = 0, maxsize[imem] = minsize[imem] = sizelist[0], avgsize[imem] = sizelist[0];
          i < NTask; i++)
        {
          if(sizelist[i] > maxsize[imem])
            {
              maxsize[imem] = sizelist[i];
              maxtask[imem] = i;
            }
          if(sizelist[i] < minsize[imem])
            {
              minsize[imem] = sizelist[i];
              mintask[imem] = i;
            }
          avgsize[imem] += sizelist[i];
        }

      free(sizelist);
    }

  if(ThisTask == 0)
    {
      printf(
          "\n-------------------------------------------------------------------------------------------------------------------------"
          "\n");
      for(imem = 0; imem < 7; imem++)
        {
          switch(imem)
            {
              case 0:
                snprintf(label, MAXLEN_PATH, "AvailMem");
                break;
              case 1:
                snprintf(label, MAXLEN_PATH, "Total Mem");
                break;
              case 2:
                snprintf(label, MAXLEN_PATH, "Committed_AS");
                break;
              case 3:
                snprintf(label, MAXLEN_PATH, "SwapTotal");
                break;
              case 4:
                snprintf(label, MAXLEN_PATH, "SwapFree");
                break;
              case 5:
                snprintf(label, MAXLEN_PATH, "AllocMem");
                break;
              case 6:
                snprintf(label, MAXLEN_PATH, "avail /dev/shm");
                break;
            }
          printf("%s:\t Largest = %10.2f Mb (on task=%4d), Smallest = %10.2f Mb (on task=%4d), Average = %10.2f Mb\n", label,
                 maxsize[imem] / (1024.0), maxtask[imem], minsize[imem] / (1024.0), mintask[imem], avgsize[imem] / (1024.0 * NTask));
        }
      printf(
          "-------------------------------------------------------------------------------------------------------------------------"
          "\n");
    }

  char name[MPI_MAX_PROCESSOR_NAME];

  if(ThisTask == maxtask[2])
    {
      int len;
      MPI_Get_processor_name(name, &len);
    }

  MPI_Bcast(name, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, maxtask[2], Communicator);

  if(ThisTask == 0)
    {
      printf("Task=%d has the maximum commited memory and is host: %s\n", ThisTask, name);
      printf(
          "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          "\n");
    }

  fflush(stdout);
}

/* the following function finds the last significant bit, as in the linux kernel */
int my_fls(unsigned int x)
{
  int r = 32;

  if(!x)
    return 0;
  if(!(x & 0xffff0000u))
    {
      x <<= 16;
      r -= 16;
    }
  if(!(x & 0xff000000u))
    {
      x <<= 8;
      r -= 8;
    }
  if(!(x & 0xf0000000u))
    {
      x <<= 4;
      r -= 4;
    }
  if(!(x & 0xc0000000u))
    {
      x <<= 2;
      r -= 2;
    }
  if(!(x & 0x80000000u))
    {
      x <<= 1;
      r -= 1;
    }
  return r;
}
