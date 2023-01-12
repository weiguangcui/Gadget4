/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/
#include "gadgetconfig.h"

#include <fcntl.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../logs/logs.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../system/system.h"

/** \file mymalloc.cc
 *
 *  \brief Manager for dynamic memory allocation
 *
 *  This module handles the dynamic memory allocation.
 *  To avoid memory allocation/dellocation overhead a big chunk of memory
 *  (which will be the maximum amount of dinamically allocatable memory)
 *  is allocated upon initialization. This chunk is then filled by the memory
 *  blocks as in a stack structure. The blocks are automatically aligned to a CACHELINESIZE bit boundary.
 *  Memory blocks come in two flavours: movable and non-movable. In non-movable
 *  blocks the starting address is fixed once the block is allocated and cannot be changed.
 *  Due to the stack structure of the dynamic memory, this implies that the last (non-movable)
 *  block allocated must be the first block to be deallocated. If this condition is not met,
 *  an abort condition is triggered. If more flexibility is needed, movable memory blocks can
 *  be used. In this case, the starting address of the block is again fixed upon allocation
 *  but the block can be shifted (therefore its initial address changes) according to needs.
 *  For a movable block to be successfully shifted it is required that all the subsequent allocated
 *  blocks are movable. Again, an abort condition is triggered if this condition is not met.
 *  Movable blocks can be deallocated in any order provided that the condition just described holds.
 *  The gap resulting form the deallocation of a block that is not in
 *  the last position will be automatically filled by shifting all the blocks coming after the
 *  deallocated block.
 */

/** \brief Initialize memory manager.
 *
 *  This function initializes the memory manager. In particular, it sets
 *  the global variables of the module to their initial value and allocates
 *  the memory for the stack.
 */
void memory::mymalloc_init(int maxmemsize, enum restart_options restartflag)
{
  BlockSize                    = (size_t *)malloc(MAXBLOCKS * sizeof(size_t));
  Table                        = (char **)malloc(MAXBLOCKS * sizeof(void *));
  MovableFlag                  = (char *)malloc(MAXBLOCKS * sizeof(char));
  GenericFlag                  = (char *)malloc(MAXBLOCKS * sizeof(char));
  BasePointers                 = (char ***)malloc(MAXBLOCKS * sizeof(void **));
  VarName                      = (char *)malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FunctionName                 = (char *)malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  ParentFileName               = (char *)malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  FileName                     = (char *)malloc(MAXBLOCKS * MAXCHARS * sizeof(char));
  LineNumber                   = (int *)malloc(MAXBLOCKS * sizeof(int));
  highmark_bufsize             = (100 + 4 * MAXCHARS) * (MAXBLOCKS + 10);
  HighMarkTabBuf               = (char *)malloc(highmark_bufsize);
  HighMarkTabBufWithoutGeneric = (char *)malloc(highmark_bufsize);

  memset(VarName, 0, MAXBLOCKS * MAXCHARS);
  memset(FunctionName, 0, MAXBLOCKS * MAXCHARS);
  memset(ParentFileName, 0, MAXBLOCKS * MAXCHARS);
  memset(FileName, 0, MAXBLOCKS * MAXCHARS);

  size_t n = maxmemsize * ((size_t)1024 * 1024);
  n        = roundup_to_multiple_of_cacheline_size(n);

  // add an extra cache line size to make sure we can guarantee that base returned from MPI_Win_allocate_shared is aligned
  n += CACHELINESIZE;

  RestartFlag = restartflag;

  MPI_Barrier(MPI_COMM_WORLD);  // wait until both regular and ghost processors are here

  double t0 = Logs.second();

  MPI_Info win_info;
  MPI_Info_create(&win_info);
  MPI_Info_set(win_info, "alloc_shared_noncontig", "true");

  if(myMPI_Win_allocate_shared(n, 1, win_info, Shmem.SharedMemComm, &Base, &Shmem.SharedMemWin) != MPI_SUCCESS)
    Terminate("Failed to allocate memory for `Base' (%d Mbytes).\n", All.MaxMemSize);

  /* we now make sure that the allocated local buffer is really aligned, not all MPI libraries guarantee this */

  int off = 0;
  if((long long)Base / CACHELINESIZE > 0)
    off = ((long long)Base / CACHELINESIZE + 1) * CACHELINESIZE - (long long)Base;

  /* allign our base */
  Base += off;

  MPI_Info_free(&win_info);

  double t1 = Logs.second();
  if(Shmem.World_ThisTask == 0)
    mpi_printf("MALLOC: Allocation of shared memory took %g sec\n", Logs.timediff(t0, t1));

  TotBytes = FreeBytes = n - CACHELINESIZE;

  AllocatedBytes                  = 0;
  Nblocks                         = 0;
  HighMarkBytes                   = 0;
  HighMarkBytesWithoutGeneric     = 0;
  OldGlobHighMarkMB               = 0;
  OldGlobHighMarkMBWithoutGeneric = 0;

  char mode[2], buf[MAXLEN_PATH_EXTRA];

  if(All.RestartFlag == RST_BEGIN)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

  MPI_Bcast(All.OutputDir, sizeof(All.OutputDir), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(Shmem.GhostRank == 0)
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "memory.txt");
  else
    snprintf(buf, MAXLEN_PATH_EXTRA, "%s%s", All.OutputDir, "memory_ghostranks.txt");

  if(!(FdMemory = fopen(buf, mode)))
    Terminate("error in opening file '%s'\n", buf);

  /* tell also the ghost ranks about the total size of the simulation partition */
  MPI_Bcast(&Shmem.Sim_NTask, 1, MPI_INT, 0, MPI_COMM_WORLD);

  Shmem.GetGhostRankForSimulCommRank = (int *)mymalloc("GetGhostRankForSimulCommRank", Shmem.Sim_NTask * sizeof(int));
  Shmem.GetShmRankForSimulCommRank   = (int *)mymalloc("GetShmRankForSimulCommRank", Shmem.Sim_NTask * sizeof(int));
  Shmem.GetNodeIDForSimulCommRank    = (int *)mymalloc("GetNodeIDForSimulCommRank", Shmem.Sim_NTask * sizeof(int));

  if(Shmem.GhostRank == 0)
    {
      MPI_Allgather(&Shmem.MyShmRankInGlobal, 1, MPI_INT, Shmem.GetGhostRankForSimulCommRank, 1, MPI_INT, Shmem.SimulationComm);
      MPI_Allgather(&Shmem.Island_ThisTask, 1, MPI_INT, Shmem.GetShmRankForSimulCommRank, 1, MPI_INT, Shmem.SimulationComm);
      MPI_Allgather(&Shmem.Island_Smallest_WorldTask, 1, MPI_INT, Shmem.GetNodeIDForSimulCommRank, 1, MPI_INT, Shmem.SimulationComm);
    }

  // to make sure that also the ghost processors have this table
  MPI_Bcast(Shmem.GetGhostRankForSimulCommRank, Shmem.Sim_NTask, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(Shmem.GetShmRankForSimulCommRank, Shmem.Sim_NTask, MPI_INT, 0, MPI_COMM_WORLD);

  /* we also need the base offsets of the other MPI ranks in the same shared memory island */
  Shmem.SharedMemBaseAddr = (void **)mymalloc("SharedMemBaseAddr", Shmem.Island_NTask * sizeof(void *));

  for(int i = 0; i < Shmem.Island_NTask; i++)
    {
      MPI_Aint size;
      int disp_unit;
      myMPI_Win_shared_query(Shmem.SharedMemWin, i, &size, &disp_unit, &Shmem.SharedMemBaseAddr[i]);
    }

  // now propagte the alignment correction also to the base addresses that all the other processes see
  int *off_list = (int *)Mem.mymalloc("off_list", Shmem.Island_NTask * sizeof(int));

  MPI_Allgather(&off, 1, MPI_INT, off_list, 1, MPI_INT, Shmem.SharedMemComm);

  for(int i = 0; i < Shmem.Island_NTask; i++)
    Shmem.SharedMemBaseAddr[i] = (char *)Shmem.SharedMemBaseAddr[i] + off_list[i];

  Mem.myfree(off_list);
}

int memory::myMPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win)
{
#ifndef ALLOCATE_SHARED_MEMORY_VIA_POSIX
  return MPI_Win_allocate_shared(size, disp_unit, info, comm, baseptr, win);
#else

  char shmpath[NAME_MAX];

  /* Base offsets of the other MPI ranks in our shared memory mapping */
  Shmem.SharedMemBaseAddrRaw = (char **)malloc(Shmem.Island_NTask * sizeof(char *));

  long long *size_list = (long long *)malloc(Shmem.Island_NTask * sizeof(long long));

  long long loc_bytes = size;
  MPI_Allgather(&loc_bytes, 1, MPI_LONG_LONG, size_list, 1, MPI_LONG_LONG, comm);

  long long tot_bytes = 0;
  for(int i = 0; i < Shmem.Island_NTask; i++)
    tot_bytes += size_list[i];

  if(Shmem.Island_ThisTask == 0)
    {
      snprintf(shmpath, NAME_MAX, "/G4-%lld.dat", (long long)getpid());

      int fd = shm_open(shmpath, O_CREAT | O_EXCL | O_RDWR, S_IRUSR | S_IWUSR);
      if(fd == -1)
        Terminate("shm_open failed in creation");

      if(ftruncate(fd, tot_bytes) == -1)
        Terminate("ftruncate failed");

      /* Map the object into the caller's address space. */

      void *buf = mmap(NULL, tot_bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
      if(buf == MAP_FAILED)
        Terminate("mmap failed");

      Shmem.SharedMemBaseAddrRaw[0] = (char *)buf;
    }

  MPI_Bcast(shmpath, NAME_MAX, MPI_BYTE, 0, comm);

  if(Shmem.Island_ThisTask != 0)
    {
      int fd = shm_open(shmpath, O_RDWR, 0);
      if(fd == -1)
        Terminate("shm open failed in access");

      void *buf = mmap(NULL, tot_bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
      if(buf == MAP_FAILED)
        Terminate("mmap failed");

      Shmem.SharedMemBaseAddrRaw[0] = (char *)buf;
    }

  for(int i = 1; i < Shmem.Island_NTask; i++)
    Shmem.SharedMemBaseAddrRaw[i] = (char *)Shmem.SharedMemBaseAddrRaw[i - 1] + size_list[i - 1];

  char **p = (char **)baseptr;

  *p = Shmem.SharedMemBaseAddrRaw[Shmem.Island_ThisTask];

  free(size_list);

  MPI_Barrier(comm);

  if(Shmem.Island_ThisTask == 0)
    shm_unlink(shmpath);

  return MPI_SUCCESS;
#endif
}

int memory::myMPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr)
{
#ifndef ALLOCATE_SHARED_MEMORY_VIA_POSIX
  return MPI_Win_shared_query(win, rank, size, disp_unit, baseptr);
#else

  char **p = (char **)baseptr;

  *p = Shmem.SharedMemBaseAddrRaw[rank];

  return MPI_SUCCESS;
#endif
}

void memory::report_memory_usage(int rank, char *tabbuf)
{
  int thistask;
  MPI_Comm_rank(Communicator, &thistask);

  if(thistask == rank)
    {
      int bufsize = (100 + 4 * MAXCHARS) * (Nblocks + 10);
      char *buf   = (char *)mymalloc("buf", bufsize);
      int cc      = 0;
      cc += snprintf(buf + cc, bufsize - cc,
                     "\nMEMORY:  Largest Allocation = %g Mbyte  |  Largest Allocation Without Generic = %g Mbyte\n\n",
                     OldGlobHighMarkMB, OldGlobHighMarkMBWithoutGeneric);

      cc += snprintf(buf + cc, bufsize - cc, "%s", tabbuf);
      if(thistask == 0)
        {
          if(RestartFlag == RST_BEGIN || RestartFlag == RST_RESUME || RestartFlag == RST_STARTFROMSNAP)
            {
              fprintf(FdMemory, "%s", buf);
              fflush(FdMemory);
            }
        }
      else
        {
          MPI_Send(&cc, 1, MPI_INT, 0, TAG_N, Communicator);
          MPI_Send(buf, cc + 1, MPI_BYTE, 0, TAG_PDATA, Communicator);
        }
      myfree(buf);
    }

  if(thistask == 0 && rank > 0)
    {
      int cc;
      MPI_Recv(&cc, 1, MPI_INT, rank, TAG_N, Communicator, MPI_STATUS_IGNORE);
      char *buf = (char *)mymalloc("buf", cc + 1);
      MPI_Recv(buf, cc + 1, MPI_BYTE, rank, TAG_PDATA, Communicator, MPI_STATUS_IGNORE);
      if(RestartFlag == RST_BEGIN || RestartFlag == RST_RESUME || RestartFlag == RST_STARTFROMSNAP)
        {
          fprintf(FdMemory, "%s", buf);
          fflush(FdMemory);
        }
      myfree(buf);
    }
}

/** \brief Output memory usage for the task with the greatest amount of memory allocated.
 *
 */
void memory::report_detailed_memory_usage_of_largest_task(void)
{
  int flag = 0;
  int thistask;
  MPI_Comm_rank(Communicator, &thistask);

  struct
  {
    double mem;
    int rank;
  } local, global;

  local.mem  = HighMarkBytes * TO_MBYTE_FAC;
  local.rank = thistask;

  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  if(global.mem >= 1.05 * OldGlobHighMarkMB)
    {
      OldGlobHighMarkMB = global.mem;
      flag |= 1;
    }

  local.mem  = HighMarkBytesWithoutGeneric * TO_MBYTE_FAC;
  local.rank = thistask;

  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Communicator);

  if(global.mem >= 1.05 * OldGlobHighMarkMBWithoutGeneric)
    {
      OldGlobHighMarkMBWithoutGeneric = global.mem;
      flag |= 2;
    }

  if(flag & 2)
    report_memory_usage(global.rank, HighMarkTabBufWithoutGeneric);

  if(flag & 1)
    report_memory_usage(global.rank, HighMarkTabBuf);
}

/** \brief Dump the buffer where the memory information is stored to the standard output.
 *
 */
void memory::dump_memory_table(void)
{
  int bufsize = 200 * (Nblocks + 10);
  char *buf   = (char *)malloc(bufsize);
  dump_memory_table_buffer(buf, bufsize);
  printf("%s", buf);
  free(buf);
}

/** \brief Fill the output buffer with the memory log.
 *
 *  \param p output buffer
 *  \return the number of characters written to p
 */
int memory::dump_memory_table_buffer(char *p, int bufsize)
{
  int cc              = 0;
  size_t totBlocksize = 0;
  int thistask;
  MPI_Comm_rank(Communicator, &thistask);

  cc += snprintf(p + cc, bufsize - cc, "-------------------------- Allocated Memory Blocks---- ( Step %8d )------------------\n",
                 All.NumCurrentTiStep);
  cc += snprintf(p + cc, bufsize - cc, "Task    Nr F                  Variable      MBytes   Cumulative  Function|File|Linenumber\n");
  cc += snprintf(p + cc, bufsize - cc, "------------------------------------------------------------------------------------------\n");
  for(int i = 0; i < Nblocks; i++)
    {
      totBlocksize += BlockSize[i];

      cc += snprintf(p + cc, bufsize - cc, "%4d %5d %d %40s  %10.4f   %10.4f  %s%s()|%s|%d\n", thistask, i, MovableFlag[i],
                     VarName + i * MAXCHARS, BlockSize[i] * TO_MBYTE_FAC, totBlocksize * TO_MBYTE_FAC, ParentFileName + i * MAXCHARS,
                     FunctionName + i * MAXCHARS, FileName + i * MAXCHARS, LineNumber[i]);
    }
  cc += snprintf(p + cc, bufsize - cc, "------------------------------------------------------------------------------------------\n");

  return cc;
}

/** \brief Allocate a movable memory block and store the relative information.
 *
 *  \param ptr pointer to the initial memory address of the block
 *  \param varname name of the variable to be stored in the allocated block
 *  \param n size of the memory block in bytes
 *  \param func name of function that has called the allocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the allocation routine resides (usually given by the __FILE__ macro)
 *  \param line line number of file where the allocation routine was called (usually given by the __LINE__ macro)
 *  \return a pointer to the beginning of the allocated memory block
 */
void *memory::mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line,
                                        int movable_flag, int clear_flag, char *callorigin)
{
  if((n % CACHELINESIZE) > 0)
    n = (n / CACHELINESIZE + 1) * CACHELINESIZE;

  if(n < CACHELINESIZE)
    n = CACHELINESIZE;

  if(Nblocks >= MAXBLOCKS)
    Terminate("No blocks left in mymalloc_fullinfo() at %s()/%s/line %d. MAXBLOCKS=%d\n", func, file, line, MAXBLOCKS);

  if(n > FreeBytes)
    {
      dump_memory_table();
      Terminate(
          "\nNot enough memory in mymalloc_fullinfo() to allocate %g MB for variable '%s' at %s()/%s/line %d (FreeBytes=%g MB).\n",
          n * TO_MBYTE_FAC, varname, func, file, line, FreeBytes * TO_MBYTE_FAC);
    }
  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  FreeBytes -= n;

  strncpy(VarName + Nblocks * MAXCHARS, varname, MAXCHARS - 1);
  if(callorigin)
    {
      strncpy(ParentFileName + Nblocks * MAXCHARS, callorigin, MAXCHARS - 1);
      GenericFlag[Nblocks] = 1;
      AllocatedBytesGeneric += n;
    }
  else
    {
      memset(ParentFileName + Nblocks * MAXCHARS, 0, MAXCHARS);
      GenericFlag[Nblocks] = 0;
    }
  strncpy(FunctionName + Nblocks * MAXCHARS, func, MAXCHARS - 1);
  strncpy(FileName + Nblocks * MAXCHARS, file, MAXCHARS - 1);
  LineNumber[Nblocks] = line;

  AllocatedBytes += n;
  BlockSize[Nblocks]    = n;
  MovableFlag[Nblocks]  = movable_flag;
  BasePointers[Nblocks] = (char **)ptr;

  Nblocks += 1;

  if(AllocatedBytes - AllocatedBytesGeneric > HighMarkBytesWithoutGeneric)
    {
      HighMarkBytesWithoutGeneric = AllocatedBytes - AllocatedBytesGeneric;
      dump_memory_table_buffer(HighMarkTabBufWithoutGeneric, highmark_bufsize);
    }

  if(AllocatedBytes > HighMarkBytes)
    {
      HighMarkBytes = AllocatedBytes;
      dump_memory_table_buffer(HighMarkTabBuf, highmark_bufsize);
    }

  if(clear_flag)
    memset(Table[Nblocks - 1], 0, n);

  return Table[Nblocks - 1];
}

size_t memory::roundup_to_multiple_of_cacheline_size(size_t n)
{
  if((n % CACHELINESIZE) > 0)
    n = (n / CACHELINESIZE + 1) * CACHELINESIZE;

  return n;
}

void *memory::myfree_query_last_block(void)
{
  if(Nblocks == 0)
    Terminate("no allocated blocks that could be returned");

  return Table[Nblocks - 1];
}

/** \brief Deallocate a movable memory block.
 *
 *  For this operation to be successful all the blocks allocated after the block that has to be freed must be of movable type.
 *
 *  \param p pointer to the memory block to be deallocated
 *  \param func name of function that has called the deallocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the deallocation routine resides (usually given by the __FILE__ macro)
 *  \param line line number of file where the deallocation routine was called (usually given by the __LINE__ macro)
 */
void memory::myfree_movable_fullinfo(void *p, const char *func, const char *file, int line, int movable_flag)
{
  if(Nblocks == 0)
    Terminate("no allocated blocks that could be freed");

  /* first, let's find the block */
  int nr;

  if(movable_flag)
    {
      for(nr = Nblocks - 1; nr >= 0; nr--)
        if(p == Table[nr])
          break;

      if(nr < 0)
        {
          dump_memory_table();
          Terminate("Wrong call of myfree_movable() from %s()/%s/line %d - this block has not been allocated!\n", func, file, line);
        }

      if(nr < Nblocks - 1) /* the block is not the last allocated block */
        {
          /* check that all subsequent blocks are actually movable */
          for(int i = nr + 1; i < Nblocks; i++)
            if(MovableFlag[i] == 0)
              {
                dump_memory_table();
                myflush(stdout);
                Terminate(
                    "Wrong call of myfree_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable allocated "
                    "blocks\n",
                    func, file, line, nr);
              }
        }
    }
  else
    {
      nr = Nblocks - 1;
      if(p != Table[nr])
        {
          dump_memory_table();
          Terminate("Wrong call of myfree() at %s()/%s/line %d: not the last allocated block!\n", func, file, line);
        }
    }

  if(GenericFlag[nr])
    AllocatedBytesGeneric -= BlockSize[nr];

  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  if(BasePointers[nr])
    *BasePointers[nr] = NULL;

  if(movable_flag)
    {
      ptrdiff_t offset = -BlockSize[nr];
      size_t length    = 0;

      for(int i = nr + 1; i < Nblocks; i++)
        length += BlockSize[i];

      if(nr < Nblocks - 1)
        memmove(Table[nr + 1] + offset, Table[nr + 1], length);

      for(int i = nr + 1; i < Nblocks; i++)
        {
          Table[i] += offset;
          *BasePointers[i] = *BasePointers[i] + offset;
        }

      for(int i = nr + 1; i < Nblocks; i++)
        {
          Table[i - 1]        = Table[i];
          BasePointers[i - 1] = BasePointers[i];
          BlockSize[i - 1]    = BlockSize[i];
          MovableFlag[i - 1]  = MovableFlag[i];
          GenericFlag[i - 1]  = GenericFlag[i];

          memmove(VarName + (i - 1) * MAXCHARS, VarName + i * MAXCHARS, MAXCHARS);
          memmove(FunctionName + (i - 1) * MAXCHARS, FunctionName + i * MAXCHARS, MAXCHARS);
          memmove(ParentFileName + (i - 1) * MAXCHARS, ParentFileName + i * MAXCHARS, MAXCHARS);
          memmove(FileName + (i - 1) * MAXCHARS, FileName + i * MAXCHARS, MAXCHARS);
          LineNumber[i - 1] = LineNumber[i];
        }
    }

  Nblocks -= 1;
}

/** \brief Reallocate an existing movable memory block.
 *
 *  For this operation to be successful all the blocks allocated after the block that has to be reallocated must be of movable type.
 *
 *  \param p pointer to the existing memory block to be reallocated
 *  \param n the new size of the memory block in bytes
 *  \param func name of function that has called the reallocation routine (usually given by the __FUNCTION__ macro)
 *  \param file file where the function that has called the reallocation routine resides (usually given by the __FILE__ macro)
 *  \param line line number of file where the reallocation routine was called (usually given by the __LINE__ macro)
 *  \return a pointer to the beginning of the newly allocated memory block
 */
void *memory::myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line, int movable_flag)
{
  if((n % CACHELINESIZE) > 0)
    n = (n / CACHELINESIZE + 1) * CACHELINESIZE;

  if(n < CACHELINESIZE)
    n = CACHELINESIZE;

  if(Nblocks == 0)
    Terminate("no allocated blocks that could be reallocated");

  /* first, let's find the block */
  int nr;

  if(movable_flag)
    {
      for(nr = Nblocks - 1; nr >= 0; nr--)
        if(p == Table[nr])
          break;

      if(nr < 0)
        {
          dump_memory_table();
          Terminate("Wrong call of myrealloc_movable() from %s()/%s/line %d - this block has not been allocated!\n", func, file, line);
        }

      if(nr < Nblocks - 1) /* the block is not the last allocated block */
        {
          /* check that all subsequent blocks are actually movable */
          for(int i = nr + 1; i < Nblocks; i++)
            if(MovableFlag[i] == 0)
              {
                dump_memory_table();
                Terminate(
                    "Wrong call of myrealloc_movable() from %s()/%s/line %d - behind block=%d there are subsequent non-movable "
                    "allocated blocks\n",
                    func, file, line, nr);
              }
        }
    }
  else
    {
      nr = Nblocks - 1;

      if(p != Table[nr])
        {
          dump_memory_table();
          Terminate("Wrong call of myrealloc() at %s()/%s/line %d - not the last allocated block!\n", func, file, line);
        }
    }

  if(GenericFlag[nr])
    Terminate("unexpected");

  AllocatedBytes -= BlockSize[nr];
  FreeBytes += BlockSize[nr];

  if(n > FreeBytes)
    {
      dump_memory_table();
      Terminate("At %s()/%s/line %d: Not enough memory in myremalloc_movable(n=%g MB). previous=%g FreeBytes=%g MB\n", func, file,
                line, n * TO_MBYTE_FAC, BlockSize[nr] * TO_MBYTE_FAC, FreeBytes * TO_MBYTE_FAC);
    }

  ptrdiff_t offset = n - BlockSize[nr];
  size_t length    = 0;

  for(int i = nr + 1; i < Nblocks; i++)
    length += BlockSize[i];

  if(nr < Nblocks - 1)
    memmove(Table[nr + 1] + offset, Table[nr + 1], length);

  for(int i = nr + 1; i < Nblocks; i++)
    {
      Table[i] += offset;

      *BasePointers[i] = *BasePointers[i] + offset;
    }

  FreeBytes -= n;
  AllocatedBytes += n;
  BlockSize[nr] = n;

  if(AllocatedBytes > HighMarkBytes)
    {
      HighMarkBytes = AllocatedBytes;
      dump_memory_table_buffer(HighMarkTabBuf, highmark_bufsize);
    }

  return Table[nr];
}

void memory::check_maxmemsize_setting(int maxmemsize)
{
  int errflag = 0, errflag_tot;

  if(maxmemsize > (MemoryOnNode / 1024.0 / TasksInThisNode) && RankInThisNode == 0)
    {
      char name[MPI_MAX_PROCESSOR_NAME];
      int len;
      MPI_Get_processor_name(name, &len);

      printf("On node '%s', we have %d MPI ranks and at most %g MB available. This is not enough space for MaxMemSize = %g MB\n", name,
             TasksInThisNode, MemoryOnNode / 1024.0, (double)maxmemsize);
      errflag = 1;
      fflush(stdout);
    }

  if(maxmemsize > (SharedMemoryOnNode / 1024.0 / TasksInThisNode) && RankInThisNode == 0)
    {
      char name[MPI_MAX_PROCESSOR_NAME];
      int len;
      MPI_Get_processor_name(name, &len);

      printf(
          "On node '%s', we have %d MPI ranks and at most %g MB of *shared* memory available. This is not enough space for MaxMemSize "
          "= %g MB\n",
          name, TasksInThisNode, SharedMemoryOnNode / 1024.0, (double)maxmemsize);
      errflag = 1;
      fflush(stdout);
    }

  MPI_Allreduce(&errflag, &errflag_tot, 1, MPI_INT, MPI_MAX, Communicator);

  if(errflag_tot)
    Terminate("At least one node has insufficient memory");
}
