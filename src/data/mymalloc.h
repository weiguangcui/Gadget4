/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file mymalloc.h
 *
 *  \brief declares class that organizes the dynamic memory allocation
 */

#ifndef MYMALLOC_H
#define MYMALLOC_H

#include "gadgetconfig.h"

#include <stdio.h>

#define CACHELINESIZE 64

#define MAXBLOCKS 5000
#define MAXCHARS 40

#define LOC __FUNCTION__, __FILE__, __LINE__
#define MMM(x, y) (x, #x, y, __FUNCTION__, __FILE__, __LINE__)
#define DDD(x) (x, __FUNCTION__, __FILE__, __LINE__)

#define mymalloc(x, y) mymalloc_movable_fullinfo(NULL, x, y, __FUNCTION__, __FILE__, __LINE__, 0, 0, NULL)
#define mymalloc_clear(x, y) mymalloc_movable_fullinfo(NULL, x, y, __FUNCTION__, __FILE__, __LINE__, 0, 1, NULL)
#define mymalloc_g(x, y) mymalloc_movable_fullinfo(NULL, x, y, __FUNCTION__, __FILE__, __LINE__, 0, 0, callorigin)

#define mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__, 1, 0, NULL)
#define mymalloc_movable_clear(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__, 1, 1, NULL)
#define mymalloc_movable_g(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__, 1, 0, callorigin)

#define myfree(x) myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__, 0)
#define myfree_movable(x) myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__, 1)

#define myrealloc(x, y) myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__, 0)
#define myrealloc_movable(x, y) myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__, 1)

#include "../mpi_utils/setcomm.h"

class memory : public setcomm
{
 public:
  memory() : setcomm("delayed init") {}

  size_t AllocatedBytes;
  size_t FreeBytes;
  char *Base; /**< Base pointer (initial memory address) of the stack. */

  void mymalloc_init(int maxmemsize, enum restart_options restartflag);

  void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line,
                                  int movable_flag, int clear_flag, char *originflag);

  void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line, int movable_flag);

  void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line, int movable_flag);

  void *myfree_query_last_block(void);

  size_t roundup_to_multiple_of_cacheline_size(size_t n);

  void report_detailed_memory_usage_of_largest_task(void);

  void check_maxmemsize_setting(int maxmemsize);

  int myMPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
  int myMPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);

  inline double getAllocatedBytesInMB(void) { return AllocatedBytes * TO_MBYTE_FAC; }

  template <typename T>
  inline T *alloc(T *&ptr, const char *varname, size_t n, const char *func, const char *file, int linenr)
  {
    return static_cast<T *>(mymalloc_movable_fullinfo(&ptr, varname, n * sizeof(T), func, file, linenr, 0, 0, NULL));
  }

  template <typename T>
  inline T *alloc_movable(T *&ptr, const char *varname, size_t n, const char *func, const char *file, int linenr)
  {
    return static_cast<T *>(mymalloc_movable_fullinfo(&ptr, varname, n * sizeof(T), func, file, linenr, 1, 0, NULL));
  }

  template <typename T>
  inline T *realloc(T *&ptr, const char *varname, size_t n, const char *func, const char *file, int linenr)
  {
    return static_cast<T *>(myrealloc_movable_fullinfo(&ptr, n * sizeof(T), func, file, linenr, 0));
  }

  template <typename T>
  inline T *realloc_movable(T *&ptr, const char *varname, size_t n, const char *func, const char *file, int linenr)
  {
    return static_cast<T *>(myrealloc_movable_fullinfo(&ptr, n * sizeof(T), func, file, linenr, 1));
  }

  void dealloc(void *ptr, const char *func, const char *file, int linenr) { myfree_movable_fullinfo(ptr, func, file, linenr, 0); }

  void dealloc_movable(void *ptr, const char *func, const char *file, int linenr)
  {
    myfree_movable_fullinfo(ptr, func, file, linenr, 1);
  }

  void dump_memory_table(void);

 private:
  size_t AllocatedBytesGeneric;

  size_t HighMarkBytes;
  size_t HighMarkBytesWithoutGeneric;

  double OldGlobHighMarkMB;
  double OldGlobHighMarkMBWithoutGeneric;

  FILE *FdMemory; /**< file handle for memory.txt log-file. */

  size_t TotBytes; /**< The total dimension (in bytes) of dynamic memory available to the current task. */
  int Nblocks;     /**< The current number of allocated memory blocks. */

  char **Table;         /**< Table containing the initial addresses of the allocated memory blocks.*/
  size_t *BlockSize;    /**< Array containing the size (in bytes) of all the allocated memory blocks. */
  char *MovableFlag;    /**< Identifies whether a block is movable. */
  char *GenericFlag;    /**< Identifies whether a block has been identified in the generic allocation routines. */
  char ***BasePointers; /**< Base pointers containing the initial addresses of movable memory blocks */
  char *VarName;        /**< The name of the variable with which the block has been allocated. */
  char *FunctionName;   /**< The function name that has allocated the memory block. */
  char *ParentFileName; /**< The location from which the generich routines were called */
  char *FileName;       /**< The file name where the function that has allocated the block is called. */
  int *LineNumber;      /**< The line number in FileName where the function that allocated the block has been called. */
  char *HighMarkTabBuf; /**< This is a buffer that holds the log-file output corresponding to the largest memory use that has occurred
                           on this task */
  char *HighMarkTabBufWithoutGeneric; /**< This is a buffer that holds the log-file output corresponding to the largest memory use that
                                         has occurred on this task */
  enum restart_options RestartFlag;

  int highmark_bufsize;

  int dump_memory_table_buffer(char *p, int bufsize);

  void report_memory_usage(int rank, char *tabbuf);
};

extern memory Mem;

#endif
