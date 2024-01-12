/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file dtypes.h
 *
 *  \brief defines some custom data types used by the code
 */

#ifndef DTYPES_H
#define DTYPES_H

#include "gadgetconfig.h"

#include <stdint.h>
#include <atomic>
#include <cstddef>
#ifdef EXPLICIT_VECTORIZATION
#include "../vectorclass/vectorclass.h"
#endif

#if !defined(POSITIONS_IN_32BIT) && !defined(POSITIONS_IN_64BIT) && !defined(POSITIONS_IN_128BIT)
/* ok, nothing has been chosen as part of the configuration, then use a default value */
#ifndef DOUBLEPRECISION
#define POSITIONS_IN_32BIT
#else
#define POSITIONS_IN_64BIT
#endif
#endif

/* Exactly one of the symbols POSITIONS_IN_32BIT, POSITIONS_IN_64BIT or POSITIONS_IN_128BIT need to be fined, otherwise
 * it is desirable to get a compile time error
 */

#ifdef POSITIONS_IN_32BIT
typedef uint32_t MyIntPosType;
typedef int32_t MySignedIntPosType;
#define BITS_FOR_POSITIONS 32
#ifdef EXPLICIT_VECTORIZATION
typedef Vec4ui Vec4MyIntPosType;
typedef Vec4i Vec4MySignedIntPosType;
#endif
#endif

#ifdef POSITIONS_IN_64BIT
typedef uint64_t MyIntPosType;
typedef int64_t MySignedIntPosType;
#define BITS_FOR_POSITIONS 64
#ifdef EXPLICIT_VECTORIZATION
typedef Vec4uq Vec4MyIntPosType;
typedef Vec4q Vec4MySignedIntPosType;
#endif
#endif

#ifdef POSITIONS_IN_128BIT
typedef uint128_t MyIntPosType;
typedef int128_t MySignedIntPosType;
#define BITS_FOR_POSITIONS 128
#ifdef EXPLICIT_VECTORIZATION
#error "EXPLICIT_VECTORIZATION and POSITIONS_IN_128BIT do not work together"
#endif
#endif

#if !defined(IDS_32BIT) && !defined(IDS_48BIT) && !defined(IDS_64BIT)
#define IDS_32BIT
#endif

#ifdef IDS_32BIT
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif

#ifdef FOF_ALLOW_HUGE_GROUPLENGTH
typedef long long MyLenType;
#else
typedef int MyLenType;
#endif

#ifdef USE_SINGLEPRECISION_INTERNALLY
typedef float MyReal;
#else
typedef double MyReal;
#endif

#ifndef DOUBLEPRECISION /* default is single-precision */
typedef float MyFloat;
typedef float MyDouble;
typedef float MyNgbTreeFloat;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_FLOAT
#define H5T_NATIVE_MYFLOAT H5T_NATIVE_FLOAT
#define H5T_NATIVE_MYDOUBLE H5T_NATIVE_FLOAT
#else
#if(DOUBLEPRECISION == 2) /* mixed precision */
typedef float MyFloat;
typedef double MyDouble;
typedef float MyNgbTreeFloat;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_DOUBLE
#define H5T_NATIVE_MYFLOAT H5T_NATIVE_FLOAT
#define H5T_NATIVE_MYDOUBLE H5T_NATIVE_DOUBLE
#else /* everything double-precision */
typedef double MyFloat;
typedef double MyDouble;
typedef double MyNgbTreeFloat;
#define MPI_MYFLOAT MPI_DOUBLE
#define MPI_MYDOUBLE MPI_DOUBLE
#define H5T_NATIVE_MYFLOAT H5T_NATIVE_DOUBLE
#define H5T_NATIVE_MYDOUBLE H5T_NATIVE_DOUBLE
#endif
#endif

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
typedef long long integertime;
#define TIMEBINS 60
#define TIMEBASE                                                                                           \
  (((long long)1) << TIMEBINS) /* The simulated timespan is mapped onto the integer interval [0,TIMESPAN], \
                                *  where TIMESPAN needs to be a power of 2. */
#else
typedef int integertime;
#define TIMEBINS 29
#define TIMEBASE (1 << TIMEBINS)
#endif

#ifndef NUMBER_OF_MPI_LISTENERS_PER_NODE
#define NUMBER_OF_MPI_LISTENERS_PER_NODE 1
#endif

#ifndef MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY
#define MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY 64
#endif

#if MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY <= 32
typedef uint32_t node_bit_field;
#elif MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY <= 64
typedef uint64_t node_bit_field;
#else
#error "unsupported MAX_NUMBER_OF_RANKS_WITH_SHARED_MEMORY setting"
#endif

struct offset_tuple
{
  char n[3];

  offset_tuple() {} /* constructor */

  offset_tuple(const char x) /* constructor  */
  {
    n[0] = x;
    n[1] = x;
    n[2] = x;
  }

  offset_tuple(const char x, const char y, const char z) /* constructor  */
  {
    n[0] = x;
    n[1] = y;
    n[2] = z;
  }
};

struct location
{
  int task;
  int index;
};

inline bool operator==(const location &left, const location &right) { return left.task == right.task && left.index == right.index; }

inline bool operator!=(const location &left, const location &right) { return left.task != right.task || left.index != right.index; }

inline bool operator<(const location &left, const location &right)
{
  if(left.task < right.task)
    return true;
  else if(left.task == right.task)
    {
      if(left.index < right.index)
        return true;
      else
        return false;
    }
  else
    return false;
}

struct halotrees_table
{
  int HaloCount;
  long long FirstHalo;
  long long TreeID;
};

struct parttrees_table
{
  int ParticleCount;
  long long ParticleFirst;
  long long TreeID;
};

struct times_catalogue
{
  double Time;
  double Redshift;
};

class peanokey
{
 public:
  MyIntPosType hs, is, ls; /* 'hs'-high significance, 'is'-intermediate, 'ls'-low significance bits */
};

inline bool operator>=(const peanokey &a, const peanokey &b)
{
  if(a.hs < b.hs)
    return false;
  else if(a.hs > b.hs)
    return true;
  else if(a.is < b.is)
    return false;
  else if(a.is > b.is)
    return true;
  else if(a.ls < b.ls)
    return false;
  else
    return true;
}

inline bool operator<(const peanokey &a, const peanokey &b)
{
  if(a.hs < b.hs)
    return true;
  else if(a.hs > b.hs)
    return false;
  else if(a.is < b.is)
    return true;
  else if(a.is > b.is)
    return false;
  else if(a.ls < b.ls)
    return true;
  else
    return false;
}

inline peanokey operator+(const peanokey &a, const peanokey &b)
{
  peanokey c;

  c.ls = a.ls + b.ls;
  c.is = a.is + b.is;
  c.hs = a.hs + b.hs;

  if(c.is < a.is || c.is < b.is) /* overflow has occurred */
    {
      c.hs += 1;
    }

  if(c.ls < a.ls || c.ls < b.ls) /* overflow has occurred */
    {
      c.is += 1;
      if(c.is == 0) /* overflown again */
        c.hs += 1;
    }

  /* note: for hs we don't check for overflow explicitly as this would not be represented in the type anyhow */

  return c;
}

inline peanokey get_peanokey_offset(unsigned int j, int bits) /* this returns the peanokey for which  j << bits */
{
  peanokey key = {j, j, j};

  if(bits < BITS_FOR_POSITIONS)
    key.ls <<= bits;
  else
    key.ls = 0;

  int is_bits = bits - BITS_FOR_POSITIONS;

  if(is_bits <= -BITS_FOR_POSITIONS)
    key.is = 0;
  else if(is_bits <= 0)
    key.is >>= -is_bits;
  else if(is_bits < BITS_FOR_POSITIONS)
    key.is <<= is_bits;
  else
    key.is = 0;

  int hs_bits = bits - 2 * BITS_FOR_POSITIONS;

  if(hs_bits <= -BITS_FOR_POSITIONS)
    key.hs = 0;
  else if(hs_bits <= 0)
    key.hs >>= -hs_bits;
  else if(hs_bits < BITS_FOR_POSITIONS)
    key.hs <<= hs_bits;
  else
    key.hs = 0;

  return key;
}

enum mysnaptype
{
  NORMAL_SNAPSHOT,
  MOST_BOUND_PARTICLE_SNAPHOT,
  MOST_BOUND_PARTICLE_SNAPHOT_REORDERED
};

enum restart_options
{
  RST_BEGIN,
  RST_RESUME,
  RST_STARTFROMSNAP,
  RST_FOF,
  RST_POWERSPEC,
  RST_CONVERTSNAP,
  RST_CREATEICS,
  RST_CALCDESC,
  RST_MAKETREES,
  RST_IOBANDWIDTH,
  RST_LCREARRANGE,
  RST_SNPREARRANGE
};

struct data_partlist
{
  int Task;  /** The task the item was exported to. */
  int Index; /** The particle index of the item on the sending task. */
};

struct thread_data
{
  int Nexport;
  int NexportNodes;

  double Interactions; /*!< The total cost of the particles/nodes processed by each thread */

  double Ewaldcount; /*!< The total cost for the Ewald correction per thread */
  int FirstExec;     /*!< Keeps track, if a given thread executes the gravity_primary_loop() for the first time */

  size_t ExportSpace;
  size_t InitialSpace;
  size_t ItemSize;

  int *P_CostCount;
  int *TreePoints_CostCount;
  int *Node_CostCount;

  data_partlist *PartList;
  int *Ngblist;
  int *Shmranklist;
  int *Exportflag;
};

template <typename T>
struct copyable_atomic : std::atomic<T>
{
  using std::atomic<T>::atomic;

  copyable_atomic(const copyable_atomic &ca) noexcept : std::atomic<T>(ca.load()) {}

  using std::atomic<T>::operator=;

  copyable_atomic &operator=(const copyable_atomic &other) noexcept
  {
    this->store(other.load());
    return *this;
  }
};

struct copyable_atomic_flag : std::atomic_flag
{
  using std::atomic_flag::atomic_flag;

  copyable_atomic_flag(const copyable_atomic_flag &ca) noexcept { this->clear(); }

  copyable_atomic_flag &operator=(const copyable_atomic_flag &other) noexcept
  {
    this->clear();
    return *this;
  }
};

#ifdef LONG_X_BITS
#define LONG_X (1 << (LONG_X_BITS))
#define MAX_LONG_X_BITS LONG_X_BITS
#else
#define LONG_X 1
#define MAX_LONG_X_BITS 0
#endif

#ifdef LONG_Y_BITS
#define LONG_Y (1 << (LONG_Y_BITS))
#define MAX_LONG_Y_BITS LONG_Y_BITS
#else
#define LONG_Y 1
#define MAX_LONG_Y_BITS 0
#endif

#ifdef LONG_Z_BITS
#define LONG_Z (1 << (LONG_Z_BITS))
#define MAX_LONG_Z_BITS LONG_Z_BITS
#else
#define LONG_Z 1
#define MAX_LONG_Z_BITS 0
#endif

#define LONG_BITS_MAX(A, B) (((A) > (B)) ? (A) : (B))

#define LEVEL_ALWAYS_OPEN LONG_BITS_MAX(MAX_LONG_X_BITS, LONG_BITS_MAX(MAX_LONG_Y_BITS, MAX_LONG_Z_BITS))

#ifdef GRAVITY_TALLBOX

#if(GRAVITY_TALLBOX == 0)
#define DBX 2
#define DBX_EXTRA 6
#define BOXX (1.0 / LONG_Y)
#define BOXY (1.0 / LONG_Z)
#else
#define DBX 1
#define DBX_EXTRA 0
#endif

#if(GRAVITY_TALLBOX == 1)
#define DBY 2
#define DBY_EXTRA 6
#define BOXX (1.0 / LONG_X)
#define BOXY (1.0 / LONG_Z)
#else
#define DBY 1
#define DBY_EXTRA 0
#endif

#if(GRAVITY_TALLBOX == 2)
#define DBZ 2
#define DBZ_EXTRA 6
#define BOXX (1.0 / LONG_X)
#define BOXY (1.0 / LONG_Y)
#else
#define DBZ 1
#define DBZ_EXTRA 0
#endif

#else

#define DBX 1
#define DBY 1
#define DBZ 1
#define DBX_EXTRA 0
#define DBY_EXTRA 0
#define DBZ_EXTRA 0
#endif

#endif
