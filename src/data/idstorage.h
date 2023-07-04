/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file idstorage.h
 *
 *  \brief defines a class that we use to store and manipulate the particle IDs
 */

#ifndef IDSTORAGE_H
#define IDSTORAGE_H

#include "gadgetconfig.h"

#include <climits>

#if !defined(IDS_48BIT)
#define ID_MSB ((MyIDType)(~((MyIDType)(~((MyIDType)0)) >> ((MyIDType)1))))
#define ID_MSK ((MyIDType)(((MyIDType)(~((MyIDType)0)) >> ((MyIDType)1))))
#define HALONR_MAX ((MyIDType)(((MyIDType)(~((MyIDType)0)) >> ((MyIDType)1))))
#define ID_MAX ID_MSK
#else
#define ID_MSB ((unsigned short)(~((unsigned short)(~((unsigned short)0)) >> ((unsigned short)1))))
#define ID_MSK ((unsigned short)(((unsigned short)(~((unsigned short)0)) >> ((unsigned short)1))))
#define HALONR_MAX ((MyIDType)(((MyIDType)(~((MyIDType)0)) >> ((MyIDType)17))))
#define ID_MAX (((long long)ID_MSK) << 32LL)
#endif

/* used to store a subhalo len in an approximate (quite accurate) way in just two bytes */
struct approxlen
{
#define ALEN_MAX 1000000000.0
#define ALEN_MIN 10.0

 private:
  unsigned short alen;

 public:
  inline void set(long long size)
  {
    double l = log(size / ALEN_MIN) / log(ALEN_MAX / ALEN_MIN) * (USHRT_MAX - 1) + 1;

    if(l < 1)
      alen = 0;
    else if(l > USHRT_MAX)
      alen = USHRT_MAX;
    else
      alen = (unsigned short)(l + 0.5);
  }

  inline long long get(void)
  {
    // relatice accuracy of this encoding is about ~0.00012 for particle numbers between 10 to 10^9
    if(alen == 0)
      return 0;
    else
      {
        return (long long)(ALEN_MIN * exp((alen - 1.0) / (USHRT_MAX - 1) * log(ALEN_MAX / ALEN_MIN)) + 0.5);
      }
  }
};

struct compactrank_t
{
 private:
  unsigned char rank;

 public:
  inline void set(MyLenType nr)
  {
    if(nr > UCHAR_MAX)
      nr = UCHAR_MAX;
    rank = nr;
  }

  inline unsigned char get(void) { return rank; }
};

class MyIDStorage
{
 private:
#if !defined(IDS_48BIT)
  MyIDType id;
#else
  unsigned short id[3];
#endif

 public:
  inline MyIDType get(void) const
  {
#if !defined(IDS_48BIT)
    return id & ID_MSK;
#else
    return (((MyIDType)(id[0] & ID_MSK)) << 32) + (((MyIDType)id[1]) << 16) + id[2];
#endif
  }

  inline void set(MyIDType ID)
  {
#if !defined(IDS_48BIT)
    id = ID;
#else
    id[2] = (unsigned short)ID;
    id[1] = (unsigned short)(ID >> 16);
    id[0] = (unsigned short)(ID >> 32);
#endif
  }

  inline void mark_as_formerly_most_bound(void)
  {
    /* we set the most significant bit */
#if !defined(IDS_48BIT)
    id |= ID_MSB;
#else
    id[0] |= ID_MSB;
#endif
  }

  inline bool is_previously_most_bound(void)
  {
    /* we set the most significant bit */
#if defined(IDS_48BIT)
    if(id[0] & ID_MSB)
      return true;
#else
    if(id & ID_MSB)
      return true;
#endif
    return false;
  }
};

class MyHaloNrType : public MyIDStorage
{
 public:
  inline MyHaloNrType &operator+=(const long long &x)
  {
    set(get() + x);
    return *this;
  }
};

inline bool operator<(const MyHaloNrType &left, const MyHaloNrType &right) { return left.get() < right.get(); }

inline bool operator>(const MyHaloNrType &left, const MyHaloNrType &right) { return left.get() > right.get(); }

inline bool operator!=(const MyHaloNrType &left, const MyHaloNrType &right) { return left.get() != right.get(); }

inline bool operator==(const MyHaloNrType &left, const MyHaloNrType &right) { return left.get() == right.get(); }

#endif /* IDSTORAGE_H */
