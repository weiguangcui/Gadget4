/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file intposconvert.h
 *
 *  \brief defines a class to convert between integer coordinates and floating point positions
 */

#ifndef CONVERT_H
#define CONVERT_H

#include "gadgetconfig.h"

#include <cmath>

#include "../data/allvars.h"
#include "../data/dtypes.h"

#define MSB ((MyIntPosType)(~((MyIntPosType)(~((MyIntPosType)0)) >> ((MyIntPosType)1))))

#if defined(LONG_X_BITS)
#define HBITS_X ((MyIntPosType)(~((MyIntPosType)(~((MyIntPosType)0)) >> ((MyIntPosType)LONG_X_BITS))))
#endif

#if defined(LONG_Y_BITS)
#define HBITS_Y ((MyIntPosType)(~((MyIntPosType)(~((MyIntPosType)0)) >> ((MyIntPosType)LONG_Y_BITS))))
#endif

#if defined(LONG_Z_BITS)
#define HBITS_Z ((MyIntPosType)(~((MyIntPosType)(~((MyIntPosType)0)) >> ((MyIntPosType)LONG_Z_BITS))))
#endif

class intposconvert
{
 public:
  MyReal FacIntToCoord;
  MyReal FacCoordToInt;
  MyReal RegionLen;

#ifndef PERIODIC
  MyReal RegionCorner[3];
  MyReal RegionCenter[3];
#endif

#ifdef RANDOMIZE_DOMAINCENTER
  MyIntPosType CurrentShiftVector[3];
#endif

#ifdef EXPLICIT_VECTORIZATION
  inline Vec4d nearest_image_intpos_to_doublepos_vectorial(MyIntPosType const &a, MyIntPosType const &b0, MyIntPosType const &b1,
                                                           MyIntPosType const &b2, MyIntPosType const &b3)
  {
#if defined(GRAVITY_TALLBOX) || defined(LONG_X_BITS) || defined(LONG_Y_BITS) || defined(LONG_Z_BITS)
    Terminate("not working in this combination");
#endif

    Vec4MyIntPosType delta = a - Vec4MyIntPosType(b0, b1, b2, b3);

    Vec4MySignedIntPosType intpos = (Vec4MySignedIntPosType)delta;

    Vec4d pos = to_double(intpos);

    return pos * FacIntToCoord;
  }
#endif

  inline void constrain_intpos(MyIntPosType *pos)
  {
#ifdef PERIODIC /* restrict the position to the primary range */

#if defined(LONG_X_BITS)
    pos[0] = (pos[0] << LONG_X_BITS) >> LONG_X_BITS;
#endif

#if defined(LONG_Y_BITS)
    pos[1] = (pos[1] << LONG_Y_BITS) >> LONG_Y_BITS;
#endif

#if defined(LONG_Z_BITS)
    pos[2] = (pos[2] << LONG_Z_BITS) >> LONG_Z_BITS;
#endif

#endif
  }

  /* function to determine the nearest periodic image distance vector in MyReal, exploiting integer wrap around */
  template <typename T>
  inline void diff_intpos_to_pos(MyIntPosType *a, MyIntPosType *b, T *posdiff, offset_tuple off = 0)
  {
    if(a[0] > b[0])
      posdiff[0] = (a[0] - b[0]) * FacIntToCoord + All.BoxSize * off.n[0];
    else
      posdiff[0] = (b[0] - a[0]) * (-FacIntToCoord) + All.BoxSize * off.n[0];

    if(a[1] > b[1])
      posdiff[1] = (a[1] - b[1]) * FacIntToCoord + All.BoxSize * off.n[1];
    else
      posdiff[1] = (b[1] - a[1]) * (-FacIntToCoord) + All.BoxSize * off.n[1];

    if(a[2] > b[2])
      posdiff[2] = (a[2] - b[2]) * FacIntToCoord + All.BoxSize * off.n[2];
    else
      posdiff[2] = (b[2] - a[2]) * (-FacIntToCoord) + All.BoxSize * off.n[2];
  }

  inline MyIntPosType nearest_image_intpos_to_intpos_X(const MyIntPosType a, const MyIntPosType b)
  {
#if defined(LONG_X_BITS)
    MyIntPosType delta = (a << LONG_X_BITS) - (b << LONG_X_BITS);

    if(delta & MSB) /* tests MSB */
      {
        delta >>= LONG_X_BITS;
        delta |= HBITS_X;
      }
    else
      delta >>= LONG_X_BITS;

    return delta;
#else
    return a - b;
#endif
  }

  inline MyIntPosType nearest_image_intpos_to_intpos_Y(const MyIntPosType a, const MyIntPosType b)
  {
#if defined(LONG_Y_BITS)
    MyIntPosType delta = (a << LONG_Y_BITS) - (b << LONG_Y_BITS);

    if(delta & MSB) /* tests MSB */
      {
        delta >>= LONG_Y_BITS;
        delta |= HBITS_Y;
      }
    else
      delta >>= LONG_Y_BITS;

    return delta;
#else
    return a - b;
#endif
  }

  inline MyIntPosType nearest_image_intpos_to_intpos_Z(const MyIntPosType a, const MyIntPosType b)
  {
#if defined(LONG_Z_BITS)
    MyIntPosType delta = (a << LONG_Z_BITS) - (b << LONG_Z_BITS);

    if(delta & MSB) /* tests MSB */
      {
        delta >>= LONG_Z_BITS;
        delta |= HBITS_Z;
      }
    else
      delta >>= LONG_Z_BITS;

    return delta;
#else
    return a - b;
#endif
  }

  /* function to determine the nearest periodic image distance vector in T, exploiting integer wrap around */
  template <typename T>
  inline void nearest_image_intpos_to_pos(const MyIntPosType *const a, const MyIntPosType *const b, T *posdiff)
  {
    MyIntPosType delta[3];
    MySignedIntPosType *intpos = (MySignedIntPosType *)delta;

    /* we use all these casts here to prevent that implicit type promotions can mess this up for types shorter than int, such as
     * unsigned char */

#if defined(GRAVITY_TALLBOX) && (GRAVITY_TALLBOX == 0)
    if(a[0] >= b[0])
      {
        delta[0]   = a[0] - b[0];
        posdiff[0] = delta[0] * FacIntToCoord;
      }
    else
      {
        delta[0]   = b[0] - a[0];
        posdiff[0] = delta[0] * (-FacIntToCoord);
      }
#else

#if defined(LONG_X_BITS)
    delta[0] = (a[0] << LONG_X_BITS) - (b[0] << LONG_X_BITS);

    if(delta[0] & MSB) /* tests MSB */
      {
        delta[0] >>= LONG_X_BITS;
        delta[0] |= HBITS_X;
      }
    else
      delta[0] >>= LONG_X_BITS;
#else
    delta[0] = a[0] - b[0];
#endif
    posdiff[0] = intpos[0] * FacIntToCoord;

#endif

      /** --- **/

#if defined(GRAVITY_TALLBOX) && (GRAVITY_TALLBOX == 1)
    if(a[1] >= b[1])
      {
        delta[1]   = a[1] - b[1];
        posdiff[1] = delta[1] * FacIntToCoord;
      }
    else
      {
        delta[1]   = b[1] - a[1];
        posdiff[1] = delta[1] * (-FacIntToCoord);
      }
#else

#if defined(LONG_Y_BITS)
    delta[1]   = (a[1] << LONG_Y_BITS) - (b[1] << LONG_Y_BITS);

    if(delta[1] & MSB) /* tests MSB */
      {
        delta[1] >>= LONG_Y_BITS;
        delta[1] |= HBITS_Y;
      }
    else
      delta[1] >>= LONG_Y_BITS;
#else
    delta[1] = a[1] - b[1];
#endif

    posdiff[1] = intpos[1] * FacIntToCoord;
#endif

      /** --- **/

#if defined(GRAVITY_TALLBOX) && (GRAVITY_TALLBOX == 2)
    if(a[2] >= b[2])
      {
        delta[2]   = a[2] - b[2];
        posdiff[2] = delta[2] * FacIntToCoord;
      }
    else
      {
        delta[2]   = b[2] - a[2];
        posdiff[2] = delta[2] * (-FacIntToCoord);
      }
#else

#if defined(LONG_Z_BITS)
    delta[2]   = (a[2] << LONG_Z_BITS) - (b[2] << LONG_Z_BITS);

    if(delta[2] & MSB) /* tests MSB */
      {
        delta[2] >>= LONG_Z_BITS;
        delta[2] |= HBITS_Z;
      }
    else
      delta[2] >>= LONG_Z_BITS;
#else
    delta[2] = a[2] - b[2];
#endif
    posdiff[2] = intpos[2] * FacIntToCoord;

#endif
  }

  inline void nearest_image_intpos_to_absolute_intdist(const MyIntPosType *a, const MyIntPosType *b, MyIntPosType *delta)
  {
#if defined(LONG_X_BITS)
    delta[0] = (a[0] << LONG_X_BITS) - (b[0] << LONG_X_BITS);

    if(delta[0] & (~((~((MyIntPosType)0)) >> 1))) /* tests MSB */
      {
        delta[0] >>= LONG_X_BITS;
        delta[0] |= (~((~((MyIntPosType)0)) >> LONG_X_BITS));
      }
    else
      delta[0] >>= LONG_X_BITS;
#else
    delta[0]   = a[0] - b[0];
#endif

#if defined(LONG_Y_BITS)
    delta[1] = (a[1] << LONG_Y_BITS) - (b[1] << LONG_Y_BITS);

    if(delta[1] & (~((~((MyIntPosType)0)) >> 1))) /* tests MSB */
      {
        delta[1] >>= LONG_Y_BITS;
        delta[1] |= (~((~((MyIntPosType)0)) >> LONG_Y_BITS));
      }
    else
      delta[1] >>= LONG_Y_BITS;
#else
    delta[1]   = a[1] - b[1];
#endif

#if defined(LONG_Z_BITS)
    delta[2] = (a[2] << LONG_Z_BITS) - (b[2] << LONG_Z_BITS);

    if(delta[2] & (~((~((MyIntPosType)0)) >> 1))) /* tests MSB */
      {
        delta[2] >>= LONG_Z_BITS;
        delta[2] |= (~((~((MyIntPosType)0)) >> LONG_Z_BITS));
      }
    else
      delta[2] >>= LONG_Z_BITS;
#else
    delta[2]   = a[2] - b[2];
#endif

    if(delta[0] & (~((~((MyIntPosType)0)) >> 1))) /* tests MSB, i.e. negative if interpreted as signed int */
      delta[0] = -delta[0];

    if(delta[1] & (~((~((MyIntPosType)0)) >> 1))) /* tests MSB, i.e. negative if interpreted as signed int */
      delta[1] = -delta[1];

    if(delta[2] & (~((~((MyIntPosType)0)) >> 1))) /* tests MSB, i.e. negative if interpreted as signed int */
      delta[2] = -delta[2];
  }

#if defined(POSITIONS_IN_32BIT)
  template <typename T>
  inline MySignedIntPosType pos_to_signedintpos(T posdiff)
  {
    return std::lrint(posdiff * FacCoordToInt);
  }

  template <typename T>
  inline void pos_to_signedintpos(T *posdiff, MySignedIntPosType *intpos)
  {
    for(int j = 0; j < 3; j++)
      intpos[j] = std::lrint(posdiff[j] * FacCoordToInt);
  }
#elif defined(POSITIONS_IN_64BIT)
  template <typename T>
  inline MySignedIntPosType pos_to_signedintpos(T posdiff)
  {
    return std::llrint(posdiff * FacCoordToInt);
  }

  template <typename T>
  inline void pos_to_signedintpos(T *posdiff, MySignedIntPosType *intpos)
  {
    for(int j = 0; j < 3; j++)
      intpos[j] = std::llrint(posdiff[j] * FacCoordToInt);
  }
#else
  template <typename T>
  inline MySignedIntPosType pos_to_signedintpos(T posdiff)
  {
    return static_cast<MySignedIntPosType>(posdiff * FacCoordToInt);
  }

  template <typename T>
  inline void pos_to_signedintpos(T *posdiff, MySignedIntPosType *intpos)
  {
    for(int j = 0; j < 3; j++)
      intpos[j] = static_cast<MySignedIntPosType>(posdiff[j] * FacCoordToInt);
  }
#endif

  template <typename T>
  inline void signedintpos_to_pos(MySignedIntPosType *intpos, T *pos)
  {
    for(int j = 0; j < 3; j++)
      pos[j] = intpos[j] * FacIntToCoord;
  }

  inline double signedintpos_to_distanceorigin(MySignedIntPosType *intpos)
  {
    double pos[3];
    for(int j = 0; j < 3; j++)
      pos[j] = intpos[j] * FacIntToCoord;
    return sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  }

  template <typename T>
  inline void intpos_to_pos(MyIntPosType *intpos, T *posdiff)
  {
#ifdef RANDOMIZE_DOMAINCENTER
#ifdef PERIODIC
    MyIntPosType xyz[3];
    for(int j = 0; j < 3; j++)
      xyz[j] = intpos[j] - CurrentShiftVector[j];
    constrain_intpos(xyz);
    for(int j = 0; j < 3; j++)
      posdiff[j] = xyz[j] * FacIntToCoord;
#else
    for(int j = 0; j < 3; j++)
      posdiff[j] = (intpos[j] - CurrentShiftVector[j]) * FacIntToCoord + RegionCorner[j];
#endif
#else
#ifdef PERIODIC
    for(int j = 0; j < 3; j++)
      posdiff[j] = intpos[j] * FacIntToCoord;
#else
    for(int j = 0; j < 3; j++)
      posdiff[j] = intpos[j] * FacIntToCoord + RegionCorner[j];
#endif
#endif
  }

  inline void intpos_to_intpos(MyIntPosType *intpos, MyIntPosType *xyz)
  {
#ifdef RANDOMIZE_DOMAINCENTER
    for(int j = 0; j < 3; j++)
      xyz[j] = intpos[j] - CurrentShiftVector[j];
#else
    for(int j = 0; j < 3; j++)
      xyz[j] = intpos[j];
#endif

#ifdef PERIODIC
    constrain_intpos(xyz);
#else
    Terminate("integer coordinate output not defined if PERIODIC is not no");
#endif
  }

  template <typename T>
  inline T constrain_pos(T pos)
  {
    int rep = 0;

    while(pos < 0)
      {
        pos += ((T)(~((MyIntPosType)0)) + static_cast<T>(1.0));

        rep++;
        if(rep > MAXITER)
          Terminate("rep > MAX_ITER");
      }

    while(pos >= ((T)(~((MyIntPosType)0)) + static_cast<T>(1.0)))
      {
        pos -= ((T)(~((MyIntPosType)0)) + static_cast<T>(1.0));

        rep++;
        if(rep > MAXITER)
          Terminate("rep > MAX_ITER");
      }

    return pos;
  }

  template <typename T>
  inline void pos_to_intpos(T *posdiff, MyIntPosType *intpos)
  {
#ifdef RANDOMIZE_DOMAINCENTER
#ifdef PERIODIC
    for(int j = 0; j < 3; j++)
      {
        intpos[j] = constrain_pos(posdiff[j] * FacCoordToInt);
        intpos[j] += CurrentShiftVector[j];
      }
    constrain_intpos(intpos);
#else
    for(int j = 0; j < 3; j++)
      {
        intpos[j] = constrain_pos((posdiff[j] - RegionCorner[j]) * FacCoordToInt);
        intpos[j] += CurrentShiftVector[j];
      }
#endif
#else
#ifdef PERIODIC
    for(int j = 0; j < 3; j++)
      intpos[j] = constrain_pos(posdiff[j] * FacCoordToInt);

    constrain_intpos(intpos);
#else
    for(int j = 0; j < 3; j++)
      intpos[j] = constrain_pos((posdiff[j] - RegionCorner[j]) * FacCoordToInt);
#endif
#endif
  }
};

#endif
