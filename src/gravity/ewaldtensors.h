/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file ewaldtensors.h
 *
 *  \brief defines derivative tensors with cubic symmetry for Ewald correction
 */

#ifndef GRAVITY_EWALDTENSORS_H
#define GRAVITY_EWALDTENSORS_H

#include "gadgetconfig.h"

#include "../data/symtensors.h"

// derivative tensors for Ewald correction - they have few independent elements due to cubic symmetry

template <typename T>
class ewaldtensor0
{
 public:
  T x0;

  // constructor
  ewaldtensor0() {}

  // constructor
  ewaldtensor0(const T x) { x0 = x; }

  inline ewaldtensor0 &operator+=(const ewaldtensor0 &right)
  {
    x0 += right.x0;

    return *this;
  }
};

template <typename T>
class ewaldtensor2
{
 public:
  T XX;

  // constructor
  ewaldtensor2() {}

  // constructor
  ewaldtensor2(const T x) { XX = x; }

  inline ewaldtensor2 &operator+=(const ewaldtensor2 &right)
  {
    XX += right.XX;

    return *this;
  }
};

template <typename T>
class ewaldtensor4
{
 public:
  T XXXX;
  T XXYY;

  // constructor
  ewaldtensor4() {}

  // constructor
  ewaldtensor4(const T x)
  {
    XXXX = x;
    XXYY = x;
  }

  inline ewaldtensor4 &operator+=(const ewaldtensor4 &right)
  {
    XXXX += right.XXXX;
    XXYY += right.XXYY;

    return *this;
  }
};

template <typename T>
class ewaldtensor6
{
 public:
  T XXXXXX;
  T XXXXYY;
  T XXYYZZ;

  // constructor
  ewaldtensor6() {}

  // constructor
  ewaldtensor6(const T x)
  {
    XXXXXX = x;
    XXXXYY = x;
    XXYYZZ = x;
  }

  inline ewaldtensor6 &operator+=(const ewaldtensor6 &right)
  {
    XXXXXX += right.XXXXXX;
    XXXXYY += right.XXXXYY;
    XXYYZZ += right.XXYYZZ;

    return *this;
  }
};

template <typename T>
class ewaldtensor8
{
 public:
  T XXXXXXXX;
  T XXXXXXYY;
  T XXXXYYYY;
  T XXXXYYZZ;

  // constructor
  ewaldtensor8() {}

  // constructor
  ewaldtensor8(const T x)
  {
    XXXXXXXX = x;
    XXXXXXYY = x;
    XXXXYYYY = x;
    XXXXYYZZ = x;
  }

  inline ewaldtensor8 &operator+=(const ewaldtensor8 &right)
  {
    XXXXXXXX += right.XXXXXXXX;
    XXXXXXYY += right.XXXXXXYY;
    XXXXYYYY += right.XXXXYYYY;
    XXXXYYZZ += right.XXXXYYZZ;

    return *this;
  }
};

template <typename T>
class ewaldtensor10
{
 public:
  T XXXXXXXXXX;
  T XXXXXXXXYY;
  T XXXXXXYYYY;
  T XXXXXXYYZZ;
  T XXXXYYYYZZ;

  // constructor
  ewaldtensor10() {}

  // constructor
  ewaldtensor10(const T x)
  {
    XXXXXXXXXX = x;
    XXXXXXXXYY = x;
    XXXXXXYYYY = x;
    XXXXXXYYZZ = x;
    XXXXYYYYZZ = x;
  }

  inline ewaldtensor10 &operator+=(const ewaldtensor10 &right)
  {
    XXXXXXXXXX += right.XXXXXXXXXX;
    XXXXXXXXYY += right.XXXXXXXXYY;
    XXXXXXYYYY += right.XXXXXXYYYY;
    XXXXXXYYZZ += right.XXXXXXYYZZ;
    XXXXYYYYZZ += right.XXXXYYYYZZ;

    return *this;
  }
};

// multiply with a scalar factor
template <typename T>
inline ewaldtensor6<T> operator*(const double fac, const ewaldtensor6<T> &S)
{
  ewaldtensor6<T> res;

  res.XXXXXX = fac * S.XXXXXX;
  res.XXXXYY = fac * S.XXXXYY;
  res.XXYYZZ = fac * S.XXYYZZ;

  return res;
}

// multiply with a scalar factor
template <typename T>
inline ewaldtensor8<T> operator*(const double fac, const ewaldtensor8<T> &S)
{
  ewaldtensor8<T> res;

  res.XXXXXXXX = fac * S.XXXXXXXX;
  res.XXXXXXYY = fac * S.XXXXXXYY;
  res.XXXXYYYY = fac * S.XXXXYYYY;
  res.XXXXYYZZ = fac * S.XXXXYYZZ;

  return res;
}

// multiply with a scalar factor
template <typename T>
inline ewaldtensor10<T> operator*(const double fac, const ewaldtensor10<T> &S)
{
  ewaldtensor10<T> res;

  res.XXXXXXXXXX = fac * S.XXXXXXXXXX;
  res.XXXXXXXXYY = fac * S.XXXXXXXXYY;
  res.XXXXXXYYYY = fac * S.XXXXXXYYYY;
  res.XXXXXXYYZZ = fac * S.XXXXXXYYZZ;
  res.XXXXYYYYZZ = fac * S.XXXXYYYYZZ;

  return res;
}

// contract a 2-ewaldtensor with a symmetric 2-tensor to yield a scalar
template <typename T>
inline T operator*(const ewaldtensor2<T> &S, const symtensor2<T> &D)
{
  T res = S.XX * D.da[qZZ] + S.XX * D.da[qYY] + S.XX * D.da[qXX];

  return res;
}

// contract a 4-ewaldtensor with a symmetric 4-tensor to yield a scalar
template <typename T>
inline T operator*(const ewaldtensor4<T> &S, const symtensor4<T> &D)
{
  T res = S.XXXX * D.da[sZZZZ] + 6.0 * S.XXYY * D.da[sYYZZ] + 6.0 * S.XXYY * D.da[sXXZZ] + S.XXXX * D.da[sYYYY] +
          6.0 * S.XXYY * D.da[sXXYY] + S.XXXX * D.da[sXXXX];

  return res;
}

// contract a 4-ewaldtensor with a symmetric 3-tensor to yield a vector
template <typename T>
inline vector<T> operator*(const ewaldtensor4<T> &S, const symtensor3<T> &D)
{
  vector<T> res;

  res.da[0] = S.XXXX * D.da[dZZZ] + 3.0 * S.XXYY * D.da[dYYZ] + 3.0 * S.XXYY * D.da[dXXZ];

  res.da[1] = 3.0 * S.XXYY * D.da[dYZZ] + S.XXXX * D.da[dYYY] + 3.0 * S.XXYY * D.da[dXXY];

  res.da[2] = 3.0 * S.XXYY * D.da[dXZZ] + 3.0 * S.XXYY * D.da[dXYY] + S.XXXX * D.da[dXXX];

  return res;
}

// contract a 6-ewaldtensor with a symmetric 5-tensor to yield a vector
template <typename T>
inline vector<T> operator*(const ewaldtensor6<T> &S, const symtensor5<T> &D)
{
  vector<T> res;

  res.da[0] = S.XXXXXX * D.da[rZZZZZ] + 10.0 * S.XXXXYY * D.da[rYYZZZ] + 10.0 * S.XXXXYY * D.da[rXXZZZ] +
              5.0 * S.XXXXYY * D.da[rYYYYZ] + 30.0 * S.XXYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXYY * D.da[rXXXXZ];

  res.da[1] = 5.0 * S.XXXXYY * D.da[rYZZZZ] + 10.0 * S.XXXXYY * D.da[rYYYZZ] + 30.0 * S.XXYYZZ * D.da[rXXYZZ] +
              S.XXXXXX * D.da[rYYYYY] + 10.0 * S.XXXXYY * D.da[rXXYYY] + 5.0 * S.XXXXYY * D.da[rXXXXY];

  res.da[2] = 5.0 * S.XXXXYY * D.da[rXZZZZ] + 30.0 * S.XXYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXYY * D.da[rXXXZZ] +
              5.0 * S.XXXXYY * D.da[rXYYYY] + 10.0 * S.XXXXYY * D.da[rXXXYY] + S.XXXXXX * D.da[rXXXXX];

  return res;
}

// contract a 2-ewaldtensor with a 0-tensor to yield a 2-tensor
template <typename T>
inline symtensor2<T> operator*(const ewaldtensor2<T> &S, const T &Q0)
{
  symtensor2<T> res;

  res.da[qXX] = S.XX * Q0;

  res.da[qXY] = 0.0;

  res.da[qXZ] = 0.0;

  res.da[qYY] = S.XX * Q0;

  res.da[qYZ] = 0.0;

  res.da[qZZ] = S.XX * Q0;

  return res;
}

// contract a 4-ewaldtensor with a symmetric 2-tensor to yield a 2-tensor
template <typename T>
inline symtensor2<T> operator*(const ewaldtensor4<T> &S, const symtensor2<T> &D)
{
  symtensor2<T> res;

  res.da[qXX] = S.XXXX * D.da[qZZ] + S.XXYY * D.da[qYY] + S.XXYY * D.da[qXX];

  res.da[qXY] = 2.0 * S.XXYY * D.da[qYZ];

  res.da[qXZ] = 2.0 * S.XXYY * D.da[qXZ];

  res.da[qYY] = S.XXYY * D.da[qZZ] + S.XXXX * D.da[qYY] + S.XXYY * D.da[qXX];

  res.da[qYZ] = 2.0 * S.XXYY * D.da[qXY];

  res.da[qZZ] = S.XXYY * D.da[qZZ] + S.XXYY * D.da[qYY] + S.XXXX * D.da[qXX];

  return res;
}

// contract a 6-ewaldtensor with a symmetric 4-tensor to yield a 3-tensor
template <typename T>
inline symtensor2<T> operator*(const ewaldtensor6<T> &S, const symtensor4<T> &D)
{
  symtensor2<T> res;

  res.da[qXX] = S.XXXXXX * D.da[sZZZZ] + 6.0 * S.XXXXYY * D.da[sYYZZ] + 6.0 * S.XXXXYY * D.da[sXXZZ] + S.XXXXYY * D.da[sYYYY] +
                6.0 * S.XXYYZZ * D.da[sXXYY] + S.XXXXYY * D.da[sXXXX];

  res.da[qXY] = 4.0 * S.XXXXYY * D.da[sYZZZ] + 4.0 * S.XXXXYY * D.da[sYYYZ] + 12.0 * S.XXYYZZ * D.da[sXXYZ];

  res.da[qXZ] = 4.0 * S.XXXXYY * D.da[sXZZZ] + 12.0 * S.XXYYZZ * D.da[sXYYZ] + 4.0 * S.XXXXYY * D.da[sXXXZ];

  res.da[qYY] = S.XXXXYY * D.da[sZZZZ] + 6.0 * S.XXXXYY * D.da[sYYZZ] + 6.0 * S.XXYYZZ * D.da[sXXZZ] + S.XXXXXX * D.da[sYYYY] +
                6.0 * S.XXXXYY * D.da[sXXYY] + S.XXXXYY * D.da[sXXXX];

  res.da[qYZ] = 12.0 * S.XXYYZZ * D.da[sXYZZ] + 4.0 * S.XXXXYY * D.da[sXYYY] + 4.0 * S.XXXXYY * D.da[sXXXY];

  res.da[qZZ] = S.XXXXYY * D.da[sZZZZ] + 6.0 * S.XXYYZZ * D.da[sYYZZ] + 6.0 * S.XXXXYY * D.da[sXXZZ] + S.XXXXYY * D.da[sYYYY] +
                6.0 * S.XXXXYY * D.da[sXXYY] + S.XXXXXX * D.da[sXXXX];

  return res;
}

// contract a 6-ewaldtensor with a symmetric 3-tensor to yield a 3-tensor
template <typename T>
inline symtensor3<T> operator*(const ewaldtensor6<T> &S, const symtensor3<T> &D)
{
  symtensor3<T> res;

  res.da[dXXX] = S.XXXXXX * D.da[dZZZ] + 3.0 * S.XXXXYY * D.da[dYYZ] + 3.0 * S.XXXXYY * D.da[dXXZ];

  res.da[dXXY] = 3.0 * S.XXXXYY * D.da[dYZZ] + S.XXXXYY * D.da[dYYY] + 3.0 * S.XXYYZZ * D.da[dXXY];

  res.da[dXXZ] = 3.0 * S.XXXXYY * D.da[dXZZ] + 3.0 * S.XXYYZZ * D.da[dXYY] + S.XXXXYY * D.da[dXXX];

  res.da[dXYY] = S.XXXXYY * D.da[dZZZ] + 3.0 * S.XXXXYY * D.da[dYYZ] + 3.0 * S.XXYYZZ * D.da[dXXZ];

  res.da[dXYZ] = 6.0 * S.XXYYZZ * D.da[dXYZ];

  res.da[dXZZ] = S.XXXXYY * D.da[dZZZ] + 3.0 * S.XXYYZZ * D.da[dYYZ] + 3.0 * S.XXXXYY * D.da[dXXZ];

  res.da[dYYY] = 3.0 * S.XXXXYY * D.da[dYZZ] + S.XXXXXX * D.da[dYYY] + 3.0 * S.XXXXYY * D.da[dXXY];

  res.da[dYYZ] = 3.0 * S.XXYYZZ * D.da[dXZZ] + 3.0 * S.XXXXYY * D.da[dXYY] + S.XXXXYY * D.da[dXXX];

  res.da[dYZZ] = 3.0 * S.XXYYZZ * D.da[dYZZ] + S.XXXXYY * D.da[dYYY] + 3.0 * S.XXXXYY * D.da[dXXY];

  res.da[dZZZ] = 3.0 * S.XXXXYY * D.da[dXZZ] + 3.0 * S.XXXXYY * D.da[dXYY] + S.XXXXXX * D.da[dXXX];

  return res;
}

// contract a 8-ewaldtensor with a symmetric 5-tensor to yield a 3-tensor
template <typename T>
inline symtensor3<T> operator*(const ewaldtensor8<T> &S, const symtensor5<T> &D)
{
  symtensor3<T> res;

  res.da[dXXX] = S.XXXXXXXX * D.da[rZZZZZ] + 10.0 * S.XXXXXXYY * D.da[rYYZZZ] + 10.0 * S.XXXXXXYY * D.da[rXXZZZ] +
                 5.0 * S.XXXXYYYY * D.da[rYYYYZ] + 30.0 * S.XXXXYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXYYYY * D.da[rXXXXZ];

  res.da[dXXY] = 5.0 * S.XXXXXXYY * D.da[rYZZZZ] + 10.0 * S.XXXXYYYY * D.da[rYYYZZ] + 30.0 * S.XXXXYYZZ * D.da[rXXYZZ] +
                 S.XXXXXXYY * D.da[rYYYYY] + 10.0 * S.XXXXYYZZ * D.da[rXXYYY] + 5.0 * S.XXXXYYZZ * D.da[rXXXXY];

  res.da[dXXZ] = 5.0 * S.XXXXXXYY * D.da[rXZZZZ] + 30.0 * S.XXXXYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXYYYY * D.da[rXXXZZ] +
                 5.0 * S.XXXXYYZZ * D.da[rXYYYY] + 10.0 * S.XXXXYYZZ * D.da[rXXXYY] + S.XXXXXXYY * D.da[rXXXXX];

  res.da[dXYY] = S.XXXXXXYY * D.da[rZZZZZ] + 10.0 * S.XXXXYYYY * D.da[rYYZZZ] + 10.0 * S.XXXXYYZZ * D.da[rXXZZZ] +
                 5.0 * S.XXXXXXYY * D.da[rYYYYZ] + 30.0 * S.XXXXYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXYYZZ * D.da[rXXXXZ];

  res.da[dXYZ] = 20.0 * S.XXXXYYZZ * D.da[rXYZZZ] + 20.0 * S.XXXXYYZZ * D.da[rXYYYZ] + 20.0 * S.XXXXYYZZ * D.da[rXXXYZ];

  res.da[dXZZ] = S.XXXXXXYY * D.da[rZZZZZ] + 10.0 * S.XXXXYYZZ * D.da[rYYZZZ] + 10.0 * S.XXXXYYYY * D.da[rXXZZZ] +
                 5.0 * S.XXXXYYZZ * D.da[rYYYYZ] + 30.0 * S.XXXXYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXXXYY * D.da[rXXXXZ];

  res.da[dYYY] = 5.0 * S.XXXXYYYY * D.da[rYZZZZ] + 10.0 * S.XXXXXXYY * D.da[rYYYZZ] + 30.0 * S.XXXXYYZZ * D.da[rXXYZZ] +
                 S.XXXXXXXX * D.da[rYYYYY] + 10.0 * S.XXXXXXYY * D.da[rXXYYY] + 5.0 * S.XXXXYYYY * D.da[rXXXXY];

  res.da[dYYZ] = 5.0 * S.XXXXYYZZ * D.da[rXZZZZ] + 30.0 * S.XXXXYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXYYZZ * D.da[rXXXZZ] +
                 5.0 * S.XXXXXXYY * D.da[rXYYYY] + 10.0 * S.XXXXYYYY * D.da[rXXXYY] + S.XXXXXXYY * D.da[rXXXXX];

  res.da[dYZZ] = 5.0 * S.XXXXYYZZ * D.da[rYZZZZ] + 10.0 * S.XXXXYYZZ * D.da[rYYYZZ] + 30.0 * S.XXXXYYZZ * D.da[rXXYZZ] +
                 S.XXXXXXYY * D.da[rYYYYY] + 10.0 * S.XXXXYYYY * D.da[rXXYYY] + 5.0 * S.XXXXXXYY * D.da[rXXXXY];

  res.da[dZZZ] = 5.0 * S.XXXXYYYY * D.da[rXZZZZ] + 30.0 * S.XXXXYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXXXYY * D.da[rXXXZZ] +
                 5.0 * S.XXXXYYYY * D.da[rXYYYY] + 10.0 * S.XXXXXXYY * D.da[rXXXYY] + S.XXXXXXXX * D.da[rXXXXX];

  return res;
}

// contract a 4-ewaldtensor with a symmetric 0-tensor to yield a 4-tensor
template <typename T>
inline symtensor4<T> operator*(const ewaldtensor4<T> &S, const T &Q0)
{
  symtensor4<T> res;

  res.da[sXXXX] = S.XXXX * Q0;

  res.da[sXXXY] = 0.0;

  res.da[sXXXZ] = 0.0;

  res.da[sXXYY] = S.XXYY * Q0;

  res.da[sXXYZ] = 0.0;

  res.da[sXXZZ] = S.XXYY * Q0;

  res.da[sXYYY] = 0.0;

  res.da[sXYYZ] = 0.0;

  res.da[sXYZZ] = 0.0;

  res.da[sXZZZ] = 0.0;

  res.da[sYYYY] = S.XXXX * Q0;

  res.da[sYYYZ] = 0.0;

  res.da[sYYZZ] = S.XXYY * Q0;

  res.da[sYZZZ] = 0.0;

  res.da[sZZZZ] = S.XXXX * Q0;

  return res;
}

// contract a 6-ewaldtensor with a symmetric 2-tensor to yield a 4-tensor
template <typename T>
inline symtensor4<T> operator*(const ewaldtensor6<T> &S, const symtensor2<T> &D)
{
  symtensor4<T> res;

  res.da[sXXXX] = S.XXXXXX * D.da[qZZ] + S.XXXXYY * D.da[qYY] + S.XXXXYY * D.da[qXX];

  res.da[sXXXY] = 2.0 * S.XXXXYY * D.da[qYZ];

  res.da[sXXXZ] = 2.0 * S.XXXXYY * D.da[qXZ];

  res.da[sXXYY] = S.XXXXYY * D.da[qZZ] + S.XXXXYY * D.da[qYY] + S.XXYYZZ * D.da[qXX];

  res.da[sXXYZ] = 2.0 * S.XXYYZZ * D.da[qXY];

  res.da[sXXZZ] = S.XXXXYY * D.da[qZZ] + S.XXYYZZ * D.da[qYY] + S.XXXXYY * D.da[qXX];

  res.da[sXYYY] = 2.0 * S.XXXXYY * D.da[qYZ];

  res.da[sXYYZ] = 2.0 * S.XXYYZZ * D.da[qXZ];

  res.da[sXYZZ] = 2.0 * S.XXYYZZ * D.da[qYZ];

  res.da[sXZZZ] = 2.0 * S.XXXXYY * D.da[qXZ];

  res.da[sYYYY] = S.XXXXYY * D.da[qZZ] + S.XXXXXX * D.da[qYY] + S.XXXXYY * D.da[qXX];

  res.da[sYYYZ] = 2.0 * S.XXXXYY * D.da[qXY];

  res.da[sYYZZ] = S.XXYYZZ * D.da[qZZ] + S.XXXXYY * D.da[qYY] + S.XXXXYY * D.da[qXX];

  res.da[sYZZZ] = 2.0 * S.XXXXYY * D.da[qXY];

  res.da[sZZZZ] = S.XXXXYY * D.da[qZZ] + S.XXXXYY * D.da[qYY] + S.XXXXXX * D.da[qXX];
  return res;
}

// contract a 8-ewaldtensor with a symmetric 4-tensor to yield a 4-tensor
template <typename T>
inline symtensor4<T> operator*(const ewaldtensor8<T> &S, const symtensor4<T> &D)
{
  symtensor4<T> res;

  res.da[sXXXX] = S.XXXXXXXX * D.da[sZZZZ] + 6.0 * S.XXXXXXYY * D.da[sYYZZ] + 6.0 * S.XXXXXXYY * D.da[sXXZZ] +
                  S.XXXXYYYY * D.da[sYYYY] + 6.0 * S.XXXXYYZZ * D.da[sXXYY] + S.XXXXYYYY * D.da[sXXXX];

  res.da[sXXXY] = 4.0 * S.XXXXXXYY * D.da[sYZZZ] + 4.0 * S.XXXXYYYY * D.da[sYYYZ] + 12.0 * S.XXXXYYZZ * D.da[sXXYZ];

  res.da[sXXXZ] = 4.0 * S.XXXXXXYY * D.da[sXZZZ] + 12.0 * S.XXXXYYZZ * D.da[sXYYZ] + 4.0 * S.XXXXYYYY * D.da[sXXXZ];

  res.da[sXXYY] = S.XXXXXXYY * D.da[sZZZZ] + 6.0 * S.XXXXYYYY * D.da[sYYZZ] + 6.0 * S.XXXXYYZZ * D.da[sXXZZ] +
                  S.XXXXXXYY * D.da[sYYYY] + 6.0 * S.XXXXYYZZ * D.da[sXXYY] + S.XXXXYYZZ * D.da[sXXXX];

  res.da[sXXYZ] = 12.0 * S.XXXXYYZZ * D.da[sXYZZ] + 4.0 * S.XXXXYYZZ * D.da[sXYYY] + 4.0 * S.XXXXYYZZ * D.da[sXXXY];

  res.da[sXXZZ] = S.XXXXXXYY * D.da[sZZZZ] + 6.0 * S.XXXXYYZZ * D.da[sYYZZ] + 6.0 * S.XXXXYYYY * D.da[sXXZZ] +
                  S.XXXXYYZZ * D.da[sYYYY] + 6.0 * S.XXXXYYZZ * D.da[sXXYY] + S.XXXXXXYY * D.da[sXXXX];

  res.da[sXYYY] = 4.0 * S.XXXXYYYY * D.da[sYZZZ] + 4.0 * S.XXXXXXYY * D.da[sYYYZ] + 12.0 * S.XXXXYYZZ * D.da[sXXYZ];

  res.da[sXYYZ] = 4.0 * S.XXXXYYZZ * D.da[sXZZZ] + 12.0 * S.XXXXYYZZ * D.da[sXYYZ] + 4.0 * S.XXXXYYZZ * D.da[sXXXZ];

  res.da[sXYZZ] = 4.0 * S.XXXXYYZZ * D.da[sYZZZ] + 4.0 * S.XXXXYYZZ * D.da[sYYYZ] + 12.0 * S.XXXXYYZZ * D.da[sXXYZ];

  res.da[sXZZZ] = 4.0 * S.XXXXYYYY * D.da[sXZZZ] + 12.0 * S.XXXXYYZZ * D.da[sXYYZ] + 4.0 * S.XXXXXXYY * D.da[sXXXZ];

  res.da[sYYYY] = S.XXXXYYYY * D.da[sZZZZ] + 6.0 * S.XXXXXXYY * D.da[sYYZZ] + 6.0 * S.XXXXYYZZ * D.da[sXXZZ] +
                  S.XXXXXXXX * D.da[sYYYY] + 6.0 * S.XXXXXXYY * D.da[sXXYY] + S.XXXXYYYY * D.da[sXXXX];

  res.da[sYYYZ] = 12.0 * S.XXXXYYZZ * D.da[sXYZZ] + 4.0 * S.XXXXXXYY * D.da[sXYYY] + 4.0 * S.XXXXYYYY * D.da[sXXXY];

  res.da[sYYZZ] = S.XXXXYYZZ * D.da[sZZZZ] + 6.0 * S.XXXXYYZZ * D.da[sYYZZ] + 6.0 * S.XXXXYYZZ * D.da[sXXZZ] +
                  S.XXXXXXYY * D.da[sYYYY] + 6.0 * S.XXXXYYYY * D.da[sXXYY] + S.XXXXXXYY * D.da[sXXXX];

  res.da[sYZZZ] = 12.0 * S.XXXXYYZZ * D.da[sXYZZ] + 4.0 * S.XXXXYYYY * D.da[sXYYY] + 4.0 * S.XXXXXXYY * D.da[sXXXY];

  res.da[sZZZZ] = S.XXXXYYYY * D.da[sZZZZ] + 6.0 * S.XXXXYYZZ * D.da[sYYZZ] + 6.0 * S.XXXXXXYY * D.da[sXXZZ] +
                  S.XXXXYYYY * D.da[sYYYY] + 6.0 * S.XXXXXXYY * D.da[sXXYY] + S.XXXXXXXX * D.da[sXXXX];

  return res;
}

// contract a 8-ewaldtensor with a symmetric 3-tensor to yield a 5-tensor
template <typename T>
inline symtensor5<T> operator*(const ewaldtensor8<T> &S, const symtensor3<T> &D)
{
  symtensor5<T> res;

  res.da[rXXXXX] = S.XXXXXXXX * D.da[dZZZ] + 3.0 * S.XXXXXXYY * D.da[dYYZ] + 3.0 * S.XXXXXXYY * D.da[dXXZ];

  res.da[rXXXXY] = 3.0 * S.XXXXXXYY * D.da[dYZZ] + S.XXXXYYYY * D.da[dYYY] + 3.0 * S.XXXXYYZZ * D.da[dXXY];

  res.da[rXXXXZ] = 3.0 * S.XXXXXXYY * D.da[dXZZ] + 3.0 * S.XXXXYYZZ * D.da[dXYY] + S.XXXXYYYY * D.da[dXXX];

  res.da[rXXXYY] = S.XXXXXXYY * D.da[dZZZ] + 3.0 * S.XXXXYYYY * D.da[dYYZ] + 3.0 * S.XXXXYYZZ * D.da[dXXZ];

  res.da[rXXXYZ] = 6.0 * S.XXXXYYZZ * D.da[dXYZ];

  res.da[rXXXZZ] = S.XXXXXXYY * D.da[dZZZ] + 3.0 * S.XXXXYYZZ * D.da[dYYZ] + 3.0 * S.XXXXYYYY * D.da[dXXZ];

  res.da[rXXYYY] = 3.0 * S.XXXXYYYY * D.da[dYZZ] + S.XXXXXXYY * D.da[dYYY] + 3.0 * S.XXXXYYZZ * D.da[dXXY];

  res.da[rXXYYZ] = 3.0 * S.XXXXYYZZ * D.da[dXZZ] + 3.0 * S.XXXXYYZZ * D.da[dXYY] + S.XXXXYYZZ * D.da[dXXX];

  res.da[rXXYZZ] = 3.0 * S.XXXXYYZZ * D.da[dYZZ] + S.XXXXYYZZ * D.da[dYYY] + 3.0 * S.XXXXYYZZ * D.da[dXXY];

  res.da[rXXZZZ] = 3.0 * S.XXXXYYYY * D.da[dXZZ] + 3.0 * S.XXXXYYZZ * D.da[dXYY] + S.XXXXXXYY * D.da[dXXX];

  res.da[rXYYYY] = S.XXXXYYYY * D.da[dZZZ] + 3.0 * S.XXXXXXYY * D.da[dYYZ] + 3.0 * S.XXXXYYZZ * D.da[dXXZ];

  res.da[rXYYYZ] = 6.0 * S.XXXXYYZZ * D.da[dXYZ];

  res.da[rXYYZZ] = S.XXXXYYZZ * D.da[dZZZ] + 3.0 * S.XXXXYYZZ * D.da[dYYZ] + 3.0 * S.XXXXYYZZ * D.da[dXXZ];

  res.da[rXYZZZ] = 6.0 * S.XXXXYYZZ * D.da[dXYZ];

  res.da[rXZZZZ] = S.XXXXYYYY * D.da[dZZZ] + 3.0 * S.XXXXYYZZ * D.da[dYYZ] + 3.0 * S.XXXXXXYY * D.da[dXXZ];

  res.da[rYYYYY] = 3.0 * S.XXXXXXYY * D.da[dYZZ] + S.XXXXXXXX * D.da[dYYY] + 3.0 * S.XXXXXXYY * D.da[dXXY];

  res.da[rYYYYZ] = 3.0 * S.XXXXYYZZ * D.da[dXZZ] + 3.0 * S.XXXXXXYY * D.da[dXYY] + S.XXXXYYYY * D.da[dXXX];

  res.da[rYYYZZ] = 3.0 * S.XXXXYYZZ * D.da[dYZZ] + S.XXXXXXYY * D.da[dYYY] + 3.0 * S.XXXXYYYY * D.da[dXXY];

  res.da[rYYZZZ] = 3.0 * S.XXXXYYZZ * D.da[dXZZ] + 3.0 * S.XXXXYYYY * D.da[dXYY] + S.XXXXXXYY * D.da[dXXX];

  res.da[rYZZZZ] = 3.0 * S.XXXXYYZZ * D.da[dYZZ] + S.XXXXYYYY * D.da[dYYY] + 3.0 * S.XXXXXXYY * D.da[dXXY];

  res.da[rZZZZZ] = 3.0 * S.XXXXXXYY * D.da[dXZZ] + 3.0 * S.XXXXXXYY * D.da[dXYY] + S.XXXXXXXX * D.da[dXXX];

  return res;
}

// contract a a 10-ewaldtensor with a symmetric 5-tensor to yield a 5-tensor
template <typename T>
inline symtensor5<T> operator*(const ewaldtensor10<T> &S, const symtensor5<T> &D)
{
  symtensor5<T> res;

  res.da[rXXXXX] = S.XXXXXXXXXX * D.da[rZZZZZ] + 10.0 * S.XXXXXXXXYY * D.da[rYYZZZ] + 10.0 * S.XXXXXXXXYY * D.da[rXXZZZ] +
                   5.0 * S.XXXXXXYYYY * D.da[rYYYYZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXXXYYYY * D.da[rXXXXZ];

  res.da[rXXXXY] = 5.0 * S.XXXXXXXXYY * D.da[rYZZZZ] + 10.0 * S.XXXXXXYYYY * D.da[rYYYZZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXXYZZ] +
                   S.XXXXXXYYYY * D.da[rYYYYY] + 10.0 * S.XXXXYYYYZZ * D.da[rXXYYY] + 5.0 * S.XXXXYYYYZZ * D.da[rXXXXY];

  res.da[rXXXXZ] = 5.0 * S.XXXXXXXXYY * D.da[rXZZZZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXXXYYYY * D.da[rXXXZZ] +
                   5.0 * S.XXXXYYYYZZ * D.da[rXYYYY] + 10.0 * S.XXXXYYYYZZ * D.da[rXXXYY] + S.XXXXXXYYYY * D.da[rXXXXX];

  res.da[rXXXYY] = S.XXXXXXXXYY * D.da[rZZZZZ] + 10.0 * S.XXXXXXYYYY * D.da[rYYZZZ] + 10.0 * S.XXXXXXYYZZ * D.da[rXXZZZ] +
                   5.0 * S.XXXXXXYYYY * D.da[rYYYYZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXYYYYZZ * D.da[rXXXXZ];

  res.da[rXXXYZ] = 20.0 * S.XXXXXXYYZZ * D.da[rXYZZZ] + 20.0 * S.XXXXYYYYZZ * D.da[rXYYYZ] + 20.0 * S.XXXXYYYYZZ * D.da[rXXXYZ];

  res.da[rXXXZZ] = S.XXXXXXXXYY * D.da[rZZZZZ] + 10.0 * S.XXXXXXYYZZ * D.da[rYYZZZ] + 10.0 * S.XXXXXXYYYY * D.da[rXXZZZ] +
                   5.0 * S.XXXXYYYYZZ * D.da[rYYYYZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXXXYYYY * D.da[rXXXXZ];

  res.da[rXXYYY] = 5.0 * S.XXXXXXYYYY * D.da[rYZZZZ] + 10.0 * S.XXXXXXYYYY * D.da[rYYYZZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXXYZZ] +
                   S.XXXXXXXXYY * D.da[rYYYYY] + 10.0 * S.XXXXXXYYZZ * D.da[rXXYYY] + 5.0 * S.XXXXYYYYZZ * D.da[rXXXXY];

  res.da[rXXYYZ] = 5.0 * S.XXXXXXYYZZ * D.da[rXZZZZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXYYYYZZ * D.da[rXXXZZ] +
                   5.0 * S.XXXXXXYYZZ * D.da[rXYYYY] + 10.0 * S.XXXXYYYYZZ * D.da[rXXXYY] + S.XXXXXXYYZZ * D.da[rXXXXX];

  res.da[rXXYZZ] = 5.0 * S.XXXXXXYYZZ * D.da[rYZZZZ] + 10.0 * S.XXXXYYYYZZ * D.da[rYYYZZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXXYZZ] +
                   S.XXXXXXYYZZ * D.da[rYYYYY] + 10.0 * S.XXXXYYYYZZ * D.da[rXXYYY] + 5.0 * S.XXXXXXYYZZ * D.da[rXXXXY];

  res.da[rXXZZZ] = 5.0 * S.XXXXXXYYYY * D.da[rXZZZZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXXXYYYY * D.da[rXXXZZ] +
                   5.0 * S.XXXXYYYYZZ * D.da[rXYYYY] + 10.0 * S.XXXXXXYYZZ * D.da[rXXXYY] + S.XXXXXXXXYY * D.da[rXXXXX];

  res.da[rXYYYY] = S.XXXXXXYYYY * D.da[rZZZZZ] + 10.0 * S.XXXXXXYYYY * D.da[rYYZZZ] + 10.0 * S.XXXXYYYYZZ * D.da[rXXZZZ] +
                   5.0 * S.XXXXXXXXYY * D.da[rYYYYZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXYYYYZZ * D.da[rXXXXZ];

  res.da[rXYYYZ] = 20.0 * S.XXXXYYYYZZ * D.da[rXYZZZ] + 20.0 * S.XXXXXXYYZZ * D.da[rXYYYZ] + 20.0 * S.XXXXYYYYZZ * D.da[rXXXYZ];

  res.da[rXYYZZ] = S.XXXXXXYYZZ * D.da[rZZZZZ] + 10.0 * S.XXXXYYYYZZ * D.da[rYYZZZ] + 10.0 * S.XXXXYYYYZZ * D.da[rXXZZZ] +
                   5.0 * S.XXXXXXYYZZ * D.da[rYYYYZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXXXYYZZ * D.da[rXXXXZ];

  res.da[rXYZZZ] = 20.0 * S.XXXXYYYYZZ * D.da[rXYZZZ] + 20.0 * S.XXXXYYYYZZ * D.da[rXYYYZ] + 20.0 * S.XXXXXXYYZZ * D.da[rXXXYZ];

  res.da[rXZZZZ] = S.XXXXXXYYYY * D.da[rZZZZZ] + 10.0 * S.XXXXYYYYZZ * D.da[rYYZZZ] + 10.0 * S.XXXXXXYYYY * D.da[rXXZZZ] +
                   5.0 * S.XXXXYYYYZZ * D.da[rYYYYZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXXYYZ] + 5.0 * S.XXXXXXXXYY * D.da[rXXXXZ];

  res.da[rYYYYY] = 5.0 * S.XXXXXXYYYY * D.da[rYZZZZ] + 10.0 * S.XXXXXXXXYY * D.da[rYYYZZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXXYZZ] +
                   S.XXXXXXXXXX * D.da[rYYYYY] + 10.0 * S.XXXXXXXXYY * D.da[rXXYYY] + 5.0 * S.XXXXXXYYYY * D.da[rXXXXY];

  res.da[rYYYYZ] = 5.0 * S.XXXXYYYYZZ * D.da[rXZZZZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXYYYYZZ * D.da[rXXXZZ] +
                   5.0 * S.XXXXXXXXYY * D.da[rXYYYY] + 10.0 * S.XXXXXXYYYY * D.da[rXXXYY] + S.XXXXXXYYYY * D.da[rXXXXX];

  res.da[rYYYZZ] = 5.0 * S.XXXXYYYYZZ * D.da[rYZZZZ] + 10.0 * S.XXXXXXYYZZ * D.da[rYYYZZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXXYZZ] +
                   S.XXXXXXXXYY * D.da[rYYYYY] + 10.0 * S.XXXXXXYYYY * D.da[rXXYYY] + 5.0 * S.XXXXXXYYYY * D.da[rXXXXY];

  res.da[rYYZZZ] = 5.0 * S.XXXXYYYYZZ * D.da[rXZZZZ] + 30.0 * S.XXXXYYYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXXXYYZZ * D.da[rXXXZZ] +
                   5.0 * S.XXXXXXYYYY * D.da[rXYYYY] + 10.0 * S.XXXXXXYYYY * D.da[rXXXYY] + S.XXXXXXXXYY * D.da[rXXXXX];

  res.da[rYZZZZ] = 5.0 * S.XXXXYYYYZZ * D.da[rYZZZZ] + 10.0 * S.XXXXYYYYZZ * D.da[rYYYZZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXXYZZ] +
                   S.XXXXXXYYYY * D.da[rYYYYY] + 10.0 * S.XXXXXXYYYY * D.da[rXXYYY] + 5.0 * S.XXXXXXXXYY * D.da[rXXXXY];

  res.da[rZZZZZ] = 5.0 * S.XXXXXXYYYY * D.da[rXZZZZ] + 30.0 * S.XXXXXXYYZZ * D.da[rXYYZZ] + 10.0 * S.XXXXXXXXYY * D.da[rXXXZZ] +
                   5.0 * S.XXXXXXYYYY * D.da[rXYYYY] + 10.0 * S.XXXXXXXXYY * D.da[rXXXYY] + S.XXXXXXXXXX * D.da[rXXXXX];

  return res;
}

#endif
