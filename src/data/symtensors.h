/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file symtensors.h
 *
 *  \brief defines symmetric tensors of different rank and basic operations for them
 */

#ifndef SYMTENSORS_H
#define SYMTENSORS_H

#include "gadgetconfig.h"

#include "symtensor_indices.h"

void symtensor_test(void);

template <typename T1, typename T2>
struct which_return;

template <typename T>
struct which_return<T, T>
{
  typedef T type;
};

template <>
struct which_return<float, double>
{
  typedef double type;
};

template <>
struct which_return<double, float>
{
  typedef double type;
};

template <>
struct which_return<int, float>
{
  typedef float type;
};

template <>
struct which_return<float, int>
{
  typedef float type;
};

template <>
struct which_return<int, double>
{
  typedef double type;
};

template <>
struct which_return<double, int>
{
  typedef double type;
};

/* with the above construction, the expression
 *
 *    typename which_return<T1, T2>::type
 *
 * gives us now the more accurate type if mixed precisions are used for T1 and T2
 */

template <typename T>
struct compl_return;

template <typename T>
struct compl_return
{
  typedef double type;
};

template <>
struct compl_return<float>
{
  typedef double type;
};

template <>
struct compl_return<double>
{
  typedef float type;
};

/* with the above construction, the expression
 *
 *    typename compl_return<T>::type
 *
 * gives us the type float of T is double, and the type double if T is type float (i.e. the complementary type)
 * If T is another type, we'll get the type double
 * We'll use this to define some implicit type casts.
 */

// vector
template <typename T>
class vector
{
 public:
  T da[3];

  vector() {} /* constructor */

  inline vector(const T x) /* constructor  */
  {
    da[0] = x;
    da[1] = x;
    da[2] = x;
  }

  inline vector(const T x, const T y, const T z) /* constructor  */
  {
    da[0] = x;
    da[1] = y;
    da[2] = z;
  }

  inline vector(const float *x) /* constructor  */
  {
    da[0] = x[0];
    da[1] = x[1];
    da[2] = x[2];
  }

  inline vector(const double *x) /* constructor  */
  {
    da[0] = x[0];
    da[1] = x[1];
    da[2] = x[2];
  }

  /* implicit conversion operator to type double or float as needed */
  typedef typename compl_return<T>::type float_or_double;
  operator vector<float_or_double>() const { return vector<float_or_double>(da); }

  inline vector &operator+=(const vector<double> &right)
  {
    da[0] += right.da[0];
    da[1] += right.da[1];
    da[2] += right.da[2];

    return *this;
  }

  inline vector &operator+=(const vector<float> &right)
  {
    da[0] += right.da[0];
    da[1] += right.da[1];
    da[2] += right.da[2];

    return *this;
  }

  inline vector &operator-=(const vector<double> &right)
  {
    da[0] -= right.da[0];
    da[1] -= right.da[1];
    da[2] -= right.da[2];

    return *this;
  }

  inline vector &operator-=(const vector<float> &right)
  {
    da[0] -= right.da[0];
    da[1] -= right.da[1];
    da[2] -= right.da[2];

    return *this;
  }

  inline vector &operator*=(const T fac)
  {
    da[0] *= fac;
    da[1] *= fac;
    da[2] *= fac;

    return *this;
  }

  inline T r2(void) { return da[0] * da[0] + da[1] * da[1] + da[2] * da[2]; }

  inline double norm(void) { return sqrt(da[0] * da[0] + da[1] * da[1] + da[2] * da[2]); }

  inline T &operator[](const size_t index) { return da[index]; }
};

// fully symmetric 2-tensor (i.e. symmetic 3x3 matrix)
template <typename T>
class symtensor2
{
 public:
  T da[6];

  symtensor2() {} /* constructor */

  inline symtensor2(const T x) /* constructor  */
  {
    da[0] = x;
    da[1] = x;
    da[2] = x;
    da[3] = x;
    da[4] = x;
    da[5] = x;
  }

  inline symtensor2(const float *x) /* constructor  */
  {
    da[0] = x[0];
    da[1] = x[1];
    da[2] = x[2];
    da[3] = x[3];
    da[4] = x[4];
    da[5] = x[5];
  }

  inline symtensor2(const double *x) /* constructor  */
  {
    da[0] = x[0];
    da[1] = x[1];
    da[2] = x[2];
    da[3] = x[3];
    da[4] = x[4];
    da[5] = x[5];
  }

  /* implicit conversion operator to type float or double as needed */
  typedef typename compl_return<T>::type float_or_double;
  operator symtensor2<float_or_double>() const { return symtensor2<float_or_double>(da); }

  inline symtensor2(const vector<T> &v, const vector<T> &w) /* constructor based on the outer product of two vectors */
  {
    da[qXX] = v.da[0] * w.da[0];
    da[qYY] = v.da[1] * w.da[1];
    da[qZZ] = v.da[2] * w.da[2];
    da[qXY] = v.da[0] * w.da[1];
    da[qXZ] = v.da[0] * w.da[2];
    da[qYZ] = v.da[1] * w.da[2];
  }

  inline symtensor2 &operator+=(const symtensor2 &right)
  {
    da[0] += right.da[0];
    da[1] += right.da[1];
    da[2] += right.da[2];
    da[3] += right.da[3];
    da[4] += right.da[4];
    da[5] += right.da[5];

    return *this;
  }

  inline symtensor2 &operator-=(const symtensor2 &right)
  {
    da[0] -= right.da[0];
    da[1] -= right.da[1];
    da[2] -= right.da[2];
    da[3] -= right.da[3];
    da[4] -= right.da[4];
    da[5] -= right.da[5];

    return *this;
  }

  inline symtensor2 &operator*=(const T fac)
  {
    da[0] *= fac;
    da[1] *= fac;
    da[2] *= fac;
    da[3] *= fac;
    da[4] *= fac;
    da[5] *= fac;

    return *this;
  }

  inline T &operator[](const size_t index) { return da[index]; }

  inline T trace(void) { return da[qXX] + da[qYY] + da[qZZ]; }

  inline double norm(void)
  {
    double sum2 = 0;
    for(int i = 0; i < 6; i++)
      sum2 += da[i] * da[i];

    return sqrt(sum2 / 6);
  }
};

// fully symmetric 3-tensor (3x3x3)
template <typename T>
class symtensor3
{
 public:
  T da[10];

  symtensor3() {} /* constructor */

  inline symtensor3(const T x) /* constructor  */
  {
    da[0] = x;
    da[1] = x;
    da[2] = x;
    da[3] = x;
    da[4] = x;
    da[5] = x;
    da[6] = x;
    da[7] = x;
    da[8] = x;
    da[9] = x;
  }

  inline symtensor3(const float *x) /* constructor  */
  {
    da[0] = x[0];
    da[1] = x[1];
    da[2] = x[2];
    da[3] = x[3];
    da[4] = x[4];
    da[5] = x[5];
    da[6] = x[6];
    da[7] = x[7];
    da[8] = x[8];
    da[9] = x[9];
  }

  inline symtensor3(const double *x) /* constructor  */
  {
    da[0] = x[0];
    da[1] = x[1];
    da[2] = x[2];
    da[3] = x[3];
    da[4] = x[4];
    da[5] = x[5];
    da[6] = x[6];
    da[7] = x[7];
    da[8] = x[8];
    da[9] = x[9];
  }

  /* implicit conversion operator to type float or double as needed */
  typedef typename compl_return<T>::type float_or_double;
  operator symtensor3<float_or_double>() const { return symtensor3<float_or_double>(da); }

  /* constructor based on the outer product  */
  inline symtensor3(const vector<T> &v, const symtensor2<T> &D)
  {
    da[dXXX] = D.da[qXX] * v.da[0];
    da[dXXY] = D.da[qXX] * v.da[1];
    da[dXXZ] = D.da[qXX] * v.da[2];
    da[dXYY] = D.da[qXY] * v.da[1];
    da[dXYZ] = D.da[qXY] * v.da[2];
    da[dXZZ] = D.da[qXZ] * v.da[2];
    da[dYYY] = D.da[qYY] * v.da[1];
    da[dYYZ] = D.da[qYY] * v.da[2];
    da[dYZZ] = D.da[qYZ] * v.da[2];
    da[dZZZ] = D.da[qZZ] * v.da[2];
  }

  inline symtensor3 &operator+=(const symtensor3 &right)
  {
    da[0] += right.da[0];
    da[1] += right.da[1];
    da[2] += right.da[2];
    da[3] += right.da[3];
    da[4] += right.da[4];
    da[5] += right.da[5];
    da[6] += right.da[6];
    da[7] += right.da[7];
    da[8] += right.da[8];
    da[9] += right.da[9];

    return *this;
  }

  inline symtensor3 &operator-=(const symtensor3 &right)
  {
    da[0] -= right.da[0];
    da[1] -= right.da[1];
    da[2] -= right.da[2];
    da[3] -= right.da[3];
    da[4] -= right.da[4];
    da[5] -= right.da[5];
    da[6] -= right.da[6];
    da[7] -= right.da[7];
    da[8] -= right.da[8];
    da[9] -= right.da[9];

    return *this;
  }

  inline symtensor3 &operator*=(const T fac)
  {
    da[0] *= fac;
    da[1] *= fac;
    da[2] *= fac;
    da[3] *= fac;
    da[4] *= fac;
    da[5] *= fac;
    da[6] *= fac;
    da[7] *= fac;
    da[8] *= fac;
    da[9] *= fac;

    return *this;
  }

  inline T &operator[](const size_t index) { return da[index]; }

  inline double norm(void)
  {
    double sum2 = 0;
    for(int i = 0; i < 10; i++)
      sum2 += da[i] * da[i];

    return sqrt(sum2 / 10);
  }
};

// fully symmetric 4-tensor (3x3x3x3)
template <typename T>
class symtensor4
{
 public:
  T da[15];

  symtensor4() {} /* constructor */

  inline symtensor4(const T x) /* constructor  */
  {
    for(int i = 0; i < 15; i++)
      da[i] = x;
  }

  inline symtensor4(const float *x) /* constructor  */
  {
    for(int i = 0; i < 15; i++)
      da[i] = x[i];
  }

  inline symtensor4(const double *x) /* constructor  */
  {
    for(int i = 0; i < 15; i++)
      da[i] = x[i];
  }

  /* implicit conversion operator to type float or double as needed */
  typedef typename compl_return<T>::type float_or_double;
  operator symtensor4<float_or_double>() const { return symtensor4<float_or_double>(da); }

  /* constructor based on an outer product  */
  inline symtensor4(const vector<T> &v, const symtensor3<T> &D)
  {
    da[sXXXX] = D.da[dXXX] * v.da[0];
    da[sXXXY] = D.da[dXXX] * v.da[1];
    da[sXXXZ] = D.da[dXXX] * v.da[2];
    da[sXXYY] = D.da[dXXY] * v.da[1];
    da[sXXYZ] = D.da[dXXY] * v.da[2];
    da[sXXZZ] = D.da[dXXZ] * v.da[2];
    da[sXYYY] = D.da[dXYY] * v.da[1];
    da[sXYYZ] = D.da[dXYY] * v.da[2];
    da[sXYZZ] = D.da[dXYZ] * v.da[2];
    da[sXZZZ] = D.da[dXZZ] * v.da[2];
    da[sYYYY] = D.da[dYYY] * v.da[1];
    da[sYYYZ] = D.da[dYYY] * v.da[2];
    da[sYYZZ] = D.da[dYYZ] * v.da[2];
    da[sYZZZ] = D.da[dYZZ] * v.da[2];
    da[sZZZZ] = D.da[dZZZ] * v.da[2];
  }

  inline symtensor4 &operator+=(const symtensor4 &right)
  {
    for(int i = 0; i < 15; i++)
      da[i] += right.da[i];

    return *this;
  }

  inline symtensor4 &operator-=(const symtensor4 &right)
  {
    for(int i = 0; i < 15; i++)
      da[i] -= right.da[i];

    return *this;
  }

  inline symtensor4 &operator*=(const T fac)
  {
    for(int i = 0; i < 15; i++)
      da[i] *= fac;

    return *this;
  }

  inline T &operator[](const size_t index) { return da[index]; }

  inline double norm(void)
  {
    double sum2 = 0;
    for(int i = 0; i < 15; i++)
      sum2 += da[i] * da[i];

    return sqrt(sum2 / 15);
  }
};

// fully symmetric 5-tensor (3x3x3x3x3)
template <typename T>
class symtensor5
{
 public:
  T da[21];

  symtensor5() {} /* constructor */

  inline symtensor5(const T x) /* constructor  */
  {
    for(int i = 0; i < 21; i++)
      da[i] = x;
  }

  inline symtensor5(const float *x) /* constructor  */
  {
    for(int i = 0; i < 21; i++)
      da[i] = x[i];
  }

  inline symtensor5(const double *x) /* constructor  */
  {
    for(int i = 0; i < 21; i++)
      da[i] = x[i];
  }

  /* implicit conversion operator to type float or double as needed */
  typedef typename compl_return<T>::type float_or_double;
  operator symtensor5<float_or_double>() const { return symtensor5<float_or_double>(da); }

  /* constructor based on an outer product  */
  inline symtensor5(const vector<T> &v, const symtensor4<T> &D)
  {
    da[rXXXXX] = D.da[sXXXX] * v.da[0];
    da[rXXXXY] = D.da[sXXXX] * v.da[1];
    da[rXXXXZ] = D.da[sXXXX] * v.da[2];
    da[rXXXYY] = D.da[sXXXY] * v.da[1];
    da[rXXXYZ] = D.da[sXXXY] * v.da[2];
    da[rXXXZZ] = D.da[sXXXZ] * v.da[2];
    da[rXXYYY] = D.da[sXXYY] * v.da[1];
    da[rXXYYZ] = D.da[sXXYY] * v.da[2];
    da[rXXYZZ] = D.da[sXXYZ] * v.da[2];
    da[rXXZZZ] = D.da[sXXZZ] * v.da[2];
    da[rXYYYY] = D.da[sXYYY] * v.da[1];
    da[rXYYYZ] = D.da[sXYYY] * v.da[2];
    da[rXYYZZ] = D.da[sXYYZ] * v.da[2];
    da[rXYZZZ] = D.da[sXYZZ] * v.da[2];
    da[rXZZZZ] = D.da[sXZZZ] * v.da[2];
    da[rYYYYY] = D.da[sYYYY] * v.da[1];
    da[rYYYYZ] = D.da[sYYYY] * v.da[2];
    da[rYYYZZ] = D.da[sYYYZ] * v.da[2];
    da[rYYZZZ] = D.da[sYYZZ] * v.da[2];
    da[rYZZZZ] = D.da[sYZZZ] * v.da[2];
    da[rZZZZZ] = D.da[sZZZZ] * v.da[2];
  }

  inline symtensor5 &operator+=(const symtensor5 &right)
  {
    for(int i = 0; i < 21; i++)
      da[i] += right.da[i];

    return *this;
  }

  inline symtensor5 &operator-=(const symtensor5 &right)
  {
    for(int i = 0; i < 21; i++)
      da[i] -= right.da[i];

    return *this;
  }

  inline symtensor5 &operator*=(const T fac)
  {
    for(int i = 0; i < 21; i++)
      da[i] *= fac;

    return *this;
  }

  inline T &operator[](const size_t index) { return da[index]; }

  inline double norm(void)
  {
    double sum2 = 0;
    for(int i = 0; i < 21; i++)
      sum2 += da[i] * da[i];

    return sqrt(sum2 / 21);
  }
};

// fully symmetric 6-tensor (3x3x3x3x3x3)
template <typename T>
class symtensor6
{
 public:
  T da[28];

  symtensor6() {} /* constructor */

  inline symtensor6(const T x) /* constructor  */
  {
    for(int i = 0; i < 28; i++)
      da[i] = x;
  }

  inline symtensor6(const float *x) /* constructor  */
  {
    for(int i = 0; i < 28; i++)
      da[i] = x[i];
  }

  inline symtensor6(const double *x) /* constructor  */
  {
    for(int i = 0; i < 28; i++)
      da[i] = x[i];
  }

  /* implicit conversion operator to type float or double as needed */
  typedef typename compl_return<T>::type float_or_double;
  operator symtensor6<float_or_double>() const { return symtensor6<float_or_double>(da); }

  /* constructor based on an outer product  */
  inline symtensor6(const vector<T> &v, const symtensor5<T> &D)
  {
    da[pXXXXXX] = D.da[rXXXXX] * v.da[0];
    da[pXXXXXY] = D.da[rXXXXX] * v.da[1];
    da[pXXXXXZ] = D.da[rXXXXX] * v.da[2];
    da[pXXXXYY] = D.da[rXXXXY] * v.da[1];
    da[pXXXXYZ] = D.da[rXXXXY] * v.da[2];
    da[pXXXXZZ] = D.da[rXXXXZ] * v.da[2];
    da[pXXXYYY] = D.da[rXXXYY] * v.da[1];
    da[pXXXYYZ] = D.da[rXXXYY] * v.da[2];
    da[pXXXYZZ] = D.da[rXXXYZ] * v.da[2];
    da[pXXXZZZ] = D.da[rXXXZZ] * v.da[2];
    da[pXXYYYY] = D.da[rXXYYY] * v.da[1];
    da[pXXYYYZ] = D.da[rXXYYY] * v.da[2];
    da[pXXYYZZ] = D.da[rXXYYZ] * v.da[2];
    da[pXXYZZZ] = D.da[rXXYZZ] * v.da[2];
    da[pXXZZZZ] = D.da[rXXZZZ] * v.da[2];
    da[pXYYYYY] = D.da[rXYYYY] * v.da[1];
    da[pXYYYYZ] = D.da[rXYYYY] * v.da[2];
    da[pXYYYZZ] = D.da[rXYYYZ] * v.da[2];
    da[pXYYZZZ] = D.da[rXYYZZ] * v.da[2];
    da[pXYZZZZ] = D.da[rXYZZZ] * v.da[2];
    da[pXZZZZZ] = D.da[rXZZZZ] * v.da[2];
    da[pYYYYYY] = D.da[rYYYYY] * v.da[1];
    da[pYYYYYZ] = D.da[rYYYYY] * v.da[2];
    da[pYYYYZZ] = D.da[rYYYYZ] * v.da[2];
    da[pYYYZZZ] = D.da[rYYYZZ] * v.da[2];
    da[pYYZZZZ] = D.da[rYYZZZ] * v.da[2];
    da[pYZZZZZ] = D.da[rYZZZZ] * v.da[2];
    da[pZZZZZZ] = D.da[rZZZZZ] * v.da[2];
  }

  inline symtensor6 &operator+=(const symtensor6 &right)
  {
    for(int i = 0; i < 28; i++)
      da[i] += right.da[i];

    return *this;
  }

  inline symtensor6 &operator-=(const symtensor6 &right)
  {
    for(int i = 0; i < 28; i++)
      da[i] -= right.da[i];

    return *this;
  }

  inline symtensor6 &operator*=(const T fac)
  {
    for(int i = 0; i < 28; i++)
      da[i] *= fac;

    return *this;
  }

  inline T &operator[](const size_t index) { return da[index]; }

  inline double norm(void)
  {
    double sum2 = 0;
    for(int i = 0; i < 28; i++)
      sum2 += da[i] * da[i];

    return sqrt(sum2 / 28);
  }
};

// fully symmetric 7-tensor (3x3x3x3x3x3x3)
template <typename T>
class symtensor7
{
 public:
  T da[36];

  symtensor7() {} /* constructor */

  inline symtensor7(const T x) /* constructor  */
  {
    for(int i = 0; i < 36; i++)
      da[i] = x;
  }

  inline symtensor7(const float *x) /* constructor  */
  {
    for(int i = 0; i < 36; i++)
      da[i] = x[i];
  }

  inline symtensor7(const double *x) /* constructor  */
  {
    for(int i = 0; i < 36; i++)
      da[i] = x[i];
  }

  /* implicit conversion operator to type float or double as needed */
  typedef typename compl_return<T>::type float_or_double;
  operator symtensor7<float_or_double>() const { return symtensor7<float_or_double>(da); }

  inline symtensor7 &operator+=(const symtensor7 &right)
  {
    for(int i = 0; i < 36; i++)
      da[i] += right.da[i];

    return *this;
  }

  inline symtensor7 &operator-=(const symtensor7 &right)
  {
    for(int i = 0; i < 36; i++)
      da[i] -= right.da[i];

    return *this;
  }

  inline symtensor7 &operator*=(const T fac)
  {
    for(int i = 0; i < 36; i++)
      da[i] *= fac;

    return *this;
  }

  inline T &operator[](const size_t index) { return da[index]; }

  inline double norm(void)
  {
    double sum2 = 0;
    for(int i = 0; i < 36; i++)
      sum2 += da[i] * da[i];

    return sqrt(sum2 / 36);
  }
};

//-------------  let's defined additions of these tensors

// add two vectors
template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> operator+(const vector<T1> &left, const vector<T2> &right)
{
  vector<typename which_return<T1, T2>::type> res;

  res.da[0] = left.da[0] + right.da[0];
  res.da[1] = left.da[1] + right.da[1];
  res.da[2] = left.da[2] + right.da[2];

  return res;
}

// add two 2-tensors
template <typename T1, typename T2>
inline symtensor2<typename which_return<T1, T2>::type> operator+(const symtensor2<T1> &left, const symtensor2<T2> &right)
{
  symtensor2<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 6; i++)
    res.da[i] = left.da[i] + right.da[i];

  return res;
}

// add two 3-tensors
template <typename T1, typename T2>
inline symtensor3<typename which_return<T1, T2>::type> operator+(const symtensor3<T1> &left, const symtensor3<T2> &right)
{
  symtensor3<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 10; i++)
    res.da[i] = left.da[i] + right.da[i];

  return res;
}

// add two 4-tensors
template <typename T1, typename T2>
inline symtensor4<typename which_return<T1, T2>::type> operator+(const symtensor4<T1> &left, const symtensor4<T2> &right)
{
  symtensor4<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 15; i++)
    res.da[i] = left.da[i] + right.da[i];

  return res;
}

// add two 5-tensors
template <typename T1, typename T2>
inline symtensor5<typename which_return<T1, T2>::type> operator+(const symtensor5<T1> &left, const symtensor5<T2> &right)
{
  symtensor5<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 21; i++)
    res.da[i] = left.da[i] + right.da[i];

  return res;
}

// add two 6-tensors
template <typename T1, typename T2>
inline symtensor6<typename which_return<T1, T2>::type> operator+(const symtensor6<T1> &left, const symtensor6<T2> &right)
{
  symtensor6<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 28; i++)
    res.da[i] = left.da[i] + right.da[i];

  return res;
}

// add two 7-tensors
template <typename T1, typename T2>
inline symtensor7<typename which_return<T1, T2>::type> operator+(const symtensor7<T1> &left, const symtensor7<T2> &right)
{
  symtensor7<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 36; i++)
    res.da[i] = left.da[i] + right.da[i];

  return res;
}

//-------------  let's defined subtractions of these tensors

// subtract two vectors
template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> operator-(const vector<T1> &left, const vector<T2> &right)
{
  vector<typename which_return<T1, T2>::type> res;

  res.da[0] = left.da[0] - right.da[0];
  res.da[1] = left.da[1] - right.da[1];
  res.da[2] = left.da[2] - right.da[2];

  return res;
}

// subtract two 2-tensors
template <typename T1, typename T2>
inline symtensor2<typename which_return<T1, T2>::type> operator-(const symtensor2<T1> &left, const symtensor2<T2> &right)
{
  symtensor2<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 6; i++)
    res.da[i] = left.da[i] - right.da[i];

  return res;
}

// subtract two 3-tensors
template <typename T1, typename T2>
inline symtensor3<typename which_return<T1, T2>::type> operator-(const symtensor3<T1> &left, const symtensor3<T2> &right)
{
  symtensor3<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 10; i++)
    res.da[i] = left.da[i] - right.da[i];

  return res;
}

// subtract two 4-tensors
template <typename T1, typename T2>
inline symtensor4<typename which_return<T1, T2>::type> operator-(const symtensor4<T1> &left, const symtensor4<T2> &right)
{
  symtensor4<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 15; i++)
    res.da[i] = left.da[i] - right.da[i];

  return res;
}

// subtract two 5-tensors
template <typename T1, typename T2>
inline symtensor5<typename which_return<T1, T2>::type> operator-(const symtensor5<T1> &left, const symtensor5<T2> &right)
{
  symtensor5<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 21; i++)
    res.da[i] = left.da[i] - right.da[i];

  return res;
}

// subtract two 6-tensors
template <typename T1, typename T2>
inline symtensor6<typename which_return<T1, T2>::type> operator-(const symtensor6<T1> &left, const symtensor6<T2> &right)
{
  symtensor6<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 28; i++)
    res.da[i] = left.da[i] - right.da[i];

  return res;
}

// subtract two 7-tensors
template <typename T1, typename T2>
inline symtensor7<typename which_return<T1, T2>::type> operator-(const symtensor7<T1> &left, const symtensor7<T2> &right)
{
  symtensor7<typename which_return<T1, T2>::type> res;

  for(int i = 0; i < 36; i++)
    res.da[i] = left.da[i] - right.da[i];

  return res;
}

//-------------  let's define multiplications with a scalar

// scalar times vector multiply
template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> operator*(const T1 fac, const vector<T2> &v)
{
  vector<typename which_return<T1, T2>::type> res;
  for(int i = 0; i < 3; i++)
    res.da[i] = fac * v.da[i];

  return res;
}

// scalar times 2-tensor multiply
template <typename T1, typename T2>
inline symtensor2<typename which_return<T1, T2>::type> operator*(const T1 fac, const symtensor2<T2> &v)
{
  symtensor2<typename which_return<T1, T2>::type> res;
  for(int i = 0; i < 6; i++)
    res.da[i] = fac * v.da[i];

  return res;
}

// scalar times 3-tensor multiply
template <typename T1, typename T2>
inline symtensor3<typename which_return<T1, T2>::type> operator*(const T1 fac, const symtensor3<T2> &v)
{
  symtensor3<typename which_return<T1, T2>::type> res;
  for(int i = 0; i < 10; i++)
    res.da[i] = fac * v.da[i];

  return res;
}

// scalar times 4-tensor multiply
template <typename T1, typename T2>
inline symtensor4<typename which_return<T1, T2>::type> operator*(const T1 fac, const symtensor4<T2> &v)
{
  symtensor4<typename which_return<T1, T2>::type> res;
  for(int i = 0; i < 15; i++)
    res.da[i] = fac * v.da[i];

  return res;
}

// scalar times 5-tensor multiply
template <typename T1, typename T2>
inline symtensor5<typename which_return<T1, T2>::type> operator*(const T1 fac, const symtensor5<T2> &v)
{
  symtensor5<typename which_return<T1, T2>::type> res;
  for(int i = 0; i < 21; i++)
    res.da[i] = fac * v.da[i];

  return res;
}

// scalar times 6-tensor multiply
template <typename T1, typename T2>
inline symtensor6<typename which_return<T1, T2>::type> operator*(const T1 fac, const symtensor6<T2> &v)
{
  symtensor6<typename which_return<T1, T2>::type> res;
  for(int i = 0; i < 28; i++)
    res.da[i] = fac * v.da[i];

  return res;
}

// scalar 7-tensor multiply
template <typename T1, typename T2>
inline symtensor7<typename which_return<T1, T2>::type> operator*(const T1 fac, const symtensor7<T2> &v)
{
  symtensor7<typename which_return<T1, T2>::type> res;
  for(int i = 0; i < 36; i++)
    res.da[i] = fac * v.da[i];

  return res;
}

//-------------  let's defined contractions of these tensors

// 2-tensor contraction with a vector (ordinary matrix-vector multiplication)
template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> operator*(const symtensor2<T1> &D, const vector<T2> &v)
{
  vector<typename which_return<T1, T2>::type> res;

  res.da[0] = D.da[qXX] * v.da[0] + D.da[qXY] * v.da[1] + D.da[qXZ] * v.da[2];
  res.da[1] = D.da[qYX] * v.da[0] + D.da[qYY] * v.da[1] + D.da[qYZ] * v.da[2];
  res.da[2] = D.da[qZX] * v.da[0] + D.da[qZY] * v.da[1] + D.da[qZZ] * v.da[2];

  return res;
}

// scalar product of two vectors
template <typename T1, typename T2>
inline typename which_return<T1, T2>::type operator*(const vector<T1> &v, const vector<T2> &w)
{
  return v.da[0] * w.da[0] + v.da[1] * w.da[1] + v.da[2] * w.da[2];
}

// contract two 2-tensors to a scalar
template <typename T1, typename T2>
inline typename which_return<T1, T2>::type operator*(const symtensor2<T1> &D, const symtensor2<T2> &Q)
{
  return (D.da[qXX] * Q.da[qXX] + D.da[qYY] * Q.da[qYY] + D.da[qZZ] * Q.da[qZZ]) +
         2 * (D.da[qXY] * Q.da[qXY] + D.da[qXZ] * Q.da[qXZ] + D.da[qYZ] * Q.da[qYZ]);
}

// contract two 3-tensors to yield a scalar
template <typename T1, typename T2>
inline typename which_return<T1, T2>::type operator*(const symtensor3<T1> &D, const symtensor3<T2> &Q)
{
  return D.da[dXXX] * Q.da[dXXX] + D.da[dYYY] * Q.da[dYYY] + D.da[dZZZ] * Q.da[dZZZ] +
         3 * (D.da[dYZZ] * Q.da[dYZZ] + D.da[dYYZ] * Q.da[dYYZ] + D.da[dXZZ] * Q.da[dXZZ] + D.da[dXYY] * Q.da[dXYY] +
              D.da[dXXY] * Q.da[dXXY] + D.da[dXXZ] * Q.da[dXXZ]) +
         6 * D.da[dXYZ] * Q.da[dXYZ];
  ;
}

// contract a 4-tensor with a 4-tensor to yield a scalar
template <typename T1, typename T2>
inline typename which_return<T1, T2>::type operator*(const symtensor4<T1> &D, const symtensor4<T2> &Q)  // checked
{
  return D.da[sZZZZ] * Q.da[sZZZZ] + D.da[sYYYY] * Q.da[sYYYY] + D.da[sXXXX] * Q.da[sXXXX] +
         6 * (D.da[sYYZZ] * Q.da[sYYZZ] + D.da[sXXZZ] * Q.da[sXXZZ] + D.da[sXXYY] * Q.da[sXXYY]) +
         4 * (D.da[sYZZZ] * Q.da[sYZZZ] + D.da[sYYYZ] * Q.da[sYYYZ] + D.da[sXZZZ] * Q.da[sXZZZ] + D.da[sXYYY] * Q.da[sXYYY] +
              D.da[sXXXZ] * Q.da[sXXXZ] + D.da[sXXXY] * Q.da[sXXXY]) +
         12 * (D.da[sXYZZ] * Q.da[sXYZZ] + D.da[sXYYZ] * Q.da[sXYYZ] + D.da[sXXYZ] * Q.da[sXXYZ]);
}

// contract a 3-tensor with a vector to yield a 2-tensor
template <typename T1, typename T2>
inline symtensor2<typename which_return<T1, T2>::type> operator*(const symtensor3<T1> &D, const vector<T2> &v)
{
  symtensor2<typename which_return<T1, T2>::type> res;

  res.da[qXX] = D.da[dXXX] * v.da[0] + D.da[dXXY] * v.da[1] + D.da[dXXZ] * v.da[2];
  res.da[qYY] = D.da[dYYX] * v.da[0] + D.da[dYYY] * v.da[1] + D.da[dYYZ] * v.da[2];
  res.da[qZZ] = D.da[dZZX] * v.da[0] + D.da[dZZY] * v.da[1] + D.da[dZZZ] * v.da[2];
  res.da[qXY] = D.da[dXYX] * v.da[0] + D.da[dXYY] * v.da[1] + D.da[dXYZ] * v.da[2];
  res.da[qXZ] = D.da[dXZX] * v.da[0] + D.da[dXZY] * v.da[1] + D.da[dXZZ] * v.da[2];
  res.da[qYZ] = D.da[dYZX] * v.da[0] + D.da[dYZY] * v.da[1] + D.da[dYZZ] * v.da[2];

  return res;
}

// contract a 3-tensor with a 2-tensor to yield a vector
template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> operator*(const symtensor3<T1> &D, const symtensor2<T2> &Q)
{
  vector<typename which_return<T1, T2>::type> res;

  res.da[0] = (D.da[dXXX] * Q.da[qXX] + D.da[dXYY] * Q.da[qYY] + D.da[dXZZ] * Q.da[qZZ]) +
              2 * (D.da[dXXY] * Q.da[qXY] + D.da[dXXZ] * Q.da[qXZ] + D.da[dXYZ] * Q.da[qYZ]);
  res.da[1] = (D.da[dYXX] * Q.da[qXX] + D.da[dYYY] * Q.da[qYY] + D.da[dYZZ] * Q.da[qZZ]) +
              2 * (D.da[dYXY] * Q.da[qXY] + D.da[dYXZ] * Q.da[qXZ] + D.da[dYYZ] * Q.da[qYZ]);
  res.da[2] = (D.da[dZXX] * Q.da[qXX] + D.da[dZYY] * Q.da[qYY] + D.da[dZZZ] * Q.da[qZZ]) +
              2 * (D.da[dZXY] * Q.da[qXY] + D.da[dZXZ] * Q.da[qXZ] + D.da[dZYZ] * Q.da[qYZ]);

  return res;
}

// contract a 4-tensor with a 3-tensor to yield a vector
template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> operator*(const symtensor4<T1> &D, const symtensor3<T2> &Q)
{
  vector<typename which_return<T1, T2>::type> res;

  res[0] = D.da[sXXXX] * Q.da[dXXX] + D.da[sXYYY] * Q.da[dYYY] + D.da[sXZZZ] * Q.da[dZZZ] +
           3 * (D.da[sXYZZ] * Q.da[dYZZ] + D.da[sXYYZ] * Q.da[dYYZ] + D.da[sXXZZ] * Q.da[dXZZ] + D.da[sXXYY] * Q.da[dXYY] +
                D.da[sXXXY] * Q.da[dXXY] + D.da[sXXXZ] * Q.da[dXXZ]) +
           6 * D.da[sXXYZ] * Q.da[dXYZ];

  res[1] = D.da[sYXXX] * Q.da[dXXX] + D.da[sYYYY] * Q.da[dYYY] + D.da[sYZZZ] * Q.da[dZZZ] +
           3 * (D.da[sYYZZ] * Q.da[dYZZ] + D.da[sYYYZ] * Q.da[dYYZ] + D.da[sYXZZ] * Q.da[dXZZ] + D.da[sYXYY] * Q.da[dXYY] +
                D.da[sYXXY] * Q.da[dXXY] + D.da[sYXXZ] * Q.da[dXXZ]) +
           6 * D.da[sYXYZ] * Q.da[dXYZ];

  res[2] = D.da[sZXXX] * Q.da[dXXX] + D.da[sZYYY] * Q.da[dYYY] + D.da[sZZZZ] * Q.da[dZZZ] +
           3 * (D.da[sZYZZ] * Q.da[dYZZ] + D.da[sZYYZ] * Q.da[dYYZ] + D.da[sZXZZ] * Q.da[dXZZ] + D.da[sZXYY] * Q.da[dXYY] +
                D.da[sZXXY] * Q.da[dXXY] + D.da[sZXXZ] * Q.da[dXXZ]) +
           6 * D.da[sZXYZ] * Q.da[dXYZ];

  return res;
}

// contract a 4-tensor with a vector to yield a 3-tensor
template <typename T1, typename T2>
inline symtensor3<typename which_return<T1, T2>::type> operator*(const symtensor4<T1> &D, const vector<T2> &v)  // checked
{
  symtensor3<typename which_return<T1, T2>::type> res;

  res.da[dXXX] = D.da[sXXXX] * v.da[0] + D.da[sXXXY] * v.da[1] + D.da[sXXXZ] * v.da[2];
  res.da[dYYY] = D.da[sYYYX] * v.da[0] + D.da[sYYYY] * v.da[1] + D.da[sYYYZ] * v.da[2];
  res.da[dZZZ] = D.da[sZZZX] * v.da[0] + D.da[sZZZY] * v.da[1] + D.da[sZZZZ] * v.da[2];
  res.da[dXYY] = D.da[sXYYX] * v.da[0] + D.da[sXYYY] * v.da[1] + D.da[sXYYZ] * v.da[2];
  res.da[dXZZ] = D.da[sXZZX] * v.da[0] + D.da[sXZZY] * v.da[1] + D.da[sXZZZ] * v.da[2];
  res.da[dYXX] = D.da[sYXXX] * v.da[0] + D.da[sYXXY] * v.da[1] + D.da[sYXXZ] * v.da[2];
  res.da[dYZZ] = D.da[sYZZX] * v.da[0] + D.da[sYZZY] * v.da[1] + D.da[sYZZZ] * v.da[2];
  res.da[dZXX] = D.da[sZXXX] * v.da[0] + D.da[sZXXY] * v.da[1] + D.da[sZXXZ] * v.da[2];
  res.da[dZYY] = D.da[sZYYX] * v.da[0] + D.da[sZYYY] * v.da[1] + D.da[sZYYZ] * v.da[2];
  res.da[dXYZ] = D.da[sXYZX] * v.da[0] + D.da[sXYZY] * v.da[1] + D.da[sXYZZ] * v.da[2];

  return res;
}

// contract a 5-tensor with a 4-tensor to yield a vector
template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> operator*(const symtensor5<T1> &D, const symtensor4<T2> &Q)
{
  vector<typename which_return<T1, T2>::type> res;

  res.da[0] = D.da[rXZZZZ] * Q.da[sZZZZ] + D.da[rXYYYY] * Q.da[sYYYY] + D.da[rXXXXX] * Q.da[sXXXX] +
              6 * (D.da[rXYYZZ] * Q.da[sYYZZ] + D.da[rXXXZZ] * Q.da[sXXZZ] + D.da[rXXXYY] * Q.da[sXXYY]) +
              4 * (D.da[rXYZZZ] * Q.da[sYZZZ] + D.da[rXYYYZ] * Q.da[sYYYZ] + D.da[rXXZZZ] * Q.da[sXZZZ] + D.da[rXXYYY] * Q.da[sXYYY] +
                   D.da[rXXXXZ] * Q.da[sXXXZ] + D.da[rXXXXY] * Q.da[sXXXY]) +
              12 * (D.da[rXXYZZ] * Q.da[sXYZZ] + D.da[rXXYYZ] * Q.da[sXYYZ] + D.da[rXXXYZ] * Q.da[sXXYZ]);

  res.da[1] = D.da[rYZZZZ] * Q.da[sZZZZ] + D.da[rYYYYY] * Q.da[sYYYY] + D.da[rYXXXX] * Q.da[sXXXX] +
              6 * (D.da[rYYYZZ] * Q.da[sYYZZ] + D.da[rYXXZZ] * Q.da[sXXZZ] + D.da[rYXXYY] * Q.da[sXXYY]) +
              4 * (D.da[rYYZZZ] * Q.da[sYZZZ] + D.da[rYYYYZ] * Q.da[sYYYZ] + D.da[rYXZZZ] * Q.da[sXZZZ] + D.da[rYXYYY] * Q.da[sXYYY] +
                   D.da[rYXXXZ] * Q.da[sXXXZ] + D.da[rYXXXY] * Q.da[sXXXY]) +
              12 * (D.da[rYXYZZ] * Q.da[sXYZZ] + D.da[rYXYYZ] * Q.da[sXYYZ] + D.da[rYXXYZ] * Q.da[sXXYZ]);

  res.da[2] = D.da[rZZZZZ] * Q.da[sZZZZ] + D.da[rZYYYY] * Q.da[sYYYY] + D.da[rZXXXX] * Q.da[sXXXX] +
              6 * (D.da[rZYYZZ] * Q.da[sYYZZ] + D.da[rZXXZZ] * Q.da[sXXZZ] + D.da[rZXXYY] * Q.da[sXXYY]) +
              4 * (D.da[rZYZZZ] * Q.da[sYZZZ] + D.da[rZYYYZ] * Q.da[sYYYZ] + D.da[rZXZZZ] * Q.da[sXZZZ] + D.da[rZXYYY] * Q.da[sXYYY] +
                   D.da[rZXXXZ] * Q.da[sXXXZ] + D.da[rZXXXY] * Q.da[sXXXY]) +
              12 * (D.da[rZXYZZ] * Q.da[sXYZZ] + D.da[rZXYYZ] * Q.da[sXYYZ] + D.da[rZXXYZ] * Q.da[sXXYZ]);

  return res;
}

// contract a 4-tensor with a 2-tensor to yield a 2-tensor
template <typename T1, typename T2>
inline symtensor2<typename which_return<T1, T2>::type> operator*(const symtensor4<T1> &D, const symtensor2<T2> &Q)  // checked
{
  symtensor2<typename which_return<T1, T2>::type> res;

  res.da[qXX] = D.da[sXXXX] * Q.da[qXX] + D.da[sXXYY] * Q.da[qYY] + D.da[sXXZZ] * Q.da[qZZ] +
                2 * (D.da[sXXXY] * Q.da[qXY] + D.da[sXXXZ] * Q.da[qXZ] + D.da[sXXYZ] * Q.da[qYZ]);
  res.da[qYY] = D.da[sYYXX] * Q.da[qXX] + D.da[sYYYY] * Q.da[qYY] + D.da[sYYZZ] * Q.da[qZZ] +
                2 * (D.da[sYYXY] * Q.da[qXY] + D.da[sYYXZ] * Q.da[qXZ] + D.da[sYYYZ] * Q.da[qYZ]);
  res.da[qZZ] = D.da[sZZXX] * Q.da[qXX] + D.da[sZZYY] * Q.da[qYY] + D.da[sZZZZ] * Q.da[qZZ] +
                2 * (D.da[sZZXY] * Q.da[qXY] + D.da[sZZXZ] * Q.da[qXZ] + D.da[sZZYZ] * Q.da[qYZ]);
  res.da[qXY] = D.da[sXYXX] * Q.da[qXX] + D.da[sXYYY] * Q.da[qYY] + D.da[sXYZZ] * Q.da[qZZ] +
                2 * (D.da[sXYXY] * Q.da[qXY] + D.da[sXYXZ] * Q.da[qXZ] + D.da[sXYYZ] * Q.da[qYZ]);
  res.da[qXZ] = D.da[sXZXX] * Q.da[qXX] + D.da[sXZYY] * Q.da[qYY] + D.da[sXZZZ] * Q.da[qZZ] +
                2 * (D.da[sXZXY] * Q.da[qXY] + D.da[sXZXZ] * Q.da[qXZ] + D.da[sXZYZ] * Q.da[qYZ]);
  res.da[qYZ] = D.da[sYZXX] * Q.da[qXX] + D.da[sYZYY] * Q.da[qYY] + D.da[sYZZZ] * Q.da[qZZ] +
                2 * (D.da[sYZXY] * Q.da[qXY] + D.da[sYZXZ] * Q.da[qXZ] + D.da[sYZYZ] * Q.da[qYZ]);

  return res;
}

// contract a 5-tensor with a 3-tensor to yield a 2-tensor
template <typename T1, typename T2>
inline symtensor2<typename which_return<T1, T2>::type> operator*(const symtensor5<T1> &D, const symtensor3<T2> &Q)
{
  symtensor2<typename which_return<T1, T2>::type> res;

  res.da[qXX] = D.da[rXXXXX] * Q.da[dXXX] + D.da[rXXYYY] * Q.da[dYYY] + D.da[rXXZZZ] * Q.da[dZZZ] +
                3 * (D.da[rXXYZZ] * Q.da[dYZZ] + D.da[rXXYYZ] * Q.da[dYYZ] + D.da[rXXXZZ] * Q.da[dXZZ] + D.da[rXXXYY] * Q.da[dXYY] +
                     D.da[rXXXXY] * Q.da[dXXY] + D.da[rXXXXZ] * Q.da[dXXZ]) +
                6 * D.da[rXXXYZ] * Q.da[dXYZ];

  res.da[qYY] = D.da[rYYXXX] * Q.da[dXXX] + D.da[rYYYYY] * Q.da[dYYY] + D.da[rYYZZZ] * Q.da[dZZZ] +
                3 * (D.da[rYYYZZ] * Q.da[dYZZ] + D.da[rYYYYZ] * Q.da[dYYZ] + D.da[rYYXZZ] * Q.da[dXZZ] + D.da[rYYXYY] * Q.da[dXYY] +
                     D.da[rYYXXY] * Q.da[dXXY] + D.da[rYYXXZ] * Q.da[dXXZ]) +
                6 * D.da[rYYXYZ] * Q.da[dXYZ];

  res.da[qZZ] = D.da[rZZXXX] * Q.da[dXXX] + D.da[rZZYYY] * Q.da[dYYY] + D.da[rZZZZZ] * Q.da[dZZZ] +
                3 * (D.da[rZZYZZ] * Q.da[dYZZ] + D.da[rZZYYZ] * Q.da[dYYZ] + D.da[rZZXZZ] * Q.da[dXZZ] + D.da[rZZXYY] * Q.da[dXYY] +
                     D.da[rZZXXY] * Q.da[dXXY] + D.da[rZZXXZ] * Q.da[dXXZ]) +
                6 * D.da[rZZXYZ] * Q.da[dXYZ];

  res.da[qXY] = D.da[rXYXXX] * Q.da[dXXX] + D.da[rXYYYY] * Q.da[dYYY] + D.da[rXYZZZ] * Q.da[dZZZ] +
                3 * (D.da[rXYYZZ] * Q.da[dYZZ] + D.da[rXYYYZ] * Q.da[dYYZ] + D.da[rXYXZZ] * Q.da[dXZZ] + D.da[rXYXYY] * Q.da[dXYY] +
                     D.da[rXYXXY] * Q.da[dXXY] + D.da[rXYXXZ] * Q.da[dXXZ]) +
                6 * D.da[rXYXYZ] * Q.da[dXYZ];

  res.da[qXZ] = D.da[rXZXXX] * Q.da[dXXX] + D.da[rXZYYY] * Q.da[dYYY] + D.da[rXZZZZ] * Q.da[dZZZ] +
                3 * (D.da[rXZYZZ] * Q.da[dYZZ] + D.da[rXZYYZ] * Q.da[dYYZ] + D.da[rXZXZZ] * Q.da[dXZZ] + D.da[rXZXYY] * Q.da[dXYY] +
                     D.da[rXZXXY] * Q.da[dXXY] + D.da[rXZXXZ] * Q.da[dXXZ]) +
                6 * D.da[rXZXYZ] * Q.da[dXYZ];

  res.da[qYZ] = D.da[rYZXXX] * Q.da[dXXX] + D.da[rYZYYY] * Q.da[dYYY] + D.da[rYZZZZ] * Q.da[dZZZ] +
                3 * (D.da[rYZYZZ] * Q.da[dYZZ] + D.da[rYZYYZ] * Q.da[dYYZ] + D.da[rYZXZZ] * Q.da[dXZZ] + D.da[rYZXYY] * Q.da[dXYY] +
                     D.da[rYZXXY] * Q.da[dXXY] + D.da[rYZXXZ] * Q.da[dXXZ]) +
                6 * D.da[rYZXYZ] * Q.da[dXYZ];

  return res;
}

// contract a 5-tensor with a 2-tensor to yield a 3-tensor
template <typename T1, typename T2>
inline symtensor3<typename which_return<T1, T2>::type> operator*(const symtensor5<T1> &D, const symtensor2<T2> &Q)
{
  symtensor3<typename which_return<T1, T2>::type> res;

  res.da[dXXX] = (D.da[rXXXXX] * Q.da[qXX] + D.da[rXXXYY] * Q.da[qYY] + D.da[rXXXZZ] * Q.da[qZZ]) +
                 2 * (D.da[rXXXXY] * Q.da[qXY] + D.da[rXXXXZ] * Q.da[qXZ] + D.da[rXXXYZ] * Q.da[qYZ]);

  res.da[dYYY] = (D.da[rYYYXX] * Q.da[qXX] + D.da[rYYYYY] * Q.da[qYY] + D.da[rYYYZZ] * Q.da[qZZ]) +
                 2 * (D.da[rYYYXY] * Q.da[qXY] + D.da[rYYYXZ] * Q.da[qXZ] + D.da[rYYYYZ] * Q.da[qYZ]);

  res.da[dZZZ] = (D.da[rZZZXX] * Q.da[qXX] + D.da[rZZZYY] * Q.da[qYY] + D.da[rZZZZZ] * Q.da[qZZ]) +
                 2 * (D.da[rZZZXY] * Q.da[qXY] + D.da[rZZZXZ] * Q.da[qXZ] + D.da[rZZZYZ] * Q.da[qYZ]);

  res.da[dXYY] = (D.da[rXYYXX] * Q.da[qXX] + D.da[rXYYYY] * Q.da[qYY] + D.da[rXYYZZ] * Q.da[qZZ]) +
                 2 * (D.da[rXYYXY] * Q.da[qXY] + D.da[rXYYXZ] * Q.da[qXZ] + D.da[rXYYYZ] * Q.da[qYZ]);

  res.da[dXZZ] = (D.da[rXZZXX] * Q.da[qXX] + D.da[rXZZYY] * Q.da[qYY] + D.da[rXZZZZ] * Q.da[qZZ]) +
                 2 * (D.da[rXZZXY] * Q.da[qXY] + D.da[rXZZXZ] * Q.da[qXZ] + D.da[rXZZYZ] * Q.da[qYZ]);

  res.da[dYXX] = (D.da[rYXXXX] * Q.da[qXX] + D.da[rYXXYY] * Q.da[qYY] + D.da[rYXXZZ] * Q.da[qZZ]) +
                 2 * (D.da[rYXXXY] * Q.da[qXY] + D.da[rYXXXZ] * Q.da[qXZ] + D.da[rYXXYZ] * Q.da[qYZ]);

  res.da[dYZZ] = (D.da[rYZZXX] * Q.da[qXX] + D.da[rYZZYY] * Q.da[qYY] + D.da[rYZZZZ] * Q.da[qZZ]) +
                 2 * (D.da[rYZZXY] * Q.da[qXY] + D.da[rYZZXZ] * Q.da[qXZ] + D.da[rYZZYZ] * Q.da[qYZ]);

  res.da[dZXX] = (D.da[rZXXXX] * Q.da[qXX] + D.da[rZXXYY] * Q.da[qYY] + D.da[rZXXZZ] * Q.da[qZZ]) +
                 2 * (D.da[rZXXXY] * Q.da[qXY] + D.da[rZXXXZ] * Q.da[qXZ] + D.da[rZXXYZ] * Q.da[qYZ]);

  res.da[dZYY] = (D.da[rZYYXX] * Q.da[qXX] + D.da[rZYYYY] * Q.da[qYY] + D.da[rZYYZZ] * Q.da[qZZ]) +
                 2 * (D.da[rZYYXY] * Q.da[qXY] + D.da[rZYYXZ] * Q.da[qXZ] + D.da[rZYYYZ] * Q.da[qYZ]);

  res.da[dXYZ] = (D.da[rXYZXX] * Q.da[qXX] + D.da[rXYZYY] * Q.da[qYY] + D.da[rXYZZZ] * Q.da[qZZ]) +
                 2 * (D.da[rXYZXY] * Q.da[qXY] + D.da[rXYZXZ] * Q.da[qXZ] + D.da[rXYZYZ] * Q.da[qYZ]);

  return res;
}

// contract a 5-tensor with a vector to yield a 4-tensor
template <typename T1, typename T2>
inline symtensor4<typename which_return<T1, T2>::type> operator*(const symtensor5<T1> &D, const vector<T2> &v)
{
  symtensor4<typename which_return<T1, T2>::type> res;

  res.da[sXXXX] = D.da[rXXXXX] * v.da[0] + D.da[rXXXXY] * v.da[1] + D.da[rXXXXZ] * v.da[2];
  res.da[sXXXY] = D.da[rXXXYX] * v.da[0] + D.da[rXXXYY] * v.da[1] + D.da[rXXXYZ] * v.da[2];
  res.da[sXXXZ] = D.da[rXXXZX] * v.da[0] + D.da[rXXXZY] * v.da[1] + D.da[rXXXZZ] * v.da[2];
  res.da[sXXYY] = D.da[rXXYYX] * v.da[0] + D.da[rXXYYY] * v.da[1] + D.da[rXXYYZ] * v.da[2];
  res.da[sXXYZ] = D.da[rXXYZX] * v.da[0] + D.da[rXXYZY] * v.da[1] + D.da[rXXYZZ] * v.da[2];
  res.da[sXXZZ] = D.da[rXXZZX] * v.da[0] + D.da[rXXZZY] * v.da[1] + D.da[rXXZZZ] * v.da[2];
  res.da[sXYYY] = D.da[rXYYYX] * v.da[0] + D.da[rXYYYY] * v.da[1] + D.da[rXYYYZ] * v.da[2];
  res.da[sXYYZ] = D.da[rXYYZX] * v.da[0] + D.da[rXYYZY] * v.da[1] + D.da[rXYYZZ] * v.da[2];
  res.da[sXYZZ] = D.da[rXYZZX] * v.da[0] + D.da[rXYZZY] * v.da[1] + D.da[rXYZZZ] * v.da[2];
  res.da[sXZZZ] = D.da[rXZZZX] * v.da[0] + D.da[rXZZZY] * v.da[1] + D.da[rXZZZZ] * v.da[2];
  res.da[sYYYY] = D.da[rYYYYX] * v.da[0] + D.da[rYYYYY] * v.da[1] + D.da[rYYYYZ] * v.da[2];
  res.da[sYYYZ] = D.da[rYYYZX] * v.da[0] + D.da[rYYYZY] * v.da[1] + D.da[rYYYZZ] * v.da[2];
  res.da[sYYZZ] = D.da[rYYZZX] * v.da[0] + D.da[rYYZZY] * v.da[1] + D.da[rYYZZZ] * v.da[2];
  res.da[sYZZZ] = D.da[rYZZZX] * v.da[0] + D.da[rYZZZY] * v.da[1] + D.da[rYZZZZ] * v.da[2];
  res.da[sZZZZ] = D.da[rZZZZX] * v.da[0] + D.da[rZZZZY] * v.da[1] + D.da[rZZZZZ] * v.da[2];

  return res;
}

// contract a 5-tensor with a 5-tensor to yield a scalar
template <typename T1, typename T2>
inline typename which_return<T1, T2>::type operator*(const symtensor5<T1> &D, const symtensor5<T2> &Q)
{
  return D.da[rXXXXX] * Q.da[rXXXXX] + D.da[rYYYYY] * Q.da[rYYYYY] + D.da[rZZZZZ] * Q.da[rZZZZZ] +
         5 * (D.da[rYZZZZ] * Q.da[rYZZZZ] + D.da[rXYYYY] * Q.da[rXYYYY] + D.da[rYYYYZ] * Q.da[rYYYYZ] + D.da[rXXXXZ] * Q.da[rXXXXZ] +
              D.da[rXXXXY] * Q.da[rXXXXY] + D.da[rXZZZZ] * Q.da[rXZZZZ]) +
         10 * (D.da[rXXZZZ] * Q.da[rXXZZZ] + D.da[rYYZZZ] * Q.da[rYYZZZ] + D.da[rYYYZZ] * Q.da[rYYYZZ] + D.da[rXXXYY] * Q.da[rXXXYY] +
               D.da[rXXYYY] * Q.da[rXXYYY] + D.da[rXXXZZ] * Q.da[rXXXZZ]) +
         20 * (D.da[rXYYYZ] * Q.da[rXYYYZ] + D.da[rXYZZZ] * Q.da[rXYZZZ] + D.da[rXXXYZ] * Q.da[rXXXYZ]) +
         30 * (D.da[rXYYZZ] * Q.da[rXYYZZ] + D.da[rXXYYZ] * Q.da[rXXYYZ] + D.da[rXXYZZ] * Q.da[rXXYZZ]);
}

template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> contract_twice(const symtensor3<T1> &D, const vector<T2> &v)
{
  typedef typename which_return<T1, T2>::type T;

  symtensor2<T> Dv = D * v;

  vector<T> res = Dv * v;

  return res;
}

template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> contract_thrice(const symtensor4<T1> &D, const vector<T2> &v)
{
  typedef typename which_return<T1, T2>::type T;

  symtensor3<T> Dv  = D * v;
  symtensor2<T> Dvv = Dv * v;

  vector<T> res = Dvv * v;

  return res;
}

template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> contract_fourtimes(const symtensor5<T1> &D, const vector<T2> &v)
{
  typedef typename which_return<T1, T2>::type T;

  symtensor4<T> Dv   = D * v;
  symtensor3<T> Dvv  = Dv * v;
  symtensor2<T> Dvvv = Dvv * v;

  vector<T> res = Dvvv * v;

  return res;
}

// contract a 6-tensor with a vector to yield a 5-tensor
template <typename T1, typename T2>
inline symtensor5<typename which_return<T1, T2>::type> operator*(const symtensor6<T1> &D, const vector<T2> &v)  // checked
{
  symtensor5<typename which_return<T1, T2>::type> res;

  res.da[rXXXXX] = D.da[pXXXXXX] * v.da[0] + D.da[pXXXXXY] * v.da[1] + D.da[pXXXXXZ] * v.da[2];
  res.da[rXXXXY] = D.da[pXXXXXY] * v.da[0] + D.da[pXXXXYY] * v.da[1] + D.da[pXXXXYZ] * v.da[2];
  res.da[rXXXXZ] = D.da[pXXXXXZ] * v.da[0] + D.da[pXXXXYZ] * v.da[1] + D.da[pXXXXZZ] * v.da[2];
  res.da[rXXXYY] = D.da[pXXXXYY] * v.da[0] + D.da[pXXXYYY] * v.da[1] + D.da[pXXXYYZ] * v.da[2];
  res.da[rXXXYZ] = D.da[pXXXXYZ] * v.da[0] + D.da[pXXXYYZ] * v.da[1] + D.da[pXXXYZZ] * v.da[2];
  res.da[rXXXZZ] = D.da[pXXXXZZ] * v.da[0] + D.da[pXXXYZZ] * v.da[1] + D.da[pXXXZZZ] * v.da[2];
  res.da[rXXYYY] = D.da[pXXXYYY] * v.da[0] + D.da[pXXYYYY] * v.da[1] + D.da[pXXYYYZ] * v.da[2];
  res.da[rXXYYZ] = D.da[pXXXYYZ] * v.da[0] + D.da[pXXYYYZ] * v.da[1] + D.da[pXXYYZZ] * v.da[2];
  res.da[rXXYZZ] = D.da[pXXXYZZ] * v.da[0] + D.da[pXXYYZZ] * v.da[1] + D.da[pXXYZZZ] * v.da[2];
  res.da[rXXZZZ] = D.da[pXXXZZZ] * v.da[0] + D.da[pXXYZZZ] * v.da[1] + D.da[pXXZZZZ] * v.da[2];
  res.da[rXYYYY] = D.da[pXXYYYY] * v.da[0] + D.da[pXYYYYY] * v.da[1] + D.da[pXYYYYZ] * v.da[2];
  res.da[rXYYYZ] = D.da[pXXYYYZ] * v.da[0] + D.da[pXYYYYZ] * v.da[1] + D.da[pXYYYZZ] * v.da[2];
  res.da[rXYYZZ] = D.da[pXXYYZZ] * v.da[0] + D.da[pXYYYZZ] * v.da[1] + D.da[pXYYZZZ] * v.da[2];
  res.da[rXYZZZ] = D.da[pXXYZZZ] * v.da[0] + D.da[pXYYZZZ] * v.da[1] + D.da[pXYZZZZ] * v.da[2];
  res.da[rXZZZZ] = D.da[pXXZZZZ] * v.da[0] + D.da[pXYZZZZ] * v.da[1] + D.da[pXZZZZZ] * v.da[2];
  res.da[rYYYYY] = D.da[pXYYYYY] * v.da[0] + D.da[pYYYYYY] * v.da[1] + D.da[pYYYYYZ] * v.da[2];
  res.da[rYYYYZ] = D.da[pXYYYYZ] * v.da[0] + D.da[pYYYYYZ] * v.da[1] + D.da[pYYYYZZ] * v.da[2];
  res.da[rYYYZZ] = D.da[pXYYYZZ] * v.da[0] + D.da[pYYYYZZ] * v.da[1] + D.da[pYYYZZZ] * v.da[2];
  res.da[rYYZZZ] = D.da[pXYYZZZ] * v.da[0] + D.da[pYYYZZZ] * v.da[1] + D.da[pYYZZZZ] * v.da[2];
  res.da[rYZZZZ] = D.da[pXYZZZZ] * v.da[0] + D.da[pYYZZZZ] * v.da[1] + D.da[pYZZZZZ] * v.da[2];
  res.da[rZZZZZ] = D.da[pXZZZZZ] * v.da[0] + D.da[pYZZZZZ] * v.da[1] + D.da[pZZZZZZ] * v.da[2];

  return res;
}

// contract a 7-tensor with a vector to yield a 6-tensor
template <typename T1, typename T2>
inline symtensor6<typename which_return<T1, T2>::type> operator*(const symtensor7<T1> &D, const vector<T2> &v)  // checked
{
  symtensor6<typename which_return<T1, T2>::type> res;

  res.da[pXXXXXX] = D.da[tXXXXXXX] * v.da[0] + D.da[tXXXXXXY] * v.da[1] + D.da[tXXXXXXZ] * v.da[2];
  res.da[pXXXXXY] = D.da[tXXXXXXY] * v.da[0] + D.da[tXXXXXYY] * v.da[1] + D.da[tXXXXXYZ] * v.da[2];
  res.da[pXXXXXZ] = D.da[tXXXXXXZ] * v.da[0] + D.da[tXXXXXYZ] * v.da[1] + D.da[tXXXXXZZ] * v.da[2];
  res.da[pXXXXYY] = D.da[tXXXXXYY] * v.da[0] + D.da[tXXXXYYY] * v.da[1] + D.da[tXXXXYYZ] * v.da[2];
  res.da[pXXXXYZ] = D.da[tXXXXXYZ] * v.da[0] + D.da[tXXXXYYZ] * v.da[1] + D.da[tXXXXYZZ] * v.da[2];
  res.da[pXXXXZZ] = D.da[tXXXXXZZ] * v.da[0] + D.da[tXXXXYZZ] * v.da[1] + D.da[tXXXXZZZ] * v.da[2];
  res.da[pXXXYYY] = D.da[tXXXXYYY] * v.da[0] + D.da[tXXXYYYY] * v.da[1] + D.da[tXXXYYYZ] * v.da[2];
  res.da[pXXXYYZ] = D.da[tXXXXYYZ] * v.da[0] + D.da[tXXXYYYZ] * v.da[1] + D.da[tXXXYYZZ] * v.da[2];
  res.da[pXXXYZZ] = D.da[tXXXXYZZ] * v.da[0] + D.da[tXXXYYZZ] * v.da[1] + D.da[tXXXYZZZ] * v.da[2];
  res.da[pXXXZZZ] = D.da[tXXXXZZZ] * v.da[0] + D.da[tXXXYZZZ] * v.da[1] + D.da[tXXXZZZZ] * v.da[2];
  res.da[pXXYYYY] = D.da[tXXXYYYY] * v.da[0] + D.da[tXXYYYYY] * v.da[1] + D.da[tXXYYYYZ] * v.da[2];
  res.da[pXXYYYZ] = D.da[tXXXYYYZ] * v.da[0] + D.da[tXXYYYYZ] * v.da[1] + D.da[tXXYYYZZ] * v.da[2];
  res.da[pXXYYZZ] = D.da[tXXXYYZZ] * v.da[0] + D.da[tXXYYYZZ] * v.da[1] + D.da[tXXYYZZZ] * v.da[2];
  res.da[pXXYZZZ] = D.da[tXXXYZZZ] * v.da[0] + D.da[tXXYYZZZ] * v.da[1] + D.da[tXXYZZZZ] * v.da[2];
  res.da[pXXZZZZ] = D.da[tXXXZZZZ] * v.da[0] + D.da[tXXYZZZZ] * v.da[1] + D.da[tXXZZZZZ] * v.da[2];
  res.da[pXYYYYY] = D.da[tXXYYYYY] * v.da[0] + D.da[tXYYYYYY] * v.da[1] + D.da[tXYYYYYZ] * v.da[2];
  res.da[pXYYYYZ] = D.da[tXXYYYYZ] * v.da[0] + D.da[tXYYYYYZ] * v.da[1] + D.da[tXYYYYZZ] * v.da[2];
  res.da[pXYYYZZ] = D.da[tXXYYYZZ] * v.da[0] + D.da[tXYYYYZZ] * v.da[1] + D.da[tXYYYZZZ] * v.da[2];
  res.da[pXYYZZZ] = D.da[tXXYYZZZ] * v.da[0] + D.da[tXYYYZZZ] * v.da[1] + D.da[tXYYZZZZ] * v.da[2];
  res.da[pXYZZZZ] = D.da[tXXYZZZZ] * v.da[0] + D.da[tXYYZZZZ] * v.da[1] + D.da[tXYZZZZZ] * v.da[2];
  res.da[pXZZZZZ] = D.da[tXXZZZZZ] * v.da[0] + D.da[tXYZZZZZ] * v.da[1] + D.da[tXZZZZZZ] * v.da[2];
  res.da[pYYYYYY] = D.da[tXYYYYYY] * v.da[0] + D.da[tYYYYYYY] * v.da[1] + D.da[tYYYYYYZ] * v.da[2];
  res.da[pYYYYYZ] = D.da[tXYYYYYZ] * v.da[0] + D.da[tYYYYYYZ] * v.da[1] + D.da[tYYYYYZZ] * v.da[2];
  res.da[pYYYYZZ] = D.da[tXYYYYZZ] * v.da[0] + D.da[tYYYYYZZ] * v.da[1] + D.da[tYYYYZZZ] * v.da[2];
  res.da[pYYYZZZ] = D.da[tXYYYZZZ] * v.da[0] + D.da[tYYYYZZZ] * v.da[1] + D.da[tYYYZZZZ] * v.da[2];
  res.da[pYYZZZZ] = D.da[tXYYZZZZ] * v.da[0] + D.da[tYYYZZZZ] * v.da[1] + D.da[tYYZZZZZ] * v.da[2];
  res.da[pYZZZZZ] = D.da[tXYZZZZZ] * v.da[0] + D.da[tYYZZZZZ] * v.da[1] + D.da[tYZZZZZZ] * v.da[2];
  res.da[pZZZZZZ] = D.da[tXZZZZZZ] * v.da[0] + D.da[tYZZZZZZ] * v.da[1] + D.da[tZZZZZZZ] * v.da[2];

  return res;
}

//-------------  let's define some outer products

/* produce a vector as the cross product of two vectors */
template <typename T1, typename T2>
inline vector<typename which_return<T1, T2>::type> operator^(const vector<T1> &v, const vector<T2> &w)
{
  vector<typename which_return<T1, T2>::type> res;

  res.da[0] = v.da[1] * w.da[2] - v.da[2] * w.da[1];
  res.da[1] = v.da[2] * w.da[0] - v.da[0] * w.da[2];
  res.da[2] = v.da[0] * w.da[1] - v.da[1] * w.da[0];

  return res;
}

/* produce a 2-tensor from the outer product of two vectors */
template <typename T1, typename T2>
inline symtensor2<typename which_return<T1, T2>::type> operator%(const vector<T1> &v, const vector<T2> &w)
{
  symtensor2<typename which_return<T1, T2>::type> res;

  res.da[qXX] = v.da[0] * w.da[0];
  res.da[qYY] = v.da[1] * w.da[1];
  res.da[qZZ] = v.da[2] * w.da[2];
  res.da[qXY] = v.da[0] * w.da[1];
  res.da[qXZ] = v.da[0] * w.da[2];
  res.da[qYZ] = v.da[1] * w.da[2];

  return res;
}

/* produce a 3-tensor from the outer product of a vector and a 2-tensor */
template <typename T1, typename T2>
inline symtensor3<typename which_return<T1, T2>::type> operator%(const vector<T1> &v, const symtensor2<T2> &D)
{
  symtensor3<typename which_return<T1, T2>::type> res;

  res.da[dXXX] = D.da[qXX] * v.da[0];
  res.da[dXXY] = D.da[qXX] * v.da[1];
  res.da[dXXZ] = D.da[qXX] * v.da[2];
  res.da[dXYY] = D.da[qXY] * v.da[1];
  res.da[dXYZ] = D.da[qXY] * v.da[2];
  res.da[dXZZ] = D.da[qXZ] * v.da[2];
  res.da[dYYY] = D.da[qYY] * v.da[1];
  res.da[dYYZ] = D.da[qYY] * v.da[2];
  res.da[dYZZ] = D.da[qYZ] * v.da[2];
  res.da[dZZZ] = D.da[qZZ] * v.da[2];

  return res;
}

/* produce a 4-tensor from the outer product of a vector and a 3-tensor */
template <typename T1, typename T2>
inline symtensor4<typename which_return<T1, T2>::type> operator%(const vector<T1> &v, const symtensor3<T2> &D)
{
  symtensor4<typename which_return<T1, T2>::type> res;

  res.da[sXXXX] = D.da[dXXX] * v.da[0];
  res.da[sXXXY] = D.da[dXXX] * v.da[1];
  res.da[sXXXZ] = D.da[dXXX] * v.da[2];
  res.da[sXXYY] = D.da[dXXY] * v.da[1];
  res.da[sXXYZ] = D.da[dXXY] * v.da[2];
  res.da[sXXZZ] = D.da[dXXZ] * v.da[2];
  res.da[sXYYY] = D.da[dXYY] * v.da[1];
  res.da[sXYYZ] = D.da[dXYY] * v.da[2];
  res.da[sXYZZ] = D.da[dXYZ] * v.da[2];
  res.da[sXZZZ] = D.da[dXZZ] * v.da[2];
  res.da[sYYYY] = D.da[dYYY] * v.da[1];
  res.da[sYYYZ] = D.da[dYYY] * v.da[2];
  res.da[sYYZZ] = D.da[dYYZ] * v.da[2];
  res.da[sYZZZ] = D.da[dYZZ] * v.da[2];
  res.da[sZZZZ] = D.da[dZZZ] * v.da[2];

  return res;
}

/* produce a 5-tensor from the outer product of a vector and a 4-tensor */
template <typename T1, typename T2>
inline symtensor5<typename which_return<T1, T2>::type> operator%(const vector<T1> &v, const symtensor4<T2> &D)
{
  symtensor5<typename which_return<T1, T2>::type> res;

  res.da[rXXXXX] = D.da[sXXXX] * v.da[0];
  res.da[rXXXXY] = D.da[sXXXX] * v.da[1];
  res.da[rXXXXZ] = D.da[sXXXX] * v.da[2];
  res.da[rXXXYY] = D.da[sXXXY] * v.da[1];
  res.da[rXXXYZ] = D.da[sXXXY] * v.da[2];
  res.da[rXXXZZ] = D.da[sXXXZ] * v.da[2];
  res.da[rXXYYY] = D.da[sXXYY] * v.da[1];
  res.da[rXXYYZ] = D.da[sXXYY] * v.da[2];
  res.da[rXXYZZ] = D.da[sXXYZ] * v.da[2];
  res.da[rXXZZZ] = D.da[sXXZZ] * v.da[2];
  res.da[rXYYYY] = D.da[sXYYY] * v.da[1];
  res.da[rXYYYZ] = D.da[sXYYY] * v.da[2];
  res.da[rXYYZZ] = D.da[sXYYZ] * v.da[2];
  res.da[rXYZZZ] = D.da[sXYZZ] * v.da[2];
  res.da[rXZZZZ] = D.da[sXZZZ] * v.da[2];
  res.da[rYYYYY] = D.da[sYYYY] * v.da[1];
  res.da[rYYYYZ] = D.da[sYYYY] * v.da[2];
  res.da[rYYYZZ] = D.da[sYYYZ] * v.da[2];
  res.da[rYYZZZ] = D.da[sYYZZ] * v.da[2];
  res.da[rYZZZZ] = D.da[sYZZZ] * v.da[2];
  res.da[rZZZZZ] = D.da[sZZZZ] * v.da[2];

  return res;
}

/* produce a 6-tensor from the outer product of a vector and a 5-tensor */
template <typename T1, typename T2>
inline symtensor6<typename which_return<T1, T2>::type> operator%(const vector<T1> &v, const symtensor5<T2> &D)
{
  symtensor6<typename which_return<T1, T2>::type> res;

  res.da[pXXXXXX] = D.da[rXXXXX] * v.da[0];
  res.da[pXXXXXY] = D.da[rXXXXX] * v.da[1];
  res.da[pXXXXXZ] = D.da[rXXXXX] * v.da[2];
  res.da[pXXXXYY] = D.da[rXXXXY] * v.da[1];
  res.da[pXXXXYZ] = D.da[rXXXXY] * v.da[2];
  res.da[pXXXXZZ] = D.da[rXXXXZ] * v.da[2];
  res.da[pXXXYYY] = D.da[rXXXYY] * v.da[1];
  res.da[pXXXYYZ] = D.da[rXXXYY] * v.da[2];
  res.da[pXXXYZZ] = D.da[rXXXYZ] * v.da[2];
  res.da[pXXXZZZ] = D.da[rXXXZZ] * v.da[2];
  res.da[pXXYYYY] = D.da[rXXYYY] * v.da[1];
  res.da[pXXYYYZ] = D.da[rXXYYY] * v.da[2];
  res.da[pXXYYZZ] = D.da[rXXYYZ] * v.da[2];
  res.da[pXXYZZZ] = D.da[rXXYZZ] * v.da[2];
  res.da[pXXZZZZ] = D.da[rXXZZZ] * v.da[2];
  res.da[pXYYYYY] = D.da[rXYYYY] * v.da[1];
  res.da[pXYYYYZ] = D.da[rXYYYY] * v.da[2];
  res.da[pXYYYZZ] = D.da[rXYYYZ] * v.da[2];
  res.da[pXYYZZZ] = D.da[rXYYZZ] * v.da[2];
  res.da[pXYZZZZ] = D.da[rXYZZZ] * v.da[2];
  res.da[pXZZZZZ] = D.da[rXZZZZ] * v.da[2];
  res.da[pYYYYYY] = D.da[rYYYYY] * v.da[1];
  res.da[pYYYYYZ] = D.da[rYYYYY] * v.da[2];
  res.da[pYYYYZZ] = D.da[rYYYYZ] * v.da[2];
  res.da[pYYYZZZ] = D.da[rYYYZZ] * v.da[2];
  res.da[pYYZZZZ] = D.da[rYYZZZ] * v.da[2];
  res.da[pYZZZZZ] = D.da[rYZZZZ] * v.da[2];
  res.da[pZZZZZZ] = D.da[rZZZZZ] * v.da[2];

  return res;
}

/* produce a 7-tensor from the outer product of a vector and a 6-tensor */
template <typename T1, typename T2>
inline symtensor7<typename which_return<T1, T2>::type> operator%(const vector<T1> &v, const symtensor6<T2> &D)
{
  symtensor7<typename which_return<T1, T2>::type> res;

  res.da[tXXXXXXX] = D.da[pXXXXXX] * v.da[0];
  res.da[tXXXXXXY] = D.da[pXXXXXX] * v.da[1];
  res.da[tXXXXXXZ] = D.da[pXXXXXX] * v.da[2];
  res.da[tXXXXXYY] = D.da[pXXXXXY] * v.da[1];
  res.da[tXXXXXYZ] = D.da[pXXXXXY] * v.da[2];
  res.da[tXXXXXZZ] = D.da[pXXXXXZ] * v.da[2];
  res.da[tXXXXYYY] = D.da[pXXXXYY] * v.da[1];
  res.da[tXXXXYYZ] = D.da[pXXXXYY] * v.da[2];
  res.da[tXXXXYZZ] = D.da[pXXXXYZ] * v.da[2];
  res.da[tXXXXZZZ] = D.da[pXXXXZZ] * v.da[2];
  res.da[tXXXYYYY] = D.da[pXXXYYY] * v.da[1];
  res.da[tXXXYYYZ] = D.da[pXXXYYY] * v.da[2];
  res.da[tXXXYYZZ] = D.da[pXXXYYZ] * v.da[2];
  res.da[tXXXYZZZ] = D.da[pXXXYZZ] * v.da[2];
  res.da[tXXXZZZZ] = D.da[pXXXZZZ] * v.da[2];
  res.da[tXXYYYYY] = D.da[pXXYYYY] * v.da[1];
  res.da[tXXYYYYZ] = D.da[pXXYYYY] * v.da[2];
  res.da[tXXYYYZZ] = D.da[pXXYYYZ] * v.da[2];
  res.da[tXXYYZZZ] = D.da[pXXYYZZ] * v.da[2];
  res.da[tXXYZZZZ] = D.da[pXXYZZZ] * v.da[2];
  res.da[tXXZZZZZ] = D.da[pXXZZZZ] * v.da[2];
  res.da[tXYYYYYY] = D.da[pXYYYYY] * v.da[1];
  res.da[tXYYYYYZ] = D.da[pXYYYYY] * v.da[2];
  res.da[tXYYYYZZ] = D.da[pXYYYYZ] * v.da[2];
  res.da[tXYYYZZZ] = D.da[pXYYYZZ] * v.da[2];
  res.da[tXYYZZZZ] = D.da[pXYYZZZ] * v.da[2];
  res.da[tXYZZZZZ] = D.da[pXYZZZZ] * v.da[2];
  res.da[tXZZZZZZ] = D.da[pXZZZZZ] * v.da[2];
  res.da[tYYYYYYY] = D.da[pYYYYYY] * v.da[1];
  res.da[tYYYYYYZ] = D.da[pYYYYYY] * v.da[2];
  res.da[tYYYYYZZ] = D.da[pYYYYYZ] * v.da[2];
  res.da[tYYYYZZZ] = D.da[pYYYYZZ] * v.da[2];
  res.da[tYYYZZZZ] = D.da[pYYYZZZ] * v.da[2];
  res.da[tYYZZZZZ] = D.da[pYYZZZZ] * v.da[2];
  res.da[tYZZZZZZ] = D.da[pYZZZZZ] * v.da[2];
  res.da[tZZZZZZZ] = D.da[pZZZZZZ] * v.da[2];

  return res;
}

/* compute the sum of the three possible outer products of a 2-tensor and a vector, yielding a symmetric 3-tensor */
template <typename T1, typename T2>
inline symtensor3<typename which_return<T1, T2>::type> outer_prod_sum(const symtensor2<T1> &D, const vector<T2> &v)
{
  symtensor3<typename which_return<T1, T2>::type> res;

  res.da[dXXX] = 3 * D.da[qXX] * v.da[0];
  res.da[dYYY] = 3 * D.da[qYY] * v.da[1];
  res.da[dZZZ] = 3 * D.da[qZZ] * v.da[2];
  res.da[dXYY] = 2 * D.da[qXY] * v.da[1] + D.da[qYY] * v.da[0];
  res.da[dXZZ] = 2 * D.da[qXZ] * v.da[2] + D.da[qZZ] * v.da[0];
  res.da[dYXX] = 2 * D.da[qXY] * v.da[0] + D.da[qXX] * v.da[1];
  res.da[dYZZ] = 2 * D.da[qYZ] * v.da[2] + D.da[qZZ] * v.da[1];
  res.da[dZXX] = 2 * D.da[qXZ] * v.da[0] + D.da[qXX] * v.da[2];
  res.da[dZYY] = 2 * D.da[qYZ] * v.da[1] + D.da[qYY] * v.da[2];
  res.da[dXYZ] = D.da[qXY] * v.da[2] + D.da[qYZ] * v.da[0] + D.da[qZX] * v.da[1];

  return res;
}

/* compute the sum of the four possible outer products of a 3-tensor and a vector, yielding a symmetric 4-tensor */
template <typename T1, typename T2>
inline symtensor4<typename which_return<T1, T2>::type> outer_prod_sum(const symtensor3<T1> &D, const vector<T2> &v)
{
  symtensor4<typename which_return<T1, T2>::type> res;

  res.da[sXXXX] = 4 * D.da[dXXX] * v.da[0];
  res.da[sZZZZ] = 4 * D.da[dZZZ] * v.da[2];
  res.da[sYYYY] = 4 * D.da[dYYY] * v.da[1];
  res.da[sXXXY] = 3 * D.da[dXXY] * v.da[0] + D.da[dXXX] * v.da[1];
  res.da[sXXXZ] = 3 * D.da[dXXZ] * v.da[0] + D.da[dXXX] * v.da[2];
  res.da[sXYYY] = 3 * D.da[dXYY] * v.da[1] + D.da[dYYY] * v.da[0];
  res.da[sYYYZ] = 3 * D.da[dYYZ] * v.da[1] + D.da[dYYY] * v.da[2];
  res.da[sXZZZ] = 3 * D.da[dXZZ] * v.da[2] + D.da[dZZZ] * v.da[0];
  res.da[sYZZZ] = 3 * D.da[dYZZ] * v.da[2] + D.da[dZZZ] * v.da[1];
  res.da[sXXYY] = 2 * D.da[dXXY] * v.da[1] + 2 * D.da[dXYY] * v.da[0];
  res.da[sXXZZ] = 2 * D.da[dXXZ] * v.da[2] + 2 * D.da[dXZZ] * v.da[0];
  res.da[sYYZZ] = 2 * D.da[dYYZ] * v.da[2] + 2 * D.da[dYZZ] * v.da[1];
  res.da[sXXYZ] = 2 * D.da[dXYZ] * v.da[0] + D.da[dXXY] * v.da[2] + D.da[dXXZ] * v.da[1];
  res.da[sXYYZ] = 2 * D.da[dXYZ] * v.da[1] + D.da[dXYY] * v.da[2] + D.da[dYYZ] * v.da[0];
  res.da[sXYZZ] = 2 * D.da[dXYZ] * v.da[2] + D.da[dXZZ] * v.da[1] + D.da[dYZZ] * v.da[0];

  return res;
}

/* compute the sum of the six possible outer products of a 2-tensor with another 2-tensor, yielding a symmetric 4-tensor */
template <typename T1, typename T2>
inline symtensor4<typename which_return<T1, T2>::type> outer_prod_sum(const symtensor2<T1> &D, const symtensor2<T2> &S)
{
  symtensor4<typename which_return<T1, T2>::type> res;

  res.da[sXXXX] = 6 * D.da[qXX] * S.da[qXX];
  res.da[sYYYY] = 6 * D.da[qYY] * S.da[qYY];
  res.da[sZZZZ] = 6 * D.da[qZZ] * S.da[qZZ];
  res.da[sXXXY] = 3 * D.da[qXX] * S.da[qXY] + 3 * D.da[qXY] * S.da[qXX];
  res.da[sXXXZ] = 3 * D.da[qXX] * S.da[qXZ] + 3 * D.da[qXZ] * S.da[qXX];
  res.da[sXYYY] = 3 * D.da[qXY] * S.da[qYY] + 3 * D.da[qYY] * S.da[qXY];
  res.da[sXZZZ] = 3 * D.da[qZZ] * S.da[qXZ] + 3 * D.da[qXZ] * S.da[qZZ];
  res.da[sYYYZ] = 3 * D.da[qYY] * S.da[qYZ] + 3 * D.da[qYZ] * S.da[qYY];
  res.da[sYZZZ] = 3 * D.da[qZZ] * S.da[qYZ] + 3 * D.da[qYZ] * S.da[qZZ];
  res.da[sXXYY] = D.da[qXX] * S.da[qYY] + D.da[qYY] * S.da[qXX] + 4 * D.da[qXY] * S.da[qXY];
  res.da[sXXZZ] = D.da[qXX] * S.da[qZZ] + D.da[qZZ] * S.da[qXX] + 4 * D.da[qXZ] * S.da[qXZ];
  res.da[sYYZZ] = D.da[qYY] * S.da[qZZ] + D.da[qZZ] * S.da[qYY] + 4 * D.da[qYZ] * S.da[qYZ];
  res.da[sXXYZ] = D.da[qXX] * S.da[qYZ] + 2 * D.da[qXY] * S.da[qXZ] + 2 * D.da[qXZ] * S.da[qXY] + D.da[qYZ] * S.da[qXX];
  res.da[sXYZZ] = D.da[qZZ] * S.da[qXY] + 2 * D.da[qXZ] * S.da[qYZ] + 2 * D.da[qYZ] * S.da[qXZ] + D.da[qXY] * S.da[qZZ];
  res.da[sXYYZ] = D.da[qYY] * S.da[qXZ] + 2 * D.da[qXY] * S.da[qYZ] + 2 * D.da[qYZ] * S.da[qXY] + D.da[qXZ] * S.da[qYY];

  return res;
}

/* compute the sum of the five possible outer products of a 4-tensor and a vector, yielding a symmetric 5-tensor */
template <typename T1, typename T2>
inline symtensor5<typename which_return<T1, T2>::type> outer_prod_sum(const symtensor4<T1> &D, const vector<T2> &v)
{
  symtensor5<typename which_return<T1, T2>::type> res;

  res.da[rXXXXX] = 5 * D.da[sXXXX] * v.da[0];
  res.da[rXXXXY] = 4 * D.da[sXXXY] * v.da[0] + D.da[sXXXX] * v.da[1];
  res.da[rXXXXZ] = 4 * D.da[sXXXZ] * v.da[0] + D.da[sXXXX] * v.da[2];
  res.da[rXXXYX] = 3 * D.da[sXXYX] * v.da[0] + D.da[sXXXX] * v.da[1] + D.da[sXXXY] * v.da[0];
  res.da[rXXXYY] = 3 * D.da[sXXYY] * v.da[0] + 2 * D.da[sXXXY] * v.da[1];
  res.da[rXXXYZ] = 3 * D.da[sXXYZ] * v.da[0] + D.da[sXXXZ] * v.da[1] + D.da[sXXXY] * v.da[2];
  res.da[rXXXZX] = 3 * D.da[sXXZX] * v.da[0] + D.da[sXXXX] * v.da[2] + D.da[sXXXZ] * v.da[0];
  res.da[rXXXZY] = 3 * D.da[sXXZY] * v.da[0] + D.da[sXXXY] * v.da[2] + D.da[sXXXZ] * v.da[1];
  res.da[rXXXZZ] = 3 * D.da[sXXZZ] * v.da[0] + 2 * D.da[sXXXZ] * v.da[2];
  res.da[rXXYYY] = 2 * D.da[sXYYY] * v.da[0] + 3 * D.da[sXXYY] * v.da[1];
  res.da[rXXYYZ] = 2 * D.da[sXYYZ] * v.da[0] + 2 * D.da[sXXYZ] * v.da[1] + D.da[sXXYY] * v.da[2];
  res.da[rXXYZY] = 2 * D.da[sXYZY] * v.da[0] + D.da[sXXZY] * v.da[1] + D.da[sXXYY] * v.da[2] + D.da[sXXYZ] * v.da[1];
  res.da[rXXYZZ] = 2 * D.da[sXYZZ] * v.da[0] + D.da[sXXZZ] * v.da[1] + 2 * D.da[sXXYZ] * v.da[2];
  res.da[rXXZZZ] = 2 * D.da[sXZZZ] * v.da[0] + 3 * D.da[sXXZZ] * v.da[2];
  res.da[rXYYYY] = D.da[sYYYY] * v.da[0] + 4 * D.da[sXYYY] * v.da[1];
  res.da[rXYYYZ] = D.da[sYYYZ] * v.da[0] + 3 * D.da[sXYYZ] * v.da[1] + D.da[sXYYY] * v.da[2];
  res.da[rXYYZY] = D.da[sYYZY] * v.da[0] + 2 * D.da[sXYZY] * v.da[1] + D.da[sXYYY] * v.da[2] + D.da[sXYYZ] * v.da[1];
  res.da[rXYYZZ] = D.da[sYYZZ] * v.da[0] + 2 * D.da[sXYZZ] * v.da[1] + 2 * D.da[sXYYZ] * v.da[2];
  res.da[rXYZZZ] = D.da[sYZZZ] * v.da[0] + D.da[sXZZZ] * v.da[1] + 3 * D.da[sXYZZ] * v.da[2];
  res.da[rXZZZZ] = D.da[sZZZZ] * v.da[0] + 4 * D.da[sXZZZ] * v.da[2];
  res.da[rYYYYY] = 5 * D.da[sYYYY] * v.da[1];
  res.da[rYYYYZ] = 4 * D.da[sYYYZ] * v.da[1] + D.da[sYYYY] * v.da[2];
  res.da[rYYYZY] = 3 * D.da[sYYZY] * v.da[1] + D.da[sYYYY] * v.da[2] + D.da[sYYYZ] * v.da[1];
  res.da[rYYYZZ] = 3 * D.da[sYYZZ] * v.da[1] + 2 * D.da[sYYYZ] * v.da[2];
  res.da[rYYZZZ] = 2 * D.da[sYZZZ] * v.da[1] + 3 * D.da[sYYZZ] * v.da[2];
  res.da[rYZZZZ] = D.da[sZZZZ] * v.da[1] + 4 * D.da[sYZZZ] * v.da[2];
  res.da[rZZZZZ] = 5 * D.da[sZZZZ] * v.da[2];

  return res;
}

/* compute the sum of the 10 possible outer products of a 3-tensor with another 2-tensor, yielding a symmetric 5-tensor */
template <typename T1, typename T2>
inline symtensor5<typename which_return<T1, T2>::type> outer_prod_sum(const symtensor3<T1> &D, const symtensor2<T2> &S)
{
  symtensor5<typename which_return<T1, T2>::type> res;

  res.da[rXXXXX] = 10 * D.da[dXXX] * S.da[qXX];
  res.da[rXXXXY] = 6 * D.da[dXXY] * S.da[qXX] + 4 * D.da[dXXX] * S.da[qXY];
  res.da[rXXXXZ] = 6 * D.da[dXXZ] * S.da[qXX] + 4 * D.da[dXXX] * S.da[qXZ];
  res.da[rXXXYX] = 3 * D.da[dXYX] * S.da[qXX] + 3 * D.da[dXXX] * S.da[qXY] + 3 * D.da[dXXY] * S.da[qXX] + D.da[dXXX] * S.da[qYX];
  res.da[rXXXYY] = 3 * D.da[dXYY] * S.da[qXX] + 6 * D.da[dXXY] * S.da[qXY] + D.da[dXXX] * S.da[qYY];
  res.da[rXXXYZ] = 3 * D.da[dXYZ] * S.da[qXX] + 3 * D.da[dXXZ] * S.da[qXY] + 3 * D.da[dXXY] * S.da[qXZ] + D.da[dXXX] * S.da[qYZ];
  res.da[rXXXZX] = 3 * D.da[dXZX] * S.da[qXX] + 3 * D.da[dXXX] * S.da[qXZ] + 3 * D.da[dXXZ] * S.da[qXX] + D.da[dXXX] * S.da[qZX];
  res.da[rXXXZY] = 3 * D.da[dXZY] * S.da[qXX] + 3 * D.da[dXXY] * S.da[qXZ] + 3 * D.da[dXXZ] * S.da[qXY] + D.da[dXXX] * S.da[qZY];
  res.da[rXXXZZ] = 3 * D.da[dXZZ] * S.da[qXX] + 6 * D.da[dXXZ] * S.da[qXZ] + D.da[dXXX] * S.da[qZZ];
  res.da[rXXYYY] = D.da[dYYY] * S.da[qXX] + 6 * D.da[dXYY] * S.da[qXY] + 3 * D.da[dXXY] * S.da[qYY];
  res.da[rXXYYZ] = D.da[dYYZ] * S.da[qXX] + 4 * D.da[dXYZ] * S.da[qXY] + 2 * D.da[dXYY] * S.da[qXZ] + D.da[dXXZ] * S.da[qYY] +
                   2 * D.da[dXXY] * S.da[qYZ];
  res.da[rXXYZY] = D.da[dYZY] * S.da[qXX] + 2 * D.da[dXZY] * S.da[qXY] + 2 * D.da[dXYY] * S.da[qXZ] + 2 * D.da[dXYZ] * S.da[qXY] +
                   D.da[dXXY] * S.da[qYZ] + D.da[dXXZ] * S.da[qYY] + D.da[dXXY] * S.da[qZY];
  res.da[rXXYZZ] = D.da[dYZZ] * S.da[qXX] + 2 * D.da[dXZZ] * S.da[qXY] + 4 * D.da[dXYZ] * S.da[qXZ] + 2 * D.da[dXXZ] * S.da[qYZ] +
                   D.da[dXXY] * S.da[qZZ];
  res.da[rXXZZZ] = D.da[dZZZ] * S.da[qXX] + 6 * D.da[dXZZ] * S.da[qXZ] + 3 * D.da[dXXZ] * S.da[qZZ];
  res.da[rXYYYY] = 4 * D.da[dYYY] * S.da[qXY] + 6 * D.da[dXYY] * S.da[qYY];
  res.da[rXYYYZ] = 3 * D.da[dYYZ] * S.da[qXY] + D.da[dYYY] * S.da[qXZ] + 3 * D.da[dXYZ] * S.da[qYY] + 3 * D.da[dXYY] * S.da[qYZ];
  res.da[rXYYZY] = 2 * D.da[dYZY] * S.da[qXY] + D.da[dYYY] * S.da[qXZ] + D.da[dYYZ] * S.da[qXY] + D.da[dXZY] * S.da[qYY] +
                   2 * D.da[dXYY] * S.da[qYZ] + 2 * D.da[dXYZ] * S.da[qYY] + D.da[dXYY] * S.da[qZY];
  res.da[rXYYZZ] = 2 * D.da[dYZZ] * S.da[qXY] + 2 * D.da[dYYZ] * S.da[qXZ] + D.da[dXZZ] * S.da[qYY] + 4 * D.da[dXYZ] * S.da[qYZ] +
                   D.da[dXYY] * S.da[qZZ];
  res.da[rXYZZZ] = D.da[dZZZ] * S.da[qXY] + 3 * D.da[dYZZ] * S.da[qXZ] + 3 * D.da[dXZZ] * S.da[qYZ] + 3 * D.da[dXYZ] * S.da[qZZ];
  res.da[rXZZZZ] = 4 * D.da[dZZZ] * S.da[qXZ] + 6 * D.da[dXZZ] * S.da[qZZ];
  res.da[rYYYYY] = 10 * D.da[dYYY] * S.da[qYY];
  res.da[rYYYYZ] = 6 * D.da[dYYZ] * S.da[qYY] + 4 * D.da[dYYY] * S.da[qYZ];
  res.da[rYYYZY] = 3 * D.da[dYZY] * S.da[qYY] + 3 * D.da[dYYY] * S.da[qYZ] + 3 * D.da[dYYZ] * S.da[qYY] + D.da[dYYY] * S.da[qZY];
  res.da[rYYYZZ] = 3 * D.da[dYZZ] * S.da[qYY] + 6 * D.da[dYYZ] * S.da[qYZ] + D.da[dYYY] * S.da[qZZ];
  res.da[rYYZZZ] = D.da[dZZZ] * S.da[qYY] + 6 * D.da[dYZZ] * S.da[qYZ] + 3 * D.da[dYYZ] * S.da[qZZ];
  res.da[rYZZZZ] = 4 * D.da[dZZZ] * S.da[qYZ] + 6 * D.da[dYZZ] * S.da[qZZ];
  res.da[rZZZZZ] = 10 * D.da[dZZZ] * S.da[qZZ];

  return res;
}

enum setup_options
{
  INIT,
  ADD
};

template <typename T, typename TypeGfac>
inline void setup_D3(enum setup_options opt, symtensor3<T> &D3, vector<T> &dxyz, symtensor2<T> &aux2, symtensor3<T> &aux3, TypeGfac g2,
                     TypeGfac g3)
{
  // Note: dxyz, aux2  are input parameters, whereas aux3 is an output parameter!

  aux3 = dxyz % aux2;  // construct outer product of the two vectors
  if(opt == INIT)
    D3 = g3 * aux3;
  else
    D3 += g3 * aux3;

  vector<T> g2_dxyz = g2 * dxyz;

  D3[dXXX] += 3 * g2_dxyz[0];
  D3[dYYY] += 3 * g2_dxyz[1];
  D3[dZZZ] += 3 * g2_dxyz[2];

  D3[dXXY] += g2_dxyz[1];
  D3[dXXZ] += g2_dxyz[2];
  D3[dXYY] += g2_dxyz[0];
  D3[dXZZ] += g2_dxyz[0];
  D3[dYYZ] += g2_dxyz[2];
  D3[dYZZ] += g2_dxyz[1];
}

template <typename T, typename TypeGfac>
inline void setup_D4(enum setup_options opt, symtensor4<T> &D4, vector<T> &dxyz, symtensor2<T> &aux2, symtensor3<T> &aux3,
                     symtensor4<T> &aux4, TypeGfac g2, TypeGfac g3, TypeGfac g4)
{
  // Note: dxyz, aux2, and aux3 are input parameters, whereas aux4 is an output parameter!

  aux4 = dxyz % aux3;  // construct outer product
  if(opt == INIT)
    D4 = g4 * aux4;
  else
    D4 += g4 * aux4;

  D4[sXXXX] += 3 * g2;
  D4[sYYYY] += 3 * g2;
  D4[sZZZZ] += 3 * g2;
  D4[sXXYY] += g2;
  D4[sXXZZ] += g2;
  D4[sYYZZ] += g2;

  symtensor2<T> g3aux2 = g3 * aux2;

  D4[sXXXX] += 6 * g3aux2[qXX];
  D4[sYYYY] += 6 * g3aux2[qYY];
  D4[sZZZZ] += 6 * g3aux2[qZZ];

  D4[sXXXY] += 3 * g3aux2[qXY];
  D4[sXYYY] += 3 * g3aux2[qXY];
  D4[sXXXZ] += 3 * g3aux2[qXZ];
  D4[sXZZZ] += 3 * g3aux2[qXZ];
  D4[sYYYZ] += 3 * g3aux2[qYZ];
  D4[sYZZZ] += 3 * g3aux2[qYZ];

  D4[sXXYY] += g3aux2[qXX] + g3aux2[qYY];
  D4[sXXZZ] += g3aux2[qXX] + g3aux2[qZZ];
  D4[sYYZZ] += g3aux2[qYY] + g3aux2[qZZ];

  D4[sXXYZ] += g3aux2[qYZ];
  D4[sXYYZ] += g3aux2[qXZ];
  D4[sXYZZ] += g3aux2[qXY];
}

template <typename T, typename TypeGfac>
inline void setup_D5(enum setup_options opt, symtensor5<T> &D5, vector<T> &dxyz, symtensor3<T> &aux3, symtensor4<T> &aux4,
                     symtensor5<T> &aux5, TypeGfac g3, TypeGfac g4, TypeGfac g5)
{
  // Note: dxyz, aux3, and aux4 are input parameters, whereas aux5 is an output parameter!

  aux5 = dxyz % aux4;  // construct outer product
  if(opt == INIT)
    D5 = g5 * aux5;
  else
    D5 += g5 * aux5;

  vector<T> g3_dxyz = g3 * dxyz;

  D5[rXXXXX] += 15 * g3_dxyz[0];
  D5[rYYYYY] += 15 * g3_dxyz[1];
  D5[rZZZZZ] += 15 * g3_dxyz[2];

  D5[rXXXXY] += 3 * g3_dxyz[1];
  D5[rXXXXZ] += 3 * g3_dxyz[2];
  D5[rXYYYY] += 3 * g3_dxyz[0];
  D5[rXZZZZ] += 3 * g3_dxyz[0];
  D5[rYYYYZ] += 3 * g3_dxyz[2];
  D5[rYZZZZ] += 3 * g3_dxyz[1];

  D5[rXXXYY] += 3 * g3_dxyz[0];
  D5[rXXXZZ] += 3 * g3_dxyz[0];
  D5[rXXYYY] += 3 * g3_dxyz[1];
  D5[rXXZZZ] += 3 * g3_dxyz[2];
  D5[rYYYZZ] += 3 * g3_dxyz[1];
  D5[rYYZZZ] += 3 * g3_dxyz[2];

  D5[rXXYZZ] += g3_dxyz[1];
  D5[rXXYYZ] += g3_dxyz[2];
  D5[rXYYZZ] += g3_dxyz[0];

  D5[rXXXYZ] += 0;
  D5[rXYYYZ] += 0;
  D5[rXYZZZ] += 0;

  ///// ----
  symtensor3<T> g4aux3 = g4 * aux3;

  D5[rXXXXX] += 10 * g4aux3[dXXX];
  D5[rYYYYY] += 10 * g4aux3[dYYY];
  D5[rZZZZZ] += 10 * g4aux3[dZZZ];

  D5[rXXXXY] += 6 * g4aux3[dXXY];
  D5[rXXXXZ] += 6 * g4aux3[dXXZ];
  D5[rXYYYY] += 6 * g4aux3[dYYX];
  D5[rXZZZZ] += 6 * g4aux3[dZZX];
  D5[rYYYYZ] += 6 * g4aux3[dYYZ];
  D5[rYZZZZ] += 6 * g4aux3[dZZY];

  D5[rXXXYY] += g4aux3[dXXX] + 3 * g4aux3[dXYY];
  D5[rXXXZZ] += g4aux3[dXXX] + 3 * g4aux3[dXZZ];
  D5[rXXYYY] += g4aux3[dYYY] + 3 * g4aux3[dYXX];
  D5[rXXZZZ] += g4aux3[dZZZ] + 3 * g4aux3[dZXX];
  D5[rYYYZZ] += g4aux3[dYYY] + 3 * g4aux3[dYZZ];
  D5[rYYZZZ] += g4aux3[dZZZ] + 3 * g4aux3[dZYY];

  D5[rXXYZZ] += g4aux3[dYZZ] + g4aux3[dXXY];
  D5[rXXYYZ] += g4aux3[dYYZ] + g4aux3[dXXZ];
  D5[rXYYZZ] += g4aux3[dXZZ] + g4aux3[dXYY];

  D5[rXXXYZ] += 3 * g4aux3[dXYZ];
  D5[rXYYYZ] += 3 * g4aux3[dXYZ];
  D5[rXYZZZ] += 3 * g4aux3[dXYZ];
}

template <typename T, typename TypeGfac>
inline void setup_D6(enum setup_options opt, symtensor6<T> &D6, vector<T> &dxyz, TypeGfac g3, TypeGfac g4, TypeGfac g5, TypeGfac g6)
{
#define X 0
#define Y 1
#define Z 2

  if(opt == INIT)
    D6 = static_cast<T>(0);

  D6[pXXXXXX] += 15 * g3 + g4 * (45 * dxyz[X] * dxyz[X]) + g5 * (15 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X];

  D6[pXXXXXY] += g4 * (15 * dxyz[X] * dxyz[Y]) + g5 * (10 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y];

  D6[pXXXXXZ] += g4 * (15 * dxyz[X] * dxyz[Z]) + g5 * (10 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z];

  D6[pXXXXYY] += 3 * g3 + g4 * (6 * dxyz[X] * dxyz[X] + 3 * dxyz[Y] * dxyz[Y]) +
                 g5 * (6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y];

  D6[pXXXXYZ] += g4 * (3 * dxyz[Y] * dxyz[Z]) + g5 * (6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z];

  D6[pXXXXZZ] += 3 * g3 + g4 * (6 * dxyz[X] * dxyz[X] + 3 * dxyz[Z] * dxyz[Z]) +
                 g5 * (6 * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z];

  D6[pXXXYYY] += g4 * (9 * dxyz[X] * dxyz[Y]) +
                 g5 * (3 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] + 3 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D6[pXXXYYZ] += g4 * (3 * dxyz[X] * dxyz[Z]) +
                 g5 * (3 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D6[pXXXYZZ] += g4 * (3 * dxyz[X] * dxyz[Y]) +
                 g5 * (3 * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D6[pXXXZZZ] += g4 * (9 * dxyz[X] * dxyz[Z]) +
                 g5 * (3 * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 3 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pXXYYYY] += 3 * g3 + g4 * (3 * dxyz[X] * dxyz[X] + 6 * dxyz[Y] * dxyz[Y]) +
                 g5 * (dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] + 6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D6[pXXYYYZ] += g4 * (3 * dxyz[Y] * dxyz[Z]) +
                 g5 * (dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] + 3 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D6[pXXYYZZ] +=
      1 * g3 + g4 * (dxyz[X] * dxyz[X] + dxyz[Y] * dxyz[Y] + dxyz[Z] * dxyz[Z]) +
      g5 * (dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y]) +
      g6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D6[pXXYZZZ] += g4 * (3 * dxyz[Y] * dxyz[Z]) +
                 g5 * (dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 3 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pXXZZZZ] += 3 * g3 + g4 * (3 * dxyz[X] * dxyz[X] + 6 * dxyz[Z] * dxyz[Z]) +
                 g5 * (dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 6 * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pXYYYYY] += g4 * (15 * dxyz[X] * dxyz[Y]) + g5 * (10 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                 g6 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D6[pXYYYYZ] += g4 * (3 * dxyz[X] * dxyz[Z]) + g5 * (6 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D6[pXYYYZZ] += g4 * (3 * dxyz[X] * dxyz[Y]) +
                 g5 * (3 * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                 g6 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D6[pXYYZZZ] += g4 * (3 * dxyz[X] * dxyz[Z]) +
                 g5 * (dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 3 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pXYZZZZ] += g4 * (3 * dxyz[X] * dxyz[Y]) + g5 * (6 * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pXZZZZZ] += g4 * (15 * dxyz[X] * dxyz[Z]) + g5 * (10 * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                 g6 * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pYYYYYY] += 15 * g3 + g4 * (45 * dxyz[Y] * dxyz[Y]) + g5 * (15 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                 g6 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D6[pYYYYYZ] += g4 * (15 * dxyz[Y] * dxyz[Z]) + g5 * (10 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                 g6 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D6[pYYYYZZ] += 3 * g3 + g4 * (6 * dxyz[Y] * dxyz[Y] + 3 * dxyz[Z] * dxyz[Z]) +
                 g5 * (6 * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] + dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                 g6 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D6[pYYYZZZ] += g4 * (9 * dxyz[Y] * dxyz[Z]) +
                 g5 * (3 * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 3 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                 g6 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pYYZZZZ] += 3 * g3 + g4 * (3 * dxyz[Y] * dxyz[Y] + 6 * dxyz[Z] * dxyz[Z]) +
                 g5 * (dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 6 * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                 g6 * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pYZZZZZ] += g4 * (15 * dxyz[Y] * dxyz[Z]) + g5 * (10 * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                 g6 * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D6[pZZZZZZ] += 15 * g3 + g4 * (45 * dxyz[Z] * dxyz[Z]) + g5 * (15 * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                 g6 * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

#undef X
#undef Y
#undef Z
}

template <typename T, typename TypeGfac>
inline void setup_D7(enum setup_options opt, symtensor7<T> &D7, vector<T> &dxyz, TypeGfac g4, TypeGfac g5, TypeGfac g6, TypeGfac g7)
{
#define X 0
#define Y 1
#define Z 2

  if(opt == INIT)
    D7 = static_cast<T>(0);

  D7[tXXXXXXX] += g4 * (105 * dxyz[X]) + g5 * (105 * dxyz[X] * dxyz[X] * dxyz[X]) +
                  g6 * (21 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X];

  D7[tXXXXXXY] += g4 * (15 * dxyz[Y]) + g5 * (45 * dxyz[X] * dxyz[X] * dxyz[Y]) +
                  g6 * (15 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y];

  D7[tXXXXXXZ] += g4 * (15 * dxyz[Z]) + g5 * (45 * dxyz[X] * dxyz[X] * dxyz[Z]) +
                  g6 * (15 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z];

  D7[tXXXXXYY] += g4 * (15 * dxyz[X]) + g5 * (10 * dxyz[X] * dxyz[X] * dxyz[X] + 15 * dxyz[X] * dxyz[Y] * dxyz[Y]) +
                  g6 * (10 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y];

  D7[tXXXXXYZ] += g5 * (15 * dxyz[X] * dxyz[Y] * dxyz[Z]) + g6 * (10 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z];

  D7[tXXXXXZZ] += g4 * (15 * dxyz[X]) + g5 * (10 * dxyz[X] * dxyz[X] * dxyz[X] + 15 * dxyz[X] * dxyz[Z] * dxyz[Z]) +
                  g6 * (10 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z];

  D7[tXXXXYYY] += g4 * (9 * dxyz[Y]) + g5 * (18 * dxyz[X] * dxyz[X] * dxyz[Y] + 3 * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g6 * (6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] + 3 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D7[tXXXXYYZ] += g4 * (3 * dxyz[Z]) + g5 * (6 * dxyz[X] * dxyz[X] * dxyz[Z] + 3 * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g6 * (6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D7[tXXXXYZZ] += g4 * (3 * dxyz[Y]) + g5 * (6 * dxyz[X] * dxyz[X] * dxyz[Y] + 3 * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g6 * (6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D7[tXXXXZZZ] += g4 * (9 * dxyz[Z]) + g5 * (18 * dxyz[X] * dxyz[X] * dxyz[Z] + 3 * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g6 * (6 * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 3 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXXXYYYY] += g4 * (9 * dxyz[X]) + g5 * (3 * dxyz[X] * dxyz[X] * dxyz[X] + 18 * dxyz[X] * dxyz[Y] * dxyz[Y]) +
                  g6 * (3 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] + 6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D7[tXXXYYYZ] += g5 * (9 * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                  g6 * (3 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] + 3 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D7[tXXXYYZZ] += g4 * (3 * dxyz[X]) +
                  g5 * (dxyz[X] * dxyz[X] * dxyz[X] + 3 * dxyz[X] * dxyz[Y] * dxyz[Y] + 3 * dxyz[X] * dxyz[Z] * dxyz[Z]) +
                  g6 * (3 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] +
                        dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D7[tXXXYZZZ] += g5 * (9 * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                  g6 * (3 * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 3 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXXXZZZZ] += g4 * (9 * dxyz[X]) + g5 * (3 * dxyz[X] * dxyz[X] * dxyz[X] + 18 * dxyz[X] * dxyz[Z] * dxyz[Z]) +
                  g6 * (3 * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 6 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXXYYYYY] += g4 * (15 * dxyz[Y]) + g5 * (15 * dxyz[X] * dxyz[X] * dxyz[Y] + 10 * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g6 * (dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] + 10 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D7[tXXYYYYZ] += g4 * (3 * dxyz[Z]) + g5 * (3 * dxyz[X] * dxyz[X] * dxyz[Z] + 6 * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g6 * (dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] + 6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D7[tXXYYYZZ] += g4 * (3 * dxyz[Y]) +
                  g5 * (3 * dxyz[X] * dxyz[X] * dxyz[Y] + dxyz[Y] * dxyz[Y] * dxyz[Y] + 3 * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g6 * (dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] + 3 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] +
                        dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D7[tXXYYZZZ] += g4 * (3 * dxyz[Z]) +
                  g5 * (3 * dxyz[X] * dxyz[X] * dxyz[Z] + 3 * dxyz[Y] * dxyz[Y] * dxyz[Z] + dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g6 * (dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] +
                        3 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXXYZZZZ] += g4 * (3 * dxyz[Y]) + g5 * (3 * dxyz[X] * dxyz[X] * dxyz[Y] + 6 * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g6 * (dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 6 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXXZZZZZ] += g4 * (15 * dxyz[Z]) + g5 * (15 * dxyz[X] * dxyz[X] * dxyz[Z] + 10 * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g6 * (dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 10 * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXYYYYYY] += g4 * (15 * dxyz[X]) + g5 * (45 * dxyz[X] * dxyz[Y] * dxyz[Y]) +
                  g6 * (15 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D7[tXYYYYYZ] += g5 * (15 * dxyz[X] * dxyz[Y] * dxyz[Z]) + g6 * (10 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D7[tXYYYYZZ] += g4 * (3 * dxyz[X]) + g5 * (6 * dxyz[X] * dxyz[Y] * dxyz[Y] + 3 * dxyz[X] * dxyz[Z] * dxyz[Z]) +
                  g6 * (6 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] + dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g7 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D7[tXYYYZZZ] += g5 * (9 * dxyz[X] * dxyz[Y] * dxyz[Z]) +
                  g6 * (3 * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 3 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXYYZZZZ] += g4 * (3 * dxyz[X]) + g5 * (3 * dxyz[X] * dxyz[Y] * dxyz[Y] + 6 * dxyz[X] * dxyz[Z] * dxyz[Z]) +
                  g6 * (dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 6 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXYZZZZZ] += g5 * (15 * dxyz[X] * dxyz[Y] * dxyz[Z]) + g6 * (10 * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tXZZZZZZ] += g4 * (15 * dxyz[X]) + g5 * (45 * dxyz[X] * dxyz[Z] * dxyz[Z]) +
                  g6 * (15 * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[X] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tYYYYYYY] += g4 * (105 * dxyz[Y]) + g5 * (105 * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g6 * (21 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g7 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y];

  D7[tYYYYYYZ] += g4 * (15 * dxyz[Z]) + g5 * (45 * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g6 * (15 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z];

  D7[tYYYYYZZ] += g4 * (15 * dxyz[Y]) + g5 * (10 * dxyz[Y] * dxyz[Y] * dxyz[Y] + 15 * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g6 * (10 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] + dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y]) +
                  g7 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z];

  D7[tYYYYZZZ] += g4 * (9 * dxyz[Z]) + g5 * (18 * dxyz[Y] * dxyz[Y] * dxyz[Z] + 3 * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g6 * (6 * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 3 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z]) +
                  g7 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tYYYZZZZ] += g4 * (9 * dxyz[Y]) + g5 * (3 * dxyz[Y] * dxyz[Y] * dxyz[Y] + 18 * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g6 * (3 * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 6 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[Y] * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tYYZZZZZ] += g4 * (15 * dxyz[Z]) + g5 * (15 * dxyz[Y] * dxyz[Y] * dxyz[Z] + 10 * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g6 * (dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] + 10 * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[Y] * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tYZZZZZZ] += g4 * (15 * dxyz[Y]) + g5 * (45 * dxyz[Y] * dxyz[Z] * dxyz[Z]) +
                  g6 * (15 * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[Y] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

  D7[tZZZZZZZ] += g4 * (105 * dxyz[Z]) + g5 * (105 * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g6 * (21 * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z]) +
                  g7 * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z] * dxyz[Z];

#undef X
#undef Y
#undef Z
}

#endif
