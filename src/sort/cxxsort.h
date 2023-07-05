/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  cxxsort.h
 *  \brief various sort routines
 */

#ifndef GADGET4_CXXSORT_H
#define GADGET4_CXXSORT_H

#include "gadgetconfig.h"

#include <algorithm>

#include "../data/allvars.h"
#include "../data/mymalloc.h"
#include "../logs/logs.h"

template <typename T, typename Tcomp>
void mycxxsort_internal_serial(T *begin, T *end, T *buf, bool res_into_buf, Tcomp comp)
{
  std::size_t n = end - begin;
  if(n <= 1)
    {
      if((n == 1) && res_into_buf)
        *buf = *begin;
      return;
    }

  mycxxsort_internal_serial(begin, begin + n / 2, buf, !res_into_buf, comp);
  mycxxsort_internal_serial(begin + n / 2, end, buf + n / 2, !res_into_buf, comp);

  res_into_buf ? std::merge(begin, begin + n / 2, begin + n / 2, begin + n, buf, comp)
               : std::merge(buf, buf + n / 2, buf + n / 2, buf + n, begin, comp);
}

template <typename T, typename Tcomp>
double mycxxsort(T *begin, T *end, Tcomp comp)
{
  if(end - begin <= 1)
    return 0.;

  double t0 = Logs.second();

  T *buf = (T *)Mem.mymalloc("buf", (end - begin) * sizeof(T));

  mycxxsort_internal_serial(begin, end, buf, false, comp);

  Mem.myfree(buf);

  return Logs.timediff(t0, Logs.second());
}

#endif
