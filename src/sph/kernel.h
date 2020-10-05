/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  kernel.h
 *
 *  \brief collects definitions of different SPH kernels
 */

#ifndef KERNEL_H
#define KERNEL_H

struct kernel_density
{
  double dpos[NUMDIMS];
  double r;
  double dv[NUMDIMS];
  double wk, dwk;
  double hinv, hinv3, hinv4;
  double mj_wk, mj_dwk_r;
};

struct kernel_hydra
{
  double dx, dy, dz;
  double r, vsig, sound_i, sound_j;
  double dvx, dvy, dvz, vdotr2;
  double wk_i, wk_j, dwk_i, dwk_j;
  double h_i, h_j, dwk_ij, rho_ij_inv;
};

#if !defined(CUBIC_SPLINE_KERNEL) && !defined(WENDLAND_C2_KERNEL) && !defined(WENDLAND_C4_KERNEL) && !defined(WENDLAND_C6_KERNEL)
#define CUBIC_SPLINE_KERNEL /* fall back to cubic spline kernel */
#endif

/* fall back to three dimensions */
#if !defined(TWODIMS) && !defined(ONEDIMS)
#define THREEDIMS
#endif

/* Norms */
#ifdef CUBIC_SPLINE_KERNEL

#ifdef THREEDIMS
#define NORM (8.0 / M_PI) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define NORM (40.0 / (7.0 * M_PI)) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define NORM (4.0 / 3.0) /*!< For 1D-normalized kernel */
#endif

#endif /* CUBIC_SPLINE_KERNEL */

#ifdef WENDLAND_C2_KERNEL

#ifdef THREEDIMS
#define NORM (21.0 / (2.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define NORM (7.0 / M_PI) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define NORM (5.0 / 4.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C2_KERNEL */

#ifdef WENDLAND_C4_KERNEL

#ifdef THREEDIMS
#define NORM (495.0 / (32.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define NORM (9.0 / M_PI) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define NORM (3.0 / 2.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL

#ifdef THREEDIMS
#define NORM (1365.0 / (64.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define NORM (78.0 / (7.0 * M_PI)) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define NORM (55.0 / 32.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C6_KERNEL */

#define COMPUTE_WK -1
#define COMPUTE_WK_AND_DWK 0
#define COMPUTE_DWK 1

static inline void kernel_hinv(double h, double *hinv, double *hinv3, double *hinv4)
{
  *hinv = 1.0 / h;

#ifdef THREEDIMS
  *hinv3 = *hinv * *hinv * *hinv;
#endif

#ifdef TWODIMS
  *hinv3 = *hinv * *hinv;
#endif

#ifdef ONEDIMS
  *hinv3 = *hinv;
#endif

  *hinv4 = *hinv3 * *hinv;
}

/* Attention: Here we assume that kernel is only called
   with range 0..1 for u as done in hydra or density !!
   Call with mode COMPUTE_WK_AND_DWK to calculate dwk and wk
   Call with mode COMPUTE_WK to calculate only wk
   Call with mode COMPUTE_DWK to calculate only dwk */
static inline void kernel_main(double u, double hinv3, double hinv4, double *wk, double *dwk, int mode)
{
#ifdef CUBIC_SPLINE_KERNEL
#if defined(WENDLAND_C2_KERNEL) || defined(WENDLAND_C4_KERNEL) || defined(WENDLAND_C6_KERNEL)
#error "Only one SPH kernel can be used"
#endif
  if(u < 0.5)
    {
      if(mode >= COMPUTE_WK_AND_DWK)
        *dwk = u * (18.0 * u - 12.0);
      if(mode <= COMPUTE_WK_AND_DWK)
        *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
  else
    {
      double t1 = (1.0 - u);
      double t2 = t1 * t1;
      if(mode >= COMPUTE_WK_AND_DWK)
        *dwk = -6.0 * t2;
      if(mode <= COMPUTE_WK_AND_DWK)
        *wk = 2.0 * t2 * t1;
    }
#endif

#ifdef WENDLAND_C2_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);

  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -12.0 * u * t2;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t2 * t1 * (1.0 + u * 3.0);

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -20.0 * u * t2 * t1;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t4 * (1.0 + u * 4.0);

#endif
#endif /* WENDLAND_C2_KERNEL */

#ifdef WENDLAND_C4_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = t1 * t1;
  double t4 = t2 * t2;
  double t5 = t4 * t1;

  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -14.0 * t4 * (4.0 * u + 1) * u;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t5 * (1.0 + u * (5.0 + 8.0 * u));

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t6 = t2 * t2 * t2;
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -56.0 / 3.0 * u * t4 * t1 * (5.0 * u + 1);
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t6 * (1.0 + u * (6.0 + 35.0 / 3.0 * u));

#endif
#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t6 = t4 * t2;
  double t7 = t4 * t2 * t1;
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -6.0 * u * t6 * (3.0 + u * (18.0 + 35.0 * u));
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t7 * (1.0 + u * (7.0 + u * (19.0 + 21.0 * u)));

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t7 = t4 * t2 * t1;
  double t8 = t4 * t4;
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk = -22.0 * u * (1.0 + u * (7.0 + 16.0 * u)) * t7;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk = t8 * (1.0 + u * (8.0 + u * (25.0 + 32.0 * u)));

#endif
#endif /* WENDLAND_C6_KERNEL */
  if(mode >= COMPUTE_WK_AND_DWK)
    *dwk *= NORM * hinv4;
  if(mode <= COMPUTE_WK_AND_DWK)
    *wk *= NORM * hinv3;
}

#if defined(WENDLAND_BIAS_CORRECTION) && (!(defined(WENDLAND_C2_KERNEL) || defined(WENDLAND_C4_KERNEL) || defined(WENDLAND_C6_KERNEL)))
#error "WENDLAND_BIAS_CORRECTION only works with a Wendland kernel"
#endif

#if defined(WENDLAND_BIAS_CORRECTION) && (defined(WENDLAND_C2_KERNEL) || defined(WENDLAND_C4_KERNEL) || defined(WENDLAND_C6_KERNEL))

#if defined(ONEDIMS) || defined(TWODIMS)
#error "WENDLAND_BIAS_CORRECTION is only implemented for 3D"
#endif

static inline void get_bias_correction_parameters(double *alpha, double *eps100)
{
#ifdef WENDLAND_C2_KERNEL
  *eps100 = 0.0294;
  *alpha  = 0.977;
#endif
#ifdef WENDLAND_C4_KERNEL
  *eps100 = 0.01342;
  *alpha  = 1.579;
#endif
#ifdef WENDLAND_C6_KERNEL
  *eps100 = 0.0116;
  *alpha  = 2.236;
#endif
}
static inline double get_density_bias(double hsml, double mass, int DesNumNgb)
{
  kernel_density kernel;
  kernel_hinv(hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
  kernel_main(0, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, COMPUTE_WK);
  double alpha  = 0;
  double eps100 = 0;
  get_bias_correction_parameters(&alpha, &eps100);
  double wc_correction = eps100 * pow(DesNumNgb * 0.01, -alpha) * mass * kernel.wk;
  return wc_correction;
}
#endif

#ifdef EXPLICIT_VECTORIZATION
static inline void kernel_main_vector(Vec4d u, Vec4d hinv3, Vec4d hinv4, Vec4d *wk, Vec4d *dwk)
{
#ifdef CUBIC_SPLINE_KERNEL
  Vec4d ucompl    = u - 1.0;
  Vec4db decision = (u < 0.5);
  Vec4d wksub     = 6.0 * u;

  Vec4d ucompsq = ucompl * ucompl;
  Vec4d uucomp  = u * ucompl;
  Vec4d wk1     = 1.0 + uucomp * wksub;
  Vec4d dwk1    = wksub + uucomp * 18.0;
  Vec4d wk2     = ucompsq * -2.0 * ucompl;
  Vec4d dwk2    = ucompsq * -6.0;
  *wk           = select(decision, wk1, wk2);
  *dwk          = select(decision, dwk1, dwk2);
#endif

#ifdef WENDLAND_C2_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  Vec4d t1 = 1.0 - u;

  Vec4d t2 = t1 * t1;

  *dwk = -12.0 * u * t2;
  *wk  = t2 * t1 * (1.0 + u * 3.0);

#else /* 2d or 3d */
  Vec4d t1 = (1.0 - u);
  Vec4d t2 = (t1 * t1);

  *dwk     = -20.0 * u * t2 * t1;
  *wk      = t2 * t2 * (1.0 + u * 4.0);
#endif
#endif /* WENDLAND_C2_KERNEL */

#ifdef WENDLAND_C4_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  Vec4d t1 = (1.0 - u);
  Vec4d t2 = t1 * t1;
  Vec4d t4 = t2 * t2;

  *dwk = -14.0 * t4 * (4.0 * u + 1) * u;
  *wk  = t1 * t4 * (1.0 + u * (5.0 + 8.0 * u));

#else /* 2d or 3d */
  Vec4d t1 = (1.0 - u);
  Vec4d t2 = (t1 * t1);
  Vec4d t4 = t2 * t2;
  Vec4d t6 = t2 * t2 * t2;
  *dwk     = -56.0 / 3.0 * u * t4 * t1 * (5.0 * u + 1);

  *wk      = t6 * (1.0 + u * (6.0 + 35.0 / 3.0 * u));

#endif
#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  Vec4d t1 = (1.0 - u);
  Vec4d t2 = (t1 * t1);
  Vec4d t4 = t2 * t2;
  Vec4d t6 = t4 * t2;
  Vec4d t7 = t4 * t2 * t1;
  *dwk     = -6.0 * u * t6 * (3.0 + u * (18.0 + 35.0 * u));
  *wk      = t7 * (1.0 + u * (7.0 + u * (19.0 + 21.0 * u)));

#else /* 2d or 3d */
  Vec4d t1 = (1.0 - u);
  Vec4d t2 = (t1 * t1);
  Vec4d t4 = t2 * t2;
  Vec4d t7 = t4 * t2 * t1;
  Vec4d t8 = t4 * t4;
  *dwk     = -22.0 * u * (1.0 + u * (7.0 + 16.0 * u)) * t7;
  *wk      = t8 * (1.0 + u * (8.0 + u * (25.0 + 32.0 * u)));

#endif
#endif /* WENDLAND_C6_KERNEL */
  *dwk *= hinv4;
  *wk *= hinv3;
}
#endif /* EXPLICIT_VECTORIZATION */
#endif
