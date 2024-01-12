/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  ngenic.h
 *
 *  \brief definition of a class for the construction of cosmological initial conditions
 */

#ifndef NGENIC_H
#define NGENIC_H

#include "gadgetconfig.h"

#ifdef NGENIC

#ifndef PERIODIC
#error NGENIC requires PERIODIC
#endif

#include <fftw3.h>

#ifdef DOUBLEPRECISION_FFTW
typedef double fft_real;
typedef fftw_complex fft_complex;
#else
typedef float fft_real;
typedef fftwf_complex fft_complex;
#endif

#include "../data/simparticles.h"
#include "../pm/pm_mpi_fft.h"

class ngenic : public pm_mpi_fft
{
 private:
  simparticles *Sp;

 public:
  ngenic(MPI_Comm comm, simparticles *Sp_ptr) : setcomm(comm), pm_mpi_fft(comm) /* constructor */ { Sp = Sp_ptr; }

 public:
  void ngenic_displace_particles(void);

  void create_grid(void);

 private:
  double ngenic_power_spec(double k);
  double ngenic_f1_omega(double a);
  double ngenic_f2_omega(double a);
  double ngenic_growth_factor(double astart, double aend);
  void ngenic_initialize_powerspectrum(void);
  void free_power_table(void);

  double Dplus;

  unsigned int *seedtable;

  fft_plan myplan;
  size_t maxfftsize;

  struct partbuf
  {
    MyIntPosType IntPos[3];
  };
  partbuf *partin, *partout;

  size_t nimport, nexport;

  size_t *Sndpm_count, *Sndpm_offset;
  size_t *Rcvpm_count, *Rcvpm_offset;

  gsl_rng *rnd_generator;
  gsl_rng *rnd_generator_conjugate;

  struct disp_data
  {
    fft_real deltapos[3];
  };

  disp_data *Pdisp;

  void ngenic_distribute_particles();
  void ngenic_setup_modes_in_kspace(fft_complex *fft_of_grid);
  void ngenic_readout_disp(fft_real *grid, int axis, double pfac, double vfac);
  void ngenic_initialize_ffts(void);
  void ngenic_get_derivate_from_fourier_field(int axes1, int axes2, fft_complex *fft_of_grid);
  void ngenic_compute_transform_of_source_potential(fft_real *pot);
  void print_spec(void);

  double R8;

  double AA, BB, CC;
  double nu;
  double Norm;

  int NPowerTable;

  struct pow_table
  {
    double logk, logD;
    bool operator<(const pow_table &other) const { return logk < other.logk; }
  };
  pow_table *PowerTable;

  double ngenic_powerspec_tabulated(double k);
  double ngenic_powerspec_efstathiou(double k);
  double ngenic_powerspec_eh(double k);
  double ngenic_tophat_sigma2(double R);
  double ngenic_tk_eh(double k);
  void read_power_table(void);

  double fnl(double x, double A, double B, double alpha, double beta, double V, double gf) /* Peacock & Dodds formula */
  {
    return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) / (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)),
                   1 / beta);
  }

  struct myparams
  {
    double R;
    ngenic *obj;
  };

  static double sigma2_int(double lnk, void *param)
  {
    myparams *par   = (myparams *)param;
    double r_tophat = par->R;

    double k   = exp(lnk);
    double kr  = r_tophat * k;
    double kr2 = kr * kr;
    double kr3 = kr2 * kr;

    if(kr < 1e-8)
      return 0;

    double w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
    double x = 4 * M_PI * k * k * w * w * par->obj->ngenic_power_spec(k);

    return k * x;
  }
};

#endif
#endif
