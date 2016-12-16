#include "ClikPar.h"

#include <math.h>

ClikPar::ClikPar() : m_free_param_amt(21) {
  // Initialise defaults
  for (int i = 0; i < TOTAL_PARAMS; ++i) {
    m_is_gaussian[i] = false;
    m_min[i] = -9999;
    m_max[i] = -9999;
    m_mean[i] = -9999;
    m_stddev[i] = -9999;
  }

  // Set flat priors. All parameters must have a flat prior
  // LCDM parameters first
  // Comments specify up to what +/- sigma the range covers
  m_min[omega_b] = 0.021;           m_max[omega_b] = 0.023;
  // 3-sigma
  m_min[omega_cdm] = 0.108;         m_max[omega_cdm] = 0.130;
  // 4-sigma
  m_min[hundredxtheta_s] = 1.0394;  m_max[hundredxtheta_s] = 1.0423;
  // 3-sigma
  m_min[tau_reio] = 0.04;           m_max[tau_reio] = 0.116;
  // 2-sigma
  m_min[n_s] = 0.94;                m_max[n_s] = 0.99;
  // 4-sigma
  m_min[ln10_10_A_s] = 2.981;       m_max[ln10_10_A_s] = 3.197;
  // 3-sigma

  // PLC nuisance parameters second
  m_min[A_cib_217] = 0.0;         m_max[A_cib_217] = 200.0;
  m_min[cib_index] = -1.3;        m_max[cib_index] = -1.3;
  m_min[xi_sz_cib] = 0.0;         m_max[xi_sz_cib] = 1.0;
  m_min[A_sz] = 0.0;              m_max[A_sz] = 10.0;
  m_min[ps_A_100_100] = 0.0;      m_max[ps_A_100_100] = 400.0;
  m_min[ps_A_143_143] = 0.0;      m_max[ps_A_143_143] = 400.0;
  m_min[ps_A_143_217] = 0.0;      m_max[ps_A_143_217] = 400.0;
  m_min[ps_A_217_217] = 0.0;      m_max[ps_A_217_217] = 400.0;
  m_min[ksz_norm] = 0.0;          m_max[ksz_norm] = 10.0;
  m_min[gal545_A_100] = 0.0;      m_max[gal545_A_100] = 50.0;
  m_min[gal545_A_143] = 0.0;      m_max[gal545_A_143] = 50.0;
  m_min[gal545_A_143_217] = 0.0;  m_max[gal545_A_143_217] = 100.0;
  m_min[gal545_A_217] = 0.0;      m_max[gal545_A_217] = 400.0;
  m_min[calib_100T] = 0.0;        m_max[calib_100T] = 3.0;
  m_min[calib_217T] = 0.0;        m_max[calib_217T] = 3.0;
  m_min[A_planck] = 0.9;          m_max[A_planck] = 1.1;

  // Check if derived parameter min and max values are the same
  // If not, throw an error. They should be the same!
  for (int i = m_free_param_amt; i < TOTAL_PARAMS; ++i) {
    if (m_min[i] != m_max[i]) {
      /* throw error */
    }
  }

  // Set Gaussian parameters
  // m_is_gaussian[tau_reio] = true;
  m_is_gaussian[gal545_A_100] = true;
  m_is_gaussian[gal545_A_143] = true;
  m_is_gaussian[gal545_A_143_217] = true;
  m_is_gaussian[gal545_A_217] = true;
  m_is_gaussian[calib_100T] = true;
  m_is_gaussian[calib_217T] = true;
  m_is_gaussian[A_planck] = true;

  for (int i = 0; i < TOTAL_PARAMS; ++i) {
    m_gaussian_param_amt += m_is_gaussian[i];
  }

  // Set Gaussian priors
  m_mean[tau_reio] = 0.07;
  m_stddev[tau_reio] = 0.02;
  m_mean[gal545_A_100] = 7.0;
  m_stddev[gal545_A_100] = 2.0;
  m_mean[gal545_A_143] = 9.0;
  m_stddev[gal545_A_143] = 2.0;
  m_mean[gal545_A_143_217] = 21.0;
  m_stddev[gal545_A_143_217] = 8.5;
  m_mean[gal545_A_217] = 80.0;
  m_stddev[gal545_A_217] = 20.0;
  m_mean[calib_100T] = 0.9990004;
  m_stddev[calib_100T] = 0.001;
  m_mean[calib_217T] = 0.99501;
  m_stddev[calib_217T] = 0.002;
  m_mean[A_planck] = 1.0;
  m_stddev[A_planck] = 0.0025;
}

ClikPar::~ClikPar() {
  /* just destructor things */
}

int ClikPar::get_free_param_amt() const {
  return m_free_param_amt;
}

int ClikPar::get_gaussian_param_amt() const {
  return m_gaussian_param_amt;
}

void ClikPar::scale_Cube(double* Cube) {
  for (int i = 0; i < TOTAL_PARAMS; i++) {
    if (i < m_free_param_amt) { // free parameter
      Cube[i] = m_min[i] + (m_max[i] - m_min[i]) * Cube[i];
    }
    else { // derived parameter
      Cube[i] = m_min[i];
    }
  }
}

double ClikPar::calculate_gaussian_priors(double* Cube) const {
  double loglike = 0.0;
  for (int i = 0; i < TOTAL_PARAMS; ++i) {
    // Only use Gaussian priors for variables that are Gaussian
    if (m_is_gaussian[i]) {
      loglike += pow(Cube[i] - m_mean[i], 2.0) / (2.0 * pow(m_stddev[i], 2.0));
    }
  }

  return loglike;
}

double ClikPar::calculate_misc_priors(double* Cube) const {
  double loglike = 0.0;

  // Calculate SZ degeneracy prior
  double SZ_val = Cube[ksz_norm] + 1.6 * Cube[A_sz];
  double SZ_mean = 9.4;
  double SZ_stddev = 1.4;

  loglike += pow(SZ_val - SZ_mean, 2.0) / (2.0 * pow(SZ_stddev, 2.0));

  return loglike;
}