#include "ClikPar.h"

#include <math.h> // for pow()
#include <iostream> // for std::cerr
#include <exception> // for std::exception

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
  // Using full Planck priors
  m_min[omega_b] = 0.005;           m_max[omega_b] = 0.1;
  m_min[omega_cdm] = 0.001;         m_max[omega_cdm] = 0.99;
  m_min[hundredxtheta_s] = 0.5;  m_max[hundredxtheta_s] = 10.0;
  m_min[tau_reio] = 0.01;           m_max[tau_reio] = 0.8;
  m_min[n_s] = 0.8;                m_max[n_s] = 1.2;
  m_min[ln10_10_A_s] = 2.0;       m_max[ln10_10_A_s] = 4.0;

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
  for (int param = m_free_param_amt; param < TOTAL_PARAMS; ++param) {
    if (m_min[param] != m_max[param]) {
      std::cerr << "[ERROR]: min and max values are different "
                << "for parameter #"
                << param
                << ", a derived parameter. Consider changing "
                << "max to min value. Unexpected behaviour "
                << "will occur."
                << std::endl;
      throw std::exception();
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
  for (int param = 0; param < TOTAL_PARAMS; param++) {
    if (param < m_free_param_amt) { // free parameter
      Cube[param] = m_min[param] + (m_max[param] - m_min[param]) * Cube[param];
    }
    else { // derived parameter
      Cube[param] = m_min[param];
    }
  }
}

/* Private Methods */

// *** The returned loglike value must always be positive ***
double ClikPar::calculate_extra_priors(double* Cube, ClassEngine* class_engine) const {
  double loglike = 0.0;

  // Calculate Gaussian priors
  for (int param = 0; param < TOTAL_PARAMS; ++param) {
    // Only use Gaussian priors for variables that are Gaussian
    if (m_is_gaussian[param]) {
      loglike += pow(Cube[param] - m_mean[param], 2.0) / (2.0 * pow(m_stddev[param], 2.0));
    }
  }

  // Calculate SZ degeneracy prior
  double SZ_val = Cube[ksz_norm] + 1.6 * Cube[A_sz];
  double SZ_mean = 9.5;
  double SZ_stddev = 3.0;

  loglike += pow(SZ_val - SZ_mean, 2.0) / (2.0 * pow(SZ_stddev, 2.0));

  // Calculate flat [20,100] prior on H0
  double H0 = class_engine->get_H0();
  double H0_min = 20.0, H0_max = 100.0;

  double delta = 0.0;
  if (H0 < H0_min or H0 > H0_max) {
    delta = 1E90;
  }

  loglike += delta;

  // Calculate any other priors needed ...

  return loglike;
}