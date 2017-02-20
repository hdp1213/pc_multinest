#include "ClikPar.h"

#include <math.h> // for pow()
#include <iostream> // for std::cerr
#include <exception> // for std::exception

ClikPar::ClikPar(int free_param_amt, int fixed_param_amt) : m_free_param_amt(free_param_amt), m_fixed_param_amt(free_param_amt + fixed_param_amt), m_class_engine(0) {
  // Initialise defaults
  for (int param = 0; param < TOTAL_PARAMS; ++param) {
    m_min[param] = -9999;
    m_max[param] = -9999;
    m_is_gaussian[param] = false;
    m_mean[param] = -9999;
    m_stddev[param] = -9999;
  }

  // Set flat priors. All parameters must have a flat prior
  // LCDM parameters first
  // Using workable priors as a mix between v0.10 and v0.11
  // priors
  // Just to make sure on the next run we get the whole bell
  // curve
  m_min[omega_b] = 0.016;          m_max[omega_b] = 0.028;
  m_min[omega_cdm] = 0.108;        m_max[omega_cdm] = 0.130;
  m_min[hundredxtheta_s] = 1.039;  m_max[hundredxtheta_s] = 1.043;
  m_min[tau_reio] = 0.01;          m_max[tau_reio] = 0.15;
  m_min[ln10_10_A_s] = 2.98;       m_max[ln10_10_A_s] = 3.20;
  m_min[n_s] = 0.92;               m_max[n_s] = 1.04;

  // PLC nuisance parameters second
  m_min[A_planck] = 0.9;           m_max[A_planck] = 1.1;
  m_min[A_cib_217] = 0.0;          m_max[A_cib_217] = 200.0;
  m_min[xi_sz_cib] = 0.0;          m_max[xi_sz_cib] = 1.0;
  m_min[A_sz] = 0.0;               m_max[A_sz] = 10.0;
  m_min[ps_A_100_100] = 0.0;       m_max[ps_A_100_100] = 400.0;
  m_min[ps_A_143_143] = 0.0;       m_max[ps_A_143_143] = 400.0;
  m_min[ps_A_143_217] = 0.0;       m_max[ps_A_143_217] = 400.0;
  m_min[ps_A_217_217] = 0.0;       m_max[ps_A_217_217] = 400.0;
  m_min[ksz_norm] = 0.0;           m_max[ksz_norm] = 10.0;
  m_min[gal545_A_100] = 0.0;       m_max[gal545_A_100] = 50.0;
  m_min[gal545_A_143] = 0.0;       m_max[gal545_A_143] = 50.0;
  m_min[gal545_A_143_217] = 0.0;   m_max[gal545_A_143_217] = 100.0;
  m_min[gal545_A_217] = 0.0;       m_max[gal545_A_217] = 400.0;
  m_min[calib_100T] = 0.0;         m_max[calib_100T] = 2.0;
  m_min[calib_217T] = 0.0;         m_max[calib_217T] = 2.0;

  // Fixed parameters
  m_min[cib_index] = -1.3;         m_max[cib_index] = -1.3;

  // There is a difference between a "derived" parameter and
  // a parameter whose value should be fixed.

  for (int param = m_free_param_amt; param < m_fixed_param_amt; ++param) {
    if (m_min[param] != m_max[param]) {
      std::cerr << "[ERROR]: min and max values are different "
                << "for parameter #"
                << param
                << ", a fixed parameter. Consider changing "
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

  for (int param = 0; param < TOTAL_PARAMS; ++param) {
    m_gaussian_param_amt += m_is_gaussian[param];
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
  delete m_class_engine;
}

void ClikPar::initialise_CLASS(int max_l) {
  // Prepare CLASS engine
  ClassParams default_params;

  // MultiNest variables
  default_params.add("omega_b", 0.022032);
  default_params.add("omega_cdm", 0.12038);
  default_params.add("100*theta_s", 1.042143);
  default_params.add("tau_reio", 0.0925);
  default_params.add("ln10^{10}A_s", 3.0980);
  default_params.add("n_s", 0.9619);

  // Neutrino values to set
  default_params.add("N_ur", 2.0328);
  default_params.add("N_ncdm", 1);
  default_params.add("m_ncdm", 0.06); // MeV
  // default_params.add("T_ncdm", 0.71611);

  // Perturbation options
  // default_params.add("P_k_max_h/Mpc", 1.);
  // default_params.add("k_pivot", 0.05); // Mpc-1
  
  // Spectra output options
  default_params.add("output", "tCl,pCl,lCl");
  default_params.add("lensing", true);   //note boolean
  default_params.add("l_max_scalars", max_l);
  default_params.add("format", "camb");

  // Initialise CLASS engine
  try {
    m_class_engine = new ClassEngine(default_params);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS initialisation in ClikPar "
              << "failed, throwing exception "
              << e.what()
              << std::endl;
    throw e;
  }
}

void ClikPar::scale_Cube(double* Cube) {
  // Set free, fixed and derived parameter values
  for (int param = 0; param < TOTAL_PARAMS; param++) {
    if (param < m_free_param_amt) { // free parameter
      Cube[param] = m_min[param] + (m_max[param] - m_min[param]) * Cube[param];
    }
    else if (param < m_fixed_param_amt) { // fixed parameter
      Cube[param] = m_min[param];
    }
    else { // derived parameter
      // Must be set to something in case CLASS fails
      Cube[param] = m_min[param];
    }
  }
}

void ClikPar::set_derived_params(double* Cube) {
  // Set derived parameter values in Cube
  // Non-standard CLASS routines
  Cube[H0] = m_class_engine->get_H0();
  Cube[Omega_b] = m_class_engine->get_Omega_b();
  Cube[Omega_cdm] = m_class_engine->get_Omega_cdm();
  Cube[Omega_L] = m_class_engine->get_Omega_L();
  Cube[Omega_g] = m_class_engine->get_Omega_g();
  // Cube[Omega_k] = m_class_engine->get_Omega_k();
  Cube[sigma8] = m_class_engine->get_sigma8();
  Cube[age] = m_class_engine->get_age();
  Cube[conf_age] = m_class_engine->get_conf_age();
  // Cube[K] = m_class_engine->get_K();

  // Standard CLASS routines
  Cube[z_drag] = m_class_engine->z_drag();
  Cube[rs_drag] = m_class_engine->rs_drag();
}

ClassEngine* ClikPar::get_CLASS() {
  return m_class_engine;
}

/* Private Methods */

// *** The returned loglike value must always be positive ***
double ClikPar::calculate_extra_priors(double* Cube) const {
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
  double H0 = m_class_engine->get_H0();
  double H0_min = 20.0, H0_max = 100.0;

  double delta = 0.0;
  if (H0 < H0_min or H0 > H0_max) {
    delta = 1E90;
  }

  loglike += delta;

  // Calculate any other priors needed ...

  return loglike;
}
