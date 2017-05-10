#include "ClikPar.h"

#include <math.h> // for pow()
#include <iostream> // for std::cerr
#include <exception> // for std::exception

// Wrapped methods for scale_Cube()
inline double pow10(double x) { return pow(10., x); }
inline double self(double x) { return x; }

ClikPar::ClikPar() : m_class_engine(0) {
  // Initialise defaults
  for (int param = 0; param < TOTAL_PARAMS; ++param) {
    m_min[param] = -9999;
    m_max[param] = -9999;
    m_has_gaussian_prior[param] = false;
    m_is_log10[param] = false;
    m_mean[param] = -9999;
    m_stddev[param] = -9999;
  }

  /***************************/
  /*  UNIFORM PRIOR SETTING  */
  /***************************/

  // Fix LCDM parameters to TTTEEE+lowP best fit values (2015)
  m_min[omega_b] = 0.022252;          m_max[omega_b] = 0.022252;
  m_min[omega_cdm] = 0.11987;         m_max[omega_cdm] = 0.11987;
  m_min[hundredxtheta_s] = 1.040778;  m_max[hundredxtheta_s] = 1.040778;
  m_min[tau_reio] = 0.0789;           m_max[tau_reio] = 0.0789;
  m_min[ln10_10_A_s] = 3.0929;        m_max[ln10_10_A_s] = 3.0929;
  m_min[n_s] = 0.96475;               m_max[n_s] = 0.96475;

  // Non-standard LCDM parameters
  m_is_log10[pbh_frac] = true;
  m_min[pbh_frac] = -20.0;          m_max[pbh_frac] = 0.0;

  // Nuisance parameters (PLC TT & TTTEEE)
  // Using thinner bounds from base_plikTTTEEE_lowTEB results
  m_min[A_planck] = 0.9;           m_max[A_planck] = 1.1;
#ifndef LITE_HI_L
  m_min[A_cib_217] = 40.0;         m_max[A_cib_217] = 90.0;
  m_min[xi_sz_cib] = 0.0;          m_max[xi_sz_cib] = 1.0;
  m_min[A_sz] = 0.0;               m_max[A_sz] = 10.0;
  m_min[ps_A_100_100] = 150.0;     m_max[ps_A_100_100] = 370.0;
  m_min[ps_A_143_143] = 15.0;      m_max[ps_A_143_143] = 75.0;
  m_min[ps_A_143_217] = 0.0;       m_max[ps_A_143_217] = 80.0;
  m_min[ps_A_217_217] = 50.0;      m_max[ps_A_217_217] = 140.0;
  m_min[ksz_norm] = 0.0;           m_max[ksz_norm] = 10.0;
  m_min[gal545_A_100] = 0.0;       m_max[gal545_A_100] = 20.0;
  m_min[gal545_A_143] = 0.0;       m_max[gal545_A_143] = 20.0;
  m_min[gal545_A_143_217] = 0.0;   m_max[gal545_A_143_217] = 40.0;
  m_min[gal545_A_217] = 50.0;       m_max[gal545_A_217] = 110.0;
  m_min[calib_100T] = 0.0;         m_max[calib_100T] = 2.0;
  m_min[calib_217T] = 0.0;         m_max[calib_217T] = 2.0;

  // Nuisance parameters (PLC TTTEEE)
  // Using approximate bounds on base_plikTTTEEE_lowTEB results
  m_min[galf_EE_A_100] = 0.04;
  m_max[galf_EE_A_100] = 0.12;
  m_min[galf_EE_A_100_143] = 0.01;
  m_max[galf_EE_A_100_143] = 0.09;
  m_min[galf_EE_A_100_217] = 0.00;
  m_max[galf_EE_A_100_217] = 0.26;
  m_min[galf_EE_A_143] = 0.060;
  m_max[galf_EE_A_143] = 0.150;
  m_min[galf_EE_A_143_217] = 0.0;
  m_max[galf_EE_A_143_217] = 0.6;
  m_min[galf_EE_A_217] = 0.00;
  m_max[galf_EE_A_217] = 1.50;
  m_min[galf_TE_A_100] = 0.00;
  m_max[galf_TE_A_100] = 0.40;
  m_min[galf_TE_A_100_143] = 0.00;
  m_max[galf_TE_A_100_143] = 0.30;
  m_min[galf_TE_A_100_217] = 0.0;
  m_max[galf_TE_A_100_217] = 0.8;
  m_min[galf_TE_A_143] = 0.0;
  m_max[galf_TE_A_143] = 0.5;
  m_min[galf_TE_A_143_217] = 0.00;
  m_max[galf_TE_A_143_217] = 0.75;
  m_min[galf_TE_A_217] = 0.0;
  m_max[galf_TE_A_217] = 4.5;

  // Fixed parameters (PLC TT & TTTEEE)
  m_min[cib_index] = -1.3;
  m_max[cib_index] = -1.3;
  // Fixed parameters (PLC TTTEEE)
  m_min[galf_EE_index] = -2.4;
  m_max[galf_EE_index] = -2.4;
  m_min[galf_TE_index] = -2.4;
  m_max[galf_TE_index] = -2.4;
  m_min[bleak_epsilon_0_0T_0E] = 0.0;
  m_max[bleak_epsilon_0_0T_0E] = 0.0;
  m_min[bleak_epsilon_1_0T_0E] = 0.0;
  m_max[bleak_epsilon_1_0T_0E] = 0.0;
  m_min[bleak_epsilon_2_0T_0E] = 0.0;
  m_max[bleak_epsilon_2_0T_0E] = 0.0;
  m_min[bleak_epsilon_3_0T_0E] = 0.0;
  m_max[bleak_epsilon_3_0T_0E] = 0.0;
  m_min[bleak_epsilon_4_0T_0E] = 0.0;
  m_max[bleak_epsilon_4_0T_0E] = 0.0;
  m_min[bleak_epsilon_0_0T_1E] = 0.0;
  m_max[bleak_epsilon_0_0T_1E] = 0.0;
  m_min[bleak_epsilon_1_0T_1E] = 0.0;
  m_max[bleak_epsilon_1_0T_1E] = 0.0;
  m_min[bleak_epsilon_2_0T_1E] = 0.0;
  m_max[bleak_epsilon_2_0T_1E] = 0.0;
  m_min[bleak_epsilon_3_0T_1E] = 0.0;
  m_max[bleak_epsilon_3_0T_1E] = 0.0;
  m_min[bleak_epsilon_4_0T_1E] = 0.0;
  m_max[bleak_epsilon_4_0T_1E] = 0.0;
  m_min[bleak_epsilon_0_0T_2E] = 0.0;
  m_max[bleak_epsilon_0_0T_2E] = 0.0;
  m_min[bleak_epsilon_1_0T_2E] = 0.0;
  m_max[bleak_epsilon_1_0T_2E] = 0.0;
  m_min[bleak_epsilon_2_0T_2E] = 0.0;
  m_max[bleak_epsilon_2_0T_2E] = 0.0;
  m_min[bleak_epsilon_3_0T_2E] = 0.0;
  m_max[bleak_epsilon_3_0T_2E] = 0.0;
  m_min[bleak_epsilon_4_0T_2E] = 0.0;
  m_max[bleak_epsilon_4_0T_2E] = 0.0;
  m_min[bleak_epsilon_0_1T_1E] = 0.0;
  m_max[bleak_epsilon_0_1T_1E] = 0.0;
  m_min[bleak_epsilon_1_1T_1E] = 0.0;
  m_max[bleak_epsilon_1_1T_1E] = 0.0;
  m_min[bleak_epsilon_2_1T_1E] = 0.0;
  m_max[bleak_epsilon_2_1T_1E] = 0.0;
  m_min[bleak_epsilon_3_1T_1E] = 0.0;
  m_max[bleak_epsilon_3_1T_1E] = 0.0;
  m_min[bleak_epsilon_4_1T_1E] = 0.0;
  m_max[bleak_epsilon_4_1T_1E] = 0.0;
  m_min[bleak_epsilon_0_1T_2E] = 0.0;
  m_max[bleak_epsilon_0_1T_2E] = 0.0;
  m_min[bleak_epsilon_1_1T_2E] = 0.0;
  m_max[bleak_epsilon_1_1T_2E] = 0.0;
  m_min[bleak_epsilon_2_1T_2E] = 0.0;
  m_max[bleak_epsilon_2_1T_2E] = 0.0;
  m_min[bleak_epsilon_3_1T_2E] = 0.0;
  m_max[bleak_epsilon_3_1T_2E] = 0.0;
  m_min[bleak_epsilon_4_1T_2E] = 0.0;
  m_max[bleak_epsilon_4_1T_2E] = 0.0;
  m_min[bleak_epsilon_0_2T_2E] = 0.0;
  m_max[bleak_epsilon_0_2T_2E] = 0.0;
  m_min[bleak_epsilon_1_2T_2E] = 0.0;
  m_max[bleak_epsilon_1_2T_2E] = 0.0;
  m_min[bleak_epsilon_2_2T_2E] = 0.0;
  m_max[bleak_epsilon_2_2T_2E] = 0.0;
  m_min[bleak_epsilon_3_2T_2E] = 0.0;
  m_max[bleak_epsilon_3_2T_2E] = 0.0;
  m_min[bleak_epsilon_4_2T_2E] = 0.0;
  m_max[bleak_epsilon_4_2T_2E] = 0.0;
  m_min[bleak_epsilon_0_0E_0E] = 0.0;
  m_max[bleak_epsilon_0_0E_0E] = 0.0;
  m_min[bleak_epsilon_1_0E_0E] = 0.0;
  m_max[bleak_epsilon_1_0E_0E] = 0.0;
  m_min[bleak_epsilon_2_0E_0E] = 0.0;
  m_max[bleak_epsilon_2_0E_0E] = 0.0;
  m_min[bleak_epsilon_3_0E_0E] = 0.0;
  m_max[bleak_epsilon_3_0E_0E] = 0.0;
  m_min[bleak_epsilon_4_0E_0E] = 0.0;
  m_max[bleak_epsilon_4_0E_0E] = 0.0;
  m_min[bleak_epsilon_0_0E_1E] = 0.0;
  m_max[bleak_epsilon_0_0E_1E] = 0.0;
  m_min[bleak_epsilon_1_0E_1E] = 0.0;
  m_max[bleak_epsilon_1_0E_1E] = 0.0;
  m_min[bleak_epsilon_2_0E_1E] = 0.0;
  m_max[bleak_epsilon_2_0E_1E] = 0.0;
  m_min[bleak_epsilon_3_0E_1E] = 0.0;
  m_max[bleak_epsilon_3_0E_1E] = 0.0;
  m_min[bleak_epsilon_4_0E_1E] = 0.0;
  m_max[bleak_epsilon_4_0E_1E] = 0.0;
  m_min[bleak_epsilon_0_0E_2E] = 0.0;
  m_max[bleak_epsilon_0_0E_2E] = 0.0;
  m_min[bleak_epsilon_1_0E_2E] = 0.0;
  m_max[bleak_epsilon_1_0E_2E] = 0.0;
  m_min[bleak_epsilon_2_0E_2E] = 0.0;
  m_max[bleak_epsilon_2_0E_2E] = 0.0;
  m_min[bleak_epsilon_3_0E_2E] = 0.0;
  m_max[bleak_epsilon_3_0E_2E] = 0.0;
  m_min[bleak_epsilon_4_0E_2E] = 0.0;
  m_max[bleak_epsilon_4_0E_2E] = 0.0;
  m_min[bleak_epsilon_0_1E_1E] = 0.0;
  m_max[bleak_epsilon_0_1E_1E] = 0.0;
  m_min[bleak_epsilon_1_1E_1E] = 0.0;
  m_max[bleak_epsilon_1_1E_1E] = 0.0;
  m_min[bleak_epsilon_2_1E_1E] = 0.0;
  m_max[bleak_epsilon_2_1E_1E] = 0.0;
  m_min[bleak_epsilon_3_1E_1E] = 0.0;
  m_max[bleak_epsilon_3_1E_1E] = 0.0;
  m_min[bleak_epsilon_4_1E_1E] = 0.0;
  m_max[bleak_epsilon_4_1E_1E] = 0.0;
  m_min[bleak_epsilon_0_1E_2E] = 0.0;
  m_max[bleak_epsilon_0_1E_2E] = 0.0;
  m_min[bleak_epsilon_1_1E_2E] = 0.0;
  m_max[bleak_epsilon_1_1E_2E] = 0.0;
  m_min[bleak_epsilon_2_1E_2E] = 0.0;
  m_max[bleak_epsilon_2_1E_2E] = 0.0;
  m_min[bleak_epsilon_3_1E_2E] = 0.0;
  m_max[bleak_epsilon_3_1E_2E] = 0.0;
  m_min[bleak_epsilon_4_1E_2E] = 0.0;
  m_max[bleak_epsilon_4_1E_2E] = 0.0;
  m_min[bleak_epsilon_0_2E_2E] = 0.0;
  m_max[bleak_epsilon_0_2E_2E] = 0.0;
  m_min[bleak_epsilon_1_2E_2E] = 0.0;
  m_max[bleak_epsilon_1_2E_2E] = 0.0;
  m_min[bleak_epsilon_2_2E_2E] = 0.0;
  m_max[bleak_epsilon_2_2E_2E] = 0.0;
  m_min[bleak_epsilon_3_2E_2E] = 0.0;
  m_max[bleak_epsilon_3_2E_2E] = 0.0;
  m_min[bleak_epsilon_4_2E_2E] = 0.0;
  m_max[bleak_epsilon_4_2E_2E] = 0.0;
  m_min[calib_100P] = 1.0;
  m_max[calib_100P] = 1.0;
  m_min[calib_143P] = 1.0;
  m_max[calib_143P] = 1.0;
  m_min[calib_217P] = 1.0;
  m_max[calib_217P] = 1.0;
  m_min[A_pol] = 1.0;
  m_max[A_pol] = 1.0;
#endif

  // There is a difference between a "derived" parameter and
  // a parameter whose value should be fixed.

  // Do checks on min and max values for fixed parameters
  for (int param = FREE_PARAMS; param < FIXED_PARAMS; ++param) {
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

  // Check if parameter priors are in log10 space.
  // If so, set transformation function to the power of 10 such that Cube[param] is set to the correct value
  for (int param = 0; param < TOTAL_PARAMS; ++param) {
    m_transform[param] = m_is_log10[param] ? &pow10 : &self;
  }

  /****************************/
  /*  GAUSSIAN PRIOR SETTING  */
  /****************************/

#ifdef GAUSS_TAU
  m_has_gaussian_prior[tau_reio] = true;
#endif
#ifndef LITE_HI_L
  // Nuisance parameters (PLC TT & TTTEEE)
  m_has_gaussian_prior[gal545_A_100] = true;
  m_has_gaussian_prior[gal545_A_143] = true;
  m_has_gaussian_prior[gal545_A_143_217] = true;
  m_has_gaussian_prior[gal545_A_217] = true;
  m_has_gaussian_prior[calib_100T] = true;
  m_has_gaussian_prior[calib_217T] = true;

  // Nuisance parameters (PLC TTTEEE)
  m_has_gaussian_prior[galf_EE_A_100] = true;
  m_has_gaussian_prior[galf_EE_A_100_143] = true;
  m_has_gaussian_prior[galf_EE_A_100_217] = true;
  m_has_gaussian_prior[galf_EE_A_143] = true;
  m_has_gaussian_prior[galf_EE_A_143_217] = true;
  m_has_gaussian_prior[galf_EE_A_217] = true;
  m_has_gaussian_prior[galf_TE_A_100] = true;
  m_has_gaussian_prior[galf_TE_A_100_143] = true;
  m_has_gaussian_prior[galf_TE_A_100_217] = true;
  m_has_gaussian_prior[galf_TE_A_143] = true;
  m_has_gaussian_prior[galf_TE_A_143_217] = true;
  m_has_gaussian_prior[galf_TE_A_217] = true;
#endif
  m_has_gaussian_prior[A_planck] = true;

  // Set Gaussian priors
#ifdef GAUSS_TAU
  m_mean[tau_reio] = 0.07;
  m_stddev[tau_reio] = 0.02;
#endif
#ifndef LITE_HI_L
  // Nuisance parameters (PLC TT & TTTEEE)
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

  // Nuisance parameters (PLC TTTEEE)
  m_mean[galf_EE_A_100] = 0.060;
  m_stddev[galf_EE_A_100] = 0.012;
  m_mean[galf_EE_A_100_143] = 0.050;
  m_stddev[galf_EE_A_100_143] = 0.015;
  m_mean[galf_EE_A_100_217] = 0.110;
  m_stddev[galf_EE_A_100_217] = 0.033;
  m_mean[galf_EE_A_143] = 0.10;
  m_stddev[galf_EE_A_143] = 0.02;
  m_mean[galf_EE_A_143_217] = 0.240;
  m_stddev[galf_EE_A_143_217] = 0.048;
  m_mean[galf_EE_A_217] = 0.72;
  m_stddev[galf_EE_A_217] = 0.14;
  m_mean[galf_TE_A_100] = 0.140;
  m_stddev[galf_TE_A_100] = 0.042;
  m_mean[galf_TE_A_100_143] = 0.120;
  m_stddev[galf_TE_A_100_143] = 0.036;
  m_mean[galf_TE_A_100_217] = 0.30;
  m_stddev[galf_TE_A_100_217] = 0.09;
  m_mean[galf_TE_A_143] = 0.240;
  m_stddev[galf_TE_A_143] = 0.072;
  m_mean[galf_TE_A_143_217] = 0.600;
  m_stddev[galf_TE_A_143_217] = 0.018;
  m_mean[galf_TE_A_217] = 1.80;
  m_stddev[galf_TE_A_217] = 0.54;
#endif
  m_mean[A_planck] = 1.0;
  m_stddev[A_planck] = 0.0025;

  /****************************/
  /*  FINISH PARAMETER ALLOC  */
  /****************************/

  std::cout << "Out of " << TOTAL_PARAMS << " total parameters, we have:\n\t"
            << FREE_PARAMS << " free parameters,\n\t"
            << (FIXED_PARAMS - FREE_PARAMS) << " fixed parameters, and\n\t"
            << (TOTAL_PARAMS - FIXED_PARAMS) << " derived parameters."
            << std::endl;

#ifdef BAO_LIKE
  // Set BAO values
  // BAO_FUDGE = 1.0275; // from Planck 2013
  BAO_FUDGE = 1.0262; // from SDSS III [arXiv:1312.4877]

  // 6DF gives value of D_v divided by rs_rescale = 1.0268
  sixDF_z = 0.106;
  sixDF_mean = 0.327;
  sixDF_stddev = 0.015;

  // BOSS LOWZ gives value of D_v/rs_fid
  LOWZ_z = 0.32;
  LOWZ_mean = 8.47;
  LOWZ_stddev = 0.17;

  // BOSS CMASS gives value of D_v/rs_fid. Measures anisotropy
  CMASS_z = 0.57;
  CMASS_mean = 13.77;
  CMASS_stddev = 0.13;

  // SDSS DR7 MGS gives value of D_v/rs_fid
  MGS_z = 0.15;
  MGS_mean = 4.47;
  MGS_stddev = 0.16;
#endif
}

ClikPar::~ClikPar() {
  delete m_class_engine;
}

// Prepare CLASS engine
void ClikPar::initialise_CLASS(int max_l, struct pbh_external* pbh_info) {
  ClassParams default_params;

  /*** FREE PARAMETERS ***/

  // PBH DM
  default_params.add("Omega_pbh_ratio", 1.E-7);

  /*** CONSTANT PARAMETERS ***/

  // LCDM variables set to best-fit TTTEEE+lowP (2015)
  default_params.add("omega_b", 0.022252);
  default_params.add("omega_cdm", 0.11987);
  default_params.add("100*theta_s", 1.040778);
  default_params.add("tau_reio", 0.0789);
  default_params.add("ln10^{10}A_s", 3.0929);
  default_params.add("n_s", 0.96475);

  // PBH DM
  default_params.add("pbh_mass_dist", "pbh_delta");
  default_params.add("pbh_mass_mean", 1.E6);
  // default_params.add("pbh_mass_width", 1.E1);
  default_params.add("read pbh splines", false); // very important!!

  // Neutrino values
  default_params.add("N_ur", 2.0328);
  default_params.add("N_ncdm", 1);
  default_params.add("m_ncdm", 0.06); // MeV
  // default_params.add("T_ncdm", 0.71611);

  // Perturbation options for matter perturbation spectrum mPk
  // default_params.add("P_k_max_h/Mpc", 1.);
  // default_params.add("k_pivot", 0.05); // Mpc-1
  
  // Spectra output options
  default_params.add("output", "tCl,pCl,lCl"); // mPk
  default_params.add("lensing", true);
  default_params.add("l_max_scalars", max_l);
  default_params.add("format", "camb");

  // Initialise CLASS engine with external PBH info
  try {
    m_class_engine = new ClassEngine(default_params, pbh_info);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS initialisation in ClikPar "
              << "failed, throwing exception "
              << e.what()
              << std::endl;
    throw e;
  }
}

// Set free, fixed and derived parameter values
void ClikPar::scale_Cube(double* Cube) {
  // Function to transform variable appropriately
  double (*func)(double);

  for (int param = 0; param < TOTAL_PARAMS; param++) {
    func = m_transform[param];

    if (param < FREE_PARAMS) { // free parameter
      Cube[param] = func(m_min[param] + (m_max[param] - m_min[param]) * Cube[param]);
    }
    else if (param < FIXED_PARAMS) { // fixed parameter
      Cube[param] = func(m_min[param]);
    }
    else { // derived parameter
      // Must be set to something in case CLASS fails
      Cube[param] = func(m_min[param]);
    }
  }
}

// Set derived parameter values in Cube
void ClikPar::set_derived_params(double* Cube) {
  // Non-standard CLASS routines
  Cube[H0] = m_class_engine->get_H0();
  Cube[Omega_b] = m_class_engine->get_Omega_b();
  Cube[Omega_cdm] = m_class_engine->get_Omega_cdm();
  Cube[Omega_L] = m_class_engine->get_Omega_L();
  Cube[Omega_g] = m_class_engine->get_Omega_g();
  Cube[sigma8] = m_class_engine->get_sigma8();
  Cube[age] = m_class_engine->get_age();
  Cube[conf_age] = m_class_engine->get_conf_age();

  // Standard CLASS routines
  Cube[z_drag] = m_class_engine->z_drag();
  Cube[rs_drag] = m_class_engine->rs_drag();
}

ClassEngine* ClikPar::get_CLASS() {
  return m_class_engine;
}

#ifdef BAO_LIKE
// *** The returned loglike value must always be positive ***
// Using values from CosmoMC
double ClikPar::calculate_BAO_likelihood() const {
  double loglike = 0.0;

  // Need rs_fid as computed by CLASS, too.
  // Using fudge factor to match Eisenstein & Hu (1998)
  double rs_fid = m_class_engine->rs_drag() * BAO_FUDGE;

  // The actual CLASS variables to compare against
  double sixDF_Dv = m_class_engine->get_Dv(sixDF_z);
  double LOWZ_Dv_rs = m_class_engine->get_Dv(LOWZ_z)/rs_fid;
  double CMASS_Dv_rs = m_class_engine->get_Dv(CMASS_z)/rs_fid;
  double MGS_Dv_rs = m_class_engine->get_Dv(MGS_z)/rs_fid;

  // Compute Gaussian likelihoods
  loglike += pow(sixDF_Dv - sixDF_mean, 2.0) / (2.0 * pow(sixDF_stddev, 2.0));
  loglike += pow(LOWZ_Dv_rs - LOWZ_mean, 2.0) / (2.0 * pow(LOWZ_stddev, 2.0));
  loglike += pow(CMASS_Dv_rs - CMASS_mean, 2.0) / (2.0 * pow(CMASS_stddev, 2.0));
  loglike += pow(MGS_Dv_rs - MGS_mean, 2.0) / (2.0 * pow(MGS_stddev, 2.0));

  return loglike;
}
#endif

// *** The returned loglike value must always be positive ***
double ClikPar::calculate_extra_priors(double* Cube) const {
  double loglike = 0.0;

  // Calculate Gaussian priors
  for (int param = 0; param < TOTAL_PARAMS; ++param) {
    // Only use Gaussian priors for variables that are Gaussian
    if (m_has_gaussian_prior[param]) {
      loglike += pow(Cube[param] - m_mean[param], 2.0) / (2.0 * pow(m_stddev[param], 2.0));
    }
  }

  // Calculate SZ degeneracy prior
#ifndef LITE_HI_L
  double SZ_val = Cube[ksz_norm] + 1.6 * Cube[A_sz];
  double SZ_mean = 9.5;
  double SZ_stddev = 3.0;

  loglike += pow(SZ_val - SZ_mean, 2.0) / (2.0 * pow(SZ_stddev, 2.0));
#endif

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
