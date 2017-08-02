#include "ClikPar.h"

#include <math.h> // for pow()
#include <iostream> // for std::cerr
#include <exception> // for std::exception

// Wrapped methods for scale_free_params()
inline double pow10(double x) { return pow(10., x); }
inline double self(double x) { return x; }

ClikPar::ClikPar() : m_class_engine(0) {
  // Initialise defaults
  for (int param = 0; param < TOTAL_PARAM_AMT; ++param) {
    m_has_gaussian_prior[param] = false;
    m_is_log10[param] = false;
    m_mean[param] = -9999.0;
    m_stddev[param] = -9999.0;
  }

  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    m_min[param] = -9999.0;
    m_max[param] = -9999.0;
  }

  for (int param = 0; param < FIXED_PARAM_AMT; ++param) {
    m_value[param] = -9999.0;
  }

// Define flat priors
// #include "TTTEEE+lowP_pbh_fixedLCDM-flat.cc"
#include "TTTEEE+lowP-flat.cc"

  // Do checks on min and max values for free parameters
  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    if (m_min[param] > m_max[param]) {
      std::cerr << "[ERROR]: min value larger than max value "
                << "for parameter #"
                << param
                << ", a free parameter. Have you got the "
                << "correct flat priors for this parameter?"
                << std::endl;
      throw std::exception();
    }

    else if (m_min[param] == m_max[param]) {
      std::cerr << "[ERROR]: min and max values are identical "
                << "for parameter #"
                << param
                << ", a free parameter. Have you got the "
                << "correct flat priors for this parameter?"
                << std::endl;
      throw std::exception();
    }
  }

  // Check if parameter priors are in log10 space.
  // If so, set transformation function to the power of 10 such that Cube[param] is set to the correct value
  for (int param = 0; param < TOTAL_PARAM_AMT; ++param) {
    m_transform[param] = m_is_log10[param] ? &pow10 : &self;
  }

// Define Gaussian priors
// #include "TTTEEE+lowP_pbh_fixedLCDM-gauss.cc"
#include "TTTEEE+lowP-gauss.cc"

  /****************************/
  /*  FINISH PARAMETER ALLOC  */
  /****************************/

  std::cout << "Out of " << TOTAL_PARAM_AMT << " total parameters, "
            << "we have:\n\t"
            << FREE_PARAM_AMT << " free parameters,\n\t"
            << FIXED_PARAM_AMT << " fixed parameters, and\n\t"
            << DERIVED_PARAM_AMT << " derived parameters."
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
  // default_params.add("Omega_pbh_ratio", 1.E-7);

  /*** FREE/FIXED PARAMETERS ***/

  // LCDM variables set to best-fit TTTEEE+lowP (2015)
  default_params.add("omega_b", 0.022252);
  default_params.add("omega_cdm", 0.11987);
  default_params.add("100*theta_s", 1.040778);
  default_params.add("tau_reio", 0.0789);
  default_params.add("ln10^{10}A_s", 3.0929);
  default_params.add("n_s", 0.96475);

  // PBH DM
  default_params.add("pbh_mass_dist", "pbh_none");
  default_params.add("read pbh splines", false);
  /*
  default_params.add("pbh_mass_dist", "pbh_delta");
  default_params.add("pbh_mass_mean", 1.E6);
  // default_params.add("pbh_mass_width", 1.E1);
  default_params.add("read pbh splines", false); // very important!!
  */

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

// MultiNest parameters used in CLASS
// Must be pushed in the same order as initialised in initialise_CLASS()
void ClikPar::run_CLASS(double* free_params) {
  std::vector<double> class_params;

  // class_params.push_back(get_par(pbh_frac, free_params));
  class_params.push_back(get_par(omega_b, free_params));
  class_params.push_back(get_par(omega_cdm, free_params));
  class_params.push_back(get_par(hundredxtheta_s, free_params));
  class_params.push_back(get_par(tau_reio, free_params));
  class_params.push_back(get_par(ln10_10_A_s, free_params));
  class_params.push_back(get_par(n_s, free_params)); // k_0 = 0.05 Mpc^-1 by default

  m_class_engine->updateParValues(class_params);
}

void ClikPar::get_CLASS_spectra(std::vector<unsigned>& cl_ls,
      std::vector<double>& cl_tt,
      std::vector<double>& cl_te,
      std::vector<double>& cl_ee,
      std::vector<double>& cl_bb) {
  m_class_engine->getCls(cl_ls, cl_tt, cl_te, cl_ee, cl_bb);
}

// Set free, fixed and derived parameter values
void ClikPar::scale_free_params(double* free_params) {
  // Function to transform variable appropriately
  double (*func)(double);

  for (int param = 0; param < FREE_PARAM_AMT; param++) {
    func = m_transform[param];

    free_params[param] = func(m_min[param] + (m_max[param] - m_min[param]) * free_params[param]);
  }
}

// Set derived parameter values in all_params
// Should be called before any likelihood functions are
void ClikPar::set_derived_params(double* all_params) {
  // Non-standard CLASS routines
  all_params[H0] = m_class_engine->get_H0();
  all_params[Omega_b] = m_class_engine->get_Omega_b();
  all_params[Omega_cdm] = m_class_engine->get_Omega_cdm();
  all_params[Omega_L] = m_class_engine->get_Omega_L();
  all_params[Omega_g] = m_class_engine->get_Omega_g();
  all_params[sigma8] = m_class_engine->get_sigma8();
  all_params[age] = m_class_engine->get_age();
  all_params[conf_age] = m_class_engine->get_conf_age();

  // Standard CLASS routines
  all_params[z_drag] = m_class_engine->z_drag();
  all_params[rs_drag] = m_class_engine->rs_drag();
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
  loglike += pow((sixDF_Dv - sixDF_mean)/sixDF_stddev, 2.0) / 2.0;
  loglike += pow((LOWZ_Dv_rs - LOWZ_mean)/LOWZ_stddev, 2.0) / 2.0;
  loglike += pow((CMASS_Dv_rs - CMASS_mean)/CMASS_stddev, 2.0) / 2.0;
  loglike += pow((MGS_Dv_rs - MGS_mean)/MGS_stddev, 2.0) / 2.0;

  return loglike;
}
#endif

// *** The returned loglike value must always be positive ***
double ClikPar::calculate_extra_likelihoods(double* in_params) const {
  double loglike = 0.0;

  // Calculate Gaussian priors
  for (int param = 0; param < TOTAL_PARAM_AMT; ++param) {
    // Only use Gaussian priors for variables that are Gaussian
    if (m_has_gaussian_prior[param]) {
      loglike += pow((in_params[param] - m_mean[param])/m_stddev[param], 2.0) / 2.0;
    }
  }

  // Calculate SZ degeneracy prior
#ifndef LITE_HI_L
  double SZ_val = in_params[ksz_norm] + 1.6 * in_params[A_sz];
  double SZ_mean = 9.5;
  double SZ_stddev = 3.0;

  loglike += pow((SZ_val - SZ_mean)/SZ_stddev, 2.0) / 2.0;
#endif

  // Calculate flat [20,100] prior on H0
  double H0_min = 20.0, H0_max = 100.0;

  if (in_params[H0] < H0_min or in_params[H0] > H0_max) {
    return 1E90;
  }

  // Calculate any other priors needed ...

  return loglike;
}

double ClikPar::calculate_flat_prior() const {
  double res = 1.0;

  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    res *= (m_max[param] - m_min[param]);
  }

  return 1.0/res;
}

double* ClikPar::get_lower_bounds() {
  return m_min;
}

double* ClikPar::get_upper_bounds() {
  return m_max;
}

double* ClikPar::get_fixed_params() {
  return m_value;
}

double ClikPar::get_par(param_t param, double* free_params) const {
  if (param < UP_TO_FREE_PARAMS) { // if free
    return free_params[param];
  }
  else if (param < UP_TO_FIXED_PARAMS) { // if fixed
    return m_value[param - UP_TO_FREE_PARAMS];
  }
  else {
    std::cerr << "[ERROR]: You have tried getting the value "
              << "of parameter #"
              << param
              << ", a derived parameter. This cannot be done "
              << "with get_par()."
              << std::endl;
    throw std::exception();
  }
}
