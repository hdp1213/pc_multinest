#ifndef PC_LOGLIKE
#define PC_LOGLIKE

#include "init_plc.h"
#include "init_params.h"

#include <exception>
#include <iostream>
#include <vector>

#ifndef PARAM_SET
#error no parameter set included!
#endif

/* "Actual" constants */
const double MIN_LOGLIKE = -1E90;
const double ZERO_CL = 0.0;
const int CLASS_CL_AMT = 4;

/* Misc. gaussian likelihood values */
const double SZ_MEAN = 9.5;
const double SZ_STDDEV = 3.0;

/* Constant locations of derived parameters in raw_params */
const int H0_LOC = static_cast<int>(H0 - FIXED_PARAM_AMT);
const int Omega_b_LOC = static_cast<int>(Omega_b - FIXED_PARAM_AMT);
const int Omega_cdm_LOC = static_cast<int>(Omega_cdm - FIXED_PARAM_AMT);
const int Omega_L_LOC = static_cast<int>(Omega_L - FIXED_PARAM_AMT);
const int Omega_g_LOC = static_cast<int>(Omega_g - FIXED_PARAM_AMT);
const int sigma8_LOC = static_cast<int>(sigma8 - FIXED_PARAM_AMT);
const int age_LOC = static_cast<int>(age - FIXED_PARAM_AMT);
const int conf_age_LOC = static_cast<int>(conf_age - FIXED_PARAM_AMT);
const int z_drag_LOC = static_cast<int>(z_drag - FIXED_PARAM_AMT);
const int rs_drag_LOC = static_cast<int>(rs_drag - FIXED_PARAM_AMT);

double calculate_PLC_likelihood(clik_struct*& clik_obj,
                                double in_vals[],
                                std::vector<std::vector<double> >& class_spectra);
// implicitly uses global arrays
double calculate_extra_likelihoods(double in_vals[],
                                   double raw_params[]);
// raw_params[] contains only free and derived parameters
void set_derived_params(double raw_params[],
                        ClassEngine*& class_engine);



/**
 * What do we need to compute the log likelihood for plc_class?
 * - the initialised .clik objects
 * - a CLASS engine
 * - the actual parameter values being considered for the particular
 *   point in parameter space. This contains just free parameters though
 * - the fixed values, which are always constant and never change
 * - a way of telling which parameter is free or fixed
 * - any gaussian priors!!!
 * 
 * Things that have already happened prior to this:
 * - .clik objects have been initialised
 * - CLASS engine has been initialised
 *
 * I feel like we need a clik struct that contains all of the
 * useful information already extracted.
**/
double pc_loglike(std::vector<clik_struct*>& clik_objs,
                  ClassEngine*& class_engine,
                  std::vector<unsigned>& cl_ls,
                  double raw_params[],
                  double in_vals[]) {

  std::vector<double> class_params;
  std::vector<double> cl_tt, cl_ee, cl_bb, cl_te;
  std::vector<std::vector<double> > class_cls;

  double loglike = 0.0;
  class_cls.reserve(CLASS_CL_AMT);

  class_params.push_back(in_vals[pbh_frac]);
  class_params.push_back(in_vals[omega_b]);
  class_params.push_back(in_vals[omega_cdm]);
  class_params.push_back(in_vals[hundredxtheta_s]);
  class_params.push_back(in_vals[tau_reio]);
  class_params.push_back(in_vals[ln10_10_A_s]);
  class_params.push_back(in_vals[n_s]);

  // Run CLASS
  try {
    class_engine->updateParValues(class_params);
  }
  catch (std::exception const& e) {
    std::cerr << "[ERROR] CLASS failed, throwing exception "
              << e.what()
              << std::endl;

    return MIN_LOGLIKE;
  }

  // Extract CLASS spectra
  try {
    class_engine->getCls(cl_ls, cl_tt, cl_te, cl_ee, cl_bb);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] Spectra extraction unsuccessful, CLASS threw "
              << e.what()
              << std::endl;

    return MIN_LOGLIKE;
  }

  // Save derived parameters into raw_params array
  set_derived_params(raw_params, class_engine);

  // Collate CLASS spectra
  class_cls.push_back(cl_tt);
  class_cls.push_back(cl_ee);
  class_cls.push_back(cl_bb);
  class_cls.push_back(cl_te);

  // Calculate total likelihoods from PLC
  for (std::vector<clik_struct*>::iterator clik_obj_it = clik_objs.begin(); clik_obj_it != clik_objs.end(); ++clik_obj_it) {
    loglike += calculate_PLC_likelihood(*clik_obj_it,
                                        in_vals,
                                        class_cls);
  }

  // Calculate extra likelihoods
  // Some could depend on derived parameters, so also need raw_params
  loglike += calculate_extra_likelihoods(in_vals, raw_params);

  return loglike;
}

double calculate_PLC_likelihood(clik_struct*& clik_obj,
                                double in_vals[],
                                std::vector<std::vector<double> >& class_spectra) {
  // First must create cl_and_pars array
  double cl_and_pars[clik_obj->cap_size];

  int cap_ind = 0;
  for (int cl_ind = 0; cl_ind < CL_AMT; ++cl_ind) {
    int old_cl_ind = cl_ind;

    if (clik_obj->has_cl[cl_ind]) { // if spectra exists
      for (int l = 0; l <= clik_obj->max_l; ++l) {
        // Correct for PLC wanting l=0,1 multipoles
        // These cannot be computed by CLASS and so are set to zero
        if (l < CLASS_MIN_L) {
          cl_and_pars[cap_ind] = ZERO_CL;
        }
        else {
          cl_and_pars[cap_ind] = class_spectra[cl_ind][l - CLASS_MIN_L];
        }

        cap_ind++;
      }
    }
  }

  // std::cout << "[calculate_PLC_likelihood] Adding nuisance parameters..." << std::endl;

  // Then add nuisance parameters at the end
  // Nuisance parameters need to be stored in the same order they
  //  initialise the clik_struct with
  for (std::vector<param_t>::iterator param_it = clik_obj->nuis_pars.begin(); param_it != clik_obj->nuis_pars.end(); ++param_it) {

    cl_and_pars[cap_ind] = in_vals[*param_it];
    // std::cout << *param_it << " = " << in_vals[*param_it] << std::endl;

    cap_ind++;
  }

  // Now must use that with the clik_id to compute log likelihood
  error* _err = initError();
  error** err = &_err;

  double loglike = MIN_LOGLIKE;
  loglike = clik_compute(clik_obj->clik_id, cl_and_pars, err);
  quitOnError(*err, __LINE__, stderr);

#ifdef DBUG
  std::cout << "[calculate_PLC_likelihood] Calculated loglike of " << loglike << std::endl;
#endif

  free(_err);

  return loglike;
}

// loglike calculated here is added to the total loglike
// so watch your signs!
double calculate_extra_likelihoods(double in_vals[],
                                   double raw_params[]) {
  double loglike = 0.0;

  // Calculate flat [20,100] prior on H0
  double H0_min = 20.0, H0_max = 100.0;

  if (raw_params[H0_LOC] < H0_min or raw_params[H0_LOC] > H0_max) {
    return MIN_LOGLIKE;
  }

  // Calculate Gaussian priors
  for (int param_ind = 0; param_ind < UP_TO_FIXED_PARAMS; ++param_ind) {
    // Only use Gaussian priors for variables that are Gaussian
    if (m_has_gaussian_prior[param_ind]) {
      loglike -= pow((in_vals[param_ind] - m_mean[param_ind])/m_stddev[param_ind], 2.0) / 2.0;
    }
  }

  // Calculate SZ degeneracy prior
#ifndef LITE_HI_L
  double SZ_val = in_vals[ksz_norm] + 1.6 * in_vals[A_sz];

  loglike -= pow((SZ_val - SZ_MEAN)/SZ_STDDEV, 2.0) / 2.0;
#endif

  // Calculate any other priors needed ...

  return loglike;
}

void set_derived_params(double raw_params[],
                        ClassEngine*& class_engine) {
  // Non-standard CLASS routines
  raw_params[H0_LOC] = class_engine->get_H0();
  raw_params[Omega_b_LOC] = class_engine->get_Omega_b();
  raw_params[Omega_cdm_LOC] = class_engine->get_Omega_cdm();
  raw_params[Omega_L_LOC] = class_engine->get_Omega_L();
  raw_params[Omega_g_LOC] = class_engine->get_Omega_g();
  raw_params[sigma8_LOC] = class_engine->get_sigma8();
  raw_params[age_LOC] = class_engine->get_age();
  raw_params[conf_age_LOC] = class_engine->get_conf_age();

  // Standard CLASS routines
  raw_params[z_drag_LOC] = class_engine->z_drag();
  raw_params[rs_drag_LOC] = class_engine->rs_drag();
}


#endif
