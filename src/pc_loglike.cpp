#include "pc_loglike.hpp"

#include <exception>
#include <iostream>
#include <iomanip>

double pc_loglike(std::vector<clik_struct*>& clik_objs,
                  ClassEngine* class_engine,
                  const std::vector<unsigned>& cl_ls,
                  double raw_params[],
                  double in_vals[]) {
  std::vector<double> class_params;
  std::vector<double> cl_tt, cl_ee, cl_bb, cl_te;
  std::vector<std::vector<double> > class_cls;

  double loglike = 0.0;
  class_cls.reserve(CLASS_CL_AMT);

  class_params.push_back(in_vals[pbh_frac]);
  class_params.push_back(in_vals[pbh_mass]);

  std::cout << std::scientific << std::setprecision(16)
            << in_vals[pbh_frac] << " "
            << in_vals[pbh_mass] << std::endl;

  // Run CLASS
  try {
    class_engine->update_parameters(class_params);
  }
  catch (std::exception const& e) {
    std::cerr << "[ERROR] CLASS failed, throwing exception "
              << e.what()
              << std::endl;

    return MIN_LOGLIKE;
  }

  // Extract CLASS spectra
  try {
    class_engine->get_Cls(cl_ls, cl_tt, cl_te, cl_ee, cl_bb);
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
  for (std::vector<clik_struct*>::iterator clik_obj_it = clik_objs.begin();
       clik_obj_it != clik_objs.end(); ++clik_obj_it) {
    loglike += calculate_PLC_likelihood(*clik_obj_it,
                                        in_vals,
                                        class_cls);
  }

  // Calculate extra likelihoods
  // Some could depend on derived parameters, so also need raw_params
  loglike += calculate_extra_likelihoods(in_vals, raw_params);

  return loglike;
}

double calculate_PLC_likelihood(clik_struct* clik_obj,
                                double in_vals[],
                                std::vector<std::vector<double> >& class_spectra) {
  // First must create cl_and_pars array
  double cl_and_pars[clik_obj->cap_size];

  int cap_ind = 0;
  for (int cl_ind = 0; cl_ind < CL_AMT; ++cl_ind) {
    if (clik_obj->has_cl[cl_ind]) {  // if spectra exists
      for (int l = 0; l <= clik_obj->max_l; ++l) {
        // Correct for PLC wanting l=0,1 multipoles
        // These cannot be computed by CLASS and so are set to zero
        if (l < CLASS_MIN_L) {
          cl_and_pars[cap_ind] = ZERO_CL;
        } else {
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
  for (std::vector<param_t>::iterator param_it = clik_obj->nuis_pars.begin();
       param_it != clik_obj->nuis_pars.end(); ++param_it) {
    cl_and_pars[cap_ind] = in_vals[*param_it];

    cap_ind++;
  }

  // Now must use that with the clik_id to compute log likelihood
  error* err = initError();

  double loglike = MIN_LOGLIKE;
  loglike = clik_compute(clik_obj->clik_id, cl_and_pars, &err);
  quitOnError(err, __LINE__, stderr);

#ifdef DBUG
  std::cout << "[calculate_PLC_likelihood] Calculated loglike of "
            << loglike << std::endl;
#endif

  free(err);

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

#ifdef DBUG
  std::cout << "[calculate_extra_likelihoods] Calculated loglike of " << loglike << std::endl;
#endif

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
  // raw_params[sigma8_LOC] = class_engine->get_sigma8();
  raw_params[age_LOC] = class_engine->get_age();
  raw_params[conf_age_LOC] = class_engine->get_conf_age();

  // Standard CLASS routines
  raw_params[z_drag_LOC] = class_engine->z_drag();
  raw_params[rs_drag_LOC] = class_engine->rs_drag();
}
