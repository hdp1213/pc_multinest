#include "loglike.hpp"

#include <iostream>

using namespace std;

double
calculate_full_likelihood(double free_params[], likelihood_context* bundle) {
  double loglike = 0.0;

  vector<double> lcdm_params(&free_params[0], &free_params[0] + CLASS_PARAM_AMT);

  // Print lcdm parameter values
  cout.setf(ios::scientific);
  auto precis = cout.precision(16);

  for (auto param : lcdm_params) {
    cout << param << " ";
  }
  cout << endl;

  cout.unsetf(ios::scientific);
  cout.precision(precis);

  vector<double> class_cls[CLASS_CL_AMT];

  try {
    bundle->engine->update_parameters(lcdm_params);
    bundle->engine->get_Cls(bundle->cl_ls, class_cls[0], class_cls[3], class_cls[1], class_cls[2]);
  }
  catch (exception const& e) {
    cerr << "[ERROR] " << e.what() << endl;
    return MIN_LOGLIKE;
  }

  // Extract derived parameters from CLASS engine
  set_derived_params(free_params, bundle->engine);

  for (auto clik_obj : bundle->clik_objs) {
    loglike += clik_obj->calculate_likelihood(free_params, class_cls);
  }

  loglike += calculate_extra_likelihood(free_params);

  return loglike;
}

double
calculate_extra_likelihood(double free_params[]) {
  // Calculate flat [20,100] prior on H0
  if (free_params[H0_LOC] < H0_MIN or free_params[H0_LOC] > H0_MAX) {
    return MIN_LOGLIKE;
  }

  double loglike = 0.0;

  // Calculate Gaussian priors
  for (auto gauss : g_gauss) {
    double param_value = get_value(gauss.param, free_params);
    loglike -= pow((param_value - gauss.mean)/gauss.stddev, 2.0) / 2.0;
  }

  // Calculate SZ degeneracy prior
#ifndef LITE_HI_L
  double SZ_val = get_value(ksz_norm, free_params) + 1.6 * get_value(A_sz, free_params);

  loglike -= pow((SZ_val - SZ_MEAN)/SZ_STDDEV, 2.0) / 2.0;
#endif

#ifdef DBUG
  cout << "extra loglike = " << loglike << endl;
#endif

  return loglike;
}

void
set_derived_params(double free_params[], ClassEngine*& class_engine) {
  // Non-standard CLASS routines
  free_params[H0_LOC] = class_engine->get_H0();
  free_params[Omega_b_LOC] = class_engine->get_Omega_b();
  free_params[Omega_cdm_LOC] = class_engine->get_Omega_cdm();
  free_params[Omega_L_LOC] = class_engine->get_Omega_L();
  free_params[Omega_g_LOC] = class_engine->get_Omega_g();
  free_params[sigma8_LOC] = class_engine->get_sigma8();
  free_params[age_LOC] = class_engine->get_age();
  free_params[conf_age_LOC] = class_engine->get_conf_age();

  // Standard CLASS routines
  free_params[z_drag_LOC] = class_engine->z_drag();
  free_params[rs_drag_LOC] = class_engine->rs_drag();
}
