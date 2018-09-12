#ifndef LOGLIKE_H
#define LOGLIKE_H

#include "ClassEngine.hpp"
#include "ClikObject.hpp"

#include <vector>

const double MIN_LOGLIKE = -1E90;

// Extra prior values
const double SZ_MEAN = 9.5;
const double SZ_STDDEV = 3.0;
const double H0_MIN = 20;
const double H0_MAX = 100;

// Constant locations of derived parameters in free_params 
const int H0_LOC = static_cast<int>(H0 - FIXED_PARAM_AMT);
const int Omega_b_LOC = static_cast<int>(Omega_b - FIXED_PARAM_AMT);
const int Omega_cdm_LOC = static_cast<int>(Omega_cdm - FIXED_PARAM_AMT);
const int Omega_L_LOC = static_cast<int>(Omega_L - FIXED_PARAM_AMT);
const int Omega_g_LOC = static_cast<int>(Omega_g - FIXED_PARAM_AMT);
const int sigma8_LOC = static_cast<int>(sigma8_p - FIXED_PARAM_AMT);
const int age_LOC = static_cast<int>(age - FIXED_PARAM_AMT);
const int conf_age_LOC = static_cast<int>(conf_age - FIXED_PARAM_AMT);
const int z_drag_LOC = static_cast<int>(z_drag - FIXED_PARAM_AMT);
const int rs_drag_LOC = static_cast<int>(rs_drag - FIXED_PARAM_AMT);

typedef struct {
  ClassEngine* engine;
  std::vector<ClikObject*> clik_objs;
  std::vector<unsigned> cl_ls;
} likelihood_context;

/**
 * Calculate full likelihood
 *
 * free_params : array of values needed for likelihood calculation. These are:
 *              * free CLASS parameters
 *              * free nuisance parameters
 */
double calculate_full_likelihood(double free_params[], likelihood_context* bundle);

double calculate_extra_likelihood(double free_params[]);

void set_derived_params(double free_params[], ClassEngine*& class_engine);

#endif
