#ifndef PARAMETER_H
#define PARAMETER_H

#include <assert.h>
#include <forward_list>

// Parameter enumeration
typedef enum
{
  // Free parameters (CLASS)
  pbh_frac,
  pbh_mass,
  // omega_b,
  // omega_cdm,
  // hundredxtheta_s,
  // tau_reio,
  // ln10_10_A_s,
  // n_s,
  UP_TO_CLASS_PARAMS,
  // Free parameters (PLC TT nuisance)
  A_planck = UP_TO_CLASS_PARAMS,
#ifndef LITE_HI_L
  A_cib_217,
  xi_sz_cib,
  A_sz,
  ps_A_100_100,
  ps_A_143_143,
  ps_A_143_217,
  ps_A_217_217,
  ksz_norm,
  gal545_A_100,
  gal545_A_143,
  gal545_A_143_217,
  gal545_A_217,
  // Free parameters (PLC TTTEEE nuisance)
  galf_EE_A_100,
  galf_EE_A_100_143,
  galf_EE_A_100_217,
  galf_EE_A_143,
  galf_EE_A_143_217,
  galf_EE_A_217,
  galf_TE_A_100,
  galf_TE_A_100_143,
  galf_TE_A_100_217,
  galf_TE_A_143,
  galf_TE_A_143_217,
  galf_TE_A_217,
  calib_100T,
  calib_217T,
#endif
  UP_TO_FREE_PARAMS,
#ifndef LITE_HI_L
  // Fixed parameters (PLC TT nuisance)
  cib_index = UP_TO_FREE_PARAMS,
  // Fixed parameters (PLC TTTEEE nuisance)
  galf_EE_index,
  galf_TE_index,
  bleak_epsilon_0_0T_0E,
  bleak_epsilon_1_0T_0E,
  bleak_epsilon_2_0T_0E,
  bleak_epsilon_3_0T_0E,
  bleak_epsilon_4_0T_0E,
  bleak_epsilon_0_0T_1E,
  bleak_epsilon_1_0T_1E,
  bleak_epsilon_2_0T_1E,
  bleak_epsilon_3_0T_1E,
  bleak_epsilon_4_0T_1E,
  bleak_epsilon_0_0T_2E,
  bleak_epsilon_1_0T_2E,
  bleak_epsilon_2_0T_2E,
  bleak_epsilon_3_0T_2E,
  bleak_epsilon_4_0T_2E,
  bleak_epsilon_0_1T_1E,
  bleak_epsilon_1_1T_1E,
  bleak_epsilon_2_1T_1E,
  bleak_epsilon_3_1T_1E,
  bleak_epsilon_4_1T_1E,
  bleak_epsilon_0_1T_2E,
  bleak_epsilon_1_1T_2E,
  bleak_epsilon_2_1T_2E,
  bleak_epsilon_3_1T_2E,
  bleak_epsilon_4_1T_2E,
  bleak_epsilon_0_2T_2E,
  bleak_epsilon_1_2T_2E,
  bleak_epsilon_2_2T_2E,
  bleak_epsilon_3_2T_2E,
  bleak_epsilon_4_2T_2E,
  bleak_epsilon_0_0E_0E,
  bleak_epsilon_1_0E_0E,
  bleak_epsilon_2_0E_0E,
  bleak_epsilon_3_0E_0E,
  bleak_epsilon_4_0E_0E,
  bleak_epsilon_0_0E_1E,
  bleak_epsilon_1_0E_1E,
  bleak_epsilon_2_0E_1E,
  bleak_epsilon_3_0E_1E,
  bleak_epsilon_4_0E_1E,
  bleak_epsilon_0_0E_2E,
  bleak_epsilon_1_0E_2E,
  bleak_epsilon_2_0E_2E,
  bleak_epsilon_3_0E_2E,
  bleak_epsilon_4_0E_2E,
  bleak_epsilon_0_1E_1E,
  bleak_epsilon_1_1E_1E,
  bleak_epsilon_2_1E_1E,
  bleak_epsilon_3_1E_1E,
  bleak_epsilon_4_1E_1E,
  bleak_epsilon_0_1E_2E,
  bleak_epsilon_1_1E_2E,
  bleak_epsilon_2_1E_2E,
  bleak_epsilon_3_1E_2E,
  bleak_epsilon_4_1E_2E,
  bleak_epsilon_0_2E_2E,
  bleak_epsilon_1_2E_2E,
  bleak_epsilon_2_2E_2E,
  bleak_epsilon_3_2E_2E,
  bleak_epsilon_4_2E_2E,
  calib_100P,
  calib_143P,
  calib_217P,
  A_pol,
  UP_TO_FIXED_PARAMS,
#else // remove all nuisance params for LITE_HI_L except A_planck
  UP_TO_FIXED_PARAMS = UP_TO_FREE_PARAMS,
#endif
  // Derived parameters (LCDM)
  H0 = UP_TO_FIXED_PARAMS,
  Omega_b,
  Omega_cdm,
  Omega_L,
  Omega_g,
  sigma8_p,
  age,
  conf_age,
  z_drag,
  rs_drag,
  // Final total number of parameters
  TOTAL_PARAM_AMT,
  // Summary constants
  CLASS_PARAM_AMT = UP_TO_CLASS_PARAMS,
  FREE_PARAM_AMT = UP_TO_FREE_PARAMS,
  FIXED_PARAM_AMT = UP_TO_FIXED_PARAMS - UP_TO_FREE_PARAMS,
  DERIVED_PARAM_AMT = TOTAL_PARAM_AMT - UP_TO_FIXED_PARAMS
} param_t;

// Transformation function typedef
typedef double (*trans_t)(double);
inline double identity(double x) {return x;};

// Free flat priors
extern double g_min[FREE_PARAM_AMT];
extern double g_max[FREE_PARAM_AMT];

// Fixed values
extern double g_value[FIXED_PARAM_AMT];

// Transforms
extern trans_t g_transform[FREE_PARAM_AMT];

// Struct for storing Gaussian priors
typedef struct {
  param_t param;
  double mean;
  double stddev;
} gauss_t;

// Gaussian priors
extern std::forward_list<gauss_t> g_gauss;

// Actual functions we use
void initialise_param_arrays();

inline void add_fixed_param(const param_t param, const double value) {
  assert(param >= UP_TO_FREE_PARAMS);
  g_value[param - UP_TO_FREE_PARAMS] = value;
}

inline void add_free_param(const param_t param, const double min, const double max) {
  assert(min < max);
  g_min[param] = min; g_max[param] = max;
}

inline void add_gauss_prior(const param_t param, const double mean, const double stddev) {
  assert(stddev > 0.0);
  gauss_t gauss_vars = {.param = param, .mean = mean, .stddev = stddev};
  g_gauss.push_front(gauss_vars);
}

inline void add_transform(const param_t param, const trans_t trans) {
  g_transform[param] = trans;
}

#endif
