#ifndef PARAM_SET
#define PARAM_SET 0
// This enum follows Planck ordering (to a degree)
enum param_t {
  // Free parameters (CLASS)
  omega_b,
  omega_cdm,
  hundredxtheta_s,
  tau_reio,
  ln10_10_A_s,
  n_s,
  // Free parameters (PLC TT nuisance)
  A_planck,
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
#endif
  UP_TO_FREE_PARAMS,
#ifndef LITE_HI_L
  // Fixed parameters (PLC TT nuisance)
  cib_index = UP_TO_FREE_PARAMS,
  UP_TO_FIXED_PARAMS,
#else // remove all nuisance params for LITE_HI_L except A_planck
  // Derived parameters (LCDM)
  UP_TO_FIXED_PARAMS = UP_TO_FREE_PARAMS,
#endif
  // Derived parameters (LCDM)
  H0 = UP_TO_FIXED_PARAMS,
  Omega_b,
  Omega_cdm,
  Omega_L,
  Omega_g,
  sigma8,
  age,
  conf_age,
  z_drag,
  rs_drag,
  // Final total number of parameters
  TOTAL_PARAM_AMT,
  // Summary constants
  FREE_PARAM_AMT = UP_TO_FREE_PARAMS,
  FIXED_PARAM_AMT = UP_TO_FIXED_PARAMS - UP_TO_FREE_PARAMS,
  DERIVED_PARAM_AMT = TOTAL_PARAM_AMT - UP_TO_FIXED_PARAMS
};
#endif
