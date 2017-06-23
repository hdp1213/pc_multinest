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
  FREE_PARAMS,
  // Fixed parameters (PLC TT nuisance)
  cib_index = FREE_PARAMS,
  FIXED_PARAMS,
  // Derived parameters (LCDM)
  H0 = FIXED_PARAMS,
#else // remove all nuisance params for LITE_HI_L except A_planck
  FREE_PARAMS,
  // Derived parameters (LCDM)
  H0 = FREE_PARAMS,
#endif
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
  TOTAL_PARAMS
};
