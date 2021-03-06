#ifndef PARAMETER_ENUM_HPP
#define PARAMETER_ENUM_HPP

typedef enum
{
  // Free parameters (CLASS)
  omega_b,
  omega_cdm,
  hundredxtheta_s,
  tau_reio,
  ln10_10_A_s,
  n_s,
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
#else
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

#else
#error "Including multiple parameter enumerations"
#endif
