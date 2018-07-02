/**********************************/
/*  UNIFORM (FLAT) PRIOR SETTING  */
/**********************************/

#if PARAM_SET == 10

#ifndef FLAT_PRIOR_SET
#define FLAT_PRIOR_SET

// Free LCDM parameters according to TTTEEE+lowP.pars
m_min[omega_b] = 0.016;          m_max[omega_b] = 0.028;
m_min[omega_cdm] = 0.108;        m_max[omega_cdm] = 0.130;
m_min[hundredxtheta_s] = 1.039;  m_max[hundredxtheta_s] = 1.043;
m_min[tau_reio] = 0.01;          m_max[tau_reio] = 0.15;
m_min[ln10_10_A_s] = 2.98;       m_max[ln10_10_A_s] = 3.20;
m_min[n_s] = 0.92;               m_max[n_s] = 1.04;

// Nuisance parameters (PLC TT & TTTEEE)
// Using thinner bounds from base_plikTTTEEE_lowTEB results
m_min[A_planck] = 0.9;           m_max[A_planck] = 1.1;

#ifndef LITE_HI_L
// Free -> fixed nuisance parameters (PLC TT) to own best-fit values
m_value[A_cib_217 - UP_TO_FREE_PARAMS] = 58.71;
m_value[xi_sz_cib - UP_TO_FREE_PARAMS] = 0.5222;
m_value[A_sz - UP_TO_FREE_PARAMS] = 5.114;
m_value[ps_A_100_100 - UP_TO_FREE_PARAMS] = 283.05;
m_value[ps_A_143_143 - UP_TO_FREE_PARAMS] = 48.38;
m_value[ps_A_143_217 - UP_TO_FREE_PARAMS] = 43.81;
m_value[ps_A_217_217 - UP_TO_FREE_PARAMS] = 105.6;
m_value[ksz_norm - UP_TO_FREE_PARAMS] = 2.378;
m_value[gal545_A_100 - UP_TO_FREE_PARAMS] = 7.2895;
m_value[gal545_A_143 - UP_TO_FREE_PARAMS] = 7.802;
m_value[gal545_A_143_217 - UP_TO_FREE_PARAMS] = 19.71;
m_value[gal545_A_217 - UP_TO_FREE_PARAMS] = 89.72;
m_value[calib_100T - UP_TO_FREE_PARAMS] = 0.07905;
m_value[calib_217T - UP_TO_FREE_PARAMS] = 0.04773;

// Free -> fixed nuisance parameters (PLC TTTEEE) to own best-fit values
m_value[galf_EE_A_100 - UP_TO_FREE_PARAMS] = 0.13276;
m_value[galf_EE_A_100_143 - UP_TO_FREE_PARAMS] = 0.10263;
m_value[galf_EE_A_100_217 - UP_TO_FREE_PARAMS] = 0.20231;
m_value[galf_EE_A_143 - UP_TO_FREE_PARAMS] = 0.65523;
m_value[galf_EE_A_143_217 - UP_TO_FREE_PARAMS] = 0.14404;
m_value[galf_EE_A_217 - UP_TO_FREE_PARAMS] = 0.15898;
m_value[galf_TE_A_100 - UP_TO_FREE_PARAMS] = 0.36603;
m_value[galf_TE_A_100_143 - UP_TO_FREE_PARAMS] = 0.21018;
m_value[galf_TE_A_100_217 - UP_TO_FREE_PARAMS] = 0.57801;
m_value[galf_TE_A_143 - UP_TO_FREE_PARAMS] = 1.79775;
m_value[galf_TE_A_143_217 - UP_TO_FREE_PARAMS] = 0.99885;
m_value[galf_TE_A_217 - UP_TO_FREE_PARAMS] = 0.99676;

// Fixed parameters (PLC TT & TTTEEE)
m_value[cib_index - UP_TO_FREE_PARAMS] = -1.3;
// Fixed parameters (PLC TTTEEE)
m_value[galf_EE_index - UP_TO_FREE_PARAMS] = -2.4;
m_value[galf_TE_index - UP_TO_FREE_PARAMS] = -2.4;
m_value[bleak_epsilon_0_0T_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_0T_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_0T_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_0T_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_0T_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_0T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_0T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_0T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_0T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_0T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_0T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_0T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_0T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_0T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_0T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_1T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_1T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_1T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_1T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_1T_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_1T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_1T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_1T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_1T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_1T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_2T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_2T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_2T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_2T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_2T_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_0E_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_0E_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_0E_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_0E_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_0E_0E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_0E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_0E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_0E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_0E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_0E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_0E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_0E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_0E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_0E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_0E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_1E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_1E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_1E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_1E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_1E_1E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_1E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_1E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_1E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_1E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_1E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_0_2E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_1_2E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_2_2E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_3_2E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[bleak_epsilon_4_2E_2E - UP_TO_FREE_PARAMS] = 0.0;
m_value[calib_100P - UP_TO_FREE_PARAMS] = 1.0;
m_value[calib_143P - UP_TO_FREE_PARAMS] = 1.0;
m_value[calib_217P - UP_TO_FREE_PARAMS] = 1.0;
m_value[A_pol - UP_TO_FREE_PARAMS] = 1.0;
#endif

#endif // FLAT_PRIOR_SET

#else
#error incorrect parameter set included
#endif // PARAM_SET == 10
