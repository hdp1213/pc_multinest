/**********************************/
/*  UNIFORM (FLAT) PRIOR SETTING  */
/**********************************/

#if PARAM_SET == 4

#ifndef FLAT_PRIOR_SET
#define FLAT_PRIOR_SET

// Non-standard LCDM parameters
m_is_log10[pbh_frac] = true;
m_min[pbh_frac] = -8.0;          m_max[pbh_frac] = 0.0;
// m_min[pbh_frac] = 0.0;           m_max[pbh_frac] = 1.E0;
m_is_log10[pbh_mass] = true;
m_min[pbh_mass] = 5;             m_max[pbh_mass] = 7;
// m_min[pbh_mass] = 1.E5;          m_max[pbh_mass] = 1.E7;

// Nuisance parameters (PLC TT & TTTEEE)
// Using thinner bounds from base_plikTTTEEE_lowTEB results
m_value[A_planck - UP_TO_FREE_PARAMS] = 1.00013;
#ifndef LITE_HI_L
m_min[A_cib_217] = 30.0;         m_max[A_cib_217] = 100.0;
m_min[xi_sz_cib] = 0.0;          m_max[xi_sz_cib] = 1.0;
m_min[A_sz] = 0.0;               m_max[A_sz] = 10.0;
m_min[ps_A_100_100] = 140.0;     m_max[ps_A_100_100] = 370.0;
m_min[ps_A_143_143] = 10.0;      m_max[ps_A_143_143] = 80.0;
m_min[ps_A_143_217] = 0.0;       m_max[ps_A_143_217] = 80.0;
m_min[ps_A_217_217] = 50.0;      m_max[ps_A_217_217] = 140.0;
m_min[ksz_norm] = 0.0;           m_max[ksz_norm] = 10.0;
m_min[gal545_A_100] = 0.0;       m_max[gal545_A_100] = 20.0;
m_min[gal545_A_143] = 0.0;       m_max[gal545_A_143] = 20.0;
m_min[gal545_A_143_217] = 0.0;   m_max[gal545_A_143_217] = 40.0;
m_min[gal545_A_217] = 50.0;      m_max[gal545_A_217] = 110.0;
m_min[calib_100T] = 0.0;         m_max[calib_100T] = 2.0;
m_min[calib_217T] = 0.0;         m_max[calib_217T] = 2.0;

// Nuisance parameters (PLC TTTEEE)
// Using approximate bounds on base_plikTTTEEE_lowTEB results
m_min[galf_EE_A_100] = 0.04;
m_max[galf_EE_A_100] = 0.12;
m_min[galf_EE_A_100_143] = 0.01;
m_max[galf_EE_A_100_143] = 0.09;
m_min[galf_EE_A_100_217] = 0.00;
m_max[galf_EE_A_100_217] = 0.26;
m_min[galf_EE_A_143] = 0.060;
m_max[galf_EE_A_143] = 0.150;
m_min[galf_EE_A_143_217] = 0.0;
m_max[galf_EE_A_143_217] = 0.6;
m_min[galf_EE_A_217] = 0.0;
m_max[galf_EE_A_217] = 1.5;
m_min[galf_TE_A_100] = 0.0;
m_max[galf_TE_A_100] = 0.4;
m_min[galf_TE_A_100_143] = 0.0;
m_max[galf_TE_A_100_143] = 0.3;
m_min[galf_TE_A_100_217] = 0.0;
m_max[galf_TE_A_100_217] = 0.8;
m_min[galf_TE_A_143] = 0.0;
m_max[galf_TE_A_143] = 0.5;
m_min[galf_TE_A_143_217] = 0.0;
m_max[galf_TE_A_143_217] = 1.0;
m_min[galf_TE_A_217] = 0.0;
m_max[galf_TE_A_217] = 4.5;

// Fixed LCDM parameters to TTTEEE+lowP best-fit values
m_value[omega_b - UP_TO_FREE_PARAMS] = 0.022250;
m_value[omega_cdm - UP_TO_FREE_PARAMS] = 0.11978;
m_value[hundredxtheta_s - UP_TO_FREE_PARAMS] = 1.041706;
m_value[tau_reio - UP_TO_FREE_PARAMS] = 0.0781;
m_value[ln10_10_A_s - UP_TO_FREE_PARAMS] = 3.0906;
m_value[n_s - UP_TO_FREE_PARAMS] = 0.96417;


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
#endif // PARAM_SET == 4
