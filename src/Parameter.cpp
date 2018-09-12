#include "Parameter.hpp"

using namespace std;

double g_value[FIXED_PARAM_AMT];
double g_min[FREE_PARAM_AMT];
double g_max[FREE_PARAM_AMT];

forward_list<gauss_t> g_gauss;

trans_t g_transform[FREE_PARAM_AMT];

void
initialise_param_arrays() {
  // Free CLASS parameter flat priors
  add_free_param(omega_b, 0.016, 0.028);
  add_free_param(omega_cdm, 0.108, 0.130);
  add_free_param(hundredxtheta_s, 1.039, 1.043);
  add_free_param(tau_reio, 0.01, 0.15);
  add_free_param(ln10_10_A_s, 2.98, 3.20);
  add_free_param(n_s, 0.92, 1.04);

  // Free PLC parameter flat priors
  add_free_param(A_planck, 0.9, 1.1);
#ifndef LITE_HI_L
  add_free_param(A_cib_217, 30.0, 100.0);
  add_free_param(xi_sz_cib, 0.0, 1.0);
  add_free_param(A_sz, 0.0, 10.0);
  add_free_param(ps_A_100_100, 140.0, 370.0);
  add_free_param(ps_A_143_143, 10.0, 80.0);
  add_free_param(ps_A_143_217, 0.0, 80.0);
  add_free_param(ps_A_217_217, 50.0, 140.0);
  add_free_param(ksz_norm, 0.0, 10.0);
  add_free_param(gal545_A_100, 0.0, 20.0);
  add_free_param(gal545_A_143, 0.0, 20.0);
  add_free_param(gal545_A_143_217, 0.0, 40.0);
  add_free_param(gal545_A_217, 50.0, 110.0);
  add_free_param(calib_100T, 0.0, 2.0);
  add_free_param(calib_217T, 0.0, 2.0);
  add_free_param(galf_EE_A_100, 0.04, 0.12);
  add_free_param(galf_EE_A_100_143, 0.01, 0.09);
  add_free_param(galf_EE_A_100_217, 0.00, 0.26);
  add_free_param(galf_EE_A_143, 0.060, 0.150);
  add_free_param(galf_EE_A_143_217, 0.0, 0.6);
  add_free_param(galf_EE_A_217, 0.0, 1.5);
  add_free_param(galf_TE_A_100, 0.0, 0.4);
  add_free_param(galf_TE_A_100_143, 0.0, 0.3);
  add_free_param(galf_TE_A_100_217, 0.0, 0.8);
  add_free_param(galf_TE_A_143, 0.0, 0.5);
  add_free_param(galf_TE_A_143_217, 0.0, 1.0);
  add_free_param(galf_TE_A_217, 0.0, 4.5);

  // Fixed parameters
  add_fixed_param(cib_index, -1.3);
  add_fixed_param(galf_EE_index, -2.4);
  add_fixed_param(galf_TE_index, -2.4);
  add_fixed_param(bleak_epsilon_0_0T_0E, 0.0);
  add_fixed_param(bleak_epsilon_1_0T_0E, 0.0);
  add_fixed_param(bleak_epsilon_2_0T_0E, 0.0);
  add_fixed_param(bleak_epsilon_3_0T_0E, 0.0);
  add_fixed_param(bleak_epsilon_4_0T_0E, 0.0);
  add_fixed_param(bleak_epsilon_0_0T_1E, 0.0);
  add_fixed_param(bleak_epsilon_1_0T_1E, 0.0);
  add_fixed_param(bleak_epsilon_2_0T_1E, 0.0);
  add_fixed_param(bleak_epsilon_3_0T_1E, 0.0);
  add_fixed_param(bleak_epsilon_4_0T_1E, 0.0);
  add_fixed_param(bleak_epsilon_0_0T_2E, 0.0);
  add_fixed_param(bleak_epsilon_1_0T_2E, 0.0);
  add_fixed_param(bleak_epsilon_2_0T_2E, 0.0);
  add_fixed_param(bleak_epsilon_3_0T_2E, 0.0);
  add_fixed_param(bleak_epsilon_4_0T_2E, 0.0);
  add_fixed_param(bleak_epsilon_0_1T_1E, 0.0);
  add_fixed_param(bleak_epsilon_1_1T_1E, 0.0);
  add_fixed_param(bleak_epsilon_2_1T_1E, 0.0);
  add_fixed_param(bleak_epsilon_3_1T_1E, 0.0);
  add_fixed_param(bleak_epsilon_4_1T_1E, 0.0);
  add_fixed_param(bleak_epsilon_0_1T_2E, 0.0);
  add_fixed_param(bleak_epsilon_1_1T_2E, 0.0);
  add_fixed_param(bleak_epsilon_2_1T_2E, 0.0);
  add_fixed_param(bleak_epsilon_3_1T_2E, 0.0);
  add_fixed_param(bleak_epsilon_4_1T_2E, 0.0);
  add_fixed_param(bleak_epsilon_0_2T_2E, 0.0);
  add_fixed_param(bleak_epsilon_1_2T_2E, 0.0);
  add_fixed_param(bleak_epsilon_2_2T_2E, 0.0);
  add_fixed_param(bleak_epsilon_3_2T_2E, 0.0);
  add_fixed_param(bleak_epsilon_4_2T_2E, 0.0);
  add_fixed_param(bleak_epsilon_0_0E_0E, 0.0);
  add_fixed_param(bleak_epsilon_1_0E_0E, 0.0);
  add_fixed_param(bleak_epsilon_2_0E_0E, 0.0);
  add_fixed_param(bleak_epsilon_3_0E_0E, 0.0);
  add_fixed_param(bleak_epsilon_4_0E_0E, 0.0);
  add_fixed_param(bleak_epsilon_0_0E_1E, 0.0);
  add_fixed_param(bleak_epsilon_1_0E_1E, 0.0);
  add_fixed_param(bleak_epsilon_2_0E_1E, 0.0);
  add_fixed_param(bleak_epsilon_3_0E_1E, 0.0);
  add_fixed_param(bleak_epsilon_4_0E_1E, 0.0);
  add_fixed_param(bleak_epsilon_0_0E_2E, 0.0);
  add_fixed_param(bleak_epsilon_1_0E_2E, 0.0);
  add_fixed_param(bleak_epsilon_2_0E_2E, 0.0);
  add_fixed_param(bleak_epsilon_3_0E_2E, 0.0);
  add_fixed_param(bleak_epsilon_4_0E_2E, 0.0);
  add_fixed_param(bleak_epsilon_0_1E_1E, 0.0);
  add_fixed_param(bleak_epsilon_1_1E_1E, 0.0);
  add_fixed_param(bleak_epsilon_2_1E_1E, 0.0);
  add_fixed_param(bleak_epsilon_3_1E_1E, 0.0);
  add_fixed_param(bleak_epsilon_4_1E_1E, 0.0);
  add_fixed_param(bleak_epsilon_0_1E_2E, 0.0);
  add_fixed_param(bleak_epsilon_1_1E_2E, 0.0);
  add_fixed_param(bleak_epsilon_2_1E_2E, 0.0);
  add_fixed_param(bleak_epsilon_3_1E_2E, 0.0);
  add_fixed_param(bleak_epsilon_4_1E_2E, 0.0);
  add_fixed_param(bleak_epsilon_0_2E_2E, 0.0);
  add_fixed_param(bleak_epsilon_1_2E_2E, 0.0);
  add_fixed_param(bleak_epsilon_2_2E_2E, 0.0);
  add_fixed_param(bleak_epsilon_3_2E_2E, 0.0);
  add_fixed_param(bleak_epsilon_4_2E_2E, 0.0);
  add_fixed_param(calib_100P, 1.0);
  add_fixed_param(calib_143P, 1.0);
  add_fixed_param(calib_217P, 1.0);
  add_fixed_param(A_pol, 1.0);
#endif  // LITE_HI_L

  // Gaussian priors
#ifdef GAUSS_TAU
  add_gauss_prior(tau_reio, 0.07, 0.02);
#endif
#ifndef LITE_HI_L
  // Nuisance parameters (PLC TT & TTTEEE)
  add_gauss_prior(gal545_A_100, 7.0, 2.0);
  add_gauss_prior(gal545_A_143, 9.0, 2.0);
  add_gauss_prior(gal545_A_143_217, 21.0, 8.5);
  add_gauss_prior(gal545_A_217, 80.0, 20.0);
  add_gauss_prior(calib_100T, 0.9990004, 0.001);
  add_gauss_prior(calib_217T, 0.99501, 0.002);

  // Nuisance parameters (PLC TTTEEE)
  add_gauss_prior(galf_EE_A_100, 0.060, 0.012);
  add_gauss_prior(galf_EE_A_100_143, 0.050, 0.015);
  add_gauss_prior(galf_EE_A_100_217, 0.110, 0.033);
  add_gauss_prior(galf_EE_A_143, 0.10, 0.02);
  add_gauss_prior(galf_EE_A_143_217, 0.240, 0.048);
  add_gauss_prior(galf_EE_A_217, 0.72, 0.14);
  add_gauss_prior(galf_TE_A_100, 0.140, 0.042);
  add_gauss_prior(galf_TE_A_100_143, 0.120, 0.036);
  add_gauss_prior(galf_TE_A_100_217, 0.30, 0.09);
  add_gauss_prior(galf_TE_A_143, 0.240, 0.072);
  add_gauss_prior(galf_TE_A_143_217, 0.60, 0.18);
  add_gauss_prior(galf_TE_A_217, 1.80, 0.54);
#endif
  add_gauss_prior(A_planck, 1.0, 0.0025);

  // Initialise all transforms to the identity function
  for (unsigned param = 0; param < FREE_PARAM_AMT; ++param) {
    add_transform((param_t) param, identity);
  }

  // Add any non-identity parameter transforms...
}
