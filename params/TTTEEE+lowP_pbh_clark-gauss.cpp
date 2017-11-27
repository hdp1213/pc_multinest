/****************************/
/*  GAUSSIAN PRIOR SETTING  */
/****************************/

#if PARAM_SET == 4

#ifndef GAUSS_PRIOR_SET
#define GAUSS_PRIOR_SET

#ifdef GAUSS_TAU
m_has_gaussian_prior[tau_reio] = true;
#endif

m_has_gaussian_prior[A_planck] = true;
#ifndef LITE_HI_L
// Nuisance parameters (PLC TT & TTTEEE)
m_has_gaussian_prior[gal545_A_100] = true;
m_has_gaussian_prior[gal545_A_143] = true;
m_has_gaussian_prior[gal545_A_143_217] = true;
m_has_gaussian_prior[gal545_A_217] = true;
m_has_gaussian_prior[calib_100T] = true;
m_has_gaussian_prior[calib_217T] = true;

// Nuisance parameters (PLC TTTEEE)
m_has_gaussian_prior[galf_EE_A_100] = true;
m_has_gaussian_prior[galf_EE_A_100_143] = true;
m_has_gaussian_prior[galf_EE_A_100_217] = true;
m_has_gaussian_prior[galf_EE_A_143] = true;
m_has_gaussian_prior[galf_EE_A_143_217] = true;
m_has_gaussian_prior[galf_EE_A_217] = true;
m_has_gaussian_prior[galf_TE_A_100] = true;
m_has_gaussian_prior[galf_TE_A_100_143] = true;
m_has_gaussian_prior[galf_TE_A_100_217] = true;
m_has_gaussian_prior[galf_TE_A_143] = true;
m_has_gaussian_prior[galf_TE_A_143_217] = true;
m_has_gaussian_prior[galf_TE_A_217] = true;
#endif

// Set Gaussian priors
#ifdef GAUSS_TAU
m_mean[tau_reio] = 0.07;
m_stddev[tau_reio] = 0.02;
#endif

m_mean[A_planck] = 1.0;
m_stddev[A_planck] = 0.0025;
#ifndef LITE_HI_L
// Nuisance parameters (PLC TT & TTTEEE)
m_mean[gal545_A_100] = 7.0;
m_stddev[gal545_A_100] = 2.0;
m_mean[gal545_A_143] = 9.0;
m_stddev[gal545_A_143] = 2.0;
m_mean[gal545_A_143_217] = 21.0;
m_stddev[gal545_A_143_217] = 8.5;
m_mean[gal545_A_217] = 80.0;
m_stddev[gal545_A_217] = 20.0;
m_mean[calib_100T] = 0.9990004;
m_stddev[calib_100T] = 0.001;
m_mean[calib_217T] = 0.99501;
m_stddev[calib_217T] = 0.002;

// Nuisance parameters (PLC TTTEEE)
m_mean[galf_EE_A_100] = 0.060;
m_stddev[galf_EE_A_100] = 0.012;
m_mean[galf_EE_A_100_143] = 0.050;
m_stddev[galf_EE_A_100_143] = 0.015;
m_mean[galf_EE_A_100_217] = 0.110;
m_stddev[galf_EE_A_100_217] = 0.033;
m_mean[galf_EE_A_143] = 0.10;
m_stddev[galf_EE_A_143] = 0.02;
m_mean[galf_EE_A_143_217] = 0.240;
m_stddev[galf_EE_A_143_217] = 0.048;
m_mean[galf_EE_A_217] = 0.72;
m_stddev[galf_EE_A_217] = 0.14;
m_mean[galf_TE_A_100] = 0.140;
m_stddev[galf_TE_A_100] = 0.042;
m_mean[galf_TE_A_100_143] = 0.120;
m_stddev[galf_TE_A_100_143] = 0.036;
m_mean[galf_TE_A_100_217] = 0.30;
m_stddev[galf_TE_A_100_217] = 0.09;
m_mean[galf_TE_A_143] = 0.240;
m_stddev[galf_TE_A_143] = 0.072;
m_mean[galf_TE_A_143_217] = 0.60;
m_stddev[galf_TE_A_143_217] = 0.18;
m_mean[galf_TE_A_217] = 1.80;
m_stddev[galf_TE_A_217] = 0.54;
#endif

#endif // GAUSS_PRIOR_SET

#else
#error incorrect parameter set included
#endif // PARAM_SET == 4
