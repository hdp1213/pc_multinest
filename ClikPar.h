#ifndef CLIKPAR_H
#define CLIKPAR_H

#include "ClassEngine.hh"

class ClikPar {

public:
  ClikPar();
  ~ClikPar();

  // This enum follows Planck ordering (to a degree)
  enum param_t
  {
    // Free parameters (CLASS)
    pbh_frac,
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
    FREE_PARAMS,
    // Fixed parameters (CLASS)
    omega_b = FREE_PARAMS,
    omega_cdm,
    hundredxtheta_s,
    tau_reio,
    ln10_10_A_s,
    n_s,
    // Fixed parameters (PLC TT nuisance)
    cib_index,
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

  // CLASS functions
  void initialise_CLASS(int max_l, struct pbh_external* pbh_info);
  void scale_Cube(double* Cube);
  void set_derived_params(double* Cube);
  ClassEngine* get_CLASS();

  // Likelihood functions
#ifdef BAO_LIKE
  double calculate_BAO_likelihood() const; // uses CLASS
#endif
  double calculate_extra_priors(double* Cube) const;
  

private:
  bool m_has_gaussian_prior[TOTAL_PARAMS];
  bool m_is_log10[TOTAL_PARAMS];
  double m_min[TOTAL_PARAMS], m_max[TOTAL_PARAMS];
  double m_mean[TOTAL_PARAMS], m_stddev[TOTAL_PARAMS];

  // Array of transforming functions
  typedef double (*trans_t)(double);
  trans_t m_transform[TOTAL_PARAMS];

  // BAO variables
#ifdef BAO_LIKE
  double sixDF_z, sixDF_mean, sixDF_stddev;
  double LOWZ_z, LOWZ_mean, LOWZ_stddev;
  double CMASS_z, CMASS_mean, CMASS_stddev;
  double MGS_z, MGS_mean, MGS_stddev;
  double BAO_FUDGE;
#endif

  ClassEngine* m_class_engine;
};

#endif
