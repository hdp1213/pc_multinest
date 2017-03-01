#ifndef CLIKPAR_H
#define CLIKPAR_H

#include "ClassEngine.hh"

class ClikPar {

public:
  ClikPar(int free_param_amt, int fixed_param_amt);
  ~ClikPar();

  // This enum follows Planck ordering (to a degree)
  enum param_t
  {
    // Free parameters (CLASS)
    omega_b,
    omega_cdm,
    hundredxtheta_s,
    tau_reio,
    ln10_10_A_s,
    n_s,
    // Free parameters (PLC)
    A_planck,
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
    calib_100T,
    calib_217T,
    // Fixed parameters (PLC)
    cib_index, // = -1.3
    // Derived parameters (LCDM)
    H0,
    Omega_b,
    Omega_cdm,
    Omega_L,
    Omega_g,
    // Omega_k, // an input parameter!
    sigma8,
    age,
    conf_age,
    // K, // presumably derived from Omega_k
    z_drag,
    rs_drag,
    // Final total number of parameters
    TOTAL_PARAMS
  };

  // CLASS functions
  void initialise_CLASS(int max_l);
  void scale_Cube(double* Cube);
  void set_derived_params(double* Cube);
  ClassEngine* get_CLASS();

  // Likelihood functions
  double calculate_BAO_likelihood() const; // uses CLASS
  double calculate_extra_priors(double* Cube) const;
  

private:
  bool m_is_gaussian[TOTAL_PARAMS];
  double m_min[TOTAL_PARAMS], m_max[TOTAL_PARAMS];
  double m_mean[TOTAL_PARAMS], m_stddev[TOTAL_PARAMS];
  int m_free_param_amt, m_gaussian_param_amt;
  int m_fixed_param_amt;

  // BAO variables
  double sixDF_z, sixDF_mean, sixDF_stddev;
  double LOWZ_z, LOWZ_mean, LOWZ_stddev;
  double CMASS_z, CMASS_mean, CMASS_stddev;
  double MGS_z, MGS_mean, MGS_stddev;
  double BAO_FUDGE;

  ClassEngine* m_class_engine;
};

#endif