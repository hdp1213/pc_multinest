#ifndef CLIKPAR_H
#define CLIKPAR_H

#include "ClassEngine.hh"
#include <vector>


class ClikPar {

public:
  ClikPar();
  ~ClikPar();

// Include commands for various parameter sets
// #include "TTTEEE+lowP_pbh_fixedLCDM.h"
// #include "TTTEEE+lowP_pbh.h"
#include "TTTEEE+lowP.h"

  // CLASS functions
  void initialise_CLASS(int max_l, struct pbh_external* pbh_info);
  void scale_free_params(double* free_params);
  void set_derived_params(double* all_params);

  void run_CLASS(double* free_params);
  void get_CLASS_spectra(std::vector<unsigned>& cl_ls,
      std::vector<double>& cl_tt,
      std::vector<double>& cl_te,
      std::vector<double>& cl_ee,
      std::vector<double>& cl_bb);

  // Likelihood functions
#ifdef BAO_LIKE
  double calculate_BAO_likelihood() const; // uses CLASS
#endif
  double calculate_extra_likelihoods(double* in_params) const;

  double calculate_flat_prior() const;

  // Not const because non-const m_min and m_max
  double* get_lower_bounds();
  double* get_upper_bounds();

  double* get_fixed_params();

private:
  bool m_is_log10[TOTAL_PARAM_AMT];

  bool m_has_gaussian_prior[TOTAL_PARAM_AMT];
  double m_mean[TOTAL_PARAM_AMT], m_stddev[TOTAL_PARAM_AMT];

  // Only needed for free parameters
  double m_min[FREE_PARAM_AMT], m_max[FREE_PARAM_AMT];

  // Only needed for fixed parameters
  double m_value[FIXED_PARAM_AMT];

  // Array of transforming functions
  typedef double (*trans_t)(double);
  trans_t m_transform[TOTAL_PARAM_AMT];

  // BAO variables
#ifdef BAO_LIKE
  double sixDF_z, sixDF_mean, sixDF_stddev;
  double LOWZ_z, LOWZ_mean, LOWZ_stddev;
  double CMASS_z, CMASS_mean, CMASS_stddev;
  double MGS_z, MGS_mean, MGS_stddev;
  double BAO_FUDGE;
#endif

  ClassEngine* m_class_engine;

  double get_par(param_t param, double* free_params) const;
};

#endif
