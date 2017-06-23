#ifndef CLIKPAR_H
#define CLIKPAR_H

#include "ClassEngine.hh"


class ClikPar {

public:
  ClikPar();
  ~ClikPar();

// Include commands for various parameter sets
#include "TTTEEE+lowP_pbh_fixedLCDM.h"
// #include "TTTEEE+lowP_pbh.h"
// #include "TTTEEE+lowP.h"

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

  double* get_lower_bounds() const;
  double* get_upper_bounds() const;


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
