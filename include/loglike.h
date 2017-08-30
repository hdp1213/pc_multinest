#ifndef PC_LOGLIKE
#define PC_LOGLIKE

#include "init_plc.h"

#include <vector>

#ifndef PARAM_SET
#error no parameter set included!
#endif

/* "Actual" constants */
const double MIN_LOGLIKE = -1E90;
const double ZERO_CL = 0.0;
const int CLASS_CL_AMT = 4;

/* Misc. gaussian likelihood values */
const double SZ_MEAN = 9.5;
const double SZ_STDDEV = 3.0;

/* Constant locations of derived parameters in raw_params */
const int H0_LOC = static_cast<int>(H0 - FIXED_PARAM_AMT);
const int Omega_b_LOC = static_cast<int>(Omega_b - FIXED_PARAM_AMT);
const int Omega_cdm_LOC = static_cast<int>(Omega_cdm - FIXED_PARAM_AMT);
const int Omega_L_LOC = static_cast<int>(Omega_L - FIXED_PARAM_AMT);
const int Omega_g_LOC = static_cast<int>(Omega_g - FIXED_PARAM_AMT);
const int sigma8_LOC = static_cast<int>(sigma8 - FIXED_PARAM_AMT);
const int age_LOC = static_cast<int>(age - FIXED_PARAM_AMT);
const int conf_age_LOC = static_cast<int>(conf_age - FIXED_PARAM_AMT);
const int z_drag_LOC = static_cast<int>(z_drag - FIXED_PARAM_AMT);
const int rs_drag_LOC = static_cast<int>(rs_drag - FIXED_PARAM_AMT);

extern double* m_mean;
extern double* m_has_gaussian_prior;
extern double* m_stddev;

/**
 * What do we need to compute the log likelihood for plc_class?
 * - the initialised .clik objects
 * - a CLASS engine
 * - the actual parameter values being considered for the particular
 *   point in parameter space. This contains just free parameters though
 * - the fixed values, which are always constant and never change
 * - a way of telling which parameter is free or fixed
 * - any gaussian priors!!!
 * 
 * Things that have already happened prior to this:
 * - .clik objects have been initialised
 * - CLASS engine has been initialised
 *
 * I feel like we need a clik struct that contains all of the
 * useful information already extracted.
**/
double pc_loglike(std::vector<clik_struct*>& clik_objs,
                  ClassEngine*& class_engine,
                  std::vector<unsigned>& cl_ls,
                  double raw_params[],
                  double in_vals[]);

double calculate_PLC_likelihood(clik_struct*& clik_obj,
                                double in_vals[],
                                std::vector<std::vector<double> >& class_spectra);
// implicitly uses global arrays
double calculate_extra_likelihoods(double in_vals[],
                                   double raw_params[]);
// raw_params[] contains only free and derived parameters
void set_derived_params(double raw_params[],
                        ClassEngine*& class_engine);

#endif
