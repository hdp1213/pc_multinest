#ifndef PC_INIT_PARAMS
#define PC_INIT_PARAMS

double m_min[FREE_PARAM_AMT], m_max[FREE_PARAM_AMT];
double m_value[FIXED_PARAM_AMT];

bool m_has_gaussian_prior[TOTAL_PARAM_AMT];
double m_mean[TOTAL_PARAM_AMT], m_stddev[TOTAL_PARAM_AMT];

void initialise_params() {
// #include "TTTEEE+lowP_pbh_fixedLCDM-flat.cc"
#include "TTTEEE+lowP-flat.cc"
// #include "TTTEEE+lowP_pbh_fixedLCDM-gauss.cc"
#include "TTTEEE+lowP-gauss.cc"
}

#endif
