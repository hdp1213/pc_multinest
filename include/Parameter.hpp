#ifndef PARAMETER_H
#define PARAMETER_H

#include <assert.h>
#include <forward_list>

#include "params/tau_degen_fixed_lcdm_nuisance.hpp"

// Transformation function typedef
typedef double (*trans_t)(double);
inline double identity(double x) {return x;};

// Free flat priors
extern double g_min[FREE_PARAM_AMT];
extern double g_max[FREE_PARAM_AMT];

// Fixed values
extern double g_value[FIXED_PARAM_AMT];

// Transforms
extern trans_t g_transform[FREE_PARAM_AMT];

// Struct for storing Gaussian priors
typedef struct {
  param_t param;
  double mean;
  double stddev;
} gauss_t;

// Gaussian priors
extern std::forward_list<gauss_t> g_gauss;

// Actual functions we use
void initialise_param_arrays();

inline void add_fixed_param(const param_t param, const double value) {
  assert(param >= UP_TO_FREE_PARAMS);
  g_value[param - UP_TO_FREE_PARAMS] = value;
}

inline void add_free_param(const param_t param, const double min, const double max) {
  assert(min < max);
  g_min[param] = min; g_max[param] = max;
}

inline void add_gauss_prior(const param_t param, const double mean, const double stddev) {
  assert(stddev > 0.0);
  gauss_t gauss_vars = {.param = param, .mean = mean, .stddev = stddev};
  g_gauss.push_front(gauss_vars);
}

inline void add_transform(const param_t param, const trans_t trans) {
  g_transform[param] = trans;
}

inline double get_value(param_t param, double free_params[]) {
  return (param < UP_TO_FREE_PARAMS) ? free_params[param] : g_value[param - UP_TO_FREE_PARAMS];
}

#endif
