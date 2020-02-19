#include "Parameter.hpp"

#include <cmath>

using namespace std;

double g_value[FIXED_PARAM_AMT];
double g_min[FREE_PARAM_AMT];
double g_max[FREE_PARAM_AMT];

forward_list<gauss_t> g_gauss;

trans_t g_transform[FREE_PARAM_AMT];

double pow10(double x) {return pow(10, x);};

void
initialise_param_arrays() {
  // Initialise all transforms to the identity function
  for (unsigned param = 0; param < FREE_PARAM_AMT; ++param) {
    add_transform((param_t) param, identity);
  }

  #include "params/tau_degen_fixed_lcdm_nuisance.cpp"
}
