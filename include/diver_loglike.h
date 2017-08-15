#ifndef DIVER_LOGLIKE
#define DIVER_LOGLIKE

#include "loglike.h"

#include <algorithm> // for std::copy
#include <limits>
#include <vector>

struct diver_bundle {
  std::vector<clik_struct*> clik_objs;
  ClassEngine* engine;
  std::vector<unsigned> cl_ls;
};

// Function to be minimized.  Corresponds to -ln(Likelihood).  Redirects to the target of context pointer.
// params contains only free and derived parameters
// at the end of likelihood computation, the latter part of the array should have derived parameters set as appropriate
double diver_loglike(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector, void*& context) {

  double loglike = std::numeric_limits<double>::max();

  // If params[] is within the parameter bounds, evaluate its likelihood
  if (validvector) {
    diver_bundle* plc_pack;
    double in_vals[UP_TO_FIXED_PARAMS];

    plc_pack = static_cast<diver_bundle*>(context);

    // Copy free and fixed parameters into in_vals, the input parameters
    std::copy(params, params + FREE_PARAM_AMT, in_vals);
    std::copy(m_value, m_value + FIXED_PARAM_AMT, in_vals + FREE_PARAM_AMT);

    loglike = pc_loglike(plc_pack->clik_objs, plc_pack->engine,
                         plc_pack->cl_ls, params, in_vals);
    fcall++;
  }

  quit = false;
  
  std::cout << "[diver_loglike] Calculated loglike of "
            << -loglike
            << std::endl;

  return -loglike;
}

double diver_prior(const double real_params[], const int real_param_dim, void*& context) {
  double res = 1.0;

  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    res *= (m_max[param] - m_min[param]);
  }

  return 1.0/res;
}


#endif
