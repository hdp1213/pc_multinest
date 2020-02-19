#include "multinest_loglike.hpp"
#include "loglike.hpp"

#include <iostream>

using namespace std;

void
multinest_loglike(double* Cube, int& ndim, int& npars, double& lnew, void* context) {
  // Transform values in Cube to physical space
  trans_t func;

  for (unsigned param = 0; param < ndim; ++param) {
    func = g_transform[param];
    Cube[param] = func(g_min[param] + Cube[param] * (g_max[param] - g_min[param]));
  }

  // Calculate log likelihood
  lnew = calculate_full_likelihood(Cube, static_cast<likelihood_context*>(context));

  cout << "[multinest_loglike] Calculated loglike of " << lnew << endl;
}
