#ifndef MULTINEST_LOGLIKE
#define MULTINEST_LOGLIKE

#include "loglike.h"

#include <algorithm> // for std::copy

// Function to be minimized.  Corresponds to ln(Likelihood).  Redirects to the target of context pointer.
// Cube contains only free and derived parameters. To begin, free params
//  are in [0,1] range
// at the end of likelihood computation, the latter part of the array has derived parameters set as appropriate
void multinest_loglike(double *Cube, int &ndim, int &npars, double &lnew, void *context) {
  plc_bundle* plc_pack;
  double in_vals[UP_TO_FIXED_PARAMS];

  plc_pack = static_cast<plc_bundle*>(context);

  // Transform free Cube parameters from unit hypercube to physical space
  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    Cube[param] = m_min[param] + Cube[param] * (m_max[param] - m_min[param]);
  }

  // Copy free and fixed parameters into in_vals, the input parameters
  std::copy(Cube, Cube + FREE_PARAM_AMT, in_vals);
  std::copy(m_value, m_value + FIXED_PARAM_AMT, in_vals + FREE_PARAM_AMT);

  lnew = pc_loglike(plc_pack->clik_objs, plc_pack->engine,
                    plc_pack->cl_ls, Cube, in_vals);

  std::cout << "[multinest_loglike] Calculated loglike of "
            << lnew
            << std::endl;
}

void multinest_dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context) {
  // blank
}

#endif
