#ifndef MULTINEST_LOGLIKE
#define MULTINEST_LOGLIKE

#include <vector>

#include "loglike.h"

extern double m_min[];
extern double m_max[];
extern double m_value[];

extern trans_t m_transform[];

// Function to be minimized.  Corresponds to ln(Likelihood).  Redirects to the target of context pointer.
// Cube contains only free and derived parameters. To begin, free params
//  are in [0,1] range
// at the end of likelihood computation, the latter part of the array has derived parameters set as appropriate
void multinest_loglike(double *Cube, int &ndim, int &npars, double &lnew, void *context);

void multinest_loglike_pt(std::vector<double>& new_vals, int &ndim, int &npars, double &lnew, void *context);

void multinest_dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context);

#endif
