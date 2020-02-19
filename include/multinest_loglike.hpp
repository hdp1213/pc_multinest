#ifndef MULTINEST_LOGLIKE_H
#define MULTINEST_LOGLIKE_H

void multinest_loglike(double* Cube, int& ndim, int& npars, double& lnew, void* context);

inline void multinest_dumper(int& nSamples, int& nlive, int& nPar, double** physLive, double** posterior, double** paramConstr, double& maxLogLike, double& logZ, double& INSlogZ, double& logZerr, void* context) {};

#endif
