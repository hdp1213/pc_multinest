#ifndef PC_MULTINEST
#define PC_MULTINEST

#include <string>

struct multinest_settings {
  bool IS;           // do Nested Importance Sampling?
  bool mmodal;       // do mode separation?
  bool ceff;         // run in constant efficiency mode?
                    // Bad for evidence calculations

  int nlive;        // number of live points
  double efr;       // set the required efficiency. 0.8
                    // for parameter estimation, 0.3 for
                    // evidence
  double tol;       // tol, defines the stopping criteria
                    // 0.5 should give enough accuracy

  int ndims;        // dimensionality (no. of free
                    // parameters)
  int nPar;         // total no. of parameters including
                    // free & derived parameters
  int nClsPar;      // no. of parameters to do mode
                    // separation on

  int updInt;       // after how many iterations feedback
                    // is required & the output files
                    // should be updated
                    // note: posterior files are updated
                    // & dumper routine is called after
                    // every updInt*10 iterations
  double Ztol;      // all the modes with logZ < Ztol are
                    // ignored
  int maxModes;     // expected max no. of modes (used
                    // only for memory allocation)
  int* pWrap;       // which parameters to have periodic
                    // boundary conditions? dim of ndims

  std::string root; // root for output files
  int seed;         // random no. generator seed, if < 0
                    // then take the seed from system
                    // clock
  bool fb;          // need feedback on standard output?
  bool resume;      // resume from a previous job?
  bool outfile;     // write output files?
  bool initMPI;     // initialize MPI routines?, relevant
                    // only if compiling with MPI
                    // set it to F if you want your main
                    // program to handle MPI
                    // initialization
  double logZero;   // points with loglike < logZero will
                    // be ignored by MultiNest
  int maxiter;      // max no. of iterations, a
                    // non-positive value means infinity.
                    // MultiNest will terminate if
                    // either it has done max no. of
                    // iterations or convergence
                    // criterion (defined through tol)
                    // has been satisfied
};

// pc_multinest: a pc_loglike wrapper for nested::run. What is nested::run, you ask? This behemoth:
// run(bool IS, bool mmodal, bool ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar, int maxModes,int updInt, double Ztol, const std::string & root, int seed, int *pWrap, bool fb, bool resume, bool outfile, bool initMPI, double logZero, int maxiter, void (*LogLike)(double *Cube, int &n_dim, int &n_par, double &lnew, void *),void (*dumper)(int &, int &, int &, double **, double **, double **, double &, double &, double &, double &, void *), void *context);

void pc_multinest(void (*loglike)(double*,
                                  int&,
                                  int&,
                                  double&,
                                  void*),
                  void (*dumper)(int&,
                                 int&,
                                 int&,
                                 double**,
                                 double**,
                                 double**,
                                 double&,
                                 double&,
                                 double&,
                                 double&,
                                 void*),
                  multinest_settings& s,
                  void*& context);

#endif
