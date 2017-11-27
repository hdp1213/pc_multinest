#ifndef PC_DIVER
#define PC_DIVER

#include <string>

struct diver_settings {
  int     nPar;               // Dimensionality of the parameter space
                              // (free) (D in literature)
  double* lowerbounds;        // Lower boundaries of parameter space
  double* upperbounds;        // Upper boundaries of parameter space
  std::string  path;          // Path to save samples, resume files,
                              // etc
  int     nDerived;           // Number of derived quantities to output
  int     nDiscrete;          // Number of parameters that are to be
                              // treated as discrete
  int*    discrete;           // Indices of discrete parameters,
                              // Fortran style, i.e. starting at 1!!
  bool    partitionDiscrete;  // Split the population evenly amongst
                              // discrete parameters and evolve
                              // separately

  int     maxciv;             // Maximum number of civilisations
  int     maxgen;             // Maximum number of generations per
                              // civilisation
  int     NP;                 // Population size (individuals per
                              // generation)
  int     nF;                 // Size of the array indicating scale
                              // factors
  double F[1] = {0.6};        // Scale factor(s).  Note that this must
                              // be entered as an array.
  double  Cr;                 // Crossover factor
  double  lambda;             // Mixing factor between best and
                              // rand/current
  bool    current;            // Use current vector for mutation
  bool    expon;              // Use exponential crossover
  int     bndry;              // Boundary constraint: 1=brick wall,
                              // 2=random re-initialization,
                              // 3=reflection
  bool    jDE;                // Use self-adaptive choices for
                              // rand/1/bin parameters as per Brest
                              // et al 2006
  bool    lambdajDE;          // Use self-adaptive rand-to-best/1/bin
                              // parameters; based on Brest et al 2006

  double  convthresh;         // Threshold for gen-level convergence:
                              // smoothed fractional improvement in the
                              // mean population value
  int     convsteps;          // Number of steps to smooth over when
                              // checking convergence
  bool    removeDuplicates;   // Weed out duplicate vectors within a
                              // single generation

  bool    doBayesian;         // Calculate approximate log evidence and
                              // posterior weightings
  double  maxNodePop;         // Population at which node is
                              // partitioned in binary space
                              // partitioning for posterior
  double  Ztolerance;         // Input tolerance in log-evidence

  int     savecount;          // Save progress every savecount
                              // generations
  bool    resume;             // Restart from a previous run
  bool    outputSamples;      // Write output .raw and .sam (if
                              // nDerived != 0) files
  int     init_pop_strategy;  // Initialisation strategy: 0=one shot,
                              // 1=n-shot, 2=n-shot with error if no
                              // valid vectors found.
  int     max_init_attempts;  // Maximum number of times to try to find
                              // a valid vector for each slot in the
                              // initial population.
  double  max_acceptable_val; // Maximum fitness to accept for the
                              // initial generation if
                              // init_population_strategy > 0.
  int     verbose;            // Output verbosity: 0=only error
                              // messages, 1=basic info, 2=civ-level
                              // info, 3+=population info
};

// pc_diver: a pc_multinest wrapper for cdiver. What is cdiver, you ask? This behemoth:
// cdiver(double (*obj_func)(double[] params, const int param_dim, int* fcall, bool* quit, const bool validvector, void** context), int nPar, const double[] lowerbounds, const double[] upperbounds, const char[] path, int nDerived, int nDiscrete, const int[] discrete, bool partitionDiscrete, const int maxciv, const int maxgen, int NP, int nF, const double[] F, double Cr, double lambda, bool current, bool expon, int bndry, bool jDE, bool lambdajDE, double convthresh, int convsteps, bool removeDuplicates, bool doBayesian, double(*prior_func)(const double[] real_params, const int real_param_dim, void** context), double maxNodePop, double Ztolerance, int savecount, bool resume, bool outputSamples, int init_pop_strategy, int max_init_attempts, double max_acceptable_val, void** context, int verbose);
void pc_diver(double (*obj_func)(double[],
                                 const int,
                                 int&,
                                 bool&,
                                 const bool,
                                 void*&),
              double (*prior_func)(const double[],
                                   const int,
                                   void*&),
              diver_settings& s,
              void*& context);

#endif
