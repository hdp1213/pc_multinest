#include "PLCPack.h"
#include "ClikObject.h"
#include "ClikPar.h"
#include "diver.hpp"

#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <limits>

#include <iostream>

#include "loglike.h"

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

// Function to be minimized.  Corresponds to -ln(Likelihood).  Redirects to the target of context pointer.
double objective(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector, void*& context)
{
  double loglike = 0.0;

  // Not needed inside LogLike, so set to zero!
  int free_dim = 0;
  int param_dim_free = param_dim;

  // If params[] is within the parameter bounds, evaluate its likelihood
  if (validvector) {
    LogLike(params, free_dim, param_dim_free, loglike, context);
    fcall++;
  }

  return -loglike;
}

// pc_multinest prior function
double pc_prior(const double real_params[], const int real_param_dim, void*& context)
{
  PLCPack* plc_pack;
  plc_pack = static_cast<PLCPack*>(context);

  return plc_pack->calculate_prior();
}

// pc_diver: a pc_multinest wrapper for cdiver. What is cdiver, you ask? This behemoth:
// cdiver(double (*obj_func)(double[] params, const int param_dim, int* fcall, bool* quit, const bool validvector, void** context), int nPar, const double[] lowerbounds, const double[] upperbounds, const char[] path, int nDerived, int nDiscrete, const int[] discrete, bool partitionDiscrete, const int maxciv, const int maxgen, int NP, int nF, const double[] F, double Cr, double lambda, bool current, bool expon, int bndry, bool jDE, bool lambdajDE, double convthresh, int convsteps, bool removeDuplicates, bool doBayesian, double(*prior_func)(const double[] real_params, const int real_param_dim, void** context), double maxNodePop, double Ztolerance, int savecount, bool resume, bool outputSamples, int init_pop_strategy, int max_init_attempts, double max_acceptable_val, void** context, int verbose);

void pc_diver(double (*obj_func)(double[], const int, int&, bool&, const bool, void*&), double (*prior_func)(const double[], const int, void*&), diver_settings& s, void* context) {

  const int maxciv = s.maxciv;
  const int maxgen = s.maxgen;

  cdiver(obj_func,
         s.nPar,
         const_cast<double*>(s.lowerbounds),
         const_cast<double*>(s.upperbounds),
         s.path.c_str(),
         s.nDerived,
         s.nDiscrete,
         const_cast<int*>(s.discrete),
         s.partitionDiscrete,
         maxciv,
         maxgen,
         s.NP,
         s.nF,
         const_cast<double*>(s.F),
         s.Cr,
         s.lambda,
         s.current,
         s.expon,
         s.bndry,
         s.jDE,
         s.lambdajDE,
         s.convthresh,
         s.convsteps,
         s.removeDuplicates,
         s.doBayesian,
         prior_func,
         s.maxNodePop,
         s.Ztolerance,
         s.savecount,
         s.resume,
         s.outputSamples,
         s.init_pop_strategy,
         s.max_init_attempts,
         s.max_acceptable_val,
         context,
         s.verbose);

}

int main(int argc, char** argv)
{
  ////// Diver settings //////

  diver_settings settings;

  // settings.nPar               = 2;
  // settings.lowerbounds        = {-6.0, -6.0};
  // settings.upperbounds        = { 6.0,  6.0};
  // settings.path               = "output/pc_diver-";
  // settings.nDerived           = 0;
  settings.nDiscrete          = 0;
  settings.discrete           = {none};
  settings.partitionDiscrete  = false;

  settings.maxciv             = 100;
  settings.maxgen             = 100;
  // settings.NP                 = 200;
  settings.nF                 = 1;
  // settings.F                  = {0.6};
  settings.Cr                 = 0.9;
  settings.lambda             = 0.8;
  settings.current            = false;
  settings.expon              = false;
  settings.bndry              = 3;
  settings.jDE                = true;
  settings.lambdajDE          = true;

  settings.convthresh         = 1.e-3;
  settings.convsteps          = 10;
  settings.removeDuplicates   = true;

  settings.doBayesian         = true;
  settings.maxNodePop         = 1.9;
  settings.Ztolerance         = 0.1;

  settings.savecount          = 1;
  settings.resume             = false; // watch out for me please!
  settings.outputSamples      = true;
  settings.init_pop_strategy  = 2; // fatal n-shot
  settings.max_init_attempts  = 1000;
  settings.max_acceptable_val = 1e7;
  settings.verbose            = 1;


  ////// pc_mulinest variables //////

  void* context = 0;

  // High l full likelihood variables
  char hi_l_clik_path[255];
  strcpy(hi_l_clik_path, PLIK_HI_L_FILE_DIR);
#ifdef LITE_HI_L
  strcat(hi_l_clik_path, "/plik_lite_v18_TTTEEE.clik/");
#else
  strcat(hi_l_clik_path, "/plik_dx11dr2_HM_v18_TTTEEE.clik/");
#endif
  ClikObject* hi_l_clik(0);
  std::vector<ClikPar::param_t> hi_l_nuis_enums;

  // Low l likelihood variables  
  char lo_l_clik_path[255];
  strcpy(lo_l_clik_path, PLIK_LOW_L_FILE_DIR);
  strcat(lo_l_clik_path, "/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/");
  ClikObject* lo_l_clik(0);
  std::vector<ClikPar::param_t> lo_l_nuis_enums;

  // Clik Par variable
  ClikPar* clik_par(0);

  // PLCPack variable
  PLCPack* plc_pack(0);

  // PBH variable
  std::string pbh_file_root = std::string(CLASS_PBH_FILE_DIR) + "/pbh_bspline_";

  // Use the first command line argument as a non-default
  // output root
  // That is, if an argument is even specified
  if (argc > 1) {
    settings.path = argv[1];
  }
  else {
    settings.path = "output/pc_diver-";
  }

  std::cout << "Printing results to " << settings.path << std::endl;

  // Read in parameters and set Diver settings
  clik_par = new ClikPar();

  settings.nPar = ClikPar::FREE_PARAM_AMT;
  settings.NP = 10*settings.nPar;
  settings.lowerbounds = clik_par->get_lower_bounds();
  settings.upperbounds = clik_par->get_upper_bounds();
  settings.nDerived = ClikPar::DERIVED_PARAM_AMT;

  std::cout << "Opening " << hi_l_clik_path << std::endl;

  // Create new clik object for high l likelihood and get nuisance parameters
  hi_l_clik = new ClikObject(hi_l_clik_path);

#ifndef LITE_HI_L
// TT & TTTEEE
  hi_l_nuis_enums.push_back(ClikPar::A_cib_217);
  hi_l_nuis_enums.push_back(ClikPar::cib_index);
  hi_l_nuis_enums.push_back(ClikPar::xi_sz_cib);
  hi_l_nuis_enums.push_back(ClikPar::A_sz);
  hi_l_nuis_enums.push_back(ClikPar::ps_A_100_100);
  hi_l_nuis_enums.push_back(ClikPar::ps_A_143_143);
  hi_l_nuis_enums.push_back(ClikPar::ps_A_143_217);
  hi_l_nuis_enums.push_back(ClikPar::ps_A_217_217);
  hi_l_nuis_enums.push_back(ClikPar::ksz_norm);
  hi_l_nuis_enums.push_back(ClikPar::gal545_A_100);
  hi_l_nuis_enums.push_back(ClikPar::gal545_A_143);
  hi_l_nuis_enums.push_back(ClikPar::gal545_A_143_217);
  hi_l_nuis_enums.push_back(ClikPar::gal545_A_217);
// TTTEEE
  hi_l_nuis_enums.push_back(ClikPar::galf_EE_A_100);
  hi_l_nuis_enums.push_back(ClikPar::galf_EE_A_100_143);
  hi_l_nuis_enums.push_back(ClikPar::galf_EE_A_100_217);
  hi_l_nuis_enums.push_back(ClikPar::galf_EE_A_143);
  hi_l_nuis_enums.push_back(ClikPar::galf_EE_A_143_217);
  hi_l_nuis_enums.push_back(ClikPar::galf_EE_A_217);
  hi_l_nuis_enums.push_back(ClikPar::galf_EE_index);
  hi_l_nuis_enums.push_back(ClikPar::galf_TE_A_100);
  hi_l_nuis_enums.push_back(ClikPar::galf_TE_A_100_143);
  hi_l_nuis_enums.push_back(ClikPar::galf_TE_A_100_217);
  hi_l_nuis_enums.push_back(ClikPar::galf_TE_A_143);
  hi_l_nuis_enums.push_back(ClikPar::galf_TE_A_143_217);
  hi_l_nuis_enums.push_back(ClikPar::galf_TE_A_217);
  hi_l_nuis_enums.push_back(ClikPar::galf_TE_index);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_0T_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_0T_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_0T_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_0T_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_0T_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_0T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_0T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_0T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_0T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_0T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_0T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_0T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_0T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_0T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_0T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_1T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_1T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_1T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_1T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_1T_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_1T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_1T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_1T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_1T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_1T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_2T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_2T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_2T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_2T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_2T_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_0E_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_0E_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_0E_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_0E_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_0E_0E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_0E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_0E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_0E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_0E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_0E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_0E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_0E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_0E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_0E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_0E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_1E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_1E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_1E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_1E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_1E_1E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_1E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_1E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_1E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_1E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_1E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_0_2E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_1_2E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_2_2E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_3_2E_2E);
  hi_l_nuis_enums.push_back(ClikPar::bleak_epsilon_4_2E_2E);
// TT & TTTEEE
  hi_l_nuis_enums.push_back(ClikPar::calib_100T);
  hi_l_nuis_enums.push_back(ClikPar::calib_217T);
// TTTEEE
  hi_l_nuis_enums.push_back(ClikPar::calib_100P);
  hi_l_nuis_enums.push_back(ClikPar::calib_143P);
  hi_l_nuis_enums.push_back(ClikPar::calib_217P);
  hi_l_nuis_enums.push_back(ClikPar::A_pol);
#endif
  hi_l_nuis_enums.push_back(ClikPar::A_planck);

  hi_l_clik->set_nuisance_param_enums(hi_l_nuis_enums);

  // Print out nuisance parameter names for the world to see
  std::cout << *hi_l_clik;

  std::cout << "Opening " << lo_l_clik_path << std::endl;

  // Create new clik object for low l likelihood and get nuisance parameters
  lo_l_clik = new ClikObject(lo_l_clik_path);

  lo_l_nuis_enums.push_back(ClikPar::A_planck);

  lo_l_clik->set_nuisance_param_enums(lo_l_nuis_enums);

  // Print out nuisance parameter names for the world to see
  std::cout << *lo_l_clik;


  // Package it all together
  plc_pack = new PLCPack();
  plc_pack->set_clik_params(clik_par);
  plc_pack->add_clik_object(hi_l_clik);
  plc_pack->add_clik_object(lo_l_clik);

  // Initialise CLASS before runing MultiNest
  // plc_pack->read_pbh_files(pbh_file_root);
  plc_pack->initialise_CLASS();

  context = plc_pack;

  // Calling Diver

  pc_diver(objective, pc_prior, settings, context);

  delete plc_pack;

  return 0;
  
}
