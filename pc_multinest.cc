#include "clik.h"
#include "ClassEngine.hh"
#include "multinest.h"

#include <cstdio> // for stderr

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <math.h>

static void print_usage(std::string name);
std::vector<unsigned> create_l_vec(int max_l);
void convert_Dl_to_Cl(std::vector<double> *dl_vec);
void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context);

/*
  Uses CLASS (wherever that will be with GAMBIT) to produce spectra
  which are used in the log likelihood using PLC to derive the log likelihood. Then MultiNest does the rest.

  5/10/16 - tested plc_class against Multinest for no free params.
            Shows identical agreement against plc_class for identical
            nuisance parameter input and CLASS initial conditions.
            Cannot find out how to use *context, though.
*/
int main(int argc, char *argv[]) {

  // MultiNest variables
  int IS = 0;             // do Nested Importance Sampling?
  int mmodal = 1;         // do mode separation?
  int ceff = 0;           // run in constant efficiency mode?
  int nlive = 500;        // number of live points
  double efr = 0.8;       // set the required efficiency
  double tol = 0.1;       // tol, defines the stopping criteria
  int ndims = 15;         // dimensionality (no. of free parameters)
  int nPar = 22;          // total no. of parameters including free &
                          // derived parameters
  int nClsPar = 15;       // no. of parameters to do mode separation on
  int updInt = 1000;      // after how many iterations feedback is
                          // required & the output files should be updated
                          // note: posterior files are updated & dumper
                          // routine is called after every updInt*10
                          // iterations
  double Ztol = -1E90;    // all the modes with logZ < Ztol are ignored
  int maxModes = 100;     // expected max no. of modes (used only for
                          // memory allocation)
  int pWrap[ndims];       // which parameters to have periodic boundary
                          // conditions?
  for(int i = 0; i < ndims; i++) pWrap[i] = 0;
  char root[100] = "output/pc_multinest_500_0.1_nuis-";      // root for output files
  int seed = -1;          // random no. generator seed, if < 0 then take
                          // the seed from system clock
  int fb = 1;             // need feedback on standard output?
  int resume = 1;         // resume from a previous job?
  int outfile = 1;        // write output files?
  int initMPI = 1;        // initialize MPI routines?, relevant only if
                          // compiling with MPI
                          // set it to F if you want your main program to
                          // handle MPI initialization
  double logZero = -1E90; // points with loglike < logZero will be ignored
                          // by MultiNest
  int maxiter = 0;        // max no. of iterations, a non-positive value
                          // means infinity. MultiNest will terminate if
                          // either it has done max no. of iterations or
                          // convergence criterion (defined through tol)
                          // has been satisfied
  void *context = 0;      // not required by MultiNest, any additional
                          // information user wants to pass


  // My own variables
  const int CL_AMT = 6;
  int do_lensing = 0; // input flags
  std::vector<std::string> spectra_names;

  // Variables for PLC manipulations
  char *clik_path = "/data/harryp/pc_multinest/plik_dx11dr2_HM_v18_TT.clik/";
  error *_err, **err;
  clik_object *clik_id;
  // clik_lensing_object *clik_lens_id;
  int is_lensed;
  int cl_flags_plc[CL_AMT], l_maxes[CL_AMT], max_l, prev_max_l;
  parname *param_names;
  int param_amts;

  // Initialise error (the PLC way)
  _err = initError();
  err = &_err;

  // Check if likelihood uses lensing (it should if it is physical!)
  // TODO: find if lensing is physical I guess (21/9/16)
  is_lensed = clik_try_lensing(clik_path, err);
  quitOnError(*err, __LINE__, stderr);

  if (is_lensed != do_lensing) {
    std::cerr << "[ERROR] Inconsistent lensing requirement for .clik "
              << "file: '"
              << clik_path
              << "'. Check your command line arguments or .clik file!"
              << std::endl;
    throw std::exception();
  }

  // Initialise Planck likelihood
  clik_id = clik_init(clik_path, err);
  quitOnError(*err, __LINE__, stderr);
  // clik_lens_id = clik_lensing_init(clik_path, err);
  // What do lensed spectra have as their cl_flags_plc??

  // cl_flags_plc is array of six integers (0 or 1) indicating presence of
  // spectra in following order:
  //   TT EE BB TE TB EB
  // spectra_names corresponds with cl_flags_plc
  clik_get_has_cl(clik_id, cl_flags_plc, err);
  quitOnError(*err, __LINE__, stderr);

  // l_maxes is a similar array containing maximum multipole of each
  // spectra present in the .clik file
  clik_get_lmax(clik_id, l_maxes, err);
  quitOnError(*err, __LINE__, stderr);

  // Check if the Cls are the same length
  max_l = -1;
  prev_max_l = -1;
  bool first = true; // flag for first spectra being checked
  for (int i = 0; i < CL_AMT; ++i) {
    if (cl_flags_plc[i] == 1) {
      max_l = l_maxes[i];
      
      if ( (max_l != prev_max_l) and (not first) ) {
        std::cerr << "[ERROR] The .clik file has uneven spectra lengths. "
                  << "This program cannot currently handle this!"
                  << std::endl;
        throw std::exception();
      }

      if (first) {
        first = false;
      }

      prev_max_l = max_l;
    }
  }

  // Check if CLASS can generate given spectra
  for (int i = 0; i < CL_AMT; ++i) {
    if (cl_flags_plc[i] == 1) {
      // i > 3 correspond to spectra which cannot be generated by CLASS
      if (i > 3) {
        std::cout << "CLASS cannot create this spectra!" << std::endl;
        std::cerr << "[ERROR] Invalid CLASS spectra, "
                  << spectra_names.at(i)
                  << "."
                  << std::endl;
        throw std::exception();
      }

      std::cout << std::endl;
    }
  }

  // Print out nuisance parameter names for the world to see
  param_amts = clik_get_extra_parameter_names(clik_id, &param_names, err);

  std::cout << "This .clik file requires "
            << param_amts
            << " nuisance parameters to be marginalised over:\n";
  for (int i = 0; i < param_amts; ++i) {
    std::cout << '\t' << param_names[i] << std::endl;
  }
  std::cout << std::endl;


  // context = clik_path;
  context = clik_id;
  
  // Calling MultiNest

  nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar,
              maxModes, updInt, Ztol, root, seed, pWrap, fb,
              resume, outfile, initMPI, logZero, maxiter, LogLike, dumper,
              context);

  clik_cleanup(&clik_id);

  return 0;
}


/*** Secondary Functions ***/

// This method is not actually used, haha!
static void print_usage(std::string name) {
  std::cerr << "Usage: " << name << " [-h -v -l] PLIK_FILE\n"
            << "Options:\n"
            << "\t-h,--help\t\tShow this help message\n"
            << "\t-v,--verbose\t\tVerbose output\n"
            << "\t-l,--do-lensing\t\tDo a lensing likelihood computation"
            // << "\t-c,--class-arg CLASS_INI\t\tGive .ini file CLASS_INI for CLASS to initialise with"
            << std::endl;
}

// Create vector of multipoles to pass to class_engine.getCls()
std::vector<unsigned> create_l_vec(int max_l) {
  std::vector<unsigned> l_vec;
  for (int l = 2; l <= max_l; ++l) {
    l_vec.push_back(l);
  }
  return l_vec;
}

// Converts Dl=l(l+1)Cl/2pi to Cl again
// Not really needed as of now (23/9/16)
void convert_Dl_to_Cl(std::vector<double> *dl_vec) {
  int i;
  double l;

  for (i = 0; i < dl_vec->size(); ++i) {
    l = double(i) + 2.0;
    // dl_vec->at(i) = dl_vec->at(i)*l*(l+1)/(2.0*M_PI); // one for 
    // dl_vec->at(i) = dl_vec->at(i)*2.0*M_PI/(l*(l+1)); // the other
  }
}

/***********************/
/*      MultiNest      */
/***********************/

// Input arguments
// ndim             = dimensionality (total number of free parameters) of the problem
// npars            = total number of free plus derived parameters
// context            void pointer, any additional information
//
// Input/Output arguments
// Cube[npars]      = on entry has the ndim parameters as [0,1]
//                    on exit, the physical parameters plus copy any
//                    derived parameters you want to store with the free
//                    parameters
//   
// Output arguments
// lnew             = loglikelihood
void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
  // ndim     = # of CLASS parameters + nuisance parameters
  //          = 6 + 16 = 22
  // npars    = total dimension of Cube
  //          = 22 + # free parameters
  // Cube     = can use nuisance parameter part to construct cl_and_pars
  // context  = PLC object. has to be! In the future for now.
  
  enum param_t
  {
    // Free parameters (PLC)
    _A_cib_217,
    _xi_sz_cib,
    _A_sz,
    _ps_A_100_100,
    _ps_A_143_143,
    _ps_A_143_217,
    _ps_A_217_217,
    _ksz_norm,
    _gal545_A_100,
    _gal545_A_143,
    _gal545_A_143_217,
    _gal545_A_217,
    _calib_100T,
    _calib_217T,
    _A_planck,
    // Free parameters (CLASS)
    _omega_b,
    _omega_cdm,
    _hundredxtheta_s,
    _tau_reio,
    _n_s,
    _ln10_10_A_s,
    // Derived parameters
    _cib_index // cib_index = -1.3 is constant
  };


  // My own variables
  const int CL_AMT = 6;
  double par_min[npars], par_max[npars];
  // par_min = new double[ndim];
  // par_max = new double[ndim];

  // Variables for PLC manipulations
  error *_err, **err;
  clik_object *clik_id;
  double *cl_and_pars;
  int cap_size;
  int cl_flags_plc[CL_AMT], l_maxes[CL_AMT], max_l;
  parname *param_names;
  int param_amts;

  // CLASS variables
  ClassParams class_params;
  ClassEngine *class_engine(0);
  std::vector<unsigned> l_vec;
  std::vector<double> cl_tt, cl_ee, cl_bb, cl_te, nuisance_pars;
  std::vector<std::vector<double> > class_cls;

  // Initialise the nuisance parameter priors using results outlined in
  //  Planck 2015 results. XI. CMB power spectra, ... . p 21,41

  // Input parameters for CLASS
  par_min[_omega_b] = 0.005;        par_max[_omega_b] = 0.1;
  par_min[_omega_cdm] = 0.001;      par_max[_omega_cdm] = 0.99;
  par_min[_hundredxtheta_s] = 0.5;  par_max[_hundredxtheta_s] = 10.0;
  par_min[_tau_reio] = 0.01;        par_max[_tau_reio] = 0.8;
  par_min[_n_s] = 0.8;              par_max[_n_s] = 1.2;
  par_min[_ln10_10_A_s] = 2.0;      par_max[_ln10_10_A_s] = 4.0;

  // Nuisance parameters for PLC
  par_min[_A_cib_217] = 0.0;        par_max[_A_cib_217] = 100.0;
  par_min[_cib_index] = -1.3;       par_max[_cib_index] = -1.3;
  par_min[_xi_sz_cib] = 0.0;        par_max[_xi_sz_cib] = 0.5;
  par_min[_A_sz] = 5.0;             par_max[_A_sz] = 10.0;
  par_min[_ps_A_100_100] = 200.0;   par_max[_ps_A_100_100] = 300.0;
  par_min[_ps_A_143_143] = 0.0;     par_max[_ps_A_143_143] = 100.0;
  par_min[_ps_A_143_217] = 0.0;     par_max[_ps_A_143_217] = 100.0;
  par_min[_ps_A_217_217] = 50.0;    par_max[_ps_A_217_217] = 150.0;
  par_min[_ksz_norm] = 0.0;         par_max[_ksz_norm] = 7.0;
  par_min[_gal545_A_100] = 0.0;     par_max[_gal545_A_100] = 20.0;
  par_min[_gal545_A_143] = 0.0;     par_max[_gal545_A_143] = 20.0;
  par_min[_gal545_A_143_217] = 0.0; par_max[_gal545_A_143_217] = 30.0;
  par_min[_gal545_A_217] = 60.0;     par_max[_gal545_A_217] = 100.0;
  par_min[_calib_100T] = 0.0;       par_max[_calib_100T] = 2.0;
  par_min[_calib_217T] = 0.0;       par_max[_calib_217T] = 2.0;
  par_min[_A_planck] = 0.9;         par_max[_A_planck] = 1.1;

  // Scale parameters
  for (int i = 0; i < npars; i++) {
    if (i < ndim) { // free parameter
      Cube[i] = par_min[i] + (par_max[i] - par_min[i]) * Cube[i];
    }
    else { // derived parameter
      Cube[i] = (par_min[i] + par_max[i]) / 2.0;
    }
  }

  // This is a set parameter
  Cube[_cib_index] = -1.3;

  // Setting to Planck best fit for base_plikHM_TT_lowTEB
  // Just as a test for Gaussianity of A_cib_217
  // Cube[_A_cib_217] = 66.6;
  // Cube[_xi_sz_cib] = 0.05;
  // Cube[_A_sz] = 7.14;
  // Cube[_ps_A_100_100] = 251.8;
  // Cube[_ps_A_143_143] = 39.2;
  // Cube[_ps_A_143_217] = 33.6;
  // Cube[_ps_A_217_217] = 97.8;
  // Cube[_ksz_norm] = 0.0;
  // Cube[_gal545_A_100] = 7.41;
  // Cube[_gal545_A_143] = 8.98;
  // Cube[_gal545_A_143_217] = 17.53;
  // Cube[_gal545_A_217] = 82.0;
  // Cube[_calib_100T] = 0.99789;
  // Cube[_calib_217T] = 0.99593;
  // Cube[_A_planck] = 1.00030;

  
  /* Old code begins now */


  /******************************/
  /*   Planck Likelihood Code   */
  /******************************/

  // Initialise error (the PLC way)
  _err = initError();
  err = &_err;

  clik_id = static_cast<clik_object*>(context);

  // cl_flags_plc is array of six integers (0 or 1) indicating presence of
  // spectra in following order:
  //   TT EE BB TE TB EB
  // spectra_names corresponds with cl_flags_plc
  // std::cout << "Help me again\n";
  clik_get_has_cl(clik_id, cl_flags_plc, err);
  quitOnError(*err, __LINE__, stderr);
  // std::cout << "Please\n";

  // l_maxes is a similar array containing maximum multipole of each
  // spectra present in the .clik file
  clik_get_lmax(clik_id, l_maxes, err);
  quitOnError(*err, __LINE__, stderr);

  // Get maximum l from .clik file
  // Already performed checks, already know each spectra has same length
  max_l = l_maxes[0];

  // Get number of nuisance parameters
  param_amts = clik_get_extra_parameter_names(clik_id, &param_names, err);

  // Determine size of cls_and_pars array
  cap_size = 0;
  for (int i = 0; i < CL_AMT; ++i) {
    if (cl_flags_plc[i] == 1) {
      cap_size += l_maxes[i] + 1; // +1 for l=0 case
    }
  }
  cap_size += param_amts;


  /******************************/
  /*           CLASS            */
  /******************************/

  // MultiNest parameters to sweep over (be careful here!)
  // class_params.add("omega_b", Cube[_omega_b]);
  // class_params.add("omega_cdm", Cube[_omega_cdm]);
  // class_params.add("100*theta_s", Cube[_hundredxtheta_s]);
  // class_params.add("tau_reio", Cube[_tau_reio]);
  // class_params.add("n_s", Cube[_n_s]); // k_0 = 0.05 Mpc^-1 by default
  // class_params.add("ln10^{10}A_s", Cube[_ln10_10_A_s]);

  // Set variables to Planck best fit
  class_params.add("omega_b", 0.022242);
  class_params.add("omega_cdm", 0.11977);
  class_params.add("100*theta_s", 1.040862);
  class_params.add("tau_reio", 0.0781);
  class_params.add("n_s", 0.9658); // k_0 = 0.05 Mpc^-1 by default
  class_params.add("ln10^{10}A_s", 3.0904);

  // Options to set for spectra output
  class_params.add("output", "tCl,pCl"); // pCl, lCl for lensed spectra
  class_params.add("l_max_scalars", max_l);
  class_params.add("format", "camb");

  try {
    class_engine = new ClassEngine(class_params);
  }
  catch (std::exception &e) {
    std::cerr << "[ERROR] CLASS failed, throwing exception " << e.what() << std::endl;
    throw e;
  }

  // Create vector of multipoles to calculate Cls for
  l_vec = create_l_vec(max_l);

  // Calculate Cls from CLASS
  try {
    class_engine->getCls(l_vec, cl_tt, cl_te, cl_ee, cl_bb);
  }
  catch (std::exception &e) {
    std::cerr << "[ERROR] Spectra extraction unsuccessful, threw " << e.what() << std::endl;
    throw e;
  }


  /******************************/
  /*   Likelihood Calculation   */
  /******************************/

  // Push all spectra into a single matrix
  class_cls.push_back(cl_tt);
  class_cls.push_back(cl_ee);
  class_cls.push_back(cl_bb);
  class_cls.push_back(cl_te);

  // Create parameter array
  nuisance_pars.push_back(Cube[_A_cib_217]);
  nuisance_pars.push_back(Cube[_cib_index]);
  nuisance_pars.push_back(Cube[_xi_sz_cib]);
  nuisance_pars.push_back(Cube[_A_sz]);
  nuisance_pars.push_back(Cube[_ps_A_100_100]);
  nuisance_pars.push_back(Cube[_ps_A_143_143]);
  nuisance_pars.push_back(Cube[_ps_A_143_217]);
  nuisance_pars.push_back(Cube[_ps_A_217_217]);
  nuisance_pars.push_back(Cube[_ksz_norm]);
  nuisance_pars.push_back(Cube[_gal545_A_100]);
  nuisance_pars.push_back(Cube[_gal545_A_143]);
  nuisance_pars.push_back(Cube[_gal545_A_143_217]);
  nuisance_pars.push_back(Cube[_gal545_A_217]);
  nuisance_pars.push_back(Cube[_calib_100T]);
  nuisance_pars.push_back(Cube[_calib_217T]);
  nuisance_pars.push_back(Cube[_A_planck]);

  // Check if we have done enough initialising
  if (nuisance_pars.size() > param_amts) {
    std::cerr << "[ERROR] More nuisance paramters initialised than needed!"
              << std::endl;
    throw std::exception();
  }
  else if (nuisance_pars.size() < param_amts) {
    std::cerr << "[ERROR] Not enough nuisance parameters have been initialised!"
              << std::endl;
    throw std::exception();
  }

  // Construct the cl_and_pars array that PLC wants
  // First transpose spectra
  cl_and_pars = new double[cap_size];
  int cap_ind = 0;
  for (int cl_ind = 0; cl_ind < CL_AMT; ++cl_ind) {
    if (cl_flags_plc[cl_ind] == 1) { // if spectra exists
      for (int l = 0; l <= l_maxes[cl_ind]; ++l) {
        // Correct for PLC wanting l=0,1 multipoles
        // These cannot be computed by CLASS and so are set to zero
        if (l < 2) {
          cl_and_pars[cap_ind] = 0.0;
        }
        else {
          // CLASS Cls start from l=2, hence the l-2
          cl_and_pars[cap_ind] = class_cls[cl_ind][l-2];
        }

        cap_ind++;
      }
    }
  }

  // Then transpose nuisance parameters
  for (int i = 0; i < param_amts; ++i) {
    cl_and_pars[cap_ind] = nuisance_pars[i];
    cap_ind++;
  }

  // Compute the log likelihood using PLC
  lnew = clik_compute(clik_id, cl_and_pars, err);
  quitOnError(*err, __LINE__, stderr);

  // This will get old pretty fast
  /*
  std::cout << "[plc_class] Calculated log likelihood of "
            << lnew
            << std::endl;
  */

  // Clean up after running
  delete class_engine;
  delete[] cl_and_pars;
}

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants

// Arguments:
//
// nSamples         = total number of samples in posterior distribution
// nlive            = total number of live points
// nPar             = total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)]      = 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)]      = posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1]  = mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1]   = standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike       = maximum loglikelihood value
// logZ             = log evidence value from the default (non-INS) mode
// INSlogZ          = log evidence value from the INS mode
// logZerr          = error on log evidence value
// context          = void pointer, any additional information

// Default example from eggbox.cc
void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
  // convert the 2D Fortran arrays to C++ arrays
  
  // the posterior distribution
  // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
  
  int i, j;
  
  double postdist[nSamples][nPar + 2];
  for( i = 0; i < nPar + 2; i++ )
    for( j = 0; j < nSamples; j++ )
      postdist[j][i] = posterior[0][i * nSamples + j];
  
  // last set of live points
  // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
  
  double pLivePts[nlive][nPar + 1];
  for( i = 0; i < nPar + 1; i++ )
    for( j = 0; j < nlive; j++ )
      pLivePts[j][i] = physLive[0][i * nlive + j];
}

/*
nested::run(int IS, int mmodal, int ceff, int nlive, double tol,
            double efr, int ndims, int nPar, int nClsPar,
            int maxModes, int updInt, double Ztol, char root[100],
            int seed, int pWrap[ndims], int fb, int resume,
            int outfile, int initMPI, double logZero, int maxiter,
            LogLike, dumper, void *context);

int IS = 0;             // do Nested Importance Sampling?

int mmodal = 1;         // do mode separation?

int ceff = 0;           // run in constant efficiency mode?

int nlive = 1000;       // number of live points

double efr = 0.8;       // set the required efficiency

double tol = 0.5;       // tol, defines the stopping criteria

int ndims = 2;          // dimensionality (no. of free parameters)

int nPar = 2;           // total no. of parameters including free & derived parameters

int nClsPar = 2;        // no. of parameters to do mode separation on

int updInt = 1000;      // after how many iterations feedback is required
                           & the output files should be updated
                        // note: posterior files are updated & dumper
                           routine is called after every updInt*10 iterations

double Ztol = -1E90;    // all the modes with logZ < Ztol are ignored

int maxModes = 100;     // expected max no. of modes (used only for
                           memory allocation)

int pWrap[ndims];       // which parameters to have periodic boundary
                           conditions?
for(int i = 0; i < ndims; i++) pWrap[i] = 0;

char root[100] = "chains/eggboxCC-";      // root for output files

int seed = -1;          // random no. generator seed, if < 0 then take
                           the seed from system clock

int fb = 1;             // need feedback on standard output?

int resume = 1;         // resume from a previous job?

int outfile = 1;        // write output files?

int initMPI = 1;        // initialize MPI routines?, relevant only if
                           compiling with MPI
                        // set it to F if you want your main program to
                           handle MPI initialization

double logZero = -1E90; // points with loglike < logZero will be ignored
                           by MultiNest

int maxiter = 0;        // max no. of iterations, a non-positive value
                           means infinity. MultiNest will terminate if
                           either it 
                        // has done max no. of iterations or convergence
                           criterion (defined through tol) has been
                           satisfied

void *context = 0;      // not required by MultiNest, any additional
                           information user wants to pass

*/