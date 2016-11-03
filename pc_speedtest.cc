#include "ClassEngine.hh"
#include "PLCPack.h"

#include <cstdio> // for stderr
#include <ctime> // for clock_t, clock()
#include <sys/time.h> // for gettimeofday()

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept> // gives various standard exceptions
#include <math.h>

#include <fstream> // for writing loop times out to file
#include <iomanip> // setw, setprecision

std::vector<unsigned> create_l_vec(int max_l);
void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);

double get_CPU_time();
double get_wall_time();
double get_accurate_wall_time();

/*
  Uses CLASS (wherever that will be with GAMBIT) to produce spectra
  which are used in the log likelihood using PLC to derive the log likelihood. Then MultiNest does the rest.

  5/10/16 - tested plc_class against Multinest for no free params.
            Shows identical agreement against plc_class for identical
            nuisance parameter input and CLASS initial conditions.
            Cannot find out how to use *context, though.
*/
int main(int argc, char *argv[]) {
  // MultiNest-like variables to replicate its behaviour
  int ndims = 1;
  int nPar = 22;
  double Cube[nPar];
  double lnew;
  void *context = 0;

  // Looping variables
  int LOOP_AMT = 1000;
  double step_size = 1.0/double(LOOP_AMT);
  double CPU_times[LOOP_AMT], wall_times[LOOP_AMT],
         acc_wall_times[LOOP_AMT], A_cib_217_values[LOOP_AMT],
         likelihood_values[LOOP_AMT];
  double CPU_time_before, wall_time_before, acc_wall_time_before;

  // Timing results printing variables
  std::ofstream time_file;
  std::string time_file_dir;
  const int DOUBLE_PRECIS = 6;

  // My own variables
  const int CL_AMT = 6;
  int do_lensing = 0; // input flags
  PLCPack *plc_pack = 0;

  // Variables for PLC manipulations
  char *clik_path = "/home/harry/plc_2.0/hi_l/plik/plik_dx11dr2_HM_v18_TT.clik/";
  error *_err, **err;
  clik_object *clik_id;
  int is_lensed;
  parname *param_names;
  int param_amts;

  // Extract out file location
  if (argc == 1) { // default directory
    time_file_dir = "output/A_cib_217_times.txt";
  }
  else if (argc == 2) { // custom directory
    time_file_dir = argv[1];
  }
  else {
    std::cerr << "Too many args given!" << std::endl;
    throw std::exception();
  }

  std::cout << "Printing results to " << time_file_dir << std::endl;

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

  // Create new PLCPack object for the clik_object to live in
  plc_pack = new PLCPack(clik_id);

  // Print out nuisance parameter names for the world to see
  param_amts = clik_get_extra_parameter_names(clik_id, &param_names, err);

  std::cout << "This .clik file requires "
            << param_amts
            << " nuisance parameters to be marginalised over:\n";
  for (int i = 0; i < param_amts; ++i) {
    std::cout << '\t' << param_names[i] << std::endl;
  }
  std::cout << std::endl;

  context = plc_pack;

  // Manually looping over LogLike method
  for (int i = 0; i < LOOP_AMT; ++i) {
    // Initialise The Cube
    for (int ipar = 0; ipar < nPar; ++ipar) {
      Cube[ipar] = double(i)*step_size;
    }

    std::cout << "------------------------"<< std::endl;
    std::cout << "Loop number " << i << std::endl;

    // Time the LogLike method
    // CPU_time_before = get_CPU_time();
    // wall_time_before = get_wall_time();
    // acc_wall_time_before = get_accurate_wall_time();

    LogLike(Cube, ndims, nPar, lnew, context);

    // LogLike times in seconds
    // CPU_times[i] = get_CPU_time() - CPU_time_before;
    // wall_times[i] = get_wall_time() - wall_time_before;
    // acc_wall_times[i] = get_accurate_wall_time() - acc_wall_time_before;

    // Also extract value of A_cib_217 and likelihood
    // A_cib_217_values[i] = Cube[0];
    // likelihood_values[i] = lnew;

    /*
    std::cout << "CPU time: "
              << std::setprecision(DOUBLE_PRECIS) << CPU_times[i]
              << " s." << std::endl;
    std::cout << "Wall time: " 
              << std::setprecision(DOUBLE_PRECIS) << wall_times[i]
              << " s." << std::endl;
    std::cout << "Accurate wall time: "
              << std::setprecision(DOUBLE_PRECIS) << acc_wall_times[i] 
              << " s." << std::endl;
    std::cout << "A_cib_217: " << A_cib_217_values[i] << std::endl;
    std::cout << "Likelihood: " << likelihood_values[i] << std::endl;
    //*/
  }

  // Print contents of CPU_times to file
  // Changed to tab delimiter 23/10/16
  /*
  time_file.open(time_file_dir.c_str());
  time_file << "A_cib_217 value" << '\t'
            << "Log likelihood" << '\t'
            << "CPU time (s)" << '\t'
            << "Wall time (s)" << '\t'
            << "Accurate wall time (s)"
            << std::endl;
  for (int i = 0; i < LOOP_AMT; ++i) {
    time_file << A_cib_217_values[i] << '\t'
              << likelihood_values[i] << '\t'
              << std::setprecision(DOUBLE_PRECIS) << CPU_times[i]
              << '\t'
              << std::setprecision(DOUBLE_PRECIS) << wall_times[i]
              << '\t'
              << std::setprecision(DOUBLE_PRECIS) << acc_wall_times[i]
              << std::endl;
  }
  time_file.close();
  //*/

  delete plc_pack;

  return 0;
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
  PLCPack *plc_pack;
  double *cl_and_pars;
  int cap_size;
  int *cl_flags_plc, max_l;
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
  par_min[_A_cib_217] = 0.0;        par_max[_A_cib_217] = 200.0;
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
  Cube[_xi_sz_cib] = 0.05;
  Cube[_A_sz] = 7.14;
  Cube[_ps_A_100_100] = 251.8;
  Cube[_ps_A_143_143] = 39.2;
  Cube[_ps_A_143_217] = 33.6;
  Cube[_ps_A_217_217] = 97.8;
  Cube[_ksz_norm] = 0.0;
  Cube[_gal545_A_100] = 7.41;
  Cube[_gal545_A_143] = 8.98;
  Cube[_gal545_A_143_217] = 17.53;
  Cube[_gal545_A_217] = 82.0;
  Cube[_calib_100T] = 0.99789;
  Cube[_calib_217T] = 0.99593;
  Cube[_A_planck] = 1.00030;

  
  /* plc_class code begins */


  /******************************/
  /*   Planck Likelihood Code   */
  /******************************/

  plc_pack = static_cast<PLCPack*>(context);

  // cl_flags_plc is array of six integers (0 or 1) indicating presence of
  // spectra in following order:
  //   TT EE BB TE TB EB
  cl_flags_plc = plc_pack->get_cl_flags();

  // Get maximum l from .clik file
  max_l = plc_pack->get_max_l();

  // Get number of nuisance parameters
  param_amts = plc_pack->get_param_amt();

  // Get size of cls_and_pars array
  cap_size = plc_pack->get_cap_size();


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
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS failed, throwing exception " << e.what() << std::endl;
    throw e;
  }

  // Create vector of multipoles to calculate Cls for
  l_vec = plc_pack->get_class_l_vec();

  // Calculate Cls from CLASS
  try {
    class_engine->getCls(l_vec, cl_tt, cl_te, cl_ee, cl_bb);
  }
  catch (std::exception const &e) {
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
    std::cerr << "[ERROR] More nuisance parameters initialised than needed!"
              << std::endl;
    throw std::exception();
  }
  else if (nuisance_pars.size() < param_amts) {
    std::cerr << "[ERROR] Not enough nuisance parameters have been initialised!"
              << std::endl;
    throw std::exception();
  }

  // Construct the cl_and_pars array that PLC wants
  // First add spectra
  cl_and_pars = new double[cap_size];
  int cap_ind = 0;
  for (int cl_ind = 0; cl_ind < CL_AMT; ++cl_ind) {
    if (cl_flags_plc[cl_ind] == 1) { // if spectra exists
      for (int l = 0; l <= max_l; ++l) {
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

  // Then add nuisance parameters
  for (int i = 0; i < param_amts; ++i) {
    cl_and_pars[cap_ind] = nuisance_pars[i];
    cap_ind++;
  }
  // cl_and_pars constructed
  // cap_ind should equal cap_size

  // Compute the log likelihood using PLC
  lnew = plc_pack->get_clik_likelihood(cl_and_pars);

  /*
  std::cout << "[plc_class] Calculated log likelihood of "
            << lnew
            << std::endl;
  */

  // Clean up after running
  delete class_engine;
  delete[] cl_and_pars;
}


/***********************/
/*      SpeedTest      */
/***********************/

// Returns the sum of user + system time
double get_CPU_time() {
  return clock() / (double) CLOCKS_PER_SEC;
}

// Might be prone to errors. POSIX.1-2008 recommends against using this
// for accurate timing due to associated hardware dependencies.
double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    std::cerr << "ERROR: gettimeofday() threw an error!"
              << std::endl;
    return 0;
  }
  return (double) time.tv_sec + (double) time.tv_usec * 0.000001;
}

// Uses a different method of getting walltime that shouldn't depend
// on how many processes are running nor over how many cores.
// Has nanosecond precision.
double get_accurate_wall_time() {
  struct timespec time;
  if (clock_gettime(CLOCK_MONOTONIC, &time)) {
    std::cerr << "ERROR: clock_gettime() threw an error!"
              << std::endl;
    return 0;
  }
  return (double) time.tv_sec + (double) time.tv_nsec * 0.000000001;  
}