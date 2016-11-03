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

#include "loglike.h"

std::vector<unsigned> create_l_vec(int max_l);

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