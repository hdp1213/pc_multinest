#include "ClassEngine.hh"
#include "PLCPack.h"
#include "ClikObject.h"

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
  int LOOP_AMT = 20;
  double step_size = 1.0/double(LOOP_AMT);
  double CPU_times[LOOP_AMT], wall_times[LOOP_AMT],
         acc_wall_times[LOOP_AMT], A_cib_217_values[LOOP_AMT],
         likelihood_values[LOOP_AMT];
  double CPU_time_before, wall_time_before, acc_wall_time_before;

  // Timing results printing variables
  std::ofstream time_file;
  std::string time_file_dir;
  const int DOUBLE_PRECIS = 6;

  // High l likelihood variables
  char *hi_l_clik_path = "/data/harryp/pc_multinest/plik_dx11dr2_HM_v18_TT.clik/";
  ClikObject *hi_l_clik = 0;
  std::vector<ClikPar::param_t> hi_l_nuis_enums;

  // Low l likelihood variables  
  ClikObject *lo_l_clik = 0;
  char *lo_l_clik_path = "/data/harryp/pc_multinest/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/";
  std::vector<ClikPar::param_t> lo_l_nuis_enums;

  int param_amts;
  parname *param_names;
  PLCPack *plc_pack = 0;

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

  // Create new clik object for high l likelihood
  hi_l_clik = new ClikObject(hi_l_clik_path);

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
  hi_l_nuis_enums.push_back(ClikPar::calib_100T);
  hi_l_nuis_enums.push_back(ClikPar::calib_217T);
  hi_l_nuis_enums.push_back(ClikPar::A_planck);

  hi_l_clik->set_nuisance_param_enums(hi_l_nuis_enums);

  // Print out nuisance parameter names for the world to see
  param_amts = hi_l_clik->get_param_amt();
  param_names = hi_l_clik->get_param_names();

  std::cout << "The high-l .clik file requires "
            << param_amts
            << " nuisance parameters to be marginalised over:\n";
  for (int i = 0; i < param_amts; ++i) {
    std::cout << '\t' << param_names[i] << std::endl;
  }
  std::cout << std::endl;


  // Create new clik object for low l likelihood
  lo_l_clik = new ClikObject(lo_l_clik_path);

  lo_l_nuis_enums.push_back(ClikPar::A_planck);

  lo_l_clik->set_nuisance_param_enums(lo_l_nuis_enums);

  // Print out nuisance parameter names for the world to see
  param_amts = lo_l_clik->get_param_amt();
  param_names = lo_l_clik->get_param_names();

  std::cout << "The low-l .clik file requires "
            << param_amts
            << " nuisance parameters to be marginalised over:\n";
  for (int i = 0; i < param_amts; ++i) {
    std::cout << '\t' << param_names[i] << std::endl;
  }
  std::cout << std::endl;


  // Package it all together
  plc_pack = new PLCPack();
  plc_pack->add_clik_object(hi_l_clik);
  plc_pack->add_clik_object(lo_l_clik);

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
    CPU_time_before = get_CPU_time();
    wall_time_before = get_wall_time();
    acc_wall_time_before = get_accurate_wall_time();

    LogLike(Cube, ndims, nPar, lnew, context);

    // LogLike times in seconds
    CPU_times[i] = get_CPU_time() - CPU_time_before;
    wall_times[i] = get_wall_time() - wall_time_before;
    acc_wall_times[i] = get_accurate_wall_time() - acc_wall_time_before;

    // Also extract value of A_cib_217 and likelihood
    A_cib_217_values[i] = Cube[0];
    likelihood_values[i] = lnew;

    //*
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
  //*
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