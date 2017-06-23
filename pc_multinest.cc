#include "PLCPack.h"
#include "ClikObject.h"
#include "ClikPar.h"
#include "multinest.h"

#include <cstdio> // for stderr

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <math.h>
#include <iomanip> // for std::setprecision()

#include "loglike.h"

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context);

/*
  Uses CLASS (wherever that will be with GAMBIT) to produce
  spectra which are used in the log likelihood using PLC to
  derive the log likelihood. Then MultiNest does the rest.

  5/10/16 - tested plc_class against Multinest for no free
            params.
            Shows identical agreement against plc_class for
            identical nuisance parameter input and CLASS
            initial conditions.
            Cannot find out how to use *context, though.
*/
int main(int argc, char* argv[]) {

  // MultiNest variables
  int IS = 0;             // do Nested Importance Sampling?
  int mmodal = 0;         // do mode separation?
  int ceff = 1;           // run in constant efficiency mode?
                          // Bad for evidence calculations
  int nlive = 1000;       // number of live points
  double efr = 0.8;       // set the required efficiency. 0.8
                          // for parameter estimation, 0.3 for
                          // evidence
  double tol = 1E-1;      // tol, defines the stopping criteria
                          // 0.5 should give enough accuracy
  int ndims = ClikPar::FREE_PARAMS;
                          // dimensionality (no. of free
                          // parameters)
  int nPar = ClikPar::TOTAL_PARAMS;
                          // total no. of parameters including
                          // free & derived parameters
  int nClsPar = 0;        // no. of parameters to do mode
                          // separation on
  int updInt = 1000;      // after how many iterations feedback
                          // is required & the output files
                          // should be updated
                          // note: posterior files are updated
                          // & dumper routine is called after
                          // every updInt*10 iterations
  double Ztol = -1E90;    // all the modes with logZ < Ztol are
                          // ignored
  int maxModes = 100;     // expected max no. of modes (used
                          // only for memory allocation)
  int pWrap[ndims];       // which parameters to have periodic
                          // boundary conditions?
  for(int i = 0; i < ndims; i++) pWrap[i] = 0;
  char* root;             // root for output files
  int seed = -1;          // random no. generator seed, if < 0
                          // then take the seed from system
                          // clock
  int fb = 1;             // need feedback on standard output?
  int resume = 1;         // resume from a previous job?
  int outfile = 1;        // write output files?
  int initMPI = 1;        // initialize MPI routines?, relevant
                          // only if compiling with MPI
                          // set it to F if you want your main
                          // program to handle MPI
                          // initialization
  double logZero = -1E90; // points with loglike < logZero will
                          // be ignored by MultiNest
  int maxiter = 0;        // max no. of iterations, a
                          // non-positive value means infinity.
                          // MultiNest will terminate if
                          // either it has done max no. of
                          // iterations or convergence
                          // criterion (defined through tol)
                          // has been satisfied

  void* context = 0;      // not required by MultiNest, any
                          // additional information user wants
                          // to pass

  // High l full likelihood variables
#ifdef LITE_HI_L
  char* hi_l_clik_path = "/home/harry/plc_2.0/hi_l/plik_lite/plik_lite_v18_TTTEEE.clik/";
#else
  char* hi_l_clik_path = "/home/harry/plc_2.0/hi_l/plik/plik_dx11dr2_HM_v18_TTTEEE.clik/";
#endif
  ClikObject* hi_l_clik(0);
  std::vector<ClikPar::param_t> hi_l_nuis_enums;

  // Low l likelihood variables  
  char* lo_l_clik_path = "/home/harry/plc_2.0/low_l/bflike/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/";
  ClikObject* lo_l_clik(0);
  std::vector<ClikPar::param_t> lo_l_nuis_enums;

  // Clik Pars variable
  ClikPar* clik_par = new ClikPar();

  PLCPack* plc_pack(0);

  // PBH variable
  std::string pbh_file_root = "/home/harry/class/pbh/pbh_bspline_";

  // Use the first command line argument as a non-default
  // output root
  // That is, if an argument is even specified
  if (argc == 2) {
    root = argv[1];
  }
  else {
    root = "output/pc_multinest_-";
  }

  std::cout << "Printing results to " << root << std::endl;

  // Create new clik object for high l likelihood
  hi_l_clik = new ClikObject(hi_l_clik_path);

  // Push nuisance parameters in the order they appear in cl_and_pars
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

  // Create new clik object for low l likelihood
  lo_l_clik = new ClikObject(lo_l_clik_path);

  lo_l_nuis_enums.push_back(ClikPar::A_planck);

  lo_l_clik->set_nuisance_param_enums(lo_l_nuis_enums);

  // Print out nuisance parameter names for the world to see
  std::cout << *lo_l_clik;


  // Package it all together
  plc_pack = new PLCPack();
  plc_pack->add_clik_object(hi_l_clik);
  plc_pack->add_clik_object(lo_l_clik);
  plc_pack->set_clik_params(clik_par);
  plc_pack->read_pbh_files(pbh_file_root);

  // Initialise CLASS before runing MultiNest
  plc_pack->initialise_CLASS();

  context = plc_pack;


  // Calling MultiNest

  nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar,
              maxModes, updInt, Ztol, root, seed, pWrap, fb,
              resume, outfile, initMPI, logZero, maxiter, LogLike,
              dumper, context);

  delete plc_pack;

  return 0;
}


/*** Secondary Functions ***/

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
