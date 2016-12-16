#include "ClassEngine.hh"
#include "PLCPack.h"
#include "ClikObject.h"
#include "multinest.h"

#include <cstdio> // for stderr

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <math.h>

#include "loglike.h"
#include "ClikPar.h"

void convert_Dl_to_Cl(std::vector<double> *dl_vec);
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
  int ndims = 21;         // dimensionality (no. of free parameters)
  int nPar = 22;          // total no. of parameters including free &
                          // derived parameters
  int nClsPar = 21;       // no. of parameters to do mode separation on
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
  char *root;             // root for output files
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


  // High l likelihood variables
  char *hi_l_clik_path = "/home/a1648400/plc_2.0/hi_l/plik/plik_dx11dr2_HM_v18_TT.clik/";
  ClikObject *hi_l_clik = 0;
  std::vector<ClikPar::param_t> hi_l_nuis_enums;

  // Low l likelihood variables  
  ClikObject *lo_l_clik = 0;
  char *lo_l_clik_path = "/home/a1648400/plc_2.0/low_l/bflike/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/";
  std::vector<ClikPar::param_t> lo_l_nuis_enums;

  // Clik Pars variable
  ClikPar *clik_par = new ClikPar();

  int param_amts;
  parname *param_names;
  PLCPack *plc_pack = 0;

  // Use the first command line argument as a non-default output root
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
  plc_pack->set_pars(clik_par);

  context = plc_pack;
  
  // Calling MultiNest

  nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar,
              maxModes, updInt, Ztol, root, seed, pWrap, fb,
              resume, outfile, initMPI, logZero, maxiter, LogLike, dumper,
              context);

  delete plc_pack;

  return 0;
}


/*** Secondary Functions ***/

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