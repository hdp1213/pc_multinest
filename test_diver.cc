#include "PLCPack.h"
#include "ClikObject.h"
#include "ClikPar.h"

#include "loglike.h"

#include <limits>

// Function to be minimized.  Corresponds to -ln(Likelihood).  Redirects to the target of context pointer.
double objective(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector, void*& context)
{
  double loglike = std::numeric_limits<double>::max();

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

int main(int argc, char** argv)
{
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
  ClikPar* clik_par = new ClikPar();

  // PLCPack variable
  PLCPack* plc_pack(0);

  // PBH variable
  std::string pbh_file_root = std::string(CLASS_PBH_FILE_DIR) + "/pbh_bspline_";

  //*
  std::cout << "Opening " << hi_l_clik_path << std::endl;

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
  //*/

  //*
  std::cout << "Opening " << lo_l_clik_path << std::endl;

  // Create new clik object for low l likelihood
  lo_l_clik = new ClikObject(lo_l_clik_path);

  lo_l_nuis_enums.push_back(ClikPar::A_planck);

  lo_l_clik->set_nuisance_param_enums(lo_l_nuis_enums);

  // Print out nuisance parameter names for the world to see
  std::cout << *lo_l_clik;
  //*/


  //*
  // Package it all together
  plc_pack = new PLCPack();
  plc_pack->add_clik_object(hi_l_clik);
  plc_pack->add_clik_object(lo_l_clik);
  plc_pack->set_clik_params(clik_par);
  // plc_pack->read_pbh_files(pbh_file_root);

  // Initialise CLASS before runing MultiNest
  plc_pack->initialise_CLASS();

  context = plc_pack;
  //*/

  // Diver parameters
  int nPar = ClikPar::FREE_PARAM_AMT;
  int NP = 10*nPar;
  double* lowerbounds = clik_par->get_lower_bounds();
  double* upperbounds = clik_par->get_upper_bounds();
  int nDerived = ClikPar::DERIVED_PARAM_AMT;

  // Objective function parameters
  const int param_dim = nPar + nDerived;
  double params[param_dim] = {
    0.022252,
    0.11987,
    1.040778,
    0.0789,
    3.0929,
    0.96475,
    1.00029,
    66.4,
    0.13,
    7.17,
    255.0,
    40.1,
    36.4,
    98.7,
    0.00,
    7.34,
    8.97,
    17.56,
    81.9,
    0.0813,
    0.0488,
    0.0995,
    0.1002,
    0.2236,
    0.645,
    0.1417,
    0.1321,
    0.307,
    0.155,
    0.338,
    1.667,
    0.99818,
    0.99598
  };
  int fcall = 0;
  bool quit = false;
  const bool validvector = true;

  double loglike = 0.0;

  double norm_params[param_dim];

  for (int i = 0; i < nPar; ++i) {
    norm_params[i] = (params[i] - lowerbounds[i])/(upperbounds[i] - lowerbounds[i]);
    std::cout << params[i] << " = " << lowerbounds[i] << " + " << norm_params[i] << " * ( " << upperbounds[i] << " - " << lowerbounds[i] << " )" << std::endl;
  }

  for (int i = 0; i < nDerived; ++i) {
    norm_params[nPar+i] = params[nPar+i];
  }

  std::cout << "Running loglikelihood function..." << std::endl;

  loglike = objective(norm_params, param_dim, fcall, quit, validvector, context); 

  delete plc_pack;

  return 0;
}
