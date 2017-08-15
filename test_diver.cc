#include <limits>

#include <vector>
#include <iostream> // you know, for kids

// double diver_loglike(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector, void*& context);

#include "diver_loglike.h"

int main(int argc, char** argv)
{
  void* context = 0;
  int total_max_l = -1;

  std::vector<clik_struct*> clik_objects;
  diver_bundle* plc_pack = new diver_bundle();

  // High l full likelihood variables
  std::string hi_l_clik_path = std::string(PLIK_HI_L_FILE_DIR);
#ifdef LITE_HI_L
  hi_l_clik_path += "/plik_lite_v18_TTTEEE.clik/";
#else
  hi_l_clik_path += "/plik_dx11dr2_HM_v18_TTTEEE.clik/";
#endif
  std::vector<param_t> hi_l_nuis_enums;
  clik_struct* hi_l_clik;

  // Low l likelihood variables  
  std::string lo_l_clik_path = std::string(PLIK_LOW_L_FILE_DIR) \
    + "/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/";
  std::vector<param_t> lo_l_nuis_enums;
  clik_struct* lo_l_clik;


  //*
  // Push nuisance parameters in the order they appear in cl_and_pars
#ifndef LITE_HI_L
// TT & TTTEEE
  hi_l_nuis_enums.push_back(A_cib_217);
  hi_l_nuis_enums.push_back(cib_index);
  hi_l_nuis_enums.push_back(xi_sz_cib);
  hi_l_nuis_enums.push_back(A_sz);
  hi_l_nuis_enums.push_back(ps_A_100_100);
  hi_l_nuis_enums.push_back(ps_A_143_143);
  hi_l_nuis_enums.push_back(ps_A_143_217);
  hi_l_nuis_enums.push_back(ps_A_217_217);
  hi_l_nuis_enums.push_back(ksz_norm);
  hi_l_nuis_enums.push_back(gal545_A_100);
  hi_l_nuis_enums.push_back(gal545_A_143);
  hi_l_nuis_enums.push_back(gal545_A_143_217);
  hi_l_nuis_enums.push_back(gal545_A_217);
// TTTEEE
  hi_l_nuis_enums.push_back(galf_EE_A_100);
  hi_l_nuis_enums.push_back(galf_EE_A_100_143);
  hi_l_nuis_enums.push_back(galf_EE_A_100_217);
  hi_l_nuis_enums.push_back(galf_EE_A_143);
  hi_l_nuis_enums.push_back(galf_EE_A_143_217);
  hi_l_nuis_enums.push_back(galf_EE_A_217);
  hi_l_nuis_enums.push_back(galf_EE_index);
  hi_l_nuis_enums.push_back(galf_TE_A_100);
  hi_l_nuis_enums.push_back(galf_TE_A_100_143);
  hi_l_nuis_enums.push_back(galf_TE_A_100_217);
  hi_l_nuis_enums.push_back(galf_TE_A_143);
  hi_l_nuis_enums.push_back(galf_TE_A_143_217);
  hi_l_nuis_enums.push_back(galf_TE_A_217);
  hi_l_nuis_enums.push_back(galf_TE_index);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_0T_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_0T_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_0T_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_0T_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_0T_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_0T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_0T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_0T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_0T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_0T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_0T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_0T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_0T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_0T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_0T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_1T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_1T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_1T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_1T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_1T_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_1T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_1T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_1T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_1T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_1T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_2T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_2T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_2T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_2T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_2T_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_0E_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_0E_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_0E_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_0E_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_0E_0E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_0E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_0E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_0E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_0E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_0E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_0E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_0E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_0E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_0E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_0E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_1E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_1E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_1E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_1E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_1E_1E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_1E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_1E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_1E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_1E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_1E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_0_2E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_1_2E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_2_2E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_3_2E_2E);
  hi_l_nuis_enums.push_back(bleak_epsilon_4_2E_2E);
// TT & TTTEEE
  hi_l_nuis_enums.push_back(calib_100T);
  hi_l_nuis_enums.push_back(calib_217T);
// TTTEEE
  hi_l_nuis_enums.push_back(calib_100P);
  hi_l_nuis_enums.push_back(calib_143P);
  hi_l_nuis_enums.push_back(calib_217P);
  hi_l_nuis_enums.push_back(A_pol);
#endif
  hi_l_nuis_enums.push_back(A_planck);

  std::cout << "Opening " << hi_l_clik_path << std::endl;

  // Create new clik object for high l likelihood
  hi_l_clik = initialise_clik_struct(hi_l_clik_path,
                                     hi_l_nuis_enums,
                                     total_max_l);
  plc_pack->clik_objs.push_back(hi_l_clik);
  //*/

  //*
  lo_l_nuis_enums.push_back(A_planck);

  std::cout << "Opening " << lo_l_clik_path << std::endl;

  // Create new clik object for low l likelihood
  lo_l_clik = initialise_clik_struct(lo_l_clik_path,
                                     lo_l_nuis_enums,
                                     total_max_l);
  plc_pack->clik_objs.push_back(lo_l_clik);
  //*/


  // Create cl_ls vector of l values!!!
  for (int l = CLASS_MIN_L; l <= total_max_l; ++l) {
    plc_pack->cl_ls.push_back(l);
  }

  //*
  // Initialise CLASS before runing MultiNest
  initialise_CLASS_engine(plc_pack->engine, total_max_l);

  context = plc_pack;
  //*/

  // Initialise m_min, m_max and the rest
  initialise_params();

  // Diver parameters
  int nPar = FREE_PARAM_AMT;
  int NP = 10*nPar;
  double* lowerbounds = m_min;
  double* upperbounds = m_max;
  int nDerived = DERIVED_PARAM_AMT;

  // Objective function parameters
  const int param_dim = nPar + nDerived;
  double params[FREE_PARAM_AMT + DERIVED_PARAM_AMT] = {
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

  double loglike;

  std::cout << "Running loglikelihood function ..."
            << std::endl;

  loglike = diver_loglike(params, param_dim, fcall, quit,
                          validvector, context);

  /*
  std::cout << "[test_diver] Calculated diver_loglike of "
            << loglike
            << std::endl;
  */

  // Deallocate memory
  for (std::vector<clik_struct*>::iterator clik_struct_it = plc_pack->clik_objs.begin(); clik_struct_it != plc_pack->clik_objs.end(); ++clik_struct_it) {
    free((*clik_struct_it)->clik_id);
    delete *clik_struct_it;
  }

  delete plc_pack->engine;
  delete plc_pack;

  return 0;
}
