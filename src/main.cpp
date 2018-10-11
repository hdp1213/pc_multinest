#include "main.hpp"
#include "multinest_loglike.hpp"
#include "loglike.hpp"
#include "external.hpp"

#include "hyrec_io.hpp"
#include "io_params.h"

#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using namespace std;

int main(int argc, char const *argv[])
{
  ////// MultiNest settings //////

  multinest_settings settings;

  settings.IS = false;
  settings.mmodal = false;
  settings.ceff = true;

  settings.nlive = 1000;
  settings.efr = 0.05;
  settings.tol = 1E-1;

  settings.ndims = FREE_PARAM_AMT;
  settings.nPar = FREE_PARAM_AMT + DERIVED_PARAM_AMT;
  settings.nClsPar = 0;

  settings.updInt = 1000;
  settings.Ztol = MIN_LOGLIKE;
  settings.maxModes = 4;
  settings.pWrap = new int[settings.ndims];
  for(int i = 0; i < settings.ndims; i++) settings.pWrap[i] = 0;

  settings.root = "output/base_planck-";
  settings.seed = -1;
  settings.fb = true;
  settings.resume = true;
  settings.outfile = true;
  settings.initMPI = true;
  settings.logZero = MIN_LOGLIKE;
  settings.maxiter = 0;

  if (argc == 2) {
    settings.root = argv[1];
  }

  cout << "[INFO] Writing output to '" << settings.root << "'..." << endl;

  // Low-l clik parameters
  const string lo_l_clik_file = string(PLIK_LOW_L_FILE_DIR) \
  + "/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/";
  vector<param_t> lo_params = {A_planck};
  ClikObject* lo_l_clik(0);

  // High-l clik parameters
#ifdef LITE_HI_L
  const std::string hi_l_clik_file = std::string(PLIK_HI_L_FILE_DIR) \
  + "/plik_lite_v18_TTTEEE.clik/";
#else
  const std::string hi_l_clik_file = std::string(PLIK_HI_L_FILE_DIR) \
  + "/plik_dx11dr2_HM_v18_TTTEEE.clik/";
#endif
  vector<param_t> hi_params = {
#ifndef LITE_HI_L
    A_cib_217,
    cib_index,
    xi_sz_cib,
    A_sz,
    ps_A_100_100,
    ps_A_143_143,
    ps_A_143_217,
    ps_A_217_217,
    ksz_norm,
    gal545_A_100,
    gal545_A_143,
    gal545_A_143_217,
    gal545_A_217,
    galf_EE_A_100,
    galf_EE_A_100_143,
    galf_EE_A_100_217,
    galf_EE_A_143,
    galf_EE_A_143_217,
    galf_EE_A_217,
    galf_EE_index,
    galf_TE_A_100,
    galf_TE_A_100_143,
    galf_TE_A_100_217,
    galf_TE_A_143,
    galf_TE_A_143_217,
    galf_TE_A_217,
    galf_TE_index,
    bleak_epsilon_0_0T_0E,
    bleak_epsilon_1_0T_0E,
    bleak_epsilon_2_0T_0E,
    bleak_epsilon_3_0T_0E,
    bleak_epsilon_4_0T_0E,
    bleak_epsilon_0_0T_1E,
    bleak_epsilon_1_0T_1E,
    bleak_epsilon_2_0T_1E,
    bleak_epsilon_3_0T_1E,
    bleak_epsilon_4_0T_1E,
    bleak_epsilon_0_0T_2E,
    bleak_epsilon_1_0T_2E,
    bleak_epsilon_2_0T_2E,
    bleak_epsilon_3_0T_2E,
    bleak_epsilon_4_0T_2E,
    bleak_epsilon_0_1T_1E,
    bleak_epsilon_1_1T_1E,
    bleak_epsilon_2_1T_1E,
    bleak_epsilon_3_1T_1E,
    bleak_epsilon_4_1T_1E,
    bleak_epsilon_0_1T_2E,
    bleak_epsilon_1_1T_2E,
    bleak_epsilon_2_1T_2E,
    bleak_epsilon_3_1T_2E,
    bleak_epsilon_4_1T_2E,
    bleak_epsilon_0_2T_2E,
    bleak_epsilon_1_2T_2E,
    bleak_epsilon_2_2T_2E,
    bleak_epsilon_3_2T_2E,
    bleak_epsilon_4_2T_2E,
    bleak_epsilon_0_0E_0E,
    bleak_epsilon_1_0E_0E,
    bleak_epsilon_2_0E_0E,
    bleak_epsilon_3_0E_0E,
    bleak_epsilon_4_0E_0E,
    bleak_epsilon_0_0E_1E,
    bleak_epsilon_1_0E_1E,
    bleak_epsilon_2_0E_1E,
    bleak_epsilon_3_0E_1E,
    bleak_epsilon_4_0E_1E,
    bleak_epsilon_0_0E_2E,
    bleak_epsilon_1_0E_2E,
    bleak_epsilon_2_0E_2E,
    bleak_epsilon_3_0E_2E,
    bleak_epsilon_4_0E_2E,
    bleak_epsilon_0_1E_1E,
    bleak_epsilon_1_1E_1E,
    bleak_epsilon_2_1E_1E,
    bleak_epsilon_3_1E_1E,
    bleak_epsilon_4_1E_1E,
    bleak_epsilon_0_1E_2E,
    bleak_epsilon_1_1E_2E,
    bleak_epsilon_2_1E_2E,
    bleak_epsilon_3_1E_2E,
    bleak_epsilon_4_1E_2E,
    bleak_epsilon_0_2E_2E,
    bleak_epsilon_1_2E_2E,
    bleak_epsilon_2_2E_2E,
    bleak_epsilon_3_2E_2E,
    bleak_epsilon_4_2E_2E,
    calib_100T,
    calib_217T,
    calib_100P,
    calib_143P,
    calib_217P,
    A_pol,
#endif
    A_planck
  };
  ClikObject* hi_l_clik(0);

  cout << "[INFO] Initialising ClikObjects..." << endl;

  // Clik object initialisation
  try {
    lo_l_clik = new ClikObject(lo_l_clik_file, lo_params);
    hi_l_clik = new ClikObject(hi_l_clik_file, hi_params);
  }
  catch (exception& e) {
    cerr << "[ERROR] " << e.what() << endl;
    return -1;
  }

  // Initialise global parameter arrays
  initialise_param_arrays();

  // Initialise external info
  string pbh_root = string(CLASS_PBH_FILE_DIR) + "/";
  string hyrec_root = string(HYREC_FILE_DIR) + "/";
  external_info* info = initialise_external_info(pbh_root, hyrec_root);

  // Allocate bundle
  likelihood_context* bundle = new likelihood_context();

  bundle->clik_objs.push_back(lo_l_clik);
  bundle->clik_objs.push_back(hi_l_clik);

  // Calculate lmax
  unsigned lmax = 0;

  for (auto clik_obj : bundle->clik_objs) {
    lmax = max(lmax, clik_obj->get_lmax());
  }

  // Initialise CLASS
  ClassParams params;

  params.add("Omega_pbh_ratio", 1e-90);
  params.add("pbh_mass_mean", 1e6);

  params.add("omega_b", 0.022231);
  params.add("omega_cdm", 0.12003);
  params.add("100*theta_s", 1.041740);
  params.add("tau_reio", 0.0807);
  params.add("ln10^{10}A_s", 3.0933);
  params.add("n_s", 0.96749);

  params.add("N_ur", 2.0328);
  params.add("N_ncdm", 1);
  params.add("m_ncdm", 0.06); // MeV

  params.add("recombination", "HyRec");

  params.add("pbh_mass_dist", "pbh_delta");
  params.add("read external files", false);

  params.add("output", "tCl,pCl,lCl");
  params.add("lensing", true);
  params.add("l_max_scalars", lmax);
  params.add("format", "camb");

  cout << "[INFO] Initialising CLASS engine..." << endl;

  try {
    bundle->engine = new ClassEngine(params, info);
  }
  catch (exception const& e) {
    cerr << "[ERROR] " << e.what() << endl;
    return -1;
  }

  // Construct Cl multipoles
  bundle->cl_ls.resize(lmax - CLASS_MIN_L + 1);
  iota(bundle->cl_ls.begin(), bundle->cl_ls.end(), CLASS_MIN_L);

  // Run MultiNest
  pc_multinest(multinest_loglike, multinest_dumper, settings, (void*) bundle);

  cout << "[INFO] Cleaning and exiting..." << endl;

  delete[] settings.pWrap;

  delete bundle->engine;
  for (auto clik_obj : bundle->clik_objs) {
    delete clik_obj;
  }

  delete bundle;

  hyrec_free_2D_array(NTM, info->logAlpha_tab[0]);
  hyrec_free_2D_array(NTM, info->logAlpha_tab[1]);

  delete info;

  return 0;
}
