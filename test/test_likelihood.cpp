#include "ClassEngine.hpp"
#include "ClikObject.hpp"
#include "loglike.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

const double LIKE_VALUE = -6.47490895946E3;
const double TOL = 5e-12;

using namespace std;

int
main(int argc, char const *argv[])
{
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

  initialise_param_arrays();

  // Clik object initialisation
  try {
    lo_l_clik = new ClikObject(lo_l_clik_file, lo_params);
    hi_l_clik = new ClikObject(hi_l_clik_file, hi_params);
  }
  catch (exception& e) {
    cerr << "[ERROR] " << e.what() << endl;
    return -1;
  }

  // Allocate bundle
  likelihood_context* bundle = new likelihood_context();

  bundle->clik_objs.push_back(lo_l_clik);
  bundle->clik_objs.push_back(hi_l_clik);

  // Calculate lmax
  unsigned lmax = 0;

  for (auto clik_obj : bundle->clik_objs) {
    lmax = max(lmax, clik_obj->get_lmax());
  }

  // CLASS parameters
  ClassParams params;

  params.add("omega_b", 0.022231);
  params.add("omega_cdm", 0.12003);
  params.add("100*theta_s", 1.041740);
  params.add("tau_reio", 0.0807);
  params.add("ln10^{10}A_s", 3.0933);
  params.add("n_s", 0.96749);

  params.add("N_ur", 2.0328);
  params.add("N_ncdm", 1);
  params.add("m_ncdm", 0.06); // MeV

  params.add("output", "tCl,pCl,lCl");
  params.add("lensing", true);
  params.add("l_max_scalars", lmax);
  params.add("format", "camb");

  params.add("pbh_root", "/home/a1648400/class/pbh/pbh_bspline_");

  // Initialise Cuuube
  double Cube[] = {
    // CLASS free variables
    0.022252,
    0.11987,
    1.040778,
    0.0789,
    3.0929,
    0.96475,
    // PLC free values
    1.00029,
#ifndef LITE_HI_L
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
    0.99598,
#endif
    // Derived variables (to set)
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };

  try {
    bundle->engine = new ClassEngine(params);
  }
  catch (exception const& e) {
    cerr << "[ERROR] " << e.what() << endl;
    return -1;
  }

  bundle->cl_ls.resize(lmax - CLASS_MIN_L + 1);
  iota(bundle->cl_ls.begin(), bundle->cl_ls.end(), CLASS_MIN_L);

  double full_like = calculate_full_likelihood(Cube, bundle);

  if (full_like != full_like) {
    cerr << "[ERROR] Your likelihood value is NaN!" << endl;
  }

  double diff = fabs((LIKE_VALUE - full_like)/LIKE_VALUE);
  cout << "Diff: " << diff << endl;

  if (diff > TOL) {
    cerr << "[ERROR] Difference above tolerance" << endl;
  } 

  cout << "H0 = " << Cube[H0_LOC] << endl;
  cout << "Omega_b = " << Cube[Omega_b_LOC] << endl;
  cout << "Omega_cdm = " << Cube[Omega_cdm_LOC] << endl;
  cout << "Omega_L = " << Cube[Omega_L_LOC] << endl;
  cout << "Omega_g = " << Cube[Omega_g_LOC] << endl;
  cout << "sigma8 = " << Cube[sigma8_LOC] << endl;
  cout << "age = " << Cube[age_LOC] << endl;
  cout << "conf_age = " << Cube[conf_age_LOC] << endl;
  cout << "z_drag = " << Cube[z_drag_LOC] << endl;
  cout << "rs_drag = " << Cube[rs_drag_LOC] << endl;

  delete bundle->engine;
  for (auto clik_obj : bundle->clik_objs) {
    delete clik_obj;
  }

  delete bundle;

  return 0;
}
