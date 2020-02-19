#include "profile.hpp"

#include "loglike.hpp"
#include "external.hpp"

#include "hyrec_io.hpp"
#include "io_params.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>  // iota
#include <stdexcept>
#include <vector>

using namespace std;

void
profile_likelihood(string& root, param_t profiled_param,
                   pair<double, double>& range, double fixed_values[],
                   unsigned num_samples, void* context) {
  ofstream fout(root.c_str(), ofstream::out);
  fout.setf(ios::scientific);
  fout.precision(16);

  for (unsigned i = 0; i < num_samples; ++i) {
    trans_t func = g_transform[profiled_param];
    double profile_value = func(range.first + (range.second - range.first) * (double) i/(num_samples - 1.0));

    fixed_values[(unsigned) profiled_param] = profile_value;

    double lnew = calculate_full_likelihood(fixed_values, static_cast<likelihood_context*>(context));

    cout << "[profile_likelihood] Calculated loglike of " << lnew << endl;

    fout << setw(26) << profile_value
         << setw(26) << lnew << endl;
  }

  fout.close();
}

int main(int argc, char const *argv[])
{
  string output_file = "output/profile.dat";
  pair<double, double> Mpbh_range(5.0, 7.0);
  unsigned num_samples = 101;

  if (argc == 2) {
    output_file = argv[1];
  }

  cout << "[INFO] Writing output to '" << output_file << "'..." << endl;

  // Low-l clik parameters
  const string lo_l_clik_file = string(PLIK_LOW_L_FILE_DIR) \
  + "/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/";
  vector<param_t> lo_params = {A_planck};
  ClikObject* lo_l_clik(0);

  // High-l clik parameters
#ifdef LITE_HI_L
  const string hi_l_clik_file = string(PLIK_HI_L_FILE_DIR) \
  + "/plik_lite_v18_TTTEEE.clik/";
#else
  const string hi_l_clik_file = string(PLIK_HI_L_FILE_DIR) \
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

  params.add("pbh_mass_mean", 1e6);
  params.add("Omega_pbh_ratio", 1e-7);
  params.add("pbh_mass_width", 1e-2);

  params.add("omega_b", 0.022220);
  params.add("omega_cdm", 0.11983);
  params.add("100*theta_s", 1.041788);
  params.add("tau_reio", 0.0803);
  params.add("ln10^{10}A_s", 3.0957);
  params.add("n_s", 0.96352);

  params.add("N_ur", 2.0328);
  params.add("N_ncdm", 1);
  params.add("m_ncdm", 0.06); // MeV

  params.add("recombination", "HyRec");

  params.add("pbh_mass_dist", "pbh_log_norm");
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

  // Run profiler
  cout << "[INFO] Starting profiler..." << endl;
  profile_likelihood(output_file, pbh_mass, Mpbh_range, (double[]){0.0}, num_samples, (void*) bundle);

  cout << "[INFO] Cleaning and exiting..." << endl;

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
