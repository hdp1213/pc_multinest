#include "multinest_loglike.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

double m_min[FREE_PARAM_AMT], m_max[FREE_PARAM_AMT];
double m_value[FIXED_PARAM_AMT];

bool m_has_gaussian_prior[TOTAL_PARAM_AMT];
double m_mean[TOTAL_PARAM_AMT], m_stddev[TOTAL_PARAM_AMT];

bool m_is_log10[FREE_PARAM_AMT];
trans_t m_transform[FREE_PARAM_AMT];

int main(int argc, char* argv[]) {

  ////// plc_class variables //////

  void* context = 0;
  int total_max_l = -1;

  double log_pbh_mass, pbh_mass, log_pbh_width, pbh_width;

  std::vector<clik_struct*> clik_objects;
  plc_bundle* plc_pack = new plc_bundle();

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

  // PBH variables
  std::string pbh_file_root = std::string(CLASS_PBH_FILE_DIR);
  std::string hyrec_file_root = std::string(HYREC_FILE_DIR) + "/";
  external_info* info;

  // Read in pbh width and mass
  if (argc == 3) {
    std::istringstream arg1stream(argv[1]);
    arg1stream >> log_pbh_mass;

    std::istringstream arg2stream(argv[2]);
    arg2stream >> log_pbh_width;
  }
  else {
    throw std::invalid_argument("Not enough parameters given");
  }

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

  lo_l_nuis_enums.push_back(A_planck);

  std::cout << "Opening " << lo_l_clik_path << std::endl;

  // Create new clik object for low l likelihood
  lo_l_clik = initialise_clik_struct(lo_l_clik_path,
                                     lo_l_nuis_enums,
                                     total_max_l);
  plc_pack->clik_objs.push_back(lo_l_clik);

  // Read in external PBH files
  info = initialise_external_info(pbh_file_root, hyrec_file_root);

  // Create cl_ls vector of l values
  for (int l = CLASS_MIN_L; l <= total_max_l; ++l) {
    plc_pack->cl_ls.push_back(l);
  }

  // Initialise CLASS
  initialise_CLASS_engine(plc_pack->engine, total_max_l, info);

  context = plc_pack;

  initialise_param_arrays();

  /* Run the dang thing */

  std::vector<double> new_pars;

  double Aplanck = 1.00029;
  pbh_mass = pow(10., log_pbh_mass);
  pbh_width = pow(10., log_pbh_width);

  new_pars.push_back(pbh_mass);
  new_pars.push_back(pbh_width);
  new_pars.push_back(Aplanck);

  int ndim = FREE_PARAM_AMT;
  int npar = FREE_PARAM_AMT + DERIVED_PARAM_AMT;
  double lnew = 1E-99;

  multinest_loglike_pt(new_pars, ndim, npar, lnew, context);

  /* Write the dang output to a dang file */
  std::ostringstream out_name;

  out_name << "point_results/"
           << std::fixed << std::setprecision(3) << log_pbh_mass
           << "_"
           << std::fixed << std::setprecision(3) << log_pbh_width << ".res";

  std::cout << "Printing result to " << out_name.str() << "..." << std::endl;

  std::ofstream file_out(out_name.str().c_str(), std::ofstream::out);

  if (!file_out.is_open()) {
    std::cerr << "[ERROR]: file output error" << std::endl;
    return -1;
  }

  file_out << std::fixed << std::setw(16) << std::setprecision(2)
           << log_pbh_mass
           << std::setw(16) << std::setprecision(2)
           << log_pbh_width
           << std::scientific << std::setw(16) << std::setprecision(6)
           << lnew
           << std::endl;

  file_out.close();

  return 0;
}
