#include "ClikObject.hpp"

#include <iostream>
#include <stdexcept>

using namespace std;

int main(int argc, char const *argv[])
{
  const string lo_l_clik_file = string(PLIK_LOW_L_FILE_DIR) \
  + "/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik/";
  vector<param_t> lo_params = {A_planck};
  ClikObject* lo_l_clik(0);

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

  try {
    lo_l_clik = new ClikObject(lo_l_clik_file, lo_params);
    hi_l_clik = new ClikObject(hi_l_clik_file, hi_params);
  }
  catch (exception& e) {
    cerr << "[ERROR] " << e.what() << endl;
    return -1;
  }

  delete lo_l_clik;
  delete hi_l_clik;

  return 0;
}
