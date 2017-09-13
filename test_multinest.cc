#include "pc_multinest.h"

#include <iostream>
#include <vector>

double m_min[FREE_PARAM_AMT], m_max[FREE_PARAM_AMT];
double m_value[FIXED_PARAM_AMT];

bool m_has_gaussian_prior[TOTAL_PARAM_AMT];
double m_mean[TOTAL_PARAM_AMT], m_stddev[TOTAL_PARAM_AMT];

bool m_is_log10[FREE_PARAM_AMT];
trans_t m_transform[FREE_PARAM_AMT];

void init_CLASS(ClassEngine*& class_engine, int max_l, pbh_external* pbh_info);

int main(int argc, char const *argv[])
{
  ////// plc_class variables //////

  int total_max_l = 2508;

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
  std::string pbh_file_root = std::string(CLASS_PBH_FILE_DIR) + "/";
  pbh_external* pbh_info;

  /// Actual Program ///

  // Read in external PBH files
  pbh_info = initialise_pbh_external(pbh_file_root);

  // Create cl_ls vector of l values!!!
  for (int l = CLASS_MIN_L; l <= total_max_l; ++l) {
    plc_pack->cl_ls.push_back(l);
  }

  init_CLASS(plc_pack->engine, total_max_l, pbh_info);

  return 0;
}

void init_CLASS(ClassEngine*& class_engine, int max_l, pbh_external* pbh_info) {
  ClassParams default_params;

  /*** FREE PARAMETERS ***/

  // PBH DM
  default_params.add("Omega_pbh_ratio", 5.68187e-06);

  /*** FREE/FIXED PARAMETERS ***/

  // LCDM variables
  default_params.add("omega_b", 0.0194805);
  default_params.add("omega_cdm", 0.128124);
  default_params.add("100*theta_s", 1.04008);
  default_params.add("tau_reio", 0.0290041);
  default_params.add("ln10^{10}A_s", 3.04091);
  default_params.add("n_s", 1.00551);

  // PBH DM
  default_params.add("pbh_mass_dist", "pbh_delta");
  default_params.add("pbh_mass_mean", 1E5);
  // default_params.add("pbh_mass_width", 1.E1);
  default_params.add("read pbh splines", false); // very important!!

  // Neutrino values
  default_params.add("N_ur", 2.0328);
  default_params.add("N_ncdm", 1);
  default_params.add("m_ncdm", 0.06); // MeV
  // default_params.add("T_ncdm", 0.71611);

  // Perturbation options for matter perturbation spectrum mPk
  // default_params.add("P_k_max_h/Mpc", 1.);
  // default_params.add("k_pivot", 0.05); // Mpc-1

  // Spectra output options
  default_params.add("output", "tCl,pCl,lCl"); // mPk
  default_params.add("lensing", true);
  default_params.add("l_max_scalars", max_l);
  default_params.add("format", "camb");

#ifdef DBUG
  default_params.add("input_verbose", 1);
  default_params.add("background_verbose", 1);
  default_params.add("thermodynamics_verbose", 1);
  default_params.add("perturbations_verbose", 1);
  default_params.add("transfer_verbose", 1);
  default_params.add("primordial_verbose", 1);
  default_params.add("spectra_verbose", 1);
  default_params.add("nonlinear_verbose", 1);
  default_params.add("lensing_verbose", 1);
  default_params.add("output_verbose", 1);
#endif

  // Initialise CLASS engine with external PBH info
  try {
    class_engine = new ClassEngine(default_params, pbh_info);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS initialisation "
              << "failed, throwing exception "
              << e.what()
              << std::endl;
    throw e;
  }
}
