#include "init_plc.h"

#include <array>
#include <iostream>
#include <vector>

double m_min[FREE_PARAM_AMT], m_max[FREE_PARAM_AMT];
double m_value[FIXED_PARAM_AMT];

bool m_has_gaussian_prior[TOTAL_PARAM_AMT];
double m_mean[TOTAL_PARAM_AMT], m_stddev[TOTAL_PARAM_AMT];

bool m_is_log10[FREE_PARAM_AMT];
trans_t m_transform[FREE_PARAM_AMT];

void init_CLASS(ClassEngine*& class_engine, int max_l, external_info* info);
void update_CLASS(ClassEngine*& class_engine);

void free_bspline(struct bspline_2d* spline);

int main(int argc, char const *argv[])
{
  ////// plc_class variables //////

  int total_max_l = 2508;

  std::vector<clik_struct*> clik_objects;
  plc_bundle* plc_pack = new plc_bundle();

  // PBH variables
  std::string pbh_file_root = std::string(CLASS_PBH_FILE_DIR) + "/";
  std::string hyrec_file_root = std::string(HYREC_FILE_DIR) + "/";
  external_info* info;

  /// Actual Program ///

  // Read in external PBH files
  info = initialise_external_info(pbh_file_root, hyrec_file_root);

  // Create cl_ls vector of l values!!!
  for (int l = CLASS_MIN_L; l <= total_max_l; ++l) {
    plc_pack->cl_ls.push_back(l);
  }

  try {
    init_CLASS(plc_pack->engine, total_max_l, info);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS initialisation "
              << "failed, throwing exception "
              << e.what()
              << std::endl;
    return -1;
  }

  update_CLASS(plc_pack->engine);

  /* Free all heap-allocated memory at the end */
  delete plc_pack->engine;
  delete plc_pack;

  delete[] info->masses;
  delete[] info->z_deps;
  free_bspline(info->hion);
  free_bspline(info->excite);
  free_bspline(info->heat);

  hyrec_free_2D_array(NTM, info->logAlpha_tab[0]);
  hyrec_free_2D_array(NTM, info->logAlpha_tab[1]);

  delete info;

  std::cout << "CLASS tested, exiting..." << std::endl;

  return 0;
}

void init_CLASS(ClassEngine*& class_engine, int max_l, external_info* info) {
  ClassParams default_params;

  /*** FREE PARAMETERS ***/

  // PBH DM
  default_params.add("Omega_pbh_ratio", 1.E-7);
  default_params.add("pbh_mass_mean", 1.E6);

  /*** FREE/FIXED PARAMETERS ***/

  // LCDM variables
  default_params.add("omega_b", 0.022252);
  default_params.add("omega_cdm", 0.11987);
  default_params.add("100*theta_s", 1.040778);
  default_params.add("tau_reio", 0.0789);
  default_params.add("ln10^{10}A_s", 3.0929);
  default_params.add("n_s", 0.96475);

  // PBH DM
  default_params.add("pbh_mass_dist", "pbh_delta");
  // default_params.add("pbh_mass_width", 1.E1);
  default_params.add("read external files", false); // very important!!

  // Neutrino values
  default_params.add("N_ur", 2.0328);
  default_params.add("N_ncdm", 1);
  default_params.add("m_ncdm", 0.06); // MeV
  // default_params.add("T_ncdm", 0.71611);

  // Use HyRec for recombination
  default_params.add("recombination", "HyRec");

  // Perturbation options for matter perturbation spectrum mPk
  // default_params.add("P_k_max_h/Mpc", 1.);
  // default_params.add("k_pivot", 0.05); // Mpc-1

  // Spectra output options
  default_params.add("output", "tCl,pCl,lCl"); // mPk
  default_params.add("lensing", true);
  default_params.add("l_max_scalars", max_l);
  default_params.add("format", "camb");

#ifdef DBUG
  default_params.add("input_verbose", 3);
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

  std::cout << "Intialising CLASS..." << std::endl;

  // Initialise CLASS engine with external PBH info
  class_engine = new ClassEngine(default_params, info);
}

/* Just updates CLASS to test values of pbh_frac and pbh_mass that break CLASS in preparation for leak checking */
void update_CLASS(ClassEngine*& class_engine) {
  std::vector<std::array<double, 2> > param_vector;

  /* TR > TR_MAX in hydrogen.c:interpolate_rates() */
  param_vector.push_back({0.940561, 643503});
  param_vector.push_back({0.454753, 300062});

  /* reionization exceeds max iterations while searching for optical depth */
  param_vector.push_back({0.592319, 3.33219e+06});
  param_vector.push_back({0.345466, 863279});

  for (size_t i = 0; i < param_vector.size(); ++i) {
    std::array<double, 2> vals = param_vector[i];
    std::vector<double> class_params(&vals[0], &vals[0]+2);

    std::cout << "Updating CLASS with following parameters:\n"
              << "pbh_frac = " << class_params[0] << " "
              << "pbh_mass = " << class_params[1]
              << std::endl;

    try {
      class_engine->update_parameters(class_params);
    }
    catch (std::exception const& e) {
      std::cerr << "[ERROR] CLASS failed, throwing exception "
                << e.what()
                << std::endl;
    }
  }
}

void free_bspline(struct bspline_2d* spline) {
  delete[] spline->xknots;
  delete[] spline->yknots;
  delete[] spline->coeffs;
  delete spline;
}
