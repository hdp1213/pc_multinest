#include "init_plc.hpp"

#include "pbh_io.hpp"
#include "hyrec_io.hpp"

#include <cstdio> // for stderr
#include <exception> // for std::exception
#include <iostream> // you know, for kids


clik_struct* initialise_clik_struct(std::string& clik_path,
                                    std::vector<param_t>& nuis_params,
                                    int& total_max_l) {
  std::cout << "[init_plc] initialising clik_struct"
            << std::endl;

  clik_struct* res_struct = new clik_struct();

  int cl_flags[CL_AMT];
  int cl_max_ls[CL_AMT];
  parname* nuis_param_names;


  error* _err = initError();
  error** err = &_err;
  char* clik_path_c = const_cast<char*>(clik_path.c_str());

  // Initialise clik_id
  res_struct->clik_id = clik_init(clik_path_c, err);
  quitOnError(*err, __LINE__, stderr);

  // Get flags indicating what spectra are contained in .clik file:
  //   TT EE BB TE TB EB
  clik_get_has_cl(res_struct->clik_id, cl_flags, err);
  quitOnError(*err, __LINE__, stderr);

  // Get array of maximum multipole for each existing spectra
  clik_get_lmax(res_struct->clik_id, cl_max_ls, err);
  quitOnError(*err, __LINE__, stderr);

  // !!!! Assuming compatibility with CLASS !!!!

  res_struct->max_l = cl_max_ls[0];

  if (res_struct->max_l > total_max_l) {
    total_max_l = res_struct->max_l;
  }

  res_struct->nuis_par_amt = clik_get_extra_parameter_names(res_struct->clik_id, &nuis_param_names, err);
  quitOnError(*err, __LINE__, stderr);

  // Calculate size of cl_and_pars array
  res_struct->cap_size = 0;
  for (int cl = 0; cl < CL_AMT; ++cl) {
    res_struct->has_cl[cl] = (cl_flags[cl] == 1);

    if (res_struct->has_cl[cl]) {
      res_struct->cap_size += cl_max_ls[cl] + 1; // +1 for l=0 case
    }
  }
  res_struct->cap_size += res_struct->nuis_par_amt;

  if (nuis_params.size() > res_struct->nuis_par_amt) {
    std::cerr << "[ERROR] More nuisance parameters given to "
              << "initialise_clik_struct() than needed!"
              << std::endl;
    throw std::exception();
  }
  else if (nuis_params.size() < res_struct->nuis_par_amt) {
    std::cerr << "[ERROR] Not enough nuisance parameters have been "
              << "given to initialise_clik_struct()!"
              << std::endl;
    throw std::exception();
  }

  res_struct->nuis_pars.reserve(res_struct->nuis_par_amt);
  res_struct->nuis_pars = nuis_params;

  free(_err);
  free(nuis_param_names);

  return res_struct;
}

void initialise_CLASS_engine(ClassEngine*& class_engine, int max_l, external_info* info) {
  std::cout << "[init_plc] initialising CLASS engine"
            << std::endl;

  ClassParams default_params;

  /*** FREE PARAMETERS ***/

  // LCDM variables set to own best-fit TTTEEE+lowTEB
  default_params.add("omega_b", 0.022231);
  default_params.add("omega_cdm", 0.12003);
  default_params.add("100*theta_s", 1.041740);
  default_params.add("tau_reio", 0.0807);
  default_params.add("ln10^{10}A_s", 3.0933);
  default_params.add("n_s", 0.96749);

  // PBH DM
  default_params.add("pbh_mass_dist", "pbh_none");
  default_params.add("read external files", false); // very important!!

  // Neutrino values
  default_params.add("N_ur", 3.046);

  // Use HyRec for recombination
  default_params.add("recombination", "RECFAST");

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
    class_engine = new ClassEngine(default_params, info);
    // class_engine = new ClassEngine(default_params);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS initialisation "
              << "failed, throwing exception "
              << e.what()
              << std::endl;
    throw e;
  }
}

external_info* initialise_external_info(std::string& pbh_root, std::string& hyrec_root) {
  std::cout << "[init_plc] initialising external_info struct"
            << std::endl;

  /* First read in the PBH stuff */

  external_info* info = new external_info();
  info->hion = new bspline_2d();
  info->excite = new bspline_2d();
  info->heat = new bspline_2d();

  // Read in axes
  read_axes(pbh_root, info);

  // Read in hydrogen ionisation
  read_bicubic_bspline(pbh_root, "hion.dat", info->hion);

  // Read in hydrogen excitation
  read_bicubic_bspline(pbh_root, "excite.dat", info->excite);

  // Read in plasma heating
  read_bicubic_bspline(pbh_root, "heat.dat", info->heat);

  /* Then read in the HyRec stuff */

  // This is good to do, too
  initialise_temp_tables(info);

  read_alpha(hyrec_root, info);

  read_RR(hyrec_root, info);

  read_two_photon(hyrec_root, info);

  // Normalise the stuff
  normalise_atomic(info);

  return info;
}

void initialise_param_arrays() {
  std::cout << "[init_plc] initialising parameter arrays"
            << std::endl;

  // Set defaults for global arrays
  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    m_is_log10[param] = false;
  }

  // Include parameter array initialisations

  // #include "TTTEEE+lowP_pbh_fixedLCDM-flat.cpp"
  // #include "TTTEEE+lowP_pbh_clark-flat.cpp"
  // #include "TTTEEE+lowP_pbh_full-flat.cpp"
  // #include "TTTEEE+lowP_pbh_full_freeLCDM-flat.cpp"
  // #include "TTTEEE+lowP_pbh_clark_freeLCDM-flat.cpp"
  // #include "TTTEEE+lowP_pbh_dist-flat.cpp"
  // #include "TTTEEE+lowP_pbh-flat.cpp"
  #include "TTTEEE+lowP-flat.cpp"

  // #include "TTTEEE+lowP_pbh_fixedLCDM-gauss.cpp"
  // #include "TTTEEE+lowP_pbh_clark-gauss.cpp"
  // #include "TTTEEE+lowP_pbh_full-gauss.cpp"
  // #include "TTTEEE+lowP_pbh_full_freeLCDM-gauss.cpp"
  // #include "TTTEEE+lowP_pbh_clark_freeLCDM-gauss.cpp"
  // #include "TTTEEE+lowP_pbh_dist-gauss.cpp"
  // #include "TTTEEE+lowP_pbh-gauss.cpp"
  #include "TTTEEE+lowP-gauss.cpp"

  // Set values of m_transform depending on includes
  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    m_transform[param] = m_is_log10[param] ? &pow10 : &self;
  }
}
