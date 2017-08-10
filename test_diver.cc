#include "clik.h"
#include "ClassEngine.hh"

#include <limits>
#include <string>

// for .clik stuff
#include <cstdio> // for stderr
#include <exception> // for std::exception
#include <iostream> // you know, for kids

/***************************** GLOBALS ********************************/

// #include "TTTEEE+lowP_pbh_fixedLCDM.h"
// #include "TTTEEE+lowP_pbh.h"
#include "TTTEEE+lowP.h"

double m_min[FREE_PARAM_AMT], m_max[FREE_PARAM_AMT];
double m_value[FIXED_PARAM_AMT];

bool m_has_gaussian_prior[TOTAL_PARAM_AMT];
double m_mean[TOTAL_PARAM_AMT], m_stddev[TOTAL_PARAM_AMT];


/**********************************************************************/

const int CL_AMT = 6;
const int CLASS_CL_AMT = 4;
const int CLASS_MIN_L = 2;

struct clik_struct {
  clik_object* clik_id;
  int max_l;
  int cap_size;
  bool has_cl[CL_AMT];
  int nuis_par_amt;
  std::vector<param_t> nuis_pars;
};

struct diver_bundle {
  std::vector<clik_struct*> clik_objs;
  ClassEngine* engine;
  std::vector<unsigned> cl_ls;
};

double diver_loglike(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector, void*& context);

clik_struct* initialise_clik_struct(std::string& clik_path,
                                    std::vector<param_t>& nuis_params,
                                    int& total_max_l);
void initialise_CLASS_engine(ClassEngine*& class_engine, int max_l);
void initialise_params();


#include "pc_loglike.h"

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

  loglike = diver_loglike(params, param_dim, fcall, quit, validvector, context);

  std::cout << "[test_diver] Calculated diver_loglike of " << loglike << std::endl;

  for (std::vector<clik_struct*>::iterator clik_struct_it = plc_pack->clik_objs.begin(); clik_struct_it != plc_pack->clik_objs.end(); ++clik_struct_it) {
    delete *clik_struct_it;
  }

  delete plc_pack;

  return 0;
}




// Function to be minimized.  Corresponds to -ln(Likelihood).  Redirects to the target of context pointer.
// params contains both free and derived parameters
// at the end of likelihood computation, the latter part of the array should have derived parameters set as appropriate
double diver_loglike(double params[], const int param_dim, int &fcall, bool &quit, const bool validvector, void*& context) {

  double loglike = std::numeric_limits<double>::max();

  // If params[] is within the parameter bounds, evaluate its likelihood
  if (validvector) {
    diver_bundle* plc_pack;
    double in_vals[UP_TO_FIXED_PARAMS];

    plc_pack = static_cast<diver_bundle*>(context);

    // Copy free and fixed parameters into in_vals, the input parameters
    std::copy(params, params + FREE_PARAM_AMT, in_vals);
    std::copy(m_value, m_value + FIXED_PARAM_AMT, in_vals + FREE_PARAM_AMT);

    loglike = pc_loglike(plc_pack->clik_objs, plc_pack->engine,
                         plc_pack->cl_ls, params, in_vals);
    fcall++;
  }

  return -loglike;
}


clik_struct* initialise_clik_struct(std::string& clik_path,
                                    std::vector<param_t>& nuis_params,
                                    int& total_max_l) {
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

  return res_struct;
}

void initialise_CLASS_engine(ClassEngine*& class_engine, int max_l) {
  ClassParams default_params;

  /*** FREE PARAMETERS ***/

  // PBH DM
  // default_params.add("Omega_pbh_ratio", 1.E-7);

  /*** FREE/FIXED PARAMETERS ***/

  // LCDM variables set to best-fit TTTEEE+lowP (2015)
  default_params.add("omega_b", 0.022252);
  default_params.add("omega_cdm", 0.11987);
  default_params.add("100*theta_s", 1.040778);
  default_params.add("tau_reio", 0.0789);
  default_params.add("ln10^{10}A_s", 3.0929);
  default_params.add("n_s", 0.96475);

  // PBH DM
  // default_params.add("pbh_mass_dist", "pbh_none");
  // default_params.add("read pbh splines", false);
  /*
  default_params.add("pbh_mass_dist", "pbh_delta");
  default_params.add("pbh_mass_mean", 1.E6);
  // default_params.add("pbh_mass_width", 1.E1);
  default_params.add("read pbh splines", false); // very important!!
  */

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

  // default_params.add("input_verbose", 1);
  // default_params.add("background_verbose", 1);
  // default_params.add("thermodynamics_verbose", 1);
  // default_params.add("perturbations_verbose", 1);
  // default_params.add("transfer_verbose", 1);
  // default_params.add("primordial_verbose", 1);
  // default_params.add("spectra_verbose", 1);
  // default_params.add("nonlinear_verbose", 1);
  // default_params.add("lensing_verbose", 1);
  // default_params.add("output_verbose", 1);

  // Initialise CLASS engine with external PBH info
  try {
    // class_engine = new ClassEngine(default_params, pbh_info);
    class_engine = new ClassEngine(default_params);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS initialisation "
              << "failed, throwing exception "
              << e.what()
              << std::endl;
    throw e;
  }
}

void initialise_params() {
// #include "TTTEEE+lowP_pbh_fixedLCDM-flat.cc"
#include "TTTEEE+lowP-flat.cc"
// #include "TTTEEE+lowP_pbh_fixedLCDM-gauss.cc"
#include "TTTEEE+lowP-gauss.cc"
}
