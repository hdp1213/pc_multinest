#include "pc_multinest.h"

#include <cstdio> // for stderr
#include <iostream>
#include <string>
#include <vector>

#include "multinest.h"
#include "multinest_loglike.h"

/*
  Uses CLASS (wherever that will be with GAMBIT) to produce
  spectra which are used in the log likelihood using PLC to
  derive the log likelihood. Then MultiNest does the rest.

  5/10/16 - tested plc_class against Multinest for no free
            params.
            Shows identical agreement against plc_class for
            identical nuisance parameter input and CLASS
            initial conditions.
            Cannot find out how to use *context, though.
*/
int main(int argc, char* argv[]) {

  ////// MultiNest settings //////

  multinest_settings settings;

  settings.IS = false;
  settings.mmodal = false;
  settings.ceff = true;

  settings.nlive = 4000;
  settings.efr = 0.8;
  settings.tol = 1E-1;

  settings.ndims = FREE_PARAM_AMT;
  settings.nPar = FREE_PARAM_AMT + DERIVED_PARAM_AMT;
  settings.nClsPar = 0;

  settings.updInt = 1000;
  settings.Ztol = MIN_LOGLIKE;
  settings.maxModes = 4;
  settings.pWrap = new int[settings.ndims];
  for(int i = 0; i < settings.ndims; i++) settings.pWrap[i] = 0;

  // settings.root = 'test';
  settings.seed = -1;
  settings.fb = true;
  settings.resume = true;
  settings.outfile = true;
  settings.initMPI = true;
  settings.logZero = MIN_LOGLIKE;
  settings.maxiter = 0;


  ////// plc_class variables //////

  void* context = 0;
  int total_max_l = -1;

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

  // Use the first command line argument as a non-default
  // output root
  // That is, if an argument is even specified
  if (argc == 2) {
    settings.root = argv[1];
  }
  else {
    settings.root = "output/pc_multinest-";
  }

  std::cout << "Printing results to " << settings.root << std::endl;

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

  // Read in external PBH files
  pbh_info = initialise_pbh_external(pbh_file_root);


  // Create cl_ls vector of l values!!!
  for (int l = CLASS_MIN_L; l <= total_max_l; ++l) {
    plc_pack->cl_ls.push_back(l);
  }

  //*
  // Initialise CLASS before runing MultiNest
  initialise_CLASS_engine(plc_pack->engine, total_max_l, pbh_info);

  context = plc_pack;
  //*/

  // Initialise m_min, m_max and the rest
  initialise_params();

  // Initialise CLASS before runing MultiNest
  // plc_pack->read_pbh_files(pbh_file_root);
  // plc_pack->initialise_CLASS();

  // Calling MultiNest

  pc_multinest(multinest_loglike, multinest_dumper, settings, context);

  // Deallocate memory
  delete[] settings.pWrap;

  for (std::vector<clik_struct*>::iterator clik_struct_it = plc_pack->clik_objs.begin(); clik_struct_it != plc_pack->clik_objs.end(); ++clik_struct_it) {
    free((*clik_struct_it)->clik_id);
    delete *clik_struct_it;
  }

  delete plc_pack->engine;
  delete plc_pack;

  return 0;
}

void pc_multinest(void (*loglike)(double*,
                                  int&,
                                  int&,
                                  double&,
                                  void*),
                  void (*dumper)(int&,
                                 int&,
                                 int&,
                                 double**,
                                 double**,
                                 double**,
                                 double&,
                                 double&,
                                 double&,
                                 double&,
                                 void*),
                  multinest_settings& s,
                  void*& context) {

  nested::run(s.IS,
              s.mmodal,
              s.ceff,
              s.nlive,
              s.tol,
              s.efr,
              s.ndims,
              s.nPar,
              s.nClsPar,
              s.maxModes,
              s.updInt,
              s.Ztol,
              s.root,
              s.seed,
              s.pWrap,
              s.fb,
              s.resume,
              s.outfile,
              s.initMPI,
              s.logZero,
              s.maxiter,
              loglike,
              dumper,
              context);
}
