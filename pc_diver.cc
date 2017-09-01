#include "pc_diver.h"

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <limits>

#include "diver.hpp"
#include "diver_loglike.h"

int main(int argc, char** argv)
{
  ////// Diver settings //////

  diver_settings settings;

  settings.nPar               = FREE_PARAM_AMT;
  settings.lowerbounds        = m_min;
  settings.upperbounds        = m_max;
  settings.nDerived           = DERIVED_PARAM_AMT;
  // settings.path               = "output/pc_diver-";
  settings.nDiscrete          = 0;
  settings.discrete           = {none};
  settings.partitionDiscrete  = false;

  settings.maxciv             = 100;
  settings.maxgen             = 100;
  settings.NP = 10*settings.nPar;
  settings.nF                 = 1;
  // settings.F                  = {0.6};
  settings.Cr                 = 0.9;
  settings.lambda             = 0.8;
  settings.current            = false;
  settings.expon              = false;
  settings.bndry              = 3;
  settings.jDE                = true;
  settings.lambdajDE          = true;

  settings.convthresh         = 1.e-3;
  settings.convsteps          = 10;
  settings.removeDuplicates   = true;

  settings.doBayesian         = true;
  settings.maxNodePop         = 1.9;
  settings.Ztolerance         = 0.1;

  settings.savecount          = 1;
  settings.resume             = false; // watch out for me please!
  settings.outputSamples      = true;
  settings.init_pop_strategy  = 2; // fatal n-shot
  settings.max_init_attempts  = 10;
  settings.max_acceptable_val = 1e6;
  settings.verbose            = 2;


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

  if (argc == 2) {
    settings.path = argv[1];
  }
  else {
    settings.path = "output/pc_diver-";
  }

  std::cout << "Printing results to " << settings.path << std::endl;

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
  initialise_param_arrays();

  // Initialise CLASS before runing MultiNest
  // plc_pack->read_pbh_files(pbh_file_root);
  // plc_pack->initialise_CLASS();

  // Calling Diver

  pc_diver(diver_loglike, diver_prior, settings, context);

  // Deallocate memory
  for (std::vector<clik_struct*>::iterator clik_struct_it = plc_pack->clik_objs.begin(); clik_struct_it != plc_pack->clik_objs.end(); ++clik_struct_it) {
    free((*clik_struct_it)->clik_id);
    delete *clik_struct_it;
  }

  delete plc_pack->engine;
  delete plc_pack;

  return 0;
  
}

void pc_diver(double (*obj_func)(double[],
                                 const int,
                                 int&,
                                 bool&,
                                 const bool,
                                 void*&),
              double (*prior_func)(const double[],
                                   const int,
                                   void*&),
              diver_settings& s,
              void*& context) {

  const int maxciv = s.maxciv;
  const int maxgen = s.maxgen;

  cdiver(obj_func,
         s.nPar,
         const_cast<double*>(s.lowerbounds),
         const_cast<double*>(s.upperbounds),
         s.path.c_str(),
         s.nDerived,
         s.nDiscrete,
         const_cast<int*>(s.discrete),
         s.partitionDiscrete,
         maxciv,
         maxgen,
         s.NP,
         s.nF,
         const_cast<double*>(s.F),
         s.Cr,
         s.lambda,
         s.current,
         s.expon,
         s.bndry,
         s.jDE,
         s.lambdajDE,
         s.convthresh,
         s.convsteps,
         s.removeDuplicates,
         s.doBayesian,
         prior_func,
         s.maxNodePop,
         s.Ztolerance,
         s.savecount,
         s.resume,
         s.outputSamples,
         s.init_pop_strategy,
         s.max_init_attempts,
         s.max_acceptable_val,
         context,
         s.verbose);

}
