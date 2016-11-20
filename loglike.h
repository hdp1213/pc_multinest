#ifndef LOGLIKE_H
#define LOGLIKE_H

/***********************/
/*      MultiNest      */
/***********************/

// Input arguments
// ndim             = dimensionality (total number of free parameters) of the problem
// npars            = total number of free plus derived parameters
// context            void pointer, any additional information
//
// Input/Output arguments
// Cube[npars]      = on entry has the ndim parameters as [0,1]
//                    on exit, the physical parameters plus copy any
//                    derived parameters you want to store with the free
//                    parameters
//   
// Output arguments
// lnew             = loglikelihood
void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
  // ndim     = # of CLASS parameters + nuisance parameters
  //          = 6 + 16 = 22
  // npars    = total dimension of Cube
  //          = 22 + # free parameters
  // Cube     = can use nuisance parameter part to construct cl_and_pars
  // context  = PLC object. has to be! In the future for now.
  
  

  // My own variables
  const int CL_AMT = 6;
  double par_min[npars], par_max[npars];

  // Variables for PLC manipulations
  PLCPack* plc_pack;
  int max_l;

  // CLASS variables
  ClassParams class_params;
  ClassEngine* class_engine(0);
  std::vector<unsigned> l_vec;
  std::vector<double> cl_tt, cl_ee, cl_bb, cl_te, nuisance_pars;
  std::vector<std::vector<double> > class_cls;

  // Initialise the nuisance parameter priors using results outlined in
  //  Planck 2015 results. XI. CMB power spectra, ... . p 21,41

  // Input parameters for CLASS
  par_min[ClikPar::omega_b] = 0.005;
  par_max[ClikPar::omega_b] = 0.1;
  par_min[ClikPar::omega_cdm] = 0.001;
  par_max[ClikPar::omega_cdm] = 0.99;
  par_min[ClikPar::hundredxtheta_s] = 0.5;
  par_max[ClikPar::hundredxtheta_s] = 10.0;
  par_min[ClikPar::tau_reio] = 0.01;
  par_max[ClikPar::tau_reio] = 0.8;
  par_min[ClikPar::n_s] = 0.8;
  par_max[ClikPar::n_s] = 1.2;
  par_min[ClikPar::ln10_10_A_s] = 2.0;
  par_max[ClikPar::ln10_10_A_s] = 4.0;

  // Nuisance parameters for PLC
  par_min[ClikPar::A_cib_217] = 60.0;
  par_max[ClikPar::A_cib_217] = 70.0;
  par_min[ClikPar::cib_index] = -1.3;
  par_max[ClikPar::cib_index] = -1.3;
  par_min[ClikPar::xi_sz_cib] = 0.0;
  par_max[ClikPar::xi_sz_cib] = 0.5;
  par_min[ClikPar::A_sz] = 5.0;
  par_max[ClikPar::A_sz] = 10.0;
  par_min[ClikPar::ps_A_100_100] = 200.0;
  par_max[ClikPar::ps_A_100_100] = 300.0;
  par_min[ClikPar::ps_A_143_143] = 0.0;
  par_max[ClikPar::ps_A_143_143] = 100.0;
  par_min[ClikPar::ps_A_143_217] = 0.0;
  par_max[ClikPar::ps_A_143_217] = 100.0;
  par_min[ClikPar::ps_A_217_217] = 50.0;
  par_max[ClikPar::ps_A_217_217] = 150.0;
  par_min[ClikPar::ksz_norm] = 0.0;
  par_max[ClikPar::ksz_norm] = 7.0;
  par_min[ClikPar::gal545_A_100] = 0.0;
  par_max[ClikPar::gal545_A_100] = 20.0;
  par_min[ClikPar::gal545_A_143] = 0.0;
  par_max[ClikPar::gal545_A_143] = 20.0;
  par_min[ClikPar::gal545_A_143_217] = 0.0;
  par_max[ClikPar::gal545_A_143_217] = 30.0;
  par_min[ClikPar::gal545_A_217] = 60.0;
  par_max[ClikPar::gal545_A_217] = 100.0;
  par_min[ClikPar::calib_100T] = 0.0;
  par_max[ClikPar::calib_100T] = 2.0;
  par_min[ClikPar::calib_217T] = 0.0;
  par_max[ClikPar::calib_217T] = 2.0;
  par_min[ClikPar::A_planck] = 0.9;
  par_max[ClikPar::A_planck] = 1.1;

  // Scale parameters
  for (int i = 0; i < npars; i++) {
    if (i < ndim) { // free parameter
      Cube[i] = par_min[i] + (par_max[i] - par_min[i]) * Cube[i];
    }
    else { // derived parameter
      Cube[i] = (par_min[i] + par_max[i]) / 2.0;
    }
  }

  // Set any fixed parameters
  Cube[ClikPar::cib_index] = -1.3;

  // Setting to Planck best fit for base_plikHM_TT_lowTEB
  // Just as a test for Gaussianity of A_cib_217
  // Cube[ClikPar::A_cib_217] = 66.6;
  Cube[ClikPar::xi_sz_cib] = 0.05;
  Cube[ClikPar::A_sz] = 7.14;
  Cube[ClikPar::ps_A_100_100] = 251.8;
  Cube[ClikPar::ps_A_143_143] = 39.2;
  Cube[ClikPar::ps_A_143_217] = 33.6;
  Cube[ClikPar::ps_A_217_217] = 97.8;
  Cube[ClikPar::ksz_norm] = 0.0;
  Cube[ClikPar::gal545_A_100] = 7.41;
  Cube[ClikPar::gal545_A_143] = 8.98;
  Cube[ClikPar::gal545_A_143_217] = 17.53;
  Cube[ClikPar::gal545_A_217] = 82.0;
  Cube[ClikPar::calib_100T] = 0.99789;
  Cube[ClikPar::calib_217T] = 0.99593;
  Cube[ClikPar::A_planck] = 1.00030;
  

  /* plc_class code begins */


  /******************************/
  /*   Planck Likelihood Code   */
  /******************************/

  plc_pack = static_cast<PLCPack*>(context);

  // Get maximum l over all included .clik files
  max_l = plc_pack->get_largest_max_l();


  /******************************/
  /*           CLASS            */
  /******************************/

  // MultiNest parameters to sweep over (be careful here!)
  // class_params.add("omega_b", Cube[ClikPar::omega_b]);
  // class_params.add("omega_cdm", Cube[ClikPar::omega_cdm]);
  // class_params.add("100*theta_s", Cube[ClikPar::hundredxtheta_s]);
  // class_params.add("tau_reio", Cube[ClikPar::tau_reio]);
  // class_params.add("n_s", Cube[ClikPar::n_s]); // k_0 = 0.05 Mpc^-1 by default
  // class_params.add("ln10^{10}A_s", Cube[ClikPar::ln10_10_A_s]);

  // Set variables to Planck best fit
  class_params.add("omega_b", 0.022242);
  class_params.add("omega_cdm", 0.11977);
  class_params.add("100*theta_s", 1.040862);
  class_params.add("tau_reio", 0.0781);
  class_params.add("n_s", 0.9658); // k_0 = 0.05 Mpc^-1 by default
  class_params.add("ln10^{10}A_s", 3.0904);

  // Options to set for spectra output
  class_params.add("output", "tCl,pCl"); // pCl, lCl for lensed spectra
  class_params.add("l_max_scalars", max_l);
  class_params.add("format", "camb");

  try {
    class_engine = new ClassEngine(class_params);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS failed, throwing exception " << e.what() << std::endl;
    throw e;
  }

  // Create vector of multipoles to calculate Cls for
  l_vec = plc_pack->get_class_l_vec();

  // Calculate Cls from CLASS
  try {
    class_engine->getCls(l_vec, cl_tt, cl_te, cl_ee, cl_bb);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] Spectra extraction unsuccessful, threw " << e.what() << std::endl;
    throw e;
  }


  /******************************/
  /*   Likelihood Calculation   */
  /******************************/

  // Push all spectra into a single matrix
  class_cls.push_back(cl_tt);
  class_cls.push_back(cl_ee);
  class_cls.push_back(cl_bb);
  class_cls.push_back(cl_te);

  plc_pack->create_all_cl_and_pars(Cube, class_cls);

  // Compute the log likelihood using PLC
  lnew = plc_pack->calculate_likelihood();

  //*
  std::cout << "[plc_class] Calculated log likelihood of "
            << lnew
            << std::endl;
  //*/

  // Clean up after running
  delete class_engine;
}

#endif
