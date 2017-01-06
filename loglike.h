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
  bool CLASS_not_failed = true;

  plc_pack = static_cast<PLCPack*>(context);

  // Initialise the nuisance parameter priors using results outlined in
  //  Planck 2015 results. XI. CMB power spectra, ... . p 21,41
  plc_pack->get_pars()->scale_Cube(Cube);

  
  /******************************/
  /*   Planck Likelihood Code   */
  /******************************/


  // Get maximum l over all included .clik files
  max_l = plc_pack->get_largest_max_l();


  /******************************/
  /*           CLASS            */
  /******************************/

  // MultiNest parameters to sweep over (be careful here!)
  class_params.add("omega_b", Cube[ClikPar::omega_b]);
  class_params.add("omega_cdm", Cube[ClikPar::omega_cdm]);
  class_params.add("100*theta_s", Cube[ClikPar::hundredxtheta_s]);
  class_params.add("tau_reio", Cube[ClikPar::tau_reio]);
  class_params.add("n_s", Cube[ClikPar::n_s]); // k_0 = 0.05 Mpc^-1 by default
  class_params.add("ln10^{10}A_s", Cube[ClikPar::ln10_10_A_s]);

  // Options to set for spectra output
  class_params.add("output", "tCl,pCl"); // pCl, lCl for lensed spectra
  class_params.add("l_max_scalars", max_l);
  class_params.add("format", "camb");

  try {
    class_engine = new ClassEngine(class_params);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS failed, throwing exception "
              << e.what()
              << std::endl;
    
    std::cout << "[INFO] LCDM parameters used at time were:\n"
              << "\tomega_b         : "
                  << Cube[ClikPar::omega_b] << '\n'
              << "\tomega_cdm       : "
                  << Cube[ClikPar::omega_cdm] << '\n'
              << "\thundredxtheta_s : "
                  << Cube[ClikPar::hundredxtheta_s] << '\n'
              << "\ttau_reio        : "
                  << Cube[ClikPar::tau_reio] << '\n'
              << "\tn_s             : "
                  << Cube[ClikPar::n_s] << '\n'
              << "\tln10_10_A_s     : "
                  << Cube[ClikPar::ln10_10_A_s] << std::endl;

    CLASS_not_failed = false;
  }

  if (CLASS_not_failed) {
    // Create vector of multipoles to calculate Cls for
    l_vec = plc_pack->get_class_l_vec();

    // Calculate Cls from CLASS
    try {
      class_engine->getCls(l_vec, cl_tt, cl_te, cl_ee, cl_bb);
    }
    catch (std::exception const &e) {
      std::cerr << "[ERROR] Spectra extraction unsuccessful, CLASS threw "
                << e.what()
                << std::endl;
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

    // Subtract any other priors for those parameters with them
    lnew -= plc_pack->get_pars()->calculate_extra_priors(Cube, class_engine);
  }

  else { // CLASS has failed, set loglike to zero for this point
    lnew = -1E90;
  }

  /*
  std::cout << "[plc_class] Calculated log likelihood of "
            << lnew
            << std::endl;
  //*/

  // Clean up after running
  delete class_engine;
}

#endif
