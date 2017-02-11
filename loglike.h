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
  
  // Variables for PLC manipulations
  PLCPack* plc_pack;

  // CLASS variables
  std::vector<double> lcdm_params;
  std::vector<double> cl_tt, cl_ee, cl_bb, cl_te;
  std::vector<std::vector<double> > class_cls;

  plc_pack = static_cast<PLCPack*>(context);

  // Initialise the nuisance parameter priors using results outlined in
  //  Planck 2015 results. XI. CMB power spectra, ... . p 21,41
  plc_pack->scale_Cube(Cube);

  /******************************/
  /*           CLASS            */
  /******************************/

  // MultiNest parameters to sweep over (be careful here!)
  lcdm_params.push_back(Cube[ClikPar::omega_b]);
  lcdm_params.push_back(Cube[ClikPar::omega_cdm]);
  lcdm_params.push_back(Cube[ClikPar::hundredxtheta_s]);
  lcdm_params.push_back(Cube[ClikPar::tau_reio]);
  lcdm_params.push_back(Cube[ClikPar::ln10_10_A_s]);
  lcdm_params.push_back(Cube[ClikPar::n_s]); // k_0 = 0.05 Mpc^-1 by default

  try {
    plc_pack->run_CLASS(lcdm_params);
  }
  // CLASS has failed, set loglike to zero and return
  catch (std::exception const &e) {
    //* An error occurred
    std::cerr << "[ERROR] CLASS failed, throwing exception "
              << e.what()
              << std::endl;
    //*
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
    //*/

    lnew = -1E90;
    return;
  }

  // Calculate Cls from CLASS
  try {
    plc_pack->get_CLASS_spectra(cl_tt, cl_te, cl_ee, cl_bb);
  }
  // Spectra extraction has failed, set loglike to zero and return
  catch (std::exception const &e) {
    //* An error occurred
    std::cerr << "[ERROR] Spectra extraction unsuccessful, CLASS threw "
              << e.what()
              << std::endl;
    //*
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
    //*/

    lnew = -1E90;
    return;
  }

  // Set any derived parameter values in Cube
  plc_pack->set_derived_params(Cube);


  /******************************/
  /*   Likelihood Calculation   */
  /******************************/

  // Push all spectra into a single matrix
  class_cls.push_back(cl_tt);
  class_cls.push_back(cl_ee);
  class_cls.push_back(cl_bb);
  class_cls.push_back(cl_te);

  // Prepare PLC clik object(s) for likelihood calculation
  plc_pack->create_all_cl_and_pars(Cube, class_cls);

  // Compute the log likelihood using PLC
  lnew = plc_pack->calculate_likelihood();

  // Subtract any other priors for those parameters with them
  lnew -= plc_pack->calculate_extra_priors(Cube);

  /*
  std::cout << "[plc_class] Calculated log likelihood of "
            << lnew
            << std::endl;
  //*/
}

#endif
