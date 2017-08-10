#ifndef LOGLIKE_H
#define LOGLIKE_H

#include <iostream>
#include <vector>
#include <stdexcept>

// void pc_loglike(double *input_params, int &n_params)

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
  std::vector<double> cl_tt, cl_ee, cl_bb, cl_te;
  std::vector<std::vector<double> > class_cls;


  plc_pack = static_cast<PLCPack*>(context);

  // Scale the cube using uniform priors defined in ClikPar::ClikPar()
  plc_pack->scale_free_params(Cube);

  /******************************/
  /*           CLASS            */
  /******************************/

  std::cout << "RUNNING CLASS ..." << std::endl;

  try {
    plc_pack->run_CLASS(Cube);
  }
  // CLASS has failed, set loglike to zero and return
  catch (std::exception const &e) {
    //* An error occurred
    std::cerr << "[ERROR] CLASS failed, throwing exception "
              << e.what()
              << std::endl;
    /*
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

  std::cout << "EXTRACTING SPECTRA FROM CLASS ..." << std::endl;

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
    /*
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

  std::cout << "SETTING DERIVED PARAMETERS IN CUBE ..." << std::endl;

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

  std::cout << "CREATING ALL CL AND PARS ..." << std::endl;

  // Prepare PLC clik object(s) for likelihood calculation
  plc_pack->create_all_cl_and_pars(Cube, class_cls);

  std::cout << "CALCULATING LIKELIHOOD ..." << std::endl;

  // Compute the log likelihood using PLC
  lnew = plc_pack->calculate_PLC_likelihood();

#ifdef BAO_LIKE
  // Calculate BAO log likelihood
  lnew -= plc_pack->calculate_BAO_likelihood();
#endif

  // Subtract any other priors for those parameters with them
  lnew -= plc_pack->calculate_extra_likelihoods(Cube);

  //*
  std::cout << "[plc_class] Calculated log likelihood of "
            << lnew
            << std::endl;
  //*/
}

#endif
