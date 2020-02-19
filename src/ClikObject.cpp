#include "ClikObject.hpp"

#include <cstdio>  // stderr
#include <iostream>
#include <stdexcept>

using namespace std;

ClikObject::ClikObject(const string& clik_file,
                       vector<param_t>& nuis_params)
                       : m_clik_id(NULL), m_nuis_params(0), m_lmax(0), m_cap_size(0) {
  int cl_flags[CL_AMT];
  int cl_max_ls[CL_AMT];
  int nuis_param_amt;
  parname* nuis_param_names;

  char* clik_path_c = const_cast<char*>(clik_file.c_str());
  error* err = initError();

  // Open clik file
  m_clik_id = clik_init(clik_path_c, &err);
  quitOnError(err, __LINE__, stderr);

  // Get flags indicating what spectra are contained in the .clik file
  // In order:  TT EE BB TE TB EB
  clik_get_has_cl(m_clik_id, cl_flags, &err);
  quitOnError(err, __LINE__, stderr);

  // Get lmax for each spectra
  clik_get_lmax(m_clik_id, cl_max_ls, &err);
  quitOnError(err, __LINE__, stderr);

  // Get nuisance parameter names and amount
  nuis_param_amt = clik_get_extra_parameter_names(m_clik_id,
                                                  &nuis_param_names,
                                                  &err);
  quitOnError(err, __LINE__, stderr);

  // Check things
  if (nuis_params.size() > nuis_param_amt) {
    throw range_error("Too many nuisance parameters given");
  }
  else if (nuis_params.size() < nuis_param_amt) {
    throw range_error("Not enough nuisance parameters given");
  }

  m_nuis_params.reserve(nuis_param_amt);
  m_nuis_params = nuis_params;

  // Check consistency of lmax across present spectra
  // Also calculate size of cl_and_pars array
  m_lmax = cl_max_ls[0];
  for (unsigned cl = 0; cl < CL_AMT; ++cl) {
    m_has_cl[cl] = (cl_flags[cl] == 1);

    if (not m_has_cl[cl]) continue;

    // Add extra one for l=0 value
    m_cap_size += 1 + cl_max_ls[cl];

    if (cl_max_ls[cl] != m_lmax) {
      throw range_error("lmax not consistent across spectra");
    }
  }
  m_cap_size += nuis_param_amt;

#ifdef DBUG
  cout << "Nuisance parameters:" << endl;
  for (int i = 0; i < nuis_param_amt; ++i) {
    cout << "  " << m_nuis_params[i] << ": " << nuis_param_names[i] << endl;
  }
#endif

  free(err);
  free(nuis_param_names);
}

ClikObject::~ClikObject() {
  free(m_clik_id);
}

double
ClikObject::calculate_likelihood(double input_vals[],
                                 vector<double> class_cls[]) const {
  double cl_and_pars[m_cap_size];
  unsigned cap_ind = 0;

  // First add Cl values
  for (unsigned cl_ind = 0; cl_ind < CL_AMT; ++cl_ind) {
    if (not m_has_cl[cl_ind]) continue;

    for (unsigned l = 0; l <= m_lmax; ++l) {
      // PLC wants multipole values at l=0,1
      // These cannot be computed by CLASS and so are set to zero
      if (l < CLASS_MIN_L) {
        cl_and_pars[cap_ind] = 0.0;
      }
      else {
        cl_and_pars[cap_ind] = class_cls[cl_ind][l - CLASS_MIN_L];
      }

      cap_ind++;
    }
  }

  // Then add nuisance parameters
  // Note that fixed parameters already have their value stored
  for (auto param : m_nuis_params) {
    if (param < UP_TO_FREE_PARAMS) {
      cl_and_pars[cap_ind] = input_vals[param];
    }
    else {
      cl_and_pars[cap_ind] = g_value[param - UP_TO_FREE_PARAMS];
    }

    cap_ind++;
  }

  // Finally, calculate likelihood
  error* err = initError();
  double loglike = clik_compute(m_clik_id, cl_and_pars, &err);
  quitOnError(err, __LINE__, stderr);

  free(err);

  return loglike;
}
