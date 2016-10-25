#include "PLCPack.h"

#include <cstdio> // for stderr
#include <stdexcept> // for std::exception
#include <iostream> // for std::cerr

// Will this do preliminary checks on the clik_object?
// I feel like it should
PLCPack::PLCPack(clik_object *clik_id) {
  _m_err = initError();
  m_err = &_m_err;

  m_clik_id = clik_id;

  // Get flags indicating what spectra are contained in .clik file:
  //   TT EE BB TE TB EB
  clik_get_has_cl(m_clik_id, m_cl_flags, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  // Get array of maximum multipole for each existing spectra
  clik_get_lmax(clik_id, m_cl_max_ls, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  // Do checks on clik_object for compatibility with CLASS
  if (invalid_clik_object()) {
    std::cerr << "[ERROR]: invalid_clik_object() returned error code "
              << invalid_clik_object()
              << std::endl;
    throw std::exception();
  }

  // Get maximum l from .clik file
  // Already know each spectra has the same length
  m_max_l = m_cl_max_ls[0];

  // Get number of nuisance parameters
  m_param_amt = clik_get_extra_parameter_names(m_clik_id, &m_param_names, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  // Determine size of cls_and_pars (cap) array
  m_cap_size = 0;
  for (int i = 0; i < CL_AMT; ++i) {
    if (m_cl_flags[i] == 1) {
      m_cap_size += m_cl_max_ls[i] + 1; // +1 for l=0 case
    }
  }
  m_cap_size += m_param_amt;
}

PLCPack::~PLCPack() {
  // This might accidentally delete the whole thing after each
  // loop though!
  clik_cleanup(&m_clik_id);
}

int PLCPack::get_max_l() const {
  return m_max_l;
}

int PLCPack::get_cap_size() const {
  return m_cap_size;
}

int* PLCPack::get_cl_flags() {
  return m_cl_flags;
}

int PLCPack::get_param_amt() const {
  return m_param_amt;
}

/*
parname PLCPack::get_param_names() const {
  return m_param_names;
}
*/

double PLCPack::get_clik_likelihood(double* cl_and_pars) {
  double loglike;

  loglike = clik_compute(m_clik_id, cl_and_pars, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  return loglike;
}

// Returns various error codes if the validation is unsuccessful
//  0  .clik file is valid
//  1  .clik file has uneven spectra lengths
//  2  .clik file contains spectra which cannot be generated by CLASS
int PLCPack::invalid_clik_object() const {
  int max_l = -1;
  int prev_max_l = -1;
  bool is_first_checked_spectra = true;

  // Check if Cls are the same length
  for (int i = 0; i < CL_AMT; ++i) {
    if (m_cl_flags[i] == 1) {
      max_l = m_cl_max_ls[i];

      if ( (max_l != prev_max_l) and (not is_first_checked_spectra) ) {
        return 1;
      }

      if (is_first_checked_spectra) {
        is_first_checked_spectra = false;
      }

      prev_max_l = max_l;
    }
  }

  // Check if CLASS can generate the given spectra
  for (int i = 0; i < CL_AMT; ++i) {
    if (m_cl_flags[i] == 1) {
      // i > 3 correspond to spectra which cannot be generated by CLASS
      if (i > 3) {
        return 2;
      }
    }
  }

  return 0;
}