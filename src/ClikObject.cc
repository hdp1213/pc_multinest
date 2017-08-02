#include "ClikObject.h"

#include <cstdio> // for stderr
#include <exception> // for std::exception
#include <iostream> // for std::cerr

ClikObject::ClikObject(char* clik_path) : m_clik_path(clik_path) {
  _m_err = initError();
  m_err = &_m_err;

  m_clik_id = clik_init(m_clik_path, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  // Get flags indicating what spectra are contained in .clik file:
  //   TT EE BB TE TB EB
  clik_get_has_cl(m_clik_id, m_cl_flags, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  // Get array of maximum multipole for each existing spectra
  clik_get_lmax(m_clik_id, m_cl_max_ls, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  // Do checks on clik_object for compatibility with CLASS
  if (incompatible_with_class()) {
    std::cerr << "[ERROR]: incompatible_with_class() returned error code "
              << incompatible_with_class()
              << std::endl;
    throw std::exception();
  }

  // Get maximum l from .clik file
  // Already know each spectra has the same length
  m_max_l = m_cl_max_ls[0];

  // Get number of nuisance parameters
  m_nuis_param_amt = clik_get_extra_parameter_names(m_clik_id, &m_nuis_param_names, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  // Determine size of cls_and_pars (cap) array
  m_cap_size = 0;
  for (int i = 0; i < CL_AMT; ++i) {
    if (m_cl_flags[i] == 1) {
      m_cap_size += m_cl_max_ls[i] + 1; // +1 for l=0 case
    }
  }
  m_cap_size += m_nuis_param_amt;

  m_cl_and_pars = new double[m_cap_size];
}

ClikObject::~ClikObject() {
  clik_cleanup(&m_clik_id);
  delete[] m_nuis_param_names;
  delete[] m_cl_and_pars;
}

int ClikObject::get_max_l() const {
  return m_max_l;
}

void ClikObject::set_nuisance_param_enums(std::vector<ClikPar::param_t>& nuisance_params) {
  // Check if we have done enough initialising
  if (nuisance_params.size() > m_nuis_param_amt) {
    std::cerr << "[ERROR] More nuisance parameters given to "
              << "ClikObject::set_nuisance_param_enums() than needed!"
              << std::endl;
    throw std::exception();
  }
  else if (nuisance_params.size() < m_nuis_param_amt) {
    std::cerr << "[ERROR] Not enough nuisance parameters have been "
              << "given to ClikObject::set_nuisance_param_enums()!"
              << std::endl;
    throw std::exception();
  }

  // Allocate enums into free or fixed vectors
  for (std::vector<ClikPar::param_t>::iterator param_it = nuisance_params.begin(); param_it != nuisance_params.end(); ++param_it) {
    if (*param_it < ClikPar::UP_TO_FREE_PARAMS) {
      m_nuis_pars.push_back(static_cast<int>(*param_it));
      m_nuis_par_is_free.push_back(true);
    }

    else if (*param_it < ClikPar::UP_TO_FIXED_PARAMS) {
      m_nuis_pars.push_back(static_cast<int>(*param_it - ClikPar::UP_TO_FIXED_PARAMS));
      m_nuis_par_is_free.push_back(false);
    }

    else {
      std::cerr << "[ERROR] Cannot give a derived nuisance parameter "
                << "to ClikObject::set_nuisance_param_enums()!"
                << std::endl;
      throw std::exception();
    }
  }
}

void ClikObject::create_cl_and_pars(double* free_params,
      double* fixed_params,
      std::vector<std::vector<double> >& class_cls) {
  // Add spectra from CLASS
  int cap_ind = 0;
  for (int cl_ind = 0; cl_ind < CL_AMT; ++cl_ind) {
    if (m_cl_flags[cl_ind] == 1) { // if spectra exists
      for (int l = 0; l <= m_max_l; ++l) {
        // Correct for PLC wanting l=0,1 multipoles
        // These cannot be computed by CLASS and so are set to zero
        if (l < 2) {
          m_cl_and_pars[cap_ind] = 0.0;
        }
        else {
          // CLASS Cls start from l=2, hence the l-2=0 first index
          m_cl_and_pars[cap_ind] = class_cls[cl_ind][l-2];
        }

        cap_ind++;
      }
    }
  }

  // Then add nuisance parameters at the end
  // Nuisance parameters need to be stored in the same order they appear
  //  in the set_nuisance_param_enums() input vector
  for (unsigned i = 0; i < m_nuis_pars.size(); ++i) {
    int param = m_nuis_pars[i];

    if (m_nuis_par_is_free[i]) {
      m_cl_and_pars[cap_ind] = free_params[param];
    }

    else {
      m_cl_and_pars[cap_ind] = fixed_params[param];
    }

    cap_ind++;
  }
}

double ClikObject::get_likelihood() const {
  double loglike;

  loglike = clik_compute(m_clik_id, m_cl_and_pars, m_err);
  quitOnError(*m_err, __LINE__, stderr);

  return loglike;
}

std::ostream& operator<<(std::ostream& os, const ClikObject& clik_object) {
  os << "This .clik file requires "
     << clik_object.m_nuis_param_amt
     << " nuisance parameters to be marginalised over:\n";

  for (int i = 0; i < clik_object.m_nuis_param_amt; ++i) {
    os << '\t'
       << clik_object.m_nuis_param_names[i]
       << std::endl;
  }
  os << std::endl;

  return os;
}

// Returns various error codes if the validation is unsuccessful
//  0  .clik file is valid
//  1  .clik file has uneven spectra lengths
//  2  .clik file contains spectra which cannot be generated by CLASS
int ClikObject::incompatible_with_class() const {
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
