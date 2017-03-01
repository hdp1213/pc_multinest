#include "PLCPack.h"

#include <exception> // for std::exception

PLCPack::PLCPack() : m_largest_max_l(1), m_clik_par(0) {
  
}

PLCPack::~PLCPack() {
  for (int i = 0; i < m_clik_objects.size(); ++i) {
    delete m_clik_objects[i];
  }

  delete m_clik_par;
}


void PLCPack::add_clik_object(ClikObject* clik_object) {
  m_clik_objects.push_back(clik_object);

  // Create a vector of multipoles l which will expand
  // each time a new largest l is found so it will always have
  // the current largest l as its final value
  if (m_largest_max_l < clik_object->get_max_l()) {
    int old_max_l = m_largest_max_l;
    m_largest_max_l = clik_object->get_max_l();

    for (int l = old_max_l + 1; l <= m_largest_max_l; ++l) {
      m_class_l_vec.push_back(l);
    }
  }
}

// Should only be called after all clik objects have been added
void PLCPack::initialise_params(ClikPar* clik_par) {
  m_clik_par = clik_par;
  m_clik_par->initialise_CLASS(m_largest_max_l);
}

double PLCPack::calculate_extra_priors(double* Cube) const {
  m_clik_par->calculate_extra_priors(Cube);
}

void PLCPack::scale_Cube(double* Cube) {
  m_clik_par->scale_Cube(Cube);
}

void PLCPack::set_derived_params(double* Cube) {
  m_clik_par->set_derived_params(Cube);
}

void PLCPack::run_CLASS(std::vector<double> class_params) {
  try {
    m_clik_par->get_CLASS()->updateParValues(class_params);
  }
  catch (std::exception const &e) {
    throw e;
  }
}

void PLCPack::get_CLASS_spectra(std::vector<double>& cl_tt,
      std::vector<double>& cl_te, 
      std::vector<double>& cl_ee, 
      std::vector<double>& cl_bb) {
  try {
    m_clik_par->get_CLASS()->getCls(m_class_l_vec, cl_tt, cl_te, cl_ee, cl_bb);
  }
  catch (std::exception const &e) {
    throw e;
  }
}

void PLCPack::create_all_cl_and_pars(double* Cube,
      std::vector<std::vector<double> >& class_cls) {
  for (int i = 0; i < m_clik_objects.size(); ++i) {
    m_clik_objects[i]->create_cl_and_pars(Cube, class_cls);
  }
}

double PLCPack::calculate_PLC_likelihood() const {
  double loglike = 0.0;

  for (int i = 0; i < m_clik_objects.size(); ++i) {
    loglike += m_clik_objects[i]->get_likelihood();
  }

  return loglike;
}

double PLCPack::calculate_BAO_likelihood() const {
  return m_clik_par->calculate_BAO_likelihood();
}