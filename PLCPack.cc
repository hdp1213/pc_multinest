#include "PLCPack.h"

#include <cstdio> // for stderr
#include <exception> // for std::exception
#include <iostream> // for std::cerr

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

int PLCPack::get_largest_max_l() const {
  return m_largest_max_l;
}

void PLCPack::set_pars(ClikPar* clik_par) {
  m_clik_par = clik_par;
}

ClikPar* PLCPack::get_pars() const {
  return m_clik_par;
}

void PLCPack::create_all_cl_and_pars(double* Cube,
      std::vector<std::vector<double> >& class_cls) {
  for (int i = 0; i < m_clik_objects.size(); ++i) {
    m_clik_objects[i]->create_cl_and_pars(Cube, class_cls);
  }
}

double PLCPack::calculate_likelihood() const {
  double loglike = 0.0;

  for (int i = 0; i < m_clik_objects.size(); ++i) {
    loglike += m_clik_objects[i]->get_likelihood();
  }

  return loglike;
}

std::vector<unsigned> PLCPack::get_class_l_vec() const {
  return m_class_l_vec;
}