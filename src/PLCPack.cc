#include "PLCPack.h"

#include <exception> // for std::exception
#include <iostream> // for std::getline
#include <fstream> // for std::ifstream
#include <sstream> // for std::stringstream

void print_bicubic_bspline(struct bspline_2d* spline) {
  std::cout << spline->nxknots << std::endl;

  for (int i = 0; i < spline->nxknots; ++i) {
    std::cout << spline->xknots[i] << " ";
  }
  std::cout << std::endl;

  std::cout << spline->nyknots << std::endl;
  
  for (int i = 0; i < spline->nyknots; ++i) {
    std::cout << spline->yknots[i] << " ";
  }
  std::cout << std::endl;

  for (int i = 0; i < (spline->nxknots-spline->degree-1)*(spline->nyknots-spline->degree-1); ++i) {
    std::cout << spline->coeffs[i] << " ";
  }
  std::cout << std::endl;
}

PLCPack::PLCPack() : m_largest_max_l(1), m_clik_par(0) {
  
}

PLCPack::~PLCPack() {
  for (int i = 0; i < m_clik_objects.size(); ++i) {
    delete m_clik_objects[i];
  }

  delete m_clik_par;

  // Deallocate b-spline memory
  delete[] m_pbh_info.hion.xknots;
  delete[] m_pbh_info.hion.yknots;
  delete[] m_pbh_info.hion.coeffs;

  delete[] m_pbh_info.excite.xknots;
  delete[] m_pbh_info.excite.yknots;
  delete[] m_pbh_info.excite.coeffs;

  delete[] m_pbh_info.heat.xknots;
  delete[] m_pbh_info.heat.yknots;
  delete[] m_pbh_info.heat.coeffs;

  // Deallocate axes
  delete[] m_pbh_info.z_deps;
  delete[] m_pbh_info.masses;
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

void PLCPack::set_clik_params(ClikPar* clik_par) {
  m_clik_par = clik_par;
}

// Should only be called after all clik objects have been added
void PLCPack::initialise_CLASS() {
  m_clik_par->initialise_CLASS(m_largest_max_l, &m_pbh_info);
}

void PLCPack::read_pbh_files(std::string pbh_root) {
  // Read in axes
  // Variables are set internally
  read_axes(pbh_root);

  // Read in hydrogen ionisation
  read_bicubic_bspline(pbh_root, "hion.dat", &(m_pbh_info.hion));

  // Read in hydrogen excitation
  read_bicubic_bspline(pbh_root, "excite.dat", &(m_pbh_info.excite));

  // Read in plasma heating
  read_bicubic_bspline(pbh_root, "heat.dat", &(m_pbh_info.heat));
}

double PLCPack::calculate_extra_priors(double* Cube) const {
  return m_clik_par->calculate_extra_priors(Cube);
}

void PLCPack::scale_Cube(double* Cube) {
  m_clik_par->scale_Cube(Cube);
}

void PLCPack::set_derived_params(double* Cube) {
  m_clik_par->set_derived_params(Cube);
}

// Throws std::exception on failure
void PLCPack::run_CLASS(std::vector<double> class_params) {
  m_clik_par->get_CLASS()->updateParValues(class_params);
}

// Throws std::exception on failure
void PLCPack::get_CLASS_spectra(std::vector<double>& cl_tt,
      std::vector<double>& cl_te, 
      std::vector<double>& cl_ee, 
      std::vector<double>& cl_bb) {
  m_clik_par->get_CLASS()->getCls(m_class_l_vec, cl_tt, cl_te, cl_ee, cl_bb);
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

#ifdef BAO_LIKE
double PLCPack::calculate_BAO_likelihood() const {
  return m_clik_par->calculate_BAO_likelihood();
}
#endif


void PLCPack::read_axes(std::string root) {
  std::ostringstream axes_filename;

  axes_filename << root << "axes.dat";
  std::cout << "reading in: '" << axes_filename.str() << "'" << std::endl;

  std::ifstream axes_file(axes_filename.str().c_str());

  if (axes_file.is_open()) {
    // Read in redshift axis
    read_1d_array(axes_file, &(m_pbh_info.z_deps), &(m_pbh_info.z_deps_size));

    // Read in mass axis
    read_1d_array(axes_file, &(m_pbh_info.masses), &(m_pbh_info.masses_size));
  }

  axes_file.close();
}

// Read a bicubic b-spline in from an external file
// Passes a pointer to a bspline_2d struct
void PLCPack::read_bicubic_bspline(std::string root, const char* channel, struct bspline_2d* spline) {
  int deg = 3;
  int n_coeffs;

  std::ostringstream spline_filename;

  spline_filename << root << channel;
  std::cout << "reading in: '" << spline_filename.str() << "'" << std::endl;

  spline->degree = deg;

  // Begin reading in file. Apparently no error is thrown if the file doesn't exist
  std::ifstream spline_file(spline_filename.str().c_str());
  
  if (spline_file.is_open()) {
    // Get knots along x-axis
    read_1d_array(spline_file, &(spline->xknots), &(spline->nxknots));

    // Get knots along y-axis
    read_1d_array(spline_file, &(spline->yknots), &(spline->nyknots));

    // Get coefficients
    read_1d_array(spline_file, &(spline->coeffs), &n_coeffs);
  }

  spline_file.close();
}

void PLCPack::read_1d_array(std::ifstream& file, double** array, int* array_size) {
  std::string line;

  // Get array size
  std::getline(file, line);

  *array_size = atoi(line.c_str());
  *array = new double[*array_size];

  // Get array contents
  std::getline(file, line);

  std::istringstream line_stream(line);
  std::string cell;

  int i = 0;
  while (std::getline(line_stream, cell, ',')) {
    std::istringstream cell_stream(cell);
    cell_stream >> (*array)[i];
    i++;
  }
}
