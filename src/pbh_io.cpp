#include "pbh_io.hpp"

#include <exception>
#include <iostream>
#include <sstream>

void read_axes(std::string& root, external_info* info) {
  std::ostringstream axes_filename;

  axes_filename << root << "axes.dat";
  std::cout << "reading in: '" << axes_filename.str() << "'" << std::endl;

  std::ifstream axes_file(axes_filename.str().c_str());

  if (axes_file.is_open()) {
    // Read in redshift axis
    read_1d_array(axes_file, &(info->z_deps), &(info->z_deps_size));

    // Read in mass axis
    read_1d_array(axes_file, &(info->masses), &(info->masses_size));
  }
  else {
    std::cerr << "[ERROR]: Could not read in axes file '"
              << axes_filename.str() << "'"
              << std::endl;

    throw std::exception();
  }

  axes_file.close();
}

void read_bicubic_bspline(std::string& root, const char* channel, bspline_2d* spline) {
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
  else {
    std::cerr << "[ERROR]: Could not read in spline file '"
             << spline_filename.str() << "'"
             << std::endl;

    throw std::exception();
  }

  spline_file.close();
}


void read_1d_array(std::ifstream& file, double** array, int* array_size) {
  std::string line;

  // Get array size from the first line
  std::getline(file, line);

  *array_size = atoi(line.c_str());
  *array = new double[*array_size];

  // Get array contents from the next line of the file
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
