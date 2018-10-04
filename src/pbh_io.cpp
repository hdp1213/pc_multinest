#include "pbh_io.hpp"

#include <exception>
#include <iostream>
#include <sstream>

using namespace std;

void read_axes(const string& root, external_info* info) {
  ostringstream axes_filename;

  axes_filename << root << "axes.dat";
  cout << "[INFO] Reading '" << axes_filename.str() << "'" << endl;

  ifstream axes_file(axes_filename.str().c_str());

  if (axes_file.is_open()) {
    // Read in redshift axis
    read_1d_array(axes_file, &(info->z_deps), &(info->z_deps_size));

    // Read in mass axis
    read_1d_array(axes_file, &(info->masses), &(info->masses_size));
  }
  else {
    cerr << "[ERROR] Could not read axes file '"
         << axes_filename.str() << "'"
         << endl;

    throw exception();
  }

  axes_file.close();
}

void read_bicubic_bspline(const string& root, const char* channel, bspline_2d* spline) {
  int deg = 3;
  int n_coeffs;

  ostringstream spline_filename;

  spline_filename << root << channel;
  cout << "[INFO] Reading '" << spline_filename.str() << "'" << endl;

  spline->degree = deg;

  // Begin reading in file. Apparently no error is thrown if the file doesn't exist
  ifstream spline_file(spline_filename.str().c_str());

  if (spline_file.is_open()) {
    // Get knots along x-axis
    read_1d_array(spline_file, &(spline->xknots), &(spline->nxknots));

    // Get knots along y-axis
    read_1d_array(spline_file, &(spline->yknots), &(spline->nyknots));

    // Get coefficients
    read_1d_array(spline_file, &(spline->coeffs), &n_coeffs);
  }
  else {
    cerr << "[ERROR] Could not read spline file '"
         << spline_filename.str() << "'"
         << endl;

    throw exception();
  }

  spline_file.close();
}


void read_1d_array(ifstream& file, double** array, int* array_size) {
  string line;

  // Get array size from the first line
  getline(file, line);

  *array_size = atoi(line.c_str());
  *array = new double[*array_size];

  // Get array contents from the next line of the file
  getline(file, line);

  istringstream line_stream(line);
  string cell;

  int i = 0;
  while (getline(line_stream, cell, ',')) {
    istringstream cell_stream(cell);
    cell_stream >> (*array)[i];
    i++;
  }
}
