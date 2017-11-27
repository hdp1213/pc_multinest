#ifndef PBH_IO
#define PBH_IO

#include <fstream>
#include <string>

#include "class.h"

// High-level functions used in initialise_external_info()
void read_axes(std::string& root, external_info* info);
void read_bicubic_bspline(std::string& root, const char* channel, bspline_2d* spline);

// Actual function used in reading arrays from file streams
// Allocates memory for the array in the function itself
void read_1d_array(std::ifstream& file, double** array, int* array_size);

#endif
