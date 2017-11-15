#ifndef HYREC_IO
#define HYREC_IO

#include <string>

#include "class.h"

void hyrec_free_2D_array(int n1, double** matrix);

void read_alpha(std::string& root, external_info* info);
void read_RR(std::string& root, external_info* info);
void read_two_photon(std::string& root, external_info* info);

void initialise_temp_tables(external_info* info);

void normalise_atomic(external_info* info);

#endif
