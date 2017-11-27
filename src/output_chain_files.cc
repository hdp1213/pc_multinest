#include "TTTEEE+lowP_pbh.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

// Global array declarations
double m_min[FREE_PARAM_AMT], m_max[FREE_PARAM_AMT];
double m_value[FIXED_PARAM_AMT];

char* m_name[FREE_PARAM_AMT + DERIVED_PARAM_AMT];
char* m_latex[FREE_PARAM_AMT + DERIVED_PARAM_AMT];

// Function prototypes
void initialise_arrays();
void write_paramnames(std::string root);
void write_ranges(std::string root);

// Main function
int
main(int argc, char const *argv[]) {
  std::string root;

  if (argc > 1) {
    root = std::string(argv[1]);
  }
  else {
    std::cerr << "[ERROR]: Must specify a root" << std::endl;
    return 1;
  }

  initialise_arrays();

  write_paramnames(root);
  write_ranges(root);

  return 0;
}


// Function implementations
void
initialise_arrays() {
  #include "TTTEEE+lowP_pbh-flat.cc"
  #include "TTTEEE+lowP_pbh-names.cc"
}

void
write_paramnames(std::string root) {
  std::ofstream fout;
  std::string paramnames_file = root + ".paramnames";

  fout.open(paramnames_file.c_str(), std::ofstream::out | std::ofstream::trunc);

  for (int param = 0; param < FREE_PARAM_AMT + DERIVED_PARAM_AMT; ++param) {
    fout << m_name[param]
         << '\t'
         << m_latex[param]
         << std::endl;
  }

  fout.close();
}

void
write_ranges(std::string root) {
  FILE * fout;
  std::string range_file = root + ".ranges";

  fout = fopen(range_file.c_str(), "w");

  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    fprintf(fout, "%-22s%#17.7E%#17.7E\n", m_name[param], m_min[param], m_max[param]);
  }

  fclose(fout);
}
