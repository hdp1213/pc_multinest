#include "init_plc.hpp"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

// Global array declarations
double m_min[FREE_PARAM_AMT], m_max[FREE_PARAM_AMT];
double m_value[FIXED_PARAM_AMT];
bool m_is_log10[FREE_PARAM_AMT];
trans_t m_transform[FREE_PARAM_AMT];

char* m_name[FREE_PARAM_AMT + DERIVED_PARAM_AMT];
char* m_latex[FREE_PARAM_AMT + DERIVED_PARAM_AMT];

// Function prototypes
void initialise_arrays();
void write_paramnames(const std::string& root);
void write_ranges(const std::string& root);

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
  // Set defaults for global arrays
  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    m_is_log10[param] = false;
  }

  #include "TTTEEE+lowP_pbh_clark-flat.cpp"
  #include "TTTEEE+lowP_pbh_clark-names.cpp"

  // Set values of m_transform depending on includes
  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    m_transform[param] = m_is_log10[param] ? &pow10 : &self;
  }
}

void
write_paramnames(const std::string& root) {
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
write_ranges(const std::string& root) {
  FILE * fout;
  std::string range_file = root + ".ranges";

  fout = fopen(range_file.c_str(), "w");

  for (int param = 0; param < FREE_PARAM_AMT; ++param) {
    trans_t func = m_transform[param];
    fprintf(fout, "%-22s%#17.7E%#17.7E\n", m_name[param], func(m_min[param]), func(m_max[param]));
  }

  fclose(fout);
}
