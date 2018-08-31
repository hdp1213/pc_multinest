#include "hyrec_io.hpp"

#include "io_params.h"

#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>

/* Stolen from HyRec, but updated for C++ */
double** hyrec_create_2D_array(int n1, int n2) {
  double** matrix = new double*[n1];

  for (int i = 0; i < n1; ++i) {
    matrix[i] = new double[n2];
  }

  return matrix;
}

/* Stolen from HyRec, but updated for C++ */
void hyrec_free_2D_array(int n1, double** matrix) {
  for (int i = 0; i < n1; ++i) {
    delete[] matrix[i];
  }

  delete[] matrix;
}

/* Stolen shamelessly from HyRec. Didn't even need to update it, haha */
void hyrec_maketab(double xmin, double xmax, unsigned Nx, double *xtab) {
    unsigned i;
    double h = (xmax - xmin)/((double) Nx - 1.0);

    for (i = 0; i < Nx; i++) xtab[i] = xmin + (double) i * h;
}


/*********** Effective rates *************/

void read_alpha(const std::string& root, external_info* info) {
  std::ostringstream alpha_filename;

  alpha_filename << root << "Alpha_inf_unity.dat";
  std::cout << "reading in: '" << alpha_filename.str() << "'" << std::endl;

  /* what the hell am i doing */
  info->logAlpha_tab[0] = hyrec_create_2D_array(NTM, NTR);
  info->logAlpha_tab[1] = hyrec_create_2D_array(NTM, NTR);

  // Read in Alpha file

  std::ifstream alpha_file(alpha_filename.str().c_str());

  if (alpha_file.is_open()) {
    for (int i = 0; i < NTR; ++i) {
      for (int j = 0; j < NTM; ++j) {
        alpha_file >> info->logAlpha_tab[0][j][i] >> info->logAlpha_tab[1][j][i];

        info->logAlpha_tab[0][j][i] = log(info->logAlpha_tab[0][j][i]);
        info->logAlpha_tab[1][j][i] = log(info->logAlpha_tab[1][j][i]);
      }
    }
  }
  else {
    std::cerr << "[ERROR]: Could not read in Alpha file '"
              << alpha_filename.str() << "'"
              << std::endl;

    throw std::exception();
  }

  alpha_file.close();
}

void read_RR(const std::string& root, external_info* info) {
  std::ostringstream RR_filename;

  RR_filename << root << "R_inf.dat";
  std::cout << "reading in: '" << RR_filename.str() << "'" << std::endl;

  // Read in RR file

  std::ifstream RR_file(RR_filename.str().c_str());

  if (RR_file.is_open()) {
    for (int i = 0; i < NTR; ++i) {
      RR_file >> info->logR2p2s_tab[i];
      info->logR2p2s_tab[i] = log(info->logR2p2s_tab[i]);
    }
  }
  else {
    std::cerr << "[ERROR]: Could not read in RR file '"
              << RR_filename.str() << "'"
              << std::endl;

    throw std::exception();
  }

  RR_file.close();
}


/************ Two-photon rates ************/

void read_two_photon(const std::string& root, external_info* info) {
  std::ostringstream two_phot_filename;

  two_phot_filename << root << "two_photon_tables.dat";
  std::cout << "reading in: '" << two_phot_filename.str() << "'" << std::endl;

  // Read in two photon file

  std::ifstream two_phot_file(two_phot_filename.str().c_str());

  for (int b = 0; b < NVIRT; b++) {
    two_phot_file >> info->Eb_tab[b]
                  >> info->A1s_tab[b]
                  >> info->A2s_tab[b]
                  >> info->A3s3d_tab[b]
                  >> info->A4s4d_tab[b];
  }

  two_phot_file.close();
}

void initialise_temp_tables(external_info* info) {
  /* These actually set the values of the table */
  hyrec_maketab(log(TR_MIN), log(TR_MAX), NTR, info->logTR_tab);
  hyrec_maketab(log(TM_TR_MIN), log(TM_TR_MAX), NTM, info->logTM_TR_tab);
  info->DlogTR = info->logTR_tab[1] - info->logTR_tab[0];
  info->DlogTM_TR = info->logTM_TR_tab[1] - info->logTM_TR_tab[0];
}

void normalise_atomic(external_info* info) {
  double L2s1s_current;

  /* Normalize 2s--1s differential decay rate to L2s1s (can be set by user in hydrogen.h) */
  L2s1s_current = 0.;
  for (int b = 0; b < NSUBLYA; b++) L2s1s_current += info->A2s_tab[b];
  for (int b = 0; b < NSUBLYA; b++) info->A2s_tab[b] *= L2s1s/L2s1s_current;


  /* Switches for the various effects considered in Hirata (2008) and diffusion:
      Effect A: correct 2s-->1s rate, with stimulated decays and absorptions of non-thermal photons
      Effect B: Sub-Lyman-alpha two-photon decays
      Effect C: Super-Lyman-alpha two-photon decays
      Effect D: Raman scattering */

   #if (EFFECT_A == 0)
     for (int b = 0; b < NSUBLYA; b++) info->A2s_tab[b] = 0;
   #endif
   #if (EFFECT_B == 0)
     for (int b = 0; b < NSUBLYA; b++) info->A3s3d_tab[b] = info->A4s4d_tab[b] = 0;
   #endif
   #if (EFFECT_C == 0)
      for (int b = NSUBLYA; b < NVIRT; b++) info->A3s3d_tab[b] = info->A4s4d_tab[b] = 0;
   #endif
   #if (EFFECT_D == 0)
      for (int b = NSUBLYA; b < NVIRT; b++) info->A2s_tab[b] = 0;
      for (int b = NSUBLYB; b < NVIRT; b++) info->A3s3d_tab[b] = 0;
   #endif
   #if (DIFFUSION == 0)
      for (int b = 0; b < NVIRT; b++) info->A1s_tab[b] = 0;
   #endif
}
