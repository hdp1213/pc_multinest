#include "external.hpp"

#include "pbh_io.hpp"
#include "hyrec_io.hpp"

#include <iostream>

using namespace std;

external_info* initialise_external_info(const string& pbh_root, const string& hyrec_root) {
  cout << "[INFO] Initialising external_info struct"
       << endl;

  /* First read in the PBH stuff */

  external_info* info = new external_info();

  info->hion = new bspline_2d();
  info->excite = new bspline_2d();
  info->heat = new bspline_2d();

  // Read in axes
  read_axes(pbh_root, info);

  // Read in hydrogen ionisation
  read_bicubic_bspline(pbh_root, "hion.dat", info->hion);

  // Read in hydrogen excitation
  read_bicubic_bspline(pbh_root, "excite.dat", info->excite);

  // Read in plasma heating
  read_bicubic_bspline(pbh_root, "heat.dat", info->heat);

  /* Then read in HyRec files */

  // This is good to do, too
  initialise_temp_tables(info);

  read_alpha(hyrec_root, info);

  read_RR(hyrec_root, info);

  read_two_photon(hyrec_root, info);

  // Normalise the stuff
  normalise_atomic(info);

  return info;
}
