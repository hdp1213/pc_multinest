#include "ClassEngine.hpp"

#include <iostream>
#include <numeric>
#include <stdexcept>

#include "external.hpp"

using namespace std;

const unsigned MIN_CL = 2;
const unsigned CLASS_CL_AMT = 4;

ClassEngine*
initialise_CLASS()
{
  ClassEngine* engine;

  unsigned lmax = 2000;

  // Initialise external info
  string pbh_root = string(CLASS_PBH_FILE_DIR) + "/";
  string hyrec_root = string(HYREC_FILE_DIR) + "/";
  external_info* info = initialise_external_info(pbh_root, hyrec_root);

  ClassParams params;

  params.add("omega_b", 0.022231);
  params.add("omega_cdm", 0.12003);
  params.add("100*theta_s", 1.041740);
  params.add("tau_reio", 0.0807);
  params.add("ln10^{10}A_s", 3.0933);
  params.add("n_s", 0.96749);

  // params.add("N_ur", 3.046);
  params.add("N_ur", 2.0328);
  params.add("N_ncdm", 1);
  params.add("m_ncdm", 0.06); // MeV

  params.add("Omega_pbh_ratio", 1e-90);
  params.add("pbh_mass_mean", 1e7);

  params.add("pbh_mass_dist", "pbh_delta");
  params.add("read external files", false);

  params.add("output", "tCl,pCl,lCl"); // mPk
  params.add("lensing", true);
  params.add("l_max_scalars", lmax);
  params.add("format", "camb");

  params.add("input_verbose", 3);

  try {
    engine = new ClassEngine(params, info);
  }
  catch (exception const& e) {
    cerr << "[ERROR] " << e.what() << endl;
  }

  return engine;
}

vector<vector<double> >
extract_spectra(ClassEngine* engine, vector<unsigned>& ls)
{
  vector<vector<double> > cls = {
    vector<double>(0),
    vector<double>(0),
    vector<double>(0),
    vector<double>(0)
  };

  try {
    engine->get_Cls(ls, cls[0], cls[3], cls[1], cls[2]);
  }
  catch (exception const& e) {
    cerr << "[ERROR] " << e.what() << endl;
  }

  return cls;
}

void
update_CLASS(ClassEngine* engine, vector<double>& params)
{
  try {
    engine->update_parameters(params);
  }
  catch (exception const& e) {
    cerr << "[ERROR] " << e.what() << endl;
  }
}

int
main(int argc, char const *argv[])
{
  cout << "initialising CLASS" << endl;
  ClassEngine* engine = initialise_CLASS();
  unsigned lmax = 2000;

  vector<unsigned> cl_ls(lmax - MIN_CL + 1);
  iota(cl_ls.begin(), cl_ls.end(), MIN_CL);

  vector<double> params = {
    2.1884959457017052e-02,
    1.2495854312444295e-01,
    1.0414785351882943e+00,
    4.9241400412509184e-02,
    3.0760186012521729e+00,
    9.4614236669140450e-01
  };

  cout << "updating CLASS" << endl;
  update_CLASS(engine, params);

  delete engine;

  return 0;
}
