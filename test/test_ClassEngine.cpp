#include "ClassEngine.hpp"

#include <iostream>
#include <numeric>
#include <stdexcept>

using namespace std;

const unsigned MIN_CL = 2;

int main(int argc, char const *argv[])
{
  ClassEngine* engine;

  unsigned lmax = 2000;
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

  params.add("output", "tCl,pCl,lCl"); // mPk
  params.add("lensing", true);
  params.add("l_max_scalars", lmax);
  params.add("format", "camb");

  try {
    engine = new ClassEngine(params);
  }
  catch (exception const& e) {
    cerr << "[ERROR] " << e.what() << endl;
    return -1;
  }

  vector<unsigned> cl_ls(lmax - MIN_CL + 1);
  iota(cl_ls.begin(), cl_ls.end(), MIN_CL);

  vector<double> cl_tt, cl_te, cl_ee, cl_bb;

  try {
    engine->get_Cls(cl_ls, cl_tt, cl_te, cl_ee, cl_bb);
  }
  catch (exception const& e) {
    cerr << "[ERROR] " << e.what() << endl;
    return -1;
  }

  return 0;
}
