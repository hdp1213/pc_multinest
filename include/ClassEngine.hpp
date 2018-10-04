//--------------------------------------------------------------------------
//
// Description:
//  class ClassEngine :
// encapsulation of class calls
//
//
// Author List:
//  Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//  Harry Poulter (harry.poulter@adelaide.edu.au)
//
// History (add to end):
//  creation:      Fri Nov  4 11:02:20 CET 2011
//  pc_multinest:  Fri Sep 29 12:23:00 AEST 2017
//
//-----------------------------------------------------------------------

#ifndef ClassEngine_hh
#define ClassEngine_hh

#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "class.h"
#include "Engine.hpp"

// General utility to safely convert numerical types to string
template<typename T> std::string str(const T &x);
// Specialisations
template<> std::string str(const float &x);
template<> std::string str(const double &x);
template<> std::string str(const bool &x);  // "yes" or "no"
template<> std::string str(const std::string &x);

std::string str(const char* x);

////////////////////////////////////////////////////////////////////////////
// Class to encapsulate CLASS parameters from any type (numerical or string)
class ClassParams {
 public:
  ClassParams() {}
  ClassParams(const ClassParams& o):m_pars(o.m_pars) {}

  // Use this to add a CLASS variable
  template<typename T> unsigned add(const std::string& key, const T& val) {
    m_pars.push_back(std::make_pair(key, str(val)));
    return m_pars.size();
  }

  // Accessors
  unsigned size() const { return m_pars.size(); }
  std::string key(const unsigned& i) const { return m_pars[i].first; }
  std::string value(const unsigned& i) const { return m_pars[i].second; }


 private:
  std::vector<std::pair<std::string, std::string> > m_pars;
};
///////////////////////////////////////////////////////////////////////////

class ClassEngine : public Engine {
  friend class ClassParams;

 public:
  // Default constructor
  explicit ClassEngine(const ClassParams& pars);
  // With a CLASS .pre file
  ClassEngine(const ClassParams& pars, const std::string& precision_file);
  // From parameters with PBH splines
  ClassEngine(const ClassParams& pars, struct external_info* info);

  // Destructor
  ~ClassEngine();

  // Modfiers: false returned if CLASS fails
  void update_parameters(const std::vector<double>& par);

  // Get value at l ( 2<l<lmax): in units = (micro-K)^2
  // Don't call if FAILURE returned previously
  // Throws std::execption if fails
  void get_Cls(const std::vector<unsigned>& lVec,  // input
               std::vector<double>& cltt,
               std::vector<double>& clte,
               std::vector<double>& clee,
               std::vector<double>& clbb);

  bool get_lensing_Cls(const std::vector<unsigned>& lVec, // input
                       std::vector<double>& clphiphi,
                       std::vector<double>& cltphi,
                       std::vector<double>& clephi);

  double get_Cl_value_at(const long &l, Engine::cltype t);

  // For BAO
  double z_drag() const { return m_th.z_d; }
  double rs_drag() const { return m_th.rs_d; }

  double get_Dv(double z);

  double get_Da(double z);
  double get_sigma8(double z);
  double get_f(double z);

  double get_Fz(double z);
  double get_Hz(double z);
  double get_Az(double z);

  double get_tau_reio() const { return m_th.tau_reio; }

  // May need these
  int num_Cls() const { return m_sp.ct_size; }
  double get_T_cmb() const { return m_ba.T_cmb; }

  int l_max_scalars() const { return m_lmax; }

  // Derived parameters
  // Get H0 in km s^-1 Mpc^-1
  double get_H0() const { return m_ba.H0 * _c_ / 1000.0; }
  // Get Omega_b
  double get_Omega_b() const { return m_ba.Omega0_b; }
  // Get Omega_cdm
  double get_Omega_cdm() const { return m_ba.Omega0_cdm; }
  // Get Omega_Lambda
  double get_Omega_L() const { return m_ba.Omega0_lambda; }
  // Get Omega_g
  double get_Omega_g() const { return m_ba.Omega0_g; }
  // Get Omega_k
  double get_Omega_k() const { return m_ba.Omega0_k; }
  // Get sigma8 perturbation parameter
  double get_sigma8() const { return m_sp.sigma8; }
  // Get age in giga years
  double get_age() const { return m_ba.age; }
  // Get conformal age in Mpc
  double get_conf_age() const { return m_ba.conformal_age; }
  // Get curvature parameter K
  double get_K() const { return m_ba.K; }

  // Print content of file_content
  void print_FC();

 private:
  // Common CLASS structures
  struct file_content m_fc;
  struct precision m_pr;        /* for precision parameters */
  struct background m_ba;       /* for cosmological background */
  struct thermo m_th;           /* for thermodynamics */
  struct perturbs m_pt;         /* for source functions */
  struct transfers m_tr;        /* for transfer functions */
  struct primordial m_pm;       /* for primordial spectra */
  struct spectra m_sp;          /* for output spectra */
  struct nonlinear m_nl;        /* for non-linear spectra */
  struct lensing m_le;          /* for lensed spectra */
  struct output m_op;           /* for output files */

  ErrorMsg m_errmsg;            /* for error messages */
  double * m_cl;

  // Helpers
  bool m_do_free;
  int free_structs();

  // Call once per model
  int compute_Cls();

  // Internal methods for intialisations
  void write_pars_to_fc(const ClassParams& pars, struct file_content* fc);

  int class_main(struct file_content *pfc,
                 struct precision * ppr,
                 struct background * pba,
                 struct thermo * pth,
                 struct perturbs * ppt,
                 struct transfers * ptr,
                 struct primordial * ppm,
                 struct spectra * psp,
                 struct nonlinear * pnl,
                 struct lensing * ple,
                 struct output * pop,
                 ErrorMsg errmsg);

  // Parameter names
  std::vector<std::string> m_parnames;
  struct external_info* m_info;
};

#endif
