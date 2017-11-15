//--------------------------------------------------------------------------
//
// Description:
// 	class ClassEngine :
// encapsulation of class calls
//
//
// Author List:
//	Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//
// History (add to end):
//	creation:   ven. nov. 4 11:02:20 CET 2011
//
//-----------------------------------------------------------------------

#ifndef ClassEngine_hh
#define ClassEngine_hh

//CLASS
#include "class.h"

#include "Engine.hh"
//STD
#include <string>
#include <vector>
#include <utility>
#include <ostream>

using std::string;

//general utility to convert safely numerical types to string
template<typename T> std::string str(const T &x);
//specialisations
template<> std::string str (const float &x);
template<> std::string str (const double &x);
template<> std::string str (const bool &x); //"yes" or "no"
template<> std::string str (const std::string &x);

std::string str(const char* x);
//////////////////////////////////////////////////////////////////////////
//class to encapsulate CLASS parameters from any type (numerical or string)
class ClassParams{
public:

  ClassParams(){};
  ClassParams( const ClassParams& o):pars(o.pars){};

  //use this to add a CLASS variable
  template<typename T> unsigned add(const string& key,const T& val){
  pars.push_back(make_pair(key,str(val)));
  return pars.size();
  }

  //accesors
  inline unsigned size() const {return pars.size();}
  inline string key(const unsigned& i) const {return pars[i].first;}
  inline string value(const unsigned& i) const {return pars[i].second;}


private:
  std::vector<std::pair<string,string> > pars;
};

///////////////////////////////////////////////////////////////////////////
class ClassEngine : public Engine
{

  friend class ClassParams;

public:
  //constructors
  ClassEngine(const ClassParams& pars);
  //with a class .pre file
  ClassEngine(const ClassParams& pars,const string& precision_file);
  //with a set of PBH splines
  ClassEngine(const ClassParams& pars, struct external_info* info);
  //from a CLASS .ini file with PBH splines
  ClassEngine(const string& init_file, int l_max, struct external_info* info);

  // destructor
  ~ClassEngine();

  //modfiers: _FAILURE_ returned if CLASS pb:
  bool updateParValues(const std::vector<double>& par);


  //get value at l ( 2<l<lmax): in units = (micro-K)^2
  //don't call if FAILURE returned previously
  //throws std::execption if pb

  double getCl(Engine::cltype t,const long &l);
  void getCls(const std::vector<unsigned>& lVec, //input
	      std::vector<double>& cltt,
	      std::vector<double>& clte,
	      std::vector<double>& clee,
	      std::vector<double>& clbb);


  bool getLensing(const std::vector<unsigned>& lVec, //input
	      std::vector<double>& clphiphi,
	      std::vector<double>& cltphi,
	      std::vector<double>& clephi);

 //for BAO
  inline double z_drag() const {return th.z_d;}
  inline double rs_drag() const {return th.rs_d;}
  double get_Dv(double z);

  double get_Da(double z);
  double get_sigma8(double z);
  double get_f(double z);

  double get_Fz(double z);
  double get_Hz(double z);
  double get_Az(double z);

  inline double getTauReio() const {return th.tau_reio;}

  //may need that
  inline int numCls() const {return sp.ct_size;}
  inline double Tcmb() const {return ba.T_cmb;}

  inline int l_max_scalars() const {return _lmax;}

  // Derived parameters
  // Get H0 in km s^-1 Mpc^-1
  inline double get_H0() const {return ba.H0 * _c_ / 1000.0;}
  // Get Omega_b
  inline double get_Omega_b() const {return ba.Omega0_b;}
  // Get Omega_cdm
  inline double get_Omega_cdm() const {return ba.Omega0_cdm;}
  // Get Omega_Lambda
  inline double get_Omega_L() const {return ba.Omega0_lambda;}
  // Get Omega_g
  inline double get_Omega_g() const {return ba.Omega0_g;}
  // Get Omega_k
  inline double get_Omega_k() const {return ba.Omega0_k;}
  // Get sigma8 perturbation parameter
  inline double get_sigma8() const {return sp.sigma8;}
  // Get age in giga years
  inline double get_age() const {return ba.age;}
  // Get conformal age in Mpc
  inline double get_conf_age() const {return ba.conformal_age;}
  // Get curvature parameter K
  inline double get_K() const {return ba.K;}

  //print content of file_content
  void printFC();

private:
  //structures class en commun
  struct file_content fc;
  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */

  ErrorMsg _errmsg;            /* for error messages */
  double * cl;

  //helpers
  bool dofree;
  int freeStructs();

  //call once /model
  int computeCls();

  int class_main(
		 struct file_content *pfc,
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
  //parnames
  std::vector<std::string> parNames;

  struct external_info* m_info;

protected:


};


;
#endif
