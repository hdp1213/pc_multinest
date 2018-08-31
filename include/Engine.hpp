//--------------------------------------------------------------------------
//
// Description:
//  class Engine :
// base class for Boltzmann code
//
//
// Author List:
//  Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//  Harry Poulter (harry.poulter@adelaide.edu.au)
//
// History (add to end):
//  creation:      Tue Mar 13 15:28:50 CET 2012
//  pc_multinest:  Fri Sep 29 12:22:00 AEST 2017
//
//------------------------------------------------------------------------

#ifndef Engine_hh
#define Engine_hh

#include <ostream>
#include <vector>

class Engine {
 public:
  // P stands for lensing potential phi
  enum cltype {TT = 0, EE, TE, BB, PP, TP, EP};

  // Constructors
  Engine();
  Engine(int lmax);

  // Pure virtual methods
  virtual void update_parameters(const std::vector<double>& cosmopars) = 0;

  // Throws std::exception on failure
  // units = (micro-K)^2
  virtual void get_Cls(const std::vector<unsigned>& lVec, // input
                       std::vector<double>& cltt,
                       std::vector<double>& clte,
                       std::vector<double>& clee,
                       std::vector<double>& clbb) = 0;


  virtual bool get_lensing_Cls(const std::vector<unsigned>& lVec, // input
                               std::vector<double>& clpp,
                               std::vector<double>& cltp,
                               std::vector<double>& clep) = 0;


  virtual double z_drag() const = 0;
  virtual double rs_drag() const = 0;

  virtual double get_Dv(double z) = 0;

  virtual double get_Da(double z) = 0;
  virtual double get_sigma8(double z) = 0;
  virtual double get_f(double z) = 0;
  virtual double get_Fz(double z) = 0;
  virtual double get_Az(double z) = 0;
  virtual double get_Hz(double z) = 0;

  virtual double get_tau_reio() const = 0;

  // destructor
  virtual ~Engine() {}

  // write Cl model+lensing in ostream
  virtual void write_Cls(std::ostream &o);
  inline int lmax() { return m_lmax; }

 protected:
  int m_lmax;
};

#endif
