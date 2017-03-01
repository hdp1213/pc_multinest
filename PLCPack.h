// Header file for plc_pack class
#ifndef PLCPACK_H
#define PLCPACK_H

#include "ClassEngine.hh"
#include "ClikObject.h"
#include "ClikPar.h"
#include <vector>


class PLCPack {

public:
  PLCPack();
  ~PLCPack();

  // Pack "building" methods
  void add_clik_object(ClikObject* clik_object);
  void initialise_params(ClikPar* clik_par);

  // Higher level ClikPar methods
  double calculate_extra_priors(double* Cube) const;
  void scale_Cube(double* Cube);
  void set_derived_params(double* Cube);
  
  // Higher level ClassEngine methods (through ClikPar)
  void run_CLASS(std::vector<double> class_params);
  void get_CLASS_spectra(std::vector<double>& cltt, 
      std::vector<double>& clte, 
      std::vector<double>& clee, 
      std::vector<double>& clbb);

  void create_all_cl_and_pars(double* Cube,
      std::vector<std::vector<double> >& class_cls);

  // Likelihood functions
  double calculate_PLC_likelihood() const;
  double calculate_BAO_likelihood() const;
  
private:
  int m_largest_max_l;

  // CLASS bits
  std::vector<unsigned> m_class_l_vec;

  std::vector<ClikObject*> m_clik_objects;
  ClikPar* m_clik_par;
};

#endif
