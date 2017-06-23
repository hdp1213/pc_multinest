// Header file for plc_pack class
#ifndef CLIKOBJECT_H
#define CLIKOBJECT_H

#include "clik.h"
#include "ClikPar.h"
#include <vector>
#include <string>


class ClikObject {

public:
  ClikObject(char* clik_path);
  ~ClikObject();

  int get_max_l() const;
  int get_cap_size() const;
  int* get_cl_flags();
  int get_param_amt() const;
  parname* get_param_names() const;

  // void set_cl_and_pars(double* cl_and_pars);
  void set_nuisance_param_enums(std::vector<ClikPar::param_t>& nuisance_params);
  void create_cl_and_pars(double* Cube,
        std::vector<std::vector<double> >& class_cls);

  double get_likelihood() const;


private:
  int incompatible_with_class() const;

  // Fixed max number of spectra that PLC can store in .clik
  // files
  static const int CL_AMT = 6;
  
  // PLC bits
  error* _m_err;
  error** m_err;
  clik_object* m_clik_id;

  char* m_clik_path;

  // Other things to store from clik object
  int m_max_l;
  int m_cap_size; // Size of cl_and_pars array
  int m_cl_max_ls[CL_AMT];
  int m_cl_flags[CL_AMT];
  int m_param_amt;
  parname* m_param_names;
  double* m_cl_and_pars;

  std::vector<ClikPar::param_t> m_nuis_par_enums;

};

#endif
