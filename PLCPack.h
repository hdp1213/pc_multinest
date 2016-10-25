// Header file for plc_pack class
#include "clik.h"


class PLCPack {

public:
  PLCPack(clik_object *clik_id);
  ~PLCPack();

  int get_max_l() const;
  int get_cap_size() const;
  int* get_cl_flags();
  int get_param_amt() const;
  // parname get_param_names() const;

  double get_clik_likelihood(double *cl_and_pars);

private:
  int invalid_clik_object() const;

  // Fixed max number of spectra that PLC can store in .clik files
  static const int CL_AMT = 6;

  // PLC bits
  error *_m_err, **m_err;
  clik_object *m_clik_id;

  // Other things to store from clik object
  int m_max_l;
  int m_cap_size; // Size of cl_and_pars array
  int m_cl_max_ls[CL_AMT];
  int m_cl_flags[CL_AMT];

  int m_param_amt;
  parname *m_param_names;

};