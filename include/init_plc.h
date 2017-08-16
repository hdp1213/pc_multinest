#ifndef PC_INIT
#define PC_INIT

#include "clik.h"
#include "ClassEngine.hh"

// #include "TTTEEE+lowP_pbh_fixedLCDM.h"
// #include "TTTEEE+lowP_pbh.h"
#include "TTTEEE+lowP.h"

#include <string>
#include <vector>

const int CL_AMT = 6;
const int CLASS_MIN_L = 2;

struct clik_struct {
  clik_object* clik_id;
  int max_l;
  int cap_size;
  bool has_cl[CL_AMT];
  int nuis_par_amt;
  std::vector<param_t> nuis_pars;
};

struct plc_bundle {
  std::vector<clik_struct*> clik_objs;
  ClassEngine* engine;
  std::vector<unsigned> cl_ls;
};

clik_struct* initialise_clik_struct(std::string& clik_path,
                                    std::vector<param_t>& nuis_params,
                                    int& total_max_l);

void initialise_CLASS_engine(ClassEngine*& class_engine, int max_l);

#endif
