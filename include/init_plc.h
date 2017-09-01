#ifndef PC_INIT
#define PC_INIT

#include "clik.h"
#include "ClassEngine.hh"

#include "pbh_io.h"

// #include "TTTEEE+lowP_pbh_fixedLCDM.h"
#include "TTTEEE+lowP_pbh.h"
// #include "TTTEEE+lowP.h"

#include <math.h>
#include <string>
#include <vector>

typedef double (*trans_t)(double);

const int CL_AMT = 6;
const int CLASS_MIN_L = 2;

extern double m_min[];
extern double m_max[];
extern double m_value[];
extern bool m_has_gaussian_prior[];
extern double m_mean[];
extern double m_stddev[];

extern bool m_is_log10[];
extern trans_t m_transform[];

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

void initialise_CLASS_engine(ClassEngine*& class_engine, int max_l, pbh_external* pbh_info);

pbh_external* initialise_pbh_external(std::string& pbh_root);

inline double pow10(double x) { return pow(10., x); }
inline double self(double x) { return x; }

void initialise_param_arrays();

#endif
