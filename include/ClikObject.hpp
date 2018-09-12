#ifndef CLIK_OBJECT_H
#define CLIK_OBJECT_H

#include "clik.h"

#include "Parameter.hpp"

#include <string>
#include <utility>
#include <vector>

const unsigned CL_AMT = 6;
const unsigned CLASS_MIN_L = 2;
const unsigned CLASS_CL_AMT = 4;

class ClikObject
{
public:
  ClikObject(const std::string& clik_file,
             std::vector<param_t>& nuis_params);
  ~ClikObject();

  unsigned get_lmax() const { return m_lmax; };

  double calculate_likelihood(double input_vals[],
                              std::vector<double> class_cls[]) const;

private:
  clik_object* m_clik_id;
  std::vector<param_t> m_nuis_params;
  unsigned m_lmax, m_cap_size;
  bool m_has_cl[CL_AMT];
};

#endif
