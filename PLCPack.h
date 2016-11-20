// Header file for plc_pack class
#ifndef PLCPACK_H
#define PLCPACK_H

#include "ClikObject.h"
#include <vector>


class PLCPack {

public:
  PLCPack();
  ~PLCPack();

  void add_clik_object(ClikObject *clik_object);
  int get_largest_max_l() const;

  void create_all_cl_and_pars(double* Cube,
      std::vector<std::vector<double> >& class_cls);

  double calculate_likelihood() const;
  
  std::vector<unsigned> get_class_l_vec() const;

private:
  int m_largest_max_l;

  // CLASS bits
  std::vector<unsigned> m_class_l_vec;

  std::vector<ClikObject*> m_clik_objects;
};

#endif
