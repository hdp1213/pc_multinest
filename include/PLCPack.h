// Header file for plc_pack class
#ifndef PLCPACK_H
#define PLCPACK_H

#include "ClassEngine.hh"
#include "ClikObject.h"
#include "ClikPar.h"

#include <vector>
#include <string>


class PLCPack {

public:
  PLCPack();
  ~PLCPack();

  // Pack "building" methods
  void add_clik_object(ClikObject* clik_object);
  void read_pbh_files(std::string pbh_root);
  void set_clik_params(ClikPar* clik_par);
  void initialise_CLASS();

  // Wrapper methods for ClikPar object
  double calculate_extra_likelihoods(double* Cube) const;
  void scale_Cube(double* Cube);
  void set_derived_params(double* Cube);
  double calculate_prior() const;
  
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
#ifdef BAO_LIKE
  double calculate_BAO_likelihood() const;
#endif
  
private:
  int m_largest_max_l;

  // CLASS bits
  std::vector<unsigned> m_class_l_vec;

  std::vector<ClikObject*> m_clik_objects;
  ClikPar* m_clik_par;

  // PBH bits
  struct pbh_external m_pbh_info;

  // PBH methods
  void read_axes(std::string root);
  void read_bicubic_bspline(std::string root, const char* channel, struct bspline_2d* spline);
  void read_1d_array(std::ifstream& file, double** array, int* array_size);
};

#endif
