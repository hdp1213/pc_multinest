//--------------------------------------------------------------------------
//
// Description:
//  class ClassEngine : see header file (ClassEngine.hpp) for description.
//
//------------------------------------------------------------------------
//-----------------------
// This Class's Header --
//-----------------------
#include "ClassEngine.hpp"
// ------------------------
// Collaborating classes --
//-------------------------
// C++
//--------------------
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

template<typename T> std::string str(const T &x) {
  std::ostringstream os;
  os << x;
  return os.str();
}

template<> std::string str (const float &x) {
  std::ostringstream os;
  os << std::setprecision(8) << x;
  return os.str();
}

template<> std::string str (const double &x) {
  std::ostringstream os;
  os << std::setprecision(16) << x;
  return os.str();
}

template<> std::string str (const bool &x) {
  return x ? "yes" : "no";
}

template<> std::string str (const std::string &x) {
  return x;
}

std::string str (const char* s) {
  return std::string(s);
}

// instanciations
template std::string str(const int &x);
template std::string str(const signed char &x);
template std::string str(const unsigned char &x);
template std::string str(const short &x);
template std::string str(const unsigned short &x);
template std::string str(const unsigned int &x);
template std::string str(const long &x);
template std::string str(const unsigned long &x);
template std::string str(const long long &x);
template std::string str(const unsigned long long &x);

//---------------
// Constructors --
//----------------
ClassEngine::ClassEngine(const ClassParams& pars): m_cl(0), m_do_free(true), m_info(NULL) {
  write_pars_to_fc(pars, &m_fc);

  // Initialise input
  if (input_init(&m_fc, &m_pr, &m_ba, &m_th, &m_pt, &m_tr, &m_pm, &m_sp, &m_nl, &m_le, &m_op, m_info, m_errmsg) == _FAILURE_) {
    throw std::invalid_argument(m_errmsg);
  }

  // Protection against invalid parameters that haven't been read
  for (std::size_t i = 0; i < pars.size(); i++) {
    if (m_fc.read[i] != _TRUE_) {
      throw std::invalid_argument(std::string("Invalid CLASS parameter: ") + m_fc.name[i]);
    }
  }

  // Run CLASS
  int status;
  status = compute_Cls();

  if (status == _FAILURE_) {
    throw std::out_of_range(m_errmsg);
  }

#ifdef DBUG
  std::cout << "Creating " << m_sp.ct_size << " arrays" << std::endl;
#endif
  m_cl = new double[m_sp.ct_size];

#ifdef DBUG
  print_FC();
#endif
}


ClassEngine::ClassEngine(const ClassParams& pars, const std::string & precision_file): m_cl(0), m_do_free(true), m_info(NULL) {

  // Decode precision structure
  struct file_content fc_precision;
  fc_precision.size = 0;
  if (parser_read_file(const_cast<char*>(precision_file.c_str()), &fc_precision, m_errmsg) == _FAILURE_) {
    throw std::invalid_argument(m_errmsg);
  }

  // Read input params into file_content struct
  struct file_content fc_input;
  fc_input.size = 0;
  fc_input.filename = new char[1];
  write_pars_to_fc(pars, &fc_input);

  // concatenate bom_th
  if (parser_cat(&fc_input, &fc_precision, &m_fc, m_errmsg) == _FAILURE_) {
    throw std::invalid_argument(m_errmsg);
  }

  parser_free(&fc_input);
  parser_free(&fc_precision);

  // Initialise input
  if (input_init(&m_fc, &m_pr, &m_ba, &m_th, &m_pt, &m_tr, &m_pm, &m_sp, &m_nl, &m_le, &m_op, m_info, m_errmsg) == _FAILURE_) {
    throw std::invalid_argument(m_errmsg);
  }

  // Protection against invalid parameters that haven't been read
  for (std::size_t i = 0; i < pars.size(); i++) {
    if (m_fc.read[i] != _TRUE_) {
      throw std::invalid_argument(std::string("Invalid CLASS parameter: ") + m_fc.name[i]);
    }
  }

  // Run CLASS
  int status;
  status = compute_Cls();

  if (status == _FAILURE_) {
    throw std::out_of_range(m_errmsg);
  }

#ifdef DBUG
  std::cout << "Creating " << m_sp.ct_size << " arrays" << std::endl;
#endif
  m_cl = new double[m_sp.ct_size];

#ifdef DBUG
  print_FC();
#endif
}

ClassEngine::ClassEngine(const ClassParams& pars, struct external_info* info): m_cl(0), m_do_free(true), m_info(info) {
  // Prepare fp structure
  write_pars_to_fc(pars, &m_fc);

  // Initialise input
  if (input_init(&m_fc, &m_pr, &m_ba, &m_th, &m_pt, &m_tr, &m_pm, &m_sp, &m_nl, &m_le, &m_op, m_info, m_errmsg) == _FAILURE_) {
    throw std::invalid_argument(m_errmsg);
  }

  // Protection against invalid parameters that haven't been read
  for (std::size_t i = 0; i < pars.size(); i++) {
    if (m_fc.read[i] != _TRUE_) {
      throw std::invalid_argument(std::string("invalid CLASS parameter: ") + m_fc.name[i]);
    }
  }

  // Run CLASS
  int status;
  status = compute_Cls();

  if (status == _FAILURE_) {
    throw std::out_of_range(m_errmsg);
  }

#ifdef DBUG
  std::cout << "Creating " << m_sp.ct_size << " arrays" << std::endl;
#endif
  m_cl = new double[m_sp.ct_size];

#ifdef DBUG
  print_FC();
#endif
}

ClassEngine::ClassEngine(const std::string& init_file, int l_max, struct external_info* info): Engine(l_max), m_cl(0), m_do_free(true), m_info(info) {

  // Variables
  std::size_t i;
  std::string lmax_str = str(l_max);
  bool found_lmax = false;

  // Pars
  m_fc.size = 0;
  // Decode init structure
  if (parser_read_file(const_cast<char*>(init_file.c_str()), &m_fc, m_errmsg) == _FAILURE_) {
    throw std::invalid_argument(m_errmsg);
  }

  // Config
  for (i = 0; i < m_fc.size; i++) {
    // Store
    m_parnames.push_back(m_fc.name[i]);
    // Identify if lmax is given in init_file, and override
    if (strcmp(m_fc.name[i], "l_max_scalars") == 0) {
      strcpy(m_fc.value[i], lmax_str.c_str());
      found_lmax = true;
      // istringstream strstrm(fc.value[i]);
      // strstrm >> _lmax;
    }
    std::cout << i << ": " << m_fc.name[i] << "\t" << m_fc.value[i] << std::endl;
  }

  if (found_lmax) {
    std::cout << "Overriding l_max_scalars value found in init file..." << std::endl;
  }
  else {
    strcpy(m_fc.name[i], "l_max_scalars");
    strcpy(m_fc.value[i], lmax_str.c_str());
  }

  std::cout << __FILE__ << " : using lmax = " << m_lmax << std::endl;
  assert(m_lmax>0);

  // Initialise input
  if (input_init(&m_fc, &m_pr, &m_ba, &m_th, &m_pt, &m_tr, &m_pm, &m_sp, &m_nl, &m_le, &m_op, m_info, m_errmsg) == _FAILURE_) {
    throw std::invalid_argument(m_errmsg);
  }

  // Protection against invalid parameters that haven't been read
  for (i = 0; i < m_fc.size; i++) {
    if (m_fc.read[i] != _TRUE_) {
      throw std::invalid_argument(std::string("invalid CLASS parameter: ") + m_fc.name[i]);
    }
  }

  // Run CLASS
  int status;
  status = compute_Cls();

  if (status == _FAILURE_) {
    throw std::out_of_range(m_errmsg);
  }

#ifdef DBUG
  std::cout << "Creating " << m_sp.ct_size << " arrays" << std::endl;
#endif
  m_cl = new double[m_sp.ct_size];

#ifdef DBUG
  print_FC();
#endif
}



//--------------
// Destructor --
//--------------
ClassEngine::~ClassEngine() {
#ifdef DBUG
  print_FC();
#endif
  std::cout << "Deleting CLASS..." << std::endl;
  m_do_free && free_structs();

  parser_free(&m_fc);

  delete[] m_cl;
}

//-----------------
// Member functions --
//-----------------
void
ClassEngine::update_parameters(const std::vector<double>& par) {
  m_do_free && free_structs();

  for (std::size_t i = 0; i < par.size(); i++) {
    double val = par[i];
    strcpy(m_fc.value[i], str(val).c_str());
#ifdef DBUG
    std::cout << "Update par values " << m_parnames[i] << "\t" <<  val << "\t" << str(val).c_str() << std::endl;
#endif
  }

  int status = compute_Cls();

#ifdef DBUG
  std::cout << "update_parameters status = " << status << ", success = " << _SUCCESS_ << std::endl;
#endif

  if (status == _FAILURE_) {
    throw std::out_of_range(m_errmsg);
  }
}

// Print content of file_content
void
ClassEngine::print_FC() {
  std::cout << "FILE_CONTENT SIZE = " << m_fc.size << std::endl;
  for (int i = 0; i < m_fc.size; i++) {
    printf("%d : %s = %s\n", i, m_fc.name[i], m_fc.value[i]);
  }
}

int
ClassEngine::class_main(struct file_content *pfc,
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
                        ErrorMsg errmsg) {


  if (input_init(pfc, ppr, pba, pth, ppt, ptr, ppm, psp, pnl, ple, pop, m_info, errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    m_do_free = false;
    return _FAILURE_;
  }

  if (background_init(ppr, pba) == _FAILURE_) {
    strcpy(m_errmsg, pba->error_message);
    background_free(&m_ba);
    m_do_free = false;
    return _FAILURE_;
  }

  if (thermodynamics_init(ppr, pba, pth, m_info) == _FAILURE_) {
    strcpy(m_errmsg, pth->error_message);
    thermodynamics_free(&m_th);
    background_free(&m_ba);
    m_do_free = false;
    return _FAILURE_;
  }

  if (perturb_init(ppr, pba, pth, ppt) == _FAILURE_) {
    strcpy(m_errmsg, ppt->error_message);
    perturb_free(&m_pt);
    thermodynamics_free(&m_th);
    background_free(&m_ba);
    m_do_free = false;
    return _FAILURE_;
  }

  if (primordial_init(ppr, ppt, ppm) == _FAILURE_) {
    strcpy(m_errmsg, ppm->error_message);
    primordial_free(&m_pm);
    perturb_free(&m_pt);
    thermodynamics_free(&m_th);
    background_free(&m_ba);
    m_do_free = false;
    return _FAILURE_;
  }

  if (nonlinear_init(ppr, pba, pth, ppt, ppm, pnl) == _FAILURE_)  {
    strcpy(m_errmsg, pnl->error_message);
    nonlinear_free(&m_nl);
    primordial_free(&m_pm);
    perturb_free(&m_pt);
    thermodynamics_free(&m_th);
    background_free(&m_ba);
    m_do_free = false;
    return _FAILURE_;
  }

  if (transfer_init(ppr, pba, pth, ppt, pnl, ptr) == _FAILURE_) {
    strcpy(m_errmsg, ptr->error_message);
    transfer_free(&m_tr);
    nonlinear_free(&m_nl);
    primordial_free(&m_pm);
    perturb_free(&m_pt);
    thermodynamics_free(&m_th);
    background_free(&m_ba);
    m_do_free = false;
    return _FAILURE_;
  }

  if (spectra_init(ppr, pba, ppt, ppm, pnl, ptr, psp) == _FAILURE_) {
    strcpy(m_errmsg, psp->error_message);
    spectra_free(&m_sp);
    transfer_free(&m_tr);
    nonlinear_free(&m_nl);
    primordial_free(&m_pm);
    perturb_free(&m_pt);
    thermodynamics_free(&m_th);
    background_free(&m_ba);
    m_do_free = false;
    return _FAILURE_;
  }

  if (lensing_init(ppr, ppt, psp, pnl, ple) == _FAILURE_) {
    strcpy(m_errmsg, ple->error_message);
    lensing_free(&m_le);
    spectra_free(&m_sp);
    transfer_free(&m_tr);
    nonlinear_free(&m_nl);
    primordial_free(&m_pm);
    perturb_free(&m_pt);
    thermodynamics_free(&m_th);
    background_free(&m_ba);
    m_do_free = false;
    return _FAILURE_;
  }


  m_do_free = true;
  return _SUCCESS_;
}


int
ClassEngine::compute_Cls() {
#ifdef DBUG
  std::cout << "call compute_Cls" << std::endl;
  print_FC();
#endif

  int status = this->class_main(&m_fc, &m_pr, &m_ba, &m_th, &m_pt, &m_tr, &m_pm, &m_sp, &m_nl, &m_le, &m_op, m_errmsg);

#ifdef DBUG
  std::cout << "status = " << status << std::endl;
#endif
  return status;
}

int
ClassEngine::free_structs() {
  if (lensing_free(&m_le) == _FAILURE_) {
    strcpy(m_errmsg, m_le.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&m_nl) == _FAILURE_) {
    strcpy(m_errmsg, m_nl.error_message);
    return _FAILURE_;
  }

  if (spectra_free(&m_sp) == _FAILURE_) {
    strcpy(m_errmsg, m_sp.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&m_pm) == _FAILURE_) {
    strcpy(m_errmsg, m_pm.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&m_tr) == _FAILURE_) {
    strcpy(m_errmsg, m_tr.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&m_pt) == _FAILURE_) {
    strcpy(m_errmsg, m_pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&m_th) == _FAILURE_) {
    strcpy(m_errmsg, m_th.error_message);
    return _FAILURE_;
  }

  if (background_free(&m_ba) == _FAILURE_) {
    strcpy(m_errmsg, m_ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;
}

// Also populates m_parnames
void
ClassEngine::write_pars_to_fc(const ClassParams& pars, struct file_content* fc) {
  // Prepare file content structure
  std::size_t n = pars.size();
  parser_init(&m_fc, n, "pipo", m_errmsg);

  // Populate m_fc with parameters from pars
  for (std::size_t i = 0; i < pars.size(); i++) {
    strcpy(m_fc.name[i], pars.key(i).c_str());
    strcpy(m_fc.value[i], pars.value(i).c_str());

    // Store names separately
    m_parnames.push_back(pars.key(i));

    std::cout << pars.key(i) << "\t" << pars.value(i) << std::endl;

    // Identify lmax
    if (pars.key(i) == "l_max_scalars") {
      std::istringstream strstrm(pars.value(i));
      strstrm >> m_lmax;
    }
  }

  std::cout << __FILE__ << " : using lmax = " << m_lmax << std::endl;
  assert(m_lmax>0);
}

double
ClassEngine::get_Cl_value_at(const long &l, Engine::cltype t) {

  if (!m_do_free) {
    throw std::out_of_range("No Cl available because CLASS failed");
  }

  if (output_total_cl_at_l(&m_sp, &m_le, &m_op, static_cast<double>(l), m_cl) == _FAILURE_) {
    std::cerr << "[ERROR] get_Cl_value_at(): Failed getting Cl type " << t << " at l = " << l << std::endl;
    throw std::out_of_range(m_sp.error_message);
  }

  double zecl = -1;

  double T0_muK = 1e6 * get_T_cmb();
  double T0_muK2 = T0_muK * T0_muK;

  switch(t) {
    case TT:
      (m_sp.has_tt == _TRUE_) ? zecl = T0_muK2*m_cl[m_sp.index_ct_tt] : throw std::invalid_argument("No Cl_TT available");
      break;
    case TE:
      (m_sp.has_te == _TRUE_) ? zecl = T0_muK2*m_cl[m_sp.index_ct_te] : throw std::invalid_argument("No Cl_TE available");
      break;
    case EE:
      (m_sp.has_ee == _TRUE_) ? zecl = T0_muK2*m_cl[m_sp.index_ct_ee] : throw std::invalid_argument("No Cl_EE available");
      break;
    case BB:
      (m_sp.has_bb == _TRUE_) ? zecl = T0_muK2*m_cl[m_sp.index_ct_bb] : throw std::invalid_argument("No Cl_BB available");
      break;
    case PP:
      (m_sp.has_pp == _TRUE_) ? zecl = m_cl[m_sp.index_ct_pp] : throw std::invalid_argument("No Cl_PhiPhi available");
      break;
    case TP:
      (m_sp.has_tp == _TRUE_) ? zecl = T0_muK*m_cl[m_sp.index_ct_tp] : throw std::invalid_argument("No Cl_TPhi available");
      break;
    case EP:
      (m_sp.has_ep == _TRUE_) ? zecl = T0_muK*m_cl[m_sp.index_ct_ep] : throw std::invalid_argument("No Cl_EPhi available");
      break;
    default:
      throw std::invalid_argument("invalid spectra");
    }

  return zecl;
}

void
ClassEngine::get_Cls(const std::vector<unsigned>& lvec, // input
                     std::vector<double>& cltt,
                     std::vector<double>& clte,
                     std::vector<double>& clee,
                     std::vector<double>& clbb) {

  cltt.resize(lvec.size());
  clte.resize(lvec.size());
  clee.resize(lvec.size());
  clbb.resize(lvec.size());

  for (std::size_t i = 0; i < lvec.size(); i++) {
    try {
      cltt[i] = get_Cl_value_at(lvec[i], Engine::TT);
      clte[i] = get_Cl_value_at(lvec[i], Engine::TE);
      clee[i] = get_Cl_value_at(lvec[i], Engine::EE);
      clbb[i] = get_Cl_value_at(lvec[i], Engine::BB);
    }
    catch (std::exception &e) {
      throw;
    }
  }

}

bool
ClassEngine::get_lensing_Cls(const std::vector<unsigned>& lvec, // input
                        std::vector<double>& clpp,
                        std::vector<double>& cltp,
                        std::vector<double>& clep) {

  clpp.resize(lvec.size());
  cltp.resize(lvec.size());
  clep.resize(lvec.size());

  for (std::size_t i = 0; i < lvec.size(); i++) {
    try {
      clpp[i] = get_Cl_value_at(lvec[i], Engine::PP);
      cltp[i] = get_Cl_value_at(lvec[i], Engine::TP);
      clep[i] = get_Cl_value_at(lvec[i], Engine::EP);
    }
    catch(std::exception &e) {
      std::cout << "plantage!" << std::endl;
      std::cout << __FILE__ << e.what() << std::endl;
      return false;
    }
  }
  return true;
}


double ClassEngine::get_f(double z)
{
  double tau;
  int index;
  double *pvecback;
  // Transform redshift in conformal time
  background_tau_of_z(&m_ba, z, &tau);

  // pvecback must be allocated
  pvecback = new double[m_ba.bg_size];

  // Call to fill pvecback
  background_at_tau(&m_ba, tau, m_ba.long_info, m_ba.inter_normal, &index, pvecback);

  double f_z = pvecback[m_ba.index_bg_f];

  delete[] pvecback;

#ifdef DBUG
  std::cout << "f_of_z = "<< f_z << std::endl;
#endif
  return f_z;
}


double ClassEngine::get_sigma8(double z)
{
  double tau;
  int index;
  double *pvecback;
  double sigma8 = 0.;
  // Transform redshift in conformal time
  background_tau_of_z(&m_ba, z, &tau);

  // pvecback must be allocated
  pvecback = new double[m_ba.bg_size];

  // Call to fill pvecback
  background_at_tau(&m_ba,tau,m_ba.long_info,m_ba.inter_normal, &index, pvecback);
  // background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback);
  spectra_sigma(&m_ba,&m_pm,&m_sp,8./m_ba.h,z,&sigma8);

  delete[] pvecback;

#ifdef DBUG
  std::cout << "sigma_8 = "<< sigma8 << std::endl;
#endif
  return sigma8;
}

// ATTENTION FONCTION BIDON - GET omegam ! -------------------
double ClassEngine::get_Az(double z)
{
  double Dv = get_Dv(z);
  // A(z) = 100DV(z)sqrt(~mh2)/cz
  double omega_bidon = 0.12 ;
  double Az = 100.*Dv*sqrt(omega_bidon)/(_c_ * z);

  return Az;
}
//      --------------------------

double ClassEngine::get_Dv(double z)
{
  double tau;
  int index;
  double *pvecback;
  // Transform redshift in conformal time
  background_tau_of_z(&m_ba,z,&tau);

  // pvecback must be allocated
  pvecback = new double[m_ba.bg_size];

  // Call to fill pvecback
  background_at_tau(&m_ba,tau,m_ba.long_info,m_ba.inter_normal, &index, pvecback);


  double H_z = pvecback[m_ba.index_bg_H];
  double D_ang = pvecback[m_ba.index_bg_ang_distance];

  delete[] pvecback;
#ifdef DBUG
  std::cout << "H_z = "<< H_z << std::endl;
  std::cout << "D_ang = "<< D_ang << std::endl;
#endif
  double D_v;

  D_v = pow(D_ang*(1+z),2)*z/H_z; // H_z is given in Mpc^-1 (i.e. H_z * c)
  D_v = pow(D_v,1./3.);
#ifdef DBUG
  std::cout << D_v << std::endl;
#endif
  return D_v;
}

double ClassEngine::get_Fz(double z)
{
  double tau;
  int index;
  double *pvecback;
  // Transform redshift in conformal time
  background_tau_of_z(&m_ba,z,&tau);

  // pvecback must be allocated
  pvecback = new double[m_ba.bg_size];

  // Call to fill pvecback
  background_at_tau(&m_ba,tau,m_ba.long_info,m_ba.inter_normal, &index, pvecback);


  double H_z = pvecback[m_ba.index_bg_H];
  double D_ang = pvecback[m_ba.index_bg_ang_distance];

  delete[] pvecback;

#ifdef DBUG
  std::cout << "H_z = "<< H_z << std::endl;
  std::cout << "D_ang = "<< D_ang << std::endl;
#endif
  double F_z = (1.+z) * D_ang * H_z /_c_;
  return F_z;
}

double ClassEngine::get_Hz(double z)
{
  double tau;
  int index;
  double *pvecback;
  // Transform redshift in conformal time
  background_tau_of_z(&m_ba,z,&tau);

  // pvecback must be allocated
  pvecback = new double[m_ba.bg_size];

  // Call to fill pvecback
  background_at_tau(&m_ba,tau,m_ba.long_info,m_ba.inter_normal, &index, pvecback);

  double H_z = pvecback[m_ba.index_bg_H];

  delete[] pvecback;

  return H_z;
}


double ClassEngine::get_Da(double z)
{
  double tau;
  int index;
  double *pvecback;
  // Transform redshift in conformal time
  background_tau_of_z(&m_ba,z,&tau);

  // pvecback must be allocated
  pvecback = new double[m_ba.bg_size];

  // Call to fill pvecback
  background_at_tau(&m_ba,tau,m_ba.long_info,m_ba.inter_normal, &index, pvecback);

#ifdef DBUG
  double H_z = pvecback[m_ba.index_bg_H];
#endif
  double D_ang = pvecback[m_ba.index_bg_ang_distance];

  delete[] pvecback;

#ifdef DBUG
  std::cout << "H_z = " << H_z << std::endl;
  std::cout << "D_ang = " << D_ang << std::endl;
#endif
  return D_ang;
}
