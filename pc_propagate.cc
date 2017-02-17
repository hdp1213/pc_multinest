#include "ClassEngine.hh"

#ifdef MPI
#include <mpi.h>
#endif

#include <cstdio> // for stderr

#include <iostream>
#include <vector>
#include <stdexcept> // gives various standard exceptions
// #include <math.h>

#include <fstream>
#include <sstream>
#include <iomanip>

struct ClikPar {
  enum param_t
  {
    llike,
    pmass,
    // Free parameters (CLASS)
    omega_b,
    omega_cdm,
    hundredxtheta_s,
    tau_reio,
    n_s,
    ln10_10_A_s,
    // Free parameters (PLC)
    A_cib_217,
    xi_sz_cib,
    A_sz,
    ps_A_100_100,
    ps_A_143_143,
    ps_A_143_217,
    ps_A_217_217,
    ksz_norm,
    gal545_A_100,
    gal545_A_143,
    gal545_A_143_217,
    gal545_A_217,
    calib_100T,
    calib_217T,
    A_planck,
    // Derived parameters
    cib_index, // cib_index = -1.3 is constant
    TOTAL_IN_PARAMS,
    // Derived parameters from sigma8
    H0 = TOTAL_IN_PARAMS,
    Omega_L,
    Omega_b,
    // sigma8,
    age,
    z_drag,
    rs_drag,
    // Others that aren't shared by Planck
    Omega_cdm,
    Omega_g,
    conf_age,
    TOTAL_OUT_PARAMS
  };
};

// Use CLASS to calculate derived parameters faster than Phoenix
int main(int argc, char* argv[])
{
  // Will need to use MPI
  int mpi_size = 1, mpi_rank = 0;

#ifdef MPI
  MPI_Comm mpi_comm;
  MPI_Status stat;

  MPI_Init(&argc, &argv);
  mpi_comm = MPI_COMM_WORLD;
  MPI_Comm_size(mpi_comm, &mpi_size);
  MPI_Comm_rank(mpi_comm, &mpi_rank);

#endif // MPI

  const int ROOT_RANK = 0;
  const int IN_ROWS = 379936; // from wc -l <in_file>
  //    32 or    31 nodes only!
  // 11873 or 12256 of size!
  // only 86 hours!
  // const int DERIVED_SIZE = ClikPar::TOTAL_OUT_PARAMS - ClikPar::TOTAL_IN_PARAMS;

  
  const char* phoenix_in;
  const char* derived_out;

  // Inputting variables
  std::ifstream in;
  double in_vals[ClikPar::TOTAL_IN_PARAMS];

  // Outputting variables
  std::ofstream out;
  int width = 28;
  int precis = 18;

  // Collecting and storing input variables
  int IN_SIZE = ClikPar::TOTAL_IN_PARAMS * IN_ROWS;
  double* in_params;

  // Variables for Hard Work
  int WORK_SIZE, WORK_ROWS;
  double* working_params;

  // Variables to collect Hard Work
  int WORKED_SIZE;
  double* worked_params;

  // Variables to gather Hard Work back into
  int OUT_SIZE = ClikPar::TOTAL_OUT_PARAMS * IN_ROWS;
  double* out_params;

  // CLASS initialisation variables
  ClassParams default_params;
  int max_l = 2509;
  ClassEngine* class_engine(0);
  std::vector<double>* lcdm_params(0);

  // Read in/out command line arguments
  if (argc < 1) {
    phoenix_in = "/data/harryp/phoenix-out/pc_multinest_mpi_full_lensed_precision_6-.txt";
  }
  else {
    phoenix_in = argv[1];
  }

  if (argc < 2) {
    derived_out = "/data/harryp/phoenix-out/pc_multinest_full_derived-.txt";
  }
  else {
    derived_out = argv[2];
  }

  // Run checks on divisibility of row number by mpi processes
  if (mpi_rank == ROOT_RANK) {
    if (IN_ROWS % mpi_size != 0) {
      std::cerr << "[ERROR]: You've chosen " << mpi_size 
                << " MPI processes. "
                << "Please choose one that divides "
                << IN_ROWS << std::endl;
      throw std::exception();
    }
  }

  WORK_SIZE = IN_SIZE/mpi_size;
  WORK_ROWS = WORK_SIZE/ClikPar::TOTAL_IN_PARAMS;
  WORKED_SIZE = ClikPar::TOTAL_OUT_PARAMS * WORK_ROWS;

  //*
  // Initialise CLASS engine
  // All nodes will need CLASS to compute things!
  default_params.add("omega_b", 0.022032);
  default_params.add("omega_cdm", 0.12038);
  default_params.add("100*theta_s", 1.042143);
  default_params.add("tau_reio", 0.0925);
  default_params.add("ln10^{10}A_s", 3.0980);
  default_params.add("n_s", 0.9619);

  default_params.add("N_ur", 2.0328);
  default_params.add("N_ncdm", 1);
  default_params.add("m_ncdm", 0.06); // MeV
  // default_params.add("T_ncdm", 0.71611);

  // default_params.add("P_k_max_h/Mpc", 1.);
  // default_params.add("k_pivot", 0.05); // Mpc-1
  
  default_params.add("output", "tCl");
  // default_params.add("output", "tCl,pCl,lCl,mPk");
  // default_params.add("lensing", true);   //note boolean
  default_params.add("l_max_scalars", max_l);
  // default_params.add("format", "camb");

  try {
    class_engine = new ClassEngine(default_params);
  }
  catch (std::exception const &e) {
    std::cerr << "[ERROR] CLASS initialisation failed "
              << "on MPI thread mpi_rank '" << mpi_rank
              << "', throwing exception "
              << e.what()
              << std::endl;
    throw e;
  }
  //*/

  // Only read input if you are root process
  if (mpi_rank == ROOT_RANK) {
    std::cout << "Starting reading from " << phoenix_in << std::endl;

    in_params = new double[IN_SIZE];

    // Just read in values to begin with
    in.open(phoenix_in);

    int line_num = 0;
    while (true) {
      in >> in_vals[ClikPar::pmass]
         >> in_vals[ClikPar::llike]
         >> in_vals[ClikPar::omega_b]
         >> in_vals[ClikPar::omega_cdm]
         >> in_vals[ClikPar::hundredxtheta_s]
         >> in_vals[ClikPar::tau_reio]
         >> in_vals[ClikPar::n_s]
         >> in_vals[ClikPar::ln10_10_A_s]
         >> in_vals[ClikPar::A_cib_217]
         >> in_vals[ClikPar::xi_sz_cib]
         >> in_vals[ClikPar::A_sz]
         >> in_vals[ClikPar::ps_A_100_100]
         >> in_vals[ClikPar::ps_A_143_143]
         >> in_vals[ClikPar::ps_A_143_217]
         >> in_vals[ClikPar::ps_A_217_217]
         >> in_vals[ClikPar::ksz_norm]
         >> in_vals[ClikPar::gal545_A_100]
         >> in_vals[ClikPar::gal545_A_143]
         >> in_vals[ClikPar::gal545_A_143_217]
         >> in_vals[ClikPar::gal545_A_217]
         >> in_vals[ClikPar::calib_100T]
         >> in_vals[ClikPar::calib_217T]
         >> in_vals[ClikPar::A_planck]
         >> in_vals[ClikPar::cib_index];

      if (!in.good()) break;

      // Push all parameters that have been read "in" to in_params.
      // in_params is a 1D vector, but we know that it "wraps" its
      // values every TOTAL_IN_PARAMS values
      for (int i = 0; i < ClikPar::TOTAL_IN_PARAMS; ++i) {
        in_params[line_num*ClikPar::TOTAL_IN_PARAMS + i] = in_vals[i];
      }

      line_num++;
    }

    in.close();

    std::cout << "Finished reading " << phoenix_in << std::endl;
  }

  else { // be safe
    in_params = NULL;
  }

#ifdef MPI
  // Where we store the scattered seeds
  working_params = new double[WORK_SIZE];

  // Scatter the seeds...
  MPI_Scatter(in_params,
    WORK_SIZE,
    MPI_DOUBLE,
    working_params,
    WORK_SIZE,
    MPI_DOUBLE,
    ROOT_RANK,
    mpi_comm);
#else
  working_params = in_params;
#endif // MPI

  // Where we put all the new goodies
  worked_params = new double[WORKED_SIZE];

  // Now, need to run CLASS into the fucking ground. Daddy, why
  for (int row = 0; row < WORK_ROWS; ++row) {
    int base = row * ClikPar::TOTAL_IN_PARAMS;

    // Update parameter vector (new one every time!)
    lcdm_params = new std::vector<double>();

    lcdm_params->push_back(working_params[base + ClikPar::omega_b]);
    lcdm_params->push_back(working_params[base + ClikPar::omega_cdm]);
    lcdm_params->push_back(working_params[base + ClikPar::hundredxtheta_s]);
    lcdm_params->push_back(working_params[base + ClikPar::tau_reio]);
    lcdm_params->push_back(working_params[base + ClikPar::ln10_10_A_s]);
    lcdm_params->push_back(working_params[base + ClikPar::n_s]); // k_0 = 0.05 Mpc^-1 by default

    // Run CLASS for me baby
    try {
      class_engine->updateParValues(*lcdm_params);
    }
    catch (std::exception const &e) {
      std::cerr << "Fuck!" << std::endl;
    }

    delete lcdm_params;

    // Copy over Known Parameters from working to worked
    for (int i = 0; i < ClikPar::TOTAL_IN_PARAMS; ++i) {
      worked_params[base + i] = working_params[base + i];
    }

    // Extract goodies
    worked_params[base + ClikPar::H0] = class_engine->get_H0();
    worked_params[base + ClikPar::Omega_b] = class_engine->get_Omega_b();
    worked_params[base + ClikPar::Omega_cdm] = class_engine->get_Omega_cdm();
    worked_params[base + ClikPar::Omega_L] = class_engine->get_Omega_L();
    worked_params[base + ClikPar::Omega_g] = class_engine->get_Omega_g();
    // worked_params[base + ClikPar::sigma8] = class_engine->get_sigma8();
    worked_params[base + ClikPar::age] = class_engine->get_age();
    worked_params[base + ClikPar::conf_age] = class_engine->get_conf_age();

    worked_params[base + ClikPar::z_drag] = class_engine->z_drag();
    worked_params[base + ClikPar::rs_drag] = class_engine->rs_drag();
  }

#ifdef MPI
  // It is done, thank fuck. Let's gather the results
  if (mpi_rank == ROOT_RANK) {
    out_params = new double[OUT_SIZE];
  }
  else {
    out_params = NULL;
  }

  MPI_Gather(worked_params,
    WORKED_SIZE,
    MPI_DOUBLE,
    out_params,
    WORKED_SIZE,
    MPI_DOUBLE,
    ROOT_RANK,
    mpi_comm);
#else
  out_params = worked_params;
#endif

  // Now print all of this fuckery to a file
  // Only write output if you are root process
  if (mpi_rank == ROOT_RANK) {
    std::cout << "Starting writing to " << derived_out << std::endl;

    out.open(derived_out);

    for (int i = 0; i < OUT_SIZE; ++i) {
      out << std::setw(width) << std::setprecision(precis) 
          << out_params[i];

      if (i % ClikPar::TOTAL_OUT_PARAMS == 0) {
        out << std::endl;
      }
    }

    out.close();

    std::cout << "Finished writing to " << derived_out << std::endl;
  }

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}