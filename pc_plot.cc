struct ClikPar {
  enum param_t
  {
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
    // Final total number of parameters
    TOTAL_PARAMS
  };
};

void print_mean(TH1D *hist, const char * title) {
  cout << title << " = " << hist->GetMean() << " +- " << hist->GetStdDev() << endl;
}

void pc_plot() {
  ifstream in;
  double pmass, like;

  int nparams = 22;
  const char* names[nparams];
  double params[nparams];
  double xlow[nparams];
  double xupp[nparams];
  TH1D * hists[nparams];

  const char* file_name = "output/pc_multinest_mpi_full_partitioned-.txt";

  // Fill names array with appropriate names
  names[ClikPar::omega_b] = "omega_b";
  names[ClikPar::omega_cdm] = "omega_cdm";
  names[ClikPar::hundredxtheta_s] = "hundredxtheta_s";
  names[ClikPar::tau_reio] = "tau_reio";
  names[ClikPar::n_s] = "n_s";
  names[ClikPar::ln10_10_A_s] = "ln10_10_A_s";
  names[ClikPar::A_cib_217] = "A_cib_217";
  names[ClikPar::xi_sz_cib] = "xi_sz_cib";
  names[ClikPar::A_sz] = "A_sz";
  names[ClikPar::ps_A_100_100] = "ps_A_100_100";
  names[ClikPar::ps_A_143_143] = "ps_A_143_143";
  names[ClikPar::ps_A_143_217] = "ps_A_143_217";
  names[ClikPar::ps_A_217_217] = "ps_A_217_217";
  names[ClikPar::ksz_norm] = "ksz_norm";
  names[ClikPar::gal545_A_100] = "gal545_A_100";
  names[ClikPar::gal545_A_143] = "gal545_A_143";
  names[ClikPar::gal545_A_143_217] = "gal545_A_143_217";
  names[ClikPar::gal545_A_217] = "gal545_A_217";
  names[ClikPar::calib_100T] = "calib_100T";
  names[ClikPar::calib_217T] = "calib_217T";
  names[ClikPar::A_planck] = "A_planck";
  names[ClikPar::cib_index] = "cib_index";

  // Fill lower and upper edges for bins (use priors from ClikPar.cc)
  xlow[ClikPar::omega_b] = 0.021; // 3-sigma
  xupp[ClikPar::omega_b] = 0.023;
  xlow[ClikPar::omega_cdm] = 0.108; // 4-sigma
  xupp[ClikPar::omega_cdm] = 0.130;
  xlow[ClikPar::hundredxtheta_s] = 1.0394; // 3-sigma
  xupp[ClikPar::hundredxtheta_s] = 1.0423;
  xlow[ClikPar::tau_reio] = 0.04;  // 2-sigma
  xupp[ClikPar::tau_reio] = 0.116;
  xlow[ClikPar::n_s] = 0.94; // 4-sigma
  xupp[ClikPar::n_s] = 0.99;
  xlow[ClikPar::ln10_10_A_s] = 2.981; // 3-sigma
  xupp[ClikPar::ln10_10_A_s] = 3.197;
  xlow[ClikPar::A_cib_217] = 0.0;
  xupp[ClikPar::A_cib_217] = 200.0;
  xlow[ClikPar::cib_index] = -1.5;
  xupp[ClikPar::cib_index] = -1.1;
  xlow[ClikPar::xi_sz_cib] = 0.0;
  xupp[ClikPar::xi_sz_cib] = 1.0;
  xlow[ClikPar::A_sz] = 0.0;
  xupp[ClikPar::A_sz] = 10.0;
  xlow[ClikPar::ps_A_100_100] = 0.0;
  xupp[ClikPar::ps_A_100_100] = 400.0;
  xlow[ClikPar::ps_A_143_143] = 0.0;
  xupp[ClikPar::ps_A_143_143] = 400.0;
  xlow[ClikPar::ps_A_143_217] = 0.0;
  xupp[ClikPar::ps_A_143_217] = 400.0;
  xlow[ClikPar::ps_A_217_217] = 0.0;
  xupp[ClikPar::ps_A_217_217] = 400.0;
  xlow[ClikPar::ksz_norm] = 0.0;
  xupp[ClikPar::ksz_norm] = 10.0;
  xlow[ClikPar::gal545_A_100] = 0.0;
  xupp[ClikPar::gal545_A_100] = 50.0;
  xlow[ClikPar::gal545_A_143] = 0.0;
  xupp[ClikPar::gal545_A_143] = 50.0;
  xlow[ClikPar::gal545_A_143_217] = 0.0;
  xupp[ClikPar::gal545_A_143_217] = 100.0;
  xlow[ClikPar::gal545_A_217] = 0.0;
  xupp[ClikPar::gal545_A_217] = 400.0;
  xlow[ClikPar::calib_100T] = 0.0;
  xupp[ClikPar::calib_100T] = 3.0;
  xlow[ClikPar::calib_217T] = 0.0;
  xupp[ClikPar::calib_217T] = 3.0;
  xlow[ClikPar::A_planck] = 0.9;
  xupp[ClikPar::A_planck] = 1.1;

  // Make histograms
  // 100 is number of bins
  for (int i = 0; i < nparams; ++i) {
    hists[i] = new TH1D(names[i], names[i], 500, xlow[i], xupp[i]);
  }

  in.open(file_name);
  
  while (1) {
    // Read in parameters
    in >> pmass
       >> like
       >> params[ClikPar::omega_b]
       >> params[ClikPar::omega_cdm]
       >> params[ClikPar::hundredxtheta_s]
       >> params[ClikPar::tau_reio]
       >> params[ClikPar::n_s]
       >> params[ClikPar::ln10_10_A_s]
       >> params[ClikPar::A_cib_217]
       >> params[ClikPar::xi_sz_cib]
       >> params[ClikPar::A_sz]
       >> params[ClikPar::ps_A_100_100]
       >> params[ClikPar::ps_A_143_143]
       >> params[ClikPar::ps_A_143_217]
       >> params[ClikPar::ps_A_217_217]
       >> params[ClikPar::ksz_norm]
       >> params[ClikPar::gal545_A_100]
       >> params[ClikPar::gal545_A_143]
       >> params[ClikPar::gal545_A_143_217]
       >> params[ClikPar::gal545_A_217]
       >> params[ClikPar::calib_100T]
       >> params[ClikPar::calib_217T]
       >> params[ClikPar::A_planck]
       >> params[ClikPar::cib_index];

    if (!in.good()) break;

    // Fill histograms
    for (int i = 0; i < nparams; ++i) {
      hists[i]->Fill(params[i], pmass);
    }
  }

  in.close();

  for (int i = 0; i < nparams; ++i) {
    print_mean(hists[i], names[i]);
    delete hists[i];
  }
}

/*
Output prior to reshuffling everything into arrays:

Processing pc_plot.cc...
omega_b = 0.0210639 +/- 5.77521e-05
omega_cdm = 0.120905 +/- 0.00179987
hundredxtheta_s = 0 +/- 0
tau_reio = 0.0697716 +/- 0.0157874
n_s = 0 +/- 0
ln10_10_A_s = 0 +/- 0
A_cib_217 = 74.2049 +/- 2.66137
xi_sz_cib = 0.0684683 +/- 0.0518798
A_sz = 0 +/- 0
ps_A_100_100 = 0 +/- 0
ps_A_143_143 = 0 +/- 0
ps_A_143_217 = 0 +/- 0
ps_A_217_217 = 0 +/- 0
ksz_norm = 0.0343596 +/- 0.0506788
gal545_A_100 = 0.104305 +/- 0.0517585
gal545_A_143 = 0.069962 +/- 0.0239933
gal545_A_143_217 = 0.196489 +/- 0
gal545_A_217 = 0 +/- 0
calib_100T = 0 +/- 0
calib_217T = 0 +/- 0
A_planck = 0 +/- 0
cib_index = 0 +/- 0
*/
