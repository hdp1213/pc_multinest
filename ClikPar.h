#ifndef CLIKPAR_H
#define CLIKPAR_H

struct ClikPar
{
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
    cib_index // cib_index = -1.3 is constant
  };
};

#endif