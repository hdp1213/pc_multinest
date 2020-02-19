/***********************************/
/*  INTERNAL GETDIST NAME SETTING  */
/***********************************/

#if PARAM_SET == 4

#ifndef NAME_SET
#define NAME_SET

m_name[pbh_frac] = "pbh_frac";
m_name[pbh_mass] = "pbh_mass";
m_name[A_planck] = "A_planck";
#ifndef LITE_HI_L
m_name[A_cib_217] = "A_cib_217";
m_name[xi_sz_cib] = "xi_sz_cib";
m_name[A_sz] = "A_sz";
m_name[ps_A_100_100] = "ps_A_100_100";
m_name[ps_A_143_143] = "ps_A_143_143";
m_name[ps_A_143_217] = "ps_A_143_217";
m_name[ps_A_217_217] = "ps_A_217_217";
m_name[ksz_norm] = "ksz_norm";
m_name[gal545_A_100] = "gal545_A_100";
m_name[gal545_A_143] = "gal545_A_143";
m_name[gal545_A_143_217] = "gal545_A_143_217";
m_name[gal545_A_217] = "gal545_A_217";
m_name[galf_EE_A_100] = "galf_EE_A_100";
m_name[galf_EE_A_100_143] = "galf_EE_A_100_143";
m_name[galf_EE_A_100_217] = "galf_EE_A_100_217";
m_name[galf_EE_A_143] = "galf_EE_A_143";
m_name[galf_EE_A_143_217] = "galf_EE_A_143_217";
m_name[galf_EE_A_217] = "galf_EE_A_217";
m_name[galf_TE_A_100] = "galf_TE_A_100";
m_name[galf_TE_A_100_143] = "galf_TE_A_100_143";
m_name[galf_TE_A_100_217] = "galf_TE_A_100_217";
m_name[galf_TE_A_143] = "galf_TE_A_143";
m_name[galf_TE_A_143_217] = "galf_TE_A_143_217";
m_name[galf_TE_A_217] = "galf_TE_A_217";
m_name[calib_100T] = "calib_100T";
m_name[calib_217T] = "calib_217T";
#endif
m_name[H0 - FIXED_PARAM_AMT] = "H0*";
m_name[Omega_b - FIXED_PARAM_AMT] = "Omega_b*";
m_name[Omega_cdm - FIXED_PARAM_AMT] = "Omega_cdm*";
m_name[Omega_L - FIXED_PARAM_AMT] = "Omega_L*";
m_name[Omega_g - FIXED_PARAM_AMT] = "Omega_g*";
// m_name[sigma8 - FIXED_PARAM_AMT] = "sigma_8*";
m_name[age - FIXED_PARAM_AMT] = "age*";
m_name[conf_age - FIXED_PARAM_AMT] = "conf_age*";
m_name[z_drag - FIXED_PARAM_AMT] = "z_drag*";
m_name[rs_drag - FIXED_PARAM_AMT] = "rs_drag*";

m_latex[pbh_frac] = "f_\\textrm{PBH}";
m_latex[pbh_mass] = "M_\\textrm{PBH}";
m_latex[A_planck] = "y_{\\rm cal}";
#ifndef LITE_HI_L
m_latex[A_cib_217] = "A^{CIB}_{217}";
m_latex[xi_sz_cib] = "\\xi^{tSZ-CIB}";
m_latex[A_sz] = "A^{tSZ}_{143}";
m_latex[ps_A_100_100] = "A^{PS}_{100}";
m_latex[ps_A_143_143] = "A^{PS}_{143}";
m_latex[ps_A_143_217] = "A^{PS}_{143\\times 217}";
m_latex[ps_A_217_217] = "A^{PS}_{217}";
m_latex[ksz_norm] = "A^{kSZ}";
m_latex[gal545_A_100] = "A^{{\\rm dust}TT}_{100}";
m_latex[gal545_A_143] = "A^{{\\rm dust}TT}_{143}";
m_latex[gal545_A_143_217] = "A^{{\\rm dust}TT}_{143\\times 217}";
m_latex[gal545_A_217] = "A^{{\\rm dust}TT}_{217}";
m_latex[galf_EE_A_100] = "A^{{\\rm dust} EE}_{100}";
m_latex[galf_EE_A_100_143] = "A^{{\\rm dust} EE}_{100\\times 143}";
m_latex[galf_EE_A_100_217] = "A^{{\\rm dust} EE}_{100\\times 217}";
m_latex[galf_EE_A_143] = "A^{{\\rm dust} EE}_{143}";
m_latex[galf_EE_A_143_217] = "A^{{\\rm dust} EE}_{143\\times 217}";
m_latex[galf_EE_A_217] = "A^{{\\rm dust} EE}_{217}";
m_latex[galf_TE_A_100] = "A^{{\\rm dust} TE}_{100}";
m_latex[galf_TE_A_100_143] = "A^{{\\rm dust} TE}_{100\\times 143}";
m_latex[galf_TE_A_100_217] = "A^{{\\rm dust} TE}_{100\\times 217}";
m_latex[galf_TE_A_143] = "A^{{\\rm dust} TE}_{143}";
m_latex[galf_TE_A_143_217] = "A^{{\\rm dust} TE}_{143\\times 217}";
m_latex[galf_TE_A_217] = "A^{{\\rm dust} TE}_{217}";
m_latex[calib_100T] = "c_{100}";
m_latex[calib_217T] = "c_{217}";
#endif
m_latex[H0 - FIXED_PARAM_AMT] = "H_0";
m_latex[Omega_b - FIXED_PARAM_AMT] = "\\Omega_\\textrm{b}";
m_latex[Omega_cdm - FIXED_PARAM_AMT] = "\\Omega_\\textrm{c}";
m_latex[Omega_L - FIXED_PARAM_AMT] = "\\Omega_\\Lambda";
m_latex[Omega_g - FIXED_PARAM_AMT] = "\\Omega_\\gamma";
// m_latex[sigma8 - FIXED_PARAM_AMT] = "\\sigma_8";
m_latex[age - FIXED_PARAM_AMT] = "{\\rm{Age}}/{\\rm{Gyr}}";
m_latex[conf_age - FIXED_PARAM_AMT] = "\\eta_0";
m_latex[z_drag - FIXED_PARAM_AMT] = "z_\\textrm{drag}";
m_latex[rs_drag - FIXED_PARAM_AMT] = "r_\\textrm{drag}";

#endif // NAME_SET

#else
#error incorrect parameter set included
#endif // PARAM_SET == 4
