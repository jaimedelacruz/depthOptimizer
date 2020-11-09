#ifndef __MTYPES_H__
#define __MTYPES_H__

#include <string>
#include <vector>

static const int nlevpf = 6;


/* ---------------------------------------------------------------------------- */

struct species{
  int anum, ion, nLev;
  double n_tot[nlevpf];
  double pf[nlevpf];
  double eion[nlevpf];
};

/* ---------------------------------------------------------------------------- */

struct line{
  char elem[8], label[15];
  double Jup, Jlow, Gup, Glow;
  double w0, nu0, width, eion, nu_min, nu_max;
  double gf, e_low, e_up, amass;
  double g_rad, g_str, g_vdw;
  double b_sig, b_alp, b_vbar, b_gvw;
  int anum, ion,off, idx;
  bool barklem, firsttime;
  
  // Zeeman splitting
  int nZ;
  std::vector<double> strength, splitting;
  std::vector<int> iL;
};
typedef line line_t;


/* /\* ---------------------------------------------------------------------------- *\/ */

/* struct region{ */
/*   double w0, dw, cscal; */
/*   int nw, off; */
/*   std::vector<double> wav, nu; */
/*   std::vector<int> idx; */
/*   std::string inst, ifile; */
/* }; */
/* typedef region region_t; */

/* ---------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------- */




/* ---------------------------------------------------------------------------- */

#endif
