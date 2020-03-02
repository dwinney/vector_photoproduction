// Parameterization of a resonant amplitude in the s-channel
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/baryon_resonance.hpp"

std::complex<double> baryon_resonance::hadronic_decay(std::vector<int> helicities, double zs)
{
  int lam_i = 2 * helicities[0] - helicities[1];
  int lam_f = 2 * helicities[2] - helicities[3];

  // Final momentum evaluated at the res
  std::complex<double> g = 8. * M_PI * xBR * gamRes;
  g *= (double(J) + 1.) / 6.;
  g *= mRes * mRes / pf_bar;
  g = sqrt(g);

  return g * wigner_d(J, lam_i, lam_f, zs);
}
