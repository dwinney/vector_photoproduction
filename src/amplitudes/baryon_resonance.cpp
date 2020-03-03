// Parameterization of a resonant amplitude in the s-channel
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/baryon_resonance.hpp"

// Combined amplitude as a Breit-Wigner with the residue as the prodect of hadronic and photo-couplings
std::complex<double> baryon_resonance::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  double residue = photo_coupling(helicities, s) * hadronic_decay(helicities, zs);
  residue *= threshold_factor(s, 1.5);

  return residue / (s - mRes*mRes - xi * mRes * gamRes);
};

// Ad-hoc threshold factor to kill the resonance at threshold
double baryon_resonance::threshold_factor(double s, double beta)
{
  double result = pow((s - sthPsiPro) / s, beta);
  result /= pow((mRes*mRes - sthPsiPro) / (mRes*mRes), beta);

  return result;
};

// Photoexcitation helicity amplitude for the process gamma p -> R
double baryon_resonance::photo_coupling(std::vector<int> helicities, double s)
{
  int lam_i = 2 * helicities[0] - helicities[1];

  int l_min; // lowest allowed relative angular momentum
  double P_t; // Combinatorial factor due to only transverse polarized J/psi contribute
  switch (J)
  {
    case 3:
    {
      if (P == -1)
        {l_min = 0; P_t = 2./3.;}
      else if (P == 1)
        {l_min = 1; P_t = 3./5.;}
      break;
    }
    case 5:
    {
      if (P == -1)
        {l_min = 1; P_t = 3./5.;}
      else if (P == 1)
        {l_min = 2; P_t = 1./3.;}
      break;
    }
    default:
    {
      std::cout << "\n baryon_resonance: spin-parity combination for J = " << J << "/2 and P = " << P << "not available.";
      std::cout << "Quiting... \n";
      exit(0);
    }
  }

  // A_1/2 or A_3/2 depending on ratio R_photo
  double a;
  (lam_i == 1) ? (a = R_photo) : (a = sqrt(1. - R_photo * R_photo));


  // Electromagnetic decay width given by VMD assumption
  double emGamma = (xBR * gamRes) * pow(fJpsi / mJpsi, 2.);
  emGamma *= pow(pi_bar / pf_bar, double(2*l_min +1)) * P_t;

  // Photo-coupling overall size of |A_1/2|^2 + |A_3/2|^2 is restriced from VMD
  double A_lam = emGamma * mPi * mRes * double(J + 1) / (2. * mPro * pi_bar * pi_bar);
  A_lam = sqrt(A_lam) * a;

  double result = sqrt(s) * pi_bar / mRes;
  result *= sqrt(8. * mPro * mRes / real(kinematics->initial.momentum("beam", s)));
  result *= A_lam;

  return result;
};

// Hadronic decay helicity amplitude for the R -> J/psi p process
double baryon_resonance::hadronic_decay(std::vector<int> helicities, double zs)
{
  int lam_i = 2 * helicities[0] - helicities[1];
  int lam_f = 2 * helicities[2] - helicities[3];

  // Hadronic coupling constant g, given in terms of branching ratio xBR
  double g = 8. * M_PI * xBR * gamRes;
  g *= double(J + 1) / 6.;
  g *= mRes * mRes / pf_bar;
  g = sqrt(g);

  // Check for extra phase from unnatural decays
  if (naturality == -1)
  {
    if (helicities[0] < 0) {g *= -1.;}
    if (lam_f < 0) {g *= -1.;}
  }

  return g * wigner_d(J, lam_i, lam_f, zs);
};
