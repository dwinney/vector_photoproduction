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
  int lam_i = 2 * helicities[0] - helicities[1];
  int lam_f = 2 * helicities[2] - helicities[3];

  std::complex<double> residue = 1.;
  residue  = photo_coupling(lam_i, s);
  residue *= hadronic_coupling(lam_f, s);
  residue *= threshold_factor(s, 1.5);

  residue *= wigner_d_half(J, lam_f, lam_i, zs);
  residue /= (s + xi * mRes * gamRes - mRes*mRes);

  return residue;
};

// Ad-hoc threshold factor to kill the resonance at threshold
double baryon_resonance::threshold_factor(double s, double beta)
{
  double result = pow((s - sthPsiPro) / s, beta);
  result /= pow((mRes*mRes - sthPsiPro) / (mRes*mRes), beta);

  return result;
};

// Photoexcitation helicity amplitude for the process gamma p -> R
std::complex<double> baryon_resonance::photo_coupling(int lam_i, double s)
{
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
      if (P == 1)
        {l_min = 1; P_t = 3./5.;}
      else if (P == -1)
        {l_min = 2; P_t = 1./3.;}
      break;
    }
    default:
    {
      std::cout << "\nbaryon_resonance: spin-parity combination for J = " << J << "/2 and P = " << P << " not available. ";
      std::cout << "Quiting... \n";
      exit(0);
    }
  }

  // A_1/2 or A_3/2 depending on ratio R_photo
  double a;
  (std::abs(lam_i) == 1) ? (a = R_photo) : (a = sqrt(1. - R_photo * R_photo));

  // Electromagnetic decay width given by VMD assumption
  std::complex<double> emGamma = (xBR * gamRes) * pow(fJpsi / mJpsi, 2.);
  emGamma *= pow(xr * pi_bar / pf_bar, double(2 * l_min + 1)) * P_t;

  // Photo-coupling overall size of |A_1/2|^2 + |A_3/2|^2 is restriced from VMD
  std::complex<double> A_lam = emGamma * M_PI * mRes * double(J + 1) / (2. * mPro * pi_bar * pi_bar);
  A_lam = sqrt(xr * A_lam);

  std::complex<double> result = sqrt(xr * s) * pi_bar / mRes;
  result *= sqrt(xr * 8. * mPro * mRes / kinematics->initial.momentum("beam", s));
  result *= A_lam * a;

  // FACTPR PF 4 PI SOMETIMES FACTORED OUT
  result *= sqrt(4. * M_PI * M_ALPHA);

  return result;
};

// Hadronic decay helicity amplitude for the R -> J/psi p process
std::complex<double> baryon_resonance::hadronic_coupling(int lam_f, double s)
{
  // Hadronic coupling constant g, given in terms of branching ratio xBR
  std::complex<double> g;
  g  = 8. * M_PI * xBR * gamRes;
  g *= double(J + 1) / 6.;
  g *= mRes * mRes / kinematics->final.momentum(kinematics->vector_particle, s);
  g = sqrt(xr * g);

  return g;
};
