// Parameterization of a u-channel background amplitude with couplings from the
// s-channel
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/reggeized_baryon.hpp"

//------------------------------------------------------------------------------
std::complex<double> reggeized_baryon::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  // NOTE WE ASSUME HERE THE HELICITIES ARE IN THE u CHANNEL!
  // if outputting anyting other than unpolarized x-section need to rotate
  // by crossing matrix.

  // To be implimented in the future

  // Sum over all helicities in the t - channel
  return u_channel_amplitude(helicities, s, zs);
};

// ---------------------------------------------------------------------------
// Helicity amplitude in terms of u-channel (unrotated) helicities
std::complex<double> reggeized_baryon::u_channel_amplitude(std::vector<int> helicities, double s, double zs)
{
  // Net helicities
  int lam  = 2 * helicities[0] - helicities[3];
  int lamp = 2 * helicities[2] - helicities[1];
  int M = std::max(std::abs(lam), std::abs(lamp));

  // Center of mass energy and scattering angle
  double u = kinematics->u_man(s, zs);
  std::complex<double> zu = kinematics->z_u(s, zs);

  // Product of residues
  std::complex<double> result = xr;
  result = photo_coupling(lam, u) * hadronic_coupling(lamp, u);

  // Angular momentum barrier factor
  auto pq = [&](double u)
  {
    std::complex<double> q = (u - mPro2) / sqrt(4. * u * xr);
    std::complex<double> p = Kallen(xr * u, xr * kinematics->mVec2, xr * mPro2) / sqrt(4. * xr * u);
    return 4. * p * q;
  };
  if (M == 3)
  {
    result /= pq(u);
  }

  result *= half_angle_factor(lam, lamp, zu);
  result *= regge_propagator(u);
  result *= pow(2. * s, alpha->eval(u) - double(M / 2.));
  return result;
};

//------------------------------------------------------------------------------
// Half angle factors
std::complex<double> reggeized_baryon::half_angle_factor(int lam, int lamp, std::complex<double> zu)
{
  std::complex<double> sinhalf = sqrt((xr - zu) / 2.);
  std::complex<double> coshalf = sqrt((xr + zu) / 2.);

  std::complex<double> result;
  result  = pow(sinhalf, double(std::abs(lam - lamp)) / 2.);
  result *= pow(coshalf, double(std::abs(lam + lamp)) / 2.);

  return result;
}

// ---------------------------------------------------------------------------
// Usual Reggeon Propagator
std::complex<double> reggeized_baryon::regge_propagator(double u)
{
  std::complex<double> alpha_u = alpha->eval(u) - 0.5;

  // the gamma function causes problesm for large t so
  if (std::abs(alpha_u) > 30.)
  {
    return 0.;
  }
  else
  {
    std::complex<double> result;
    result  = - alpha->slope();
    result *= 0.5 * (double(alpha->signature) + exp(-xi * M_PI * alpha_u));
    result *= cgamma(1. - alpha_u);

    return result;
  }
};
