// Axial-vector meson photoproduction proceeding through a Reggeized meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/reggeized_meson.hpp"

// ---------------------------------------------------------------------------
// Assemble the helicity amplitude by contracting the lorentz indices
std::complex<double> reggeized_meson::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  // NOTE WE ASSUME HERE THE HELICITIES ARE IN THE t CHANNEL!
  // if outputting anyting other than unpolarized x-section need to rotate
  // by crossing matrix.

  // To be implimented in the future

  // Sum over all helicities in the t - channel
  return t_channel_amplitude(helicities, s, zs);
};

// ---------------------------------------------------------------------------
// Helicity amplitude in terms of t-channel (unrotated) helicities
std::complex<double> reggeized_meson::t_channel_amplitude(std::vector<int> helicities, double s, double zs)
{
  // Net helicities
  int lam  = helicities[0] - helicities[2];
  int lamp = (helicities[1] - helicities[3]) / 2;
  int M = std::max(std::abs(lam), std::abs(lamp));

  // Double flip is assumed to be negligable
  if (M == 2)
  {
    return 0.;
  }

  // Center of mass energy and scattering angle
  double t = kinematics->t_man(s, zs);
  std::complex<double> z_t = kinematics->z_t(s,zs);

  // Product of residues
  std::complex<double> result = xr;
  result = top_residue(lam, t) * bottom_residue(lamp, t);

  // Angular momentum barrier factor
  auto pq = [&](double t)
  {
    std::complex<double> q = (t - kinematics->mVec2) / sqrt(4. * t * xr);
    std::complex<double> p = sqrt(xr * t - 4.*mPro2) / 2.;
    return 2. * p * q;
  };

  // Theres only a nontrivial barrier factor in the lam = lamp = 0 case
  if (M == 0)
  {
    result /= pq(t);
  }

  result *= half_angle_factor(lam, lamp, z_t);
  result *= regge_propagator(t);
  result *= pow(s, alpha->eval(t) - double(M));
  return result;
};

//------------------------------------------------------------------------------
// Half angle factors
std::complex<double> reggeized_meson::half_angle_factor(int lam, int lamp, std::complex<double> z_t)
{
  std::complex<double> sinhalf = sqrt((xr - z_t) / 2.);
  std::complex<double> coshalf = sqrt((xr + z_t) / 2.);

  std::complex<double> result;
  result  = pow(sinhalf, double(std::abs(lam - lamp)) / 2.);
  result *= pow(coshalf, double(std::abs(lam + lamp)) / 2.);

  return result;
}

// ---------------------------------------------------------------------------
// Usual Reggeon Propagator
std::complex<double> reggeized_meson::regge_propagator(double t)
{
  std::complex<double> alpha_t = alpha->eval(t);

  // the gamma function causes problesm for large t so
  if (std::abs(alpha_t) > 30.)
  {
    return 0.;
  }
  else
  {
    std::complex<double> result;
    result  = - alpha->slope();
    result *= 0.5 * (double(alpha->signature) + exp(-xi * M_PI * alpha_t));
    result *= cgamma(1. - alpha_t);

    return result;
  }
};
