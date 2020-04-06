// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/amplitude.hpp"

// ---------------------------------------------------------------------------
// Differential cross section dsigma / dt
double amplitude::diff_xsection(double s, double zs)
{
  // Sum all the helicity amplitudes
  double sum = 0.;
  for (int i = 0; i < 24; i++)
  {
    std::complex<double> square;
    square = helicity_amplitude(kinematics->helicities[i], s, zs);
    square *= conj(helicity_amplitude(kinematics->helicities[i], s, zs));

    sum += real(square);
  }

  double norm = (6084.375 * alpha) / pow(0.5 * (s - mPro_sqr), 2.);

  return sum * norm;
};

// ---------------------------------------------------------------------------
// Polarizatiopn asymmetry between beam and recoil proton
double amplitude::K_LL(double s, double zs)
{
  double sigmapp = 0., sigmapm = 0.;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> squarepp, squarepm;

    // Amplitudes with lam_gam = + and lam_recoil = +
    squarepp = helicity_amplitude(kinematics->helicities[2*i+1], s, zs);
    squarepp *= conj(squarepp);
    sigmapp += real(squarepp);

    // Amplitudes with lam_gam = + and lam_recoil = -
    squarepm = helicity_amplitude(kinematics->helicities[2*i], s, zs);
    squarepm *= conj(squarepm);
    sigmapm += real(squarepm);
  }

  return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}

// ---------------------------------------------------------------------------
// Polarizatiopn asymmetry between beam and target proton
double amplitude::A_LL(double s, double zs)
{
  double sigmapp = 0., sigmapm = 0.;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> squarepp, squarepm;

    // Amplitudes with lam_gam = + and lam_targ = +
    squarepp = helicity_amplitude(kinematics->helicities[i+6], s, zs);
    squarepp *= conj(squarepp);
    sigmapp += real(squarepp);

    // Amplitudes with lam_gam = + and lam_targ = -
    squarepm = helicity_amplitude(kinematics->helicities[i], s, zs);
    squarepm *= conj(squarepm);
    sigmapm += real(squarepm);
  }

  return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}
