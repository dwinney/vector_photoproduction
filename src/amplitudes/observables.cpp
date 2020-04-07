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
double amplitude::differential_xsection(double s, double zs)
{
  if (s - kinematics->sth < 0.001)
  {
    return 0.;
  }

  // Sum all the helicity amplitudes
  double sum = 0.;
  for (int i = 0; i < 24; i++)
  {
    std::complex<double> square;
    square = helicity_amplitude(kinematics->helicities[i], s, zs);
    square *= conj(helicity_amplitude(kinematics->helicities[i], s, zs));

    sum += real(square);
  }

  double norm = 64. * M_PI * s;
  norm *= real(kinematics->initial.momentum("beam", s)) * real(kinematics->initial.momentum("beam", s));

  return sum / norm;
};

// ---------------------------------------------------------------------------
// Inegrated total cross-section
double amplitude::integrated_xsection(double s)
{
  if (s - kinematics->sth < 0.1)
  {
    return 0.;
  }

  int xN = 15;
  double x[xN+1], w[xN+1];
  NR_gauleg(-1., +1., x, w, xN);

  double sum = 0.;
  for (int i = 0; i <= xN; i++)
  {
    double jacobian; // 2. k * q
    jacobian = 2. * real(kinematics->initial.momentum("beam", s));
    jacobian *= real(kinematics->final.momentum(kinematics->vector_particle, s));

    sum += w[i] * differential_xsection(s, x[i]) * jacobian;
  }

  return sum;
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
