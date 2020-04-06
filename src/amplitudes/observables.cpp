// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/amplitude.hpp"

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
