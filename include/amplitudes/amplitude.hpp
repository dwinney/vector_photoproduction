// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AMPLITUDE_
#define _AMPLITUDE_

// ---------------------------------------------------------------------------
// Abstract class to define helicity amplitudes. This will allow multiple different
// classes (for s, t, and u- channels but also multiple contibutions in each channel)
// to be added together and evaluated in observables.
//
// Any generic amplitude needs a reaction_kinematics object
// and a way to evaluate the helicity amplitude for given set of helicities,
// CoM energy and scattering angle.
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"

#include <string>

class amplitude
{
public:
  // Constructor
  amplitude(reaction_kinematics * xkinem)
  : kinematics(xkinem)
  {};

  amplitude(reaction_kinematics * xkinem, std::string id)
  : kinematics(xkinem), identifier(id)
  {};

  // Copy constructor
  amplitude(const amplitude & old)
  : kinematics(old.kinematics), identifier(old.identifier)
  {};

  // Kinematics object for thresholds and etc.
  reaction_kinematics * kinematics;

  // Some saveable string by which to identify the amplitude
  std::string identifier;

  virtual std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs) = 0;

  // ---------------------------------------------------------------------------
  // Observables

  // Differential and total cross-section
  double diff_xsection(double s, double zs);

  // Spin asymmetries
  double K_LL(double s, double zs);
  double A_LL(double s, double zs);
};

#endif
