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

class amplitude
{
public:
  // Constructor
  amplitude(reaction_kinematics * xkinem)
  : kinematics(xkinem)
  {};

  // Copy constructor
  amplitude(const amplitude & old)
  : kinematics(old.kinematics)
  {};

  reaction_kinematics * kinematics;

  virtual std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs) = 0;
};

#endif
