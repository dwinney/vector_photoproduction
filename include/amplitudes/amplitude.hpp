// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AMPLITUDE_
#define _AMPLITUDE_

#include "reaction_kinematics.hpp"

class amplitude
{
public:
  // Constructor
  amplitude(reaction_kinematics * xkinem)
  : kinematics(xkinem)
  {};

  reaction_kinematics * kinematics;

  virtual std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs) = 0;
};

#endif
