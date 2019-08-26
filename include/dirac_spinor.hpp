// Class for dirac spinors for the spin-1/2 nucleon
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _SPINOR_
#define _SPINOR_

#include <cmath>
#include <complex>

#include "two_body_state.hpp"

using std::complex;

class dirac_spinor
{
private:
  two_body_state twobody;
  const string particle;
  const double mass;
  
  const bool ANTI_PARTICLE = false;

  double cos_half(double zs);
  double sin_half(double zs);

  complex<double> momentum(int sign, double s);

public:
  dirac_spinor(two_body_state xstate, string name, bool if_anti = false)
  : twobody(xstate), particle(name), mass(xstate.get_mass(name)),
    ANTI_PARTICLE(if_anti)
  {};

  complex<double> component(int i, int lambda, double s, double zs);
};
#endif
