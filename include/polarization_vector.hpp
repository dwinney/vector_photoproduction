// Class for the polarization vector of vector particles
// coded up independently to not require ROOT to be installed
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _4VEC_
#define _4VEC_

#include <iostream>
#include "constants.hpp"
#include "two_body_state.hpp"

// ---------------------------------------------------------------------------
// Polarization vectors for vector particles
// in the s-channel center of mass frame
// ---------------------------------------------------------------------------
class polarization_vector
{
private:
    two_body_state state;
    const string particle;
    const double mass;

    const bool conj = false; // whether or not this is a complex conjugate vector

public:
  // Constructor
  polarization_vector(two_body_state xstate, string name, bool xconj = false)
    : state(xstate), particle(name), mass(state.get_mass(name)),
      conj(xconj)
  {};

  // Destructor
  ~polarization_vector(){};

  // Components
  complex<double> component(int i, int lambda, double s, double zs);
};

#endif
