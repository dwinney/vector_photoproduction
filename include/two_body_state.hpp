// Base class to work with four vectors in a scattering event
//
// All code can be easily modified to use TLorentzVector instead
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TWO_BODY_
#define _TWO_BODY_

#include <string>
#include <complex>
#include <iostream>

using std::cout;
using std::complex;
using std::string;

// ---------------------------------------------------------------------------
// The two_body_state is the base object for defining a reaction in the
// s-channel center of mass scatering frame.
//
// Two particles of mass m1 and m2 are defined with momenta opposite along the
// same axis such that the energy and momenta of both particles is entirely
// determined by the center-of-mass energy, s, and the cosing of the angle
// from the z-axis (define to be at theta = 0).
// ---------------------------------------------------------------------------

class two_body_state
{
private:
  const double m1, m2;
  const string particle1, particle2;

public:
  two_body_state(double xm1, double xm2, string name1, string name2)
  : m1(xm1), m2(xm2), particle1(name1), particle2(name2)
  {};

  double get_mass(string name);

  complex<double> momentum(string name, double s);
  complex<double> energy(string name, double s);

  complex<double> component(int i, string name, double s, double zs);
};

#endif
