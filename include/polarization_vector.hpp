// Classes for operations involving Lorentz 4-vectors
// coded up independently to not require ROOT to be installed
//
// Dependencies: None
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _4VEC_
#define _4VEC_

#include <iostream>

#include "constants.hpp"

// ---------------------------------------------------------------------------
// Polarization vectors for vector particles
// in the s-channel center of mass frame
// ---------------------------------------------------------------------------
class polarization_vector
{
private:
    const double mass;
    const bool conj = false; // whether or not this is a complex conjugate vector
    complex<double> Kallen(double x, double y, double z);

public:
  // Constructor
  polarization_vector(double xmass, bool xconj = false)
    : mass(xmass), conj(xconj)
  {};

  // Destructor
  ~polarization_vector(){};

  // vector particle momentum
  complex<double> E(double s);
  complex<double> q(double s);

  // Components
  complex<double> p0(int lambda, double s, double zs);
  complex<double> p1(int lambda, double s, double zs);
  complex<double> p2(int lambda, double s, double zs);
  complex<double> p3(int lambda, double s, double zs);
};

#endif
