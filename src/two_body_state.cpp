// Base class to work with four vectors.
//
// All code can be easily modified to use TLorentzVector instead
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "two_body_state.hpp"

double two_body_state::get_mass(string name)
{
  if (name == particle1)
  {
    return m1;
  }
  else if (name == particle2)
  {
    return m2;
  }
  else
  {
    cout << "two_body_state: Unknown particle name passed as argument. Quiting... \n";
    exit(0);
  }
};

// ---------------------------------------------------------------------------
complex<double> two_body_state::energy(double s, string name)
{
  if (name == particle1)
  {
    return (s + m1*m1 - m2*m2) / (2. * sqrt(s));
  }
  else if (name == particle2)
  {
    return (s - m1*m1 + m2*m2) / (2. * sqrt(s));
  }
  else
  {
    cout << "two_body_state: Unknown particle name passed as argument. Quiting... \n";
    exit(0);
  }
};

complex<double> two_body_state::momentum(double s, string name)
{
  if (name == particle1)
  {
    complex<double> E1 = energy(s, name);
    return sqrt(E1*E1 - m1*m1);
  }
  else if (name == particle2)
  {
    complex<double> E2 = energy(s, name);
    return sqrt(E2*E2 - m2*m2);
  }
  else
  {
    cout << "two_body_state: Unknown particle name passed as argument. Quiting... \n";
    exit(0);
  }
};

// ---------------------------------------------------------------------------
// Three momenta components in the x-z plane

complex<double> two_body_state::p0(double s, double zs, string name)
{
  return energy(s, name);
};

complex<double> two_body_state::p1(double s, double zs, string name)
{
  return momentum(s, name) * sqrt(1. - zs * zs);
};

complex<double> two_body_state::p2(double s, double zs, string name)
{
  return 0.;
};

complex<double> two_body_state::p3(double s, double zs, string name)
{
  return momentum(s, name) * zs;
};
