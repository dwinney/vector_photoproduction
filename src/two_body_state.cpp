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
complex<double> two_body_state::energy(string name, double s)
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

complex<double> two_body_state::momentum(string name, double s)
{
  if (name == particle1)
  {
    complex<double> E1 = energy(name, s);
    return sqrt(E1*E1 - m1*m1);
  }
  else if (name == particle2)
  {
    complex<double> E2 = energy(name, s);
    return -sqrt(E2*E2 - m2*m2);
  }
  else
  {
    cout << "two_body_state: Unknown particle name passed as argument. Quiting... \n";
    exit(0);
  }
};

// ---------------------------------------------------------------------------
// Three momenta components in the x-z plane
complex<double> two_body_state::component(int i, string name, double s, double zs)
{
  switch (i)
  {
    case 0: return energy(name, s);
    case 1: momentum(name, s) * sqrt(1. - zs * zs);
    case 2: return 0.;
    case 3: return momentum(name, s) * zs;
    default: cout << "two_body_state: Invalid four vector component! Quiting...\n";
             exit(0);
  }
};
