// Class for the polarization vector of vector particles
// coded up independently to not require ROOT to be installed
//
// Dependencies: None
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "polarization_vector.hpp"

// ---------------------------------------------------------------------------
// Components
complex<double> polarization_vector::p0(int lambda, double s, double zs)
{
  if (abs(mass) > 0.01 && lambda == 0)
  {
    return state.momentum(s, particle) / mass;
  }
  else
  {
    return 0.;
  }
};

complex<double> polarization_vector::p1(int lambda, double s, double zs)
{
  if (abs(lambda) > 1)
  {
    std::cout << "polarization_vector: Invalid helicity! Quitting... \n";
    exit(0);
  }

  if (abs(mass) < 0.01 && lambda == 0)
  {
    return 0.;
  }
  else if (abs(lambda) == 1)
  {
    return - double(lambda) * zs / sqrt(2.);
  }
  else
  {
    return state.energy(s, particle) / mass;
  }
};

complex<double> polarization_vector::p2(int lambda, double s, double zs)
{
  if (abs(lambda) > 1)
  {
    std::cout << "polarization_vector: Invalid helicity! Quitting... \n";
    exit(0);
  }
  if (lambda == 0)
  {
    return 0.;
  }
  else if (conj == true)
  {
    return xi / sqrt(2.);
  }
  else
  {
    return - xi / sqrt(2.);
  }
};

complex<double> polarization_vector::p3(int lambda, double s, double zs)
{
  if (abs(lambda) > 1)
  {
    std::cout << "polarization_vector: Invalid helicity! Quitting... \n";
    exit(0);
  }

  if (abs(mass) < 0.01 && lambda == 0)
  {
    return 0.;
  }
  else if (abs(lambda) == 1)
  {
    return double(lambda) * sqrt(1. - zs * zs) / sqrt(2.);
  }
  else
  {
    return state.energy(s, particle) * zs / mass;
  }
};
