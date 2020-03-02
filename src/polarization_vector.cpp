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
std::complex<double> polarization_vector::component(int i, int lambda, double s, double zs)
{
  if (abs(lambda) == 1)
  {
    switch (i)
    {
      case 0: return 0.;
      case 1: return - double(lambda) * zs / sqrt(2.);
      case 2: return - xi / sqrt(2.);
      case 3: return double(lambda) * sqrt(1. - zs*zs) / sqrt(2.);
    }
  }
  else if (lambda == 0)
  {
    if (abs(mass) < 0.01) return 0.;

    switch (i)
    {
      case 0: return state.momentum(particle, s) / mass;
      case 1: return state.energy(particle, s) * sqrt(1. - zs*zs) / mass;
      case 2: return 0.;
      case 3: return state.energy(particle, s) * zs / mass;
    }
  }
  else
  {
    std::cout << "polarization_vector: Invalid helicity! Quitting... \n";
    exit(0);
  }
};

std::complex<double> polarization_vector::conjugate_component(int i, int lambda, double s, double zs)
{
  return conj(component(i, lambda, s, zs));
};
