// Class for dirac spinors for the spin-1/2 nucleon
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dirac_spinor.hpp"

// ---------------------------------------------------------------------------
// Angular half angle factors
double dirac_spinor::cos_half(double zs)
{
  double result = (1. + zs) / 2.;
  return sqrt(result);
};

double dirac_spinor::sin_half(double zs)
{
  double result = (1. - zs) / 2.;
  return sqrt(result);
};

complex<double> dirac_spinor::momentum(int sign, double s)
{
  if (ANTI_PARTICLE)
  {
    sign *= -1;
  }

  complex<double> E = twobody.energy(particle, s);
  return sqrt(E*E + sign * mass*mass);
}
// ---------------------------------------------------------------------------

complex<double> dirac_spinor::component(int i, int lambda, double s, double zs)
{
  if (lambda == 1)
  {
    switch (i)
    {
      case 0: return momentum(+1, s) * cos_half(zs);
      case 1: return momentum(+1, s) * sin_half(zs);
      case 2: return momentum(-1, s) * cos_half(zs);
      case 3: return momentum(-1, s) * sin_half(zs);
      default : cout << "dirac_spinor: Invalid component index passed as argument. Quitting... \n";
                exit(0);
    }
  }
  else if (lambda == -1)
    {
      switch (i)
      {
        case 0: return - momentum(+1, s) * sin_half(zs);
        case 1: return momentum(+1, s) * cos_half(zs);
        case 2: return momentum(-1, s) * sin_half(zs);
        case 3: return - momentum(-1, s) * cos_half(zs);
        default : cout << "dirac_spinor: Invalid component index passed as argument. Quitting... \n";
                  exit(0);
      }
    }
  else
  {
    cout << "dirac_spinor: Invalid helicity projection passed as argument. Quitting... \n";
    exit(0);
  }
};
