// Class for dirac spinors for the spin-1/2 nucleon
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dirac_spinor.hpp"

using namespace jpacPhoto;

// ---------------------------------------------------------------------------
// Angular half angle factors
double jpacPhoto::dirac_spinor::cos_half(double zs)
{
  double result = (1. + zs) / 2.;
  return sqrt(result);
};

double jpacPhoto::dirac_spinor::sin_half(double zs)
{
  double result = (1. - zs) / 2.;
  return sqrt(result);
};

std::complex<double> jpacPhoto::dirac_spinor::momentum(int sign, double s)
{
  if (ANTI_PARTICLE)
  {
    sign *= -1;
  }

  std::complex<double> E = twobody.energy(particle, s);
  return sqrt(E + sign * mass);
}
// ---------------------------------------------------------------------------
// Components for both the regular spinor or adjoint
std::complex<double> jpacPhoto::dirac_spinor::component(int i, int lambda, double s, double zs)
{
  if (lambda == 1)
  {
    switch (i)
    {
      case 0: return - momentum(+1, s) * sin_half(zs);
      case 1: return momentum(+1, s) * cos_half(zs);
      case 2: return - momentum(-1, s) * sin_half(zs);
      case 3: return momentum(-1, s) * cos_half(zs);
      default : std::cout << "dirac_spinor: Invalid component index " << i << " passed as argument. Quitting... \n";
                exit(0);
    }
  }
  else if (lambda == -1)
    {
      switch (i)
      {
        case 0: return momentum(+1, s) * cos_half(zs);
        case 1: return momentum(+1, s) * sin_half(zs);
        case 2: return - momentum(-1, s) * cos_half(zs);
        case 3: return - momentum(-1, s) * sin_half(zs);
        default : std::cout << "dirac_spinor: Invalid component index " << i << " passed as argument. Quitting... \n";
                  exit(0);
      }
    }
  else
  {
    std::cout << "dirac_spinor: Invalid helicity projection passed as argument. Quitting... \n";
    exit(0);
  }
};

std::complex<double> jpacPhoto::dirac_spinor::adjoint_component(int i, int lambda, double s, double zs)
{
  std::complex<double> result = 0.;
  for (int j = 0; j < 4; j++)
  {
    result += conj(component(j, lambda, s, zs)) * gamma_matrices[0][j][i];
  }
  return result;
};
