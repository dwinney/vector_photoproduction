// Class for dirac spinors for the spin-1/2 nucleon
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dirac_spinor.hpp"

// ---------------------------------------------------------------------------
// Energy part
std::complex<double> jpacPhoto::dirac_spinor::omega(int sign, double s)
{
  if (ANTI_PARTICLE)
  {
    sign *= -1;
  }

  std::complex<double> E = twobody.energy(particle, s);
  return sqrt(E + sign * mass);
}

// ---------------------------------------------------------------------------
// Angular half angle factors
double jpacPhoto::dirac_spinor::xi(int lam, double zs)
{
  double result = (1. + double(lam) * zs) / 2.;
  return sqrt(result);
};

// ---------------------------------------------------------------------------
// Components for both the regular spinor or adjoint
std::complex<double> jpacPhoto::dirac_spinor::component(int i, int lambda, double s, double zs)
{
  if (abs(lambda) != 1)
  {
    std::cout << "\ndirac_spinor: Invalid helicity projection passed as argument. Quitting... \n";
    exit(0);
  }

  switch (i)
  {
    case 0: return omega(+1, s) * xi(lambda, zs);
    case 1: return double(lambda) * omega(+1, s) * xi(-lambda, zs);
    case 2: return double(lambda) * omega(-1, s) * xi(lambda, zs);
    case 3: return omega(-1, s) * xi(-lambda, zs);
    default : std::cout << "dirac_spinor: Invalid component index " << i << " passed as argument. Quitting... \n";
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
