// Classes for operations involving Lorentz 4-vectors
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
// Kinematic functions such that the polarization vector is the only class
// needed to characterize the kinematics of the vector meson
complex<double> polarization_vector::Kallen(double x, double y, double z)
{
    return x*x + y*y + z*z - 2.*(x*y + x*z + z*y);
};

// vector particle momentum
complex<double> polarization_vector::E(double s)
{
  complex<double> result = (s + pow(mPro, 2.) - pow(mass, 2.));
  result /= (2. * sqrt(s));
  return result;
}

complex<double> polarization_vector::q(double s)
{
  complex<double> result = Kallen(s, mPro_sqr, mass*mass);
  result /= 4. * s;

  return sqrt(result);
};

// ---------------------------------------------------------------------------
// Components

complex<double> polarization_vector::p0(int lambda, double s, double zs)
{
  if (abs(mass) > 0.01 && lambda == 0)
  {
    return q(s) / mass;
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
    return (sqrt(s) - E(s)) * sqrt(1. - zs * zs) / mass;
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
    return (sqrt(s) - E(s)) * zs / mass;
  }
};
