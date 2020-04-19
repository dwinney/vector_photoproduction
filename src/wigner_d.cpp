// Wigner little-d functions for half-integer spin particles
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "misc_math.hpp"

void wigner_error(int j, int lam1, int lam2)
{
  std::cout << "\n";
  std::cout << "wigner_d: Argument combination with j = " << j << "/2";
  std::cout << ", lam1 = " << lam1 << "/2";
  std::cout << ", lam2 = " << lam2 << "/2";
  std::cout << " does not exist. Quitting... \n";
  exit(0);
}

std::complex<double> wigner_d_half(int j, int lam1, int lam2, std::complex<double> z)
{
  // if (abs(z) > 1.)
  // {
  //   std::cout << "wigner_d: Angular argument outside of real region (-1 < z < 1).";
  //   std::cout << "Quitting... \n";
  //   exit(0);
  // }

  double phase = 1.;
  // If first lam argument is smaller, switch them
  if (abs(lam1) < abs(lam2))
  {
    int temp = lam1;
    lam1 = lam2;
    lam2 = temp;

    phase *= pow(-1., double(lam1 - lam2)/2.);
  };

  // If first lam is negative, smitch them
  if (lam1 < 0)
  {
    lam1 *= -1;
    lam2 *= -1;

    phase *= pow(-1., double(lam1 - lam2)/2.);
  }

  std::complex<double> result = 0.;
  switch (j)
  {
    // Spin-1/2
    case 1:
    {
      if (lam1 == 1)
      {
        if (lam2 == 1)
        {
          result = sqrt((xr + z) / 2.); break;
        }
        else
        {
          result = sqrt((xr - z) / 2.); break;
        }
      }
      else
      {
        wigner_error(j, lam1, lam2);
      }
      break;
    }

    // Spin-3/2
    case 3:
    {
      // Lambda = 3/2
      if (lam1 == 3)
      {
        switch (lam2)
        {
          case 3: result = pow((xr + z) / 2., 1.5); break;
          case 1:
          {
            result = - sqrt(3.) * (xr + z) / 2.;
            result *= sqrt((xr - z) / 2.); break;
          }
          case -1:
          {
            result = sqrt(3.) * (xr - z) / 2.;
            result *= sqrt((xr + z) / 2.); break;
          }
          case -3: result = - pow((xr - z) / 2., 1.5); break;
          default: wigner_error(j, lam1, lam2);
        }
      }

      // Lambda = 1/2
      else if (lam1 == 1)
      {
        switch (lam2)
        {
          case 1:
          {
            result = (3. * z - 1.) / 2.;
            result *= sqrt((xr + z) / 2.); break;
          }
          case -1:
          {
            result = - (3. * z + 1.) / 2.;
            result *= sqrt((xr - z) / 2.); break;
          }
          default: wigner_error(j, lam1, lam2);
        }
      }

      // Error
      else
      {
        wigner_error(j, lam1, lam2);
      }
      break;
    }

    // Spin-5/2
    case 5:
    {
      switch (lam1)
      {
        case 5: wigner_error(j, lam1, lam2);
        case 3:
        {
          switch (lam2)
          {
            case 3:
            {
              result = (-3. + 5. * z) * (1. + z) / 4.;
              result *= sqrt((xr + z)); break;
            }
            case 1:
            {
              result = (1. + z) * (-1. + 5.*z) / 4.;
              result *= - sqrt(xr - z); break;
            }
            case -1:
            {
              result = (-1. + z) * (1. + 5.*z) / 4.;
              result *= - sqrt(xr + z); break;
            }
            case -3:
            {
              result = - pow((xr - z)/ 2., 1.5) * (3. + 5.*z); break;
            }
            default: wigner_error(j, lam1, lam2);
          }
          break;
        }
        case 1:
        {
          switch (lam2)
          {
            case 1:
            {
              result = sqrt((xr + z) / 2.);
              result *= (-1. - 2.*z + 5.*z*z) / 2.; break;
            }
            case -1:
            {
              result = sqrt((xr - z) / 2.);
              result *= (1. - 2.*z - 5.*z*z) / 2.; break;
            }
            default: wigner_error(j, lam1, lam2);
          }
          break;
        }
        default: wigner_error(j, lam1, lam2);
      }
      break;
    }
    // Error
    default: wigner_error(j, lam1, lam2);
  }

  return phase * result;
};
