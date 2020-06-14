// Wigner little-d functions for half-integer spin particles
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "misc_math.hpp"

// ---------------------------------------------------------------------------
double jpacPhoto::wigner_error(int j, int lam1, int lam2, bool half)
{
  // if (half == true)
  // {
  //   std::cout << "\n";
  //   std::cout << "wigner_d: Argument combination with j = " << j << "/2";
  //   std::cout << ", lam1 = " << lam1 << "/2";
  //   std::cout << ", lam2 = " << lam2 << "/2";
  //   std::cout << " does not exist. Quitting... \n";
  // }
  // else
  // {
  //   std::cout << "\n";
  //   std::cout << "wigner_d: Argument combination with j = " << j;
  //   std::cout << ", lam1 = " << lam1;
  //   std::cout << ", lam2 = " << lam2;
  //   std::cout << " does not exist. Quitting... \n";
  // }
  //
  // exit(0);

  return 0.;
}
// ---------------------------------------------------------------------------
std::complex<double> jpacPhoto::wigner_d_half(int j, int lam1, int lam2, std::complex<double> z)
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
        wigner_error(j, lam1, lam2, true);
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
          case 3:
          {
             result = pow((xr + z) / 2., 1.5);
             break;
          }
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
          case -3:
          {
            result = - pow((xr - z) / 2., 1.5);
            break;
          }
          default: wigner_error(j, lam1, lam2, true);
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
            result *= sqrt((xr + z) / 2.);
            break;
          }
          case -1:
          {
            result = - (3. * z + 1.) / 2.;
            result *= sqrt((xr - z) / 2.);
            break;
          }
          default: wigner_error(j, lam1, lam2, true);
        }
      }

      // Error
      else
      {
        wigner_error(j, lam1, lam2, true);
      }
      break;
    }

    // Spin-5/2
    case 5:
    {
      switch (lam1)
      {
        case 5:
        {
          wigner_error(j, lam1, lam2, true);
          break;
        }
        // lam1 == 3
        case 3:
        {
          switch (lam2)
          {
            case 3:
            {
              result  = (-3. + 5. * z) * (1. + z) / 4.;
              result *= sqrt((xr + z) / 2.);
              break;
            }
            case 1:
            {
              result = (1. + z) * (-1. + 5.*z) / 4.;
              result *= - sqrt(xr - z);
              break;
            }
            case -1:
            {
              result = (-1. + z) * (1. + 5.*z) / 4.;
              result *= - sqrt(xr + z);
              break;
            }
            case -3:
            {
              result  = -(3. + 5. * z) * (1. - z) / 4.;
              result *= sqrt((xr - z) / 2.);
              break;
            }
            default: wigner_error(j, lam1, lam2, true);
          }
          break;
        }
        // lam1 == 1
        case 1:
        {
          switch (lam2)
          {
            case 1:
            {
              result  = sqrt((xr + z) / 2.);
              result *= (-1. - 2.*z + 5.*z*z) / 2.;
              break;
            }
            case -1:
            {
              result  = sqrt((xr - z) / 2.);
              result *= (1. - 2.*z - 5.*z*z) / 2.;
              break;
            }
            default: wigner_error(j, lam1, lam2, true);
          }
          break;
        }
        default: wigner_error(j, lam1, lam2, true);
      }
      break;
    }
    // Error
    default: wigner_error(j, lam1, lam2, true);
  }

  return phase * result;
};

// ---------------------------------------------------------------------------
std::complex<double> jpacPhoto::wigner_d_int(int j, int lam1, int lam2, std::complex<double> z)
{

  double phase = 1.;
  // If first lam argument is smaller, switch them
  if (abs(lam1) < abs(lam2))
  {
    int temp = lam1;
    lam1 = lam2;
    lam2 = temp;

    phase *= pow(-1., double(lam1 - lam2));
  };

  // If first lam is negative, smitch them
  if (lam1 < 0)
  {
    lam1 *= -1;
    lam2 *= -1;

    phase *= pow(-1., double(lam1 - lam2));
  }

  std::complex<double> result = 0.;
  switch (j)
  {
    // spin - 1
    case 1:
    {
      if (lam1 == 1)
      {
        switch (lam2)
        {
          case 1:
          {
            result = (1. + z) / 2.;
            break;
          }
          case 0:
          {
            result = - sqrt(xr - z*z) / sqrt(2.);
            break;
          }
          case -1:
          {
            result = (1. - z) / 2.;
            break;
          }
          default: wigner_error(j, lam1, lam2, false);
        }
      }
      else if (lam1 == 0)
      {
        result = z;
      }
      else
      {
        wigner_error(j, lam1, lam2, false);
      }
      break;
    }
  // Error
  default: wigner_error(j, lam1, lam2, false);
}

return phase * result;
};
