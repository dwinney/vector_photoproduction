// Header file for all things gamma matrix related
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "gamma_technology.hpp"

// ---------------------------------------------------------------------------
// Rank two gamma tensor
std::complex<double> jpacPhoto::sigma(int mu, int nu, int i, int j)
{
  std::complex<double> result = 0.;
  for (int k = 0; k < 4; k++)
  {
    result += gamma_matrices[mu][i][k] * gamma_matrices[nu][k][j];
    result -= gamma_matrices[nu][i][k] * gamma_matrices[mu][k][j];
  }

  result *= xr / 2.;

  return result;
};

// ---------------------------------------------------------------------------
// Four dimensional Levi-Civita symbol
double jpacPhoto::levi_civita(int mu, int alpha, int beta, int gamma)
{
  // Error check
  if ((mu > 3 || alpha > 3 || beta > 3 || gamma > 3) || (mu < 0 || alpha < 0 || beta < 0 || gamma < 0))
  {
    std::cout << " \nLeviCivita: Error! Invalid argument recieved. Quitting... \n";
    exit(0);
  }

  // Return 0 if any are equal
  if (mu == alpha || mu == beta || mu == gamma || alpha == beta || alpha == gamma || beta == gamma)
  {
    return 0.;
  }

  // Else compare with strings
  std::string input = std::to_string(mu) + std::to_string(alpha) + std::to_string(beta) + std::to_string(gamma);
  std::vector<std::string>  even_permutations =
  {
    "0123",
    "0231",
    "0312",
    "1032",
    "1203",
    "1320",
    "2013",
    "2130",
    "2301",
    "3021",
    "3102",
    "3210"
  };

  for (int i = 0; i < 12; i++)
  {
    if (input == even_permutations[i])
    {
      return 1.;
    }
  }

  std::vector<std::string> odd_permutations =
  {
    "0132",
    "0213",
    "0321",
    "1023",
    "2103",
    "3120",
    "1230",
    "1302",
    "2031",
    "2310",
    "3012",
    "3201"
  };

  for (int i = 0; i < 12; i++)
  {
    if (input == odd_permutations[i])
    {
      return -1.;
    }
  }

  return 0.;
};
