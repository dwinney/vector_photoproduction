// ---------------------------------------------------------------------------
// Analytic model for the photoproduction of chi_c1 near threshold at GlueX
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] 10.1103/PhysRevD.96.093008
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/vector_meson_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "utilities.hpp"

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

int main( int argc, char** argv )
{
  double theta = 45.;
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
  }

  // Set up kinematics for the chi_c1
  reaction_kinematics * ptr = new reaction_kinematics(3.510, "chi_c1");

  vector_meson_exchange rho(ptr, .770, "rho");
  rho.set_params({9.2E-4, 2.4, 14.6});
  exchanges.push_back(&rho);

  vector_meson_exchange omega(ptr, .780, "omega");
  omega.set_params({5.2E-4, 16., 0.});
  exchanges.push_back(&omega);

  int N = 50; // how many points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Print contributions from each exchange seperately
double zs = cos(theta * deg2rad);

for (int n = 0; n < exchanges.size(); n++)
{
  if (INTEG == false)
  {
    std::cout << "\n";
    std::cout << "Printing DXS contribution from " << exchanges[n]->identifier;
    std::cout << " exchange at " << theta << " degrees. \n";
  }
  else
  {
    std::cout << "\n";
    std::cout << "Printing integrated cross-section contribution from " << exchanges[n]->identifier;
    std::cout << " exchange. \n";
  }

  std::vector<double> s, dxs;
  for (int i = 0; i < N; i++)
  {
    double si = ptr->sth + double(i) * (100. - ptr->sth) / N;
    double dxsi;

    if (INTEG == false)
    {
      dxsi = exchanges[n]->differential_xsection(si, zs);
    }
    else
    {
      dxsi = exchanges[n]->integrated_xsection(si);
    }

    s.push_back(si);
    dxs.push_back(dxsi);
  }

  quick_print(s, dxs, exchanges[n]->identifier + "_exchange");
  quick_plot(s, dxs, exchanges[n]->identifier + "_exchange");
}
