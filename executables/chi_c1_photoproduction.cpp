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
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "utilities.hpp"

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

int main( int argc, char** argv )
{
  double theta = 45.;
  bool INTEG = false;
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-integ")==0) INTEG = true;
  }

  // Set up kinematics for the chi_c1
  reaction_kinematics * ptr = new reaction_kinematics(3.510, "chi_c1");

  // Reaction proceeds through multiple vector exchanges.
  // Which we will sum incoherently
  std::vector<amplitude*> exchanges;

  vector_exchange rho(ptr, .770, "rho");
  rho.set_params({9.2E-4, 2.4, 14.6});
  exchanges.push_back(&rho);

  vector_exchange omega(ptr, .780, "omega");
  omega.set_params({5.2E-4, 16., 0.});
  exchanges.push_back(&omega);

  vector_exchange phi(ptr, 1.10, "phi");
  phi.set_params({4.2E-4, -6.2, 2.1});
  exchanges.push_back(&phi);

  vector_exchange jpsi(ptr, 3.097, "jpsi");
  jpsi.set_params({1., 3.3E-3, 0.});
  exchanges.push_back(&jpsi);

  // The total amplitude with all the above exchanges
  amplitude_sum total(ptr, exchanges);

  int N = 20; // how many points to plot

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

// ---------------------------------------------------------------------------
// Print the total cross-section
if (INTEG == false)
{
  std::cout << "\n";
  std::cout << "Printing total DXS for " << ptr->vector_particle;
  std::cout << " photoproduction at " << theta << " degrees. \n";
}
else
{
  std::cout << "\n";
  std::cout << "Printing total integrated cross-section for " << ptr->vector_particle;
  std::cout << " photoproduction. \n";
}

std::vector<double> s, dxs;
for (int i = 0; i < N; i++)
{
  double si = ptr->sth + double(i) * (100. - ptr->sth) / N;

  double dxsi;
  if (INTEG == false)
  {
    dxsi = total.differential_xsection(si, zs);
  }
  else
  {
    dxsi = total.integrated_xsection(si);
  }

  s.push_back(si);
  dxs.push_back(dxsi);
}

quick_print(s, dxs, ptr->vector_particle + "_photoproduction");
quick_plot(s, dxs, ptr->vector_particle + "_photoproduction");
std::cout << "\n";

delete ptr;
return 1.;
};
