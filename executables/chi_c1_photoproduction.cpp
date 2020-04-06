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
#include "regge_trajectory.hpp"
#include "utilities.hpp"

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

int main( int argc, char** argv )
{
  double theta = 0.;
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
  }

  // Set up kinematics for the chi_ci
  reaction_kinematics * ptr = new reaction_kinematics(3.510, "chi_c1");

  // Assume it proceeds entirely by omega exchange
  vector_meson_exchange omega(ptr, .780, "omega");

  // Set the three couplings
  std::vector<double> omega_couplings = {1., 16, 1.};
  omega.set_params(omega_couplings);

  int N = 100; // how many points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Print .dat files of differential cross section and plot it to a .pdf

double zs = cos(theta * deg2rad);

std::cout << "\n";
std::cout << "Printing DXS for " << ptr->vector_particle;
std::cout << " production at " << theta << " degrees in center-of-mass frame.";
std::cout << "\n";

std::vector<double> s, dxs;
for (int i = 0; i < N; i++)
{
  double si = ptr->sth + double(i) * (100. - ptr->sth) / N;
  double dxsi = omega.diff_xsection(si, zs);

  s.push_back(si);
  dxs.push_back(dxsi);
}

quick_print(s, dxs, ptr->vector_particle + "_photoproduction");
quick_plot(s, dxs, ptr->vector_particle + "_photoproduction");

std::cout << "\n";

return 1.;
};
