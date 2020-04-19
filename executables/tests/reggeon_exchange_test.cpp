// ---------------------------------------------------------------------------
// Numerical test code for the new reggeon_exchange amplitude
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/reggeon_exchange.hpp"
#include "utilities.hpp"

#include <cstring>
#include <complex>
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

  // Set up kinematics for the chi_c1
  reaction_kinematics * ptr = new reaction_kinematics(3.510, "chi_c1");

  reggeon_exchange amp(ptr, .780, "omega");
  amp.set_params({5.2E-4, 16., 0.});

  int N = 50; // how many points to plot

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Print contributions from each exchange seperately
  double zs = cos(theta * deg2rad);

  std::cout << "\n";
  std::cout << "Printing DXS contribution from " << amp.identifier;
  std::cout << " exchange at " << theta << " degrees. \n";

  std::vector<double> s, dxs;
  std::vector<std::complex<double>> mom;
  for (int i = 0; i < N; i++)
  {
    double si = (ptr->sth + EPS) + double(i) * (100. - (ptr->sth + EPS)) / N;
    s.push_back(si);

    double dxsi = amp.differentia_xsection(si, zs);;
    dxs.push_back(dxsi);
  }

  quick_print(s, dxs, amp.identifier + "_exchange");
  quick_plot(s, dxs, amp.identifier + "_exchange");
}
