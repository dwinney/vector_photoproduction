// ---------------------------------------------------------------------------
// Comparing the photoproduction cross-sections of the Jpsi 1s and 2s states
// near threshold at GlueX
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] 10.1103/PhysRevD.100.034019
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/pomeron_exchange.hpp"
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

  // Set up kinematics, determined entirely by vector meson mass
  reaction_kinematics * ptr1s = new reaction_kinematics(3.097, "PSI(1s)");
  reaction_kinematics * ptr2s = new reaction_kinematics(3.686, "PSI(2s)");

  // Set up pomeron trajectory
  // Here we use (real) linear trajectory with intercept and slope only free params
  // Best fit values from [1]
  linear_traj alpha(0.941, 0.364);

  // Create amplitudes with kinematics and same trajectory
  pomeron_exchange pomeron_1s(ptr1s, &alpha);
  pomeron_exchange pomeron_2s(ptr2s, &alpha);

  // Feed in other two parameters (normalization and t-slope)
  // Best fit from [1]
  std::vector<double> params_1s = {0.379, 0.12};
  pomeron_1s.set_params(params_1s);

  // Same t-slope but coupling is scaled by 1/4
  std::vector<double> params_2s = {0.379 / 4., 0.12};
  pomeron_2s.set_params(params_2s);

  int N = 100; // how many points to plot

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Print .dat files of differential cross section and plot it to a .pdf

  std::vector<amplitude*> amp = {&pomeron_1s, &pomeron_2s};
  double zs = cos(theta * deg2rad);

  for (int n = 0; n < amp.size(); n++)
  {
    std::cout << "\n";
    std::cout << "Printing DXS for " << amp[n]->kinematics->vector_particle;
    std::cout << " production at " << theta << " degrees in center-of-mass frame.";
    std::cout << "\n";

    std::vector<double> s_n, dxs_n;

    s_n.push_back((amp[n]->kinematics->sth/mPro - mPro)/2.);
    dxs_n.push_back(0.);

    for (int i = 2; i < N; i++)
    {
      double si = amp[n]->kinematics->sth + double(i) * (30. - amp[n]->kinematics->sth) / N;
      double dxsi = amp[n]->diff_xsection(si, zs);

      s_n.push_back((si/mPro - mPro)/2.);

      dxs_n.push_back(dxsi);
    }

    quick_print(s_n, dxs_n, amp[n]->kinematics->vector_particle + "_photoproduction");
    quick_plot(s_n, dxs_n, amp[n]->kinematics->vector_particle + "_photoproduction");
  }
  std::cout << "\n";

  // ---------------------------------------------------------------------------
  // Do the same for the ratio though
  std::cout << "Printing ratio of cross-sections. \n";

  std::vector<double> s, ratio;
  s.push_back((ptr2s->sth/mPro - mPro) / 2.);
  ratio.push_back(0.);

  for (int i = 10; i < N; i++)
  {
    double si = ptr2s->sth + double(i) * (30 - ptr2s->sth) / N;
    double ratioi = pomeron_2s.diff_xsection(si, zs) / pomeron_1s.diff_xsection(si, zs);

    s.push_back((si/mPro - mPro)/2.);
    ratio.push_back(ratioi);
  }

  quick_print(s, ratio, "dxs_ratio");
  quick_plot(s, ratio, "dxs_ratio");
  std::cout << "\n";

  delete ptr1s, ptr2s;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
  return 1.;
};
