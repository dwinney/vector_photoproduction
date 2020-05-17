// ---------------------------------------------------------------------------
// Predicted sensativity to LHCb pentaquarks in double polarization Observables
// at Hall A at JLab.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] 10.1103/PhysRevD.100.034019
// [2] 10.1103/PhysRevLett.115.072001
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "utilities.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/baryon_resonance.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

int main( int argc, char** argv )
{
  double theta = 0.;
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
  }

  // Set up Kinematics
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");

  // Incoherent sum of the s and t channels
  amplitude_sum sum(ptr);

  // ---------------------------------------------------------------------------
  // S - CHANNEL

  // Two different pentaquarks
  // masses and widths from 2015 LHCb paper [2]
  baryon_resonance P_c1(3, -1, 4.45, 0.4, ptr);
  baryon_resonance P_c2(5, +1, 4.38, 0.1, ptr);

  // 2% branching fraction and equal photocouplings for both
  std::vector<double> params = {0.02, .7071};
  P_c1.set_params(params);
  P_c2.set_params(params);

  // Add them to the sum
  sum.add_amplitude(&P_c1);
  sum.add_amplitude(&P_c2);

  // ---------------------------------------------------------------------------
  // T - CHANNEL

  // Set up pomeron trajectory
  linear_trajectory alpha(1., 0.941, 0.364, "pomeron");

  // Create amplitude with kinematics and trajectory
  pomeron_exchange background(ptr, &alpha);

  // normalization and t-slope
  std::vector<double> back_params = {0.379, 0.12};
  background.set_params(back_params);

  // Add to the sum
  sum.add_amplitude(&background);

  int N = 100; // how many points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Print .dat files of differential cross section and plot it to a .pdf

double zs = cos(theta * deg2rad);
std::vector<double> s, kll, all;
for (int i = 0; i < N; i++)
{
  double si = (ptr->sth + 10.*EPS) + double(i) * (30. - (ptr->sth + 10.*EPS)) / N;
  s.push_back((si/mPro - mPro)/2.);

  double klli = sum.K_LL(si, zs);
  kll.push_back(klli);

  double alli = sum.A_LL(si, zs);
  all.push_back(alli);
}

// Print results
std::cout << "\n";
quick_print(s, kll, "pentaquark_KLL");
quick_plot(s, kll, "pentaquark_KLL");

quick_print(s, all, "pentaquark_ALL");
quick_plot(s, all, "pentaquark_ALL");
std::cout << "\n";

return 1.;
};
