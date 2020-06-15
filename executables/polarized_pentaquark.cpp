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
#include "reaction_kinematics.hpp"
#include "amplitudes/baryon_resonance.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
  double theta = 0.;
  double y[2] = {-0.1, 0.1};
  std::string filename = "polarized_pentaquark.pdf";
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-y")==0) y_range(argv[i+1], y);
  }

  // Set up Kinematics
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");

  // Incoherent sum of the s and t channels
  amplitude_sum sum(ptr);

  // ---------------------------------------------------------------------------
  // S - CHANNEL

  // Two different pentaquarks
  // masses and widths from 2015 LHCb paper [2]
  baryon_resonance P_c1(ptr, 3, -1, 4.45, 0.040, "P_{c}(4450)");
  baryon_resonance P_c2(ptr, 5, +1, 4.38, 0.205, "P_{c}(4380)");

  // 2% branching fraction and equal photocouplings for both
  std::vector<double> params = {0.01, .7071};
  P_c1.set_params(params);
  P_c2.set_params(params);

  // Add them to the sum
  sum.add_amplitude(&P_c1);
  sum.add_amplitude(&P_c2);

  // ---------------------------------------------------------------------------
  // T - CHANNEL

  // Set up pomeron trajectory
  linear_trajectory alpha(+1, 0.941, 0.364, "pomeron");

  // Create amplitude with kinematics and trajectory
  pomeron_exchange background(ptr, &alpha);

  // normalization and t-slope
  std::vector<double> back_params = {0.367, 0.12};
  background.set_params(back_params);

  // Add to the sum
  sum.add_amplitude(&background);

  int N = 200; // how many points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// Plotter objects
jpacGraph1D* plotter = new jpacGraph1D();

// ---------------------------------------------------------------------------
// Print .dat files of differential cross section and plot it to a .pdf

double zs = cos(theta * deg2rad);
std::vector<double> s, kll, all;
for (int i = 1; i < N; i++)
{
  double si = (ptr->sth + 10.*EPS) + double(i) * (30. - (ptr->sth + 10.*EPS)) / N;
  s.push_back((si/mPro - mPro)/2.);

  double klli = sum.K_LL(si, zs);
  kll.push_back(klli);

  double alli = sum.A_LL(si, zs);
  all.push_back(alli);
}

plotter->AddEntry(s, kll, "K_{LL}");
plotter->AddEntry(s, all, "A_{LL}");

plotter->SetLegend(0.7, 0.2);
plotter->SetXaxis("E_{#gamma}  (GeV)", 8.5, 12.);
plotter->SetYaxis("", y[0], y[1]);

plotter->Plot(filename);

return 1.;
};
