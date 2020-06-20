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
#include "jpacUtils.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
  // ---------------------------------------------------------------------------
  // COMMAND LINE OPTIONS
  // ---------------------------------------------------------------------------

  // Default values
  double zs = 1.;
  double y[2]; bool custom_y = false;
  double W = 4.45;
  int N = 200; // how many points to plot
  bool TENQ = false;
  double max = 5.;
  std::string filename = "polarized_5q.pdf";
  std::string observable = "dxs", ylabel;

  // Parse input string
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) zs = cos(atof(argv[i+1]) * deg2rad);
    if (std::strcmp(argv[i],"-e")==0) W = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-m")==0) max = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-o")==0) observable = argv[i+1];
    if (std::strcmp(argv[i],"-10q")==0) TENQ = true;
    if (std::strcmp(argv[i],"-y")==0)
    {
      custom_y = true;
      y_range(argv[i+1], y);
    }
  }

  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // Set up Kinematics
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");

  // ---------------------------------------------------------------------------
  // S - CHANNEL

  // Two different pentaquarks
  // masses and widths from 2015 LHCb paper [2]
  baryon_resonance P_c4450(ptr, 3, -1, 4.45, 0.040, "P_{c}(4450)");
  P_c4450.set_params({0.01, .7071}); // 2% branching fraction and equal photocouplings

  baryon_resonance P_c4380(ptr, 5, +1, 4.38, 0.205, "P_{c}(4380)");
  P_c4380.set_params({0.01, .7071}); // 2% branching fraction and equal photocouplings

  // ---------------------------------------------------------------------------
  // T - CHANNEL

  // Set up pomeron trajectory
  // Best fit values from [1]
  linear_trajectory alpha(+1, 0.941, 0.364, "pomeron");

  // Create amplitude with kinematics and trajectory
  pomeron_exchange background(ptr, &alpha, "Background");

  // normalization and t-slope
  background.set_params({0.379, 0.12});

  // ---------------------------------------------------------------------------
  // SUM
  // ---------------------------------------------------------------------------
  // if 10q then include both pentaquarks
  std::vector<amplitude*> amps {&background, &P_c4450};
  if (TENQ == true) (amps.push_back(&P_c4380));

  // Incoherent sum of the s and t channels
  amplitude_sum sum(ptr, amps, "Sum");
  amps.push_back(&sum);

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // Plotter objects
  jpacGraph1D* plotter = new jpacGraph1D();

  // ---------------------------------------------------------------------------
  // scan over energy
  for (int n = 0; n < amps.size(); n++)
  {
    // find the desired observable
    auto F = [&](double W)
    {
      if (observable == "dxs")
      {
        ylabel = "d#sigma/dt    (nb GeV^{-2})";
        return amps[n]->differential_xsection(W*W, zs);
      }
      else if (observable == "kll")
      {
        ylabel = "K_{LL}";
        return amps[n]->K_LL(W*W, zs);
      }
      else if (observable == "all")
      {
        ylabel = "A_{LL}";
        return amps[n]->A_LL(W*W, zs);
      }
      else
      {
        std::cout << "invalid observable passed. Quitting.... \n"; exit(0);
      }
    };

    std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, sqrt(ptr->sth) + 0.01, max);
    plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
  }

  plotter->SetLegend(0.2, 0.7);

  // X axis
  plotter->SetXaxis("W  (GeV)", sqrt(ptr->sth) + 0.01, max);

  // To change the range of the Y-axis or the position of the Legend change the arguments here
  (custom_y == true) ? (plotter->SetYaxis(ylabel, y[0], y[1])) : (plotter->SetYaxis(ylabel));


  plotter->Plot(filename);

return 1.;
};
