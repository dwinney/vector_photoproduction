// ---------------------------------------------------------------------------
// Prediction for the photoproduction of X(3872) or Chi_c1(3872) at the future
// EIC at low momentum transfer and high center of mass energy.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/reggeon_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "jpacGraph1D.hpp"
#include "regge_trajectory.hpp"
#include "utilities.hpp"

#include <cstring>
#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <tuple>

int main( int argc, char** argv )
{
  double theta = 0.;
  std::string filename = "X3872_photoproduction.pdf";
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
  }

  // Set up kinematics for the X(3872)
  reaction_kinematics * ptr = new reaction_kinematics(3.872, "X(3872)");

  // Linear trajectory for the rho
  linear_trajectory alpha(1, 0.5, 0.9, "EXD_linear");

  // Initialize Reggeon amplitude with the above kinematics and regge_trajectory
  reggeon_exchange regge(ptr, &alpha, "#rho");
  regge.set_params({0.20, 2.4, 14.6});
  regge.set_signature(-1);

  // Fixed spin rho-exhcange with same couplings for comparison
  vector_exchange not_regge(ptr, .770, "#rho");
  not_regge.set_params({0.20, 2.4, 14.6});

  int N = 4500; // how many points to plot

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Print contributions from each exchange seperately
  double zs = cos(theta * deg2rad);
  double max = 1.E4;

  std::vector<double> s, dxs, NRdxs;
  std::vector<std::complex<double>> mom;
  for (int i = 0; i < N; i++)
  {
    double si = (ptr->sth + 100.*EPS) + double(i) * (max - (ptr->sth + EPS)) / N;
    s.push_back(si);

    // Divide by 4 to average over helicities in the final state
    double xsi, nrxsi;
    xsi   = regge.differential_xsection(si, zs) / 4.;
    nrxsi = not_regge.differential_xsection(si, zs) / 4.;

    dxs.push_back(xsi * 200. * 1.E-3);
    NRdxs.push_back(nrxsi * 1.E-3);
  }

  // ---------------------------------------------------------------------------
  // Plotting Settings
  // ---------------------------------------------------------------------------

  // Initialize plotting object
  jpacGraph1D* plotter = new jpacGraph1D();

  // Add the data sets from above
  plotter->AddEntry(s, NRdxs, "#rho exchange");
  plotter->AddEntry(s, dxs, "200 x #alpha_{#rho}(t)");

// // Tweak the axes
  plotter->SetYaxis(ROOT_italics("d#sigma/dt") + "  (#mub GeV^{-2})", 0., 6.);

  plotter->SetXlogscale(true);
  plotter->SetXaxis(ROOT_italics("s") + "  (GeV^{2})", ptr->sth, max);

  plotter->SetLegend(0.6, 0.5);
  plotter->Plot(filename);
}
