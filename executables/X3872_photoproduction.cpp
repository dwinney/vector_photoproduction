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
#include "regge_trajectory.hpp"
#include "utilities.hpp"

#include "amplitudes/reggeized_meson_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"

#include "jpacGraph1D.hpp"


int main( int argc, char** argv )
{
  double theta = 0.;
  bool integ = false;
  bool compare = false;
  std::string filename = "X3872_photoproduction.pdf";
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-integ")==0) integ = true;
    if (std::strcmp(argv[i],"-compare")==0) compare = true;
  }

  // Set up kinematics for the X(3872)
  reaction_kinematics * ptr = new reaction_kinematics(3.872, "X(3872)");

  // Linear trajectory for the rho
  linear_trajectory alpha(1, -1, 0.5, 0.9, "EXD_linear");

  // Initialize Reggeon amplitude with the above kinematics and regge_trajectory
  reggeized_meson_exchange regge(ptr, &alpha, "#rho");
  regge.set_params({0.20, 2.4, 14.6});

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
  double max;
  if (integ == true)
  {
    max = 200.;
  }
  else
  {
    max = 1.E4;
  }

  std::vector<double> s, dxs, NRdxs;
  std::vector<std::complex<double>> mom;
  for (int i = 0; i < N; i++)
  {
    double si = (ptr->sth + 100.*EPS) + double(i) * (max - (ptr->sth + EPS)) / N;
    s.push_back(si);

    // Divide by 4 to average over helicities in the final state
    double xsi, nrxsi;
    if (integ == false)
    {
      xsi   = regge.differential_xsection(si, zs) / 4.;
      dxs.push_back(xsi * 1.E-3);

      if (compare == true)
      {
        nrxsi = not_regge.differential_xsection(si, zs) / 4.;
        NRdxs.push_back(nrxsi * 1.E-6);
      }
    }
    else
    {
      xsi   = regge.integrated_xsection(si) / 4.;
      dxs.push_back(xsi);

      if (compare == true)
      {
        nrxsi = not_regge.integrated_xsection(si) / 4.;
        NRdxs.push_back(nrxsi * 1.E-6);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Plotting Settings
  // ---------------------------------------------------------------------------

  // Initialize plotting object
  jpacGraph1D* plotter = new jpacGraph1D();

  // Add the data sets from above
  plotter->AddEntry(s, dxs, "#alpha_{#rho}(t)");

  if (compare == true)
  {
    plotter->AddEntry(s, NRdxs, "#rho exchange x 10^{-3}");
  }

// // Tweak the axes
  if (integ == false)
  {
    plotter->SetYaxis(ROOT_italics("d#sigma/dt") + "  (#mub GeV^{-2})");
    plotter->SetXlogscale(true);
  }
  else
  {
    plotter->SetYaxis("#sigma   (nb)");
  }

  plotter->SetXaxis(ROOT_italics("s") + "  (GeV^{2})", ptr->sth, max);

  plotter->SetLegend(0.6, 0.5);
  plotter->Plot(filename);
}
