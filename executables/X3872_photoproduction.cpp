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

#include "amplitudes/reggeized_meson.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>


int main( int argc, char** argv )
{
  int N = 100;
  double theta = 0.;
  bool integ = false;
  std::string filename = "X3872_photoproduction.pdf";
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-n")==0) N = atoi(argv[i+1]);
    if (std::strcmp(argv[i],"-integ")==0) integ = true;
  }

  // Set up kinematics for the X(3872)
  reaction_kinematics * ptr = new reaction_kinematics(3.872, "X(3872)");

  // Linear trajectory for the rho
  linear_trajectory alpha(-1, 0.5, 0.9, "EXD_linear");

  // Initialize Reggeon amplitude with the above kinematics and regge_trajectory
  reggeized_meson rho(ptr, &alpha, "#rho");
  rho.set_params({3.81E-3, 2.4, 14.6});

  reggeized_meson omega(ptr, &alpha, "#omega");
  omega.set_params({9.51E-3, 16, 0.});

  std::vector<amplitude*> exchanges = {&rho, &omega};
  amplitude_sum total(ptr, exchanges);

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

  std::vector<double> s, dxs;
  for (int i = 0; i < N; i++)
  {
    double si = (ptr->sth + 100.*EPS) + double(i) * (max - (ptr->sth + EPS)) / N;
    s.push_back(si);

    // Divide by 4 to average over helicities in the final state
    double xsi;
    if (integ == false)
    {
      xsi   = total.differential_xsection(si, zs) / 4.;
      dxs.push_back(xsi);
    }
    else
    {
      xsi   = total.integrated_xsection(si) / 4.;
      dxs.push_back(xsi);
    }
  }

  // Initialize plotting object
  jpacGraph1D* plotter = new jpacGraph1D();

  plotter->AddEntry(s, dxs, "Sum");

  // ---------------------------------------------------------------------------
  for (int n = 0; n < exchanges.size(); n++)
  {
    std::vector<double> s, dxs;
    for (int i = 0; i <= N; i++)
    {
      double si = (ptr->sth + EPS) + double(i) * (max - (ptr->sth + EPS)) / N;
      double dxsi;

      if (integ == false)
      {
        dxsi = exchanges[n]->differential_xsection(si, zs) / 4.;
      }
      else
      {
        dxsi = exchanges[n]->integrated_xsection(si) / 4.;
      }

      s.push_back(si);
      dxs.push_back(dxsi);
    }

    plotter->AddEntry(s, dxs, exchanges[n]->identifier);
  }

  // ---------------------------------------------------------------------------
  // Plotting Settings
  // ---------------------------------------------------------------------------

 // Tweak the axes
  if (integ == false)
  {
    plotter->SetYaxis(ROOT_italics("d#sigma/dt") + "  (nb GeV^{-2})", 0., 2.0);
    plotter->SetXlogscale(true);
  }
  else
  {
    plotter->SetYaxis("#sigma   (nb)", 0., .15);
  }

  plotter->SetXaxis(ROOT_italics("s") + "  (GeV^{2})", ptr->sth, max);

  plotter->SetLegend(0.73, 0.2);
  plotter->Plot(filename);
}
