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

#include "jpacGraph1D.hpp"

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

int main( int argc, char** argv )
{
  double theta = 0.;
  bool LAB = false;
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-lab")==0) LAB = true;
  }

  // Set up kinematics, determined entirely by vector meson mass
  reaction_kinematics * ptr1s = new reaction_kinematics(3.097, "#psi(1S)");
  reaction_kinematics * ptr2s = new reaction_kinematics(3.686, "#psi(2S)");

  // Set up pomeron trajectory
  // Here we use (real) linear trajectory with intercept and slope only free params
  // Best fit values from [1]
  linear_trajectory alpha(+1, 0.941, 0.364, "pomeron");

  // Create amplitudes with kinematics and same trajectory
  pomeron_exchange pomeron_1s(ptr1s, &alpha, "#psi(1S)");
  pomeron_exchange pomeron_2s(ptr2s, &alpha, "#psi(2S)");

  // Feed in other two parameters (normalization and t-slope)
  // Best fit from [1]
  std::vector<double> params_1s = {sqrt(4. * M_PI * M_ALPHA) * 0.379, 0.12};
  pomeron_1s.set_params(params_1s);

  // Same t-slope but coupling is scaled by 1/4
  std::vector<double> params_2s = {sqrt(4. * M_PI * M_ALPHA) * 0.379 / 4., 0.12};
  pomeron_2s.set_params(params_2s);

  int N = 100; // how many points to plot

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // Initialize two plotting objects for the ratio and the dxs
  jpacGraph1D* plotter = new jpacGraph1D();
  plotter->SetLegend(0.2, .8);

  if (LAB == true)
  {
    plotter->SetXaxis("E_{#gamma} (GeV)", 10.5, 16.);
  }
  else
  {
    plotter->SetXaxis("W (GeV)", 4.6 , 5.55);
  }

  double zs = cos(theta * deg2rad);
  std::vector<double> s, dxs1, dxs2, ratio;

  for (int i = 1; i <= N; i++)
  {
    double si = ptr2s->sth + EPS + double(i) * (30. - ptr2s->sth - EPS) / N;
    double dxsi1 = pomeron_1s.differential_xsection(si, zs) / 4.;
    double dxsi2 = pomeron_2s.differential_xsection(si, zs) * 100. / 4.;
    double ratioi = dxsi2 / dxsi1;

    //Convert center of mass energy to lab frame energy
    if (LAB == true)
    {
      s.push_back((si/mPro - mPro)/2.);
    }
    else
    {
      s.push_back(sqrt(si));
    }

    dxs1.push_back(dxsi1);
    dxs2.push_back(dxsi2);
    ratio.push_back(ratioi);
  }


  plotter->AddEntry(s, dxs1, "#psi(1S)");
  plotter->AddEntry(s, dxs2, "#psi(2S) x 100");
  plotter->SetYaxis("d#sigma/dt  (nb GeV^{-2})", 0., 15.);
  plotter->SetLegend(0.2, 0.7);
  plotter->Plot("psi_dxs.pdf");

  plotter->ClearData();
  plotter->AddEntry(s, ratio,"");
  plotter->SetLegend(false);
  plotter->SetYaxis("d#sigma(2S)/d#sigma(1S) x 100", 0., 2.);
  plotter->Plot("psi_dxs_ratio.pdf");

  // Clean up pointers
  delete ptr1s, ptr2s;
  delete plotter;

  return 1.;
};
