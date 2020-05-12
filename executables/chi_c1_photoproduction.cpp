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
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "utilities.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

int main( int argc, char** argv )
{
  std::cout << "n";

  // Default values
  double theta = 0.;
  bool INTEG = false;
  bool FEYN = false;
  double upper_limit = 7.;
  std::string filename = "chi_c1_photoproduction.pdf";
  std::string ylabel = ROOT_italics("d#sigma/dt") + "  (" + ROOT_italics("nB") + " / GeV^{2})";
  // ---------------------------------------------------------------------------
  // Parse command line arguments
  for (int i = 0; i < argc; i++)
  {
    // File name of output
    if (std::strcmp(argv[i],"-f")==0)
    {
      filename = argv[i+1];
    }
    // Scattering angle in CoM
    if (std::strcmp(argv[i],"-c")==0)
    {
      theta = atof(argv[i+1]);
      if (std::abs(theta) < 90.)
      {
        upper_limit = .5;
      }
    }
    // Plot ixs or dxs
    if (std::strcmp(argv[i],"-integ")==0)
    {
       INTEG = true;
       upper_limit = 13.;
       ylabel = ROOT_italics("#sigma") + "  (" + ROOT_italics("nB") + ")";
    }
    // Use feynman rules or analytic expressions
    if (std::strcmp(argv[i],"-feyn")==0)
    {
      FEYN = true;
    }
  }

  // ---------------------------------------------------------------------------
  // Plotter object
  jpacGraph1D* plotter = new jpacGraph1D();
  plotter->SetXaxis(ROOT_italics("s") + "  (GeV^{2})", 19.5, 100.);

  // To change the range of the Y-axis or the position of the Legend change the arguments here
  plotter->SetYaxis(ylabel , 0., upper_limit);
  plotter->SetLegend(0.2, .75);

  // Set up kinematics for the chi_c1
  reaction_kinematics * ptr = new reaction_kinematics(3.510, "chi_c1");

  // Reaction proceeds through multiple vector exchanges.
  // Which we will sum incoherently
  std::vector<amplitude*> exchanges;

  vector_exchange omega(ptr, .780, "#omega", FEYN);
  omega.set_params({5.2E-4, 16., 0.}); // hadronic and nucleon vector & tensor couplings
  exchanges.push_back(&omega); // Add to the sum vector

  vector_exchange rho(ptr, .770, "#rho", FEYN);
  rho.set_params({9.2E-4, 2.4, 14.6});
  exchanges.push_back(&rho);

  vector_exchange phi(ptr, 1.10, "#phi", FEYN);
  phi.set_params({4.2E-4, -6.2, 2.1});
  exchanges.push_back(&phi);

  vector_exchange jpsi(ptr, 3.097, "J/#psi", FEYN);
  jpsi.set_params({1., 3.3E-3, 0.});
  exchanges.push_back(&jpsi);

  // The total amplitude with all the above exchanges
  amplitude_sum total(ptr, exchanges);

  int N = 50; // how many points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------
double zs = cos(theta * deg2rad);

// ---------------------------------------------------------------------------
// Print the total cross-section
std::vector<double> s, dxs;
for (int i = 0; i <= N; i++)
{
  double si = (ptr->sth + EPS) + double(i) * (100. - ptr->sth) / N;

  double dxsi;
  if (INTEG == false)
  {
    dxsi = total.differential_xsection(si, zs);
  }
  else
  {
    dxsi = total.integrated_xsection(si);
  }

  s.push_back(si);
  dxs.push_back(dxsi);
}

plotter->AddEntry(s, dxs, "Sum");

// ---------------------------------------------------------------------------
// Print contributions from each exchange seperately
for (int n = 0; n < exchanges.size(); n++)
{
  std::vector<double> s, dxs;
  for (int i = 0; i <= N; i++)
  {
    double si = (ptr->sth + EPS) + double(i) * (100. - (ptr->sth + EPS)) / N;
    double dxsi;

    if (INTEG == false)
    {
      dxsi = exchanges[n]->differential_xsection(si, zs);
    }
    else
    {
      dxsi = exchanges[n]->integrated_xsection(si);
    }

    s.push_back(si);
    dxs.push_back(dxsi);
  }

  plotter->AddEntry(s, dxs, exchanges[n]->identifier);
}

// Output to file
plotter->Plot(filename);

delete ptr, plotter;
return 1.;
};
