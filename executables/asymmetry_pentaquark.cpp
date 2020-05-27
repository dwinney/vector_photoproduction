// ---------------------------------------------------------------------------
// Predicted sensativity to LHCb pentaquarks in beam and parity asymmetries at
// GlueX at JLab
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

int main( int argc, char** argv )
{
  double y[2] = {-0.4, 0.3};
  double W = 4.45;
  bool LAB = false, TENQ = false;
  std::string filename = "5q_beam_asymmetry.pdf";
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-e")==0) W = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-y")==0) y_range(argv[i+1], y);
    if (std::strcmp(argv[i],"-lab")==0) LAB = true;
    if (std::strcmp(argv[i],"-10q")==0) TENQ = true;
  }

  // Set up Kinematics
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");
  std::vector<amplitude*> amps;

  // ---------------------------------------------------------------------------
  // T - CHANNEL

  // Set up pomeron trajectory
  linear_trajectory alpha(+1, 0.941, 0.364);

  // Create amplitude with kinematics and trajectory
  pomeron_exchange background(ptr, &alpha, "Background");

  // normalization and t-slope
  std::vector<double> back_params = {0.367, 0.12};
  background.set_params(back_params);

  // ---------------------------------------------------------------------------
  // S - CHANNEL  // Two different pentaquarks

  // masses and widths from 2015 LHCb paper [2]
  baryon_resonance P_c4450(ptr, 3, -1, 4.45, 0.040, "P_{c}(4450)");
  P_c4450.set_params({0.01, .7071});

  // 1% branching fraction and equal photocouplings for both
  baryon_resonance P_c4380(ptr, 5, +1, 4.38, 0.205, "P_{c}(4380)");
  P_c4380.set_params({0.01, .7071});

  // Incoherent sum of the s and t channels
  amplitude_sum sum(ptr, {&background, &P_c4450, &P_c4380}, "Sum");

  // ---------------------------------------------------------------------------
  // S - CHANNEL  // 1 5q different BR scenarios

    baryon_resonance P_c1(ptr, 3, -1, 4.45, 0.040, "1%");
    P_c1.set_params({0.01, .7071});

    baryon_resonance P_c05(ptr, 3, -1, 4.45, 0.040, "0.5%");
    P_c05.set_params({0.005, .7071});

    baryon_resonance P_c01(ptr, 3, -1, 4.45, 0.040, "0.1%");
    P_c01.set_params({0.001, .7071});

    // Add to the sum
    amplitude_sum sum1(ptr, {&background, &P_c1}, "1%");
    amplitude_sum sum2(ptr, {&background, &P_c05}, "0.5%");
    amplitude_sum sum3(ptr, {&background, &P_c01}, "0.1%");

    if (TENQ == true)
    {
      amps = {&sum, &background, &P_c4450, &P_c4380};
    }
    else
    {
      amps = {&background, &sum1, &sum2, &sum3};
    }

  int N = 200; // how many points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// Plotter objects
jpacGraph1D* plotter = new jpacGraph1D();

// ---------------------------------------------------------------------------
// scan over theta

double s;
(LAB == true) ? (s = 2.*mPro* W + mPro2) : (s = W*W);

for (int n = 0; n < amps.size(); n++)
{
  std::vector<double> theta, sigma;
  for (int i = 1; i <= N; i++)
  {
    double theta_i = double(i) * 90. / N;
    theta.push_back(theta_i);

    double sigma_i = amps[n]->beam_asymmetry(s, cos(theta_i * deg2rad));
    sigma.push_back(sigma_i);
  }

  plotter->AddEntry(theta, sigma, amps[n]->identifier);
}

std::ostringstream streamObj;
streamObj << std::setprecision(4);
streamObj << W;

std::string header;
(LAB == true) ? (header = "E_{#gamma} = ") : (header = "W = ");
plotter->SetLegend(0.18, 0.7, header + streamObj.str() + " GeV");
plotter->SetXaxis("#theta  (GeV)", 0., 90.);
plotter->SetYaxis("#Sigma_{4#pi}", y[0], y[1]);

plotter->Plot(filename.c_str());

return 1.;
};
