// ---------------------------------------------------------------------------
// Test executable to print out all 24 helicity amplitude to command line.
// run with optional flag -e and -c to set the lab beam energy and cos(theta) in
// center of mass frame.
//
// example: ./test -e 12. -c 1.
// for the forward amplitude at 12 GeV photon beam
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/baryon_resonance.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using std::setw;
using std::cout;
using std::endl;

int main( int argc, char** argv )
{
  double egam = 10.;
  double zs = .7071;

  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-e")==0) egam = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-c")==0) zs = atof(argv[i+1]);
  }

  // Set up Kinematics
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");

  // Incoherent sum of the s and t channels
  amplitude_sum sum;

  // ---------------------------------------------------------------------------
  // S - CHANNEL

  // Two different pentaquarks
  baryon_resonance P_c1(3, 4.45, 0.04, ptr);
  baryon_resonance P_c2(-5, 4.38, 0.01, ptr);

  // 2% branching fraction and equal photocouplings for both
  std::vector<double> params = {0.02, .7071};
  P_c1.set_params(params);  P_c2.set_params(params);

  // Add them to the sum
  sum.add_amplitude(&P_c1);  sum.add_amplitude(&P_c2);

  // ---------------------------------------------------------------------------
  // T - CHANNEL

  // Set up pomeron trajectory
  linear_traj alpha(0.941, 0.364);

  // Create amplitude with kinematics and trajectory
  pomeron_exchange background(ptr, &alpha);

  // normalization and t-slope
  std::vector<double> back_params = {0.379, 0.12};
  background.set_params(back_params);

  // Add to the sum
  sum.add_amplitude(&background);

  // ---------------------------------------------------------------------------
  // Print helicity amplitudes to command-line

  cout << std::right << setw(5) << " ";
  cout << setw(10) << "lam_gam";
  cout << setw(10) << "lam_targ";
  cout << setw(10) << "lam_vec";
  cout << setw(10) << "lam_rec";
  cout << setw(25) << "helicity_amplitude" << endl;

  double s = mPro * (2.l * egam + mPro);
  for (int i = 0; i < 24; i++)
  {
    cout << std::right << setw(5) << i;
    cout << setw(10) << ptr->helicities[i][0];
    cout << setw(10) << ptr->helicities[i][1];
    cout << setw(10) << ptr->helicities[i][2];
    cout << setw(10) << ptr->helicities[i][3];
    cout << setw(25) << sum.helicity_amplitude(ptr->helicities[i], s, zs) << endl;
  }

  delete ptr;
  return 1.;
};
