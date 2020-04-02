// ---------------------------------------------------------------------------
// Test executable to print out all 24 helicity amplitude to command line.
// run with optional flag -e and -c to set the lab beam energy and cos(theta) in
// center of mass frame.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "regge_trajectory.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using std::setw;
using std::cout;
using std::endl;

int main( int argc, char** argv )
{
  double sqrts = 5.5;
  double zs = .7071;

  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-e")==0) sqrts = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-c")==0) zs = atof(argv[i+1]);
  }

  // Set up kinematics, determined entirely by vector meson mass
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");

  // Set up pomeron trajectory
  // here we use linear trajectory with intercept and slope only free params
  linear_traj alpha(0.941, 0.364);

  // Create amplitude with kinematics and trajectory
  pomeron_exchange t_channel(ptr, &alpha);

  // Feed in other two parameters (normalization and t-slope)
  std::vector<double> params = {0.379, 0.12};
  t_channel.set_params(params);

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Print out helicity amplitudes to cmd line

  cout << std::right << setw(5) << " ";
  cout << setw(10) << "lam_gam";
  cout << setw(10) << "lam_targ";
  cout << setw(10) << "lam_vec";
  cout << setw(10) << "lam_rec";
  cout << setw(25) << "helicity_amplitude" << endl;

  for (int i = 0; i < 24; i++)
  {
    cout << std::right << setw(5) << i;
    cout << setw(10) << ptr->helicities[i][0];
    cout << setw(10) << ptr->helicities[i][1];
    cout << setw(10) << ptr->helicities[i][2];
    cout << setw(10) << ptr->helicities[i][3];
    cout << setw(25) << t_channel.helicity_amplitude(ptr->helicities[i], sqrts*sqrts, zs) << endl;
  }
  std::cout << "\n";

  delete ptr;
  return 1.;
};
