// ---------------------------------------------------------------------------
// Test executable to print out all 24 helicity amplitude to command line.
// run with optional flag -e and -c to set the lab beam energy and cos(theta) in
// center of mass frame.
//
// example: ./test -e 12. -c 1.
// for the forward amplitude at 12 GeV photon beam
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "pomeron_exchange.hpp"

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

  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");
  pomeron_exchange  amp(ptr);

  vector<double> params = {0.379, 0.941, 0.364, 0.12};
  amp.set_params(params);


  double s = mPro * (2.l * egam + mPro);
  for (int i = 0; i < 24; i++)
  {
    cout << std::left << setw(5) << i << setw(15) << amp.helicity_amplitude(ptr->helicities[i], s, zs) << endl;
  }

  delete ptr;
  return 1.;
};
