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
  baryon_resonance  amp(3, 4.45, 0.04, ptr);

  std::vector<double> params = {0.02, .7071};
  amp.set_params(params);

  cout << std::right << setw(5) << " ";
  cout << setw(10) << "lam_gam" << setw(10) << "lam_targ" << setw(10) << "lam_vec" << setw(10) << "lam_rec";
  cout << setw(25) << "helicity_amplitude" << endl;

  double s = mPro * (2.l * egam + mPro);
  for (int i = 0; i < 24; i++)
  {
    cout << std::right << setw(5) << i;
    cout << setw(10) << ptr->helicities[i][0];
    cout << setw(10) << ptr->helicities[i][1];
    cout << setw(10) << ptr->helicities[i][2];
    cout << setw(10) << ptr->helicities[i][3];
    cout << setw(25) << amp.helicity_amplitude(ptr->helicities[i], s, zs) << endl;
  }

  delete ptr;
  return 1.;
};
