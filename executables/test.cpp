#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "pomeron_exchange.hpp"

#include <iostream>
#include <iomanip>

using std::setw;
using std::cout;
using std::endl;

int main()
{
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");
  pomeron_exchange  amp(ptr);

  vector<double> params = {0.379, 0.941, 0.364, 0.12};
  amp.set_params(params);

  double egam = 10.;
  double s = mPro * (2.l * egam + mPro);

  for (int i = 0; i < 24; i++)
  {
    cout << std::left << setw(5) << i << setw(15) << amp.helicity_amplitude(ptr->helicities[i], s, 1.) << endl;
  }

  delete ptr;
  return 1.;
};
