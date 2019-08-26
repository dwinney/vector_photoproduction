#include "constants.hpp"

#include "polarization_vector.hpp"
#include "dirac_spinor.hpp"
#include "two_body_state.hpp"
#include "reaction_kinematics.hpp"

#include <iostream>
using std::cout;
using std::endl;

int main()
{

  reaction_kinematics amp(mJpsi, "jpsi");
  return 1.;
};
