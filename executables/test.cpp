#include "constants.hpp"
#include "polarization_vector.hpp"
#include "two_body_state.hpp"
#include "reaction_kinematics.hpp"

#include <iostream>
using std::cout;
using std::endl;

int main()
{

  reaction_kinematics(mJpsi, "jpsi");
  return 1.;
};
