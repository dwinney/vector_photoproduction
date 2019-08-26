// Class to contain all relevant kinematic quantities. The kinematics of the reaction
// gamma p -> V p' is entirely determined by specifying the mass of the vector particle
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "two_body_state.hpp"
#include "polarization_vector.hpp"

#include <string>

#ifndef _KINEMATICS_
#define _KINEMATICS_

using std::string;

template <typename T>
T Kallen(T x, T y, T z)
{
  return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
};

class reaction_kinematics
{
private:
  // vector meson mass and identifier
  double mVec;
  string vector_particle;

  // inital and final state kinematics
  two_body_state initial, final;
  polarization_vector eps_vec_star, eps_gamma;
  dirac_spinor target, recoil;

public:
  reaction_kinematics(double vec_mass, string vec_name)
  : mVec(vec_mass), vector_particle(vec_name),
    initial(0., mPro, "beam", "target"),
    final(vec_mass, mPro, vec_name, "recoil"),
    eps_vec_star(final, vec_name, true),
    eps_gamma(initial, "beam"),
    target(initial, "target"), recoil(final, "recoil")
  {};

};

#endif
