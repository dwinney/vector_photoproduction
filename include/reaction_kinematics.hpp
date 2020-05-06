// Class to contain all relevant kinematic quantities. The kinematics of the reaction
// gamma p -> V p' is entirely determined by specifying the mass of the vector particle
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------


#ifndef _KINEMATICS_
#define _KINEMATICS_

#include "constants.hpp"
#include "misc_math.hpp"
#include "two_body_state.hpp"
#include "dirac_spinor.hpp"
#include "polarization_vector.hpp"

#include <vector>
#include <string>

class reaction_kinematics
{
public:
  // Constructor
  reaction_kinematics(double vec_mass, std::string vec_name)
  : mVec(vec_mass), mVec2(vec_mass*vec_mass), vector_particle(vec_name),
    initial(0., mPro, "beam", "target"),
    final(vec_mass, mPro, vec_name, "recoil"),
    eps_vec(final, vec_name),
    eps_gamma(initial, "beam"),
    target(initial, "target"), recoil(final, "recoil")
  {};

// Copy Constructor
  reaction_kinematics(const reaction_kinematics & old)
  : mVec(old.mVec), mVec2(old.mVec2), vector_particle(old.vector_particle),
    initial(old.initial), final(old.final),
    eps_vec(old.eps_vec), eps_gamma(old.eps_gamma),
    target(old.target), recoil(old.recoil)
  {};

  const std::string vector_particle;
  const double mVec, mVec2;

  // inital and final state kinematics
  const double sth = (mVec + mPro) * (mVec + mPro);
  two_body_state initial, final;
  polarization_vector eps_vec, eps_gamma;
  dirac_spinor target, recoil;

  // Invariant variables
  double t_man(double s, double zs);
  double u_man(double s, double zs);

  // Scattering angle in the t-channel
  std::complex<double> z_t(double s, double zs);

  // Cosine of crossing angles
  std::complex<double> crossing_angle(std::string particle, double s, double zs);

  // Helicity configurations
  // Photon [0], Incoming Proton [1], Vector meson [2], Outgoing Proton [3]
  std::vector< std::vector<int> > helicities =
  {
    {  1, -1,  1, -1 },
    {  1, -1,  1,  1 },
    {  1, -1,  0, -1 },
    {  1, -1,  0,  1 },
    {  1, -1, -1, -1 },
    {  1, -1, -1,  1 },
    {  1,  1,  1, -1 },
    {  1,  1,  1,  1 },
    {  1,  1,  0, -1 },
    {  1,  1,  0,  1 },
    {  1,  1, -1, -1 },
    {  1,  1, -1,  1 },
    { -1, -1,  1, -1 },
    { -1, -1,  1,  1 },
    { -1, -1,  0, -1 },
    { -1, -1,  0,  1 },
    { -1, -1, -1, -1 },
    { -1, -1, -1,  1 },
    { -1,  1,  1, -1 },
    { -1,  1,  1,  1 },
    { -1,  1,  0, -1 },
    { -1,  1,  0,  1 },
    { -1,  1, -1, -1 },
    { -1,  1, -1,  1 },
  };
};

#endif
