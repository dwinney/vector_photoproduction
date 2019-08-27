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
#include "two_body_state.hpp"
#include "dirac_spinor.hpp"
#include "polarization_vector.hpp"

#include <vector>
#include <string>

using std::vector;
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
  const double mVec;
  const string vector_particle;

public:
  // Constructor
  reaction_kinematics(double vec_mass, string vec_name)
  : mVec(vec_mass), vector_particle(vec_name),
    initial(0., mPro, "beam", "target"),
    final(vec_mass, mPro, vec_name, "recoil"),
    eps_vec_star(final, vec_name, true),
    eps_gamma(initial, "beam"),
    target(initial, "target"), recoil(final, "recoil")
  {};

// Copy Constructor
  reaction_kinematics(const reaction_kinematics & old)
  : mVec(old.mVec), vector_particle(old.vector_particle),
    initial(old.initial), final(old.final),
    eps_vec_star(old.eps_vec_star), eps_gamma(old.eps_gamma),
    target(old.target), recoil(old.recoil)
  {};

  // inital and final state kinematics
  const double sth = (mVec + mPro) * (mVec + mPro);
  two_body_state initial, final;
  polarization_vector eps_vec_star, eps_gamma;
  dirac_spinor target, recoil;

  double t_man(double s, double zs)
  {
    complex<double> ks = initial.momentum("beam", s) * final.momentum("vec_name", s);
    return (2. * mPro * mPro + mVec * mVec  - s) / 2. + 2. * abs(ks) * zs;
  };

  // Helicity configurations for Photon [0], Incoming Proton [1], vector m [2], Outgoing Proton [3]
  vector<vector<double>> helicities =
  {
    {  1., -1.,  1., -1. },
    {  1., -1.,  1.,  1. },
    {  1., -1.,  0., -1. },
    {  1., -1.,  0.,  1. },
    {  1., -1., -1., -1. },
    {  1., -1., -1.,  1. },
    {  1.,  1.,  1., -1. },
    {  1.,  1.,  1.,  1. },
    {  1.,  1.,  0., -1. },
    {  1.,  1.,  0.,  1. },
    {  1.,  1., -1., -1. },
    {  1.,  1., -1.,  1. },
    { -1., -1.,  1., -1. },
    { -1., -1.,  1.,  1. },
    { -1., -1.,  0., -1. },
    { -1., -1.,  0.,  1. },
    { -1., -1., -1., -1. },
    { -1., -1., -1.,  1. },
    { -1.,  1.,  1., -1. },
    { -1.,  1.,  1.,  1. },
    { -1.,  1.,  0., -1. },
    { -1.,  1.,  0.,  1. },
    { -1.,  1., -1., -1. },
    { -1.,  1., -1.,  1. },
  };
};

#endif
