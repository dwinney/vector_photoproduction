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
    eps_vec(final, vec_name),
    eps_gamma(initial, "beam"),
    target(initial, "target"), recoil(final, "recoil")
  {};

// Copy Constructor
  reaction_kinematics(const reaction_kinematics & old)
  : mVec(old.mVec), vector_particle(old.vector_particle),
    initial(old.initial), final(old.final),
    eps_vec(old.eps_vec), eps_gamma(old.eps_gamma),
    target(old.target), recoil(old.recoil)
  {};

  // inital and final state kinematics
  const double sth = (mVec + mPro) * (mVec + mPro);
  two_body_state initial, final;
  polarization_vector eps_vec, eps_gamma;
  dirac_spinor target, recoil;

  double t_man(double s, double zs)
  {
    complex<double> kq = initial.momentum("beam", s) * final.momentum(vector_particle, s);
    complex<double> E1E3 = initial.energy("beam", s) * final.energy(vector_particle, s);
    return mVec*mVec - 2. * abs(E1E3) + 2. * abs(kq) * zs;
  };

  // Helicity configurations
  // Photon [0], Incoming Proton [1], Vector meson [2], Outgoing Proton [3]
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
