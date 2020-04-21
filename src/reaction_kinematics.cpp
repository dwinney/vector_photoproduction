// Class to contain all relevant kinematic quantities. The kinematics of the reaction
// gamma p -> V p' is entirely determined by specifying the mass of the vector particle
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"

// ---------------------------------------------------------------------------
// Invariant variables
double reaction_kinematics::t_man(double s, double zs)
{
  std::complex<double> kq = initial.momentum("beam", s) * final.momentum(vector_particle, s);
  std::complex<double> E1E3 = initial.energy("beam", s) * final.energy(vector_particle, s);
  return mVec*mVec - 2. * abs(E1E3) + 2. * abs(kq) * zs;
};

// ---------------------------------------------------------------------------
// Scattering angle in the t-channel
std::complex<double> reaction_kinematics::z_t(double s, double zs)
{
  double t = t_man(s, zs);
  std::complex<double> p_t = sqrt(xr * t - 4. * mPro2) / 2.;
  std::complex<double> q_t = sqrt(xr * Kallen(t, mVec2, 0.)) / sqrt(xr * 4. * t);

  std::complex<double> result;
  result = 2. * s + t - 2. * mPro2 - mVec2; // s - u
  result /= 4. * p_t * q_t;

  return result;
};

// ---------------------------------------------------------------------------
// Cosine of crossing angles
std::complex<double> reaction_kinematics::crossing_angle(std::string particle, double s, double zs)
{
  double t = t_man(s, zs);

  std::complex<double> S_gp, S_vp, T_gv, T_pp;
  S_gp = sqrt(xr * Kallen(s, 0., mPro2));
  S_vp = sqrt(xr * Kallen(s, mVec2, mPro2));
  T_gv = sqrt(xr * Kallen(t, 0., mVec2));
  T_pp = sqrt(xr * Kallen(t, mPro2, mPro2));

  std::complex<double> result;
  if (particle == "beam")
  {
    result = - (s - mPro2) * (t - mVec2);
    result /= S_gp * T_gv;
  }
  else if (particle == vector_particle)
  {
    result = -(s + mVec2 - mPro2) * ( t + mVec2) + 2. * pow(mVec2, 2.);
    result /= S_vp * T_gv;
  }
  else if (particle == "target")
  {
    result = - (s + mPro2) * t + 2.*mPro2*mVec2;
    result /= S_gp * T_pp;
  }
  else if (particle == "recoil")
  {
    result = (s + mPro2 - mVec2) * t + 2.*mPro2*mVec2;
    result /= S_vp * T_pp;
  }
  else
  {
    std::cout << "\nreaction_kinematics: particle = " << particle << " not recognized in crossing_angle. Quitting... \n";
    exit(0);
  }

  return result;
};
