// Axial-vector meson photoproduction proceeding through a Reggeized meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/reggeon_exchange.hpp"

// ---------------------------------------------------------------------------
// Assemble the helicity amplitude by contracting the lorentz indices
std::complex<double> reggeon_exchange::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  // Sum over all helicities in the t - channel
  std::complex<double> result = 0.;
  for (int i = 0; i < 24; i++)
  {
    int lamp_gam = kinematics->helicities[i][0];
    int lamp_targ = kinematics->helicities[i][1];
    int lamp_vec = kinematics->helicities[i][2];
    int lamp_rec = kinematics->helicities[i][3];

    // wigner functions with appropriate crossing angles
    std::complex<double> temp = xr;
    temp *= wigner_d_int(1, lamp_gam, lam_gam, crossing_angle("beam", s, zs));
    temp *= wigner_d_int(1, lamp_vec, lam_vec, crossing_angle(kinematics->vector_particle, s, zs));
    temp *= wigner_d_half(1, lamp_targ, lam_targ, crossing_angle("target", s, zs));
    temp *= wigner_d_half(1, lamp_rec, lam_rec, crossing_angle("recoil", s, zs));
    temp *= t_channel_amplitude(kinematics->helicities[i], s, zs);

    result += temp;
  }

  return result;
};

// ---------------------------------------------------------------------------
// Helicity amplitude in terms of t-channel (unrotated) helicities
std::complex<double> reggeon_exchange::t_channel_amplitude(std::vector<int> helicities, double s, double zs)
{
  // Net helicities
  int lam  = helicities[0] - helicities[2];
  int lamp = (helicities[1] - helicities[3]) / 2.;

  if (abs(lam) == 2)
  {
    return 0.;
  }

  // Center of mass energy and scattering angle
  double t = kinematics->t_man(s, zs);
  std::complex<double> zt = z_t(s, zs);

  // Product of residues
  std::complex<double> result = xr;
  result *= top_vertex(lam, t);
  result *= bottom_vertex(lamp, t);

  // Pole piece
  result /= t - mEx2;

  // angular function
  result *= wigner_d_int(1, lamp, lam, zt);

  return result;
};

// ---------------------------------------------------------------------------
// Photon - Axial Vector - Vector vertex
std::complex<double> reggeon_exchange::top_vertex(int lam, double t)
{
  std::complex<double> result;
  switch (abs(lam))
  {
    case 2:
    {
      result = (t - kinematics->mVec2) / sqrt(4. * t * xr);
      break;
    }
    case 1:
    {
      result = (t / (2. * kinematics->mVec2) - 1.);
      break;
    }
    case 0:
    {
      return 0.;
    }
    default:
    {
      std::cout << "\nreggeon_exchange: invalid helicity flip lambda = " << lam << ". Quitting... \n";
      exit(0);
    }
  }

  return result * gGamma;
};

// ---------------------------------------------------------------------------
// Nucleon - Nucleon - Vector vertex
std::complex<double> reggeon_exchange::bottom_vertex(int lamp, double t)
{
  std::complex<double> result;
  switch (abs(lamp))
  {
    case 1:
    {
      result = - 2. * mPro2;
      break;
    }
    case 0:
    {
      result = sqrt(2. * xr * t);
      break;
    }
    default:
    {
      std::cout << "\nreggeon_exchange: invalid helicity flip lambda^prime = " << lamp << ". Quitting... \n";
      exit(0);
    }
  }

  return result * gV;
};

// ---------------------------------------------------------------------------
// Scattering angle in the t-channel
std::complex<double> reggeon_exchange::z_t(double s, double zs)
{
  double t = kinematics->t_man(s, zs);

  std::complex<double> result;
  result = 2. * s + t - 2. * mPro2 - kinematics->mVec2;
  result /= t - kinematics->mVec2;
  result *= sqrt(xr * t) / sqrt(xr * t - 4. * mPro2);

  return result;
};

// ---------------------------------------------------------------------------
// Cosine of crossing angles
std::complex<double> reggeon_exchange:: crossing_angle(std::string particle, double s, double zs)
{
  double t = kinematics->t_man(s, zs);

  std::complex<double> S_gv, S_pp, T_gp, T_vp;
  S_gv = sqrt(xr * Kallen(s, 0., kinematics->mVec2));
  S_pp = sqrt(xr * Kallen(s, mPro2, mPro2));
  T_gp = sqrt(xr * Kallen(t, 0., mPro2));
  T_vp = sqrt(xr * Kallen(t, kinematics->mVec2, mPro2));

  std::complex<double> result;
  if (particle == "beam")
  {
    result = - (s - kinematics->mVec2) * (t - mPro2);
    result /= S_gv * T_gp;
  }
  else if (particle == kinematics->vector_particle)
  {
    result = (s + kinematics->mVec2) * ( t + kinematics->mVec2 - mPro2) - 2. * pow(kinematics->mVec2, 2.);
    result /= S_gv * T_vp;
  }
  else if (particle == "target")
  {
    result = - s * (t + mPro2 - kinematics->mVec2) - 2.*mPro2*kinematics->mVec2;
    result /= S_pp * T_vp;
  }
  else if (particle == "recoil")
  {
    result = s * (t - mPro2) - 2.*mPro2*kinematics->mVec2;
    result /= S_pp * T_gp;
  }
  else
  {
    std::cout << "\nreggeon_exchange: particle = " << particle << " not recognized in crossing_angle. Quitting... \n";
    exit(0);
  }

  return result;
};
