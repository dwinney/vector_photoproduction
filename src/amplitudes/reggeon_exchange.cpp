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
  // for (int i = 0; i < 24; i++)
  // {
  //   int lamp_gam = kinematics->helicities[i][0];
  //   int lamp_targ = kinematics->helicities[i][1];
  //   int lamp_vec = kinematics->helicities[i][2];
  //   int lamp_rec = kinematics->helicities[i][3];
  //
  //   // wigner functions with appropriate crossing angles
  //   std::complex<double> temp = xr;
  //   temp *= wigner_d_int(1, lamp_gam, lam_gam, M_PI);
  //   temp *= wigner_d_int(1, lamp_vec, lam_vec, -kinematics->crossing_angle(kinematics->vector_particle, s, zs));
  //   temp *= wigner_d_half(1, lamp_targ, lam_targ, kinematics->crossing_angle("target", s, zs));
  //   temp *= wigner_d_half(1, lamp_rec, lam_rec, kinematics->crossing_angle("recoil", s, zs));
  //   temp *= t_channel_amplitude(kinematics->helicities[i], s, zs);
  //
  //   result += temp;
  // }

  result = t_channel_amplitude(helicities, s, zs);
  return result;
};

// ---------------------------------------------------------------------------
// Helicity amplitude in terms of t-channel (unrotated) helicities
std::complex<double> reggeon_exchange::t_channel_amplitude(std::vector<int> helicities, double s, double zs)
{
  // Net helicities
  int lam  = helicities[0] - helicities[2];
  int lamp = (helicities[1] - helicities[3]) / 2.;

  // Center of mass energy and scattering angle
  double t = kinematics->t_man(s, zs);
  std::complex<double> zt = kinematics->z_t(s, zs);

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
  switch (std::abs(lam))
  {
    case 2:
    {
      result = (t - kinematics->mVec2) / sqrt(4. * t * xr);
      break;
    }
    case 1:
    {
      result = (t - kinematics->mVec2) / kinematics->mVec;
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
  switch (std::abs(lamp))
  {
    case 1:
    {
      result = - 2. * mPro2;
      break;
    }
    case 0:
    {
      result = sqrt(xr * 2. * t);
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
