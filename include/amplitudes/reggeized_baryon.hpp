// Parameterization of a u-channel background amplitude with couplings from the
// s-channel
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BARYEGGEON_
#define _BARYEGGEON_

#include "baryon_resonance.hpp"
#include "regge_trajectory.hpp"
#include "misc_math.hpp"

class reggeized_baryon : public baryon_resonance
{
public:
  //Constructor
  reggeized_baryon(reaction_kinematics * xkinem, linear_trajectory * traj, int j, double mass, double width, std::string exchange = "")
  : baryon_resonance(xkinem, j, traj->signature, mass, width, exchange),
    alpha(traj)
  {};

  // Copy constructor
  reggeized_baryon(const reggeized_baryon & old)
  : baryon_resonance(old), alpha(old.alpha)
  {};

  // Assemble the helicity amplitude
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

private:
  linear_trajectory * alpha;

  // Usual reggeon propagator
  std::complex<double> regge_propagator(double t);

  // Half angle factors
  std::complex<double> half_angle_factor(int lam, int lamp, std::complex<double> z_t);

  // Helicity amplitude in terms of t-channel (unrotated) helicities
  std::complex<double> u_channel_amplitude(std::vector<int> helicities, double s, double zs);
};

#endif
