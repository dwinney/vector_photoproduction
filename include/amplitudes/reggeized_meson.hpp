// Axial-vector meson photoproduction proceeding through a Reggeized meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _REGGEON_
#define _REGGEON_

#include "vector_exchange.hpp"
#include "regge_trajectory.hpp"
#include "misc_math.hpp"

// ---------------------------------------------------------------------------
// reggeon_exchange class describes the amplitude for a reggeon exchange
// in the t-channel. Residues are given by the feynman amplitude in vector_exchange
//
// Initialization required a reaction_kinematics object, a linear_trajectory,
// and an optional string to identify the amplitude with.
//
//  Evaluation requires the same three couplings photon coupling, gGamma, and vector/tensor
//  nucleon couplings, gV and gT respectively.
//
// Set couplings with amp.set_params({gGamma, gV, gT});
// ---------------------------------------------------------------------------

class reggeized_meson : public vector_exchange
{
public:
  reggeized_meson(reaction_kinematics * xkinem, linear_trajectory * traj, std::string exchange = "")
  : vector_exchange(xkinem, -1., exchange), alpha(traj)
    {};

  // Copy constructor
  reggeized_meson(const reggeized_meson & old)
  : vector_exchange(old), alpha(old.alpha)
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
  std::complex<double> t_channel_amplitude(std::vector<int> helicities, double s, double zs);

};

#endif
