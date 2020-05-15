// Axial-vector meson photoproduction proceeding through a Reggeized meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _REGGEON_
#define _REGGEON_

#include <string>
#include <vector>
#include <iostream>
#include <tuple>
#include <iomanip>
#include <algorithm>

#include "amplitude.hpp"
#include "regge_trajectory.hpp"
#include "misc_math.hpp"

class reggeon_exchange : public amplitude
{
public:
  reggeon_exchange(reaction_kinematics * xkinem, linear_trajectory * traj, std::string exchange = "")
  : amplitude(xkinem, exchange), alpha(traj)
    {};

  // Copy constructor
  reggeon_exchange(const reggeon_exchange & old)
  : amplitude(old),
    gGam(old.gGam), gV(old.gV), gT(old.gT),
    alpha(old.alpha)
  {};

  // Setting utilities
  void set_signature(int i)
  {
    signature = i;

    if (std::abs(signature) != 1)
    {
      std::cout << "\nWarning! Invalid signature (" << signature << ")";
      std::cout << " passed to reggeon_exchange: " << amplitude::identifier << ".\n";
    }
  };

  void set_params(std::vector<double> params)
  {
    gGam    = params[0];
    gV      = params[1];
    gT      = params[2];
  };

  // Assemble the helicity amplitude
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

private:
  linear_trajectory * alpha;
  int signature = -1.;

  // Couplings to the axial-vector/photon and vector/tensor couplings to nucleon
  double gGam = 0., gV = 0., gT = 0.;

  // Photon - Axial Vector - Vector vertex
  std::complex<double> top_residue(int lam, double t);

  // Nucleon - Nucleon - Vector vertex
  std::complex<double> bottom_residue(int lamp, double t);

  // Usual reggeon propagator
  std::complex<double> regge_propagator(double t);

  // Half angle factors
  std::complex<double> half_angle_factor(int lam, int lamp, std::complex<double> z_t);

  // Helicity amplitude in terms of t-channel (unrotated) helicities
  std::complex<double> t_channel_amplitude(std::vector<int> helicities, double s, double zs);

};

#endif
