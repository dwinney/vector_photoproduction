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

#include <iomanip>

#include "amplitude.hpp"
#include "misc_math.hpp"

class reggeon_exchange : public amplitude
{
public:
  // For testing have a fixed mass, eventually this will require a trajectory object
  reggeon_exchange(reaction_kinematics * xkinem, double mass, std::string exchange)
  : amplitude(xkinem, exchange), mEx2(mass*mass)
  {};

  // Copy constructor
  reggeon_exchange(const reggeon_exchange & old)
  : amplitude(old), mEx2(old.mEx2),
    gGamma(old.gGamma), gV(old.gV), gT(old.gT)
  {};

  // Setting utility
  void set_params(std::vector<double> params)
  {
    gGamma = params[0];
    gV = params[1];
    gT = params[2];
  };

  // Assemble the helicity amplitude by contracting the lorentz indices
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

private:
  // Mass of the exchange
  double mEx2;

  // Couplings to the axial-vector/photon and vector/tensor couplings to nucleon
  double gGamma = 0., gV = 0., gT = 0.;

  // Photon - Axial Vector - Vector vertex
  std::complex<double> top_vertex(int lam, double t);

  // Nucleon - Nucleon - Vector vertex
  std::complex<double> bottom_vertex(int lamp,  double t);

  // Helicity amplitude in terms of t-channel (unrotated) helicities
  std::complex<double> t_channel_amplitude(std::vector<int> helicities, double s, double zs);

};

#endif
