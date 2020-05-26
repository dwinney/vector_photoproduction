// Spin-1/2 exchange ampltiude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PROTON_
#define _PROTON_

#include <string>
#include <vector>
#include <iostream>

#include <iomanip>

#include "amplitude.hpp"
#include "gamma_technology.hpp"

class fermion_exchange : public amplitude
{
public:
  // constructor
  fermion_exchange(reaction_kinematics * xkinem, double mass, std::string name = "")
  : amplitude(xkinem, name), mEx2(mass*mass)
  {};

  // Copy constructor
  fermion_exchange(const fermion_exchange & old)
  : amplitude(old), mEx2(old.mEx2),
    gGam(old.gGam), gVec(old.gVec)
  {};

  // Setting utility
  void set_params(std::vector<double> params)
  {
    gGam = params[0];
    gVec = params[1];
  };

  // Assemble the helicity amplitude by contracting the spinor indices
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

  // Should be exactly u_man(s, zs);
  double exchange_mass(double s, double zs);

  // debugging options to make either the photon or vector into scalars
  void set_debug(int i)
  {
    switch (i)
    {
      case 3: ScTOP = true; ScBOT = true; break;
      case 2: ScTOP = true; break;
      case 1: ScBOT = true; break;
    }
  }

private:
  // TESTING PARAMS
  bool ScTOP = false, ScBOT = false;

  // Exchange nucleon mass
  double mEx2;

  // couplings
  double gGam = 0., gVec = 0.;

  // Four-momentum of the exhange (u - channel)
  std::complex<double> exchange_momentum(int mu, double s, double zs);

  // Slashed momentumn
  std::complex<double> slashed_exchange_momentum(int i, int j, double s, double zs);

  // Slashed polarization vectors
  std::complex<double> slashed_eps(int i, int j, double lam, polarization_vector eps, bool STARRED, double s, double zs);

  // Photon - excNucleon - recNucleon vertex
  std::complex<double> top_vertex(int i, int lam_gam, int lam_rec, double s, double zs);

  // excNucleon - recNucleon - Vector vertex
  std::complex<double> bottom_vertex(int j, int lam_vec, int lam_targ, double s, double zs);

  // Spin-1/2 propagator
  std::complex<double> fermion_propagator(int i, int j, double s, double zs);
};

#endif
