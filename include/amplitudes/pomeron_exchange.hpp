// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _POMERON_
#define _POMERON_

#include "amplitude.hpp"

// ---------------------------------------------------------------------------
// The pomeron_exchange class describes the amplitude correspinding to
// vector meson photoproduction via a vector pomeron coupling.
//
// As a reggeon exchange model it is parameterized in terms three functions
// 1. top_vertex() coupling the vector meson to the incoming photon
// 2. bottom_vertex() coupling the two proton dirac spinors
// 3. regge_factor() the function describing the energy dependence of the amplitude
//
// To define a pomeron_exchange object all is needed is a pointer to
// a reaction_kinematics object containing the mass and name of the vector particle
// being produced
// ---------------------------------------------------------------------------

class pomeron_exchange : public amplitude
{
private:
  double a0, aprime; // Regge trajectory intercept and slope
  double norm, b0; // Regge factor parameters: normalization and t-slope

public:
  // Constructor
  pomeron_exchange(reaction_kinematics * xkinem)
  : amplitude(xkinem)
  {};

  // Copy constructor
  pomeron_exchange(const pomeron_exchange & old)
  : amplitude(old),
    a0(old.a0), aprime(old.aprime), norm(old.norm), b0(old.b0)
  {};

  // Usual (real) linear Regge trajectory
  double trajectory(double s)
  {
    return a0 + aprime * s;
  };

  // Setting utility
  void set_params(std::vector<double> params)
  {
    norm = params[0];
    a0 = params[1];
    aprime = params[2];
    b0 = params[3];
  };

  // Photon - J/Psi - Pomeron vertex
  std::complex<double> top_vertex(int mu, int lam_gam, int lam_vec, double s, double zs);

  // Nucleon - Nucleon - Pomeron vertex
  std::complex<double> bottom_vertex(int mu, int lam_targ, int lam_rec, double s, double zs);

  // Energy dependence from Pomeron propogator
  std::complex<double> regge_factor(double s, double zs);

  // Assemble the helicity amplitude by contracting the lorentz indices
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);
};

#endif
