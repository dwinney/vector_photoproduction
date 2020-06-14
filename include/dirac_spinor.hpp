// Class for dirac spinors for the spin-1/2 nucleon
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _SPINOR_
#define _SPINOR_

#include <cmath>
#include <complex>

#include "constants.hpp"
#include "gamma_technology.hpp"
#include "two_body_state.hpp"

// ---------------------------------------------------------------------------
// The dirac_spinor object is exactly as it sounds. :)
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class dirac_spinor
  {
  private:
    two_body_state twobody;
    const std::string particle;
    const double mass;

    //Whether its an anti-particle or not
    const bool ANTI_PARTICLE = false;

    double cos_half(double zs);
    double sin_half(double zs);
    std::complex<double> momentum(int sign, double s);

  public:
    // Constructor
    dirac_spinor(two_body_state xstate, std::string name, bool if_anti = false)
    : twobody(xstate), particle(name), mass(xstate.get_mass(name)),
      ANTI_PARTICLE(if_anti)
    {};

  	// Copy Constructor
  	dirac_spinor(const dirac_spinor & old)
  	: twobody(old.twobody), particle(old.particle), mass(old.mass),
  		ANTI_PARTICLE(old.ANTI_PARTICLE)
  	{};

    // Default destructor
    ~dirac_spinor(){};

    // Components
    std::complex<double> component(int i, int lambda, double s, double zs);
  	std::complex<double> adjoint_component(int i, int lambda, double s, double zs);
  };
};

#endif
