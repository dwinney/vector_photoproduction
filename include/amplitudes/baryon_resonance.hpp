// Parameterization of a resonant amplitude in the s-channel
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _RESONANCE_
#define _RESONANCE_

#include "amplitude.hpp"

// ---------------------------------------------------------------------------
// baryon_resonance class describes the amplitude corresponding to a narrow
// (Breit-Wigner) resonance in the channel. It is parameterized in terms of
// 3 functions:
//
// 1. Breit-Wigner pole for given mass and width
// 2. Hadronic decay coupling to J/psi p final state
// 3. Photo-excitation coupling to the gamma p initial state
// ---------------------------------------------------------------------------

class baryon_resonance : public amplitude
{
private:
  int P; double J; // Spin and parity of the resonance
  double mRes, gamRes; // Resonant mass and width
  double xBR, gamPhoto; // Branching ratio and electromagnetic decay width

public:
  // Constructor
  baryon_resonance(reaction_kinematics * xkinem)
  : amplitude(xkinem)
  {};



};

#endif
