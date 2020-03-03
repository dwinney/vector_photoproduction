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
  int J, P; // (2xSpin) and parity of the resonance
  double mRes, gamRes; // Resonant mass and width
  double xBR; // Hadronic banching fraction to j/psi p
  double R_photo; // Photocoupling ratio

  // Initial and final CoM momenta evaluated at resonance energy.
  double pi_bar, pf_bar;

public:
  // Constructor
  baryon_resonance(int j, double mass, double width, reaction_kinematics * xkinem)
  : amplitude(xkinem), mRes(mass), gamRes(width), J(abs(j))
  {
    pi_bar = - real(kinematics->initial.momentum("target", mass * mass));
    pf_bar = - real(kinematics->final.momentum("recoil", mass * mass));

    if (j < 0){P = -1;}
    else {P = 1;}
  };

  void set_params(std::vector<double> params)
  {
    xBR = params[0];
    R_photo = params[1];
  };

  // Photoexcitation helicity amplitude for the process gamma p -> R
  double photo_coupling(int lam_i, double s, double zs);

  // Hadronic decay helicity amplitude for the R -> J/psi p process
  double hadronic_decay(int lam_i, int lam_f, double zs);

  // Ad-hoc threshold factor to kill the resonance at threshold
  double threshold_factor(double s, double beta);

  // Combined total amplitude including Breit Wigner pole
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

};

#endif
