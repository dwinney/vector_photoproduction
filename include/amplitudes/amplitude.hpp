// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AMPLITUDE_
#define _AMPLITUDE_

// ---------------------------------------------------------------------------
// Abstract class to define helicity amplitudes. This will allow multiple different
// classes (for s, t, and u- channels but also multiple contibutions in each channel)
// to be added together and evaluated in observables.
//
// Any generic amplitude needs a reaction_kinematics object
// and a way to evaluate the helicity amplitude for given set of helicities,
// CoM energy and scattering angle.
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

#include <string>

class amplitude
{
public:
  // Constructor
  amplitude(reaction_kinematics * xkinem)
  : kinematics(xkinem)
  {};

  amplitude(reaction_kinematics * xkinem, std::string id)
  : kinematics(xkinem), identifier(id)
  {};

  // Copy constructor
  amplitude(const amplitude & old)
  : kinematics(old.kinematics), identifier(old.identifier)
  {};

  // Kinematics object for thresholds and etc.
  reaction_kinematics * kinematics;

  // Some saveable string by which to identify the amplitude
  std::string identifier;

  // How the calculate the helicity amplitude
  virtual std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs) = 0;

  // ---------------------------------------------------------------------------
  // Observables

  // Differential and total cross-section
  double differential_xsection(double s, double zs);
  double integrated_xsection(double s);

  // Spin asymmetries
  double K_LL(double s, double zs);
  double A_LL(double s, double zs);

  // Spin density matrix elements
  std::complex<double> SDME(int alpha, int lam, int lamp, double s, double zs);

  // Asymmetries
  double beam_asymmetry(double s, double zs);
  double parity_asymmetry(double s, double zs);

  // ---------------------------------------------------------------------------
  // Nparams error message
  void check_Nparams(int n, std::vector<double> params)
  {
    if (n != params.size())
    {
      std::cout << "\nWarning! Invalid number of parameters (" << params.size() << ") passed to " << identifier << ".\n";
    }
  }

};

#endif
