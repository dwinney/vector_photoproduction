// Abstract class for a regge trajectory. Any user-defined class for a specific
// for of the RT can be used in amplitudes (e.g. pomeron_exchange) as long as it
// is derived from this class
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _REGGE_TRAJ_
#define _REGGE_TRAJ_

#include <complex>

class regge_trajectory
{
public:
  // constructor
  regge_trajectory(){};

  // copy constructor
  regge_trajectory(const regge_trajectory & old)
  {};

  // Only need a function to evaluate the trajectory at some s
  virtual std::complex<double> eval(double s) = 0;
};


// Basic linear regge_trajectory
class linear_traj : public regge_trajectory
{
private:
  // Intercept and slope
  double a0, aprime;

public:
  // Empty constructor
  linear_traj(){};

  // Parameterized constructor
  linear_traj(double inter, double slope)
  : a0(inter), aprime(slope)
  {};

  // copy Constructor
  linear_traj(const linear_traj & old)
  : a0(old.a0), aprime(old.aprime)
  {};
  // Setting utility
  void set_params(double inter, double slope)
  {
    a0 = inter; aprime = slope;
  };

  std::complex<double> eval(double s)
  {
    return a0 + aprime * s;
  }
};

#endif
