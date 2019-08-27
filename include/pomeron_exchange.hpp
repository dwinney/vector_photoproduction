// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _POMERON_
#define _POMERON_

#include "reaction_kinematics.hpp"

class pomeron_exchange
{
private:
  reaction_kinematics * kinematics;
  double norm, a0, aprime, b0;

public:
  pomeron_exchange(reaction_kinematics * xkinem)
  : kinematics(xkinem)
  {};

  double trajectory(double s)
  {
    return a0 + aprime * s;
  };

  complex<double> top_vertex(int mu, int lam_gam, int lam_vec, double s, double zs);
  complex<double> bottom_vertex(int mu, int lam_targ, int lam_rec, double s, double zs);
  complex<double> regge_factor(double s, double zs);

  complex<double> helicity_amplitude(vector<double> helicities, double s, double zs);
  complex<double> eval(double s, double zs);
};

#endif
