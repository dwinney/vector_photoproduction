// Misc math functions that are useful to have
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _MATH_
#define _MATH_

#include "constants.hpp"
#include <iostream>

template <typename T>
T Kallen(T x, T y, T z)
{
  return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
};

// Wigner d-function for half-integer spin
void wigner_error(int j, int lam1, int lam2);
double wigner_d(int j, int lam1, int lam2, double z);


#endif
