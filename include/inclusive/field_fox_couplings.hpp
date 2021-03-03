// Phenomenological expressions for triple Regge couplings as extracted by Field and Fox
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef FF_COUPS
#define FF_COUPS

#include <cmath>

inline double G_PPP(double t)
{
    return 2.31 * exp(3.94 * t) + 0.33 * exp(1.12 * t);
};

inline double G_RRP(double t)
{
    return 26.81 * exp(7.26 * t) + 4.80 * exp(-1.83 * t);
};

inline double G_RRR(double t)
{
    return 18.1 * exp(12. * t);
};

inline double G_PPR(double t)
{
    return 0.95 * exp(-0.01 * t) + 3.47 * exp(4.41 * t);
};

#endif