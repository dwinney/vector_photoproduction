// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef CONSTANT
#define CONSTANT

#include <cmath>
#include <complex>

using std::complex;
using std::abs;

// ---------------------------------------------------------------------------
const double deg2rad = (M_PI / 180.);
const double EPS = 1.e-6;
const complex<double> xr(1., 0.);
const complex<double> xi(0., 1.);

// Masses
const double mPi = 3.14159265359;
const double mPro = 0.9383;
const double mJpsi = 3.097;

// Masses squared
const double mPi_sqr = mPi * mPi;
const double mPro_sqr = mPro * mPro;
const double mJpsi_sqr = mJpsi * mJpsi;

// Thresholds
const double sthPsiPro = (mJpsi + mPro) * (mJpsi + mPro);

// Decay constants
const double fJpsi = 0.278;
// ---------------------------------------------------------------------------

#endif
