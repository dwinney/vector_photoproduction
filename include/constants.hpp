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

// ---------------------------------------------------------------------------
const double deg2rad = (M_PI / 180.);
const double EPS = 1.e-6;
const double M_ALPHA = 1. / 137.;
const std::complex<double> xr(1., 0.);
const std::complex<double> xi(0., 1.);

// Masses
const double mPi = 3.14159265359;
const double mPro = 0.9383;
const double mJpsi = 3.097;
const double mPsi2S = 3.686;

// Masses squared
const double mPi_sqr = mPi * mPi;
const double mPro_sqr = mPro * mPro;
const double mJpsi_sqr = mJpsi * mJpsi;

// Thresholds
const double sthPsiPro = (mJpsi + mPro) * (mJpsi + mPro);

// Decay constants
const double fJpsi = 0.278;
// ---------------------------------------------------------------------------
const std::complex<double> gamma_matrices[4][4][4] =
{
//gamma0
	{ { 1., 0., 0., 0. },
    { 0., 1., 0., 0. },
    { 0., 0., -1., 0. },
    { 0., 0., 0., -1. } },
//gamma1
	{ { 0., 0., 0., 1. },
    { 0., 0., 1., 0. },
    { 0., -1., 0., 0. },
    { -1., 0., 0., 0. } },
//gamma2
	{ { 0., 0., 0., -xi },
    { 0., 0., xi, 0. },
    { 0.,  xi, 0., 0. },
    { -xi, 0., 0., 0. } },
//gamma3
	{ { 0., 0., 1., 0. },
    { 0., 0., 0., -1. },
    { -1., 0., 0., 0. },
    { 0., 1., 0., 0. } }
};

const std::complex<double> gamma_5[4][4] =
{
  { 0., 0., 1., 0. },
  { 0., 0., 0., 1. },
  { 1., 0., 0., 0. },
  { 0., 1., 0., 0. }
};

const double metric[4] = {1., -1., -1., -1.};

#endif
