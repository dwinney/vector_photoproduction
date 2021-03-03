// Energy and momentum of a two-particle state in the center of mass scattering frame
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef TWO_BODY
#define TWO_BODY

#include <string>
#include <complex>
#include <iostream>

#include "misc_math.hpp"

// ---------------------------------------------------------------------------
// The two_body_state is the base object for defining a reaction in the
// s-channel center of mass scatering frame.
//
// Two particles of mass mV and mB are defined with momenta opposite along the
// same axis such that the energy and momenta of both particles is entirely
// determined by the center-of-mass energy, s, and the cosing of the angle
// from the z-axis (define to be at theta = 0).
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class two_body_state
  {
        private:

        double _mV2; // Vector mass
        double _mB2; // Baryon mass (allowed to be float for N* or Δ)

        public:

        // Constructor
        two_body_state(double mV2, double mB2)
        : _mV2(mV2), _mB2(mB2)
        {};

        // return mass
        inline double get_mV() 
        { 
            if (_mV2 >= 0.) 
            {
                return sqrt(_mV2);
            }
            else
            {
                return sqrt(-_mV2);
            }
        };

        inline double get_mB() { return sqrt(_mB2); };

        // Return mass squared
        inline double get_mV2() { return _mV2; };
        inline double get_mB2() { return _mB2; };

        // set masses independently
        inline void set_mV2(double mV2)
        {
            _mV2 = mV2;
        };

        inline void set_mB2(double mB2)
        {
            _mB2 = mB2;
        };

        // Momenta
        // V is always particle 1 in + z direction, 
        inline std::complex<double> momentum(double s)
        {
            return sqrt( Kallen(XR * s, XR *_mV2, XR * _mB2)) / (2. * sqrt(XR * s));
        };

        // Energies
        inline std::complex<double> energy_V(double s)
        {
            return (s + _mV2 - _mB2) / (2. * sqrt(XR * s));
        };

        inline std::complex<double> energy_B(double s)
        {
            return (s - _mV2 + _mB2) / (2. * sqrt(XR * s));
        };

        // Full 4-momenta 
        std::complex<double> q(int mu, double s, double theta); // 4vector of vector, particle 1
        std::complex<double> p(int mu, double s, double theta); // 4vector of baryon, particle 2
    };
};

#endif
