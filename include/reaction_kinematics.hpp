// Class to contain all relevant kinematic quantities. The kinematics of the reaction
// gamma p -> X p' is entirely determined by specifying the mass of the vector particle.
//
// Additional options to include virtual photon and different baryons (e.g. gamma p -> X Lambda_c) 
// also available
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------


#ifndef _KINEMATICS_
#define _KINEMATICS_

#include "constants.hpp"
#include "misc_math.hpp"
#include "two_body_state.hpp"
#include "dirac_spinor.hpp"
#include "polarization_vector.hpp"
#include "helicities.hpp"

#include "TMath.h"

#include <array>
#include <vector>
#include <string>
#include <cmath>

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // The reaction kinematics object is intended to have all relevant kinematic quantities
    // forthe reaction. Here you'll find the momenta and energies of all particles,
    //  spinors for the baryons and polarization vectors for the gamma and produced meson
    // ---------------------------------------------------------------------------

    class reaction_kinematics
    {
        public: 
        // Empty constructor,
        // defaults to compton scattering: gamma p -> gamma p
        reaction_kinematics()
        {
            initial   = new two_body_state(0., mPro2);
            eps_gamma = new polarization_vector(initial);
            target    = new dirac_spinor(initial);

            final     = new two_body_state(0., mPro2);
            eps_vec   = new polarization_vector(final);
            recoil    = new dirac_spinor(final);
        };

        // Constructor with a set mX and JP
        // defaults to proton as baryon and real photon
        // string ID is deprecated but kept for backward compatibility
        reaction_kinematics(double _mX, std::string id = "")
        : mX(_mX), mX2(_mX*_mX)
        {
            initial   = new two_body_state(0., mPro2);
            eps_gamma = new polarization_vector(initial);
            target    = new dirac_spinor(initial);

            final     = new two_body_state(_mX*_mX, mPro2);
            eps_vec   = new polarization_vector(final);
            recoil    = new dirac_spinor(final);
        };


        // Constructor with a set mX and baryon mass mR
        // defaults to real photon
        reaction_kinematics(double _mX, double _mR)
        : mX(_mX), mX2(_mX*_mX),
          mR(_mR), mR2(_mR*_mR)
        {
            initial   = new two_body_state(0., mPro2);
            eps_gamma = new polarization_vector(initial);
            target    = new dirac_spinor(initial);

            final     = new two_body_state(_mX*_mX, _mR*_mR);
            eps_vec   = new polarization_vector(final);
            recoil    = new dirac_spinor(final);
        };

        // Constructor with a set mV and baryon mass mR
        // and massive incoming scatterer and target
        reaction_kinematics(double _mX, double _mR, double _mT, double _mB = 0.)
        : mX(_mX), mX2(_mX*_mX),
          mR(_mR), mR2(_mR*_mR),
          mB(_mB), mB2(_mB*_mB),
          mT(_mT), mT2(_mT*_mT)
        {
            initial   = new two_body_state(_mB*_mB, _mT*_mT);
            eps_gamma = new polarization_vector(initial);
            target    = new dirac_spinor(initial);

            final     = new two_body_state(_mX*_mX, _mR*_mR);
            eps_vec   = new polarization_vector(final);
            recoil    = new dirac_spinor(final);
        };

        // destructor
        ~reaction_kinematics()
        {
            delete initial, final;
            delete eps_gamma, eps_vec;
            delete target, recoil;
        }

        // ---------------------------------------------------------------------------
        // Masses
        
        double mB = 0., mB2 = 0.;       // mass and mass squared of the "beam" 
        double mX = 0., mX2 = 0.;       // mass and mass squared of the produced particle

        double mT = mPro, mT2 = mPro2;  // mass of the target, assumed to be proton unless overriden
        double mR = mPro, mR2 = mPro2;  // mass of the recoil baryon, assumed to be proton unless overriden

        inline double Wth(){ return (mX + mR); }; // square root of the threshold
        inline double sth(){ return Wth() * Wth(); }; // final state threshold

        // Change the meson mass
        inline void set_mX(double m)
        {
            mX  = m;
            mX2 = m*m;

            // also update the meson mass in two_body_state
            final->set_mV2(m*m);
        };

        inline void set_mX2(double m2)
        {
            mX  = sqrt(m2);
            mX2 = m2;

            // also update the meson mass in two_body_state
            final->set_mV2(m2);
        };

        // Change virtuality of the photon
        // Q2 > 0
        inline void set_Q2(double q2)
        {
            if (q2 > 0) { std::cout << "Caution! set_Q2(x) requires x > 0! \n"; }
            mB2 = -q2;
            initial->set_mV2(-q2);
        };

        // ---------------------------------------------------------------------------
        // Quantum numbers of produced meson. 
        std::array<int,2> JP{{1,1}};
        inline void set_JP(int _J, int _P)
        { 
            JP = {_J, _P};
            helicities = get_helicities(_J);
            nAmps = helicities.size();
        };

        // Helicity configurations
        // Photon [0], Incoming Proton [1], Produced meson [2], Outgoing Proton [3]
        int nAmps = 24; 
        std::vector< std::array<int, 4> > helicities = spin_one_helicities;

        //--------------------------------------------------------------------------
        two_body_state * initial,  * final;
        polarization_vector * eps_vec, * eps_gamma;
        dirac_spinor * target, * recoil;

        // Get s-channel scattering angle from invariants
        inline double z_s(double s, double t)
        {
            std::complex<double> qdotqp = initial->momentum(s) * final->momentum(s);
            std::complex<double> E1E3   = initial->energy_V(s) * final->energy_V(s);

            double result = t - mX2 - mB2 + 2.*abs(E1E3);
            result /= 2. * abs(qdotqp);

            return result;
        };

        // Scattering angle in the s-channel
        // Use TMath::ACos instead of std::acos because its safer at the end points
        inline double theta_s(double s, double t)
        {
            return TMath::ACos( z_s(s, t) );
        };

        // Invariant variables
        inline double t_man(double s, double theta)
        {
            std::complex<double> qdotqp = initial->momentum(s) * final->momentum(s);
            std::complex<double> E1E3   = initial->energy_V(s) * final->energy_V(s);

            return mX2 + mB2 - 2. * abs(E1E3) + 2. * abs(qdotqp) * cos(theta);
        };

        inline double u_man(double s, double theta)
        { 
            return mX2 + mB2 + mT2 + mR2 - s - t_man(s, theta);
        };

        // Scattering angles in t and u channel frames
        inline std::complex<double> z_t(double s, double theta)
        {
            double t = t_man(s, theta);
            std::complex<double> p_t = sqrt(xr * Kallen(t, mT2, mR2)) / sqrt(xr * 4. * t);
            std::complex<double> q_t = sqrt(xr * Kallen(t, mX2, mB2)) / sqrt(xr * 4. * t);

            std::complex<double> result;
            result = 2. * s + t - mT2 - mR2 - mX2 - mB2; // s - u
            result /= 4. * p_t * q_t;

            return result;
        };

        inline std::complex<double> z_u(double s, double theta)
        {
            // TODO: fix this
            std::cout << "z_u not implimented yet i fucked up here :p\n";
            std::cout << "if youre seeing this message email me and yell at me because i forgot\n";
            
            return xr;
        };
    };
};

#endif
