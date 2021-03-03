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


#ifndef KINEMATICS
#define KINEMATICS

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
            _initial_state   = new two_body_state(0., M2_PROTON);
            _eps_gamma       = new polarization_vector(_initial_state);
            _target          = new dirac_spinor(_initial_state);

            _final_state     = new two_body_state(0., M2_PROTON);
            _eps_vec         = new polarization_vector(_final_state);
            _recoil          = new dirac_spinor(_final_state);
        };

        // Constructor with a set mX and JP
        // defaults to proton as baryon and real photon
        // string ID is deprecated but kept for backward compatibility
        reaction_kinematics(double mX, std::string id = "")
        : _mX(mX), _mX2(mX*mX)
        {
            _initial_state   = new two_body_state(0., M2_PROTON);
            _eps_gamma       = new polarization_vector(_initial_state);
            _target          = new dirac_spinor(_initial_state);

            _final_state     = new two_body_state(mX*mX, M2_PROTON);
            _eps_vec         = new polarization_vector(_final_state);
            _recoil          = new dirac_spinor(_final_state);
        };


        // Constructor with a set mX and baryon mass mR
        // defaults to real photon
        reaction_kinematics(double mX, double mR)
        : _mX(mX), _mX2(mX*mX),
          _mR(mR), _mR2(mR*mR)
        {
            _initial_state   = new two_body_state(0., M2_PROTON);
            _eps_gamma       = new polarization_vector(_initial_state);
            _target          = new dirac_spinor(_initial_state);

            _final_state     = new two_body_state(mX*mX, mR*mR);
            _eps_vec         = new polarization_vector(_final_state);
            _recoil          = new dirac_spinor(_final_state);
        };

        // Constructor with a set mV and baryon mass mR
        // and massive incoming scatterer and target
        reaction_kinematics(double mX, double mR, double mT, double mB = 0.)
        : _mX(mX), _mX2(mX*mX),
          _mR(mR), _mR2(mR*mR),
          _mB(mB), _mB2(mB*mB),
          _mT(mT), _mT2(mT*mT)
        {
            _initial_state   = new two_body_state(mB*mB, mT*mT);
            _eps_gamma       = new polarization_vector(_initial_state);
            _target          = new dirac_spinor(_initial_state);

            _final_state     = new two_body_state(mX*mX, mR*mR);
            _eps_vec         = new polarization_vector(_final_state);
            _recoil          = new dirac_spinor(_final_state);
        };

        // destructor
        ~reaction_kinematics()
        {
            delete _initial_state;
            delete _final_state;
            delete _eps_gamma;
            delete _eps_vec;
            delete _target;
            delete _recoil;
        }

        // ---------------------------------------------------------------------------
        // Masses
        
        double _mB = 0., _mB2 = 0.;       // mass and mass squared of the "beam" 
        double _mX = 0., _mX2 = 0.;       // mass and mass squared of the produced particle

        double _mT = M_PROTON, _mT2 = M2_PROTON;  // mass of the target, assumed to be proton unless overriden
        double _mR = M_PROTON, _mR2 = M2_PROTON;  // mass of the recoil baryon, assumed to be proton unless overriden

        inline double Wth(){ return (_mX + _mR); }; // square root of the threshold
        inline double sth(){ return Wth() * Wth(); }; // final state threshold

        // Change the meson mass
        inline void set_mX(double m)
        {
            _mX  = m;
            _mX2 = m*m;

            // also update the meson mass in two_body_state
            _final_state->set_mV2(m*m);
        };

        inline void set_mX2(double m2)
        {
            _mX  = sqrt(m2);
            _mX2 = m2;

            // also update the meson mass in two_body_state
            _final_state->set_mV2(m2);
        };

        // Change virtuality of the photon
        // Q2 > 0
        inline void set_Q2(double q2)
        {
            if (q2 < 0) { std::cout << "Caution! set_Q2(x) requires x > 0! \n"; }
            _mB2 = -q2;
            _initial_state->set_mV2(-q2);
        };

        // ---------------------------------------------------------------------------
        // Quantum numbers of produced meson. 
        std::array<int,2> _jp{{1,1}};
        inline void set_JP(int J, int P)
        { 
            _jp = {J, P};
            _helicities = get_helicities(J);
            _nAmps = _helicities.size();
        };

        // Helicity configurations
        // Defaults to spin-1
        // Photon [0], Incoming Proton [1], Produced meson [2], Outgoing Proton [3]
        int _nAmps = 24; 
        std::vector< std::array<int, 4> > _helicities = SPIN_ONE_HELICITIES;

        //--------------------------------------------------------------------------
        two_body_state * _initial_state,  * _final_state;
        polarization_vector * _eps_vec, * _eps_gamma;
        dirac_spinor * _target, * _recoil;

        // Get s-channel scattering angle from invariants
        inline double z_s(double s, double t)
        {
            std::complex<double> qdotqp = _initial_state->momentum(s) * _final_state->momentum(s);
            std::complex<double> E1E3   = _initial_state->energy_V(s) * _final_state->energy_V(s);

            double result = t - _mX2 - _mB2 + 2.*abs(E1E3);
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
            std::complex<double> qdotqp = _initial_state->momentum(s) * _final_state->momentum(s);
            std::complex<double> E1E3   = _initial_state->energy_V(s) * _final_state->energy_V(s);

            return _mX2 + _mB2 - 2. * abs(E1E3) + 2. * abs(qdotqp) * cos(theta);
        };

        inline double u_man(double s, double theta)
        { 
            return _mX2 + _mB2 + _mT2 + _mR2 - s - t_man(s, theta);
        };

        // Scattering angles in t and u channel frames
        inline std::complex<double> z_t(double s, double theta)
        {
            double t = t_man(s, theta);
            std::complex<double> p_t = sqrt(XR * Kallen(t, _mT2, _mR2)) / sqrt(XR * 4. * t);
            std::complex<double> q_t = sqrt(XR * Kallen(t, _mX2, _mB2)) / sqrt(XR * 4. * t);

            std::complex<double> result;
            result = 2. * s + t - _mT2 - _mR2 - _mX2 - _mB2; // s - u
            result /= 4. * p_t * q_t;

            return result;
        };
    };
};

#endif
