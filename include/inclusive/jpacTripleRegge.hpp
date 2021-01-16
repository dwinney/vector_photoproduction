// Form of the terms following JPAC's ds
//
// Nucl. Phys. B80(1974) 367
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TRIPLE_JPAC_
#define _TRIPLE_JPAC_

#include <complex>
#include "misc_math.hpp"
#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

namespace jpacPhoto
{
    class jpacTripleRegge
    {
        public:
        jpacTripleRegge(inclusive_kinematics * xkinem, regge_trajectory * trajectory, std::array<double,3> couplings)
        : _kinematics(xkinem), _trajectory(trajectory), _couplings(couplings)
        {};
        
        regge_trajectory * _trajectory;
        std::array<double, 3> _couplings;

        inclusive_kinematics * _kinematics;

        constexpr static double _scale = 1.0;

        inline double sigmaTOT(double s)
        {
            return _couplings[1] * pow(s, _couplings[2]);
        };

        inline std::complex<double> xi(double t)
        {
            std::complex<double> signature_factor, gamma_factor;
            signature_factor = 0.5 * (1. + double(_trajectory->_signature) * exp(XI * PI * _trajectory->eval(t)));
            gamma_factor = cgamma(1. - _trajectory->eval(t));

            return signature_factor * gamma_factor;
        };

        inline double xi_squared(double t)
        {
            return real(xi(t) * std::conj(xi(t)));
        };
        
        inline double eval(double s, double t, double M2)
        {
            double nu = M2 - t - _kinematics->get_mT2();

            double result;
            result  = sigmaTOT(s);
            result *= pow((s / M2), real(2. * _trajectory->eval(t))); 
            result *= xi_squared(t);
            result *= _couplings[0] * _couplings[0];
            result *= M2 / s;

            double normalization = real(_trajectory->slope()) / (16. * PI*PI*PI);

            return normalization * result;
        };
    };
};

#endif