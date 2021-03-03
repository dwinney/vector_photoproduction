// Form of the terms following JPAC's normalization
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef TRIPLE_JPAC
#define TRIPLE_JPAC

#include <complex>
#include <tuple>
#include <functional>
#include "misc_math.hpp"
#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

namespace jpacPhoto
{
    class jpacTripleRegge
    {
        public:

        // Fully general constructor
        jpacTripleRegge(inclusive_kinematics * xkinem, regge_trajectory * trajectory, 
                        const std::function<double(double)>& coupling, const std::function<double(double)>& sigmatot)
        : _kinematics(xkinem), _trajectory(trajectory), _coupling(coupling), _sigmatot(sigmatot)
        {};
        
        inline double eval(double s, double t, double M2)
        {
            double nu = M2 - t - _kinematics->get_mT2();

            double result;
            result  = _sigmatot(s);
            result *= pow((s / M2), real(2. * _trajectory->eval(t)) - 1.); 
            result *= norm(xi(t));
            result *= _coupling(t) * _coupling(t);
            
            double normalization = real(_trajectory->slope()) / (16. * PI*PI*PI);

            return normalization * result;
        };

        protected: 
        
        inclusive_kinematics * _kinematics;
        regge_trajectory * _trajectory;

        std::function<double(double)> _coupling;
        std::function<double(double)> _sigmatot;

        constexpr static double _scale = 1.0;

        inline std::complex<double> xi(double t)
        {
            if (std::abs(t) > 40.) return 0.;

            std::complex<double> signature_factor, gamma_factor;
            signature_factor = 0.5 * (1. + double(_trajectory->_signature) * exp(XI * PI * _trajectory->eval(t)));
            gamma_factor = cgamma(double(_trajectory->_minJ) - _trajectory->eval(t));

            return signature_factor * gamma_factor;
        };
    };
};

#endif