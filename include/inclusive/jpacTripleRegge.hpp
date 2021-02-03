// Form of the terms following JPAC's normalization
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TRIPLE_JPAC_
#define _TRIPLE_JPAC_

#include <complex>
#include <tuple>
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
                        std::vector<std::tuple<int,double>> couplings, std::vector<std::array<double, 2>> sigmatot)
        : _kinematics(xkinem), _trajectory(trajectory), _couplings(couplings), _sigmaParams(sigmatot)
        {};
        
        inline double eval(double s, double t, double M2)
        {
            double nu = M2 - t - _kinematics->get_mT2();

            double result;
            result  = sigmaTOT(s);
            result *= pow((s / M2), real(2. * _trajectory->eval(t) - 1.)); 
            result *= norm(xi(t));
            result *= norm(coupling(t));

            double normalization = real(_trajectory->slope()) / (16. * PI*PI*PI);

            return normalization * result;
        };

        protected: 
        
        regge_trajectory * _trajectory;

        std::vector<std::tuple<int,double>> _couplings;
        std::vector<std::array<double, 2>> _sigmaParams;

        inclusive_kinematics * _kinematics;

        constexpr static double _scale = 1.0;

        inline double sigmaTOT(double s)
        {
            double result = 0.;
            for (int i = 0; i < _sigmaParams.size(); i++)
            {
                result += _sigmaParams[i][0] * pow(s, _sigmaParams[i][1]);
            }
            return result;
        };

        inline std::complex<double> xi(double t)
        {
            std::complex<double> signature_factor, gamma_factor;
            signature_factor = 0.5 * (1. + double(_trajectory->_signature) * exp(XI * PI * _trajectory->eval(t)));
            gamma_factor = cgamma(double(_trajectory->_minJ) - _trajectory->eval(t));

            return signature_factor * gamma_factor;
        };
        
        inline std::complex<double> coupling(double t)
        {
            std::complex<double> result = 0.; 
            for (int i = 0; i < _couplings.size(); i++)
            {
                double power = (double) std::get<0>(_couplings[i]);
                double beta  = std::get<1>(_couplings[i]);

                if (power == 0) result += beta;
                else if (std::abs(t) < 1.E-4 && power != 0) return 0.;
                else result += beta * pow(sqrt(XR * -t), power);
            };

            return result;
        };
    };
};

#endif