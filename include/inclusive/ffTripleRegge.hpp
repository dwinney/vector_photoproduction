// Form of the terms inspired by the Field and Fox parameterization
//
// Nucl. Phys. B80(1974) 367
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _FIELD_FOX_
#define _FIELD_FOX_

#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

namespace jpacPhoto
{
    class ffTripleRegge
    {
        public:
        ffTripleRegge(inclusive_kinematics * xkinem, std::array<regge_trajectory*, 3> trajectories, std::array<double, 2> couplings)
        : _kinematics(xkinem), _trajectories(trajectories), _couplings(couplings)
        {};
        
        std::array<regge_trajectory*, 3> _trajectories;
        std::array<double, 2> _couplings;

        inclusive_kinematics * _kinematics;

        constexpr static double _scale = 1.0;

        inline double coupling(double t)
        {
            return _couplings[0] * exp(_couplings[1] * t);
        };
        
        inline double eval(double s, double t, double M2)
        {
            double alpha_i  = real(_trajectories[0]->eval(t));
            double alpha_j  = real(_trajectories[1]->eval(t));
            double alpha_k0 = real(_trajectories[2]->eval(0.));

            double nu = M2 - t - _kinematics->get_mT2();

            double result;
            result  = pow((s / nu), alpha_i + alpha_j); 
            result *= pow((nu / _scale), alpha_k0);
            result *= coupling(t) / (PI * _scale * s);

            return result;
        };
    };
};

#endif