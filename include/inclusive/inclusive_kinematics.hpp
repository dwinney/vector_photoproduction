// Extension of the reaction_kinematics class to include quantities revevant
// for semi-inclusive reactions at high energies.
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef INC_KINEM
#define INC_KINEM

#include "constants.hpp"
#include "misc_math.hpp"

namespace jpacPhoto
{
    class inclusive_kinematics 
    {
        public:
        // Empty constructor
        inclusive_kinematics()
        {};

        // Construtor with produced meson mass
        inclusive_kinematics(double mX)
        : _mX(mX), _mX2(mX * mX)
        {};

        inline void set_minM2(double xM2)
        {
            _minM2 = xM2;
        };

        inline double get_mT2()
        {
            return _mT2;
        };

        double _mX, _mX2;                         // Mass of the produce (observed particle)
        double _mT = M_PROTON, _mT2 = M2_PROTON;  // Mass of the target

        double _minM2 = M2_PROTON; // Default to proton is minimum mass unobserved

        inline double M2(double s, double x)
        {
            return _minM2 + (s - _minM2) * (1. - x);
        };

        // ---------------------------------------------------------------------------
        // Center-of-mass kinematics

        inline double cosTheta_CM(double xs, double xt, double xM2)
        {
            double u = _mX2 + _mT2 + xM2 - xs - xt;

            double result;
            result  = xs * (xt - u) - _mT2 * (_mX2 - xM2);
            result /= sqrt(Kallen(xs, 0., _mT2) * Kallen(xs, _mX2, xM2));
            
            return result;
        };

        inline double t_man(double s, double costheta, double M2)
        {
            double result;
            result  = 2. * pGamma_CM(s) * pX_CM(s, M2) * costheta;
            result -= (s * (s - _mT2 - _mX2 - M2) - _mT2 * (_mX2 - M2)) / (2. * s);
            
            return result;
        };

        inline double pX_CM(double xs, double xM2)
        {
            return sqrt(Kallen(xs, _mX2, xM2)) / (2. * sqrt(xs));
        };

        inline double pGamma_CM(double xs)
        {
            return sqrt(Kallen(xs, 0., M2_PROTON)) / (2. * sqrt(xs));
        };

        inline double pPerp_CM(double xs, double xt, double xM2)
        {
            return pX_CM(xs, xM2) * cosTheta_CM(xs, xt, xM2);
        };

        inline double x_CM(double xs, double xt, double xM2)
        {
            return pPerp_CM(xs, xt, xM2) / pX_CM(xs, _minM2);
        };
    };
};

#endif