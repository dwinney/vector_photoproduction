// Discontinuity across unitarity cut for vector production via a one-loop box diagram
// Used as a container class in integration processes
// 
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BOX_DISC_
#define _BOX_DISC_

#include "constants.hpp"
#include "amplitudes/amplitude.hpp"
#include "amplitudes/reaction_kinematics.hpp"

#include "cubature.h"

namespace jpacPhoto
{
    // Wrapper class for integration with cubature
    struct disc_integrand
    {
        disc_integrand(amplitude * left, amplitude * right)
        : _initialAmp(left), _finalAmp(right)
        {};

        // External helicities and center-of-mass energy / angle
        inline void set_externals(std::array<int,4> helicities, double s, double theta)
        {
            _lam_gam = helicities[0];
            _lam_tar = helicities[1];
            _lam_vec = helicities[2];
            _lam_rec = helicities[3];

            _s = s; _theta = theta;
        };
        double _s, _theta;
        int _lam_gam, _lam_tar, _lam_vec, _lam_rec;

        // Evaluate the product of subamplitudes
        std::complex<double> eval(double theta_gam, double phi_gam);
        
        // Save the pointers to the two sub-amplitudes
        amplitude * _initialAmp;
        amplitude * _finalAmp;
    };

    // Actual discontinuity
    class box_discontinuity
    {
        public: 
        box_discontinuity(amplitude * left, amplitude * right)
        {
            _integrand = new disc_integrand(left, right);
        };

        ~box_discontinuity()
        {
            delete _integrand;
        };

        // Evaluate the discontinuity integrated over intermediate phase space
        std::complex<double> eval(double s);

        inline void set_externals(std::array<int,4> helicities, double s, double theta)
        {
            _external_s = s;
            _external_theta = theta; 
            _external_helicities = helicities;
        };

        inline double s(){return _external_s;};

        private:
        double _external_s, _external_theta;
        std::array<int,4> _external_helicities;

        // Wrapper class for integration with cubature
        disc_integrand * _integrand;
        static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
    };
};

#endif