// Vector production via a one-loop box diagram
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BOX_AMP_
#define _BOX_AMP_

#include "constants.hpp"
#include "amplitudes/amplitude.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "one_loop/box_discontinuity.hpp"

#include "cubature.h"

namespace jpacPhoto
{
    class box_amplitude : public amplitude
    {
        public: 
        // Constructor.
        // Need the parent reaction kinematics and pre-set sub-amplitudes
        box_amplitude(reaction_kinematics * xkinem, amplitude * left, amplitude * right, std::string id = "Box Amplitude")
        : amplitude(xkinem, id)
        {
            _s_thr = left->_kinematics->sth();
            _discontinuity = new box_discontinuity(left, right);
        };

        // Destructor, delete the only new pointer
        ~box_amplitude()
        {
            delete _discontinuity;
        };

        // Setter for max cutoff in dispersion relation
        inline void set_cutoff(double s_cut)
        {
            _s_cut = s_cut;
        };

        // only vector available
        inline std::vector<std::array<int,2>> allowedJP()
        {
            return { {1, -1} };
        };

        // Evaluate the helicity amplitude by dispersing
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        private:

        // The two sub-amplitudes that contribute to the discontinuity 
        box_discontinuity * _discontinuity;
        static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
    
        // External helicites
        std::array<int,4> _helicities;
        
        // Integration momentum cutoff. Defaults to 2 GeV (an arbitrary but sensible value)
        double _s_cut = 2.;
        double _s_thr; 
    };
};

#endif