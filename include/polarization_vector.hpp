// Class for the polarization vector of vector particles
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef POLVEC
#define POLVEC

#include <iostream>
#include "constants.hpp"
#include "two_body_state.hpp"

// ---------------------------------------------------------------------------
// Polarization vectors for vector particles
// in the s-channel center of mass frame
// Vector particles are always particle 1
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
    class polarization_vector
    {
        public:

        // Constructor
        polarization_vector(two_body_state * xstate)
            : _state(xstate)
        {};

        // Destructor
        ~polarization_vector(){};

        // Components
        std::complex<double> component(int i, int lambda, double s, double theta);
        std::complex<double> conjugate_component(int i, int lambda, double s, double theta);

        private:
        
        two_body_state * _state;
    };
};

#endif
