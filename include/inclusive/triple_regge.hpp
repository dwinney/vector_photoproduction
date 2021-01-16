// Form of the invariant cross-section from a triple regge interaction.
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TRIP_REGGE_
#define _TRIP_REGGE_

#include "inclusive_kinematics.hpp"
#include "ffTripleRegge.hpp"
#include "jpacTripleRegge.hpp"
#include "regge_trajectory.hpp"

#include <vector>

namespace jpacPhoto
{
    class triple_regge
    {
        public:
        // Constructor only needs a kinematics object
        triple_regge(inclusive_kinematics * xkinem, std::string id = "")
        : _kinematics(xkinem), _identifier(id)
        {};

        // Destructor to clean up pointers
        ~triple_regge()
        {
            for (int i = 0; i < _termsFF.size(); i++)
            {
                delete _termsFF[i];
            }
            for (int i = 0; i < _termsJPAC.size(); i++)
            {
                delete _termsJPAC[i];
            }
        };

        inline void add_term(std::array<regge_trajectory*, 3> trajectories, std::array<double,2> couplings)
        {
            auto new_term = new ffTripleRegge(_kinematics, trajectories, couplings);
            _termsFF.push_back(new_term);
        };

        inline void add_term(regge_trajectory* trajectory, std::array<double,3> couplings)
        {
            auto new_term = new jpacTripleRegge(_kinematics, trajectory, couplings);
            _termsJPAC.push_back(new_term);
        };

        // E d sigma / d^3p
        double invariant_xsection(double s, double t, double M2);

        std::string _identifier;
        
        protected:
        inclusive_kinematics * _kinematics;
        
        // Cross-section build out of sum of triple regge interaction terms
        std::vector<ffTripleRegge*> _termsFF;
        std::vector<jpacTripleRegge*> _termsJPAC;
    };
};

#endif 