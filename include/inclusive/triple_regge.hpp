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
#include <tuple>

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

        //--------------------------------------------------------------------
        // Methods to add field and fox like terms

        // If only one term in the coupling
        inline void add_term(std::array<regge_trajectory*, 3> trajectories, std::array<double,2> couplings)
        {
            auto new_term = new ffTripleRegge(_kinematics, trajectories, {couplings});
            _termsFF.push_back(new_term);
        };

        // If many (passing a vector)
        inline void add_term(std::array<regge_trajectory*, 3> trajectories,std::vector<std::array<double,2>> couplings)
        {
            auto new_term = new ffTripleRegge(_kinematics, trajectories, couplings);
            _termsFF.push_back(new_term);
        };
        
        //--------------------------------------------------------------------
        // Methods to add terms following Vincent's normalization

        // Only one coupling and one sigma term
        inline void add_term(regge_trajectory* trajectory, std::tuple<int,double> coupling, std::array<double,2> sigmaparams)
        {
            auto new_term = new jpacTripleRegge(_kinematics, trajectory, {coupling}, {sigmaparams});
            _termsJPAC.push_back(new_term);
        };

        // One coupling but many sigma
        inline void add_term(regge_trajectory* trajectory, std::tuple<int,double> coupling, std::vector<std::array<double,2>> sigmaparams)
        {
            auto new_term = new jpacTripleRegge(_kinematics, trajectory, {coupling}, sigmaparams);
            _termsJPAC.push_back(new_term);
        };

        // Many coupling one sigma
        inline void add_term(regge_trajectory* trajectory, std::vector<std::tuple<int,double>> coupling, std::array<double,2> sigmaparams)
        {
            auto new_term = new jpacTripleRegge(_kinematics, trajectory, coupling, {sigmaparams});
            _termsJPAC.push_back(new_term);
        };
        
        // Passing multiple coupling and sigma terms
        inline void add_term(regge_trajectory* trajectory, std::vector<std::tuple<int,double>> coupling, std::vector<std::array<double,2>> sigmaparams)
        {
            auto new_term = new jpacTripleRegge(_kinematics, trajectory, coupling, sigmaparams);
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