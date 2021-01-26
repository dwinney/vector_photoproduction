// Charmonium production via a loop of open charm exchanges
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _CLOOP_
#define _CLOOP_

#include "constants.hpp"
#include "amplitude.hpp"
#include "reaction_kinematics.hpp"
#include "vector_exchange.hpp"
#include "pseudoscalar_exchange.hpp"
#include "dirac_exchange.hpp"

namespace jpacPhoto
{
    class charm_loop : public amplitude
    {
        public:

        // Constructor
        charm_loop(reaction_kinematics * xkinem, std::string name = "charm_loop")
        : amplitude(xkinem, name)
        {
            set_nParams(2);
            check_JP(xkinem->_jp);
        };

        // Destructor
        ~charm_loop()
        {};

        // Setting utility
        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _eta  = params[0];
            _qmax = params[1];
        };

        protected:

        // Two free parameters:
        double _eta;  // formfactor parameter
        double _qmax; // integration cutoff

        // All other couplings assumed fixed:
        double _lambdaQCD  = 0.250, _e = sqrt(4. * PI * ALPHA);
        double _gGamDDstar = 0.134, _gGamDstarDstar = 0.641;
        double _gDNLam     = -4.3,  _gDstarNLam     = -13.2;
        double _gPsiDD     = 7.44,  _gPsiLamLam     = -1.4;

        inline void set_up_amplitudes()
        {
            gamDDstarEx.set_params({_gGamDDstar, _gDstarNLam, 0.});
            gamDLamEx.set_params({_e, _gDNLam, 0.});
            gamDstarDEX.set_params({_gGamDDstar, _gDNLam, 0.});
        };

        // Kinematics for the sum-process jpsi p -> Lam D(*)
        reaction_kinematics kgamD     = reaction_kinematics(M_D,     M_LAMBDAC, M_PROTON);
        reaction_kinematics kgamDstar = reaction_kinematics(M_DSTAR, M_LAMBDAC, M_PROTON);

        vector_exchange gamDDstarEx = vector_exchange( &kgamD, M_DSTAR,   "#gamma p #to D #Lambda_{c} (D^{*} exchange)"); 
        dirac_exchange  gamDLamEx   = dirac_exchange(  &kgamD, M_LAMBDAC, "#gamma p #to D #Lambda_{c} (#Lambda_{c} exchange)"); 

        pseudoscalar_exchange gamDstarDEX = pseudoscalar_exchange( &kgamDstar, M_D,       "#gamma p #to D^{*} #Lambda_{c} (D exchange)");
        vector_exchange gamDstarDstarEx   = vector_exchange(       &kgamDstar, M_DSTAR,   "#gamma p #to D^{*} #Lambda_{c} (D^{*} exchange)"); 
        dirac_exchange  gamDstarLamEx     = dirac_exchange(        &kgamDstar, M_LAMBDAC, "#gamma p #to D^{*} #Lambda_{c} (#Lambda_{c} exchange)"); 

        // Kinematics for the sum-process jpsi p -> Lam D(*)
        // reaction_kinematics kpsiD(M_D, M_LAMBDAC, M_PROTON, M_JPSI);
        // reaction_kinematics kpsiDstar(M_DSTAR, M_LAMBDAC, M_PROTON, M_JPSI);

    };
};

#endif