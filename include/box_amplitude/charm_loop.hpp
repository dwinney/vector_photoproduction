// Charmonium production via a loop of open charm exchanges
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _CLOOP_
#define _CLOOP_

#include "constants.hpp"
#include "amplitudes/amplitude.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"

#include "cubature.h"

namespace jpacPhoto
{
    class loop_discontinuity
    {
        public:
        loop_discontinuity()
        {
            kgamD.set_JP(0, -1);
            kpsiD.set_JP(0, -1);

            gamDDstarEx.set_params({_gGamDDstar, _gDstarNLam, 0.});
            psiDDstarEx.set_params({_gPsiDDstar, _gDstarNLam, 0.});
        };
        
        // Evaluate the disconitnuity in terms of the intermediate state angles
        // Sum of the D and Dstar loops
        inline std::complex<double> eval(double sp, double thetaGam, double phiGam)
        {
            return eval_d(sp, thetaGam, phiGam);
        };

        inline void pass_params(std::array<int,4> helicities, double s, double theta, double eta)
        {
            // Save energy and theta for the external reaction
            _s = s; _theta = theta;
            
            // Save helicities
            _lam_gam  = helicities[0];
            _lam_tar  = helicities[1];
            _lam_psi  = helicities[2];
            _lam_rec  = helicities[3];

            // Save the updated form-factor parameter
            _eta = eta;
            gamDDstarEx.set_formfactor(2, M_DSTAR + _eta * _lambdaQCD);
            psiDDstarEx.set_formfactor(2, M_DSTAR + _eta * _lambdaQCD);
        };

         inline double s(){ return _s; };

        //---------------------------------------------------------------------------
        private:

        // energies and momentum transfers
        double _s, _theta;
        
        // Momentum transfers
        double _tGamma, _tPsi;

        // External helicites
        int _lam_gam, _lam_tar, _lam_psi, _lam_rec;

        // Form factor parameter
        double _eta;

        // All couplings assumed fixed:
        double _lambdaQCD  = 0.250, _e = sqrt(4. * PI * ALPHA);

        // Photon couplings
        double _gGamDDstar = 0.134, _gGamDstarDstar = 0.641;

        // Lambda_C couplings
        double _gDNLam     = -4.3,  _gDstarNLam     = -13.2;
        double _gPsiLamLam = -1.4;

        // Psi couplings
        double _gPsiDD          = 7.44;
        double _gPsiDDstar      = _gPsiDD / sqrt(M_D * M_DSTAR);
        double _gPsiDstarDstart = _gPsiDD * M_DSTAR / M_D;

        // ---------------------------------------------------------------------------
        // Functions to evaluate the D meson loop contribution
        std::complex<double> eval_d(double sp, double thetaGam, double phiGam);

        // Kinematics for the sub-process gamma p -> Lam D
        reaction_kinematics kgamD = reaction_kinematics(M_D, M_LAMBDAC, M_PROTON);
        reaction_kinematics kpsiD = reaction_kinematics(M_D, M_LAMBDAC, M_PROTON, M_JPSI);

        // Dstar exchanges
        vector_exchange gamDDstarEx = vector_exchange(&kgamD, M_DSTAR); 
        vector_exchange psiDDstarEx = vector_exchange(&kpsiD, M_DSTAR);
        
        amplitude_sum gamD = amplitude_sum(&kgamD, {&gamDDstarEx});
        amplitude_sum psiD = amplitude_sum(&kpsiD, {&psiDDstarEx});

        // ---------------------------------------------------------------------------
    };

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

        // Evaluate loop by dispersing over discontinuity
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // only vector kinematics allowed
        inline std::vector<std::array<int,2>> allowedJP()
        {
            return {{1, -1}};
        };

        protected:

        // Two free parameters:
        double _eta;  // formfactor parameter
        double _qmax; // integration cutoff

        // Wrapper for cubature
        loop_discontinuity _disc = loop_discontinuity();
        static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
    };
};



#endif