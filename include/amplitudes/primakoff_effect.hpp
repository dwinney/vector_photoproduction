// Axial-vector meson photoproduction proceeding through the primakoff effect
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PRIMAKOFF_
#define _PRIMAKOFF_

#include "amplitude.hpp"
#include "gamma_technology.hpp"

namespace jpacPhoto
{
    class primakoff_effect : public amplitude
    {
        public:
        // Constructor 
        primakoff_effect(reaction_kinematics * xkinem, std::string amp_id = "")
        : amplitude(xkinem, amp_id, 4)
        {};

        void set_params(std::vector<double> params)
        {
            check_Nparams(params); 
            Z = params[0];
            R = params[1];
            a = params[2];
            g = params[3];

            calculate_norm();
        };

        inline void set_Q2(double q2)
        {
            Q2 = q2;
        };

        // individual helicity amplitudes not supported but need to provide definition for virtual class.
        inline std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t)
        {
            std::cout << "Warning! Individual helicity amplitudes not supported by primakoff_effect!\n";
            return 0.;
        }

        // instead we override the definition of differential_xsection in amplitude.hpp
        double differential_xsection(double s, double t);
        double integrated_xsection(double s);

        // Fermi model nuclear charge distribution
        inline double charge_distribution(double r)
        {
            return 1. / ( 1. + exp((r - R) / a) );
        };

        // Normalized fourier transform of the above charge_distributions 
        double form_factor(double x);
        double F_0;

        private:

        // Parameters
        int    Z     = 0 ;  // atomic number
        double R     = 0.;  // radius parameter
        double a     = 0.;  // skin thickness parameter
        double g     = 0.;  // X -> gamma gamma* coupling
        double rho_0 = 0.;  // normalizaton
        double Q2 = 0.;      // Q^2 of the external photon

        // Calculate the normalization with above parameters
        void calculate_norm();        

        // Top tensor (spin averaged tensor of X -> gamma gamma* interaction)
        std::complex<double> top_tensor(int mu, int nu);

        // Hadronic bottom tensor (gamma* -> N Nbar)
        std::complex<double> bottom_tensor(int mu, int nu);       
    };
};

#endif