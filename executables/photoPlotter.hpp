// ---------------------------------------------------------------------------
// Wrapper for the jpacGraph1D class to get rid of copying the same code snippets
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PHOTOPLOT_
#define _PHOTOPLOT_

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include "reaction_kinematics.hpp"
#include "amplitudes/amplitude.hpp"

#include <vector>
#include <string>

namespace jpacPhoto
{
    class photoPlotter : public jpacGraph1D
    {
        public:
        
        // default constructor
        photoPlotter( std::vector<amplitude*> amps_)
        : amps(amps_)
        {};
        
        int N = 20;
        bool PRINT_TO_COMMANDLINE = true;
        bool LAB_ENERGY = false;
        double xmin, xmax, ymin, ymax;
        std::string xlabel, ylabel, filename;

        void plot_integrated_xsection()
        {
            
            for (int n = 0; n < amps.size(); n++)
            {
                std::cout << std::endl << "Printing amplitude: " << amps[n]->identifier << "\n";
    
                double th;
                (LAB_ENERGY) ? (th = E_beam(amps[n]->kinematics->Wth())) : (th = amps[n]->kinematics->Wth());

                auto F = [&](double x)
                {
                    double W;
                    (LAB_ENERGY) ? (W = W_cm(x)) : (W = x);
                    return amps[n]->integrated_xsection(W*W);
                };

                std::array<std::vector<double>, 2> x_fx;
                if (xmin < th)
                {
                    x_fx = vec_fill(N, F, th + EPS, xmax, PRINT_TO_COMMANDLINE);
                }
                else
                {
                    x_fx = vec_fill(N, F, xmin, xmax, PRINT_TO_COMMANDLINE);
                }

                AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
            }

            SetXaxis(xlabel, xmin, xmax);
            SetYaxis(ylabel, ymin, ymax);
            SetLegend(false);

            // Output to file
            Plot(filename);
        };

        void plot_differential_xsection(double theta)
        {
            
            for (int n = 0; n < amps.size(); n++)
            {
                std::cout << std::endl << "Printing amplitude: " << amps[n]->identifier << "\n";
    
                double th;
                (LAB_ENERGY) ? (th = E_beam(amps[n]->kinematics->Wth())) : (th = amps[n]->kinematics->Wth());

                auto F = [&](double x)
                {
                    double W;
                    (LAB_ENERGY) ? (W = W_cm(x)) : (W = x);
                    double t = amps[n]->kinematics->t_man(W*W, theta * deg2rad);
                    return amps[n]->differential_xsection(W*W, t);
                };

                std::array<std::vector<double>, 2> x_fx;
                if (xmin < th)
                {
                    x_fx = vec_fill(N, F, th + EPS, xmax, PRINT_TO_COMMANDLINE);
                }
                else
                {
                    x_fx = vec_fill(N, F, xmin, xmax, PRINT_TO_COMMANDLINE);
                }

                AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
            }

            SetXaxis(xlabel, xmin, xmax);
            SetYaxis(ylabel, ymin, ymax);
            SetLegend(false);

            // Output to file
            Plot(filename);
        };

        private:
        std::vector<amplitude*> amps;
    }; 
};

#endif