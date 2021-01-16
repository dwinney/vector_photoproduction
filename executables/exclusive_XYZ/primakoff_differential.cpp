// --------------------------------------------------------------------------- 
// Prediction for X(3872) Primakoff production off a nuclear target
//
// USAGE:
// make primakoff_differential && ../bin/primakoff_differential
//
// OUTPUT:
// x_primakoff.pdf
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/primakoff_effect.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Amplitude
    // ---------------------------------------------------------------------------
    double Q2 = 0.5;   // fixed Q2
    double W  = 2.;    // Fixed energy per nucleon
    double mX = 3.872; // Mass of produced X

    // Uranium
    double mU = 221.6977; // GeV
    reaction_kinematics * kU = new reaction_kinematics(mX, mU, mU);
    kU->set_Q2(Q2);
    kU->set_JP(1, 1);

    primakoff_effect U(kU, "^{238}U");
    U.set_params({92, 34.48, 3.07, 3.2E-3});

    // Tin
    double mSn = 115.3924; // GeV
    reaction_kinematics * kSn = new reaction_kinematics(mX, mSn, mSn);
    kSn->set_Q2(Q2);
    kSn->set_JP(1, 1);

    primakoff_effect Sn(kSn, "^{124}Sn");
    Sn.set_params({50, 27.56, 2.73, 3.2E-3});

    // Zinc
    double mZn = 65.1202; // GeV
    reaction_kinematics * kZn = new reaction_kinematics(mX, mZn, mZn);
    kZn->set_Q2(Q2);
    kZn->set_JP(1, 1);

    primakoff_effect Zn(kZn, "^{70}Zn");
    Zn.set_params({30, 22.34, 2.954, 3.2E-3});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<primakoff_effect*> amps;
    amps.push_back(&Zn);
    amps.push_back(&Sn);
    amps.push_back(&U);

    // number of nucleons 
    double xNs[3] = {70., 124., 238.}; 

    int N = 400;
    std::string filename = "primakoff_differential.pdf";

    // x - axis params
    double  xmax = 0.1;
    std::string xlabel = "#it{-t}   [GeV^{2}]";

    // y - axis params
    double  ymin = 2.E-6;
    double  ymax = 100.;

    std::string ylabel  = "#it{d#sigma/dt} (#gamma* A #rightarrow X A)   [nb GeV^{-2}]";
    bool print_to_cmd = true;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        double s = W * W * xNs[n] * xNs[n];
        double xmin = -amps[n]->_kinematics->t_man(s, 0.);

        auto F = [&](double t)
        {
            return amps[n]->differential_xsection(s, -t);
        };

        std::array<std::vector<double>, 2> x_fx, x_fx2; 

        std::cout << std::endl << "Printing longitudinal xsection: " << amps[n]->_identifier << "\n";
        x_fx = vec_fill(N, F, xmin, xmax, print_to_cmd);
        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);

        std::cout << std::endl << "Printing tranverse xsection: " << amps[n]->_identifier << "\n";
        amps[n]->set_LT(1);
        x_fx2 = vec_fill  (N, F, xmin, xmax, print_to_cmd);
        plotter->AddDashedEntry(x_fx2[0], x_fx2[1]);
    }

    // Add a header to legend to specify the fixed energy
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "Q^{2} = " << Q2 << " GeV^{2},  W_{#gammaN} = " << W << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.52, 0.6, header);

    // Set up axes
    plotter->SetXaxis(xlabel, 0., xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);    
    plotter->SetYlogscale(true);

    // Output to file
    plotter->Plot(filename);

    // Cleanup
    delete plotter;
    delete kU, kSn, kZn;
    return 1.;
};