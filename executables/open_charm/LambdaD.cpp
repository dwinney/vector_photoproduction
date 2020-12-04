// ---------------------------------------------------------------------------
// Photoproduction cross-sections of the Lambda_c Dbar/D* final state
// by open charm exchanges 
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:2009.08345v1
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "photoPlotter.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Set up Kinematics for Dbar LambdaC in final state
    auto ptr = new reaction_kinematics(mD, mLambdaC, mPro);
    ptr->set_JP(0, -1);

    // Amplitude
    double eta = 1.;

    auto dstarEx = new vector_exchange(ptr, mDstar, "D^{*} exchange");
    dstarEx->set_params({0.134, -13.2, 0.});
    dstarEx->set_formfactor(2, mDstar + eta * 0.250);

    auto lamcEx = new dirac_exchange(ptr, mLambdaC, "#Lambda_{c} exchange");
    lamcEx->set_params({sqrt(4.*M_PI*M_ALPHA), -4.3});
    lamcEx->set_formfactor(2, mLambdaC + eta * 0.250);

    auto sum = new amplitude_sum(ptr, {dstarEx, lamcEx}, "Sum");

    // ---------------------------- -----------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(dstarEx);
    amps.push_back(lamcEx);
    amps.push_back(sum);

    auto plotter = new photoPlotter(amps);

    plotter->N = 30;
    plotter->PRINT_TO_COMMANDLINE = true;
    plotter->LAB_ENERGY = true;

    plotter->xmin = 8.5;
    plotter->xmax = 10.5;

    plotter->ymin = 0.;
    plotter->ymax = 230.;

    // plotter->SHOW_LEGEND = true;
    // plotter->xlegend = 0.2;
    // plotter->ylegend = 0.6;

    plotter->filename  = "open_charm.pdf";
    plotter->ylabel    = "#it{#sigma(#gamma p #rightarrow #bar{D} #Lambda_{c}^{+})}  [nb]";
    plotter->xlabel    = "#it{E_{#gamma}}  [GeV]";

    plotter->Plot("integrated_xsection");

    delete ptr;
    return 1;
};