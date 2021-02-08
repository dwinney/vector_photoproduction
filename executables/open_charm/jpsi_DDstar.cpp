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
#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "box_amplitude/charm_loop.hpp"

#include "photoPlotter.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Form factor parameter
    double eta = 1.;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    auto kPsi = new reaction_kinematics(M_JPSI, M_PROTON);
    kPsi->set_JP(1, -1);

    auto loop = new charm_loop(kPsi, "test");
    loop->set_params({1., 1.});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(loop);

    auto plotter = new photoPlotter(amps);

    plotter->N = 30;
    plotter->PRINT_TO_COMMANDLINE = true;
    plotter->LAB_ENERGY = true;

    plotter->xmin = 8.5;
    plotter->xmax = 10.5;

    plotter->ymin = 0.;
    plotter->ymax = 10.;

    plotter->SHOW_LEGEND = false;
    plotter->xlegend = 0.2;
    plotter->ylegend = 0.6;

    plotter->filename  = "open_charm.pdf";
    plotter->ylabel    = "#it{#sigma(#gamma p #rightarrow #bar{D} #Lambda_{c}^{+})}  [nb]";
    plotter->xlabel    = "#it{E_{#gamma}}  [GeV]";

    plotter->Plot("differential_xsection", 0.);

    return 1;
};