#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/pomeron_exchange.hpp"

#include "photoPlotter.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // ---------------------------------------------------------------------------
    // AMPLITUDES
    // ---------------------------------------------------------------------------

    // Set up Kinematics for jpsi in final state
    auto * ptr = new reaction_kinematics(mJpsi);
    ptr->set_JP(1, -1);

    // ---------------------------------------------------------------------------
    // Our amplitude

    auto alpha1 = new linear_trajectory(+1, 0.941, 0.364);
    auto background1 = new pomeron_exchange(ptr, alpha1, 0, "JPAC");

    // Normalization and t-slope
    background1->set_params({0.379, 0.12});


    // ---------------------------------------------------------------------------
    // Wang et al. amplitude

    auto alpha2 = new linear_trajectory(+1, 1. - 0.08, 0.25);
    auto background2 = new pomeron_exchange(ptr, alpha2, 2, "Wang et al.");

    // Pomeron-charm coupling and cutoff
    background2->set_params({sqrt(0.8), 1.2});


    // ---------------------------- -----------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(background1);
    amps.push_back(background2);

    auto plotter = new photoPlotter(amps);

    plotter->N = 30;
    plotter->PRINT_TO_COMMANDLINE = true;
    plotter->LAB_ENERGY = true;

    plotter->xmin = E_beam(ptr->Wth()) + EPS;
    plotter->xmax = 12.;

    plotter->ymin = 0.;
    plotter->ymax = 3.5;

    plotter->SHOW_LEGEND = true;
    plotter->xlegend = 0.2;
    plotter->ylegend = 0.6;

    plotter->filename  = "jpsi_compare.pdf";
    plotter->ylabel    = "#it{#sigma(#gamma p #rightarrow J/#psi p)}  [nb]";
    plotter->xlabel    = "#it{E_{#gamma}}  [GeV]";

    plotter->Plot("integrated_xsection");

    delete ptr, alpha1, alpha2, background1, background2, plotter;
    return 1;
};
