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
#include "amplitudes/pseudoscalar_exchange.hpp"
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
    // Form factor parameter
    double eta = 1.;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    auto kD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON);
    kD->set_JP(0, -1);

    auto d_dstarEx = new vector_exchange(kD, M_DSTAR, "D^{*} exchange");
    d_dstarEx->set_params({0.134, -13.2, 0.});
    d_dstarEx->set_formfactor(2, M_DSTAR + eta * 0.250);

    auto d_lamcEx = new dirac_exchange(kD, M_LAMBDAC, "#Lambda_{c} exchange");
    d_lamcEx->set_params({sqrt(4.* PI * ALPHA), -4.3});
    d_lamcEx->set_formfactor(2, M_LAMBDAC + eta * 0.250);

    auto d_sum = new amplitude_sum(kD, {d_dstarEx, d_lamcEx}, "Sum");

    // ---------------------------------------------------------------------------
    // D* phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for D* LambdaC in final state
    auto kDstar = new reaction_kinematics(M_DSTAR, M_LAMBDAC, M_PROTON);
    kDstar->set_JP(1, -1);

    auto dstar_dEx = new pseudoscalar_exchange(kDstar, M_D, "D exchange");
    dstar_dEx->set_params({0.134, -4.3, 0.});
    dstar_dEx->set_formfactor(2, M_D + eta * 0.250);

    auto dstar_dstarEx = new vector_exchange(kDstar, M_DSTAR, "D^{*} exchange");
    dstar_dstarEx->set_params({0.641, -13.2, 0.});
    dstar_dstarEx->set_formfactor(2, M_DSTAR + eta * 0.250);

    auto dstar_lamcEx = new dirac_exchange(kDstar, M_LAMBDAC, "#Lambda_{c} exchange");
    dstar_lamcEx->set_params({sqrt(4.* PI * ALPHA), -13.2});
    dstar_lamcEx->set_formfactor(2, M_LAMBDAC + eta * 0.250);

    auto dstar_sum = new amplitude_sum(kDstar, {dstar_dEx, dstar_dstarEx, dstar_lamcEx}, "Sum");

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(d_sum);
    amps.push_back(dstar_sum);

    auto plotter = new photoPlotter(amps);

    plotter->N = 30;
    plotter->PRINT_TO_COMMANDLINE = true;
    plotter->LAB_ENERGY = true;

    plotter->xmin = 8.5;
    plotter->xmax = 10.5;

    plotter->ymin = 0.;
    plotter->ymax = 250.;

    plotter->SHOW_LEGEND = false;
    plotter->xlegend = 0.2;
    plotter->ylegend = 0.6;

    plotter->filename  = "open_charm.pdf";
    plotter->ylabel    = "#it{#sigma(#gamma p #rightarrow #bar{D} #Lambda_{c}^{+})}  [nb]";
    plotter->xlabel    = "#it{E_{#gamma}}  [GeV]";

    plotter->Plot("integrated_xsection");

    delete kD, kDstar;
    return 1;
};