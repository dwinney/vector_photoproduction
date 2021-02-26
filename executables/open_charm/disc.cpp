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
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "one_loop/box_discontinuity.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Form factor parameter
    double eta  = 1.;
    double qmax = 1.;

    // Couplings
    double lambdaQCD       = 0.250, e = sqrt(4. * PI * ALPHA);
    double gGamDDstar      = 0.134, gGamDstarDstar = 0.641;
    double gDNLam          = -4.3,  gDstarNLam     = -13.2;
    double gPsiLamLam      = -1.4;
    double gPsiDD          = 7.44;
    double gPsiDDstar      = gPsiDD / sqrt(M_D * M_DSTAR);
    double gPsiDstarDstar  = gPsiDD * M_DSTAR / M_D;

    // ---------------------------------------------------------------------------
    // D loop 
    // ---------------------------------------------------------------------------
    // Kinematics for the sub-processes (gamma / psi) p -> Lam D
    auto kgamD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON);
    kgamD->set_JP(0, -1);  // Pseudo-scalar production

    auto kpsiD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON, M_JPSI);
    kpsiD->set_JP(0, -1);  // Pseudo-scalar production

    // Photon amplitudes
    auto gamDDstarEx = new vector_exchange(kgamD, M_DSTAR); 
    gamDDstarEx->set_params({gGamDDstar, gDstarNLam, 0.});
    gamDDstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    auto gamDLambdaEx = new dirac_exchange(kgamD, M_LAMBDAC);
    gamDLambdaEx->set_params({e, gDNLam});
    gamDLambdaEx->set_formfactor(2, M_LAMBDAC + lambdaQCD * eta);

    // Psi amplitudes
    auto psiDDstarEx = new vector_exchange(kpsiD, M_DSTAR);
    psiDDstarEx->set_params({gPsiDDstar, gDstarNLam, 0.});
    psiDDstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    auto psiDLambdaEx = new dirac_exchange(kpsiD, M_DSTAR);
    psiDLambdaEx->set_params({gPsiLamLam, gDNLam});
    psiDLambdaEx->set_formfactor(2, M_LAMBDAC + lambdaQCD * eta);

    // Sums 
    auto psiAmp = new amplitude_sum(kpsiD, {psiDDstarEx, psiDLambdaEx});
    auto gamAmp = new amplitude_sum(kgamD, {gamDDstarEx, gamDLambdaEx});

    // ---------------------------------------------------------------------------
    // Loop
    // ---------------------------------------------------------------------------

    auto disc = new box_discontinuity(gamAmp, psiAmp);
    disc->set_externals({1, 1, -1, 1}, 0.);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int Np = 100;
    double low  = E_beam(M_D + M_LAMBDAC);
    double high = 9.5;
    std::string filename = "disc.pdf";

    std::vector<double> E;
    std::vector<double> discontinuity;

    for (int i = 0; i <= Np; i++)
    {
        double Ei = low + EPS + double(i) * (high - low) / double(Np);
        E.push_back(Ei);

        double si = W_cm(Ei) * W_cm(Ei);

        double fx = disc->eval(si);
        discontinuity.push_back(fx);

        std::cout << std::left;
        std::cout << std::setw(7)  << i;
        std::cout << std::setw(15) << Ei;
        std::cout << std::setw(30) << fx;
        std::cout << std::endl;
    }

    jpacGraph1D* plotter = new jpacGraph1D();
    plotter->AddEntry(E, discontinuity, "Disc");

    plotter->SetLegend(false);

    plotter->SetXaxis("E_{#gamma}   [GeV]", low, high);
    plotter->Plot(filename);

    return 1;
};