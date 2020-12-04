#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/baryon_resonance.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

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
    auto ptr = new reaction_kinematics(mJpsi);
    ptr->set_JP(1, -1);

    // ---------------------------------------------------------------------------
    // T - CHANNEL // this is the same for all cases

    // Set up pomeron trajectory
    // best fit values from [1]
    auto alpha = new linear_trajectory(+1, 0.941, 0.364);

    // Create amplitude with kinematics and trajectory
    auto background = new pomeron_exchange(ptr, alpha, false, "Background");

    // normalization and t-slope
    // best fit values from [1]
    std::vector<double> back_params = {0.379, 0.12};
    background->set_params(back_params);

    // ---------------------------------------------------------------------------
    // S - CHANNEL 

    // Pc4312
    auto Pc4312_1m = new baryon_resonance(ptr, 1, -1, 4.3119, 9.8E-3, "P_{c}(4312)");
    Pc4312_1m->set_params({0.01, .7071});

    auto Pc4312_3m = new baryon_resonance(ptr, 3, -1, 4.3119, 9.8E-3, "P_{c}(4312)");
    Pc4312_3m->set_params({0.01, .7071});

    // Pc4440
    auto Pc4440_1m = new baryon_resonance(ptr, 1, -1, 4.4403, 20.6E-3, "P_{c}(4440)");
    Pc4440_1m->set_params({0.01, .7071});

    auto Pc4440_3m = new baryon_resonance(ptr, 3, -1, 4.4403, 20.6E-3, "P_{c}(4440)");
    Pc4440_3m->set_params({0.01, .7071});

    auto Pc4440_3p = new baryon_resonance(ptr, 3, +1, 4.4403, 20.6E-3, "P_{c}(4440)");
    Pc4440_3p->set_params({0.01, .7071});

    // Pc4457
    auto Pc4457_1m = new baryon_resonance(ptr, 1, -1, 4.4573, 6.4E-3, "P_{c}(4457)");
    Pc4457_1m->set_params({0.01, .7071});

    auto Pc4457_3m = new baryon_resonance(ptr, 3, -1, 4.4573, 6.4E-3, "P_{c}(4457)");
    Pc4457_3m->set_params({0.01, .7071});

    auto Pc4457_5m = new baryon_resonance(ptr, 5, -1, 4.4573, 6.4E-3, "P_{c}(4457)");
    Pc4457_5m->set_params({0.01, .7071});


    // Add to the sum
    auto sumA = new amplitude_sum (ptr, {background, Pc4312_1m, Pc4440_3m}, "A");
    auto sumB = new amplitude_sum (ptr, {background, Pc4312_3m, Pc4440_1m}, "B");
    auto sumC = new amplitude_sum (ptr, {background, Pc4312_3m, Pc4440_3p}, "C");

    // ---------------------------------------------------------------------------
    // Choose which scenario to plot
    std::vector<amplitude*> amps;
    amps.push_back(background);
    // amps.push_back(sumA);
    // amps.push_back(sumB);
    // amps.push_back(sumC);

    int N = 100;
    bool PRINT_TO_COMMANDLINE = true;

    double ymax =  1.;
    double ymin = -0.1;

    double Emin = E_beam(ptr->Wth()) + EPS;
    double Emax = 12.;

    // double Emin = 10.;
    // double Emax = 10.1;

    std::string filename = "sigma_integrated.pdf";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter objects
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // scan over theta
    // for (int n = 0; n < amps.size(); n++)
    // {
    //     std::cout << "\nPrinting amplitude " << amps[n]->identifier << ".\n";

    //     auto F = [&](double theta)
    //     {

    //         auto f = [&](double E)
    //         {
    //             double W = W_cm(E);
    //             double t = ptr->t_man(W*W, theta * deg2rad);
    //             return amps[n]->beam_asymmetry_4pi(W*W, t);
    //         };

    //         ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
    //         ROOT::Math::Functor1D wF(f);
    //         ig.SetFunction(wF);
            
    //         return ig.Integral(Emin, Emax);
    //     };

    //     std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, 0., 90., PRINT_TO_COMMANDLINE);
    //     plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
    // }

    // // Add a header to legend to specify the fixed energy
    // plotter->SetLegend(0.2, 0.7);

    // plotter->SetXaxis("#theta", 0., 90.);

    // // To change the range of the Y-axis or the position of the Legend change the arguments here
    // plotter->SetYaxis("#Sigma integrated over E_{#gamma} = \{th, 12.\}", ymin, ymax);


    // plotter->Plot(filename.c_str());

    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude " << amps[n]->identifier << ".\n";

        auto F = [&](double egam)
        {
            double W = W_cm(egam);
            auto f = [&](double t)
            {
                return amps[n]->beam_asymmetry_4pi(W*W, t);
            };

            ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
            ROOT::Math::Functor1D wF(f);
            ig.SetFunction(wF);
            
            double t_min = amps[n]->kinematics->t_man(W*W, 0.);
            double t_max = amps[n]->kinematics->t_man(W*W, M_PI);
            return ig.Integral(t_max, t_min);
        };

        std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, Emin, Emax, PRINT_TO_COMMANDLINE);
        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
    }

    // Add a header to legend to specify the fixed energy
    plotter->SetLegend(0.2, 0.7);

    plotter->SetXaxis("E_{#gamma}", Emin, Emax);

    // To change the range of the Y-axis or the position of the Legend change the arguments here
    plotter->SetYaxis("#Sigma integrated over t", ymin, ymax);


    plotter->Plot(filename.c_str());

    return 1.;
};
