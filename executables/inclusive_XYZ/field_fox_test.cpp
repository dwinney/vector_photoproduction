
#include "regge_trajectory.hpp"
#include "inclusive/inclusive_kinematics.hpp"
#include "inclusive/triple_regge.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Amplitude
    // ---------------------------------------------------------------------------

    auto kTest = new inclusive_kinematics(0.);
    auto field_and_fox = new triple_regge(kTest, "Field & Fox");

    // trajectories
    auto alphaP = new linear_trajectory(+1,  1., 0.37, "Pomeron");
    auto alphaR = new linear_trajectory(+1, 0.5,   1., "Reggeon");

    // Triple Reggeon
    field_and_fox->add_term({alphaR, alphaR, alphaR}, {18.1, 12.});
    
    // Reggeon-Reggeon-Pomeron
    field_and_fox->add_term({alphaR, alphaR, alphaP}, {26.81, 7.26});
    field_and_fox->add_term({alphaR, alphaR, alphaP}, {4.80, -1.83});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    amps.push_back(field_and_fox);

    int N = 100;

    double W = 50.;
    double x = 0.9;

    double  xmin = 0.;
    double  xmax = 1.;

    double  ymin = 1.E-2;
    double  ymax = 1.E2;

    std::string filename = "FF.pdf";
    std::string xlabel   = "#it{-t} [GeV^{2}]";
    std::string ylabel   = "E #frac{d#sigma}{d^{3}p}      [mb]";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        auto F = [&](double mt)
        {
            double s = W*W;
            double M2 = s * (1. - x);
            return amps[n]->invariant_xsection(s, -mt, M2);
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, xmin, xmax, true);

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);
    
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "x = " << x << ",  W_{#gammaN} = " << W << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.2, 0.75, header);
    plotter->SetLegendOffset(0.4, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
    delete alphaP;
    delete alphaR;
    delete field_and_fox;
    delete kTest;

    return 1.;
};