
#include "regge_trajectory.hpp"
#include "inclusive/inclusive_kinematics.hpp"
#include "inclusive/triple_regge.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>
#include <memory>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Amplitudes
    // ---------------------------------------------------------------------------

    // For all Z's we only have pion exchange
    auto alphaPi = new linear_trajectory(+1, -M2_PION * 0.7, 0.7, "#pi trajectory");

    auto kZc3900 = new inclusive_kinematics(M_ZC3900);
    auto Zc3900 = new triple_regge(kZc3900, "Z_{c}(3900)");


    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    // amps.push_back(field_and_fox);
    // amps.push_back(vincent);

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
    auto plotter = new jpacGraph1D();

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
    streamObj << std::setprecision(4) << "x = " << x << ",  W = " << W << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.2, 0.75, header);
    plotter->SetLegendOffset(0.4, 0.1);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};