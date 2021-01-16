// ---------------------------------------------------------------------------
// Prediction for Z_c and Z_b photoproduction based on pi reggeon exchange
// at high energies. 
//
// Reproduces right plot in FIG 3 of [1] 
// 
// USAGE:
// make Z_high && ./Z_high
//
// OUTPUT:
// Z_regge.pdf
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

  // ---------------------------------------------------------------------------
  // Preliminaries
  // ---------------------------------------------------------------------------

  double g_NN = sqrt(4. * PI * 13.81); // Nucleon coupling same for all
  double LamPi = .9;  // 900 MeV cutoff for formfactor
  double bPi = 1. / (LamPi * LamPi);

  // Zc(3900)
  reaction_kinematics * kZc = new reaction_kinematics(M_ZC3900);
  kZc->set_JP(1, 1);

  double gc_Psi = 1.91; // psi coupling before VMD scaling
  double gc_Gamma = E * F_JPSI * gc_Psi / M_JPSI;
  std::vector<double> Zc_couplings = {gc_Gamma, g_NN};

  // Zb(10610)
  reaction_kinematics * kZb = new reaction_kinematics(M_ZB10610);
  kZb->set_JP(1, 1);

  double gb_Ups1 = 0.49, gb_Ups2 = 3.30, gb_Ups3 = 9.22;
  double gb_Gamma = E * (F_UPSILON1S * gb_Ups1 / M_UPSILON1S 
                       + F_UPSILON2S * gb_Ups2 / M_UPSILON2S
                       + F_UPSILON3S * gb_Ups3 / M_UPSILON3S); 
  std::vector<double> Zb_couplings = {gb_Gamma, g_NN};

  
  // Zb(10650)
  reaction_kinematics * kZbp = new reaction_kinematics(M_ZB10650);
  kZbp->set_JP(1, 1);

  double gbp_Ups1 = 0.21, gbp_Ups2 = 1.47, gbp_Ups3 = 4.8;
  double gbp_Gamma = E * (F_UPSILON1S * gbp_Ups1 / M_UPSILON1S 
                       +  F_UPSILON2S * gbp_Ups2 / M_UPSILON2S
                       +  F_UPSILON3S * gbp_Ups3 / M_UPSILON3S);  
  std::vector<double> Zbp_couplings = {gbp_Gamma, g_NN};
  
  // Pion trajectory 
  int signature = +1;
  double alpha_prime = 0.7; // GeV^-2
  double alpha_0 =  - alpha_prime * M2_PION;
  linear_trajectory * alpha = new linear_trajectory(signature, alpha_0, alpha_prime);

  // ---------------------------------------------------------------------------
  // Reggeized amplitudes
  // ---------------------------------------------------------------------------

  pseudoscalar_exchange Zc(kZc, alpha, "#it{Z_{c}}(3900)^{+}");
  Zc.set_params(Zc_couplings);
  Zc.set_formfactor(true, bPi);

  pseudoscalar_exchange Zb(kZb, alpha,  "#it{Z_{b}}(10610)^{+}");
  Zb.set_params(Zb_couplings);
  Zb.set_formfactor(true, bPi);

  pseudoscalar_exchange Zbp(kZbp, alpha, "#it{Z_{b}}'(10650)^{+}");
  Zbp.set_params(Zbp_couplings);
  Zbp.set_formfactor(true, bPi);

  // ---------------------------------------------------------------------------
  // Plotting options
  // ---------------------------------------------------------------------------
  
  // which amps to plot
  std::vector<amplitude*> amps;
  amps.push_back(&Zc);
  amps.push_back(&Zb);
  amps.push_back(&Zbp);

  int   N = 100;
  bool PRINT_TO_COMMANDLINE     = true;

  // ---------------------------------------------------------------------------
  double  xmin = 20.;
  double  xmax = 70.;

  double  ymin = 1.E-4;
  double  ymax = 2.;
  // ---------------------------------------------------------------------------

  std::string ylabel    = "#it{#sigma(#gamma p #rightarrow Z n)}  [nb]";
  std::string xlabel    = "#it{W_{#gammap}}  [GeV]";
  std::string filename  = "Z_regge.pdf";

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // Plotter object
  jpacGraph1D* plotter = new jpacGraph1D();

  // ---------------------------------------------------------------------------
  // Print the desired observable for each amplitude
  for (int n = 0; n < amps.size(); n++)
  {
    std::cout << std::endl << "Printing amplitude: " << amps[n]->_identifier << "\n";

    auto F = [&](double x)
    {
      return amps[n]->integrated_xsection(x*x);
    };

    std::array<std::vector<double>, 2> x_fx, x_fx1;
    if (xmin < amps[n]->_kinematics->Wth())
    {
        x_fx = vec_fill(N, F, amps[n]->_kinematics->Wth() + EPS, xmax, PRINT_TO_COMMANDLINE);
        x_fx[0].insert(x_fx[0].begin(), amps[n]->_kinematics->Wth());
        x_fx[1].insert(x_fx[1].begin(), 0.);
    }
    else
    {
      x_fx = vec_fill(N, F, xmin, xmax, PRINT_TO_COMMANDLINE);
    }

    plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
  }

  plotter->SetXaxis(xlabel, xmin, xmax);
  plotter->SetYaxis(ylabel, ymin, ymax);
  plotter->SetYlogscale(1);
  plotter->SetLegend(0.2, 0.2);

  // Output to file
  plotter->Plot(filename);

  delete kZc, kZb, kZbp, plotter;

  return 1.;
}
