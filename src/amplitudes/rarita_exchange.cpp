// Spin-3/2 exchange amplitude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/rarita_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> rarita_exchange::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  std::complex<double> result = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp;
      temp  = top_vertex(i, lam_gam, lam_rec, s, zs);
      temp *= rarita_propagator(i, j, s, zs);
      temp *= bottom_vertex(j, lam_vec, lam_targ, s, zs);

      result += temp;
    }
  }

  return result;
};

//------------------------------------------------------------------------------
// rank-2 traceless tensor
std::complex<double> rarita_exchange::g_bar(int mu, int nu, double s, double zs)
{
  std::complex<double> result;
  result = exchange_momentum(mu, s, zs) * exchange_momentum(nu, s, zs) / mEx2;

  if (mu == nu)
  {
    result -= metric[mu];
  }

  return result;
};

// g_bar contracted with gamma^nu
std::complex<double> rarita_exchange::slashed_g_bar(int mu, int i, int j, double s, double zs)
{
  std::complex<double> result = 0.;

  for (int nu = 0; nu < 4; nu++)
  {
    std::complex<double> temp;
    temp  = g_bar(mu, nu, s, zs);
    temp *= metric[nu];
    temp *= gamma_matrices[nu][i][j];

    result += temp;
  }

  return result;
};

//------------------------------------------------------------------------------
// Relative momentum either entering (top vertex) or exiting (bottom vertex) the propagator
std::complex<double> rarita_exchange::relative_momentum(int mu, double s, double zs, std::string in_out)
{
  std::complex<double> q1_mu, q2_mu;

  if ((in_out == "in") || (in_out == "top") || (in_out == "initial") )
  {
    q1_mu = kinematics->initial.component(mu, "beam",   s, 1.);
    q2_mu = kinematics->initial.component(mu, "target", s, 1.);
  }
  else if ((in_out == "out") || (in_out == "bot") || (in_out == "final"))
  {
    q1_mu = kinematics->final.component(mu, kinematics->vector_particle, s, zs);
    q2_mu = kinematics->final.component(mu, "recoil",                    s, zs);
  }
  else
  {
    std::cout << "Error! Unkown parameter: " << in_out << "passed to relative_momentum. ";
    std::cout << "Quitting...";
    exit(1);
  }

  return q1_mu - q2_mu;
}

//------------------------------------------------------------------------------
// Rarita-Schwinger Propagator
std::complex<double> rarita_exchange::rarita_propagator(int i, int j, double s, double zs)
{
  std::complex<double> result = 0.;

  for (int mu = 0; mu < 4; mu++)
  {
    for(int nu = 0; nu < 4; nu++)
    {
      std::complex<double> term_1;
      term_1  = relative_momentum(mu, s, zs, "in");
      term_1 *= metric[mu];
      term_1 *= g_bar(mu, nu, s, zs);
      term_1 *= metric[nu];
      term_1 *= relative_momentum(nu, s, zs, "out");

      std::complex<double> term_2;
      term_2  = relative_momentum(mu, s, zs, "in");
      term_2 *= metric[mu];
      term_2 *= slashed_g_bar(mu, i, j, s, zs);
      term_2 *= slashed_g_bar(nu, i, j, s, zs);
      term_2 *= metric[nu];
      term_2 *= relative_momentum(nu, s, zs, "out");

      result += -term_1 + term_2 / 3.;
    }
  }

  result *= dirac_propagator(i, j, s, zs);
  
  return result;
}
