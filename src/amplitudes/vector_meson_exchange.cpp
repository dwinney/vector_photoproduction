// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/vector_meson_exchange.hpp"

// ---------------------------------------------------------------------------
// Assemble the helicity amplitude by contracting the lorentz indices
std::complex<double> vector_meson_exchange::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  std::complex<double> result = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    for(int nu = 0; nu < 4; nu++)
    {
      std::complex<double> temp;
      temp = top_vertex(mu, lam_gam, lam_vec, s, zs);
      temp *= metric[mu][nu];
      temp *= vector_propagator(mu, nu, s, zs);
      temp *= metric[mu][nu];
      temp *= bottom_vertex(nu, lam_targ, lam_rec, s, zs);

      result += temp;
    }
  }

  return result;
};

// ---------------------------------------------------------------------------
// Photon - Axial Vector - Vector vertex
std::complex<double> vector_meson_exchange::top_vertex(int mu, int lam_gam, int lam_vec, double s, double zs)
{
  // Contract with LeviCivita
  std::complex<double> result = 0.;
  for (int alpha = 0; alpha < 3; alpha++)
  {
    for (int beta = 0; beta < 3; beta++)
    {
      for (int gamma = 0; gamma < 3; gamma++)
      {
        std::complex<double> temp;
        temp = LeviCivita(mu, alpha, beta, gamma);
        temp *= exchange_momenta(alpha, s, zs);
        temp *= kinematics->eps_gamma.component(beta, lam_gam, s, 1.);
        temp *= kinematics->eps_vec.component(gamma, lam_vec, s, zs);

        result += temp;
      }
    }
  }

  // Multiply by coupling
  return result * gGamma;
};

// ---------------------------------------------------------------------------
// Nucleon - Nucleon - Vector vertex
std::complex<double> vector_meson_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec, double s, double zs)
{
  std::complex<double> vector = 0., tensor = 0.;

  // Vector coupling piece
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp;
      temp = kinematics->recoil.adjoint_component(i, lam_rec, s, zs);
      temp *= gamma_matrices[mu][i][j];
      temp *= kinematics->target.component(j, lam_targ, s, 1.);

      vector += temp;
    }
  }

  // Tensor coupling piece
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp = 0., temp2 = 0.;
      temp = kinematics->recoil.adjoint_component(i, lam_rec, s, zs);

      for (int nu = 0; nu < 4; nu++)
      {
        temp2 += Sigma(mu, nu, i, j) * exchange_momenta(nu, s, zs) / (2. * mPro);
      }
      temp *= temp2;

      temp *= kinematics->target.component(j, lam_targ, s, 1.);
    }
  }

  return gV * vector - gT * tensor;
};

// ---------------------------------------------------------------------------
std::complex<double> vector_meson_exchange::vector_propagator(int mu, int nu, double s, double zs)
{
  std::complex<double> result;
  result = - metric[mu][nu];
  result +=  exchange_momenta(mu, s, zs) * exchange_momenta(nu, s, zs) / mEx2;

  return result / (momentum_transfer(s,zs) - mEx2);
};

// ---------------------------------------------------------------------------
// Misc other functions

// ---------------------------------------------------------------------------
// Four-momentum of the exchanged meson.
// Simply the difference of the photon and axial 4-momenta
std::complex<double> vector_meson_exchange::exchange_momenta(int mu, double s, double zs)
{
  std::complex<double> qGamma_mu, qA_mu;
  qGamma_mu = kinematics->initial.component(mu, "beam", s, 1.);
  qA_mu = kinematics->final.component(mu, kinematics->vector_particle, s, zs);

  return  qGamma_mu - qA_mu;
};

// Mandelstam t momentum transfer
double vector_meson_exchange::momentum_transfer(double s, double zs)
{
  double t;
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp;
    temp = exchange_momenta(mu, s, zs);
    temp *= metric[mu][mu];
    temp *= exchange_momenta(mu, s, zs);

    t += real(temp);
  }

  return t;
};

// ---------------------------------------------------------------------------
std::complex<double> vector_meson_exchange::Sigma(int mu, int nu, int i, int j)
{
  std::complex<double> result;
  result = gamma_matrices[mu][i][j] * gamma_matrices[nu][i][j];
  result -= gamma_matrices[nu][i][j] * gamma_matrices[mu][i][j];

  return result / 2.;
};

// ---------------------------------------------------------------------------
// Four dimensional Levi-Civita symbol
double vector_meson_exchange::LeviCivita(int mu, int alpha, int beta, int gamma)
{
  // Error check
  if ((mu > 3 || alpha > 3 || beta > 3 || gamma > 3) || (mu < 0 || alpha < 0 || beta < 0 || gamma < 0))
  {
    std::cout << " \nLeviCivita: Error! Invalid argument recieved. Quitting... \n";
    exit(0);
  }

  // Return 0 if any are equal
  if (mu == alpha || mu == beta || mu == gamma || alpha == beta || alpha == gamma || beta == gamma)
  {
    return 0.;
  }

  // Else compare with strings
  std::string input = std::to_string(mu) + std::to_string(alpha) + std::to_string(beta) + std::to_string(gamma);
  std::vector<std::string>  even_permutations =
  {
    "0123",
    "0231",
    "0312",
    "1032",
    "1203",
    "1320",
    "2013",
    "2130",
    "2301",
    "3021",
    "3102",
    "3210"
  };

  for (int i = 0; i < 12; i++)
  {
    if (input == even_permutations[i])
    {
      return 1.;
    }
  }

  std::vector<std::string> odd_permutations =
  {
    "0132",
    "0213",
    "0321",
    "1023",
    "2103",
    "3120",
    "1230",
    "1302",
    "2031",
    "2310",
    "3012",
    "3201"
  };

  for (int i = 0; i < 12; i++)
  {
    if (input == odd_permutations[i])
    {
      return -1.;
    }
  }
};
