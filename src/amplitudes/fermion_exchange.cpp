// Spin-1/2 exchange ampltiude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/fermion_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> fermion_exchange::helicity_amplitude(std::vector<int> helicities, double s, double zs)
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
      std::complex<double> temp = 1.;
      temp *= top_vertex(i, lam_gam, lam_rec, s, zs);
      temp *= fermion_propagator(i, j, s, zs);
      temp *= bottom_vertex(j, lam_vec, lam_targ, s, zs);

      result += temp;
    }
  }

  return result;
};


//------------------------------------------------------------------------------
// Photon fermion fermion vertex
// (ubar epsilon-slashed)
std::complex<double> fermion_exchange::top_vertex(int i, int lam_gam, int lam_rec, double s, double zs)
{
  if (ScTOP == true)
  {
    // Scalar for testing purposes
    return gGam * kinematics->recoil.adjoint_component(i, lam_rec, s, zs);
  }

  std::complex<double> result = 0.;
  for (int k = 0; k < 4; k++)
  {
    std::complex<double> temp;
    temp  = kinematics->recoil.adjoint_component(k, lam_rec, s, zs);
    temp *= slashed_eps(k, i, lam_gam, kinematics->eps_gamma, false, s, 1.);

    result += temp;
  }

  return gGam * result;
};

//------------------------------------------------------------------------------
// Vector fermion fermion vertex
// (epsilon*-slashed u)
std::complex<double> fermion_exchange::bottom_vertex(int j, int lam_vec, int lam_targ, double s, double zs)
{
  if (ScBOT == true)
  {
    // Scalar for testing purposes
    return gVec * kinematics->target.component(j, lam_targ, s , 1.);
  }

  std::complex<double> result = 0.;
  for (int k = 0; k < 4; k++)
  {
    std::complex<double> temp;
    temp  = slashed_eps(j, k, lam_vec, kinematics->eps_vec, true, s, zs);
    temp *= kinematics->target.component(k, lam_targ, s , 1.);

    result += temp;
  }

  return gVec * result;
};

//------------------------------------------------------------------------------
double fermion_exchange::exchange_mass(double s, double zs)
{
    double result = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
      std::complex<double> temp;
      temp  = exchange_momentum(mu, s, zs);
      temp *= metric[mu];
      temp *= exchange_momentum(mu, s, zs);

      result += real(temp);
    }
    return result;
}

std::complex<double> fermion_exchange::exchange_momentum(int mu, double s, double zs)
{
  std::complex<double> qGamma_mu, qRec_mu;
  qGamma_mu   = kinematics->initial.component(mu, "beam", s, 1.);
  qRec_mu     = kinematics->final.component(mu, "recoil", s, zs);

  return (qGamma_mu - qRec_mu);
};

std::complex<double> fermion_exchange::slashed_exchange_momentum(int i, int j, double s, double zs)
{
  std::complex<double> result = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp;
    temp  = gamma_matrices[mu][i][j];
    temp *= metric[mu];
    temp *= exchange_momentum(mu, s, zs);

    result += temp;
  }

  return result;
};

//------------------------------------------------------------------------------
// Slashed polarization vectors
std::complex<double> fermion_exchange::slashed_eps(int i, int j, double lam, polarization_vector eps, bool STAR, double s, double zs)
{
  std::complex<double> result = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp;
    if (STAR == false)
    {
      temp  = eps.component(mu, lam, s, zs);
    }
    else
    {
      temp = eps.conjugate_component(mu, lam, s, zs);
    }
    temp *= metric[mu];
    temp *= gamma_matrices[mu][i][j];

    result += temp;
  }

  return result;
};


//------------------------------------------------------------------------------
std::complex<double> fermion_exchange::fermion_propagator(int i, int j, double s, double zs)
{
  std::complex<double> result;
  result = slashed_exchange_momentum(i, j, s, zs);

  if (i == j)
  {
    result += mEx2;
  }

  result /= exchange_mass(s, zs) - mEx2;

  return result;
};
