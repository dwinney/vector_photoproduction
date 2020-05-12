// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/vector_exchange.hpp"

// ---------------------------------------------------------------------------
// Assemble the helicity amplitude by contracting the lorentz indices
std::complex<double> vector_exchange::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  std::complex<double> result = 0.;
  if (FEYN == true)
  {
    for (int mu = 0; mu < 4; mu++)
    {
      for(int nu = 0; nu < 4; nu++)
      {
        std::complex<double> temp;
        temp = top_vertex(mu, lam_gam, lam_vec, s, zs);
        temp *= metric[mu];
        temp *= vector_propagator(mu, nu, s, zs);
        temp *= metric[nu];
        temp *= bottom_vertex(nu, lam_targ, lam_rec, s, zs);

        result += temp;
      }
    }
  }
  else
  {
    double t = kinematics->t_man(s, zs);
    std::complex<double> zt = kinematics->z_t(s, zs);

    int lam  = lam_gam - lam_vec;
    int lamp = (lam_targ - lam_rec) / 2.;

    // Product of residues
    result  = top_residue(lam, t);
    result *= bottom_residue(lamp, t);

    // Pole piece
    result /= t - mEx2;

    // angular function
    result *= wigner_d_int(1, lamp, lam, zt);
  }

  return result;
};

// ---------------------------------------------------------------------------
// Analytic residues
std::complex<double> vector_exchange::top_residue(int lam, double t)
{
  std::complex<double> result;
  switch (std::abs(lam))
  {
    case 0:
    {
      result = 1.;
      break;
    }
    case 1:
    {
      result = sqrt(xr * t) / kinematics->mVec;
      break;
    }
    case 2:
    {
      return gpGam * (t / kinematics->mVec2 - 1.);
    }
    default:
    {
      std::cout << "\nvector_exchange: invalid helicity flip lambda = " << lam << ". Quitting... \n";
      exit(0);
    }
  }

  std::complex<double> q = (t - kinematics->mVec2) / sqrt(4. * t * xr);
  return  result * q * gGam;
};

std::complex<double> vector_exchange::bottom_residue(int lamp, double t)
{
  std::complex<double> vector, tensor;
  switch (std::abs(lamp))
  {
    case 0:
    {
      vector =  1.;
      tensor = sqrt(xr * t) / (2. * mPro);
      break;
    }
    case 1:
    {
      vector = sqrt(2.) * sqrt(xr * t) / (2. * mPro);
      tensor = sqrt(2.);
      break;
    }
    case 2:
    {
      return 0.;
    }
    default:
    {
      std::cout << "\nreggeon_exchange: invalid helicity flip lambda^prime = " << lamp << ". Quitting... \n";
      exit(0);
    }
  }

  std::complex<double> result;
  result = gV * vector + gT * tensor * sqrt(xr * t) / (2. * mPro);
  result *= 2. * mPro;

  return result;
};

// ---------------------------------------------------------------------------
// Photon - Axial Vector - Vector vertex
// Feynman rules
std::complex<double> vector_exchange::top_vertex(int mu, int lam_gam, int lam_vec, double s, double zs)
{
  // Contract with LeviCivita
  std::complex<double> result = 0.;
  for (int alpha = 0; alpha < 4; alpha++)
  {
    for (int beta = 0; beta < 4; beta++)
    {
      for (int gamma = 0; gamma < 4; gamma++)
      {
        std::complex<double> temp;
        temp = levi_civita(mu, alpha, beta, gamma);
        temp *= metric[mu];
        temp *= kinematics->initial.component(alpha, "beam", s, 1.);
        temp *= kinematics->eps_gamma.component(beta, lam_gam, s, 1.);
        temp *= kinematics->eps_vec.component(gamma, lam_vec, s, zs);

        result += temp;
      }
    }
  }

  // Multiply by coupling
  return result * gGam;
};

// ---------------------------------------------------------------------------
// Nucleon - Nucleon - Vector vertex
std::complex<double> vector_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec, double s, double zs)
{
  // Vector coupling piece
  std::complex<double> vector = 0.;
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
  std::complex<double> tensor = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> sigma_q_ij = 0.;
      for (int nu = 0; nu < 4; nu++)
      {
        sigma_q_ij += sigma(mu, nu, i, j) * metric[nu] * exchange_momenta(nu, s, zs) / (2. * mPro);
      }

      std::complex<double> temp;
      temp = kinematics->recoil.adjoint_component(i, lam_rec, s, zs);
      temp *= sigma_q_ij;
      temp *= kinematics->target.component(j, lam_targ, s, 1.);

      tensor += temp;
    }
  }

  return gV * vector - gT * tensor;
};

// ---------------------------------------------------------------------------
// Four-momentum of the exchanged meson.
// Simply the difference of the photon and axial 4-momenta
std::complex<double> vector_exchange::exchange_momenta(int mu, double s, double zs)
{
  std::complex<double> qGamma_mu, qA_mu;
  qGamma_mu = kinematics->initial.component(mu, "beam", s, 1.);
  qA_mu = kinematics->final.component(mu, kinematics->vector_particle, s, zs);

  return (qGamma_mu - qA_mu);
};

// ---------------------------------------------------------------------------
// Propagator of a massive spin-one particle
std::complex<double> vector_exchange::vector_propagator(int mu, int nu, double s, double zs)
{
  std::complex<double> result;
  result = exchange_momenta(mu, s, zs) * exchange_momenta(nu, s, zs) / mEx2;

  if (mu == nu)
  {
    result -= metric[mu];
  }

  // pole piece (zero width here)
  result /= kinematics->t_man(s,zs) - mEx2;

  return result;
};
