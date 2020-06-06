// Charged axial-vector meson photoproduction proceeding through a pseudoscalar (pion) exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// ---------------------------------------------------------------------------

#include "amplitudes/pseudoscalar_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> pseudoscalar_exchange::helicity_amplitude(std::vector<int> helicities, double s, double zs)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  // Because its a scalar exchange we dont have any loose indices to contract
  std::complex<double> result;
  result  = top_vertex(lam_gam, lam_vec, s, zs);
  result *= scalar_propagator(s, zs);
  result *= bottom_vertex(lam_rec, lam_targ, s, zs);

  return result;
};

//------------------------------------------------------------------------------
// Pion form factors
double pseudoscalar_exchange::form_factor(double m, double s, double zs)
{
  double result;
  result  = (m*m - mEx2);
  result /= (m*m - kinematics->t_man(s,zs));
  
  return result;
};

//------------------------------------------------------------------------------
// Nucleon vertex
std::complex<double> pseudoscalar_exchange::bottom_vertex(double lam_rec, double lam_targ, double s, double zs)
{
  std::complex<double> result = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp;
      temp  = kinematics->recoil.adjoint_component(i, lam_rec, s, zs);
      temp *= gamma_5[i][j];
      temp *= kinematics->target.component(j, lam_targ, s, 1.);

      result += temp;
    }
  }

  result *= sqrt(2.) * gNN;
  result *= form_factor(LamPi, s, zs);

  return result;
};

//------------------------------------------------------------------------------
// Photon vertex
std::complex<double> pseudoscalar_exchange::top_vertex(double lam_gam, double lam_vec, double s, double zs)
{
  std::complex<double> term1 = 0., term2 = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    for (int nu = 0; nu < 4; nu++)
    {
      // (eps*_lam . eps_gam)(p . q)
      std::complex<double> temp1;
      temp1  = kinematics->eps_vec.conjugate_component(mu, lam_vec, s, zs);
      temp1 *= metric[mu];
      temp1 *= kinematics->eps_gamma.component(mu, lam_gam, s, 1.);
      temp1 *= kinematics->initial.component(nu, "beam", s, 1.);
      temp1 *= metric[nu];
      temp1 *= kinematics->final.component(nu, kinematics->vector_particle, s, zs);

      term1 += temp1;

      // (eps*_lam . p)(eps_gam . q)
      std::complex<double> temp2;
      temp2  = kinematics->eps_vec.conjugate_component(mu, lam_vec, s, zs);
      temp2 *= metric[mu];
      temp2 *= kinematics->initial.component(mu, "beam", s, 1.);
      temp2 *= kinematics->eps_gamma.component(nu, lam_gam, s, 1.);
      temp2 *= metric[nu];
      temp2 *= kinematics->final.component(nu, kinematics->vector_particle, s, zs);

      term2 += temp2;
    }
  }

  std::complex<double> result;
  result = term1 - term2;
  result *= gPsi / kinematics->mVec;
  result *= M_ALPHA / fJpsi;
  result *= form_factor(mJpsi, s, zs);

  return result;
};

//------------------------------------------------------------------------------
// Simple pole propagator
double pseudoscalar_exchange::scalar_propagator(double s, double zs)
{
  return 1. / (kinematics->t_man(s,zs) - mEx2);
};
