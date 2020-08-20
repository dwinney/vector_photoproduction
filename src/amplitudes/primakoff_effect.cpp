// Axial-vector meson photoproduction proceeding through the primakoff effect
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/primakoff_effect.hpp"

// ---------------------------------------------------------------------------
// Contract Lorentz indices
double jpacPhoto::primakoff_effect::differential_xsection(double xs, double xt)
{
    // update saved energies
    s = xs; t = xt; theta = kinematics->theta_s(xs, xt);

    // Form factor
    F_0 = form_factor(t);
    
    // Contract the Lorentz indices
    double sum = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp;
            temp  = top_tensor(mu, nu);
            temp *= metric[mu] * metric[nu];
            temp *= bottom_tensor(mu, nu);

            sum += real(temp);
        }
    }

    // photon propagator (squared)
    sum *= e * e / (t * t);

    // Normalization for dxs
    double norm = 1.;
    norm /= 64. * M_PI * s;
    norm /= real(pow(kinematics->initial->momentum(s), 2.));
    norm /= (2.56819E-6); // Convert from GeV^-2 -> nb

    return norm * sum;
};

// ---------------------------------------------------------------------------
// Inegrated total cross-section
// IN NANOBARN
double jpacPhoto::primakoff_effect::integrated_xsection(double s)
{
  auto F = [&](double t)
  {
    return differential_xsection(s, t);
  };

  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
  ROOT::Math::Functor1D wF(F);
  ig.SetFunction(wF);

  double t_min = kinematics->t_man(s, 0.);
  double t_max = kinematics->t_man(s, M_PI);

  return ig.Integral(t_max, t_min);
};

// ---------------------------------------------------------------------------
// Normalization
void jpacPhoto::primakoff_effect::calculate_norm()
{
    auto F = [&](double r)
    {
        return r * r * charge_distribution(r);
    };

    ROOT::Math::GSLIntegrator ig(   ROOT::Math::IntegrationOneDim::kADAPTIVE,
                                    ROOT::Math::Integration::kGAUSS61);

    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    // Integrate [0:inf]
    rho_0 = 1. / ig.IntegralUp(0.);
};

// ---------------------------------------------------------------------------
// Fourier transform the charge_distribution
double jpacPhoto::primakoff_effect::form_factor(double x)
{
    // momentum in the t channel
    double q = sqrt(x * (x - 4. * mPro2)) / (2. * mPro);

    auto dF = [&] (double r)
    {
        return r * r * sin(q * r) * charge_distribution(r);
    };

    ROOT::Math::GSLIntegrator ig(   ROOT::Math::IntegrationOneDim::kADAPTIVE,
                                    ROOT::Math::Integration::kGAUSS61);

    ROOT::Math::Functor1D wF(dF);
    ig.SetFunction(wF);
    
    return rho_0 * ig.IntegralUp(0.) / q;

};

// ---------------------------------------------------------------------------
// Top vertex
std::complex<double> jpacPhoto::primakoff_effect::top_tensor(int mu, int nu)
{
    std::complex<double> qq_term;
    qq_term  = t * Q2 / pow(kinematics->mVec, 4.);
    qq_term *= (1. - Q2 / kinematics->mVec2);
    qq_term *= kinematics->initial->q(mu, s, 0.);
    qq_term *= kinematics->initial->q(nu, s, 0.);

    std::complex<double> g_term;
    g_term  = (kinematics->mVec2 - Q2) * (kinematics->mVec2 + Q2) * (kinematics->mVec2 + Q2);
    g_term -= 2. * t * (kinematics->mVec2 + Q2) * (kinematics->mVec2 + Q2);
    g_term += t * t * (kinematics->mVec2 - Q2);
    g_term *= Q2 / (4. * pow(kinematics->mVec, 6.));
    g_term *= metric[mu];

    return g * g * (g_term);
};

// Bottom vertex
std::complex<double> jpacPhoto::primakoff_effect::bottom_tensor(int mu, int nu)
{
    std::complex<double> result;
    result  = 16. * M_PI;
    result *= kinematics->initial->p(mu, s, M_PI); + 0.5 * kinematics->initial->q(nu, s, 0.);
    result *= kinematics->initial->p(nu, s, M_PI); + 0.5 * kinematics->initial->q(mu, s, 0.);
    result *= Z * Z * 16. * mPro2 * mPro2 / (4. * M_PI);
    result *= F_0 * F_0;
    result /= (t - 4. * mPro2) * (t - 4. * mPro2);

    return result;
};
