// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/pomeron_exchange.hpp"

// ---------------------------------------------------------------------------
// Given a set of helicities for each particle, assemble the helicity amplitude by contracting Lorentz indicies
std::complex<double> jpacPhoto::pomeron_exchange::helicity_amplitude(std::array<int, 4> helicities, double xs, double xt)
{
    int lam_gam = helicities[0];
    int lam_targ = helicities[1];
    int lam_vec = helicities[2];
    int lam_rec = helicities[3];

    // Save energies 
    s = xs; t = xt; theta = kinematics->theta_s(xs, xt);

    std::complex<double> result = 0.;

    // IF using helicity conserving delta fuction model
    if (model == 1)
    {
        (lam_gam == lam_vec && lam_rec == lam_targ) ? (result = regge_factor()) : (result = 0.);
        return result;
    }

    // else contract indices
    for (int mu = 0; mu < 4; mu++)
    {
        std::complex<double> temp = 1.;
        temp *= regge_factor();
        temp *= top_vertex(mu, lam_gam, lam_vec);
        temp *= metric[mu];
        temp *= bottom_vertex(mu, lam_targ, lam_rec);

        result += temp;
    }

    return result;
};

// ---------------------------------------------------------------------------
// Bottom vertex coupling the target and recoil proton spinors to the vector pomeron
std::complex<double> jpacPhoto::pomeron_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec)
{
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            // Recoil oriented an angle theta + pi
            temp = kinematics->recoil->adjoint_component(i, lam_rec, s, theta + M_PI);

            // vector coupling
            temp *= gamma_matrices[mu][i][j];

            // target oriented in negative z direction
            temp *= kinematics->target->component(j, lam_targ, s, M_PI);

            result += temp;
        }
    }

    return result;
};

// ---------------------------------------------------------------------------
// Top vertex coupling the photon, pomeron, and vector meson.
std::complex<double> jpacPhoto::pomeron_exchange::top_vertex(int mu, int lam_gam, int lam_vec)
{
    std::complex<double> result = 0.;

    if (model == 0)
    {
        std::complex<double> sum1 = 0., sum2 = 0.;
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp1, temp2;

            // (q . eps_vec^*) eps_gam^mu
            temp1  = kinematics->initial->q(nu, s, 0.);
            temp1 *= metric[nu];
            temp1 *= kinematics->eps_vec->conjugate_component(nu, lam_vec, s, theta);
            sum1  += kinematics->eps_gamma->component(mu, lam_gam, s, 0.) * temp1;

            // (eps_vec^* . eps_gam) q^mu
            temp2  = kinematics->eps_gamma->component(nu, lam_gam, s, 0.);
            temp2 *= metric[nu];
            temp2 *= kinematics->eps_vec->conjugate_component(nu, lam_vec, s, theta);
            sum2  += kinematics->initial->q(mu, s, 0.) * temp2;
        }

        result = -sum1 + sum2;
    }
    else if (model == 2)
    {
        std::complex<double> sum1 = 0., sum2 = 0.;
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp1, temp2;

            // -2 * (q . eps_vec^*) eps_gam^mu
            temp1  = kinematics->initial->q(nu, s, 0.);
            temp1 *= metric[nu];
            temp1 *= kinematics->eps_vec->conjugate_component(nu, lam_vec, s, theta);
            sum1  += -2. * kinematics->eps_gamma->component(mu, lam_gam, s, 0.) * temp1;

            // (eps_vec . eps_gam) (q + q')^mu
            temp2  = kinematics->eps_vec->conjugate_component(nu, lam_vec, s, theta);
            temp2 *= metric[nu];
            temp2 *= kinematics->eps_gamma->component(nu, lam_gam, s, 0.);
            sum2  += (kinematics->initial->q(mu, s, 0.) + kinematics->final->q(mu, s, theta)) * temp2;
        }
      
        result = (sum1 + sum2);
    };

    return result;
};

// ---------------------------------------------------------------------------
// Usual Regge power law behavior, s^alpha(t) with an exponential fall from the forward direction
std::complex<double> jpacPhoto::pomeron_exchange::regge_factor()
{
    if (s < kinematics->sth())
    {
        std::cout << " \n pomeron_exchange: Trying to evaluate below threshold (sqrt(s) = " << sqrt(s) << ")! Quitting... \n";
        exit(0);
    }

    std::complex<double> result = 0.;
    
    switch (model)
    {
        case 0:
        {
            double t_min = kinematics->t_man(s, 0.); // t_min = t(theta = 0)
            result  = exp(b0 * (t - t_min));
            result *= pow(s - kinematics->sth(), pomeron_traj->eval(t));
            result *= xi * norm * e;
            result /= s;
            break;
        }
        case 1:
        {
            double t_min = kinematics->t_man(s, 0.); // t_min = t(theta = 0)
            result  = exp(b0 * (t - t_min));
            result *= pow(s - kinematics->sth(), pomeron_traj->eval(t));
            result *= xi * norm * e;
            break;
        }
        case 2:
        {
            double mX2 = kinematics->mX2;
            double th  = pow((kinematics->mT + kinematics->mR), 2.);

            double beta_0 = 2.;         // Pomeron - light quark coupling
            double beta_c = norm;       // Pomeron - charm quark coupling
            double mu2 = b0 * b0;       // cutoff parameter 
            double etaprime = std::real(pomeron_traj->slope());

            std::complex<double> F_t;
            F_t  = 3. * beta_0;
            F_t *= (th - 2.8*t);
            F_t /= (th - t) *  pow((1. - (t / 0.7)) , 2.);

            std::complex<double> G_p = -xi;
            G_p  *= pow(xr * etaprime * s, pomeron_traj->eval(t) - 1.);

            result  = - xi * 8. * beta_c * mu2 * G_p * F_t;
            result /= (mX2 - t) * (2.*mu2 + mX2 - t);
            break;
        }
        default: return 0.;
    }

    return result;
};
