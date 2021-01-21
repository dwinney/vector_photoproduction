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
std::complex<double> jpacPhoto::pseudoscalar_exchange::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    int lam_gam = helicities[0];
    int lam_targ = helicities[1];
    int lam_vec = helicities[2];
    int lam_rec = helicities[3];

    // Store the invariant energies to avoid having to pass them around 
    _s = s; _t = t, _theta = _kinematics->theta_s(s, t);

    std::complex<double> result;

    if (_useFourVecs == true)
    {
        // Because its a scalar exchange we dont have any loose indices to contract
        result  = top_vertex(lam_gam, lam_vec);
        result *= scalar_propagator();
        result *= bottom_vertex(lam_targ, lam_rec);
    }
    else
    {
        if (lam_vec != lam_gam || lam_targ != lam_rec) 
        {
            return 0.; 
        }
        else
        {
            result  = sqrt(2.) * _gNN;
            result *= _gGamma / _kinematics->_mX;
            result *= sqrt(XR * _t) / 2.;
            result *= (_kinematics->_mX2 - _t);
            result *= scalar_propagator();
        }
    }

    // Multiply by the optional expontial form factor
    if (_useFF == true)
    {
        double tprime = _t - _kinematics->t_man(s, 0.);
        result *= exp(_b * tprime);
    }

    return result;
};

//------------------------------------------------------------------------------
// Nucleon vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_vertex(double lam_targ, double lam_rec)
{
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            // ubar(recoil) * gamma_5 * u(target)
            std::complex<double> temp;
            temp  = _kinematics->_recoil->adjoint_component(i, lam_rec, _s, _theta + PI); // theta_recoil = theta + pi
            temp *= GAMMA_5[i][j];
            temp *= _kinematics->_target->component(j, lam_targ, _s, PI); // theta_target = pi

            result += temp;
        }
    }

    // Sqrt(2) from isospin considering a charged pion field
    // remove the Sqrt(2) if considering a neutral pion exchange
    result *= sqrt(2.) * _gNN;

    return result;
};

//------------------------------------------------------------------------------
// Photon vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::top_vertex(double lam_gam, double lam_vec)
{
    std::complex<double> result = 0.;

    // A - V - P
    if (_kinematics->_jp[0] == 1 && _kinematics->_jp[1] == 1)
    {
         std::complex<double> term1 = 0., term2 = 0.;
        for (int mu = 0; mu < 4; mu++)
        {
            for (int nu = 0; nu < 4; nu++)
            {
                // (eps*_lam . eps_gam)(q_vec . q_gam)
                std::complex<double> temp1;
                temp1  = _kinematics->_eps_vec->conjugate_component(mu, lam_vec, _s, _theta);
                temp1 *= METRIC[mu];
                temp1 *= _kinematics->_eps_gamma->component(mu, lam_gam, _s, 0.);
                temp1 *= _kinematics->_initial_state->q(nu, _s, 0.);
                temp1 *= METRIC[nu];
                temp1 *= _kinematics->_final_state->q(nu, _s, _theta);

                term1 += temp1;

                // (eps*_lam . q_gam)(eps_gam . q_vec)
                std::complex<double> temp2;
                temp2  = _kinematics->_eps_vec->conjugate_component(mu, lam_vec, _s, _theta);
                temp2 *= METRIC[mu];
                temp2 *= _kinematics->_initial_state->q(mu, _s, 0.);
                temp2 *= _kinematics->_eps_gamma->component(nu, lam_gam, _s, 0.);
                temp2 *= METRIC[nu];
                temp2 *= _kinematics->_final_state->q(nu, _s, _theta);

                term2 += temp2;

                result = (temp1 - temp2) / _kinematics->_mX;
            }
        }
    }

    // V - V - P
    if (_kinematics->_jp[0] == 1 && _kinematics->_jp[1] == -1)
    {
        // Contract with LeviCivita
        for (int mu = 0; mu < 4; mu ++)
        {
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    for (int gamma = 0; gamma < 4; gamma++)
                    {
                        std::complex<double> temp;
                        temp = levi_civita(mu, alpha, beta, gamma);
                        if (std::abs(temp) < 0.001) continue;
                        temp *= _kinematics->_eps_vec->conjugate_component(mu, lam_vec, _s, _theta);
                        temp *= _kinematics->_eps_gamma->field_tensor(alpha, beta, lam_gam, _s, 0.);
                        temp *= _kinematics->_final_state->q(gamma, _s, _theta) - _kinematics->t_exchange_momentum(gamma, _s, _theta);
                        result += temp;
                    }
                }
            }
        }
    }

    return _gGamma * result;
};

//------------------------------------------------------------------------------
// Simple pole propagator
std::complex<double> jpacPhoto::pseudoscalar_exchange::scalar_propagator()
{
    if (_reggeized == false)
    {
        return 1. / (_t - _mEx2);
    }
    else
    {
        std::complex<double> alpha_t = _alpha->eval(_t);

        if (std::abs(alpha_t) > 20.) return 0.;

        // Else use the regge propagator
        std::complex<double> result = 1.;
        result  = - _alpha->slope();
        result *= 0.5 * (double(_alpha->_signature) +  exp(-XI * PI * alpha_t));
        result *= cgamma(0. - alpha_t);
        result *= pow(_s, alpha_t);
        return result;
    }
};
