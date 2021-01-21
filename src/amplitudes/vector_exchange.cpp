// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/vector_exchange.hpp"

// ---------------------------------------------------------------------------
// Assemble the helicity amplitude by contracting the lorentz indices
std::complex<double> jpacPhoto::vector_exchange::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    int lam_gam = helicities[0];
    int lam_targ = helicities[1];
    int lam_vec = helicities[2];
    int lam_rec = helicities[3];

    // Update the saved energies and angles
    _s = s; _t = t;
    _theta = _kinematics->theta_s(s, t);
    _zt = real(_kinematics->z_t(s, _theta));

    // Output
    std::complex<double> result;

    // if psuedo scalar or scalar production do covariant
    if (!(_kinematics->_jp[0] == 1 && _kinematics->_jp[1] == 1) || _useCovariant == true)
    {
        result = covariant_amplitude(helicities);
    }
    else
    {
        // NOTE THIS ONLY WORKS FOR UNPOLARIZED CROSS-SECTIONS
        // NEED TO WIGNER-ROTATE HELICITES TO S CHANNEL fOR POLARIZED 

        // TODO: ADD CROSSING-MATRICES
        int lam  = lam_gam - lam_vec;
        int lamp = (lam_targ - lam_rec) / 2.;

        if (abs(lam) == 2) return 0.; // double flip forbidden!

        // Product of residues  
        result  = top_residue(lam_gam, lam_vec);
        result *= bottom_residue(lam_targ, lam_rec);

        // Pole with d function residue if fixed spin
        if (_ifReggeized == false)
        {
            result *= wigner_d_int_cos(1, lam, lamp, _zt);
            result /= t - _mEx2;
        }
        // or regge propagator if reggeized
        else
        {
            result *= regge_propagator(1, lam, lamp);
        }
    }

    // add form factor if wanted
    result *= form_factor();    

    return result;
};

double jpacPhoto::vector_exchange::form_factor()
{
    switch (_useFormFactor)
    {
        // exponential form factor
        case 1: 
        {
            return exp((_t - _kinematics->t_man(_s, 0.)) / _cutoff*_cutoff);
        };

        // monopole form factor
        case 2:
        {
            return (_cutoff*_cutoff - _mEx2) / (_cutoff*_cutoff - _t); 
        };

        default:
        {
            return 1.;
        };
    }
};

// ---------------------------------------------------------------------------
// REGGE EVALUATION
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Analytic residues for Regge form

// Photon - Axial - Vector
std::complex<double> jpacPhoto::vector_exchange::top_residue(int lam_gam, int lam_vec)
{
    int lam = lam_gam - lam_vec;

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
            result = sqrt(XR * _t) / _kinematics->_mX;
            break;
        }
        default:
        {
            std::cout << "\nvector_exchange: invalid helicity flip lambda = " << lam << "!\n";
            return 0.;
        }
    }

    std::complex<double> q = (_t - _kinematics->_mX2) / sqrt(4. * _t * XR);
    return  XI * double(lam_gam) * result * q * _gGam;
};

// Nucleon - Nucleon - Vector
std::complex<double> jpacPhoto::vector_exchange::bottom_residue(int lam_targ, int lam_rec)
{
    // TODO: Explicit phases in terms of lam_targ and lam_rec instead of difference
    int lamp = (lam_targ - lam_rec) / 2.;

    std::complex<double> vector, tensor;
    switch (std::abs(lamp))
    {
        case 0:
        {
            vector =  1.;
            tensor = sqrt(XR * _t) / (2. * M_PROTON);
            break;
        }
        case 1:
        {
            vector = sqrt(2.) * sqrt(XR * _t) / (2. * M_PROTON);
            tensor = sqrt(2.);
            break;
        }
        case 2:
        {
            return 0.;
        }
        default:
        {
            std::cout << "\nreggeon_exchange: invalid helicity flip lambda^prime = " << lamp << "!\n";
            return 0.;
        }
    }

    std::complex<double> result;
    result = _gV * vector + _gT * tensor * sqrt(XR * _t) / (2. * M_PROTON);
    result *= 2. * M_PROTON;

    return result;
};

// ---------------------------------------------------------------------------
// Reggeon Propagator
std::complex<double> jpacPhoto::vector_exchange::regge_propagator(int j, int lam, int lamp)
{
    int M = std::max(std::abs(lam), std::abs(lamp));

    if (M > j)
    {
        return 0.;
    }

    std::complex<double> alpha_t = _alpha->eval(_t);

    // the gamma function causes problesm for large t so
    if (std::abs(alpha_t) > 30.)
    {
        return 0.;
    }
    else
    {
        std::complex<double> result;
        result  = wigner_leading_coeff(j, lam, lamp);
        result /= barrier_factor(j, M);
        result *= half_angle_factor(lam, lamp);

        result *= - _alpha->slope();
        result *= 0.5 * (double(_alpha->_signature) + exp(-XI * PI * alpha_t));
        result *= cgamma(1. - alpha_t);
        result *= pow(_s, alpha_t - double(M));

        return result;
    }
};

//------------------------------------------------------------------------------
// Half angle factors
std::complex<double> jpacPhoto::vector_exchange::half_angle_factor(int lam, int lamp)
{
    std::complex<double> sinhalf = sqrt((XR - _zt) / 2.);
    std::complex<double> coshalf = sqrt((XR + _zt) / 2.);

    std::complex<double> result;
    result  = pow(sinhalf, double(std::abs(lam - lamp)));
    result *= pow(coshalf, double(std::abs(lam + lamp)));

    return result;
};

//------------------------------------------------------------------------------
// Angular momentum barrier factor
std::complex<double> jpacPhoto::vector_exchange::barrier_factor(int j, int M)
{
    std::complex<double> q = (_t - _kinematics->_mX2) / sqrt(4. * _t * XR);
    std::complex<double> p = sqrt(XR * _t - 4.* M2_PROTON) / 2.;

    std::complex<double> result = pow(2. * p * q, double(j - M));

    return result;
};

// ---------------------------------------------------------------------------
// FEYNMAN EVALUATION
// ---------------------------------------------------------------------------

std::complex<double> jpacPhoto::vector_exchange::covariant_amplitude(std::array<int, 4> helicities)
{
    int lam_gam = helicities[0];
    int lam_targ = helicities[1];
    int lam_vec = helicities[2];
    int lam_rec = helicities[3];

    std::complex<double> result = 0.;

    // Need to contract the Lorentz indices
    for (int mu = 0; mu < 4; mu++)
    {
        for(int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp;
            temp  = top_vertex(mu, lam_gam, lam_vec);
            temp *= METRIC[mu];
            temp *= vector_propagator(mu, nu);
            temp *= METRIC[nu];
            temp *= bottom_vertex(nu, lam_targ, lam_rec);

            result += temp;
        }
    }

    return result;
};

// ---------------------------------------------------------------------------
// Photon - Axial Vector - Vector vertex
std::complex<double> jpacPhoto::vector_exchange::top_vertex(int mu, int lam_gam, int lam_vec)
{
    std::complex<double> result = 0.;

    // A-V-V coupling
    if (_kinematics->_jp[0]== 1 && _kinematics->_jp[1] == 1)
    {
        // Contract with LeviCivita
        for (int alpha = 0; alpha < 4; alpha++)
        {
            for (int beta = 0; beta < 4; beta++)
            {
                for (int gamma = 0; gamma < 4; gamma++)
                {
                    std::complex<double> temp;
                    temp = levi_civita(mu, alpha, beta, gamma);
                    if (std::abs(temp) < 0.001) continue;
                
                    temp *= METRIC[mu];
                    temp *= _kinematics->_initial_state->q(alpha, _s, 0.);
                    temp *= _kinematics->_eps_gamma->component(beta, lam_gam, _s, 0.);
                    temp *= _kinematics->_eps_vec->component(gamma, lam_vec, _s, _theta);

                    result += temp;
                }
            }
        }
    }

    // V-V-V coupling
    else if (_kinematics->_jp[0] == 1 && _kinematics->_jp[1] == -1)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp = XI;
            temp *= field_tensor(mu, nu, lam_gam);
            temp *= METRIC[nu];
            temp *= _kinematics->_eps_vec->component(nu, lam_vec, _s, _theta);

            result += temp;
        }
    }

    // S-V-V coupling
    else if (_kinematics->_jp[0]== 0 && _kinematics->_jp[1] == 1)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> term1, term2;

            // (k . q) eps_gamma^mu
            term1  = exchange_momenta(nu);
            term1 *= METRIC[nu];
            term1 *= _kinematics->_initial_state->q(nu, _s, 0.);
            term1 *= _kinematics->_eps_gamma->component(mu, lam_gam, _s, 0.);

            // (eps_gam . k) q^mu
            term2  = _kinematics->_eps_gamma->component(nu, lam_gam, _s, 0.);
            term2 *= METRIC[nu];
            term2 *= exchange_momenta(nu);
            term2 *= _kinematics->_initial_state->q(mu, _s, 0.);

            result += term1 - term2;
        }

        // Dimensionless coupling requires dividing by the mX
        result /= _kinematics->_mX;
    }

    // P-V-V coupling
    if (_kinematics->_jp[0]== 0 && _kinematics->_jp[1] == -1)
    {
        // Contract with LeviCivita
        for (int alpha = 0; alpha < 4; alpha++)
        {
            for (int beta = 0; beta < 4; beta++)
            {
                for (int gamma = 0; gamma < 4; gamma++)
                {
                    std::complex<double> temp;
                    temp = levi_civita(mu, alpha, beta, gamma);
                    if (std::abs(temp) < 0.001) continue;
                
                    temp *= field_tensor(alpha, beta, lam_gam);
                    temp *= _kinematics->_final_state->q(gamma, _s, _theta) - exchange_momenta(gamma);
                    result += temp;
                }
            }
        }
    }

    // Multiply by coupling
    return result * _gGam;
};

// ---------------------------------------------------------------------------
// Nucleon - Nucleon - Vector vertex
std::complex<double> jpacPhoto::vector_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec)
{
    // Vector coupling piece
    std::complex<double> vector = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            temp  = _kinematics->_recoil->adjoint_component(i, lam_rec, _s, _theta + PI); // theta_rec = theta + pi
            temp *= GAMMA[mu][i][j];
            temp *= _kinematics->_target->component(j, lam_targ, _s, PI); // theta_targ = pi

            vector += temp;
        }
    }

    // Tensor coupling piece
    std::complex<double> tensor = 0.;
    if (abs(_gT) > 0.001)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::complex<double> sigma_q_ij = 0.;
                for (int nu = 0; nu < 4; nu++)
                {
                sigma_q_ij += sigma(mu, nu, i, j) * METRIC[nu] * exchange_momenta(nu) / (2. * M_PROTON);
                }

                std::complex<double> temp;
                temp = _kinematics->_recoil->adjoint_component(i, lam_rec, _s, _theta + PI); // theta_rec = theta + pi
                temp *= sigma_q_ij;
                temp *= _kinematics->_target->component(j, lam_targ, _s, PI); // theta_targ = pi

                tensor += temp;
            }
        }
    }

    return _gV * vector - _gT * tensor;
};

std::complex<double> jpacPhoto::vector_exchange::field_tensor(int mu, int nu, int lambda)
{
    std::complex<double> result;
    result  = _kinematics->_initial_state->q(mu, _s, 0.) * _kinematics->_eps_gamma->component(nu, lambda, _s, 0.);
    result -= _kinematics->_initial_state->q(nu, _s, 0.) * _kinematics->_eps_gamma->component(mu, lambda, _s, 0.); 

    return result;
};

// ---------------------------------------------------------------------------
// Four-momentum of the exchanged meson.
// Simply the difference of the photon and axial 4-momenta
std::complex<double> jpacPhoto::vector_exchange::exchange_momenta(int mu)
{
    std::complex<double> qGamma_mu, qA_mu;
    qGamma_mu = _kinematics->_initial_state->q(mu, _s, 0.);
    qA_mu = _kinematics->_final_state->q(mu, _s, _theta);

    return (qGamma_mu - qA_mu);
};

// ---------------------------------------------------------------------------
// Propagator of a massive spin-one particle
std::complex<double> jpacPhoto::vector_exchange::vector_propagator(int mu, int nu)
{
    // q_mu q_nu / mEx2 - g_mu nu
    std::complex<double> result;
    result = exchange_momenta(mu) * exchange_momenta(nu) / _mEx2;

    if (mu == nu)
    {
        result -= METRIC[mu];
    }

    result /= _t - _mEx2;

    return result;
};