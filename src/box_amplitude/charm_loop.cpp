// Charmonium production via a loop of open charm exchanges
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "box_amplitude/charm_loop.hpp"

// ---------------------------------------------------------------------------
// Evaluate the entire loop diagram in terms of external helicites and invariant mass / momentum transfer of the gamma p -> jpsi p process
std::complex<double> jpacPhoto::charm_loop::helicity_amplitude(std::array<int,4> helicities, double s, double t)
{
    // Store the invariant energies to avoid having to pass them around 
    _s = s; _t = t, _theta = _kinematics->theta_s(s, t);
    _disc.pass_params(helicities, _s, _theta, _eta);

    // destination for the integral result + associated errors
    double val[2], err[2];

    // bounds of integrtation for s, theta, and phi
    double s_thr  = pow(M_LAMBDAC + M_D, 2.) + EPS;
    double s_cut  = sqrt(_qmax*_qmax + M2_LAMBDAC) + sqrt(_qmax*_qmax + M2_D);
    double min[3] = {s_thr, 0,  0.};
    double max[3] = {s_cut, PI, 2.*PI};

    // Integrate the integral
    pcubature(2, wrapped_integrand, &_disc, 3, min, max, 2E3, 0, 1e-3, ERROR_INDIVIDUAL, val, err);
    std::complex<double> result = val[0] + XI * val[1];

    return result;
};

// ---------------------------------------------------------------------------
// The integrand for the triple integration (dispersion + intermediate phase space):
int jpacPhoto::charm_loop::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
{
    loop_discontinuity* disc = (loop_discontinuity *) fdata;

    double s  = disc->s();
    double sp         = in[0];
    double thetaGamma = in[1];
    double phiGamma   = in[2];

    std::complex<double> result = disc->eval(sp, thetaGamma, phiGamma) / (sp - s - IEPS);
    result /= (1. - in[0]) * (1. - in[0]); // jacobian of the sp map
    result *= sin(in[1]);
    
    fval[0] = std::real(result);
    fval[1] = std::imag(result);

    return 0;
};

// ---------------------------------------------------------------------------
// Functions for the auxiliary class loop_disconitnuity

// Evaluate the disconitnuity by multiplying tree diagrams on either side of unitarity cut
// This evaluates specifically for the intermediate D meson exchanges
std::complex<double> jpacPhoto::loop_discontinuity::eval_d(double sp, double thetaGam, double phiGam)
{
    // Calculate the sub-process momentum transfers
    _tGamma = kgamD.t_man(sp, thetaGam);

    double costhetaPsi = cos(_theta) * cos(thetaGam) + sin(_theta) * sin(thetaGam) * cos(phiGam);
    double thetaPsi = TMath::ACos(costhetaPsi);
    _tPsi   = kpsiD.t_man(sp, thetaPsi);

    std::complex<double> phase_space = kgamD._initial_state->momentum(_s) / (8.*PI*sqrt(_s));

    std::complex<double> result = 0.;

    // For D contribution only need to sum over intermediate Lambda_c helicities
    int LAMBDA_HELICITIES[2] = {-1, 1};
    for (int i = 0; i < 2; i++)
    {
        std::complex<double> ampGamma, ampPsi;
        std::array<int,4> heliGam{ {_lam_gam, _lam_tar, 0, LAMBDA_HELICITIES[i]} };
        std::array<int,4> heliPsi{ {_lam_psi, _lam_rec, 0, LAMBDA_HELICITIES[i]} }; 

        ampGamma = gamD.helicity_amplitude(heliGam, sp, _tGamma);
        ampPsi   = psiD.helicity_amplitude(heliPsi, sp, _tPsi);

        result += phase_space * ampGamma * ampPsi;
    };

    return result;
};