// Vector production via a one-loop box diagram
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "one_loop/box_amplitude.hpp"

// ---------------------------------------------------------------------------
// Evaluate the entire loop diagram in terms of external helicites and invariant mass / momentum transfer of the gamma p -> jpsi p process
std::complex<double> jpacPhoto::box_amplitude::helicity_amplitude(std::array<int,4> helicities, double s, double t)
{
    // Store external helicites
    _helicities = helicities;

    // Store the invariant energies to avoid having to pass them around 
    _s = s; _t = t, _theta = _kinematics->theta_s(s, t);

    // Pass 
    _discontinuity->set_externals(helicities, _s, _theta);

    // Integrate from threshold to cutoff
    double val[2], err[2];
    double min[1] = {_s_thr + EPS};
    double max[1] = {_s_cut};
    hcubature(2, wrapped_integrand, _discontinuity, 1, min, max, 2E4, 0, 1e-3, ERROR_INDIVIDUAL, val, err);

    return (val[0] + XI * val[1]) / PI;
};

// Static wrapper to pass all the necessary data to the integrater
int jpacPhoto::box_amplitude::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
{
    box_discontinuity * integrand = (box_discontinuity *) fdata;

    double sp = in[0];
    double s  = integrand->s();

    std::complex<double> result;
    result  = integrand->eval(sp);
    result /= (sp - s - IEPS);

    fval[0] = std::real(result);
    fval[1] = std::imag(result);

    return 0;
};

