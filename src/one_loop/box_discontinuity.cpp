// Discontinuity across unitarity cut for vector production via a one-loop box diagram
// Used as a container class in integration processes
// 
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "one_loop/box_discontinuity.hpp"

// ---------------------------------------------------------------------------
// Evaluate the integrand of the discontinuity as the product of the two
// sub amplitudes in terms of their respective t's
std::complex<double> jpacPhoto::disc_integrand::eval(double theta_gam, double phi_gam)
{
    // Check that the two amplitudes line up
    std::array<int,2> _jp_left  = _initialAmp->_kinematics->_jp;
    std::array<int,2> _jp_right = _finalAmp->_kinematics->_jp;

    if (_jp_left != _jp_right) 
    {
        std::cout << std::left;
        std::cout << std::setw(40) << "box_amplitude: Intermediate state between sub-amplitudes dont match!" << std::endl;
        std::cout << std::setw(20) << _initialAmp->_identifier << ": \t (" << _jp_left[0]  << ", " << _jp_left[1]  << ")\n";
        std::cout << std::setw(20) << _finalAmp->_identifier   << ": \t (" << _jp_right[0] << ", " << _jp_right[1] << ")\n";
        std::cout << "Returning 0!" << std::endl;

        return 0.;
    };

    std::vector< std::array<int,4> > intermediate_helicities = get_helicities(_jp_left[0]);

    // calculate the sub-process momentum transfers
    double t_gam        = _initialAmp->_kinematics->t_man(_s, theta_gam);

    double costheta_vec = cos(_theta) * cos(theta_gam) + sin(_theta) * sin(theta_gam) * cos(phi_gam);
    double theta_vec    = TMath::ACos(costheta_vec);
    double t_vec        = _finalAmp->_kinematics->t_man(_s, theta_gam);

    // 
    std::complex<double> phase_space;
    phase_space  = _initialAmp->_kinematics->_final_state->momentum(_s);
    phase_space *= 2. / sqrt(_s);

    // Sum over intermediate helicities 
    std::complex<double> result = 0.;
    for (int i = 0; i < 4*_jp_left[0]+2; i++)
    {
        int lam_meson  = intermediate_helicities[i][2];
        int lam_baryon = intermediate_helicities[i][3];

        std::complex<double> temp;
        temp  = _initialAmp->helicity_amplitude({_lam_gam, _lam_tar, lam_meson, lam_baryon}, _s, t_gam);
        temp *= _finalAmp->helicity_amplitude(  {_lam_vec, _lam_rec, lam_meson, lam_baryon}, _s, t_vec);
        
        result += temp;
    };
    
    return result * phase_space / (64. * PI*PI);
};

// ---------------------------------------------------------------------------
// Evaluate the product of sub-amplitudes integrating over intermediate phase-space
std::complex<double> jpacPhoto::box_discontinuity::eval(double s)
{
    // Pass stuff to the integrand
    _integrand->set_externals(_external_helicities, s, _external_theta);

    // Integrate over theta_gamma = [0, pi] and phi = [0, 2pi]
    double val[2], err[2];
    double min[2] = {0., 0.};
    double max[2] = {PI, 2. * PI};
    pcubature(2, wrapped_integrand, _integrand, 2, min, max, 2E5, 0, 1e-3, ERROR_INDIVIDUAL, val, err);
    
    std::complex<double> result = val[0] + XI * val[1];
    // debug(result);
    return result;
};

// Static wrapper to pass all the necessary data to the integrater
int jpacPhoto::box_discontinuity::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
{
    double theta_gam = in[0], phi_gam = in[1];
    disc_integrand * integrand = (disc_integrand *) fdata;

    std::complex<double> result;
    result  = integrand->eval(theta_gam, phi_gam);
    result *= sin(theta_gam); // Jacobian

    fval[0] = std::real(result);
    fval[1] = std::imag(result);

    return 0;
};

