// Form of the invariant cross-section from a triple regge interaction.
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "inclusive/triple_regge.hpp"

double jpacPhoto::triple_regge::invariant_xsection(double s, double t, double M2)
{
    std:double result = 0.;

    for (int i = 0; i < _termsFF.size(); i++)
    {
        result += _termsFF[i]->eval(s, t, M2);
    }
    for (int i = 0; i < _termsJPAC.size(); i++)
    {
        result += _termsJPAC[i]->eval(s, t, M2);
    }

    return result;
};