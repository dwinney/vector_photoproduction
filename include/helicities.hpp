// Arrays containing helicity combinations for indexing.
// Moved to a seperate file to not clutter up reaction_kinematics.hpp
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _HELIC_COMBO_
#define _HELIC_COMBO_

#include <iostream>
#include <vector>
#include <array>

namespace jpacPhoto
{
    const std::vector< std::array<int, 4> > SPIN_ZERO_HELICITIES =
    {
    //  {  γ,  p,  S,  p'}
        {  1, -1,  0, -1},
        {  1, -1,  0,  1},
        {  1,  1,  0, -1},
        {  1,  1,  0,  1},
        { -1, -1,  0, -1},
        { -1, -1,  0,  1},
        { -1,  1,  0, -1},
        { -1,  1,  0,  1} 
    };

    const std::vector< std::array<int, 4> > SPIN_ONE_HELICITIES =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  1, -1},
        {  1, -1,  1,  1},
        {  1, -1,  0, -1},
        {  1, -1,  0,  1},
        {  1, -1, -1, -1},
        {  1, -1, -1,  1},
        {  1,  1,  1, -1},
        {  1,  1,  1,  1},
        {  1,  1,  0, -1},
        {  1,  1,  0,  1},
        {  1,  1, -1, -1},
        {  1,  1, -1,  1},
        { -1, -1,  1, -1},
        { -1, -1,  1,  1},
        { -1, -1,  0, -1},
        { -1, -1,  0,  1},
        { -1, -1, -1, -1},
        { -1, -1, -1,  1},
        { -1,  1,  1, -1},
        { -1,  1,  1,  1},
        { -1,  1,  0, -1},
        { -1,  1,  0,  1}, 
        { -1,  1, -1, -1},
        { -1,  1, -1,  1}
    };

    const std::vector< std::array<int, 4> > SPIN_TWO_HELICITIES =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  2, -1},
        {  1, -1,  2,  1},
        {  1, -1,  1, -1},
        {  1, -1,  1,  1},
        {  1, -1,  0, -1},
        {  1, -1,  0,  1},
        {  1, -1, -1, -1},
        {  1, -1, -1,  1},
        {  1, -1, -2, -1},
        {  1, -1, -2,  1},
        {  1,  1,  2, -1},
        {  1,  1,  2,  1},
        {  1,  1,  1, -1},
        {  1,  1,  1,  1},
        {  1,  1,  0, -1},
        {  1,  1,  0,  1},
        {  1,  1, -1, -1},
        {  1,  1, -1,  1},
        {  1,  1, -2, -1},
        {  1,  1, -2,  1},
        { -1, -1,  2, -1},
        { -1, -1,  2,  1},
        { -1, -1,  1, -1},
        { -1, -1,  1,  1},
        { -1, -1,  0, -1},
        { -1, -1,  0,  1},
        { -1, -1, -1, -1},
        { -1, -1, -1,  1},
        { -1, -1, -2, -1},
        { -1, -1, -2,  1},
        { -1,  1,  2, -1},
        { -1,  1,  2,  1},
        { -1,  1,  1, -1},
        { -1,  1,  1,  1},
        { -1,  1,  0, -1},
        { -1,  1,  0,  1}, 
        { -1,  1, -1, -1},
        { -1,  1, -1,  1},
        { -1,  1, -2, -1},
        { -1,  1, -2,  1}
    };


    inline std::vector<std::array<int, 4>> get_helicities(int J)
    {
        switch (J)
        {   
            case 0: return SPIN_ZERO_HELICITIES;
            case 1: return SPIN_ONE_HELICITIES;
            case 2: return SPIN_TWO_HELICITIES;
            default:
            {
                std::cout << "Error! Amplitudes for spin J = " << J << " not yet implemented. Quitting...\n";
                exit(0);
            }
        };
        
        return {};
    };
};

#endif