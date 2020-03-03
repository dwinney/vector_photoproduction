// Class to sum up any number of generic amplitudes and build observables.
// Amplitudes are loaded up in a vector and summed incoherently
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _SUM_
#define _SUM_

#include "amplitudes/amplitude.hpp"

// ---------------------------------------------------------------------------
// The amplitude_sum class can take a vector of the above amplitude objects
// and add them together to get observables!
// ---------------------------------------------------------------------------

class amplitude_sum
{
private:
  // Store a vector of all the amplitudes you want to sum incoherently
  std::vector<amplitude*> amps;

public:
  // Empty constructor
  amplitude_sum()
  {};

  // Constructor with a vector already set up
  amplitude_sum(std::vector<amplitude*> vec)
  : amps(vec)
  {};

  // Add a new amplitude to the vector
  void add_amplitude(amplitude * new_amp)
  {
    amps.push_back(new_amp);
  };

  // Add all the members of an existing sum to a new sum
  void add_amplitude(amplitude_sum new_sum)
  {
    for (int i = 0; i < new_sum.amps.size(); i++)
    {
      amps.push_back(new_sum.amps[i]);
    }
  };

  // Evaluate the sum for given set of helicites, energy, and cos
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

};

#endif
