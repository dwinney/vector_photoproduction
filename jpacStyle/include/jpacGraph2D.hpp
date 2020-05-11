// The Wrapper for quick two-dimensional graphs using JPAC style
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _JPAC_2D_
#define _JPAC_2D_

#include "jpacPlotter.hpp"

#include <TLegend.h>
#include <TGraph2D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <tuple>

class jpacGraph2D : public jpacPlotter
{
public:
  // Default constructor
  jpacGraph2D(){};

  // Parameterized constructor for effeciency
  jpacGraph2D(std::vector<double> x, std::vector<double> y, std::vector<double> z)
  {
    SetData(x, y, z);
  };

  ~jpacGraph2D()
  {
    delete data;
  };

  // Add in the data as vectors
  void SetData(std::vector<double> x, std::vector<double> y, std::vector<double> z);

  // Set up the Legend
  void SetLegend(double xx, double yy);

  // Plot to file
  void Plot(std::string filename);

private:
  // The actual ROOT object being called
  TGraph2D * data = NULL;

  // Position the logo in to top right
  void AddLogo();
};

#endif
