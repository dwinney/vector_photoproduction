// The OFFICIAL JPAC graphing style for ROOT
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _JPAC_STYLE_
#define _JPAC_STYLE_

#include <iostream>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TError.h>
#include <TLatex.h>

class jpacPlotter
{
public:
  jpacPlotter()
  {
    SetStyle();
    Set2DPalette();
  };

  ~jpacPlotter()
  {
    delete jpacStyle;
    delete jpacBlue, jpacRed, jpacGreen;
    delete jpacOrange, jpacPurple, jpacBrown;
    delete jpacPink, jpacGold, jpacAqua, jpacGrey;
    delete canvas, logo;
  };

  static Int_t  kjpacBlue, kjpacRed, kjpacGreen,
                kjpacOrange, kjpacPurple, kjpacBrown,
                kjpacPink, kjpacGold, kjpacAqua, kjpacGrey;

  static std::vector<Int_t> jpacColors;

  static TColor *jpacBlue, *jpacRed, *jpacGreen,
                *jpacOrange, *jpacPurple, *jpacBrown,
                *jpacPink, *jpacGold, *jpacAqua, *jpacGrey;

  static std::string JPAC, JPAC_BW, JPAC_BIG;

  // Set all the style options
  TStyle* jpacStyle = new TStyle("jpacStyle", "JPAC Style");
  void SetStyle();

  // Stuff regarding the Logo
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);

  // For the 2D contour plots
  void Set2DPalette();

  // Axes Parameters
  bool xCustom = false, yCustom = false;
  std::string xLabel = "", yLabel = "";
  double xlow, xhigh, ylow, yhigh;

  // Set label and range for x axis
  void SetXaxis(std::string label, double low = 0., double high = 0.);
  void SetYaxis(std::string label, double low = 0., double high = 0.);

  // Plot all the saved entries and print to file given by filename
  // This will depend on specific implementation and needs to be overridden
  virtual void Plot(std::string filename) = 0;

  // Place the logo in the top right corner
  // Specific locations will depend on specific implementation
  TLatex* logo = NULL;
  virtual void AddLogo() = 0;
};

#endif
