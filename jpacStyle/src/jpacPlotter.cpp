// The OFFICIAL JPAC graphing style for ROOT
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "jpacPlotter.hpp"

// -----------------------------------------------------------------------------
// All other noncolor settings
void jpacPlotter::SetStyle()
{
  gErrorIgnoreLevel = kWarning;

  // remove info box
  jpacStyle->SetOptStat(0);

  // Centre title
  jpacStyle->SetTitleAlign(22);
  jpacStyle->SetTitleX(.5);
  jpacStyle->SetTitleY(.95);
  jpacStyle->SetTitleBorderSize(0);

  // set background colors to white
  jpacStyle->SetFillColor(10);
  jpacStyle->SetFrameFillColor(10);
  jpacStyle->SetCanvasColor(10);
  jpacStyle->SetPadColor(10);
  jpacStyle->SetTitleFillColor(0);
  jpacStyle->SetStatColor(10);

  // Don't put a colored frame around the plots
  jpacStyle->SetFrameBorderMode(0);
  jpacStyle->SetCanvasBorderMode(0);
  jpacStyle->SetPadBorderMode(0);

  // No border on legends
  jpacStyle->SetLegendBorderSize(0);
  jpacStyle->SetLegendTextSize(0.03);

  // Axis titles
  jpacStyle->SetNdivisions(506, "xy");
  jpacStyle->SetTitleSize(.045, "xyz");
  jpacStyle->SetTitleOffset(1.1, "xyz");

  // More space for y-axis to avoid clashing with big numbers
  jpacStyle->SetTitleOffset(1.6, "y");

  // This applies the same settings to the overall plot title
  jpacStyle->SetTitleSize(.055, "");
  jpacStyle->SetTitleOffset(.8, "");

  // Axis labels (numbering)
  jpacStyle->SetLabelSize(.035, "xyz");
  jpacStyle->SetLabelOffset(.01, "xyz");

  // Thicker lines
  jpacStyle->SetFrameLineWidth(2);

  // Set the tick mark style
  jpacStyle->SetPadTickX(1);
  jpacStyle->SetPadTickY(1);

  const int kjpacFont = 132;
  jpacStyle->SetStatFont(kjpacFont);
  jpacStyle->SetLabelFont(kjpacFont, "xyz");
  jpacStyle->SetTitleFont(kjpacFont, "xyz");
  jpacStyle->SetTitleFont(kjpacFont, ""); // Apply same setting to plot titles
  jpacStyle->SetTextFont(kjpacFont);
  jpacStyle->SetLegendFont(kjpacFont);

  // Make it the global default style
  gROOT->SetStyle("jpacStyle");
  canvas->UseCurrentStyle();
};

// -----------------------------------------------------------------------------
// Set up the color map for 2D plots
void jpacPlotter::Set2DPalette()
{
  Int_t NCont = 999;
  Double_t Red[3]   = { 0.12156862745098039, 1.0,  0.8392156862745098};
  Double_t Green[3] = { 0.4666666666666667, 1.0, 0.15294117647058825};
  Double_t Blue[3]  = { 0.7058823529411765, 1.0, 0.1568627450980392};

  Double_t Stops[3] = { 0.0, 0.5, 1.0 };
  Int_t FI = TColor::CreateGradientColorTable(3, Stops, Red, Green, Blue, NCont);

  Int_t jpac3D[NCont];
  for (int i=0; i < NCont; i++)
  {
    jpac3D[i] = FI+i;
  }

  jpacStyle->SetPalette(NCont, jpac3D);
  jpacStyle->SetNumberContours(NCont);
};

// -----------------------------------------------------------------------------
// Set the axes
void jpacPlotter::SetXaxis(std::string label, double low, double high)
{
  xLabel = label;

  if (std::abs(low) > 0.000001 || std::abs(high) > 0.000001)
  {
    xlow = low; xhigh = high;
    xCustom = true;
  }
};

void jpacPlotter::SetYaxis(std::string label, double low, double high)
{
  yLabel = label;
  if (std::abs(low) > 0.000001 || std::abs(high) > 0.000001)
  {
    ylow = low; yhigh = high;
    yCustom = true;
  }
};
