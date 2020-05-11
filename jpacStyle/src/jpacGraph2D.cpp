// The Wrapper for quick two-dimensional graphs using JPAC style
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "jpacGraph2D.hpp"

// -----------------------------------------------------------------------------
// Import the x, y, and z data to plot
void jpacGraph2D::SetData(std::vector<double> x, std::vector<double> y, std::vector<double> z)
{
  data = new TGraph2D(x.size(), &(x[0]), &(y[0]), &(z[0]));
};

// -----------------------------------------------------------------------------
// Clear saved data
void jpacGraph2D::ClearData()
{
  data = NULL;
};

// -----------------------------------------------------------------------------
// Add the J^{PAC} logo in black and white in top right corner
void jpacGraph2D::AddLogo()
{
  logo = new TLatex(.82, .89,  JPAC_BW.c_str());
  logo->SetNDC();
  logo->SetTextSize(2/30.);
  logo->SetTextAlign(32);
  logo->Draw();
};

// -----------------------------------------------------------------------------
// Plot to file
void jpacGraph2D::Plot(std::string filename)
{
  if (data == NULL)
  {
    std::cout << "\n";
    std::cout << "Warning! Trying to plot empty graph. Call to Plot() will be ignored... \n";
    std::cout << "\n";
    return;
  };

  // Force the canvas to be square
  // Also make sure to give enough room for the axes labels
  canvas->SetRightMargin(0.15);
  canvas->SetLeftMargin(0.15);
  canvas->SetBottomMargin(0.12);
  canvas->SetTopMargin(0.05);

  TAxis* xAxis = data->GetHistogram()->GetXaxis();
  xAxis->SetTitle(xLabel.c_str());
  xAxis->CenterTitle(true);
  if (xCustom == true)
  {
    xAxis->SetRangeUser(xlow, xhigh);
  }

  TAxis* yAxis = data->GetHistogram()->GetYaxis();
  yAxis->SetTitle(yLabel.c_str());
  yAxis->CenterTitle(true);
  if (yCustom == true)
  {
    yAxis->SetRangeUser(ylow, yhigh);
  }

  data->SetTitle("");
  data->Draw("COLZ");

  // Add the logo
  AddLogo();

  // Print to file
  canvas->Print(filename.c_str());
};
