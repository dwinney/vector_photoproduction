// The Wrapper for quick one-dimensional graphs using JPAC style
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "jpacGraph1D.hpp"

// -----------------------------------------------------------------------------
// Take in x and f(x) values as a vector and a legend entry
void jpacGraph1D::AddEntry(std::vector<double> xs, std::vector<double> fxs, std::string name)
{
  TGraph *g = new TGraph(xs.size(), &(xs[0]), &(fxs[0]));
  auto entry = std::make_tuple(g, name);
  entries.push_back(entry);
};

// -----------------------------------------------------------------------------
// Clear saved data
void jpacGraph1D::ClearData()
{
  for (int i = 0; i < entries.size(); i++)
  {
    delete std::get<0>(entries[i]);
  }
  entries.clear();
};

// -----------------------------------------------------------------------------
// Toggle legAdd which if false wont draw a legend at all
void jpacGraph1D::SetLegend(bool ifremove)
{
  legAdd = ifremove;
};

// Flip legCustom to true, indicate that we want manual placement of legend
void jpacGraph1D::SetLegend(double xx, double yy)
{
  xCord = xx; yCord = yy;
  legCustom = true;
};

// -----------------------------------------------------------------------------
// Add the J^{PAC} logo in appropriate colors at the top right of the plot
void jpacGraph1D::AddLogo()
{
  logo = new TLatex(.93, .89,  JPAC.c_str());
  logo->SetNDC();
  logo->SetTextSize(2/30.);
  logo->SetTextAlign(32);
  logo->Draw();
};

// -----------------------------------------------------------------------------
// Plot all the saved entries and print to file given by filename
void jpacGraph1D::Plot(std::string filename)
{
  if (entries.size() == 0)
  {
    std::cout << "\n";
    std::cout << "Error! Trying to plot empty graph. Call to Plot() will be ignored... \n";
    return;
  }
  if (entries.size() > 10)
  {
    std::cout << "\n";
    std::cout << "Warning! Number of curve greater than number of colors (9)! \n";
  }

  // Force the canvas to be square
  // Also make sure to give enough room for the axes labels
  canvas->SetTopMargin(0.05);
  canvas->SetRightMargin(0.05);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  canvas->SetFixedAspectRatio();

  // Set up the Legend
  if (legCustom == true)
  {
    legend = new TLegend(xCord, yCord, xCord + .3, yCord + .15);
  }
  else
  {
    legend = new TLegend(); // Automatic placement
  }

  legend->SetFillStyle(0);
  if (entries.size() > 5)
  {
    legend->SetNColumns(2); // Two column style if more than 5 entries
  }

  // Draw the first entry
  std::get<0>(entries[0])->UseCurrentStyle();
  std::get<0>(entries[0])->SetTitle("");
  std::get<0>(entries[0])->SetLineWidth(3);
  std::get<0>(entries[0])->SetLineColor(jpacColors[0]);

  // Set up the axes
  TAxis* xAxis = std::get<0>(entries[0])->GetXaxis();
  xAxis->SetTitle(xLabel.c_str());
  xAxis->CenterTitle(true);
  if (xCustom == true)
  {
    xAxis->SetRangeUser(xlow, xhigh);
  }

  TAxis* yAxis = std::get<0>(entries[0])->GetYaxis();
  yAxis->SetTitle(yLabel.c_str());
  yAxis->CenterTitle(true);
  if (yCustom == true)
  {
    yAxis->SetRangeUser(ylow, yhigh);
  }

  // Draw the first curve
  std::get<0>(entries[0])->Draw("AL");
  legend->AddEntry(std::get<0>(entries[0]), std::get<1>(entries[0]).c_str(), "l");

  // Add the logo
  AddLogo();

  // And any subsequent ones on the same canvas
  for (int i = 1; i < entries.size(); i++)
  {
    std::get<0>(entries[i])->SetLineWidth(3);
    std::get<0>(entries[i])->SetLineColor(jpacColors[i]);
    std::get<0>(entries[i])->Draw("same");

    legend->AddEntry(std::get<0>(entries[i]), std::get<1>(entries[i]).c_str(), "l");
  };

  if (legAdd == true)
  {
    legend->Draw();
  }
  canvas->Print(filename.c_str());
};
