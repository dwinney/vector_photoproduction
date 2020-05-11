// The Wrapper for quick one-dimensional graphs using JPAC style
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "jpacGraph1Dc.hpp"

// Utility functions that take in a vector of complex<doubles> and return vector<double>s of the same size
// containing only the real or imaginary parts.
//-----------------------------------------------------------------------------
std::vector<double> jpacGraph1Dc::vec_real(std::vector<std::complex<double>> fx)
{
  std::vector<double> result;
  for (int i = 0; i < fx.size(); i++)
  {
    result.push_back(real(fx[i]));
  }

  // Quick Error check
  if (result.size() != fx.size())
  {
    std::cout << "vec_real: ERROR Output and Input vector sizes dont match. Quitting... \n";
    std::exit(1);
  }

  return result;
};

std::vector<double> jpacGraph1Dc::vec_imag(std::vector<std::complex<double>> fx)
{
  std::vector<double> result;
  for (int i = 0; i < fx.size(); i++)
  {
    result.push_back(imag(fx[i]));
  }

  // Quick Error check
  if (result.size() != fx.size())
  {
    std::cout << "vec_imag: ERROR Output and Input vector sizes dont match. Quitting... \n";
    std::exit(1);
  }

  return result;
};

// -----------------------------------------------------------------------------
// Toggle legAdd which if false wont draw a legend at all
void jpacGraph1Dc::SetLegend(bool ifremove)
{
  legAdd = ifremove;
};

// Flip legCustom to true, indicate that we want manual placement of legend
void jpacGraph1Dc::SetLegend(double xx, double yy)
{
  xCord = xx; yCord = yy;
  legCustom = true;
};

// -----------------------------------------------------------------------------
// Set up the custom labels and ranges for the 2 y axes
void jpacGraph1Dc::SetYRealaxis(std::string label, double low, double high)
{
  yRLabel = label;
  if (std::abs(low) > 0.000001 || std::abs(high) > 0.000001)
  {
    yRlow = low; yRhigh = high;
    yRCustom = true;
  }
};

void jpacGraph1Dc::SetYImagaxis(std::string label, double low, double high)
{
  yILabel = label;
  if (std::abs(low) > 0.000001 || std::abs(high) > 0.000001)
  {
    yIlow = low; yIhigh = high;
    yICustom = true;
  }
};

// -----------------------------------------------------------------------------
// Take in x and f(x) values as a vector and a legend entry
void jpacGraph1Dc::AddEntry(std::vector<double> xs, std::vector<std::complex<double>> fxs, std::string name)
{
  std::vector<double> realx = vec_real(fxs);
  std::vector<double> imagx = vec_imag(fxs);

  TGraph* gReal = new TGraph(xs.size(), &(xs[0]), &(realx[0]));
  TGraph* gImag = new TGraph(xs.size(), &(xs[0]), &(imagx[0]));

  auto entry = std::make_tuple(gReal, gImag, name);
  entries.push_back(entry);
};

// -----------------------------------------------------------------------------
// Clear saved data
void jpacGraph1Dc::ClearData()
{
  for (int i = 0; i < entries.size(); i++)
  {
    delete std::get<0>(entries[i]);
    delete std::get<1>(entries[i]);
  }
  entries.clear();
};

// -----------------------------------------------------------------------------
// Add the J^{PAC} logo in appropriate colors at the top right of the plot
// Also add labels for Real and imaginary part
void jpacGraph1Dc::AddLogo()
{
  logo = new TLatex(.92, .87,  JPAC_BIG.c_str());
  logo->SetNDC();
  logo->SetTextSize(2/30.);
  logo->SetTextAlign(32);

  realTag->SetNDC();
  realTag->SetTextSize(2/30.);
  realTag->SetTextAlign(32);

  imagTag->SetNDC();
  imagTag->SetTextSize(2/30.);
  imagTag->SetTextAlign(32);

  canvas->cd(1);
  realTag->Draw();
  logo->Draw();
  canvas->cd(2);
  imagTag->Draw();
  logo->Draw();
};

// -----------------------------------------------------------------------------
// Plot all the saved entries and print to file given by filename
void jpacGraph1Dc::Plot(std::string filename)
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
    std::cout << "Warning! Number of curve greater than number of colors (10)! \n";
    std::cout << "\n";
  }
  if (yCustom == true || yLabel != "")
  {
    std::cout << "\n";
    std::cout << "Warning! Custom Y-axes for complex plots implemented with SetYRealaxis and SetYImagaxis! \n";
    std::cout << "Call to SetYaxis() will be ignored... \n";
    std::cout << "\n";
  }

  // Set up the margins of both canvases
  canvas->cd(1)->SetTopMargin(0.04);
  canvas->cd(1)->SetRightMargin(0.05);
  canvas->cd(1)->SetLeftMargin(0.09);
  canvas->cd(1)->SetBottomMargin(0.12);
  canvas->cd(1)->SetFixedAspectRatio();

  canvas->cd(2)->SetTopMargin(0.04);
  canvas->cd(2)->SetRightMargin(0.05);
  canvas->cd(2)->SetLeftMargin(0.09);
  canvas->cd(2)->SetBottomMargin(0.12);
  canvas->cd(2)->SetFixedAspectRatio();

  // Set up the Legend
  if (legCustom == true)
  {
    legend = new TLegend(xCord, yCord, xCord + .2, yCord + .12);
  }
  else
  {
    if (legAdd == true)
    {
      std::cout << "Warning! Automatic Legend placement may be wonky! \n";
      std::cout << "Recommned to specify coordinates using SetLegend(double, double) \n";
      std::cout << "\n";
    }

    legend = new TLegend(); // Automatic placement
  }
  //
  legend->SetFillStyle(0);
  if (entries.size() > 5)
  {
    legend->SetNColumns(2); // Two column style if more than 5 entries
    legend->SetColumnSeparation(.3);

  }

  // Plot the zeroth curve
  // Real part first
  canvas->cd(1);
  TGraph* gReal_0 = std::get<0>(entries[0]);
  gReal_0->SetTitle("");
  gReal_0->SetLineWidth(3);
  gReal_0->SetLineColor(jpacColors[0]);
  gReal_0->Draw("AL");

  // Add only the real parts to the Legend
  legend->AddEntry(gReal_0, std::get<2>(entries[0]).c_str(), "l");

  // Axes of Real Plot
  TAxis* xAxisReal = gReal_0->GetXaxis();
  xAxisReal->SetTitle(xLabel.c_str());
  xAxisReal->CenterTitle(true);
  if (xCustom == true)
  {
    xAxisReal->SetRangeUser(xlow, xhigh);
  }
  TAxis* yAxisReal = gReal_0->GetYaxis();
  yAxisReal->SetTitle(yRLabel.c_str());
  yAxisReal->CenterTitle(true);
  if (yRCustom == true)
  {
    yAxisReal->SetRangeUser(yRlow, yRhigh);
  }

  // Now the imaginary part on lower graph
  TGraph* gImag_0 = std::get<1>(entries[0]);
  canvas->cd(2);
  canvas->SetTopMargin(0.0);
  gImag_0->SetTitle("");
  gImag_0->SetLineWidth(3);
  gImag_0->SetLineColor(jpacColors[0]);
  gImag_0->Draw("AL");

  // Axes of Imag Plot
  TAxis* xAxisImag = gImag_0->GetXaxis();
  xAxisImag->SetTitle(xLabel.c_str());
  xAxisImag->CenterTitle(true);
  if (xCustom == true)
  {
    xAxisImag->SetRangeUser(xlow, xhigh);
  }
  TAxis* yAxisImag = gImag_0->GetYaxis();
  yAxisImag->SetTitle(yILabel.c_str());
  yAxisImag->CenterTitle(true);
  if (yICustom == true)
  {
    yAxisImag->SetRangeUser(yIlow, yIhigh);
  }

  // Repeat for each subsequent curve on top of the set up canvas
  for (int i = 1; i < entries.size(); i++)
  {
    canvas->cd(1);
    TGraph* gReal_i = std::get<0>(entries[i]);
    gReal_i->SetLineWidth(3);
    gReal_i->SetLineColor(jpacColors[i]);
    gReal_i->Draw("same");

    // Again only adding the real parts to legend
    legend->AddEntry(gReal_i, std::get<2>(entries[i]).c_str(), "l");

    canvas->cd(2);
    TGraph* gImag_i = std::get<1>(entries[i]);
    gImag_i->SetLineWidth(3);
    gImag_i->SetLineColor(jpacColors[i]);
    gImag_i->Draw("same");
  };

  // Add the logo to both plots!
  AddLogo();

  // Draw the legend
  if (legAdd == true)
  {
    canvas->cd(0);
    legend->Draw();
  }

  canvas->Modified();
  canvas->Print(filename.c_str());
};
