// The OFFICIAL JPAC graphing style for ROOT
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#include "jpacPlotter.hpp"

// -----------------------------------------------------------------------------
// Functions to italicize and bold
std::string ROOT_italics(std::string in)
{
  return "#font[12]{" + in + "}";
};

std::string ROOT_bold(std::string in)
{
  return "#font[22]{" + in + "}";
};

std::string ROOT_bold_italics(std::string in)
{
  return "#font[32]{" + in + "}";
};

// -----------------------------------------------------------------------------
// JPAC Color Palette
Int_t jpacPlotter::kjpacBlue = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacBlue = new TColor(kjpacBlue, 0.12156862745098039, 0.4666666666666667, 0.7058823529411765);

Int_t jpacPlotter::kjpacRed = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacRed = new TColor(kjpacRed, 0.8392156862745098, 0.15294117647058825, 0.1568627450980392);

Int_t jpacPlotter::kjpacGreen = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacGreen = new TColor(kjpacGreen, 0.0, 0.6196078431372549, 0.45098039215686275);

Int_t jpacPlotter::kjpacOrange = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacOrange = new TColor(kjpacOrange, 0.8823529411764706, 0.4980392156862745, 0.054901960784313725);

Int_t jpacPlotter::kjpacPurple = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacPurple = new TColor(kjpacPurple, 0.5803921568627451, 0.403921568627451, 0.7411764705882353);

Int_t jpacPlotter::kjpacBrown = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacBrown = new TColor(kjpacBrown, 0.5490196078431373, 0.33725490196078434, 0.29411764705882354);

Int_t jpacPlotter::kjpacPink = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacPink = new TColor(kjpacPink, 0.8901960784313725, 0.4666666666666667, 0.7607843137254902);

Int_t jpacPlotter::kjpacGold = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacGold = new TColor(kjpacGold, 0.7372549019607844, 0.7411764705882353, 0.13333333333333333);

Int_t jpacPlotter::kjpacAqua = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacAqua = new TColor(kjpacAqua, 0.09019607843137255, 0.7450980392156863, 0.8117647058823529);

Int_t jpacPlotter::kjpacGrey = TColor::GetFreeColorIndex();
TColor * jpacPlotter::jpacGrey = new TColor(kjpacGrey, 0.4980392156862745, 0.4980392156862745, 0.4980392156862745);

std::vector<Int_t> jpacPlotter::jpacColors = {kjpacBlue, kjpacRed, kjpacGreen,
                                kjpacOrange, kjpacPurple, kjpacBrown,
                                kjpacPink, kjpacGold, kjpacAqua, kjpacGrey};

std::string jpacPlotter::JPAC_BW = "#scale[1.3]{#font[32]{J}^{#scale[0.8]{#font[32]{PAC}}}}";

std::string jpacPlotter::JPAC = "#scale[1.3]{#font[32]{#color[" + std::to_string(kjpacBlue) + "]{J}}"
                            + "^{#scale[0.8]{#font[32]{" + "#color[" + std::to_string(kjpacBlue) + "]{P}"
                            + "#color[" + std::to_string(kjpacRed) + "]{A}"
                            + "#color[" + std::to_string(kjpacBlue) + "]{C}"
                            + "}}}}";

std::string jpacPlotter::JPAC_BIG = "#scale[1.8]{#font[32]{#color[" + std::to_string(kjpacBlue) + "]{J}"
                          + "^{#scale[.7]{" + "#color[" + std::to_string(kjpacBlue) + "]{P}"
                          + "#color[" + std::to_string(kjpacRed) + "]{A}"
                          + "#color[" + std::to_string(kjpacBlue) + "]{C}"
                          + "}}}}";
