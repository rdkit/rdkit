//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx) on 08/05/2020.
//

#include <GraphMol/MolDraw2D/DrawTextFTSVG.h>

namespace RDKit {

std::string DrawColourToSVG(const RDKit::DrawColour &col);

// ****************************************************************************
DrawTextFTSVG::DrawTextFTSVG(double max_fnt_sz, double min_fnt_sz,
                             const std::string &font_file, std::ostream &oss,
                             std::string &d_act_class)
    : DrawTextFT(max_fnt_sz, min_fnt_sz, font_file),
      oss_(oss),
      d_active_class_(d_act_class) {}

// ****************************************************************************
double DrawTextFTSVG::extractOutline() {
  std::string col = DrawColourToSVG(colour());

  oss_ << "<path ";
  if (!d_active_class_.empty()) {
    oss_ << " class='" << d_active_class_ << "'"
         << " d='";
  } else {
    oss_ << " d='";
  }

  double adv = DrawTextFT::extractOutline();
  oss_ << "' fill='" << col << "'/>" << std::endl;

  return adv;
}

// ****************************************************************************
int DrawTextFTSVG::MoveToFunctionImpl(const FT_Vector *to) {
  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);
  oss_ << "M " << dx << ' ' << dy << std::endl;

  return 0;
}

// ****************************************************************************
int DrawTextFTSVG::LineToFunctionImpl(const FT_Vector *to) {
  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);
  oss_ << "L " << dx << ' ' << dy << std::endl;

  return 0;
}

// ****************************************************************************
int DrawTextFTSVG::ConicToFunctionImpl(const FT_Vector *control,
                                       const FT_Vector *to) {
  double controlX, controlY;
  fontPosToDrawPos(control->x, control->y, controlX, controlY);

  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);

  oss_ << "Q " << controlX << ' ' << controlY << ", " << dx << ' ' << dy
       << std::endl;

  return 0;
}

// ****************************************************************************
int DrawTextFTSVG::CubicToFunctionImpl(const FT_Vector *controlOne,
                                       const FT_Vector *controlTwo,
                                       const FT_Vector *to) {
  double controlOneX, controlOneY;
  fontPosToDrawPos(controlOne->x, controlOne->y, controlOneX, controlOneY);
  double controlTwoX, controlTwoY;
  fontPosToDrawPos(controlTwo->x, controlTwo->y, controlTwoX, controlTwoY);

  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);

  oss_ << "C " << controlOneX << ' ' << controlOneY << ", " << controlTwoX
       << ' ' << controlTwoY << ", " << dx << ' ' << dy << std::endl;

  return 0;
}

}  // namespace RDKit
