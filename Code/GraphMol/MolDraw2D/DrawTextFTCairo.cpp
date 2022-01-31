//
//  Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx).
//

#include <GraphMol/MolDraw2D/DrawTextFTCairo.h>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawTextFTCairo::DrawTextFTCairo(double max_fnt_sz, double min_fnt_sz,
                                 const std::string &font_file, cairo_t *dp_cr)
    : DrawTextFT(max_fnt_sz, min_fnt_sz, font_file), dp_cr_(dp_cr) {}

void DrawTextFTCairo::setCairoContext(cairo_t *cr) { dp_cr_ = cr; }

// ****************************************************************************
double DrawTextFTCairo::extractOutline() {
  cairo_set_source_rgba(dp_cr_, colour().r, colour().g, colour().b, colour().a);
  double adv = DrawTextFT::extractOutline();
  cairo_fill(dp_cr_);

  return adv;
}

// ****************************************************************************
int DrawTextFTCairo::MoveToFunctionImpl(const FT_Vector *to) {
  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);
  cairo_move_to(dp_cr_, dx, dy);

  return 0;
}

// ****************************************************************************
int DrawTextFTCairo::LineToFunctionImpl(const FT_Vector *to) {
  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);
  cairo_line_to(dp_cr_, dx, dy);

  return 0;
}

// ****************************************************************************
int DrawTextFTCairo::ConicToFunctionImpl(const FT_Vector *control,
                                         const FT_Vector *to) {
  // Cairo doesn't have a quadratic bezier function, we need to promote
  // it to a cubic using formula from
  // https://lists.cairographics.org/archives/cairo/2010-April/019691.html
  // and correcting the typo on the last line.

  double x0, y0;
  cairo_get_current_point(dp_cr_, &x0, &y0);

  double x1, y1;
  fontPosToDrawPos(control->x, control->y, x1, y1);

  double x2, y2;
  fontPosToDrawPos(to->x, to->y, x2, y2);

  cairo_curve_to(dp_cr_, (2.0 / 3.0) * x1 + (1.0 / 3.0) * x0,
                 (2.0 / 3.0) * y1 + (1.0 / 3.0) * y0,
                 (2.0 / 3.0) * x1 + (1.0 / 3.0) * x2,
                 (2.0 / 3.0) * y1 + (1.0 / 3.0) * y2, x2, y2);

  return 0;
}

// ****************************************************************************
int DrawTextFTCairo::CubicToFunctionImpl(const FT_Vector *controlOne,
                                         const FT_Vector *controlTwo,
                                         const FT_Vector *to) {
  double controlOneX, controlOneY;
  fontPosToDrawPos(controlOne->x, controlOne->y, controlOneX, controlOneY);
  double controlTwoX, controlTwoY;
  fontPosToDrawPos(controlTwo->x, controlTwo->y, controlTwoX, controlTwoY);

  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);

  cairo_curve_to(dp_cr_, controlOneX, controlOneY, controlTwoX, controlTwoY, dx,
                 dy);

  return 0;
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit
