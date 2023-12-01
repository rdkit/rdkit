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
// A concrete class derived from DrawText that uses the Cairo
// toy API to draw text onto a surface.

#ifndef RDKIT_DRAWTEXTCAIRO_H
#define RDKIT_DRAWTEXTCAIRO_H

#include <cairo.h>

#include <GraphMol/MolDraw2D/DrawTextNotFT.h>

namespace RDKit {

class MolDraw2DCairo;
namespace MolDraw2D_detail {

// ****************************************************************************
class DrawTextCairo : public DrawTextNotFT {
 public:
  DrawTextCairo(double max_fnt_sz, double min_fnt_sz, cairo_t *dp_cr);
  DrawTextCairo(const DrawTextCairo &) = delete;
  DrawTextCairo(DrawTextCairo &&) = delete;
  DrawTextCairo &operator=(const DrawTextCairo &) = delete;
  DrawTextCairo &operator=(DrawTextCairo &&) = delete;
  void drawChar(char c, const Point2D &cds) override;
  void setCairoContext(cairo_t *cr);

 private:
  cairo_t *dp_cr_;

  // return a vector of StringRects, one for each char in text, with
  // super- and subscripts taken into account.  Sizes in pixel coords,
  // i.e. scaled by fontScale().
  void getStringRects(const std::string &text,
                      std::vector<std::shared_ptr<StringRect>> &rects,
                      std::vector<TextDrawType> &draw_modes,
                      std::vector<char> &draw_chars) const override;
};

}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // RDKIT_DRAWTEXTCAIRO_H
