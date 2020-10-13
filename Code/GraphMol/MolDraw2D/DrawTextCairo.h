//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx) on 29/04/2020.
//
// A concrete class derived from DrawText that uses the Cairo
// toy API to draw text onto a surface.

#ifndef RDKIT_DRAWTEXTCAIRO_H
#define RDKIT_DRAWTEXTCAIRO_H

#include <cairo.h>

#include <GraphMol/MolDraw2D/DrawText.h>

namespace RDKit {

// ****************************************************************************
class DrawTextCairo : public DrawText {
 public:
  DrawTextCairo(double max_fnt_sz, double min_fnt_sz, cairo_t *dp_cr);

#if 0
  void getStringSize(const std::string &label, double &label_width,
                     double &label_height) const override;
#endif
  void drawChar(char c, const Point2D &cds) override;

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

}  // namespace RDKit

#endif  // RDKIT_DRAWTEXTCAIRO_H
