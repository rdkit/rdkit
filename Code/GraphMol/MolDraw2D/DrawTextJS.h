//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: Greg Landrum
//
// A concrete class derived from DrawText that uses the JS Canvas
// to draw text onto a picture.
// Requires emscripten

#ifndef RDKIT_DRAWTEXTJS_H
#define RDKIT_DRAWTEXTJS_H

#include <iosfwd>

#include <emscripten.h>
#include <emscripten/val.h>
#include <GraphMol/MolDraw2D/DrawText.h>

namespace RDKit {

// ****************************************************************************

class DrawTextJS : public DrawText {
 public:
  DrawTextJS(double max_fnt_sz, double min_fnt_sz, emscripten::val &context);

  void drawChar(char c, const Point2D &cds) override;

 private:
  emscripten::val &context_;

  // fills a vector of StringRects, one for each char in text, with
  // super- and subscripts taken into account.  Sizes in pixel coords,
  // i.e. scaled by fontScale().
  void getStringRects(const std::string &text,
                      std::vector<std::shared_ptr<StringRect>> &rects,
                      std::vector<TextDrawType> &draw_modes,
                      std::vector<char> &draw_chars) const override;
};

}  // namespace RDKit

#endif  // RDKIT_DRAWTEXTSVG_H
