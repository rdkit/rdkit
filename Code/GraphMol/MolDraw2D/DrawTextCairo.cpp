//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx) on 29/04/2020.
//

#include <GraphMol/MolDraw2D/DrawTextCairo.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

namespace RDKit {

// ****************************************************************************
DrawTextCairo::DrawTextCairo(cairo_t *dp_cr) : DrawText(), dp_cr_(dp_cr) {
  cairo_select_font_face(dp_cr, "sans", CAIRO_FONT_SLANT_NORMAL,
                         CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(dp_cr, fontSize());
}

// ****************************************************************************
// using the current scale, work out the size of the label in molecule
// coordinates
void DrawTextCairo::getStringSize(const std::string &label, double &label_width,
                                  double &label_height) const {

  PRECONDITION(dp_cr_, "no draw context");
  label_width = 0.0;
  label_height = 0.0;

  DrawColour col = colour();
  cairo_set_source_rgb(dp_cr_, col.r, col.g, col.b);

  MolDraw2D::TextDrawType draw_mode = MolDraw2D::TextDrawNormal;
  // we have seen different behaviour on different OSes if the font is sized to
  // drawFontSize() / scale() which is what we really want.  Adjust at the end.
  cairo_set_font_size(dp_cr_, fontSize());

  bool had_a_super = false;

  char txt[2];
  txt[1] = 0;
  for (size_t i = 0, is = label.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == label[i] && setStringDrawMode(label, draw_mode, i)) {
      continue;
    }

    txt[0] = label[i];
    cairo_text_extents_t extents;
    cairo_text_extents(dp_cr_, txt, &extents);
    double twidth = extents.x_advance, theight = extents.height;
    label_height = std::max(label_height, theight);
    double char_width = twidth;
    if (MolDraw2D::TextDrawSubscript == draw_mode) {
      char_width *= 0.75;
    } else if (MolDraw2D::TextDrawSuperscript == draw_mode) {
      char_width *= 0.75;
      had_a_super = true;
    }
    label_width += char_width;
  }

  // subscript keeps its bottom in line with the bottom of the bit chars,
  // superscript goes above the original char top by a quarter
  if (had_a_super) {
    label_height *= 1.25;
  }
  label_height *= 1.2;  // empirical
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void DrawTextCairo::drawChar(char c, const Point2D &cds) {
  PRECONDITION(dp_cr_, "no draw context");
  char txt[2];
  txt[0] = c;
  txt[1] = 0;
  cairo_set_font_size(dp_cr_, fontSize());
  Point2D c1 = cds;
  cairo_move_to(dp_cr_, c1.x, c1.y);
  cairo_show_text(dp_cr_, txt);
  cairo_stroke(dp_cr_);
}

} // namespace RDKit
