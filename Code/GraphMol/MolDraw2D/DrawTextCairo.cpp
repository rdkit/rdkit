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

using namespace std;

namespace RDKit {

// ****************************************************************************
DrawTextCairo::DrawTextCairo(double max_fnt_sz, cairo_t *dp_cr)
    : DrawText(max_fnt_sz), dp_cr_(dp_cr) {
  cout << "DrawTextCairo  fontSize = " << fontSize() << endl;
  cairo_select_font_face(dp_cr, "sans", CAIRO_FONT_SLANT_NORMAL,
                         CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(dp_cr, fontSize());
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void DrawTextCairo::drawChar(char c, const Point2D &cds) {
  PRECONDITION(dp_cr_, "no draw context");

  cairo_set_font_size(dp_cr_, fontSize());
  DrawColour col = colour();
  cairo_set_source_rgb(dp_cr_, col.r, col.g, col.b);

  char txt[2];
  txt[0] = c;
  txt[1] = 0;
  Point2D c1 = cds;
  cairo_move_to(dp_cr_, c1.x, c1.y);
  cairo_show_text(dp_cr_, txt);
  cairo_stroke(dp_cr_);
}

// ****************************************************************************
void DrawTextCairo::getStringRects(const string &text,
                                   vector<shared_ptr<StringRect>> &rects,
                                   vector<TextDrawType> &draw_modes) const {

  cout << "DrawTextCairo::getStringRects" << endl;
  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  double running_x = 0.0;
  std::vector<char> to_draw;
  char char_str[2];
  char_str[1] = 0;

  cairo_set_font_size(dp_cr_, fontSize());
  for(size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    to_draw.push_back(text[i]);
    char_str[0] = text[i];
    cairo_text_extents_t extents;
    cairo_text_extents(dp_cr_, char_str, &extents);
    cout << extents.width << " by " << extents.height << " and "
         << extents.x_advance << endl;
    double w_mult = 1.0;
    if(draw_mode == TextDrawType::TextDrawSuperscript
       || draw_mode == TextDrawType::TextDrawSubscript) {
      w_mult  = 0.75;
    }
    double twidth = w_mult * extents.x_advance;
    double theight = extents.height;
    Point2D centre(running_x + twidth / 2, theight / 2);
    rects.push_back(shared_ptr<StringRect>(new StringRect(centre, twidth, theight)));
    draw_modes.push_back(draw_mode);
    running_x += w_mult * extents.x_advance;
  }

  adjustStringRectsForSuperSubScript(draw_modes, to_draw, rects);

}

} // namespace RDKit
