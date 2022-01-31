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

#include <GraphMol/MolDraw2D/DrawTextCairo.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

using namespace std;

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawTextCairo::DrawTextCairo(double max_fnt_sz, double min_fnt_sz,
                             cairo_t *dp_cr)
    : DrawTextNotFT(max_fnt_sz, min_fnt_sz), dp_cr_(dp_cr) {
  setCairoContext(dp_cr);
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void DrawTextCairo::drawChar(char c, const Point2D &cds) {
  PRECONDITION(dp_cr_, "no draw context");

  cairo_set_font_size(dp_cr_, fontSize());
  DrawColour col = colour();
  cairo_set_source_rgba(dp_cr_, col.r, col.g, col.b, col.a);

  char txt[2];
  txt[0] = c;
  txt[1] = 0;
  cairo_move_to(dp_cr_, cds.x, cds.y);
  cairo_show_text(dp_cr_, txt);
  cairo_stroke(dp_cr_);
}

void DrawTextCairo::setCairoContext(cairo_t *cr) {
  dp_cr_ = cr;
  if (dp_cr_) {
    cairo_select_font_face(dp_cr_, "sans", CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(dp_cr_, fontSize());
  }
}

// ****************************************************************************
void DrawTextCairo::getStringRects(const string &text,
                                   vector<shared_ptr<StringRect>> &rects,
                                   vector<TextDrawType> &draw_modes,
                                   vector<char> &draw_chars) const {
  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  double running_x = 0.0;
  char char_str[2];
  char_str[1] = 0;
  double max_y = 0.0;
  double full_fs = fontSize();
  auto *p_cr = dp_cr_;
  if (p_cr == nullptr) {
    // allocate an arbitrarily sized context temporarily
    cairo_surface_t *surf =
        cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 200, 200);
    p_cr = cairo_create(surf);
    cairo_surface_destroy(surf);  // dp_cr has a reference to this now;
    cairo_select_font_face(p_cr, "sans", CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
  }
  for (size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    draw_chars.push_back(text[i]);

    char_str[0] = text[i];
    cairo_text_extents_t extents;
    cairo_set_font_size(p_cr, selectScaleFactor(text[i], draw_mode) * full_fs);
    cairo_text_extents(p_cr, char_str, &extents);
    cairo_set_font_size(p_cr, full_fs);
    double twidth = extents.width;
    double theight = extents.height;
    Point2D offset(extents.x_bearing + twidth / 2.0, -extents.y_bearing / 2.0);
    Point2D g_centre(offset.x, -extents.y_bearing - theight / 2.0);
    rects.push_back(shared_ptr<StringRect>(
        new StringRect(offset, g_centre, twidth, theight)));
    rects.back()->trans_.x = running_x;
    draw_modes.push_back(draw_mode);
    running_x += extents.x_advance;
    max_y = max(max_y, -extents.y_bearing);
  }
  for (auto r : rects) {
    r->g_centre_.y = max_y - r->g_centre_.y;
    r->offset_.y = max_y / 2.0;
  }

  if (dp_cr_ == nullptr) {
    if (cairo_get_reference_count(p_cr) > 0) {
      cairo_destroy(p_cr);
    }
  }
  adjustStringRectsForSuperSubScript(draw_modes, rects);
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit
