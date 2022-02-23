//
//  Copyright (C) 2020-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: Greg Landrum
// A concrete class derived from DrawText that uses the JS Canvas
// to draw text onto a picture.
#ifdef __EMSCRIPTEN__
#include <sstream>
#include <boost/algorithm/string.hpp>

#include <GraphMol/MolDraw2D/DrawTextJS.h>
#include <GraphMol/MolDraw2D/MolDraw2DJS.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <boost/format.hpp>

namespace RDKit {

std::string DrawColourToSVG(const RDKit::DrawColour &col);

namespace MolDraw2D_detail {
// ****************************************************************************
DrawTextJS::DrawTextJS(double max_fnt_sz, double min_fnt_sz,
                       emscripten::val &context)
    : DrawTextNotFT(max_fnt_sz, min_fnt_sz), context_(context) {}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void DrawTextJS::drawChar(char c, const Point2D &cds) {
  std::string txt(1, c);
  std::string font = (boost::format("%fpx sans-serif") % fontSize()).str();
  std::string col = DrawColourToSVG(colour());
  context_.set("font", font);
  context_.set("fillStyle", col);
  context_.call<void>("fillText", txt, cds.x, cds.y);
}

// ****************************************************************************
void DrawTextJS::getStringRects(const std::string &text,
                                std::vector<std::shared_ptr<StringRect>> &rects,
                                std::vector<TextDrawType> &draw_modes,
                                std::vector<char> &draw_chars) const {
  double running_x = 0.0;
  double act_font_size = fontSize();
  double char_height;
  double max_width = 0.0;
  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  for (size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    draw_modes.push_back(draw_mode);
    draw_chars.push_back(text[i]);

    max_width = std::max(
        max_width,
        static_cast<double>(MolDraw2D_detail::char_widths[(int)text[i]]));
  }

  for (size_t i = 0; i < draw_chars.size(); ++i) {
    double char_width =
        0.6 * act_font_size *
        static_cast<double>(MolDraw2D_detail::char_widths[(int)draw_chars[i]]) /
        max_width;
    // Absent a proper set of font metrics (we don't know what font we'll be
    // using, for starters) this is something of an empirical bodge.
    if (draw_chars[i] == '+') {
      char_height = 0.6 * act_font_size;
    } else if (draw_chars[i] == '-') {
      char_height = 0.4 * act_font_size;
    } else {
      char_height = 0.8 * act_font_size;
    }
    double cscale = selectScaleFactor(draw_chars[i], draw_modes[i]);
    char_height *= cscale;
    char_width *= cscale;
    Point2D offset(char_width / 2, char_height / 2);
    if (draw_chars[i] == '+' || draw_chars[i] == '-') {
      offset.y /= 2.0;
    }
    Point2D g_centre(char_width / 2, char_height / 2);
    rects.push_back(std::shared_ptr<StringRect>(
        new StringRect(offset, g_centre, char_width, char_height)));
    rects.back()->trans_.x += running_x;
    // empirical spacing.
    if (draw_modes[i] != TextDrawType::TextDrawNormal) {
      running_x += char_width * 1.05;
    } else {
      running_x += char_width * 1.15;
    }
  }
  for (auto r : rects) {
    r->g_centre_.y = act_font_size - r->g_centre_.y;
    r->offset_.y = act_font_size / 2.0;
  }

  adjustStringRectsForSuperSubScript(draw_modes, rects);
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit
#endif  // __EMSCRIPTEN__
