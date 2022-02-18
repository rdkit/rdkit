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
// A concrete class derived from DrawText that uses SVG
// to draw text onto a picture.

#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/DrawTextSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>

namespace RDKit {

std::string DrawColourToSVG(const DrawColour &col);

namespace MolDraw2D_detail {
// ****************************************************************************
DrawTextSVG::DrawTextSVG(double max_fnt_sz, double min_fnt_sz,
                         std::ostream &oss, std::string &d_act_class)
    : DrawTextNotFT(max_fnt_sz, min_fnt_sz),
      oss_(oss),
      d_active_class_(d_act_class) {}

namespace {
void escape_xhtml(std::string &data) {
  boost::algorithm::replace_all(data, "&", "&amp;");
  boost::algorithm::replace_all(data, "\"", "&quot;");
  boost::algorithm::replace_all(data, "\'", "&apos;");
  boost::algorithm::replace_all(data, "<", "&lt;");
  boost::algorithm::replace_all(data, ">", "&gt;");
}
}  // namespace

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void DrawTextSVG::drawChar(char c, const Point2D &cds) {
  unsigned int fontSz = fontSize();
  std::string col = DrawColourToSVG(colour());

  oss_ << "<text";
  oss_ << " x='" << MolDraw2D_detail::formatDouble(cds.x);
  oss_ << "' y='" << MolDraw2D_detail::formatDouble(cds.y) << "'";
  if (!d_active_class_.empty()) {
    oss_ << " class='" << d_active_class_ << "'";
  }
  std::string cs;
  cs += c;
  escape_xhtml(cs);
  oss_ << " style='font-size:" << fontSz
       << "px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;"
          "font-family:sans-serif;text-anchor:start;"
       << "fill:" << col << "'";
  oss_ << " >";
  oss_ << cs;
  oss_ << "</text>" << std::endl;
}

// ****************************************************************************
void DrawTextSVG::getStringRects(
    const std::string &text, std::vector<std::shared_ptr<StringRect>> &rects,
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
