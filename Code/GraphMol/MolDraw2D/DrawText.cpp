//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx) on 29/04/2020.
//


#include <GraphMol/MolDraw2D/DrawText.h>

using namespace std;

namespace RDKit {

// ****************************************************************************
DrawText::DrawText() : font_scale_(1.0) {}

// ****************************************************************************
DrawColour const &DrawText::colour() const {
  return colour_;
}

// ****************************************************************************
void DrawText::setColour(const DrawColour &col) {
  colour_ = col;
}

// ****************************************************************************
double DrawText::fontSize() const {
  return font_scale_ * FONT_SIZE;
}

// ****************************************************************************
double DrawText::fontScale() const {
  return font_scale_;
}

// ****************************************************************************
void DrawText::setFontScale(double new_scale) { font_scale_ = new_scale; }

// ****************************************************************************
void DrawText::drawString(const string &str, const Point2D &cds,
                          TextAlignType align) {
  RDUNUSED_PARAM(align);
  cout << "Drawing string : " << str << " at " << cds.x << ", " << cds.y << endl;
  cout << "font size  " << fontSize() << " font scale : " << fontScale() << endl;

  double string_width, string_height;
  getStringSize(str, string_width, string_height);

  // FIX: this shouldn't stay
  double M_width, M_height;
  getStringSize(std::string("M"), M_width, M_height);

  double draw_x = cds.x - string_width / 2.0;
  double draw_y = cds.y + string_height / 2.0;

  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  string next_char(" ");

  for (size_t i = 0, is = str.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == str[i] && setStringDrawMode(str, draw_mode, i)) {
      continue;
    }

    char next_c = str[i];
    next_char[0] = next_c;
    double char_width, char_height;
    getStringSize(next_char, char_width, char_height);

    // these font sizes and positions work best for Qt, IMO. They may want
    // tweaking for a more general solution.
    if (TextDrawType::TextDrawSubscript == draw_mode) {
      // y goes from top to bottom, so add for a subscript!
      double oscale = fontScale();
      setFontScale(0.75 * oscale);
      char_width *= 0.5;
      drawChar(next_c,
               Point2D(draw_x, draw_y + 0.5 * char_height));
      setFontScale(oscale);
    } else if (TextDrawType::TextDrawSuperscript == draw_mode) {
      double oscale = fontScale();
      setFontScale(0.75 * oscale);
      char_width *= 0.5;
      drawChar(next_c,
               Point2D(draw_x, draw_y - .5 * M_height));
      setFontScale(oscale);
    } else {
      drawChar(next_c, Point2D(draw_x, draw_y));
    }
    draw_x += char_width;
  }
}

// ****************************************************************************
void DrawText::alignString(const std::string &str, const Point2D &in_cds,
                           TextAlignType align, Point2D &out_cds) {

  double width, height;
  getStringSize(str, width, height);

  // freetype has the origin at bottom left, we want it in the top left.
  // centring at the same time.
  // TODO - this is probably wrong sign for height / 2 for non-freetype.
  out_cds.y = in_cds.y + height / 2;

  double c_width, c_height;
  getStringSize(str.substr(0, 1), c_width, c_height);

  switch(align) {
    case TextAlignType::START:
      cout << "align = start" << endl;
      out_cds.x = in_cds.x - c_width / 2;
      break;
    case TextAlignType::MIDDLE:
      cout << "align = middle" << endl;
      out_cds.x = in_cds.x - width / 2;
      break;
    case TextAlignType::END:
      cout << "align = end" << endl;
      out_cds.x = in_cds.x - width + c_width / 2;
      break;
  }

}

// ****************************************************************************
bool setStringDrawMode(const string &instring,
                       TextDrawType &draw_mode,
                       size_t &i) {

  string bit1 = instring.substr(i, 5);
  string bit2 = instring.substr(i, 6);

  // could be markup for super- or sub-script
  if (string("<sub>") == bit1) {
    draw_mode = TextDrawType::TextDrawSubscript;
    i += 4;
    return true;
  } else if (string("<sup>") == bit1) {
    draw_mode = TextDrawType::TextDrawSuperscript;
    i += 4;
    return true;
  } else if (string("</sub>") == bit2) {
    draw_mode = TextDrawType::TextDrawNormal;
    i += 5;
    return true;
  } else if (string("</sup>") == bit2) {
    draw_mode = TextDrawType::TextDrawNormal;
    i += 5;
    return true;
  }

  return false;
}

std::ostream& operator<<(std::ostream &oss, TextAlignType tat) {
  switch(tat) {
    case TextAlignType::START: oss << "START"; break;
    case TextAlignType::MIDDLE: oss << "MIDDLE"; break;
    case TextAlignType::END: oss << "END"; break;
  }
  return oss;
}
std::ostream& operator<<(std::ostream &oss, TextDrawType tdt) {
  switch(tdt) {
    case TextDrawType::TextDrawNormal: oss << "TextDrawNormal"; break;
    case TextDrawType::TextDrawSuperscript: oss << "TextDrawSuperscript"; break;
    case TextDrawType::TextDrawSubscript: oss << "TextDrawSubscript"; break;
  }
  return oss;
}

} // namespace RDKit