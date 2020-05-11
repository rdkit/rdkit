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
                          MolDraw2D::AlignType align) {
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

  MolDraw2D::TextDrawType draw_mode = MolDraw2D::TextDrawNormal;
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
    if (MolDraw2D::TextDrawSubscript == draw_mode) {
      // y goes from top to bottom, so add for a subscript!
      double oscale = fontScale();
      setFontScale(0.75 * oscale);
      char_width *= 0.5;
      drawChar(next_c,
               Point2D(draw_x, draw_y + 0.5 * char_height));
      setFontScale(oscale);
    } else if (MolDraw2D::TextDrawSuperscript == draw_mode) {
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
                           MolDraw2D::AlignType align, Point2D &out_cds) {

  double width, height;
  getStringSize(str, width, height);

  // freetype has the origin at bottom left, we want it in the top left.
  // centring at the same time.
  out_cds.y = in_cds.y - height / 2;

  double c_width, c_height;
  getStringSize(str.substr(0, 1), c_width, c_height);

  switch(align) {
    case MolDraw2D::START:
      cout << "align = start" << endl;
      out_cds.x = in_cds.x - c_width / 2;
      break;
    case MolDraw2D::MIDDLE:
      cout << "align = middle" << endl;
      out_cds.x = in_cds.x - width / 2;
      break;
    case MolDraw2D::END:
      cout << "align = end" << endl;
      out_cds.x = in_cds.x - width + c_width / 2;
      break;
  }

}

// ****************************************************************************
bool setStringDrawMode(const string &instring,
                       MolDraw2D::TextDrawType &draw_mode,
                       size_t &i) {

  string bit1 = instring.substr(i, 5);
  string bit2 = instring.substr(i, 6);

  // could be markup for super- or sub-script
  if (string("<sub>") == bit1) {
    draw_mode = MolDraw2D::TextDrawSubscript;
    i += 4;
    return true;
  } else if (string("<sup>") == bit1) {
    draw_mode = MolDraw2D::TextDrawSuperscript;
    i += 4;
    return true;
  } else if (string("</sub>") == bit2) {
    draw_mode = MolDraw2D::TextDrawNormal;
    i += 5;
    return true;
  } else if (string("</sup>") == bit2) {
    draw_mode = MolDraw2D::TextDrawNormal;
    i += 5;
    return true;
  }

  return false;
}

} // namespace RDKit