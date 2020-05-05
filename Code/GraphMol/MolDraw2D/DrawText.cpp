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
DrawText::DrawText() : font_size_(16) {}

// ****************************************************************************
DrawColour DrawText::colour() const {
  return colour_;
}

// ****************************************************************************
void DrawText::setColour(const DrawColour &col) {
  colour_ = col;
}

// ****************************************************************************
double DrawText::fontSize() const {
  return font_size_;
}

// ****************************************************************************
void DrawText::setFontSize(double new_size) { font_size_ = new_size; }

// ****************************************************************************
void DrawText::drawString(const string &str, const Point2D &cds,
                          MolDraw2D::AlignType align) {
  RDUNUSED_PARAM(align);
  cout << "Drawing string : " << str << " at " << cds.x << ", " << cds.y << endl;
  cout << "font size  " << fontSize() << endl;

  double string_width, string_height;
  getStringSize(str, string_width, string_height);

  // FIX: this shouldn't stay
  double M_width, M_height;
  getStringSize(std::string("M"), M_width, M_height);

  double draw_x = cds.x - string_width / 2.0;
  double draw_y = cds.y + string_height / 2.0;

  double full_font_size = fontSize();
  MolDraw2D::TextDrawType draw_mode = MolDraw2D::TextDrawNormal;
  string next_char(" ");

  for (int i = 0, is = str.length(); i < is; ++i) {
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
      setFontSize(0.75 * full_font_size);
      char_width *= 0.5;
      drawChar(next_c,
               Point2D(draw_x, draw_y + 0.5 * char_height));
      setFontSize(full_font_size);
    } else if (MolDraw2D::TextDrawSuperscript == draw_mode) {
      setFontSize(0.75 * full_font_size);
      char_width *= 0.5;
      drawChar(next_c,
               Point2D(draw_x, draw_y - .5 * M_height));
      setFontSize(full_font_size);
    } else {
      setFontSize(full_font_size);
      drawChar(next_c, Point2D(draw_x, draw_y));
    }
    draw_x += char_width;
  }
}

// ****************************************************************************
bool setStringDrawMode(const string &instring,
                       MolDraw2D::TextDrawType &draw_mode,
                       int &i) {

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