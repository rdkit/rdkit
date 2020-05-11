//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Original author: David Cosgrove (CozChemIx) on 06/05/2020.
//

#include <cstdio>
#include <iostream>

#include <GraphMol/MolDraw2D/DrawTextFT.h>

using namespace std;

namespace RDKit {

// ****************************************************************************
DrawTextFT::DrawTextFT() :
    x_trans_(0), y_trans_(0), string_y_max_(0) {

  string fontfile = getFontFile();

  int err_code = FT_Init_FreeType(&library_);
  if(err_code != FT_Err_Ok) {
    throw runtime_error(string("Couldn't initialise Freetype."));
  }
  // take the first face
  err_code = FT_New_Face(library_, fontfile.c_str(), 0, &face_);
  if(err_code != FT_Err_Ok) {
    throw runtime_error(string("Font file ") + fontfile + string(" not found."));
  }

}

// ****************************************************************************
DrawTextFT::~DrawTextFT() {

  FT_Done_Face(face_);
  FT_Done_FreeType(library_);

}

// ****************************************************************************
double DrawTextFT::fontSize() const {

  return DrawText::fontSize() / POINT_SIZE;

}

// ****************************************************************************
void DrawTextFT::getStringSize(const std::string &label, double &label_width,
                               double &label_height) const {

  double x_min, y_min, x_max, y_max;
  GetStringBBox(label, x_min, y_min, x_max, y_max);

  label_width = x_max - x_min;
  label_height = y_max - y_min;

}

// ****************************************************************************
void DrawTextFT::drawString(const string &str, const Point2D &cds,
                            MolDraw2D::AlignType align) {

  cout << "DrawTextFT : " << str << " at " << cds.x << ", " << cds.y
       << " align = " << align << " fontScale = " << fontScale() << endl;

  Point2D char_cds;
  alignString(str, cds, align, char_cds);
  FT_Pos full_string_y_max = string_y_max_;
  MolDraw2D::TextDrawType draw_mode = MolDraw2D::TextDrawNormal;

  bool last_was_subs = false;
  double last_char_height = 0.0;
  FT_Pos last_x_trans = x_trans_;
  for (size_t i = 0, is = str.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == str[i] && setStringDrawMode(str, draw_mode, i)) {
      continue;
    }

    char next_c = str[i];
    cout << "Next_c : " << next_c << " : " << draw_mode
         << "last_char_height = " << last_char_height << endl;
    double x_min, y_min, x_max, y_max;
    string next_char( " ");
    next_char[0] = next_c;
    double oscale = fontScale();

    if(draw_mode == MolDraw2D::TextDrawSuperscript
       || draw_mode == MolDraw2D::TextDrawSubscript) {
      setFontScale(SUBS_SCALE * oscale);
    }
    GetStringBBox(next_char, x_min, y_min, x_max, y_max);
    string_y_max_ = full_string_y_max;
    if(draw_mode == MolDraw2D::TextDrawNormal) {
      last_x_trans = x_trans_;
      drawChar(next_c, char_cds);
      last_was_subs = false;
      last_char_height = y_max - y_min;
    } else {
      double char_height = y_max - y_min;
      cout << "XX " << next_c << " : fontScale = " << fontScale() << " char_height = "
           << char_height << endl;
      Point2D draw_cds;
      draw_cds.x = char_cds.x;
      if(draw_mode == MolDraw2D::TextDrawSubscript) {
        // if the last was also a super-/subscript put them one above
        // the other.  We really ought to step it back by width of
        // previous character, but we no longer have that.
        if(last_was_subs) {
          draw_cds.x = last_x_trans;
        }
        draw_cds.y = char_cds.y + 0.5 * last_char_height + 0.25 * char_height;
        last_was_subs = true;
      } else if(draw_mode == MolDraw2D::TextDrawSuperscript) {
        if(last_was_subs) {
          draw_cds.x = last_x_trans;
        }
        if(next_c == '-') {
          char_height *= 2;
        }
        draw_cds.y = char_cds.y - 0.5 * last_char_height + 0.5 * char_height;
        last_was_subs = true;
      }
      last_x_trans = x_trans_;
      drawChar(next_c, draw_cds);
    }
    setFontScale(oscale);
    // this somewhat inelegant thing is because drawChar updates x_trans_
    // by the glyph advance value which is where we want the next
    // character.
    char_cds.x = x_trans_;
  }

}

// ****************************************************************************
void DrawTextFT::drawChar(char c, const Point2D &cds) {

  cout << "draw " << c << " at " << cds.x << ", " << cds.y << endl;
  FT_Load_Char(face_, c, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP);
  x_trans_ = cds.x;
  y_trans_ = cds.y;
  double advance = ExtractOutline();
  x_trans_ += advance;

}

// ****************************************************************************
double DrawTextFT::FontCoordToPixelCoord(FT_Pos fc) const {

  static const double pixel_size = POINT_SIZE * RESOLUTION / 72.0;
  double pc = fontSize() * fc * pixel_size / (face_->units_per_EM);

  return pc;

}

// ****************************************************************************
void DrawTextFT::FontPosToDrawPos(FT_Pos fx, FT_Pos fy,
                                  double &dx, double &dy) const {

  dx = x_trans_ + FontCoordToPixelCoord(fx);
  // freetype has the origin at bottom left, we want it in the top left.
  FT_Pos rev_y = string_y_max_ - fy;
  dy = y_trans_ + FontCoordToPixelCoord(rev_y);

}

// ****************************************************************************
double DrawTextFT::ExtractOutline() {

  FT_Outline_Funcs callbacks;

  callbacks.move_to = MoveToFunction;
  callbacks.line_to = LineToFunction;
  callbacks.conic_to = ConicToFunction;
  callbacks.cubic_to = CubicToFunction;

  callbacks.shift = 0;
  callbacks.delta = 0;

  FT_GlyphSlot slot = face_->glyph;
  FT_Outline &outline = slot->outline;

  FT_Error error = FT_Outline_Decompose(&outline, &callbacks, this);
  if(error != FT_Err_Ok) {
    /* not sure what to do in this case */ ;
  }
  return FontCoordToPixelCoord(slot->advance.x);

}

// ****************************************************************************
string DrawTextFT::getFontFile() const {

  string ff_name = getenv("RDBASE") ? getenv("RDBASE") : "";
  if(ff_name.empty()) {
    return "NO_FONT_FILE"; // the c'tor will throw when this file isn't found.
  }
//  ff_name += "/Code/GraphMol/MolDraw2D/Roboto-Regular.ttf";
//  ff_name += "/Code/GraphMol/MolDraw2D/Roboto-Bold.ttf";
//  ff_name += "/Code/GraphMol/MolDraw2D/Roboto-Black.ttf";
//  ff_name += "/Code/GraphMol/MolDraw2D/Courier Prime Code.ttf";
//  ff_name += "/Code/GraphMol/MolDraw2D/Alatsi-Regular.ttf";
  ff_name += "/Code/GraphMol/MolDraw2D/Telex-Regular.ttf";
  return ff_name;

}

// ****************************************************************************
void DrawTextFT::GetStringBBox(const std::string &text,
                               double &x_min, double &y_min ,
                               double &x_max, double &y_max) const {

  cout << "GetStringBBox for " << text << endl;

  FT_Pos this_x_min, this_y_min, this_x_max, this_y_max, advance;
  x_min = x_max = y_max = 0;
  y_min = numeric_limits<double>::max();
  string_y_max_ = numeric_limits<FT_Pos>::min();

  MolDraw2D::TextDrawType draw_mode = MolDraw2D::TextDrawNormal;

  bool last_was_subs = false;
  for(size_t i = 0; i < text.length(); ++i) {

    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    CalcGlyphBBox(text[i], this_x_min, this_y_min, this_x_max, this_y_max,
                  advance);
    cout << "next char : " << text[i] << " last_was_subs = " << last_was_subs << endl;
    // advance is the full width of the glyph, including the space either side.
    double oscale = 1.0;
    if(draw_mode == MolDraw2D::TextDrawSubscript
       || draw_mode == MolDraw2D::TextDrawSuperscript) {
      oscale = SUBS_SCALE;
    }
    if(!i) {
      x_min = oscale * FontCoordToPixelCoord(this_x_min);
    }
    if(!last_was_subs || draw_mode == MolDraw2D::TextDrawNormal) {
      x_max += oscale * FontCoordToPixelCoord(advance);
    }
    y_min = min(y_min, oscale * FontCoordToPixelCoord(this_y_min));
    y_max = max(y_max, oscale * FontCoordToPixelCoord(this_y_max));
    string_y_max_ = max(string_y_max_, this_y_max);
    if(draw_mode == MolDraw2D::TextDrawSubscript
       || draw_mode == MolDraw2D::TextDrawSuperscript) {
      last_was_subs = true;
    } else {
      last_was_subs = false;
    }
    cout << "x_max = " << x_max << endl;
  }

  cout << text << " : fs = " << fontScale() << " : " << x_min << ", " << y_min
       << " to " << x_max << ", " << y_max << endl;

}

// ****************************************************************************
void DrawTextFT::CalcGlyphBBox(char c, FT_Pos &x_min, FT_Pos &y_min,
                               FT_Pos &x_max , FT_Pos &y_max,
                               FT_Pos &advance) const {

  FT_Load_Char(face_, c, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP);
  FT_GlyphSlot slot = face_->glyph;
  FT_Outline &outline = slot->outline;
  FT_BBox bbox;

  FT_Outline_Get_BBox(&outline, &bbox);

  x_min = bbox.xMin;
  y_min = bbox.yMin;
  x_max = bbox.xMax;
  y_max = bbox.yMax;
  advance = slot->advance.x;

  cout << "bbox for " << c << " : " << x_min << ", " << y_min
       << " and " << x_max << ", " << y_max
       << " advance = " << advance << endl;

}

// ****************************************************************************
int MoveToFunction(const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->MoveToFunctionImpl(to);

}

// ****************************************************************************
int LineToFunction(const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->LineToFunctionImpl(to);

}

// ****************************************************************************
int ConicToFunction(const FT_Vector *control, const FT_Vector *to,
                    void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->ConicToFunctionImpl(control, to);

}

// ****************************************************************************
int CubicToFunction(const FT_Vector *controlOne,
                    const FT_Vector *controlTwo,
                    const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->CubicToFunctionImpl(controlOne, controlTwo, to);

}

} // namespace RDKit