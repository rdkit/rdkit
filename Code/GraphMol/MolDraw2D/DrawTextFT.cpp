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

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/DrawTextFT.h>

using namespace std;

namespace RDKit {

// ****************************************************************************
DrawTextFT::DrawTextFT() :
    x_trans_(0), y_trans_(0), string_y_max_(0) {

  cout << "DrawTextFT : " << FONT_SIZE << endl;

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

  return fontScale() * FONT_SIZE / POINT_SIZE;

}

// ****************************************************************************
void DrawTextFT::drawChar(char c, const Point2D &cds) {

//  cout << "draw " << c << " at " << cds.x << ", " << cds.y << endl;
  FT_Load_Char(face_, c, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP);
  x_trans_ = cds.x;
  y_trans_ = cds.y;
  double advance = extractOutline();
  x_trans_ += advance;

}

// ****************************************************************************
double DrawTextFT::fontCoordToPixelCoord(FT_Pos fc) const {

  static const double pixel_size = POINT_SIZE * RESOLUTION / 72.0;
  double pc = fontSize() * fc * pixel_size / (face_->units_per_EM);

  return pc;

}

// ****************************************************************************
void DrawTextFT::fontPosToDrawPos(FT_Pos fx, FT_Pos fy,
                                  double &dx, double &dy) const {

  dx = x_trans_ + fontCoordToPixelCoord(fx);
  // freetype has the origin at bottom left, we want it in the top left.
  FT_Pos rev_y = string_y_max_ - fy;
  dy = y_trans_ + fontCoordToPixelCoord(rev_y);

}

// ****************************************************************************
double DrawTextFT::extractOutline() {

  FT_Outline_Funcs callbacks;

  callbacks.move_to = moveToFunction;
  callbacks.line_to = lineToFunction;
  callbacks.conic_to = conicToFunction;
  callbacks.cubic_to = cubicToFunction;

  callbacks.shift = 0;
  callbacks.delta = 0;

  FT_GlyphSlot slot = face_->glyph;
  FT_Outline &outline = slot->outline;

  FT_Error error = FT_Outline_Decompose(&outline, &callbacks, this);
  if(error != FT_Err_Ok) {
    /* not sure what to do in this case */ ;
  }
  return fontCoordToPixelCoord(slot->advance.x);

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
void DrawTextFT::getStringRects(const string &text,
                                vector<shared_ptr<StringRect>> &rects,
                                vector<TextDrawType> &draw_modes) const {

  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  double running_x = 0.0;
  vector<char> to_draw;
  for(size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    to_draw.push_back(text[i]);
//    cout << "Rects : " << text[i] << " : " << running_x << endl;
    FT_Pos this_x_min, this_y_min, this_x_max, this_y_max, advance;
    calcGlyphBBox(text[i], this_x_min, this_y_min, this_x_max, this_y_max,
                  advance);
    double oscale = 1.0;
    if(draw_mode == TextDrawType::TextDrawSubscript) {
      oscale = SUBS_SCALE;
    } else if(draw_mode == TextDrawType::TextDrawSuperscript) {
      if(text[i] == '-' || text[i] == '+') {
        oscale = SUBS_SCALE;
      } else {
        oscale = SUPER_SCALE;
      }
    }
    double p_x_min = oscale * fontCoordToPixelCoord(this_x_min);
    double p_y_min = oscale * fontCoordToPixelCoord(this_y_min);
    double width = oscale * fontCoordToPixelCoord(this_x_max) - p_x_min;
    double height = oscale * fontCoordToPixelCoord(this_y_max) - p_y_min;
    Point2D centre(running_x - p_x_min + 0.5 * width, 0.5 * height);
    rects.push_back(shared_ptr<StringRect>(new StringRect(centre, width, height)));
    draw_modes.push_back(draw_mode);
    // TODO: include kerning
    running_x += oscale * (fontCoordToPixelCoord(advance) - p_x_min);
  }

  adjustStringRectsForSuperSubScript(draw_modes, to_draw, rects);

//  for(size_t i = 0; i < rects.size(); ++i) {
//    cout << i << " : " << rects[i]->centre_ << " and " << rects[i]->width_
//         << " by " << rects[i]->height_ << " type " << draw_modes[i] << endl;
//  }
}

// ****************************************************************************
void DrawTextFT::calcGlyphBBox(char c, FT_Pos &x_min, FT_Pos &y_min,
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

//  cout << "bbox for " << c << " : " << x_min << ", " << y_min
//       << " and " << x_max << ", " << y_max
//       << " advance = " << advance << endl;

}

// ****************************************************************************
int moveToFunction(const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->MoveToFunctionImpl(to);

}

// ****************************************************************************
int lineToFunction(const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->LineToFunctionImpl(to);

}

// ****************************************************************************
int conicToFunction(const FT_Vector *control, const FT_Vector *to,
                    void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->ConicToFunctionImpl(control, to);

}

// ****************************************************************************
int cubicToFunction(const FT_Vector *controlOne,
                    const FT_Vector *controlTwo,
                    const FT_Vector *to, void *user) {

  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->CubicToFunctionImpl(controlOne, controlTwo, to);

}

} // namespace RDKit