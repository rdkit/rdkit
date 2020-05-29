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
DrawTextFT::DrawTextFT(double max_fnt_sz) :
    DrawText(max_fnt_sz), x_trans_(0), y_trans_(0), string_y_max_(0) {

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
  cout << "max pixel thing : " << fontCoordToDrawCoord(face_->units_per_EM) << endl;
  cout << "units per EM = " << face_->units_per_EM << endl;
  em_scale_ = 1.0 / face_->units_per_EM;

}

// ****************************************************************************
DrawTextFT::~DrawTextFT() {

  FT_Done_Face(face_);
  FT_Done_FreeType(library_);

}

// ****************************************************************************
void DrawTextFT::drawChar(char c, const Point2D &cds) {

  cout << "draw " << c << " at " << cds.x << ", " << cds.y << endl;
  FT_Load_Char(face_, c, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP);
  auto min_it = char_mins_.find(c);
  FT_Pos x_min = 0, y_min = 0;
  if(min_it != char_mins_.end()) {
    x_min = min_it->second.first;
    y_min = min_it->second.second;
  }
  x_trans_ = cds.x - fontCoordToDrawCoord(x_min);
  y_trans_ = cds.y + fontCoordToDrawCoord(y_min);
  double advance = extractOutline();
  x_trans_ += advance;

}

// ****************************************************************************
double DrawTextFT::fontCoordToDrawCoord(FT_Pos fc) const {

  double pc = fontSize() * fc * em_scale_;

  return pc;

}

// ****************************************************************************
void DrawTextFT::fontPosToDrawPos(FT_Pos fx, FT_Pos fy,
                                  double &dx, double &dy) const {

  dx = x_trans_ + fontCoordToDrawCoord(fx);
  // freetype has the origin at bottom left, we want it in the top left.
  FT_Pos rev_y = string_y_max_ - fy;
  dy = y_trans_ + fontCoordToDrawCoord(rev_y);

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
  return fontCoordToDrawCoord(slot->advance.x);

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
                                vector<TextDrawType> &draw_modes,
                                vector<char> &draw_chars) const {

  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  double running_x = 0.0;
  for(size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    draw_chars.push_back(text[i]);
//    cout << "Rects : " << text[i] << " : " << running_x << endl;
    FT_Pos this_x_min, this_y_min, this_x_max, this_y_max, advance;
    calcGlyphBBox(text[i], this_x_min, this_y_min, this_x_max, this_y_max,
                  advance);
    double oscale = selectScaleFactor(text[i], draw_mode);
    double p_x_min = oscale * fontCoordToDrawCoord(this_x_min);
    double p_y_min = oscale * fontCoordToDrawCoord(this_y_min);
    double p_x_max = oscale * fontCoordToDrawCoord(this_x_max);
    double p_y_max = oscale * fontCoordToDrawCoord(this_y_max);
    cout << "p_y_min to p_y_max : " << p_y_min << " -> " << p_y_max << endl;
    double width = rect_scale_ * (p_x_max - p_x_min);
    double height = rect_scale_ * (p_y_max - p_y_min);
    cout << "height = " << height << endl;
    if(!i) {
      running_x = -p_x_min;
    }
    Point2D centre(running_x + 0.5 * (p_x_max + p_x_min),
                   0.5 * p_y_max);
    rects.push_back(shared_ptr<StringRect>(new StringRect(centre, width, height)));
    draw_modes.push_back(draw_mode);
    running_x += oscale * fontCoordToDrawCoord(advance);
    char_mins_[text[i]] = make_pair(this_x_min, this_y_min);
  }

  adjustStringRectsForSuperSubScript(draw_modes, draw_chars, rects);

  cout << "DrawTextFT::getStringRects : " << endl;
  for(size_t i = 0; i < rects.size(); ++i) {
    cout << i << " : " << rects[i]->trans_ << " and " << rects[i]->width_
         << " by " << rects[i]->height_ << " type " << draw_modes[i] << endl;
  }
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

  cout << "bbox for " << c << " : " << x_min << ", " << y_min
       << " and " << x_max << ", " << y_max
       << " advance = " << advance << endl;

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