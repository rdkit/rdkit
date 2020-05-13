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
  getStringBBox(label, x_min, y_min, x_max, y_max);

  label_width = x_max - x_min;
  label_height = y_max - y_min;

}

// ****************************************************************************
void DrawTextFT::drawString(const string &str, const Point2D &cds,
                            TextAlignType align) {

  cout << "XXXXXXXXXXXXXXXXX" << endl;
  cout << "DrawTextFT : " << str << " at " << cds.x << ", " << cds.y
       << " align = " << align << " fontScale = " << fontScale() << endl;

  Point2D a_cds;
  alignString(str, cds, align, a_cds);

  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  getStringRects(str, rects, draw_modes);

  double full_scale = fontScale();
  double subs_scale = SUBS_SCALE * fontScale();
  double super_scale = SUPER_SCALE * fontScale();
  for(size_t i = 0, j = 0; i < str.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == str[i] && setStringDrawMode(str, draw_modes[j], i)) {
      continue;
    }

    switch(draw_modes[j]) {
      case TextDrawType::TextDrawNormal:
        setFontScale(full_scale); break;
      case TextDrawType::TextDrawSubscript:
        setFontScale(subs_scale); break;
      case TextDrawType::TextDrawSuperscript:
        if(str[i] == '-' || str[i] == '+') {
          setFontScale(subs_scale);
        } else {
          setFontScale(super_scale);
        }
        break;
    }

    Point2D draw_cds;
    draw_cds.x = a_cds.x + rects[j]->centre_.x - 0.5 * rects[j]->width_;
    draw_cds.y = a_cds.y + rects[j]->centre_.y - 0.5 * rects[j]->height_;
    cout << "Drawing " << str[i] << " at " << draw_cds << endl;
    drawChar(str[i], draw_cds);
    ++j;
  }
  setFontScale(full_scale);

}

// ****************************************************************************
void DrawTextFT::drawChar(char c, const Point2D &cds) {

  cout << "draw " << c << " at " << cds.x << ", " << cds.y << endl;
  FT_Load_Char(face_, c, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP);
  x_trans_ = cds.x;
  y_trans_ = cds.y;
  double advance = extractOutline();
  x_trans_ += advance;

}

// ****************************************************************************
void DrawTextFT::alignString(const std::string &str, const Point2D &in_cds,
                             TextAlignType align, Point2D &out_cds) {

  cout << "Align string : " << str << endl;
  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  getStringRects(str, rects, draw_modes);

  double full_width = 0.0;
  for(const auto r: rects) {
    full_width += r->width_;
  }

  double align_h = 0.0;
  double align_w = 0.0;
  if(align == TextAlignType::START || align == TextAlignType::END) {
    size_t align_char = 0;
    for (size_t i = 0; i < rects.size(); ++i) {
      if (draw_modes[i] == TextDrawType::TextDrawNormal) {
        align_char = i;
        if (align == TextAlignType::START) {
          break;
        }
      }
    }
    cout << "Align char : " << align_char << endl;
    align_h = rects[align_char]->height_;
    align_w = rects[align_char]->width_;
  } else {
    double y_min = numeric_limits<double>::max();
    double y_max = -y_min;
    for(const auto r: rects) {
      y_max = max(y_max, r->centre_.y + r->height_ / 2);
      y_min = min(y_min, r->centre_.y - r->height_ / 2);
    }
    align_h = y_max - y_min;
    align_w = full_width;
  }
  cout << "Align height : " << align_h << endl;
  out_cds.y = in_cds.y + align_h / 2;

  switch(align) {
    case TextAlignType::START:
      cout << "align = start" << endl;
      out_cds.x = in_cds.x - align_w / 2;
      break;
    case TextAlignType::MIDDLE:
      cout << "align = middle" << endl;
      out_cds.x = in_cds.x - align_w / 2;
      break;
    case TextAlignType::END:
      cout << "align = end" << endl;
      out_cds.x = in_cds.x - full_width + align_w / 2;
      break;
  }

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
void DrawTextFT::getStringBBox(const std::string &text,
                               double &x_min, double &y_min ,
                               double &x_max, double &y_max) const {

  cout << "getStringBBox for " << text << endl;
  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes; // not needed for bounding box
  getStringRects(text, rects, draw_modes);
  x_min = y_min = numeric_limits<double>::max();
  x_max = y_max = 0.0;
  for(const auto rect: rects) {
    x_min = min(x_min, rect->centre_.x - rect->width_ / 2.0);
    y_min = min(y_min, rect->centre_.y - rect->height_ / 2.0);
    x_max = max(x_max, rect->centre_.x + rect->width_ / 2.0);
    y_max = max(y_max, rect->centre_.y + rect->height_ / 2.0);
  }

  cout << text << " : fs = " << fontScale() << " : " << x_min << ", " << y_min
       << " to " << x_max << ", " << y_max << endl;

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
    cout << "Rects : " << text[i] << " : " << running_x << endl;
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

  // adjust super- and subscripts
  double last_height = -1.0;
  for(size_t i = 0; i < draw_modes.size(); ++i) {
    switch(draw_modes[i]) {
      case TextDrawType::TextDrawSuperscript:
        // superscripts may come before or after a letter.  If before,
        // spin through to first non-superscript for last_height
        if(last_height < 0.0) {
          for (size_t j = i + 1; j < draw_modes.size(); ++j) {
            if (draw_modes[j] == TextDrawType::TextDrawNormal) {
              last_height = rects[j]->height_;
              break;
            }
          }
        }
        // adjust y up by last_height / 2, unless it's a + which tends
        // to stick above.
        if(to_draw[i] == '+') {
          rects[i]->centre_.y -= 0.33 * last_height;
        } else {
          rects[i]->centre_.y -= 0.5 * last_height;
        }
        // if the last char was a subscript, remove the advance
        if (i && draw_modes[i - 1] == TextDrawType::TextDrawSubscript) {
          rects[i]->centre_.x = rects[i - 1]->centre_.x;
        }
        break;
      case TextDrawType::TextDrawSubscript:
        // adjust y down by last_height / 2
        rects[i]->centre_.y += 0.5 * last_height;
        // if the last char was a superscript, remove the advance
        if (i && draw_modes[i - 1] == TextDrawType::TextDrawSuperscript) {
          rects[i]->centre_.x = rects[i - 1]->centre_.x;
        }
        break;
      case TextDrawType::TextDrawNormal:
        last_height = rects[i]->height_;
        break;
    }
  }

  for(size_t i = 0; i < rects.size(); ++i) {
    cout << i << " : " << rects[i]->centre_ << " and " << rects[i]->width_
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