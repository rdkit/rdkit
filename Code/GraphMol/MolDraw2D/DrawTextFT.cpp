//
//  Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Original author: David Cosgrove (CozChemIx).
//

#include <cstdio>
#include <iostream>

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/DrawTextFT.h>
#include <GraphMol/MolDraw2D/Fonts/telex_regular.h>
#include <GraphMol/MolDraw2D/Fonts/roboto_regular.h>

namespace RDKit {
extern const std::string telex_ttf;

namespace MolDraw2D_detail {

// ****************************************************************************
DrawTextFT::DrawTextFT(double max_fnt_sz, double min_fnt_sz,
                       const std::string &font_file)
    : DrawText(max_fnt_sz, min_fnt_sz),
      library_(nullptr),
      face_(nullptr),
      x_trans_(0),
      y_trans_(0),
      string_y_max_(0) {
  int err_code = FT_Init_FreeType(&library_);
  if (err_code != FT_Err_Ok) {
    throw std::runtime_error(std::string("Couldn't initialise Freetype."));
  }
  setFontFile(font_file);
}

// ****************************************************************************
DrawTextFT::~DrawTextFT() {
  FT_Done_Face(face_);
  FT_Done_FreeType(library_);
}

// ****************************************************************************
void DrawTextFT::drawChar(char c, const Point2D &cds) {
  FT_Load_Char(face_, c, FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP);
  x_trans_ = cds.x;
  y_trans_ = cds.y;
  extractOutline();
}

// ****************************************************************************
double DrawTextFT::fontCoordToDrawCoord(FT_Pos fc) const {
  double pc = fontSize() * fc * em_scale_;

  return pc;
}

// ****************************************************************************
void DrawTextFT::fontPosToDrawPos(FT_Pos fx, FT_Pos fy, double &dx,
                                  double &dy) const {
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
  if (error != FT_Err_Ok) {
    /* not sure what to do in this case */;
  }
  return fontCoordToDrawCoord(slot->advance.x);
}

// ****************************************************************************
// This is currently not being used
std::string DrawTextFT::getFontFile() const {
  if (!font_file_.empty()) {
    return font_file_;
  }
  std::string ff_name = getenv("RDBASE") ? getenv("RDBASE") : "";
  if (ff_name.empty()) {
    throw std::runtime_error(
        "Freetype not initialised because RDBASE not defined.");
  }
  ff_name += "/Data/Fonts/Telex-Regular.ttf";
  return ff_name;
}

// ****************************************************************************
void DrawTextFT::setFontFile(const std::string &font_file) {
  if (face_ && font_file == font_file_) {
    return;
  }

  font_file_ = font_file;
  FT_Done_Face(face_);
  face_ = nullptr;
  // take the first face
  const std::string *font_string = nullptr;
  if (!font_file_.empty()) {
    if (font_file == "BuiltinTelexRegular") {
      font_string = &telex_regular_ttf;
    } else if (font_file == "BuiltinRobotoRegular") {
      font_string = &roboto_regular_ttf;
    } else {
      int err_code = FT_New_Face(library_, font_file_.c_str(), 0, &face_);
      if (err_code != FT_Err_Ok) {
        throw std::runtime_error(std::string("Font file ") + font_file_ +
                                 std::string(" not found."));
      }
    }
  } else {
    font_string = &telex_regular_ttf;
  }
  if (font_string) {
    int err_code = FT_New_Memory_Face(library_, (FT_Byte *)font_string->c_str(),
                                      font_string->size(), 0, &face_);
    if (err_code != FT_Err_Ok) {
      throw std::runtime_error("could not load embedded font data");
    }
  }
  em_scale_ = 1.0 / face_->units_per_EM;
}

// ****************************************************************************
void DrawTextFT::getStringRects(const std::string &text,
                                std::vector<std::shared_ptr<StringRect>> &rects,
                                std::vector<TextDrawType> &draw_modes,
                                std::vector<char> &draw_chars) const {
  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  double max_y = 0.0;
  double mean_width = 0.0;
  std::vector<double> extras;
  for (size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    draw_chars.push_back(text[i]);
    FT_Pos this_x_min, this_y_min, this_x_max, this_y_max, advance;
    calcGlyphBBox(text[i], this_x_min, this_y_min, this_x_max, this_y_max,
                  advance);
    double oscale = selectScaleFactor(text[i], draw_mode);
    double p_x_min = oscale * fontCoordToDrawCoord(this_x_min);
    double p_y_min = oscale * fontCoordToDrawCoord(this_y_min);
    double p_x_max = oscale * fontCoordToDrawCoord(this_x_max);
    double p_y_max = oscale * fontCoordToDrawCoord(this_y_max);
    double p_advance = oscale * fontCoordToDrawCoord(advance);
    double width = p_x_max - p_x_min;
    // The mean width is to define the spacing between the characters, so
    // use full size characters.
    mean_width += width / oscale;
    // reduce the horizontal offset, which is the distance
    // of the start of the glyph from the start of the char box.
    // Otherwise spacing is uneven.
    extras.push_back(p_advance - p_x_max);
    if (!this_x_max) {
      // it was a space, probably, and we want small spaces because screen
      // real estate is limited.
      width = p_advance / 3;
    }
    double height = p_y_max - p_y_min;
    Point2D offset(p_x_min + width / 2.0, p_y_max / 2.0);
    Point2D g_centre(offset.x, p_y_max - height / 2.0);
    rects.push_back(std::shared_ptr<StringRect>(
        new StringRect(offset, g_centre, width, height)));
    draw_modes.push_back(draw_mode);
    max_y = std::max(max_y, p_y_max);
  }
  // Use the mean width of the characters to define some extra space between
  // the characters.
  mean_width /= rects.size();
  for (auto i = 0u; i < rects.size(); ++i) {
    extras[i] += mean_width / 10;
    rects[i]->g_centre_.y = max_y - rects[i]->g_centre_.y;
    rects[i]->offset_.y = max_y / 2.0;
    if (i) {
      rects[i]->trans_.x = rects[i - 1]->trans_.x + rects[i - 1]->width_ / 2 +
                           rects[i]->width_ / 2 + extras[i - 1];
    }
  }

  adjustStringRectsForSuperSubScript(draw_modes, rects);
}

// ****************************************************************************
void DrawTextFT::calcGlyphBBox(char c, FT_Pos &x_min, FT_Pos &y_min,
                               FT_Pos &x_max, FT_Pos &y_max,
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
int conicToFunction(const FT_Vector *control, const FT_Vector *to, void *user) {
  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->ConicToFunctionImpl(control, to);
}

// ****************************************************************************
int cubicToFunction(const FT_Vector *controlOne, const FT_Vector *controlTwo,
                    const FT_Vector *to, void *user) {
  auto *rdft = static_cast<DrawTextFT *>(user);
  return rdft->CubicToFunctionImpl(controlOne, controlTwo, to);
}
}  // namespace MolDraw2D_detail

}  // namespace RDKit
