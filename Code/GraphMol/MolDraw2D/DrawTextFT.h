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
// This is an abstract base class derived from DrawText that does drawing
// using FreeType.

#ifndef RDKIT_DRAWTEXTFT_H
#define RDKIT_DRAWTEXTFT_H

#include <string>

#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_BBOX_H
#include FT_OUTLINE_H

#include <GraphMol/MolDraw2D/DrawText.h>

namespace RDKit {

struct StringRect;

// ****************************************************************************
class DrawTextFT : public DrawText {
 public:
  DrawTextFT(double max_fnt_sz, double min_fnt_sz,
             const std::string &font_file);
  ~DrawTextFT();

  void drawChar(char c, const Point2D &cds) override;

  virtual int MoveToFunctionImpl(const FT_Vector *to) = 0;
  virtual int LineToFunctionImpl(const FT_Vector *to) = 0;
  virtual int ConicToFunctionImpl(const FT_Vector *control,
                                  const FT_Vector *to) = 0;
  virtual int CubicToFunctionImpl(const FT_Vector *controlOne,
                                  const FT_Vector *controlTwo,
                                  const FT_Vector *to) = 0;

  // unless over-ridden by the c'tor, this will return a hard-coded
  // file from $RDBASE.
  std::string getFontFile() const override;
  void setFontFile(const std::string &font_file) override;

 protected:
  double fontCoordToDrawCoord(FT_Pos fc) const;
  void fontPosToDrawPos(FT_Pos fx, FT_Pos fy, double &dx, double &dy) const;
  // adds x_trans_ and y_trans_ to coords returns x advance distance
  virtual double extractOutline();

 private:
  FT_Library library_;
  FT_Face face_;
  std::string font_file_;  // over-rides default if not empty.
  double x_trans_, y_trans_;
  mutable FT_Pos
      string_y_max_;  // maximum y value of string drawn, for inverting y
  double em_scale_;

  // return a vector of StringRects, one for each char in text, with
  // super- and subscripts taken into account.  Sizes in pixel coords,
  // i.e. scaled by fontScale().
  void getStringRects(const std::string &text,
                      std::vector<std::shared_ptr<StringRect>> &rects,
                      std::vector<TextDrawType> &draw_modes,
                      std::vector<char> &draw_chars) const override;

  // calculate the bounding box of the glyph for c in
  // font units (0 -> face_->units_per_EM (2048 for roboto font).
  void calcGlyphBBox(char c, FT_Pos &x_min, FT_Pos &y_min, FT_Pos &x_max,
                     FT_Pos &y_max, FT_Pos &advance) const;
};

// Callbacks for FT_Outline_Decompose.  user should be a pointer to
// an instance of DrawTextFT.
int moveToFunction(const FT_Vector *to, void *user);
int lineToFunction(const FT_Vector *to, void *user);
int conicToFunction(const FT_Vector *control, const FT_Vector *to, void *user);
int cubicToFunction(const FT_Vector *controlOne, const FT_Vector *controlTwo,
                    const FT_Vector *to, void *user);

}  // namespace RDKit

#endif  // RDKIT_DRAWTEXTFT_H
