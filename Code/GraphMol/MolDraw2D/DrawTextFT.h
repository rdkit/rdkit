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

  DrawTextFT();
  ~DrawTextFT();

  double fontSize() const override;
  void getStringSize(const std::string &label, double &label_width,
                     double &label_height) const override;
  //! drawString centres the string on cds.
  void drawString(const std::string &str, const Point2D &cds,
                  TextAlignType align) override;

  void drawChar(char c, const Point2D &cds) override;

  virtual int MoveToFunctionImpl(const FT_Vector *to) = 0;
  virtual int LineToFunctionImpl(const FT_Vector *to) = 0;
  virtual int ConicToFunctionImpl(const FT_Vector *control,
                                  const FT_Vector *to) = 0;
  virtual int CubicToFunctionImpl(const FT_Vector *controlOne,
                                  const FT_Vector *controlTwo,
                                  const FT_Vector *to) = 0;

 protected:

  void alignString(const std::string &str, const Point2D &in_cds,
                   TextAlignType align, Point2D &out_cds) override;
  double fontCoordToPixelCoord(FT_Pos fc) const;
  void fontPosToDrawPos(FT_Pos fx, FT_Pos fy, double &dx, double &dy) const;
  // adds x_trans_ and y_trans_ to coords returns x advance distance
  virtual double extractOutline();

 private:

  FT_Library library_;
  FT_Face face_;
  FT_Pos x_trans_, y_trans_;
  mutable FT_Pos string_y_max_; // maximum y value of string drawn, for inverting y

  // resolution of the target display in dpi.  96 is normal for a computer
  // display, apparently, although I have roughly 160 and 224 on mine.
  constexpr static int RESOLUTION = 96;
  // define a 16 point font.  It will be scaled up and down as required.
  // A point is 1/72 of an inch, so it will only appear at exactly that
  // size on a display of RESOLUTION dpi.
  constexpr static int POINT_SIZE = 16;
  // amount to scale subscripts and superscripts by
  constexpr static double SUBS_SCALE = 0.75;
  constexpr static double SUPER_SCALE = 0.5;

  // look for a hard-coded font file under RDBase.  We may want to
  // improve on this to give the user scope for having a different
  // one.
  std::string getFontFile() const;

  void getStringBBox(const std::string &text,
                     double &x_min, double &y_min ,
                     double &x_max, double &y_max) const;
  // return a vector of StringRects, one for each char in text, with
  // super- and subscripts taken into account.  Sizes in pixel coords,
  // i.e. scaled by fontScale().
  void getStringRects(const std::string &text,
                      std::vector<std::shared_ptr<StringRect>> &rects,
                      std::vector<TextDrawType> &draw_modes) const;

    // calculate the bounding box of the currently loaded glyph in
  // font units (0 -> face_->units_per_EM (2048 for roboto font).
  void calcGlyphBBox(char c, FT_Pos &x_min, FT_Pos &y_min,
                     FT_Pos &x_max, FT_Pos &y_max, FT_Pos &advance) const;

};

// Callbacks for FT_Outline_Decompose.  user should be a pointer to
// an instance of RDFreeType.
int moveToFunction(const FT_Vector *to, void *user);
int lineToFunction(const FT_Vector *to, void *user);
int conicToFunction(const FT_Vector *control, const FT_Vector *to,
                    void *user);
int cubicToFunction(const FT_Vector *controlOne,
                    const FT_Vector *controlTwo,
                    const FT_Vector *to, void *user);

} // namespace RDKit

#endif  // RDKIT_DRAWTEXTFT_H
