//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Original author: David Cosgrove (CozChemIx) on 29/04/2020.
//
// This is an abstract base class for drawing text into a MolDraw2D
// object.

#ifndef RDKIT_DRAWTEXT_H
#define RDKIT_DRAWTEXT_H

#include <string>
#include <vector>

#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

using RDGeom::Point2D;

namespace RDKit {

// for aligning the drawing of text to the passed in coords.
enum class OrientType:unsigned char { C = 0, N, E, S, W };
enum class TextAlignType:unsigned char { START, MIDDLE, END };
enum class TextDrawType:unsigned char  {
  TextDrawNormal = 0,
  TextDrawSuperscript,
  TextDrawSubscript
};
std::ostream& operator<<(std::ostream &oss, const TextAlignType &tat);
std::ostream& operator<<(std::ostream &oss, const TextDrawType &tdt);
std::ostream& operator<<(std::ostream &oss, const OrientType &o);

// ****************************************************************************
class DrawText {

 public:

  static constexpr double FONT_SIZE = 0.4; // based on a bond length of 1

  DrawText(double max_fnt_sz);
  virtual ~DrawText() {}

  DrawColour const &colour() const;
  void setColour(const DrawColour &col);

  virtual double fontSize() const;
  double maxFontSize() const;
  void setMaxFontSize(double new_max);
  double fontScale() const;
  void setFontScale(double new_scale);

  //! using the current scale, work out the size of the label
  /*!
     Bear in mind when implementing this, that, for example, NH2 will appear as
     NH<sub>2</sub> to convey that the 2 is a subscript, and this needs to
     accounted for in the width and height.
   */
  virtual void getStringSize(const std::string &label, double &label_width,
                             double &label_height) const;
  // returns a rectangle that goes round the label, in pixel coords.  Centred
  // on origin.
  StringRect getStringRect(const std::string &label, OrientType orient) const;

  //! drawString centres the string on cds.
  virtual void drawString(const std::string &str, const Point2D &cds,
                          TextAlignType align);
  // Aligns them according to OrientType.
  void drawString(const std::string &label, const Point2D &cds,
                  OrientType orient);

  // draw the char, with the bottom left hand corner at cds
  virtual void drawChar(char c, const Point2D &cds) = 0;

 protected:
  // amount to scale subscripts and superscripts by
  constexpr static double SUBS_SCALE = 0.75;
  constexpr static double SUPER_SCALE = 0.5;

  virtual void alignString(TextAlignType align,
                           const std::vector<TextDrawType> &draw_modes,
                           std::vector<std::shared_ptr<StringRect> > &rects) const;
  // adjust the string rectangles up and down for super- and subscripts
  void adjustStringRectsForSuperSubScript(const std::vector<TextDrawType> &draw_modes,
                                          const std::vector<char> &to_draw,
                                          std::vector<std::shared_ptr<StringRect>> &rects) const;
  // return a scale factor appropriate for the character and draw type
  // (normal or super- or subscript)
  double selectScaleFactor(char c, TextDrawType draw_type) const;

 private:

  DrawColour colour_;
  double font_scale_;
  double max_font_size_;

  // return a vector of StringRects, one for each char in text, with
  // super- and subscripts taken into account.  Sizes in pixel coords,
  // i.e. scaled by fontScale().
  virtual void getStringRects(const std::string &text,
                              std::vector<std::shared_ptr<StringRect>> &rects,
                              std::vector<TextDrawType> &draw_modes,
                              std::vector<char> &draw_chars) const = 0;
  void getStringRects(const std::string &text, OrientType orient,
                      std::vector<std::shared_ptr<StringRect>> &rects,
                      std::vector<TextDrawType> &draw_modes,
                      std::vector<char> &draw_chars);
  void drawRects(const Point2D &a_cds,
                 const std::vector<std::shared_ptr<StringRect>> &rects,
                 const std::vector<TextDrawType> &draw_modes,
                 const std::vector<char> &draw_chars);
};

//! establishes whether to put string draw mode into super- or sub-script
//! mode based on contents of instring from i onwards. Increments i
//! appropriately
//! \returns true or false depending on whether it did something or not
bool setStringDrawMode(const std::string &instring,
                       TextDrawType &draw_mode, size_t &i);

// take the label for the given atom and return the individual pieces
// that need to be drawn for it.  So NH<sub>2</sub> will return
// "N", "H<sub>2</sub>".
std::vector<std::string> atomLabelToPieces(const std::string &label,
                                           OrientType orient);

}  // namespace RDKit

#endif  // RDKIT_DRAWTEXT_H
