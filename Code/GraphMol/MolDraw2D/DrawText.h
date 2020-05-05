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

#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

using RDGeom::Point2D;

namespace RDKit {

// ****************************************************************************
class DrawText {

 public:

  DrawText();
  virtual ~DrawText() {}

  DrawColour colour() const;
  void setColour(const DrawColour &col);

  double fontSize() const;
  void setFontSize(double new_size);

  //! using the current scale, work out the size of the label
  /*!
     Bear in mind when implementing this, that, for example, NH2 will appear as
     NH<sub>2</sub> to convey that the 2 is a subscript, and this needs to
     accounted for in the width and height.
   */
  virtual void getStringSize(const std::string &label, double &label_width,
                             double &label_height) const = 0;

  //! drawString centres the string on cds.
  virtual void drawString(const std::string &str, const Point2D &cds,
                          MolDraw2D::AlignType align);
  // draw the char, with the bottom left hand corner at cds
  virtual void drawChar(char c, const Point2D &cds) = 0;

 private:

  DrawColour colour_;
  double font_size_;

};

//! establishes whether to put string draw mode into super- or sub-script
//! mode based on contents of instring from i onwards. Increments i
//! appropriately
//! \returns true or false depending on whether it did something or not
bool setStringDrawMode(const std::string &instring,
                       MolDraw2D::TextDrawType &draw_mode,
                       int &i);

}  // namespace RDKit

#endif  // RDKIT_DRAWTEXT_H
