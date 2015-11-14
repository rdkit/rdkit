//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// derived from Dave Cosgrove's MolDraw2D
//
// This is a concrete class derived from MolDraw2D that uses RDKit to draw a
// molecule into an SVG file

#ifndef MOLDRAW2DSVG_H
#define MOLDRAW2DSVG_H

#include <iostream>
#include <sstream>
#include "MolDraw2D.h"

// ****************************************************************************

namespace RDKit {

class MolDraw2DSVG : public MolDraw2D {
 public:
  // initialize to use a particular ostream
  MolDraw2DSVG(int width, int height, std::ostream &os)
      : MolDraw2D(width, height), d_os(os) {
    initDrawing();
  };
  // initialize to use the internal stringstream
  MolDraw2DSVG(int width, int height) : MolDraw2D(width, height), d_os(d_ss) {
    initDrawing();
  };

  // set font size in molecule coordinate units. That's probably Angstrom for
  // RDKit. It will turned into drawing units using scale_, which might be
  // changed as a result, to make sure things still appear in the window.
  void setFontSize(double new_size);
  void setColour(const DrawColour &col);

  // not sure if this goes here or if we should do a dtor since initDrawing() is
  // called in the ctor,
  // but we'll start here
  void finishDrawing();

  void drawLine(const Point2D &cds1, const Point2D &cds2);
  void drawString(const std::string &str, const Point2D &cds);
  void drawPolygon(const std::vector<Point2D> &cds);
  void drawEllipse(const Point2D &cds1, const Point2D &cds2);
  void clearDrawing();

  // using the current scale, work out the size of the label in molecule
  // coordinates
  void getStringSize(const std::string &label, double &label_width,
                     double &label_height) const;

  // this only makes sense if the object was initialized without a stream
  std::string getDrawingText() const { return d_ss.str(); };

  void tagAtoms(const ROMol &mol);

 private:
  std::ostream &d_os;
  std::stringstream d_ss;

  void drawChar(char c, const Point2D &cds);
  void initDrawing();
};
}
#endif  // MOLDRAW2DSVG_H
