//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (AstraZeneca)
//
// This is a concrete class derived from MolDraw2D that uses RDKit to draw a
// molecule into a QPainter.

#include <RDGeneral/export.h>
#ifndef MOLDRAW2DQT_H
#define MOLDRAW2DQT_H

#include "MolDraw2D.h"

class QPainter;
class QString;

// ****************************************************************************

namespace RDKit {

class RDKIT_MOLDRAW2D_EXPORT MolDraw2DQt : public MolDraw2D {
 public:
  MolDraw2DQt(int width, int height, QPainter &qp, int panelWidth = -1,
              int panelHeight = -1);

  // set font size in molecule coordinate units. That's probably Angstrom for
  // RDKit. It will turned into drawing units using scale_, which might be
  // changed as a result, to make sure things still appear in the window.
  void setFontSize(double new_size);
  void setColour(const DrawColour &col);

  void drawLine(const Point2D &cds1, const Point2D &cds2);
  void drawChar(char c, const Point2D &cds);
  void drawPolygon(const std::vector<Point2D> &cds);
  void clearDrawing();

  // using the current scale, work out the size of the label in molecule
  // coordinates
  void getStringSize(const std::string &label, double &label_width,
                     double &label_height) const;

 private:
  QPainter &qp_;
};
}  // namespace RDKit
#endif  // MOLDRAW2DQT_H
