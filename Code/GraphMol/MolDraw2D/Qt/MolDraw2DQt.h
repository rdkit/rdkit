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

#include <GraphMol/MolDraw2D/MolDraw2D.h>

class QPainter;
class QString;

// ****************************************************************************

namespace RDKit {

RDKIT_MOLDRAW2DQT_EXPORT extern const char *rdkitQtVersion;

class RDKIT_MOLDRAW2DQT_EXPORT MolDraw2DQt : public MolDraw2D {
 public:
  MolDraw2DQt(int width, int height, QPainter *qp, int panelWidth = -1,
              int panelHeight = -1, bool noFreetype = false);

  void setColour(const DrawColour &col) override;

  void drawLine(const Point2D &cds1, const Point2D &cds2,
                bool rawCoords = false) override;
  void drawChar(char c, const Point2D &cds);
  void drawPolygon(const std::vector<Point2D> &cds,
                   bool rawCoords = false) override;
  void clearDrawing() override;
  void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                    const DrawColour &col1, const DrawColour &col2,
                    unsigned int nSegments = 16, double vertOffset = 0.05,
                    bool rawCoords = false) override;

 private:
  QPainter *d_qp;
  void initDrawing() override;
  void initTextDrawer(bool noFreetype) override;
};
}  // namespace RDKit
#endif  // MOLDRAW2DQT_H
