//
//  Copyright (C) 2020 Greg Landrum and T5 informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// derived from Dave Cosgrove's MolDraw2D
//
// This is a concrete class derived from MolDraw2D that uses RDKit to draw using
// the JS canvas. This requires emscripten and is only intended for the RDKit
// Javascript builds

#include <RDGeneral/export.h>
#ifndef MOLDRAW2DJS_H
#define MOLDRAW2DJS_H

#include <iostream>
#include <sstream>
#include <emscripten.h>
#include <emscripten/val.h>

#include <GraphMol/MolDraw2D/MolDraw2D.h>

// ****************************************************************************

namespace RDKit {

class RDKIT_MOLDRAW2D_EXPORT MolDraw2DJS : public MolDraw2D {
 public:
  // initialize to use a particular ostream
  MolDraw2DJS(int width, int height, emscripten::val &context,
              int panelWidth = -1, int panelHeight = -1,
              bool noFreetype = false)
      : MolDraw2D(width, height, panelWidth, panelHeight), d_context(context) {
    PRECONDITION(width > 0, "bad width");
    PRECONDITION(height > 0, "bad height");
    initTextDrawer(noFreetype);
  }
  MolDraw2DJS(const MolDraw2DJS &) = delete;
  MolDraw2DJS(MolDraw2DJS &&) = delete;
  MolDraw2DJS &operator=(const MolDraw2DJS &) = delete;
  MolDraw2DJS &operator=(MolDraw2DJS &&) = delete;

  void drawLine(const Point2D &cds1, const Point2D &cds2,
                bool rawCoords = false) override;
  void drawPolygon(const std::vector<Point2D> &cds,
                   bool rawCoords = false) override;
  void drawEllipse(const Point2D &cds1, const Point2D &cds2,
                   bool rawCoords = false) override;
  void clearDrawing() override;
  void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                    const DrawColour &col1, const DrawColour &col2,
                    unsigned int nSegments = 16, double vertOffset = 0.05,
                    bool rawCoords = false) override;

 private:
  emscripten::val &d_context;
  void initDrawing() override;
  void initTextDrawer(bool noFreetype) override;
};

}  // namespace RDKit
#endif  // MOLDRAW2DSVG_H
