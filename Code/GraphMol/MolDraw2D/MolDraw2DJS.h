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
    initDrawing();
    initTextDrawer(noFreetype);
  };

  void drawLine(const Point2D &cds1, const Point2D &cds2) override;
  void drawPolygon(const std::vector<Point2D> &cds) override;
  void drawEllipse(const Point2D &cds1, const Point2D &cds2) override;
  void clearDrawing() override;

 private:
  emscripten::val &d_context;
  void initDrawing() override;
  void initTextDrawer(bool noFreetype) override;
};

}  // namespace RDKit
#endif  // MOLDRAW2DSVG_H
