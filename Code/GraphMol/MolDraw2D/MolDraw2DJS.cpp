//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

#ifdef __EMSCRIPTEN__

#include <GraphMol/MolDraw2D/MolDraw2DJS.h>

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <Geometry/point.h>
#ifdef RDK_BUILD_FREETYPE_SUPPORT
#include <GraphMol/MolDraw2D/DrawTextFTJS.h>
#else
#include <GraphMol/MolDraw2D/DrawTextJS.h>
#endif

#include <boost/format.hpp>
#include <sstream>

namespace RDKit {

// ****************************************************************************
void MolDraw2DJS::initDrawing() {}

// ****************************************************************************
void MolDraw2DJS::initTextDrawer(bool noFreetype) {
  double max_fnt_sz = drawOptions().maxFontSize;
  double min_fnt_sz = drawOptions().minFontSize;

  if (noFreetype) {
    text_drawer_.reset(new DrawTextJS(max_fnt_sz, min_fnt_sz, d_context));
  } else {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    try {
      text_drawer_.reset(new DrawTextFTJS(max_fnt_sz, min_fnt_sz,
                                          drawOptions().fontFile, d_context));
    } catch (std::runtime_error &e) {
      BOOST_LOG(rdWarningLog)
          << e.what() << std::endl
          << "Falling back to native JS text handling." << std::endl;
      text_drawer_.reset(new DrawTextJS(max_fnt_sz, min_fnt_sz, d_context));
    }
#else
    text_drawer_.reset(new DrawTextJS(max_fnt_sz, min_fnt_sz, d_context));
#endif
  }
}

// ****************************************************************************
void MolDraw2DJS::drawLine(const Point2D &cds1, const Point2D &cds2) {
  Point2D c1 = getDrawCoords(cds1);
  Point2D c2 = getDrawCoords(cds2);
  // std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();
  std::string dashString = "";
  const DashPattern &dashes = dash();
}

// ****************************************************************************
void MolDraw2DJS::drawPolygon(const std::vector<Point2D> &cds) {
  PRECONDITION(cds.size() >= 3, "must have at least three points");

  // std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();
}

// ****************************************************************************
void MolDraw2DJS::drawEllipse(const Point2D &cds1, const Point2D &cds2) {
  Point2D c1 = getDrawCoords(cds1);
  Point2D c2 = getDrawCoords(cds2);
  double w = c2.x - c1.x;
  double h = c2.y - c1.y;
  double cx = c1.x + w / 2;
  double cy = c1.y + h / 2;
  w = w > 0 ? w : -1 * w;
  h = h > 0 ? h : -1 * h;

  // std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();
}

// ****************************************************************************
void MolDraw2DJS::clearDrawing() {
  // std::string col = DrawColourToSVG(drawOptions().backgroundColour);
}

}  // namespace RDKit
#endif  // __EMSCRIPTEN__
