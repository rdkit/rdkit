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
#endif
#include <GraphMol/MolDraw2D/DrawTextJS.h>

#include <boost/format.hpp>
#include <sstream>
#include <cmath>

namespace RDKit {

std::string DrawColourToSVG(const RDKit::DrawColour &col);

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
  std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();
  std::string dashString = "";
  const DashPattern &dashes = dash();

  d_context.call<void>("beginPath");
  d_context.set("lineWidth", width);
  d_context.set("strokeStyle", col);
  if (dashes.size()) {
    d_context.call<void>("setLineDash", emscripten::typed_memory_view(
                                            dashes.size(), dashes.data()));
  }
  d_context.call<void>("moveTo", std::round(c1.x), std::round(c1.y));
  d_context.call<void>("lineTo", std::round(c2.x), std::round(c2.y));
  d_context.call<void>("stroke");
  if (dashes.size()) {
    static const DashPattern nodash;
    d_context.call<void>("setLineDash", emscripten::typed_memory_view(
                                            nodash.size(), nodash.data()));
  }
}
// ****************************************************************************
void MolDraw2DJS::drawPolygon(const std::vector<Point2D> &cds) {
  PRECONDITION(cds.size() >= 3, "must have at least three points");

  std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();

  d_context.call<void>("beginPath");
  d_context.set("lineWidth", width);
  d_context.set("lineCap", std::string("butt"));
  d_context.set("lineJoin", std::string("round"));
  d_context.set("strokeStyle", col);
  Point2D c0 = getDrawCoords(cds[0]);
  d_context.call<void>("moveTo", std::round(c0.x), std::round(c0.y));
  for (unsigned int i = 1; i < cds.size(); ++i) {
    Point2D ci = getDrawCoords(cds[i]);
    d_context.call<void>("lineTo", std::round(ci.x), std::round(ci.y));
  }
  if (fillPolys()) {
    d_context.call<void>("closePath");
    d_context.set("fillStyle", col);
    d_context.call<void>("fill");
  }
  d_context.call<void>("stroke");
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

  std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();
  d_context.call<void>("beginPath");
  d_context.set("lineWidth", width);
  d_context.set("strokeStyle", col);
  d_context.call<void>("ellipse", cx, cy, w / 2, h / 2, 0, 0, 2 * M_PI);
  if (fillPolys()) {
    d_context.set("fillStyle", col);
    d_context.call<void>("fill");
  }
  d_context.call<void>("stroke");
}

// ****************************************************************************
void MolDraw2DJS::clearDrawing() {
  std::string col = DrawColourToSVG(drawOptions().backgroundColour);
  d_context.call<void>("clearRect", offset().x, offset().y, width(), height());
  d_context.set("fillStyle", col);
  d_context.call<void>("fillRect", offset().x, offset().y, width(), height());
}

}  // namespace RDKit
#endif  // __EMSCRIPTEN__
