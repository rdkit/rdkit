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
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>

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
    text_drawer_.reset(
        new MolDraw2D_detail::DrawTextJS(max_fnt_sz, min_fnt_sz, d_context));
  } else {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    try {
      text_drawer_.reset(new MolDraw2D_detail::DrawTextFTJS(
          max_fnt_sz, min_fnt_sz, drawOptions().fontFile, d_context));
    } catch (std::runtime_error &e) {
      BOOST_LOG(rdWarningLog)
          << e.what() << std::endl
          << "Falling back to native JS text handling." << std::endl;
      text_drawer_.reset(
          new MolDraw2D_detail::DrawTextJS(max_fnt_sz, min_fnt_sz, d_context));
    }
#else
    text_drawer_.reset(
        new MolDraw2D_detail::DrawTextJS(max_fnt_sz, min_fnt_sz, d_context));
#endif
  }
  if (drawOptions().baseFontSize > 0.0) {
    text_drawer_->setBaseFontSize(drawOptions().baseFontSize);
  }
}

// ****************************************************************************
void MolDraw2DJS::drawLine(const Point2D &cds1, const Point2D &cds2,
                           bool rawCoords) {
  Point2D c1 = rawCoords ? cds1 : getDrawCoords(cds1);
  Point2D c2 = rawCoords ? cds2 : getDrawCoords(cds2);
  std::string col = DrawColourToSVG(colour());
  double width = getDrawLineWidth();
  std::string dashString = "";
  const DashPattern &dashes = dash();

  d_context.call<void>("beginPath");
  d_context.set("lineWidth", width);
  d_context.set("strokeStyle", col);
  if (dashes.size()) {
    d_context.call<void>("setLineDash", emscripten::typed_memory_view(
                                            dashes.size(), dashes.data()));
  }
  d_context.call<void>("moveTo", c1.x, c1.y);
  d_context.call<void>("lineTo", c2.x, c2.y);
  d_context.call<void>("stroke");
  if (dashes.size()) {
    static const DashPattern nodash;
    d_context.call<void>("setLineDash", emscripten::typed_memory_view(
                                            nodash.size(), nodash.data()));
  }
}
// ****************************************************************************
void MolDraw2DJS::drawPolygon(const std::vector<Point2D> &cds, bool rawCoords) {
  PRECONDITION(cds.size() >= 3, "must have at least three points");

  std::string col = DrawColourToSVG(colour());
  double width = getDrawLineWidth();

  d_context.call<void>("beginPath");
  d_context.set("lineWidth", width);
  d_context.set("lineCap", std::string("butt"));
  d_context.set("lineJoin", std::string("round"));
  d_context.set("strokeStyle", col);
  Point2D c0 = rawCoords ? cds[0] : getDrawCoords(cds[0]);
  d_context.call<void>("moveTo", c0.x, c0.y);
  for (unsigned int i = 1; i < cds.size(); ++i) {
    Point2D ci = rawCoords ? cds[i] : getDrawCoords(cds[i]);
    d_context.call<void>("lineTo", ci.x, ci.y);
  }
  if (fillPolys()) {
    d_context.call<void>("closePath");
    d_context.set("fillStyle", col);
    d_context.call<void>("fill");
  }
  d_context.call<void>("stroke");
}

// ****************************************************************************
void MolDraw2DJS::drawEllipse(const Point2D &cds1, const Point2D &cds2,
                              bool rawCoords) {
  Point2D c1 = rawCoords ? cds1 : getDrawCoords(cds1);
  Point2D c2 = rawCoords ? cds2 : getDrawCoords(cds2);
  double w = c2.x - c1.x;
  double h = c2.y - c1.y;
  double cx = c1.x + w / 2;
  double cy = c1.y + h / 2;
  w = w > 0 ? w : -1 * w;
  h = h > 0 ? h : -1 * h;

  std::string col = DrawColourToSVG(colour());
  double width = getDrawLineWidth();
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
void MolDraw2DJS::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                               const DrawColour &col1, const DrawColour &,
                               unsigned int nSegments, double vertOffset,
                               bool rawCoords) {
  PRECONDITION(nSegments > 1, "too few segments");

  std::string col = DrawColourToSVG(colour());
  double width = getDrawLineWidth();

  auto segments =
      MolDraw2D_detail::getWavyLineSegments(cds1, cds2, nSegments, vertOffset);

  d_context.call<void>("beginPath");
  d_context.set("lineWidth", width);
  d_context.set("strokeStyle", col);

  auto c1 = std::get<0>(segments[0]);
  c1 = rawCoords ? c1 : getDrawCoords(c1);

  d_context.call<void>("moveTo", c1.x, c1.y);
  for (const auto &segment : segments) {
    auto cpt1 = std::get<1>(segment);
    cpt1 = rawCoords ? cpt1 : getDrawCoords(cpt1);
    auto cpt2 = std::get<2>(segment);
    cpt2 = rawCoords ? cpt2 : getDrawCoords(cpt2);
    auto segpt = std::get<3>(segment);
    segpt = rawCoords ? segpt : getDrawCoords(segpt);
    d_context.call<void>("bezierCurveTo", cpt1.x, cpt1.y, cpt2.x, cpt2.y,
                         segpt.x, segpt.y);
  }
  d_context.call<void>("stroke");
}

// ****************************************************************************
void MolDraw2DJS::clearDrawing() {
  MolDraw2D::clearDrawing();

  std::string col = DrawColourToSVG(drawOptions().backgroundColour);
  d_context.call<void>("clearRect", offset().x, offset().y, width(), height());
  d_context.set("fillStyle", col);
  d_context.call<void>("fillRect", offset().x, offset().y, width(), height());
}

}  // namespace RDKit
#endif  // __EMSCRIPTEN__
