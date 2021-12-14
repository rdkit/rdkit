//
//  Copyright (C) 2014-2021 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <GraphMol/MolDraw2D/DrawShape.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>

namespace RDKit {
// ****************************************************************************
DrawShape::DrawShape(const std::vector<Point2D> &points, int lineWidth,
                     bool scaleLineWidth, DrawColour lineColour, bool fill)
    : points_(points),
      lineWidth_(lineWidth),
      scaleLineWidth_(scaleLineWidth),
      lineColour_(lineColour),
      fill_(fill) {}

// ****************************************************************************
void DrawShape::draw(MolDraw2D &drawer) {
  const auto ocolour = drawer.colour();
  const auto olw = drawer.lineWidth();
  const auto ofill = drawer.fillPolys();

  myDraw(drawer);

  drawer.setColour(ocolour);
  drawer.setLineWidth(olw);
  drawer.setFillPolys(ofill);
}

// ****************************************************************************
void DrawShape::findExtremes(double &xmin, double &xmax,
                             double &ymin, double &ymax) const {
  for(auto &p: points_) {
    xmin = std::min(xmin, p.x);
    xmax = std::max(xmax, p.x);
    ymin = std::min(ymin, p.y);
    ymax = std::max(ymax, p.y);
  }
}

// ****************************************************************************
void DrawShape::scale(const Point2D &scale) {
  for (auto &p : points_) {
    p.x *= scale.x;
    p.y *= scale.y;
  }
}

// ****************************************************************************
void DrawShape::move(const Point2D &trans) {
  for (auto &p : points_) {
    p += trans;
  }
}

// ****************************************************************************
DrawShapeArrow::DrawShapeArrow(const std::vector<Point2D> &points,
                               int lineWidth, bool scaleLineWidth,
                               DrawColour lineColour, bool fill)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill) {
  PRECONDITION(points_.size() == 4, "arrow bad points size");
}

// ****************************************************************************
void DrawShapeArrow::myDraw(MolDraw2D &drawer) const {
  drawer.setColour(lineColour_);
  drawer.setLineWidth(lineWidth_);
  drawer.drawLine(points_[0], points_[1]);
  if (!fill_) {
    drawer.drawLine(points_[1], points_[2]);
    drawer.drawLine(points_[1], points_[3]);
  } else {
    drawer.setFillPolys(true);
    std::vector<Point2D> head(points_.begin() + 1, points_.end());
    drawer.drawPolygon(head);
  }
}

// ****************************************************************************
DrawShapeEllipse::DrawShapeEllipse(const std::vector<Point2D> &points,
                                   int lineWidth, bool scaleLineWidth,
                                   DrawColour lineColour, bool fill)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill) {
  PRECONDITION(points_.size() == 2, "ellipse wrong points");
  Point2D p1(points_[0]);
  Point2D p2(points_[1]);
  points_.clear();
  MolDraw2D_detail::arcPoints(p1, p2, points_, 0.0, 360.0);
}

// ****************************************************************************
void DrawShapeEllipse::myDraw(MolDraw2D &drawer) const {
  drawer.setColour(lineColour_);
  drawer.setLineWidth(lineWidth_);
  drawer.setFillPolys(fill_);
  drawer.drawPolygon(points_);
}

// ****************************************************************************
DrawShapePolyline::DrawShapePolyline(const std::vector<Point2D> &points,
                                     int lineWidth, bool scaleLineWidth,
                                     DrawColour lineColour, bool fill)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill){
  PRECONDITION(points_.size() > 1, "polyline not enough points");
}

// ****************************************************************************
void DrawShapePolyline::myDraw(MolDraw2D &drawer) const {
  drawer.setColour(lineColour_);
  auto lw = lineWidth_;
  if (scaleLineWidth_) {
    lw *= drawer.scale() * 0.02;
    if (lw < 0.0) {
      lw = 0.0;
    }
  }
  drawer.setLineWidth(lw);
  if (points_.size() > 2 && fill_) {
    drawer.setFillPolys(true);
  } else {
    drawer.setFillPolys(false);
  }
  if (drawer.drawOptions().comicMode) {
    auto drawPoints =
        MolDraw2D_detail::handdrawnLine(points_[0], points_[1], drawer.scale());
    for (unsigned int i = 2; i < points_.size(); ++i) {
      auto lpts = MolDraw2D_detail::handdrawnLine(
          points_[i - 1], points_[i], drawer.scale());
      std::move(lpts.begin(), lpts.end(), std::back_inserter(drawPoints));
    }
    drawer.drawPolygon(drawPoints);
  } else {
    if (points_.size() > 2) {
      drawer.drawPolygon(points_);
    } else {
      //      std::cout << "Draw line from " << points_[0] << " to " <<
      //      points_[1]
      //                << " in colour " << lineColour_.r << ", " <<
      //                lineColour_.g
      //                << ", " << lineColour_.b << "  width = " << lw <<
      //                std::endl;
      drawer.drawLine(points_[0], points_[1]);
    }
  }
}

} // namespace RDKit