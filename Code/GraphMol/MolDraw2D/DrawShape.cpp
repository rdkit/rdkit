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
                     bool scaleLineWidth, DrawColour lineColour, bool fill,
                     int atom1, int atom2, int bond)
    : points_(points),
      lineWidth_(lineWidth),
      scaleLineWidth_(scaleLineWidth),
      lineColour_(lineColour),
      fill_(fill),
      atom1_(atom1),
      atom2_(atom2),
      bond_(bond) {}

// ****************************************************************************
void DrawShape::draw(MolDraw2D &drawer) {
  // the various myDraw functions may over-ride these.
  const auto ocolour = drawer.colour();
  drawer.setColour(lineColour_);
  const auto olw = drawer.lineWidth();
  drawer.setLineWidth(lineWidth_);
  const auto ofill = drawer.fillPolys();
  drawer.setFillPolys(fill_);
  const auto lineScale = drawer.drawOptions().scaleBondWidth;
  drawer.drawOptions().scaleBondWidth = scaleLineWidth_;
  drawer.setActiveAtmIdx(atom1_, atom2_);
  drawer.setActiveBndIdx(bond_);

  myDraw(drawer);

  drawer.setActiveAtmIdx();
  drawer.setActiveBndIdx();
  drawer.setColour(ocolour);
  drawer.setLineWidth(olw);
  drawer.setFillPolys(ofill);
  drawer.drawOptions().scaleBondWidth = lineScale;
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
void DrawShape::scale(const Point2D &scale_factor) {
  for (auto &p : points_) {
    p.x *= scale_factor.x;
    p.y *= scale_factor.y;
  }
}

// ****************************************************************************
void DrawShape::move(const Point2D &trans) {
  for (auto &p : points_) {
    p += trans;
  }
}

// ****************************************************************************
bool DrawShape::doesRectClash(const StringRect &rect, double padding) const {
  return false;
}

// ****************************************************************************
DrawShapeArrow::DrawShapeArrow(const std::vector<Point2D> &points,
                               int lineWidth, bool scaleLineWidth,
                               DrawColour lineColour, bool fill,
                               double frac, double angle)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill),
      frac_(frac),
      angle_(angle) {
  PRECONDITION(points_.size() == 2, "arrow bad points size");
}

// ****************************************************************************
void DrawShapeArrow::myDraw(MolDraw2D &drawer) const {
  if (drawer.drawOptions().splitBonds) {
    drawer.setActiveAtmIdx(atom1_);
  }
  drawer.drawLine(points_[0], points_[1], lineColour_, lineColour_);
  if (drawer.drawOptions().splitBonds) {
    drawer.setActiveAtmIdx(atom2_);
  }
  if (!fill_) {
    drawer.drawLine(points_[1], points_[2], lineColour_, lineColour_);
    drawer.drawLine(points_[1], points_[3], lineColour_, lineColour_);
  } else {
    drawer.setFillPolys(true);
    std::vector<Point2D> head(points_.begin() + 1, points_.end());
    drawer.drawPolygon(head);
  }
}

// ****************************************************************************
bool DrawShapeArrow::doesRectClash(const StringRect &rect,
                                   double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  if (MolDraw2D_detail::doesLineIntersect(rect, points_[0], points_[1],
                                          padding)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersect(rect, points_[1], points_[2],
                                          padding)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersect(rect, points_[1], points_[3],
                                          padding)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersect(rect, points_[2], points_[3],
                                          padding)) {
    return true;
  }
  return false;
}

// ****************************************************************************
DrawShapeEllipse::DrawShapeEllipse(const std::vector<Point2D> &points,
                                   int lineWidth, bool scaleLineWidth,
                                   DrawColour lineColour, bool fill,
                                   int atom1)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill, atom1) {
  PRECONDITION(points_.size() == 2, "ellipse wrong points");
}

// ****************************************************************************
void DrawShapeEllipse::myDraw(MolDraw2D &drawer) const {
  if (fill_) {
    drawer.setLineWidth(1);
    drawer.drawOptions().scaleBondWidth = false;
  }
  drawer.drawEllipse(points_[0], points_[1]);
}

// ****************************************************************************
void DrawShapeEllipse::findExtremes(double &xmin, double &xmax, double &ymin,
                                    double &ymax) const {
  double wb2 = (points_[1].x - points_[0].x) / 2;
  double hb2 = (points_[1].y - points_[0].y) / 2;
  double cx = points_[0].x + wb2;
  double cy = points_[0].y + hb2;
  wb2 = wb2 > 0 ? wb2 : -1 * wb2;
  hb2 = hb2 > 0 ? hb2 : -1 * hb2;
  xmin = std::min(cx - wb2, xmin);
  xmax = std::max(cx + wb2, xmax);
  ymin = std::min(cy - hb2, ymin);
  ymax = std::max(cy + hb2, ymax);
}

// ****************************************************************************
bool DrawShapeEllipse::doesRectClash(const StringRect &rect,
                                     double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  Point2D tl, tr, br, bl;
  rect.calcCorners(tl, tr, br, bl, padding);
  double w = points_[1].x - points_[0].x;
  double h = points_[1].y - points_[0].y;
  double cx = points_[0].x + w / 2;
  double cy = points_[0].y + h / 2;
  w = w > 0 ? w : -1 * w;
  h = h > 0 ? h : -1 * h;
  Point2D centre{cx, cy};
  if (MolDraw2D_detail::doesLineIntersectEllipse(centre, w / 2.0, h / 2.0,
                                                 padding, tl, tr)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersectEllipse(centre, w / 2.0, h / 2.0,
                                                 padding, tr, br)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersectEllipse(centre, w / 2.0, h / 2.0,
                                                 padding, br, bl)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersectEllipse(centre, w / 2.0, h / 2.0,
                                                 padding, bl, tl)) {
    return true;
  }
  return false;
}

// ****************************************************************************
DrawShapeSimpleLine::DrawShapeSimpleLine(const std::vector<Point2D> &points,
                                         int lineWidth, bool scaleLineWidth,
                                         DrawColour lineColour,
                                         int atom1, int atom2, int bond,
                                         DashPattern dashPattern)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, false, atom1,
                atom2, bond),
      dashPattern_(dashPattern) {
  PRECONDITION(points_.size() == 2, "simple line wrong number of points");
}

// ****************************************************************************
void DrawShapeSimpleLine::myDraw(MolDraw2D &drawer) const {
  auto od = drawer.dash();
  auto dp = dashPattern_;
  if (dp == shortDashes) {
    // these are roughly equivalent to the original checks on scale,
    // which we don't have any more.
    double sq_len = (points_[0] - points_[1]).lengthSq();
    if (sq_len < 55.0) {
      dp[0] /= 4;
      dp[1] /= 3;
    } else if(sq_len < 900.0) {
      dp[0] /= 2;
      dp[1] /= 1.5;
    }
  }
  drawer.setDash(dp);
  drawer.setActiveAtmIdx(atom1_, atom2_);
  drawer.setActiveBndIdx(bond_);
  drawer.drawLine(points_[0], points_[1], lineColour_, lineColour_);
  drawer.setActiveAtmIdx();
  drawer.setActiveBndIdx();
  drawer.setDash(od);
}

// ****************************************************************************
bool DrawShapeSimpleLine::doesRectClash(const StringRect &rect,
                                        double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  return MolDraw2D_detail::doesLineIntersect(rect, points_[0],
                                             points_[1], padding);
}

// ****************************************************************************
DrawShapePolyLine::DrawShapePolyLine(const std::vector<Point2D> &points,
                                     int lineWidth, bool scaleLineWidth,
                                     DrawColour lineColour, bool fill,
                                     int atom1, int atom2,
                                     int bond, DashPattern dashPattern)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill,
                atom1, atom2, bond),
      dashPattern_(dashPattern) {
  PRECONDITION(points_.size() > 2, "polyline not enough points");
}

// ****************************************************************************
void DrawShapePolyLine::myDraw(MolDraw2D &drawer) const {
  drawer.drawPolygon(points_);
}

// ****************************************************************************
bool DrawShapePolyLine::doesRectClash(const StringRect &rect,
                                      double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  Point2D tl, tr, br, bl;
  rect.calcCorners(tl, tr, br, bl, padding);
  for (size_t i=1; i < points_.size(); ++i) {
    if (MolDraw2D_detail::doesLineIntersect(rect, points_[i - 1],
                                            points_[i], padding)) {
      return true;
    }
  }
  if (MolDraw2D_detail::doesLineIntersect(rect, points_.front(),
                                          points_.back(), padding)) {
    return true;
  }
  return false;
}

// ****************************************************************************
DrawShapeSolidWedge::DrawShapeSolidWedge(const std::vector<Point2D> points,
                                         const DrawColour &col1,
                                         const DrawColour &col2, bool inverted,
                                         bool splitBonds, int atom1, int atom2,
                                         int bond)
    : DrawShape(points, 1, false, col1, false, atom1, atom2, bond),
      col2_(col2),
      inverted_(inverted),
      splitBonds_(splitBonds) {
  PRECONDITION(points_.size() == 3, "solid wedge wrong points");
  buildTriangles();
}

// ****************************************************************************
void DrawShapeSolidWedge::buildTriangles() {
  if (!(lineColour_ == col2_) || splitBonds_) {
    Point2D point = points_[0];
    Point2D end1 = points_[1];
    Point2D end2 = points_[2];
    points_.clear();
    Point2D e1 = end1 - point;
    Point2D e2 = end2 - point;
    Point2D mid1 = point + e1 * 0.5;
    Point2D mid2 = point + e2 * 0.5;
    points_.emplace_back(point);
    points_.emplace_back(mid1);
    points_.emplace_back(mid2);
    points_.emplace_back(mid1);
    points_.emplace_back(end2);
    points_.emplace_back(end1);
    points_.emplace_back(mid1);
    points_.emplace_back(mid2);
    points_.emplace_back(end2);
  }
}

// ****************************************************************************
void DrawShapeSolidWedge::myDraw(MolDraw2D &drawer) const {
  drawer.setFillPolys(true);
  int at1Idx = inverted_ ? atom2_ : atom1_;
  int at2Idx = inverted_ ? atom1_ : atom2_;
  if (drawer.drawOptions().splitBonds) {
    drawer.setActiveAtmIdx(at1Idx);
  } else {
    drawer.setActiveAtmIdx(at1Idx, at2Idx);
  }
  drawer.setActiveBndIdx(bond_);
  drawer.drawTriangle(points_[0], points_[1], points_[2]);
  if (points_.size() > 3) {
    if (drawer.drawOptions().splitBonds) {
      drawer.setActiveAtmIdx(at2Idx);
    }
    drawer.setColour(col2_);
    drawer.drawTriangle(points_[3], points_[4], points_[5]);
    drawer.drawTriangle(points_[6], points_[7], points_[8]);
  }
}

// ****************************************************************************
void DrawShapeSolidWedge::scale(const Point2D &scale_factor) {
  DrawShape::scale(scale_factor);

  // don't let it get too wide.  A width of 12 pixels is similar to the
  // original 0.15 scaled by 40.
  Point2D point, end1, end2;
  point = points_[0];
  if (points_.size() == 3) {
    end1 = points_[1];
    end2 = points_[2];
  } else {
    end1 = points_[5];
    end2 = points_[4];
  }
  double widthsq = (end1 - end2).lengthSq();
  std::cout << "squared width = " << widthsq << std::endl;
  if (widthsq > 144.0) {
    points_ = calcScaledWedgePoints(points_[0], end1, end2, widthsq);
    buildTriangles();
  }

}

// ****************************************************************************
bool DrawShapeSolidWedge::doesRectClash(const StringRect &rect,
                                        double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  if (MolDraw2D_detail::doesTriangleIntersect(rect, points_[0],
                                              points_[1], points_[2], padding)) {
    return true;
  }
  if (points_.size() > 3) {
    if (MolDraw2D_detail::doesTriangleIntersect(
            rect, points_[3], points_[4], points_[5], padding)) {
      return true;
    }
    if (MolDraw2D_detail::doesTriangleIntersect(
            rect, points_[6], points_[7], points_[8], padding)) {
      return true;
    }
  }
  return false;
}

// ****************************************************************************
DrawShapeDashedWedge::DrawShapeDashedWedge(const std::vector<Point2D> points,
                                           const DrawColour &col1,
                                           const DrawColour &col2,
                                           bool inverted, int atom1, int atom2,
                                           int bond)
    : DrawShape(points, 1, false, col1, false, atom1, atom2, bond),
      col2_(col2),
      inverted_(inverted) {
  PRECONDITION(points_.size() == 3, "dashed wedge wrong points");
  at1Cds_ = points[0];
  buildLines();
}

// ****************************************************************************
void DrawShapeDashedWedge::buildLines() {
  // assumes this is the starting configuration, where the 3 points define
  // the enclosing triangle of the wedge.
  Point2D point = points_[0];
  Point2D end1 = points_[1];
  Point2D end2 = points_[2];
  points_.clear();
  lineColours_.clear();
  Point2D e1 = end1 - point;
  Point2D e2 = end2 - point;

  unsigned int nDashes = 6;
  double sideLen = e1.length();
  double dashSep = sideLen / nDashes;
  // don't have the dashes too far apart or too close together.
  if (dashSep > 20.0) {
    nDashes = sideLen / 20;
  } else if (dashSep < 5.0) {
    nDashes = sideLen / 5;
    nDashes = nDashes < 3 ? 3 : nDashes;
  }
  for (unsigned int i = 1; i < nDashes + 1; ++i) {
    Point2D e11 = point + e1 * (rdcast<double>(i) / nDashes);
    Point2D e22 = point + e2 * (rdcast<double>(i) / nDashes);
    points_.emplace_back(e11);
    points_.emplace_back(e22);
    if (i > nDashes / 2) {
      lineColours_.emplace_back(col2_);
    } else {
      lineColours_.emplace_back(lineColour_);
    }
  }
}

// ****************************************************************************
void DrawShapeDashedWedge::myDraw(MolDraw2D &drawer) const {
  drawer.setFillPolys(false);
  int at1Idx = inverted_ ? atom2_ : atom1_;
  int at2Idx = inverted_ ? atom1_ : atom2_;
  drawer.setActiveAtmIdx(at1Idx, at2Idx);
  drawer.setActiveBndIdx(bond_);
  for(size_t i=0, j=0; i < points_.size(); i+= 2, ++j) {
    if (drawer.drawOptions().splitBonds) {
      if (i < points_.size() / 2) {
        drawer.setActiveAtmIdx(at1Idx);
      } else {
        drawer.setActiveAtmIdx(at2Idx);
      }
    }
    drawer.drawLine(points_[i], points_[i + 1], lineColours_[j],
                    lineColours_[j]);
  }
}

// ****************************************************************************
void DrawShapeDashedWedge::scale(const Point2D &scale_factor) {
  DrawShape::scale(scale_factor);
  at1Cds_.x *= scale_factor.x;
  at1Cds_.y *= scale_factor.y;

  // don't let it get too wide.  A width of 12 pixels is similar to the
  // original 0.15 scaled by 40.
  Point2D point, end1, end2;
  size_t lastPoint = points_.size() - 1;
  point = at1Cds_;
  end1 = points_[lastPoint];
  end2 = points_[lastPoint -1];
  double widthsq = (end1 - end2).lengthSq();

  if (widthsq > 144.0) {
    points_ = calcScaledWedgePoints(point, end1, end2, widthsq);
  } else {
    points_ = std::vector<Point2D>{point, end1, end2};
  }
  // always rebuild the lines so we retain a sensible number of dashes.
  buildLines();
}

// ****************************************************************************
void DrawShapeDashedWedge::move(const Point2D &trans) {
  DrawShape::move(trans);
  at1Cds_ += trans;
}

// ****************************************************************************
bool DrawShapeDashedWedge::doesRectClash(const StringRect &rect,
                                         double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  size_t last_point = points_.size() - 1;
  if (MolDraw2D_detail::doesLineIntersect(rect, at1Cds_,
                                          points_[last_point], padding)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersect(rect, points_[last_point],
                                          points_[last_point - 1], padding)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersect(rect, points_[last_point - 1],
                                          at1Cds_, padding)) {
    return true;
  }
  return false;
}

// ****************************************************************************
DrawShapeWavyLine::DrawShapeWavyLine(const std::vector<Point2D> points,
                                     int lineWidth, bool scaleLineWidth,
                                     const DrawColour &col1,
                                     const DrawColour &col2, double offset,
                                     int atom1, int atom2, int bond)
    : DrawShape(points, lineWidth, scaleLineWidth, col1, false, atom1, atom2,
                bond),
      col2_(col2),
      offset_(offset) {
  PRECONDITION(points_.size() == 2, "wavy line wrong points");
}

// ****************************************************************************
void DrawShapeWavyLine::myDraw(MolDraw2D &drawer) const {
  // nSegments is 16 by default in MolDraw2D.
  drawer.drawWavyLine(points_[0], points_[1], lineColour_, col2_, 16, offset_);
}

// ****************************************************************************
bool DrawShapeWavyLine::doesRectClash(const StringRect &rect, double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  padding += offset_;
  return MolDraw2D_detail::doesLineIntersect(rect, points_[0],
                                             points_[1], padding);
}

// ****************************************************************************
std::vector<Point2D> calcScaledWedgePoints(const Point2D &point,
                                           const Point2D &end1,
                                           const Point2D &end2,
                                           double widthsq) {
  Point2D lastLine = end1 - end2;
  Point2D mid = (end1 + end2) / 2.0;
  lastLine /= sqrt(widthsq);
  lastLine *= 6;
  std::vector<Point2D> new_pts{point, mid + lastLine, mid - lastLine};
  return new_pts;
}

// ****************************************************************************
DrawShapeArc::DrawShapeArc(const std::vector<Point2D> points, double ang1,
                           double ang2, int lineWidth, bool scaleLineWidth,
                           const DrawColour &col1, bool fill, int atom1)
    : DrawShape(points, lineWidth, scaleLineWidth, col1, fill, atom1),
      ang1_(ang1),
      ang2_(ang2) {
  PRECONDITION(points_.size() == 2, "arc wrong points");
  PRECONDITION(ang1_ >= 0.0 && ang1_ <= 360.0, "ang1_ not 0-360")
  PRECONDITION(ang2_ >= 0.0 && ang2_ <= 360.0, "ang2_ not 0-360")
}

// ****************************************************************************
void DrawShapeArc::myDraw(MolDraw2D &drawer) const {
  double start_ang = ang1_ > ang2_ ? ang1_ - 360.0 : ang1_;
  drawer.drawArc(points_[0], points_[1].x, points_[1].y, start_ang, ang2_);
}

// ****************************************************************************
void DrawShapeArc::findExtremes(double &xmin, double &xmax, double &ymin, double &ymax) const {
  xmin = std::min(xmin, points_[0].x - points_[1].x);
  xmax = std::max(xmax, points_[0].x + points_[1].x);
  ymin = std::min(ymin, points_[0].y - points_[1].y);
  ymax = std::max(ymax, points_[0].y + points_[1].y);
}

// ****************************************************************************
void DrawShapeArc::move(const Point2D &trans) {
  points_[0] += trans;
}

// ****************************************************************************
bool DrawShapeArc::doesRectClash(const StringRect &rect,
                                 double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  Point2D tl, tr, br, bl;
  rect.calcCorners(tl, tr, br, bl, padding);
  if (MolDraw2D_detail::doesLineIntersectArc(points_[0], points_[1].x,
                                             points_[1].y, ang1_, ang2_,
                                             padding, tl, tr)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersectArc(points_[0], points_[1].x,
                                             points_[1].y, ang1_, ang2_,
                                             padding, tr, br)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersectArc(points_[0], points_[1].x,
                                             points_[1].y, ang1_, ang2_,
                                             padding, br, bl)) {
    return true;
  }
  if (MolDraw2D_detail::doesLineIntersectArc(points_[0], points_[1].x,
                                             points_[1].y, ang1_, ang2_,
                                             padding, bl, tl)) {
    return true;
  }
  return false;
}

} // namespace RDKit