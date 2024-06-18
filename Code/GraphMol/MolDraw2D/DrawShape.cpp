//
//  Copyright (C) 2021-2022 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <cmath>

#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/DrawShape.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/StringRect.h>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawShape::DrawShape(const std::vector<Point2D> &points, double lineWidth,
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
void DrawShape::findExtremes(double &xmin, double &xmax, double &ymin,
                             double &ymax) const {
  for (const auto &p : points_) {
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
bool DrawShape::doesRectClash(const StringRect &, double) const {
  return false;
}

// ****************************************************************************
DrawShapeArrow::DrawShapeArrow(const std::vector<Point2D> &points,
                               double lineWidth, bool scaleLineWidth,
                               DrawColour lineColour, bool fill, int atom1,
                               int atom2, int bond, double frac, double angle)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill, atom1,
                atom2, bond),
      frac_(frac),
      angle_(angle) {
  PRECONDITION(points_.size() == 2, "arrow bad points size");

  // the two ends of the arrowhead are used for collision detection so we
  // may as well store them.
  Point2D ab(points_[1]), p1, p2;
  MolDraw2D_detail::calcArrowHead(ab, p1, p2, points_[0], fill_, frac_, angle_);
  points_.push_back(p1);
  points_.push_back(p2);
}

// ****************************************************************************
void DrawShapeArrow::myDraw(MolDraw2D &drawer) const {
  if (drawer.drawOptions().splitBonds) {
    drawer.setActiveAtmIdx(atom2_);
  } else {
    drawer.setActiveAtmIdx(atom1_, atom2_);
  }
  drawer.setActiveBndIdx(bond_);
  drawer.drawArrow(points_[0], points_[1], fill_, frac_, angle_, lineColour_,
                   true);
}

// ****************************************************************************
bool DrawShapeArrow::doesRectClash(const StringRect &rect,
                                   double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  if (doesLineIntersect(rect, points_[0], points_[1], padding)) {
    return true;
  }
  if (doesLineIntersect(rect, points_[1], points_[2], padding)) {
    return true;
  }
  if (doesLineIntersect(rect, points_[1], points_[3], padding)) {
    return true;
  }
  if (doesLineIntersect(rect, points_[2], points_[3], padding)) {
    return true;
  }
  return false;
}

// ****************************************************************************
DrawShapeEllipse::DrawShapeEllipse(const std::vector<Point2D> &points,
                                   double lineWidth, bool scaleLineWidth,
                                   DrawColour lineColour, bool fill, int atom1)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill, atom1) {
  PRECONDITION(points_.size() == 2, "ellipse wrong points");
}

// ****************************************************************************
void DrawShapeEllipse::myDraw(MolDraw2D &drawer) const {
  if (fill_) {
    drawer.setLineWidth(1);
    drawer.drawOptions().scaleBondWidth = false;
  }
  auto p1 = points_[0] - points_[1];
  auto p2 = points_[0] + points_[1];
  drawer.drawEllipse(p1, p2, true);
}

// ****************************************************************************
void DrawShapeEllipse::findExtremes(double &xmin, double &xmax, double &ymin,
                                    double &ymax) const {
  // points_[0] is the centre, points_[1] the radii
  xmin = std::min(points_[0].x - points_[1].x, xmin);
  xmax = std::max(points_[0].x + points_[1].x, xmax);
  ymin = std::min(points_[0].y - points_[1].y, ymin);
  ymax = std::max(points_[0].y + points_[1].y, ymax);
}

// ****************************************************************************
void DrawShapeEllipse::move(const Point2D &trans) { points_[0] += trans; }
// ****************************************************************************
bool DrawShapeEllipse::doesRectClash(const StringRect &rect,
                                     double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  Point2D tl, tr, br, bl;
  rect.calcCorners(tl, tr, br, bl, padding);
  auto w = points_[1].x;
  auto h = points_[1].y;
  auto cx = points_[0].x;
  auto cy = points_[0].y;
  w = w > 0 ? w : -1 * w;
  h = h > 0 ? h : -1 * h;
  Point2D centre{cx, cy};
  if (doesLineIntersectEllipse(centre, w / 2.0, h / 2.0, padding, tl, tr)) {
    return true;
  }
  if (doesLineIntersectEllipse(centre, w / 2.0, h / 2.0, padding, tr, br)) {
    return true;
  }
  if (doesLineIntersectEllipse(centre, w / 2.0, h / 2.0, padding, br, bl)) {
    return true;
  }
  if (doesLineIntersectEllipse(centre, w / 2.0, h / 2.0, padding, bl, tl)) {
    return true;
  }
  return false;
}

// ****************************************************************************
DrawShapeSimpleLine::DrawShapeSimpleLine(const std::vector<Point2D> &points,
                                         double lineWidth, bool scaleLineWidth,
                                         DrawColour lineColour, int atom1,
                                         int atom2, int bond,
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
    } else if (sq_len < 900.0) {
      dp[0] /= 2;
      dp[1] /= 1.5;
    }
  }
  drawer.setDash(dp);
  drawer.setActiveAtmIdx(atom1_, atom2_);
  drawer.setActiveBndIdx(bond_);
  drawer.drawLine(points_[0], points_[1], lineColour_, lineColour_, true);
  drawer.setActiveAtmIdx();
  drawer.setActiveBndIdx();
  drawer.setDash(od);
}

// ****************************************************************************
bool DrawShapeSimpleLine::doesRectClash(const StringRect &rect,
                                        double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  return doesLineIntersect(rect, points_[0], points_[1], padding);
}

// ****************************************************************************
DrawShapePolyLine::DrawShapePolyLine(const std::vector<Point2D> &points,
                                     double lineWidth, bool scaleLineWidth,
                                     DrawColour lineColour, bool fill,
                                     int atom1, int atom2, int bond,
                                     DashPattern dashPattern)
    : DrawShape(points, lineWidth, scaleLineWidth, lineColour, fill, atom1,
                atom2, bond),
      dashPattern_(dashPattern) {
  PRECONDITION(points_.size() > 2, "polyline not enough points");
}

// ****************************************************************************
void DrawShapePolyLine::myDraw(MolDraw2D &drawer) const {
  drawer.drawPolygon(points_, true);
}

// ****************************************************************************
bool DrawShapePolyLine::doesRectClash(const StringRect &rect,
                                      double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  Point2D tl, tr, br, bl;
  rect.calcCorners(tl, tr, br, bl, padding);
  for (size_t i = 1; i < points_.size(); ++i) {
    if (doesLineIntersect(rect, points_[i - 1], points_[i], padding)) {
      return true;
    }
  }
  return doesLineIntersect(rect, points_.front(), points_.back(), padding);
}

// ****************************************************************************
DrawShapeSolidWedge::DrawShapeSolidWedge(const std::vector<Point2D> points,
                                         const DrawColour &col1,
                                         const DrawColour &col2,
                                         bool splitBonds,
                                         std::vector<Point2D> &otherBondVecs,
                                         double lineWidth, int atom1, int atom2,
                                         int bond)
    : DrawShape(points, lineWidth / 2.0, false, col1, false, atom1, atom2,
                bond),
      col2_(col2),
      splitBonds_(splitBonds),
      otherBondVecs_(otherBondVecs) {
  PRECONDITION(points_.size() == 3, "solid wedge wrong points");
  if (otherBondVecs_.size() > 2) {
    trimOtherBondVecs();
  }
  if (otherBondVecs_.size() == 2) {
    orderOtherBondVecs();
  }
  buildTriangles();
}

// ****************************************************************************
void DrawShapeSolidWedge::buildTriangles() {
  if (!(lineColour_ == col2_) || splitBonds_) {
    buildTwoColorTriangles();
  } else {
    buildSingleColorTriangles();
  }
}

// ****************************************************************************
void DrawShapeSolidWedge::buildSingleColorTriangles() {
  auto point = points_[0];
  auto end1 = points_[1];
  auto end2 = points_[2];
  auto midEnd = (end1 + end2) / 2.0;
  auto adjend1 = end1;
  auto adjend2 = end2;
  points_.clear();
  // adjust adjend1 and adjend2 to line up with otherBondVecs_.
  if (otherBondVecs_.empty()) {
    points_.push_back(point);
    points_.push_back(adjend1);
    points_.push_back(adjend2);
  } else if (otherBondVecs_.size() == 1) {
    auto side1 = (end1 - point) * 2.0;
    if (!doLinesIntersect(point, point + side1, midEnd - otherBondVecs_[0],
                          midEnd + otherBondVecs_[0], &adjend1)) {
      adjend1 = end1;
    }
    auto side2 = (end2 - point) * 2.0;
    if (!doLinesIntersect(point, point + side2, midEnd - otherBondVecs_[0],
                          midEnd + otherBondVecs_[0], &adjend2)) {
      adjend2 = end2;
    }
    points_.push_back(point);
    points_.push_back(adjend1);
    points_.push_back(adjend2);
  } else if (otherBondVecs_.size() == 2) {
    auto side1 = (end1 - point) * 2.0;
    if (!doLinesIntersect(point, point + side1, midEnd - otherBondVecs_[0],
                          midEnd + otherBondVecs_[0], &adjend1)) {
      adjend1 = end1;
    }
    points_.push_back(point);
    points_.push_back(adjend1);
    points_.push_back(midEnd);
    auto side2 = (end2 - point) * 2.0;
    if (!doLinesIntersect(point, point + side2, midEnd - otherBondVecs_[1],
                          midEnd + otherBondVecs_[1], &adjend2)) {
      adjend2 = end2;
    }
    points_.push_back(point);
    points_.push_back(midEnd);
    points_.push_back(adjend2);
  }
}

// ****************************************************************************
void DrawShapeSolidWedge::buildTwoColorTriangles() {
  auto point = points_[0];
  auto end1 = points_[1];
  auto end2 = points_[2];
  auto midEnd = (end1 + end2) / 2.0;
  auto adjend1 = end1;
  auto adjend2 = end2;
  points_.clear();
  auto e1 = end1 - point;
  auto e2 = end2 - point;
  auto mid1 = point + e1 * 0.5;
  auto mid2 = point + e2 * 0.5;
  points_.push_back(point);
  points_.push_back(mid1);
  points_.push_back(mid2);
  if (otherBondVecs_.empty()) {
    points_.push_back(mid1);
    points_.push_back(adjend2);
    points_.push_back(adjend1);
    points_.push_back(mid1);
    points_.push_back(mid2);
    points_.push_back(adjend2);
  } else if (otherBondVecs_.size() == 1) {
    auto side1 = (end1 - point) * 2.0;
    if (!doLinesIntersect(point, point + side1, midEnd - otherBondVecs_[0],
                          midEnd + otherBondVecs_[0], &adjend1)) {
      adjend1 = end1;
    }
    auto side2 = (end2 - point) * 2.0;
    if (!doLinesIntersect(point, point + side2, midEnd - otherBondVecs_[0],
                          midEnd + otherBondVecs_[0], &adjend2)) {
      adjend2 = end2;
    }
    points_.push_back(mid1);
    points_.push_back(adjend2);
    points_.push_back(adjend1);
    points_.push_back(mid1);
    points_.push_back(mid2);
    points_.push_back(adjend2);
  } else if (otherBondVecs_.size() == 2) {
    auto side1 = (end1 - point) * 2.0;
    if (!doLinesIntersect(point, point + side1, midEnd - otherBondVecs_[0],
                          midEnd + otherBondVecs_[0], &adjend1)) {
      adjend1 = end1;
    }
    auto side2 = (end2 - point) * 2.0;
    if (!doLinesIntersect(point, point + side2, midEnd - otherBondVecs_[1],
                          midEnd + otherBondVecs_[1], &adjend2)) {
      adjend2 = end2;
    }
    points_.push_back(mid1);
    points_.push_back(adjend1);
    points_.push_back(midEnd);
    points_.push_back(midEnd);
    points_.push_back(mid2);
    points_.push_back(mid1);
    points_.push_back(midEnd);
    points_.push_back(adjend2);
    points_.push_back(adjend1);
  }
}

// ****************************************************************************
void DrawShapeSolidWedge::myDraw(MolDraw2D &drawer) const {
  drawer.setFillPolys(true);
  if (drawer.drawOptions().splitBonds) {
    drawer.setActiveAtmIdx(atom1_);
  } else {
    drawer.setActiveAtmIdx(atom1_, atom2_);
  }
  drawer.setActiveBndIdx(bond_);
  drawer.drawTriangle(points_[0], points_[1], points_[2], true);
  if (points_.size() > 3) {
    if (drawer.drawOptions().splitBonds) {
      drawer.setActiveAtmIdx(atom2_);
    }
    drawer.setColour(col2_);
  }
  for (unsigned int i = 3; i < points_.size(); i += 3) {
    drawer.drawTriangle(points_[i], points_[i + 1], points_[i + 2], true);
  }
}

// ****************************************************************************
bool DrawShapeSolidWedge::doesRectClash(const StringRect &rect,
                                        double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  if (doesTriangleIntersect(rect, points_[0], points_[1], points_[2],
                            padding)) {
    return true;
  }
  if (points_.size() >= 6) {
    if (doesTriangleIntersect(rect, points_[3], points_[4], points_[5],
                              padding)) {
      return true;
    }
    if (points_.size() >= 9) {
      if (doesTriangleIntersect(rect, points_[6], points_[7], points_[8],
                                padding)) {
        return true;
      }
    }
  }
  return false;
}

// ****************************************************************************
void DrawShapeSolidWedge::trimOtherBondVecs() {
  if (otherBondVecs_.size() < 3) {
    return;
  }
  int firstVec = 0, secondVec = 1;
  double largestAng = -361.0;
  for (unsigned int i = 0; i < otherBondVecs_.size() - 1; ++i) {
    for (unsigned int j = i + 1; j < otherBondVecs_.size(); ++j) {
      auto ang = otherBondVecs_[i].angleTo(otherBondVecs_[j]);
      if (ang > largestAng) {
        firstVec = i;
        secondVec = j;
        largestAng = ang;
      }
    }
  }
  std::vector<Point2D> newVecs{otherBondVecs_[firstVec],
                               otherBondVecs_[secondVec]};
  otherBondVecs_ = newVecs;
}

// ****************************************************************************
void DrawShapeSolidWedge::orderOtherBondVecs() {
  if (otherBondVecs_.size() < 2) {
    return;
  }
  // otherBondVecs_[0] needs to be on the same side as points_[1], which
  // implies the larger angle between the 2 vectors.
  auto side1 = (points_[0] - points_[1]);
  if (side1.angleTo(otherBondVecs_[0]) < side1.angleTo(otherBondVecs_[1])) {
    std::swap(otherBondVecs_[0], otherBondVecs_[1]);
  }
}

// ****************************************************************************
DrawShapeDashedWedge::DrawShapeDashedWedge(const std::vector<Point2D> points,
                                           const DrawColour &col1,
                                           const DrawColour &col2,
                                           bool oneLessDash, double lineWidth,
                                           int atom1, int atom2, int bond)
    : DrawShape(points, lineWidth, false, col1, false, atom1, atom2, bond),
      col2_(col2),
      oneLessDash_(oneLessDash) {
  PRECONDITION(points_.size() == 3, "dashed wedge wrong points");
  at1Cds_ = points[0];
  end1Cds_ = points[1];
  end2Cds_ = points[2];
  buildLines();
}

// ****************************************************************************
void DrawShapeDashedWedge::buildLines() {
  auto midend = (end1Cds_ + end2Cds_) * 0.5;
  points_.clear();
  lineColours_.clear();
  auto e1 = at1Cds_.directionVector(end1Cds_);
  auto e2 = at1Cds_.directionVector(end2Cds_);
  // the ACS1996 rules say the dash separation should be 2.5px.  It seems
  // like a good result for all of them.
  // It appears that this means a 2.5px gap between each line, so we need
  // to take the line width into account.  Each line that the gap is
  // between will contribute half a width.
  double dashSep = 2.5 + lineWidth_;
  double centralLen = (at1Cds_ - midend).length();
  unsigned int nDashes = rdcast<unsigned int>(std::round(centralLen / dashSep));
  // There should be at least 3 dashes so we can see which way the wedge
  // is going (Github6041b).
  unsigned int numDashesNeeded = oneLessDash_ ? 4 : 3;
  if (nDashes < numDashesNeeded) {
    nDashes = numDashesNeeded;
  }
  if (!nDashes) {
    points_.push_back(end1Cds_);
    points_.push_back(end2Cds_);
    lineColours_.push_back(lineColour_);
  } else {
    // re-adjust so the last dash is on the end of the wedge.
    dashSep = centralLen / rdcast<double>(nDashes);
    // if doing one less dash, we want a shorter wedge that is just as wide
    // at the end as it would have been.
    if (oneLessDash_) {
      double endlenb2 = (end1Cds_ - end2Cds_).length() / 2.0;
      auto centralLine = at1Cds_.directionVector(midend);
      Point2D centralPerp{-centralLine.y, centralLine.x};
      Point2D newEnd1 = at1Cds_ + centralLine * (centralLen - dashSep) +
                        centralPerp * endlenb2;
      Point2D newEnd2 = at1Cds_ + centralLine * (centralLen - dashSep) -
                        centralPerp * endlenb2;
      e1 = at1Cds_.directionVector(newEnd1);
      e2 = at1Cds_.directionVector(newEnd2);
    }
    // we want the separation down the sides of the triangle, so use
    // similar triangles to scale.
    dashSep *= (end1Cds_ - at1Cds_).length() / centralLen;
    int extra = oneLessDash_ ? 0 : 1;
    for (unsigned int i = 1; i < nDashes + extra; ++i) {
      auto e11 = at1Cds_ + e1 * rdcast<double>(i) * dashSep;
      auto e22 = at1Cds_ + e2 * rdcast<double>(i) * dashSep;
      points_.push_back(e11);
      points_.push_back(e22);
      if (i > nDashes / 2) {
        lineColours_.push_back(col2_);
      } else {
        lineColours_.push_back(lineColour_);
      }
    }
  }
}

// ****************************************************************************
void DrawShapeDashedWedge::myDraw(MolDraw2D &drawer) const {
  drawer.setFillPolys(false);
  drawer.setActiveAtmIdx(atom1_, atom2_);
  drawer.setActiveBndIdx(bond_);
  for (size_t i = 0, j = 0; i < points_.size(); i += 2, ++j) {
    if (drawer.drawOptions().splitBonds) {
      if (i < points_.size() / 2) {
        drawer.setActiveAtmIdx(atom1_);
      } else {
        drawer.setActiveAtmIdx(atom2_);
      }
    }
    drawer.drawLine(points_[i], points_[i + 1], lineColours_[j],
                    lineColours_[j], true);
  }
}

// ****************************************************************************
void DrawShapeDashedWedge::scale(const Point2D &scale_factor) {
  DrawShape::scale(scale_factor);
  at1Cds_.x *= scale_factor.x;
  at1Cds_.y *= scale_factor.y;
  end1Cds_.x *= scale_factor.x;
  end1Cds_.y *= scale_factor.y;
  end2Cds_.x *= scale_factor.x;
  end2Cds_.y *= scale_factor.y;
  buildLines();
}

// ****************************************************************************
void DrawShapeDashedWedge::move(const Point2D &trans) {
  DrawShape::move(trans);
  at1Cds_ += trans;
  end1Cds_ += trans;
  end2Cds_ += trans;
}

// ****************************************************************************
void DrawShapeDashedWedge::findExtremes(double &xmin, double &xmax,
                                        double &ymin, double &ymax) const {
  xmin = std::min({at1Cds_.x, end1Cds_.x, end2Cds_.x, xmin});
  xmax = std::max({at1Cds_.x, end1Cds_.x, end2Cds_.x, xmax});
  ymin = std::min({at1Cds_.y, end1Cds_.y, end2Cds_.y, ymin});
  ymax = std::max({at1Cds_.y, end1Cds_.y, end2Cds_.y, ymax});
}

// ****************************************************************************
bool DrawShapeDashedWedge::doesRectClash(const StringRect &rect,
                                         double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  return doesTriangleIntersect(rect, at1Cds_, end1Cds_, end2Cds_, padding);
}

// ****************************************************************************
DrawShapeWavyLine::DrawShapeWavyLine(const std::vector<Point2D> points,
                                     double lineWidth, bool scaleLineWidth,
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
  // use a negative offset because of inverted y coords to make it look the
  // same as it used to.
  int nsegs = int(
      std::round((points_[0] - points_[1]).length() / (offset_ * 2.0 / 3.0)));
  drawer.drawWavyLine(points_[0], points_[1], lineColour_, col2_, nsegs,
                      -offset_ / 2.0, true);
}

// ****************************************************************************
void DrawShapeWavyLine::scale(const Point2D &scaleFactor) {
  DrawShape::scale(scaleFactor);
  offset_ *= scaleFactor.x;
}

// ****************************************************************************
bool DrawShapeWavyLine::doesRectClash(const StringRect &rect,
                                      double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  padding += offset_;
  return doesLineIntersect(rect, points_[0], points_[1], padding);
}

// ****************************************************************************
DrawShapeArc::DrawShapeArc(const std::vector<Point2D> points, double ang1,
                           double ang2, double lineWidth, bool scaleLineWidth,
                           const DrawColour &col1, bool fill, int atom1)
    : DrawShape(points, lineWidth, scaleLineWidth, col1, fill, atom1),
      ang1_(ang1),
      ang2_(ang2) {
  PRECONDITION(points_.size() == 2, "arc wrong points");
  RANGE_CHECK(0., ang1_, 360.);
  RANGE_CHECK(0., ang2_, 360.);
}

// ****************************************************************************
void DrawShapeArc::myDraw(MolDraw2D &drawer) const {
  if (fill_) {
    drawer.setLineWidth(1);
    drawer.drawOptions().scaleBondWidth = false;
  }
  double start_ang = ang1_ > ang2_ ? ang1_ - 360.0 : ang1_;
  drawer.drawArc(points_[0], points_[1].x, points_[1].y, start_ang, ang2_,
                 true);
}

// ****************************************************************************
void DrawShapeArc::findExtremes(double &xmin, double &xmax, double &ymin,
                                double &ymax) const {
  xmin = std::min(xmin, points_[0].x - points_[1].x);
  xmax = std::max(xmax, points_[0].x + points_[1].x);
  ymin = std::min(ymin, points_[0].y - points_[1].y);
  ymax = std::max(ymax, points_[0].y + points_[1].y);
}

// ****************************************************************************
void DrawShapeArc::move(const Point2D &trans) { points_[0] += trans; }

// ****************************************************************************
bool DrawShapeArc::doesRectClash(const StringRect &rect, double padding) const {
  padding = scaleLineWidth_ ? padding * lineWidth_ : padding;
  Point2D tl, tr, br, bl;
  rect.calcCorners(tl, tr, br, bl, padding);
  if (doesLineIntersectArc(points_[0], points_[1].x, points_[1].y, ang1_, ang2_,
                           padding, tl, tr)) {
    return true;
  }
  if (doesLineIntersectArc(points_[0], points_[1].x, points_[1].y, ang1_, ang2_,
                           padding, tr, br)) {
    return true;
  }
  if (doesLineIntersectArc(points_[0], points_[1].x, points_[1].y, ang1_, ang2_,
                           padding, br, bl)) {
    return true;
  }
  if (doesLineIntersectArc(points_[0], points_[1].x, points_[1].y, ang1_, ang2_,
                           padding, bl, tl)) {
    return true;
  }
  return false;
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit
