//
//  Copyright (C) 2015-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/StringRect.h>
#include <GraphMol/Chirality.h>

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ****************************************************************************

namespace RDKit {
namespace MolDraw2D_detail {
// implementation from $RDBASE/rdkit/sping/pid.py
void arcPoints(const Point2D &cds1, const Point2D &cds2,
               std::vector<Point2D> &res, float startAng, float extent) {
  // Note: this implementation is simple and not particularly efficient.
  float xScale = (cds2.x - cds1.x) / 2.0;
  float yScale = (cds2.y - cds1.y) / 2.0;
  if (xScale < 0) {
    xScale *= -1;
  }
  if (yScale < 0) {
    yScale *= -1;
  }

  float x = std::min(cds1.x, cds2.x) + xScale;
  float y = std::min(cds1.y, cds2.y) + yScale;

  int steps = std::max(static_cast<int>(extent * 2), 5);
  float step = M_PI * extent / (180 * steps);
  float angle = M_PI * startAng / 180;
  for (int i = 0; i <= steps; ++i) {
    Point2D point(x + xScale * cos(angle), y - yScale * sin(angle));
    res.emplace_back(point);
    angle += step;
  }
}

namespace {
// note, this is approximate since we're just using it for drawing
bool lineSegmentsIntersect(const Point2D &s1, const Point2D &s2,
                           const Point2D &s3, const Point2D &s4) {
  auto d1x = (s1.x - s2.x);
  auto d1y = (s1.y - s2.y);
  auto d2x = (s3.x - s4.x);
  auto d2y = (s3.y - s4.y);

  if (fabs(d1x) < 1e-4) {
    // fudge factor, since this isn't super critical
    d1x = 1e-4;
  }
  if (fabs(d2x) < 1e-4) {
    // fudge factor, since this isn't super critical
    d2x = 1e-4;
  }

  auto m1 = d1y / d1x;
  auto m2 = d2y / d2x;
  if (m1 == m2 || m1 == -m2) {
    // parallel
    return false;
  }
  auto b1 = (s1.x * s2.y - s2.x * s1.y) / d1x;
  auto b2 = (s3.x * s4.y - s4.x * s3.y) / d2x;

  auto intersectX = (b2 - b1) / (m1 - m2);
  return ((intersectX < s1.x) ^ (intersectX < s2.x)) &&
         ((intersectX < s3.x) ^ (intersectX < s4.x));
}
}  // namespace

std::vector<Point2D> getBracketPoints(
    const Point2D &p1, const Point2D &p2, const Point2D &refPt,
    const std::vector<std::pair<Point2D, Point2D>> &bondSegments,
    double bracketFrac) {
  std::vector<Point2D> res;
  auto v = p2 - p1;
  Point2D bracketDir{v.y, -v.x};
  bracketDir *= bracketFrac;

  // we'll default to use the refPt
  auto refVect = p2 - refPt;
  // but check if we intersect any of the bonds:
  for (const auto &seg : bondSegments) {
    if (lineSegmentsIntersect(p1, p2, seg.first, seg.second)) {
      refVect = p2 - seg.first;
    }
  }
  if (bracketDir.dotProduct(refVect) > 0) {
    bracketDir *= -1;
  }
  auto p0 = p1 + bracketDir;
  auto p3 = p2 + bracketDir;
  return {p0, p1, p2, p3};
}

// there are a several empirically determined constants here.
std::vector<Point2D> handdrawnLine(Point2D cds1, Point2D cds2, double scale,
                                   bool shiftBegin, bool shiftEnd,
                                   unsigned nSteps, double deviation,
                                   double endShift) {
  // std::cout << "   " << scale << " " << endShift / scale << std::endl;
  while (endShift / scale > 0.02) {
    endShift *= 0.75;
  }
  if (shiftBegin) {
    cds1.x += (std::rand() % 10 >= 5 ? endShift : -endShift) / scale;
    cds1.y += (std::rand() % 10 >= 5 ? endShift : -endShift) / scale;
  }
  if (shiftEnd) {
    cds2.x += (std::rand() % 10 >= 5 ? endShift : -endShift) / scale;
    cds2.y += (std::rand() % 10 >= 5 ? endShift : -endShift) / scale;
  }

  Point2D step = (cds2 - cds1) / nSteps;
  // make sure we aren't adding loads of wiggles to short lines
  while (step.length() < 0.2 && nSteps > 2) {
    --nSteps;
    step = (cds2 - cds1) / nSteps;
  }
  // make sure the wiggles aren't too big
  while (deviation / step.length() > 0.15 || deviation * scale > 0.70) {
    deviation *= 0.75;
  }
  Point2D perp{step.y, -step.x};
  perp.normalize();
  std::vector<Point2D> pts;
  pts.push_back(cds1);
  for (unsigned int i = 1; i < nSteps; ++i) {
    auto tgt = cds1 + step * i;
    tgt += perp * deviation * (std::rand() % 20 - 10) / 10.0;
    pts.push_back(tgt);
  }
  pts.push_back(cds2);
  return pts;
}

// ****************************************************************************
bool doesLineIntersect(const StringRect &rect, const Point2D &end1,
                       const Point2D &end2, double padding) {
  Point2D tl, tr, bl, br;
  rect.calcCorners(tl, tr, br, bl, padding);
  if (doLinesIntersect(end2, end1, tl, tr, nullptr)) {
    return true;
  }
  if (doLinesIntersect(end2, end1, tr, br, nullptr)) {
    return true;
  }
  if (doLinesIntersect(end2, end1, br, bl, nullptr)) {
    return true;
  }
  if (doLinesIntersect(end2, end1, bl, tl, nullptr)) {
    return true;
  }
  return false;
}

// ****************************************************************************
bool doesTriangleIntersect(const StringRect &rect, const Point2D &pt1,
                           const Point2D &pt2, const Point2D &pt3,
                           double padding) {
  // the quick test is for any of the triangle points inside the rectangle.
  if (rect.isPointInside(pt1, padding) || rect.isPointInside(pt2, padding) ||
      rect.isPointInside(pt3, padding)) {
    return true;
  }
  // But if the rectangle is inside the triangle, that's not enough of a test.
  Point2D tl, tr, br, bl;
  rect.calcCorners(tl, tr, br, bl, padding);
  if (isPointInTriangle(tl, pt1, pt2, pt3) ||
      isPointInTriangle(tr, pt1, pt2, pt3) ||
      isPointInTriangle(br, pt1, pt2, pt3) ||
      isPointInTriangle(bl, pt1, pt2, pt3)) {
    return true;
  }
  // and finally all the points in the rectangle can be outside the triangle,
  // but the sides can cross it.  And vice versa.  So see if any of the sides
  // intersect.
  if (doLinesIntersect(tl, tr, pt1, pt2, nullptr) ||
      doLinesIntersect(tl, tr, pt2, pt3, nullptr) ||
      doLinesIntersect(tl, tr, pt3, pt1, nullptr) ||
      doLinesIntersect(tr, br, pt1, pt2, nullptr) ||
      doLinesIntersect(tr, br, pt2, pt3, nullptr) ||
      doLinesIntersect(tr, br, pt3, pt1, nullptr) ||
      doLinesIntersect(br, bl, pt1, pt2, nullptr) ||
      doLinesIntersect(br, bl, pt2, pt3, nullptr) ||
      doLinesIntersect(br, bl, pt3, pt1, nullptr) ||
      doLinesIntersect(bl, tl, pt1, pt2, nullptr) ||
      doLinesIntersect(bl, tl, pt2, pt3, nullptr) ||
      doLinesIntersect(bl, tl, pt3, pt1, nullptr)) {
    return true;
  }
  return false;
}

// ****************************************************************************
bool doesLineIntersectEllipse(const Point2D &centre, double xradius,
                              double yradius, double padding,
                              const Point2D &end1, const Point2D &end2) {
  // using
  // https://math.stackexchange.com/questions/76457/check-if-a-point-is-within-an-ellipse
  // to see if either end1 or end2 are inside the ellipse.
  double xr2 = (xradius + padding) * (xradius + padding);
  double yr2 = (yradius + padding) * (yradius + padding);
  double xdisc = (end1.x - centre.x) * (end1.x - centre.x) / xr2;
  double ydisc = (end1.y - centre.y) * (end1.y - centre.y) / yr2;
  if (xdisc + ydisc <= 1.0) {
    return true;
  }
  xdisc = (end2.x - centre.x) * (end2.x - centre.x) / xr2;
  ydisc = (end2.y - centre.y) * (end2.y - centre.y) / yr2;
  return xdisc + ydisc <= 1.0;
}

// ****************************************************************************
bool doesLineIntersectArc(const Point2D &centre, double xradius, double yradius,
                          double start_ang, double stop_ang, double padding,
                          const Point2D &end1, const Point2D &end2) {
  double xr = xradius + padding;
  double yr = yradius + padding;
  double xr2 = xr * xr;
  double yr2 = yr * yr;

  auto pointInArc = [&](const Point2D &p) -> bool {
    double xdisc = p.x * p.x / xr2;
    double ydisc = p.y * p.y / yr2;
    // start_ang can be more than stop_ang, if, for example, the arc goes from
    // 315 to 45.
    if (xdisc + ydisc <= 1.0) {
      // end1 is inside the whole ellipse.  See if the angle it makes with the
      // x axis lies between start_and and stop_ang.
      double th = atan2(p.x, p.y) * 180.0 / M_PI;
      if (th < 0.0) {
        th += 360.0;
      }
      if (start_ang < stop_ang) {
        // it's pretty straightforward
        if (th >= start_ang && th <= stop_ang) {
          return true;
        }
      } else {
        // the arc crosses 0.
        if (th >= start_ang && th <= 360.0) {
          return true;
        }
        if (th >= 0.0 && th <= stop_ang) {
          return true;
        }
      }
    }
    return false;
  };

  Point2D p1 = end1 - centre;
  if (pointInArc(p1)) {
    return true;
  }
  Point2D p2 = end2 - centre;
  return pointInArc(p2);
}

// ****************************************************************************
bool doLinesIntersect(const Point2D &l1s, const Point2D &l1f,
                      const Point2D &l2s, const Point2D &l2f, Point2D *ip) {
  // using spell from answer 2 of
  // https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
  double s1_x = l1f.x - l1s.x;
  double s1_y = l1f.y - l1s.y;
  double s2_x = l2f.x - l2s.x;
  double s2_y = l2f.y - l2s.y;

  double d = (-s2_x * s1_y + s1_x * s2_y);
  if (d == 0.0) {
    // parallel lines.
    return false;
  }
  double s, t;
  s = (-s1_y * (l1s.x - l2s.x) + s1_x * (l1s.y - l2s.y)) / d;
  t = (s2_x * (l1s.y - l2s.y) - s2_y * (l1s.x - l2s.x)) / d;

  if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
    if (ip) {
      ip->x = l1s.x + t * s1_x;
      ip->y = l1s.y + t * s1_y;
    }
    return true;
  }

  return false;
}

// ****************************************************************************
bool isPointInTriangle(const Point2D &pt, const Point2D &t1, const Point2D &t2,
                       const Point2D &t3) {
  double d = ((t2.y - t3.y) * (t1.x - t3.x) + (t3.x - t2.x) * (t1.y - t3.y));
  double a =
      ((t2.y - t3.y) * (pt.x - t3.x) + (t3.x - t2.x) * (pt.y - t3.y)) / d;
  double b =
      ((t3.y - t1.y) * (pt.x - t3.x) + (t1.x - t3.x) * (pt.y - t3.y)) / d;
  double c = 1 - a - b;
  return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1;
}

std::vector<std::tuple<Point2D, Point2D, Point2D, Point2D>> getWavyLineSegments(
    const Point2D &p1, const Point2D &p2, unsigned int nSegments,
    double vertOffset) {
  std::vector<std::tuple<Point2D, Point2D, Point2D, Point2D>> res;

  PRECONDITION(nSegments > 1, "too few segments");

  if (nSegments % 2) {
    ++nSegments;  // we're going to assume an even number of segments
  }

  Point2D delta = (p2 - p1);
  Point2D perp(delta.y, -delta.x);
  perp.normalize();
  perp *= vertOffset;
  delta /= nSegments;

  for (unsigned int i = 0; i < nSegments; ++i) {
    Point2D startpt = p1 + delta * i;
    Point2D segpt = startpt + delta;
    Point2D cpt1 = startpt + perp * (i % 2 ? -1 : 1);
    Point2D cpt2 = segpt + perp * (i % 2 ? -1 : 1);
    res.emplace_back(startpt, cpt1, cpt2, segpt);
  }
  return res;
}

RDKIT_MOLDRAW2D_EXPORT void calcArrowHead(Point2D &arrowEnd, Point2D &arrow1,
                                          Point2D &arrow2,
                                          const Point2D &arrowBegin,
                                          bool asPolygon, double frac,
                                          double angle) {
  auto delta = arrowBegin - arrowEnd;
  double cos_angle = std::cos(angle), sin_angle = std::sin(angle);
  // to have the arrowhead a consistent fraction of the line length, we need
  // the hypotenuse
  frac /= cos_angle;
  if (asPolygon) {
    // allow for the mitring, using an empirically derived guess.
    arrowEnd += delta * 0.1;
  }
  arrow1 = arrowEnd;
  arrow1.x += frac * (delta.x * cos_angle + delta.y * sin_angle);
  arrow1.y += frac * (delta.y * cos_angle - delta.x * sin_angle);

  arrow2 = arrowEnd;
  arrow2.x += frac * (delta.x * cos_angle - delta.y * sin_angle);
  arrow2.y += frac * (delta.y * cos_angle + delta.x * sin_angle);
}

// ****************************************************************************
void adjustLineEndForEllipse(const Point2D &centre, double xradius,
                             double yradius, Point2D p1, Point2D &p2) {
  // move everything so the ellipse is centred on the origin.
  p1 -= centre;
  p2 -= centre;
  double a2 = xradius * xradius;
  double b2 = yradius * yradius;
  double A =
      (p2.x - p1.x) * (p2.x - p1.x) / a2 + (p2.y - p1.y) * (p2.y - p1.y) / b2;
  double B = 2.0 * p1.x * (p2.x - p1.x) / a2 + 2.0 * p1.y * (p2.y - p1.y) / b2;
  double C = p1.x * p1.x / a2 + p1.y * p1.y / b2 - 1.0;

  auto t_to_point = [&](double t) -> Point2D {
    Point2D ret_val;
    ret_val.x = p1.x + (p2.x - p1.x) * t + centre.x;
    ret_val.y = p1.y + (p2.y - p1.y) * t + centre.y;
    return ret_val;
  };

  double disc = B * B - 4.0 * A * C;
  if (disc < 0.0) {
    // no solutions, leave things as they are.  Bit crap, though.
    p2 += centre;
    return;
  } else if (fabs(disc) < 1.0e-6) {
    // 1 solution
    double t = -B / (2.0 * A);
    p2 = t_to_point(t);
  } else {
    // 2 solutions - take the one nearest p1.
    double disc_rt = sqrt(disc);
    double t1 = (-B + disc_rt) / (2.0 * A);
    double t2 = (-B - disc_rt) / (2.0 * A);
    double t;
    // prefer the t between 0 and 1, as that must be between the original
    // points.  If both are, prefer the lower, as that will be nearest p1,
    // so on the bit of the ellipse the line comes to first.
    bool t1_ok = (t1 >= 0.0 && t1 <= 1.0);
    bool t2_ok = (t2 >= 0.0 && t2 <= 1.0);
    if (t1_ok && !t2_ok) {
      t = t1;
    } else if (t2_ok && !t1_ok) {
      t = t2;
    } else if (t1_ok && t2_ok) {
      t = std::min(t1, t2);
    } else {
      // the intersections are both outside the line between p1 and p2
      // so don't do anything.
      p2 += centre;
      return;
    }
    p2 = t_to_point(t);
  }
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit
