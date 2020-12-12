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
#include <GraphMol/Conformer.h>
#include <GraphMol/SubstanceGroup.h>

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

void addStereoAnnotation(const ROMol &mol, bool includeRelativeCIP) {
  const auto &sgs = mol.getStereoGroups();
  std::vector<unsigned int> doneAts(mol.getNumAtoms(), 0);
  unsigned int grpid = 1;
  for (const auto &sg : sgs) {
    for (const auto atom : sg.getAtoms()) {
      if (doneAts[atom->getIdx()]) {
        BOOST_LOG(rdWarningLog) << "Warning: atom " << atom->getIdx()
                                << " is in more than one stereogroup. Only the "
                                   "label from the first group will be used."
                                << std::endl;
        continue;
      }
      std::string lab;
      std::string cip;
      if (includeRelativeCIP ||
          sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
        atom->getPropIfPresent(common_properties::_CIPCode, cip);
      }
      switch (sg.getGroupType()) {
        case StereoGroupType::STEREO_ABSOLUTE:
          lab = "abs";
          break;
        case StereoGroupType::STEREO_OR:
          lab = (boost::format("or%d") % grpid).str();
          break;
        case StereoGroupType::STEREO_AND:
          lab = (boost::format("and%d") % grpid).str();
          break;
        default:
          break;
      }
      if (!lab.empty()) {
        doneAts[atom->getIdx()] = 1;
        if (!cip.empty()) {
          lab += " (" + cip + ")";
        }
        atom->setProp(common_properties::atomNote, lab);
      }
    }
    if (sg.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
      ++grpid;
    }
  }
  for (auto atom : mol.atoms()) {
    std::string cip;
    if (!doneAts[atom->getIdx()] &&
        atom->getPropIfPresent(common_properties::_CIPCode, cip)) {
      std::string lab = "(" + cip + ")";
      atom->setProp(common_properties::atomNote, lab);
    }
  }
  for (auto bond : mol.bonds()) {
    std::string cip;
    if (!bond->getPropIfPresent(common_properties::_CIPCode, cip)) {
      if (bond->getStereo() == Bond::STEREOE) {
        cip = "E";
      } else if (bond->getStereo() == Bond::STEREOZ) {
        cip = "Z";
      }
    }
    if (!cip.empty()) {
      std::string lab = "(" + cip + ")";
      bond->setProp(common_properties::bondNote, lab);
    }
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

namespace {
void drawArrow(MolDraw2D &drawer, const MolDrawShape &shape) {
  PRECONDITION(shape.shapeType == MolDrawShapeType::Arrow, "bad shape type");
  PRECONDITION(shape.points.size() == 4, "bad points size");
  drawer.setColour(shape.lineColour);
  drawer.setLineWidth(shape.lineWidth);
  drawer.drawLine(shape.points[0], shape.points[1]);
  if (!shape.fill) {
    drawer.drawLine(shape.points[1], shape.points[2]);
    drawer.drawLine(shape.points[1], shape.points[3]);
  } else {
    drawer.setFillPolys(true);
    std::vector<Point2D> head(shape.points.begin() + 1, shape.points.end());
    drawer.drawPolygon(head);
  }
}
void drawPolyline(MolDraw2D &drawer, const MolDrawShape &shape) {
  PRECONDITION(shape.shapeType == MolDrawShapeType::Polyline, "bad shape type");
  PRECONDITION(shape.points.size() > 1, "not enough points");
  drawer.setColour(shape.lineColour);
  auto lw = shape.lineWidth;
  if (shape.scaleLineWidth) {
    lw *= drawer.scale() * 0.02;
    if (lw < 0.0) {
      lw = 0.0;
    }
  }
  drawer.setLineWidth(lw);
  if (shape.points.size() > 2 && shape.fill) {
    drawer.setFillPolys(true);
  } else {
    drawer.setFillPolys(false);
  }
  if (drawer.drawOptions().comicMode) {
    auto drawPoints =
        handdrawnLine(shape.points[0], shape.points[1], drawer.scale());
    for (unsigned int i = 2; i < shape.points.size(); ++i) {
      auto lpts = MolDraw2D_detail::handdrawnLine(
          shape.points[i - 1], shape.points[i], drawer.scale());
      std::move(lpts.begin(), lpts.end(), std::back_inserter(drawPoints));
    }
    drawer.drawPolygon(drawPoints);
  } else {
    if (shape.points.size() > 2) {
      drawer.drawPolygon(shape.points);
    } else {
      drawer.drawLine(shape.points[0], shape.points[1]);
    }
  }
}
void drawEllipse(MolDraw2D &drawer, const MolDrawShape &shape) {
  PRECONDITION(shape.shapeType == MolDrawShapeType::Ellipse, "bad shape type");
  PRECONDITION(shape.points.size() == 2, "wrong points");
  drawer.setColour(shape.lineColour);
  drawer.setLineWidth(shape.lineWidth);
  if (shape.fill) {
    drawer.setFillPolys(true);
  } else {
    drawer.setFillPolys(false);
  }
  drawer.drawEllipse(shape.points[0], shape.points[1]);
}
}  // namespace
void drawShapes(MolDraw2D &drawer, const std::vector<MolDrawShape> &shapes) {
  const auto ocolour = drawer.colour();
  const auto olw = drawer.lineWidth();
  const auto ofill = drawer.fillPolys();
  for (const auto &shape : shapes) {
    switch (shape.shapeType) {
      case MolDrawShapeType::Polyline:
        drawPolyline(drawer, shape);
        break;
      case MolDrawShapeType::Arrow:
        drawArrow(drawer, shape);
        break;
      case MolDrawShapeType::Ellipse:
        drawEllipse(drawer, shape);
        break;
      default:
        ASSERT_INVARIANT(false, "unrecognized shape type");
    }
  }
  drawer.setColour(ocolour);
  drawer.setLineWidth(olw);
  drawer.setFillPolys(ofill);
};

// there are a several empirically determined constants here.
std::vector<Point2D> handdrawnLine(Point2D cds1, Point2D cds2, double scale,
                                   bool shiftBegin, bool shiftEnd,
                                   unsigned nSteps, double deviation,
                                   double endShift) {
  // std::cerr << "   " << scale << " " << endShift / scale << endl;
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
}  // namespace MolDraw2D_detail
}  // namespace RDKit
