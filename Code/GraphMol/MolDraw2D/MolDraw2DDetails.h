//
//  Copyright (C) 2015-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDKITMOLDRAW2DDETAILS_H
#define RDKITMOLDRAW2DDETAILS_H

#include <vector>

#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>

// ****************************************************************************
using RDGeom::Point2D;

namespace RDKit {
namespace MolDraw2D_detail {
struct StringRect;

// data taken from the helvetica font info in
// $RDBASE/rdkit/sping/PDF/pdfmetrics.py
const int char_widths[] = {
    0,   0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   0,    0,   0,   278, 278, 355,  556,  556,  889, 667,  222, 333,  333,
    389, 584,  278, 333, 278, 278, 556,  556,  556,  556, 556,  556, 556,  556,
    556, 556,  278, 278, 584, 584, 584,  556,  1015, 667, 667,  722, 722,  667,
    611, 778,  722, 278, 500, 667, 556,  833,  722,  778, 667,  778, 722,  667,
    611, 722,  667, 944, 667, 667, 611,  278,  278,  278, 469,  556, 222,  556,
    556, 500,  556, 556, 278, 556, 556,  222,  222,  500, 222,  833, 556,  556,
    556, 556,  333, 500, 278, 556, 500,  722,  500,  500, 500,  334, 260,  334,
    584, 0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   0,    0,   0,   0,   0,   0,    333,  556,  556, 167,  556, 556,  556,
    556, 191,  333, 556, 333, 333, 500,  500,  0,    556, 556,  556, 278,  0,
    537, 350,  222, 333, 333, 556, 1000, 1000, 0,    611, 0,    333, 333,  333,
    333, 333,  333, 333, 333, 0,   333,  333,  0,    333, 333,  333, 1000, 0,
    0,   0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   1000, 0,   370, 0,   0,   0,    0,    556,  778, 1000, 365, 0,    0,
    0,   0,    0,   889, 0,   0,   0,    278,  0,    0,   222,  611, 944,  611,
    0,   0,    834};

// angles in degrees.
RDKIT_MOLDRAW2D_EXPORT void arcPoints(const Point2D &cds1, const Point2D &cds2,
                                      std::vector<Point2D> &res,
                                      float startAng = 0, float extent = 360);

//! add R/S, relative stereo, and E/Z annotations to atoms and bonds
RDKIT_MOLDRAW2D_EXPORT void addStereoAnnotation(
    const ROMol &mol, bool includeRelativeCIP = false);

//! add annotations with atom indices.
RDKIT_MOLDRAW2D_EXPORT inline void addAtomIndices(const ROMol &mol) {
  // we don't need this in the global set of tags since it will only be used
  // here
  if (mol.hasProp("_atomIndicesAdded")) {
    return;
  }
  bool computed = true;
  mol.setProp("_atomIndicesAdded", 1, computed);
  for (auto atom : mol.atoms()) {
    auto lab = std::to_string(atom->getIdx());
    if (atom->hasProp(common_properties::atomNote)) {
      lab += "," + atom->getProp<std::string>(common_properties::atomNote);
    }
    atom->setProp(common_properties::atomNote, lab);
  }
};

//! add annotations with bond indices.
RDKIT_MOLDRAW2D_EXPORT inline void addBondIndices(const ROMol &mol) {
  // we don't need this in the global set of tags since it will only be used
  // here
  if (mol.hasProp("_bondIndicesAdded")) {
    return;
  }
  bool computed = true;
  mol.setProp("_bondIndicesAdded", 1, computed);
  for (auto bond : mol.bonds()) {
    auto lab = std::to_string(bond->getIdx());
    if (bond->hasProp(common_properties::bondNote)) {
      lab += "," + bond->getProp<std::string>(common_properties::bondNote);
    }
    bond->setProp(common_properties::bondNote, lab);
  }
};

RDKIT_MOLDRAW2D_EXPORT std::vector<Point2D> getBracketPoints(
    const Point2D &p1, const Point2D &p2, const Point2D &refPt,
    const std::vector<std::pair<Point2D, Point2D>> &bondSegments,
    double bracketFrac = 0.1);
// there are a several empirically determined constants here.
RDKIT_MOLDRAW2D_EXPORT std::vector<Point2D> handdrawnLine(
    Point2D cds1, Point2D cds2, double scale, bool shiftBegin = false,
    bool shiftEnd = false, unsigned nSteps = 4, double deviation = 0.03,
    double endShift = 0.5);

inline std::string formatDouble(double val) {
  return boost::str(boost::format("%.1f") % val);
}

RDKIT_MOLDRAW2D_EXPORT bool doesLineIntersect(const StringRect &rect,
                                              const Point2D &end1,
                                              const Point2D &end2,
                                              double padding);
// returns true if any corner of triangle is inside the rectangle.
RDKIT_MOLDRAW2D_EXPORT bool doesTriangleIntersect(const StringRect &rect,
                                                  const Point2D &pt1,
                                                  const Point2D &pt2,
                                                  const Point2D &pt3,
                                                  double padding);
RDKIT_MOLDRAW2D_EXPORT bool doesLineIntersectEllipse(
    const Point2D &centre, double xradius, double yradius, double padding,
    const Point2D &end1, const Point2D &end2);
// angles expected in degrees, between 0 and 360.
RDKIT_MOLDRAW2D_EXPORT bool doesLineIntersectArc(
    const Point2D &centre, double xradius, double yradius, double start_ang,
    double stop_ang, double padding, const Point2D &end1, const Point2D &end2);
RDKIT_MOLDRAW2D_EXPORT bool doLinesIntersect(const Point2D &l1s,
                                             const Point2D &l1f,
                                             const Point2D &l2s,
                                             const Point2D &l2f, Point2D *ip);
// This uses the barycentric coordinate system method from
// http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
// where it notes and provides a solution for instabilities when the point
// in exactly on one of the edges of the triangle.  That refinement is not
// implemented because it seems a bit of overkill for most uses.  It is an
// issue when, for example, two triangles share an edge and the point is on that
// edge, when it might give the disappointing result that the point is in
// neither triangle.
RDKIT_MOLDRAW2D_EXPORT bool isPointInTriangle(const Point2D &pt,
                                              const Point2D &t1,
                                              const Point2D &t2,
                                              const Point2D &t3);

// returns a vector of p1,c1,c2,p2 tuples for bezier curves
RDKIT_MOLDRAW2D_EXPORT
std::vector<std::tuple<Point2D, Point2D, Point2D, Point2D>> getWavyLineSegments(
    const Point2D &p1, const Point2D &p2, unsigned int nSegments,
    double vertOffset);

// calculate the points making up the arrowhead of a DrawShapeArrow, allowing
// for the fact that in polygon mode the point can extend over the end
// of the point, because of the mitring.
RDKIT_MOLDRAW2D_EXPORT void calcArrowHead(Point2D &arrowEnd, Point2D &arrow1,
                                          Point2D &arrow2, const Point2D &arrowBegin,
                                          bool asPolygon, double frac, double angle);
}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif
