//
//  Copyright (C) 2023 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This is based on a suggestion and code from Christian Feldmann.
// It was discussion #4607.  His Python implementation (which I haven't
// followed to any great extent) is at
// https://github.com/c-feldmann/lassohighlight

#include <vector>

#include <GraphMol/RWMol.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/DrawMolMCHLasso.h>

namespace RDKit {
namespace MolDraw2D_detail {

// an empirically derived lineWidth.
int LINE_WIDTH = 3;
bool SCALE_LINE_WIDTH = true;

// ****************************************************************************
DrawMolMCHLasso::DrawMolMCHLasso(
    const ROMol &mol, const std::string &legend, int width, int height,
    MolDrawOptions &drawOptions, DrawText &textDrawer,
    const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
    const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
    const std::map<int, double> &highlight_radii,
    const std::map<int, int> &highlight_linewidth_multipliers, int confId)
    : DrawMolMCH(mol, legend, width, height, drawOptions, textDrawer,
                 highlight_atom_map, highlight_bond_map, highlight_radii,
                 highlight_linewidth_multipliers, confId) {}

// ****************************************************************************
void DrawMolMCHLasso::extractHighlights(double /* scale  */) {
  extractMCHighlights();
}

// ****************************************************************************
void DrawMolMCHLasso::extractMCHighlights() {
  std::vector<DrawColour> colours;
  std::vector<std::vector<int>> colourAtoms;
  std::vector<std::vector<int>> colourLists;
  extractAtomColourLists(colours, colourAtoms, colourLists);
  for (const auto &colourList : colourLists) {
    for (size_t i = 0U; i < colourList.size(); ++i) {
      drawLasso(i, colours[colourList[i]], colourAtoms[colourList[i]]);
    }
  }
}

// ****************************************************************************
// Get the atoms to lasso in the given colours.  Split the lists up into
// different sets that overlap.  That way, we won't have a single lasso
// that is larger than ones in a different set that doesn't overlap with it.
void DrawMolMCHLasso::extractAtomColourLists(
    std::vector<DrawColour> &colours,
    std::vector<std::vector<int>> &colourAtoms,
    std::vector<std::vector<int>> &colourLists) const {
  std::vector<boost::dynamic_bitset<>> inColourAtoms;
  for (const auto &cm : mcHighlightAtomMap_) {
    for (const auto &col : cm.second) {
      auto cpos = std::find(colours.begin(), colours.end(), col);
      if (cpos == colours.end()) {
        colours.push_back(col);
        colourAtoms.push_back(std::vector<int>(1, cm.first));
        inColourAtoms.push_back(
            boost::dynamic_bitset<>(drawMol_->getNumAtoms()));
        inColourAtoms.back().set(cm.first);
      } else {
        auto ln = std::distance(colours.begin(), cpos);
        // it's important that each atom is only in the list once - Github6749
        if (!inColourAtoms[ln][cm.first]) {
          colourAtoms[ln].push_back(cm.first);
          inColourAtoms[ln].set(cm.first);
        }
      }
    }
  }
  for (size_t i = 0U; i < colourAtoms.size(); ++i) {
    colourLists.push_back(std::vector<int>(1, i));
  }
  if (colourLists.size() < 2) {
    return;
  }
  // This is pretty inefficient, but the lists should all be small.  It doesn't
  // seem worth doing anything more sophisticated.
  auto listsIntersect =
      [](const std::vector<int> &v1, const std::vector<int> &v2,
         const std::vector<std::vector<int>> &colourAtoms) -> bool {
    for (auto &mv1 : v1) {
      for (auto ca : colourAtoms[mv1]) {
        for (auto &mv2 : v2) {
          if (std::find(colourAtoms[mv2].begin(), colourAtoms[mv2].end(), ca) !=
              colourAtoms[mv2].end()) {
            return true;
          }
        }
      }
    }
    return false;
  };
  bool didSomething = true;
  while (didSomething && colourLists.size() > 1) {
    didSomething = false;
    for (size_t i = 0U; i < colourLists.size() - 1; ++i) {
      for (size_t j = i + 1; j < colourLists.size(); ++j) {
        if (listsIntersect(colourLists[i], colourLists[j], colourAtoms)) {
          colourLists[i].insert(colourLists[i].end(), colourLists[j].begin(),
                                colourLists[j].end());
          colourLists[j].clear();
          didSomething = true;
          break;
        }
      }
      if (didSomething) {
        break;
      }
    }
    colourLists.erase(std::remove_if(colourLists.begin(), colourLists.end(),
                                     [](const std::vector<int> &v) -> bool {
                                       return v.empty();
                                     }),
                      colourLists.end());
  }
}

// ****************************************************************************
void DrawMolMCHLasso::drawLasso(size_t lassoNum, const RDKit::DrawColour &col,
                                const std::vector<int> &colAtoms) {
  // Extract the arcs and lines for the given atoms in the given colour.
  // lassoNum is the number of the lasso being done, and hence dictates
  // the radii of the arcs and the displacements of the lines.
  if (colAtoms.empty()) {
    return;
  }
  std::vector<std::unique_ptr<DrawShapeArc>> arcs;
  std::vector<std::unique_ptr<DrawShapeSimpleLine>> lines;
  std::vector<std::vector<LinePair>> atomLines(drawMol_->getNumAtoms());
  extractBondLines(lassoNum, col, colAtoms, lines, atomLines);
  extractAtomArcs(atomLines, arcs);
  addNoLineArcs(colAtoms, lassoNum, col, lines, arcs);

  for (auto &it : arcs) {
    highlights_.push_back(std::move(it));
  }
  for (auto &it : lines) {
    highlights_.emplace_back(std::move(it));
  }
}

namespace {
double getLassoWidth(const DrawMolMCH *dm, int atNum, int lassoNum) {
  PRECONDITION(dm, "Needs valid DrawMolMCH")
  double xrad, yrad;
  dm->getAtomRadius(atNum, xrad, yrad);
  // Double the area of the circles for successive lassos.
  const static double rats[] = {1.0, 1.414, 2, 2.828, 4};
  if (lassoNum > 4) {
    // It's going to look horrible, probably, but it's a lot of lassos.
    return xrad * (1 + lassoNum) * 0.75;
  } else {
    return xrad * rats[lassoNum];
  }
}
}  // namespace

namespace {
Point2D arcEnd(const DrawShapeArc &arc, double ang) {
  // angles are in degrees
  ang *= M_PI / 180.0;
  return Point2D{arc.points_[0].x + arc.points_[1].x * cos(ang),
                 arc.points_[0].y + arc.points_[1].x * sin(ang)};
};

// given 2 points pt1, pt2, assumed to be either side of the line between
// points at1 and at2, compute the angle of the line from at1 to pt1 and
// the x-axis and likewise for pt2 and the mid-point of the line between
// at1 and at2.  All angles go anti-clockwise from the x-axis and are in
// degrees ready for plugging straight into a DrawShapeArc.  If they
// had to be swapped so that ang1 is less than ang2, swapped is set true.
void calcSubtendedAngles(const Point2D &pt1, const Point2D &pt2,
                         const Point2D &at1, const Point2D &at2, double &rang1,
                         double &rang2, double &bang, bool &swapped) {
  static const Point2D index{1.0, 0.0};

  auto rad1 = at1.directionVector(pt1);
  auto rad2 = at1.directionVector(pt2);
  auto brad = at1.directionVector(at2);

  double ang1 = 360.0 - rad1.signedAngleTo(index) * 180.0 / M_PI;
  if (ang1 >= 360.0) {
    ang1 -= 360.0;
  }
  double ang2 = 360.0 - rad2.signedAngleTo(index) * 180.0 / M_PI;
  if (ang2 >= 360.0) {
    ang2 -= 360.0;
  }
  bang = 360.0 - brad.signedAngleTo(index) * 180.0 / M_PI;
  if (bang >= 360.0) {
    bang -= 360.0;
  }
  double cross = rad1.x * rad2.y - rad1.y * rad2.x;
  swapped = false;
  if (cross > 0.0) {
    std::swap(ang1, ang2);
    swapped = true;
  }
  double minAng = std::min({ang1, ang2, bang});
  if (ang1 - minAng < bang - minAng && bang - minAng < ang2 - minAng) {
    rang1 = ang1;
    rang2 = ang2;
  } else {
    rang1 = ang2;
    rang2 = ang1;
    swapped = !swapped;
  }
}

}  // namespace

// ****************************************************************************
void DrawMolMCHLasso::addNoLineArcs(
    const std::vector<int> &colAtoms, size_t lassoNum,
    const RDKit::DrawColour &col,
    const std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines,
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs) const {
  boost::dynamic_bitset<> inColAtoms(drawMol_->getNumAtoms());
  for (auto &ca : colAtoms) {
    inColAtoms.set(ca);
  }
  std::vector<boost::dynamic_bitset<>> inLines(
      drawMol_->getNumAtoms(),
      boost::dynamic_bitset<>(drawMol_->getNumAtoms()));
  for (auto &l : lines) {
    inLines[l->atom2_].set(l->atom1_);
    inLines[l->atom1_].set(l->atom2_);
  }
  std::vector<std::vector<unsigned int>> intersects(
      drawMol_->getNumAtoms(), std::vector<unsigned int>());
  for (unsigned int i = 0; i < drawMol_->getNumAtoms(); ++i) {
    if (inColAtoms[i] && inLines[i].none()) {
      intersects[i].push_back(i);
    }
  }

  for (unsigned int i = 0; i < drawMol_->getNumAtoms() - 1; ++i) {
    if (!inColAtoms[i]) {
      continue;
    }
    auto radI = getLassoWidth(this, i, lassoNum);
    for (unsigned int j = i + 1; j < drawMol_->getNumAtoms(); ++j) {
      if (!inColAtoms[j]) {
        continue;
      }
      auto radJ = getLassoWidth(this, j, lassoNum);
      if (!inLines[i][j] && (atCds_[i] - atCds_[j]).length() < radI + radJ) {
        if (intersects[i].empty()) {
          intersects[i].push_back(i);
        }
        intersects[i].push_back(j);
        if (intersects[j].empty()) {
          intersects[j].push_back(j);
        }
        intersects[j].push_back(i);
      }
    }
  }
  std::vector<DrawShapeArc *> newArcs;
  makeIntersectingArcs(intersects, lassoNum, col, arcs, newArcs);
  for (auto &a : newArcs) {
    arcs.emplace_back(a);
  }
}

namespace {
void intersectingCirclePointsSameSize(const Point2D &at1Cds,
                                      const Point2D &at2Cds, double rad,
                                      Point2D &isect1, Point2D &isect2) {
  // if the radii are the same, the intersection points are where
  // the perpendicular at the mid-point between the circles intersects
  // with either circle.
  Point2D mid = (at1Cds + at2Cds) / 2.0;
  auto perp = calcPerpendicular(at1Cds, at2Cds);
  isect1 = mid + perp * 2.0 * rad;
  adjustLineEndForEllipse(at1Cds, rad, rad, mid, isect1);
  isect2 = mid - perp * (mid - isect1).length();
}

void intersectingCirclePointsDiffSize(const Point2D &at1Cds, double rad1,
                                      const Point2D &at2Cds, double rad2,
                                      Point2D &isect1, Point2D &isect2) {
  // Using the same approach as for the same size, but having to do some
  // extra work because the perpendicular isn't on the mid-point of the
  // line between the centres.  Code adapted from
  // http://ambrnet.com/TrigoCalc/Circles2/circle2intersection/CircleCircleIntersection.htm
  // with explanation
  // https://math.stackexchange.com/questions/256100/how-can-i-find-the-points-at-which-two-circles-intersect
  double dist = (at1Cds - at2Cds).length();
  double dist2 = dist * dist;
  // We already know the circles intersect, so no need for further
  // checks.
  // Triangle area via Heron's formula
  double a1 = dist + rad1 + rad2;
  double a2 = dist + rad1 - rad2;
  double a3 = dist - rad1 + rad2;
  double a4 = rad1 + rad2 - dist;
  double area = sqrt(a1 * a2 * a3 * a4) / 4.0;
  double c = (rad1 * rad1 - rad2 * rad2) / (2.0 * dist2);
  // x values
  double val1 = (at1Cds.x + at2Cds.x) / 2.0 + (at2Cds.x - at1Cds.x) * c;
  double val2 = 2.0 * (at1Cds.y - at2Cds.y) * area / dist2;
  isect1.x = val1 + val2;
  isect2.x = val1 - val2;

  // now y
  val1 = (at1Cds.y + at2Cds.y) / 2.0 + (at2Cds.y - at1Cds.y) * c;
  val2 = 2.0 * (at1Cds.x - at2Cds.x) * area / dist2;
  isect1.y = val1 + val2;
  isect2.y = val1 - val2;

  // Finally, need to check that the 2 combinations of y for x we've taken
  // are the correct 2 by making sure isect1 is on the circle
  if (fabs((isect1 - at1Cds).length() - rad1 * rad1) > 1.0e-8) {
    std::swap(isect1.y, isect2.y);
  }
}

// Find the points of intersection of 2 circles.
void intersectingCirclePoints(const Point2D &at1Cds, double rad1,
                              const Point2D &at2Cds, double rad2,
                              Point2D &isect1, Point2D &isect2) {
  if (fabs(rad1 - rad2) < 1.0e-8) {
    intersectingCirclePointsSameSize(at1Cds, at2Cds, rad1, isect1, isect2);
  } else {
    intersectingCirclePointsDiffSize(at1Cds, rad1, at2Cds, rad2, isect1,
                                     isect2);
  }
}

// get the angles for the intersections of the 2 circles given,
// relative to the first one.  They are the angles from the x axis
// going anti-clockwise. Also returns the angle for where the line
// between the centres intersects with the circle on at1Cds (bang).
void calcIntersectingArcAngles(const Point2D &at1Cds, double rad1,
                               const Point2D &at2Cds, double rad2, double &ang1,
                               double &ang2, double &bang) {
  Point2D isect1, isect2;
  intersectingCirclePoints(at1Cds, rad1, at2Cds, rad2, isect1, isect2);
  bool swapped;
  calcSubtendedAngles(isect1, isect2, at1Cds, at2Cds, ang1, ang2, bang,
                      swapped);
}

// See if the arc defined by the centre, radius etc. has either of its ends
// inside a circle of one of the otherAtoms
bool areArcEndsInOtherCircle(const Point2D &centre, double radius, double ang1,
                             double ang2, size_t atNum, const DrawMolMCH *dm,
                             const std::vector<unsigned int> &otherAtoms,
                             int lassoNum) {
  ang1 *= M_PI / 180.0;
  Point2D end1{centre.x + radius * cos(ang1), centre.y + radius * sin(ang1)};
  ang2 *= M_PI / 180.0;
  Point2D end2{centre.x + radius * cos(ang2), centre.y + radius * sin(ang2)};

  for (size_t i = 0; i < otherAtoms.size(); ++i) {
    if (otherAtoms[i] == atNum) {
      continue;
    }
    auto otherRad = getLassoWidth(dm, otherAtoms[i], lassoNum);
    if ((end1 - dm->atCds_[otherAtoms[i]]).lengthSq() - otherRad * otherRad <
        -1.0e-4) {
      return true;
    }
    if ((end2 - dm->atCds_[otherAtoms[i]]).lengthSq() - otherRad * otherRad <
        -1.0e-4) {
      return true;
    }
  }
  return false;
}

// Find the start of any existing arcs involving atom i, feeding them into
// arcAngs and bangs so as to be able to merge in any new arcs that intersect
// them.  Removes the existing arcs if found.
void addExistingArcs(size_t i,
                     std::vector<std::unique_ptr<DrawShapeArc>> &currArcs,
                     std::vector<double> &arcAngs,
                     std::vector<std::pair<double, unsigned int>> &bangs) {
  for (auto &arc : currArcs) {
    if (arc->atom1_ == static_cast<int>(i)) {
      arcAngs.push_back(arc->ang2_);
      arcAngs.push_back(arc->ang1_);
      if (arc->ang2_ < arc->ang1_) {
        bangs.push_back(
            std::make_pair((arc->ang2_ + arc->ang1_) / 2.0, bangs.size()));
      } else {
        double bang = -360.0 + (arc->ang2_ + arc->ang1_ + 360.0) / 2.0;
        bangs.push_back(std::make_pair(bang, bangs.size()));
      }
      arc.reset();
    }
  }
  currArcs.erase(std::remove_if(currArcs.begin(), currArcs.end(),
                                [](const std::unique_ptr<DrawShapeArc> &arc) {
                                  return !arc;
                                }),
                 currArcs.end());
}
}  // namespace

// ****************************************************************************
void DrawMolMCHLasso::makeIntersectingArcs(
    const std::vector<std::vector<unsigned int>> &intersects, int lassoNum,
    const RDKit::DrawColour &col,
    std::vector<std::unique_ptr<DrawShapeArc>> &currArcs,
    std::vector<DrawShapeArc *> &arcs) const {
  for (size_t i = 0; i < intersects.size(); ++i) {
    if (intersects[i].empty()) {
      continue;
    }
    auto radI = getLassoWidth(this, i, lassoNum);
    if (intersects[i].size() == 1 && intersects[i][0] == i) {
      std::vector<Point2D> pts{atCds_[i], Point2D{radI, radI}};
      DrawShapeArc *arc = new DrawShapeArc(pts, 0.0, 360.0, LINE_WIDTH,
                                           SCALE_LINE_WIDTH, col, false, i);
      arcs.push_back(arc);
    } else {
      std::vector<double> arcAngs;
      std::vector<std::pair<double, unsigned int>> bangs;
      addExistingArcs(i, currArcs, arcAngs, bangs);
      // the first entry is the atom itself, that's how it's set up
      for (size_t j = 1; j < intersects[i].size(); ++j) {
        auto radJ = getLassoWidth(this, intersects[i][j], lassoNum);
        double ang1, ang2, bang;
        calcIntersectingArcAngles(atCds_[i], radI, atCds_[intersects[i][j]],
                                  radJ, ang1, ang2, bang);
        arcAngs.push_back(ang1);
        arcAngs.push_back(ang2);
        bangs.push_back(std::make_pair(bang, bangs.size()));
      }
      std::sort(bangs.begin(), bangs.end());
      std::vector<double> sortedArcAngs;
      for (auto &b : bangs) {
        sortedArcAngs.push_back(arcAngs[2 * b.second]);
        sortedArcAngs.push_back(arcAngs[2 * b.second + 1]);
      }
      // Now make the arcs from i.second to (i + 1).first, and finally from
      // last.first to first.second
      for (size_t j = 0; j < sortedArcAngs.size() - 2; j += 2) {
        std::vector<Point2D> pts{atCds_[i], Point2D{radI, radI}};
        // if either end of this arc we're about to make are inside another
        // of the intersecting rings, skip it.
        if (!areArcEndsInOtherCircle(atCds_[i], radI, sortedArcAngs[j + 1],
                                     sortedArcAngs[j + 2], i, this,
                                     intersects[i], lassoNum)) {
          DrawShapeArc *arc =
              new DrawShapeArc(pts, sortedArcAngs[j + 1], sortedArcAngs[j + 2],
                               LINE_WIDTH, SCALE_LINE_WIDTH, col, false, i);
          arcs.push_back(arc);
        }
      }
      // if either end of this arc we're about to make is inside another
      // of the intersecting rings, skip it.
      std::vector<Point2D> pts{atCds_[i], Point2D{radI, radI}};
      if (!areArcEndsInOtherCircle(atCds_[i], radI, sortedArcAngs.back(),
                                   sortedArcAngs.front(), i, this,
                                   intersects[i], lassoNum)) {
        DrawShapeArc *arc =
            new DrawShapeArc(pts, sortedArcAngs.back(), sortedArcAngs.front(),
                             LINE_WIDTH, SCALE_LINE_WIDTH, col, false, i);
        arcs.push_back(arc);
      }
    }
  }
}

// ****************************************************************************
void DrawMolMCHLasso::extractBondLines(
    size_t lassoNum, const RDKit::DrawColour &col,
    const std::vector<int> &colAtoms,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines,
    std::vector<std::vector<LinePair>> &atomLines) const {
  if (colAtoms.size() > 1) {
    for (size_t i = 0U; i < colAtoms.size() - 1; ++i) {
      if (colAtoms[i] < 0 ||
          static_cast<unsigned int>(colAtoms[i]) >= drawMol_->getNumAtoms()) {
        // there's an error in the colour map.
        continue;
      }
      auto lassoWidthI = getLassoWidth(this, colAtoms[i], lassoNum);
      auto dispI = lassoWidthI * 0.75;
      for (size_t j = i + 1; j < colAtoms.size(); ++j) {
        if (colAtoms[j] < 0 ||
            static_cast<unsigned int>(colAtoms[j]) >= drawMol_->getNumAtoms()) {
          // there's an error in the colour map.
          continue;
        }
        auto lassoWidthJ = getLassoWidth(this, colAtoms[j], lassoNum);
        auto dispJ = lassoWidthJ * 0.75;
        auto bond = drawMol_->getBondBetweenAtoms(colAtoms[i], colAtoms[j]);
        if (bond) {
          if (!mcHighlightBondMap_.empty()) {
            auto bndIt = mcHighlightBondMap_.find(bond->getIdx());
            if (bndIt == mcHighlightBondMap_.end()) {
              continue;
            }
          }
          auto atCdsI = atCds_[colAtoms[i]];
          auto atCdsJ = atCds_[colAtoms[j]];
          auto perp = calcPerpendicular(atCdsI, atCdsJ);
          LinePair thisPair;
          bool skip = false;
          for (auto m : {-1.0, +1.0}) {
            auto p1 = atCdsI + perp * dispI * m;
            auto p2 = atCdsJ + perp * dispJ * m;
            // if the mid point of the line is inside one of the circles
            // then it gets messy because it's then intersecting arcs to deal
            // with.  To duck the problem completely, push the line out to just
            // less than the radii of the circles (just less, so that they still
            // intersect rather than hitting on the tangent)
            auto mid = (p1 + p2) / 2.0;
            if ((atCdsI - mid).lengthSq() < lassoWidthI * lassoWidthI ||
                (atCdsI - mid).lengthSq() < lassoWidthJ * lassoWidthJ) {
              skip = true;
              continue;
            }
            adjustLineEndForEllipse(atCds_[colAtoms[i]], lassoWidthI,
                                    lassoWidthI, p2, p1);
            adjustLineEndForEllipse(atCds_[colAtoms[j]], lassoWidthJ,
                                    lassoWidthJ, p1, p2);
            DrawShapeSimpleLine *pl = new DrawShapeSimpleLine(
                {p1, p2}, LINE_WIDTH, SCALE_LINE_WIDTH, col, colAtoms[i],
                colAtoms[j], bond->getIdx(), noDash);
            lines.emplace_back(pl);
            if (m < 0.0) {
              thisPair.line1 = pl;
            } else {
              thisPair.line2 = pl;
            }
          }
          if (!skip) {
            thisPair.radius = lassoWidthI;
            thisPair.atom = colAtoms[i];
            atomLines[colAtoms[i]].push_back(thisPair);
            thisPair.radius = lassoWidthJ;
            thisPair.atom = colAtoms[j];
            atomLines[colAtoms[j]].push_back(thisPair);
          }
        }
      }
    }
  }
  orderAtomLines(atomLines);
}

namespace {
DrawShapeArc *makeArc(LinePair &lp1, LinePair &lp2, const Point2D &atCds) {
  // Make an arc between line2 of lp1 and line1 of lp2.  If the 2 lines
  // intersect outside the radius of the arc, trim them both back to the
  // point of intersection and return nullptr.
  auto adjustLine = [](DrawShapeSimpleLine *line, const Point2D &pt) -> void {
    double d1 = (line->points_[0] - pt).lengthSq();
    double d2 = (line->points_[1] - pt).lengthSq();
    if (d1 < d2) {
      line->points_[0] = pt;
    } else {
      line->points_[1] = pt;
    }
  };
  Point2D ip;
  if (doLinesIntersect(lp1.line2->points_[0], lp1.line2->points_[1],
                       lp2.line1->points_[0], lp2.line1->points_[1], &ip)) {
    adjustLine(lp1.line2, ip);
    adjustLine(lp2.line1, ip);
    return nullptr;
  }
  DrawShapeArc *arc = new DrawShapeArc(
      {atCds, {lp1.radius, lp1.radius}}, lp1.angle2, lp2.angle1,
      lp1.line1->lineWidth_, lp1.line1->scaleLineWidth_, lp1.line1->lineColour_,
      false, lp1.atom);
  return arc;
}

}  // namespace

// ****************************************************************************
void DrawMolMCHLasso::extractAtomArcs(
    std::vector<std::vector<LinePair>> &atomLines,
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs) const {
  for (size_t k = 0; k < atomLines.size(); ++k) {
    auto &atomLine = atomLines[k];
    if (atomLine.empty()) {
      continue;
    }
    if (atomLine.size() == 1) {
      DrawShapeArc *arc = new DrawShapeArc(
          {atCds_[atomLine[0].atom], {atomLine[0].radius, atomLine[0].radius}},
          atomLine[0].angle2, atomLine[0].angle1, atomLine[0].line1->lineWidth_,
          atomLine[0].line1->scaleLineWidth_, atomLine[0].line1->lineColour_,
          false, atomLine[0].atom);
      arcs.emplace_back(arc);
    } else {
      for (size_t i = 0; i < atomLine.size() - 1; ++i) {
        auto arc =
            makeArc(atomLine[i], atomLine[i + 1], atCds_[atomLine[i].atom]);
        if (arc) {
          arcs.emplace_back(arc);
        }
      }
      auto arc = makeArc(atomLine.back(), atomLine.front(),
                         atCds_[atomLine.front().atom]);
      if (arc) {
        arcs.emplace_back(arc);
      }
    }
  }
}

namespace {
Point2D *adjustLineEnd(const DrawShapeArc &arc, DrawShapeSimpleLine &line) {
  Point2D *adjEnd = nullptr;
  Point2D *fixEnd = nullptr;
  if ((arc.points_[0] - line.points_[0]).lengthSq() >
      (arc.points_[0] - line.points_[1]).lengthSq()) {
    adjEnd = &line.points_[1];
    fixEnd = &line.points_[0];
  } else {
    adjEnd = &line.points_[0];
    fixEnd = &line.points_[1];
  }
  if (fabs((arc.points_[0] - *adjEnd).lengthSq() -
           arc.points_[1].x * arc.points_[1].x) < 1.0e-4) {
    return nullptr;
  }
  adjustLineEndForEllipse(arc.points_[0], arc.points_[1].x, arc.points_[1].x,
                          *fixEnd, *adjEnd);
  return adjEnd;
}

// Calculate the angles of the 2 points  around the centre, measured from the
// X axis, guaranteeing a consistent rotation around the z axis i.e. always
// clockwise or always anti-clockwise.
void calcAnglesFromXAxis(Point2D &centre, Point2D &end1, Point2D &end2,
                         double &ang1, double &ang2) {
  static const Point2D index{1.0, 0.0};
  Point2D rad1, rad2;
  rad1 = centre.directionVector(end1);
  ang1 = 360.0 - rad1.signedAngleTo(index) * 180.0 / M_PI;
  rad2 = centre.directionVector(end2);
  ang2 = 360.0 - rad2.signedAngleTo(index) * 180.0 / M_PI;
  // make sure they're going round in the same direction
  auto crossZ = rad1.x * rad2.y - rad1.y * rad2.x;
  if (crossZ > 0.0) {
    std::swap(ang1, ang2);
  }
}
}  // namespace

// ****************************************************************************
void DrawMolMCHLasso::orderAtomLines(
    std::vector<std::vector<LinePair>> &atomLines) const {
  for (size_t i = 0; i < drawMol_->getNumAtoms(); ++i) {
    if (atomLines[i].empty()) {
      continue;
    }
    std::vector<std::pair<double, size_t>> bondAngles;
    double minBondAngle = 720.0;
    for (size_t j = 0; j < atomLines[i].size(); ++j) {
      // because the same DrawShapeSimpleLine is used for both atoms, the
      // points_[0] may not be the nearer point to atom i.
      int oatom = atomLines[i][j].line1->atom1_ == static_cast<int>(i)
                      ? atomLines[i][j].line1->atom2_
                      : atomLines[i][j].line1->atom1_;
      int pt = 0;
      if ((atCds_[i] - atomLines[i][j].line1->points_[0]).lengthSq() >
          (atCds_[i] - atomLines[i][j].line1->points_[1]).lengthSq()) {
        pt = 1;
      }
      double bang;
      bool swapped;
      calcSubtendedAngles(atomLines[i][j].line1->points_[pt],
                          atomLines[i][j].line2->points_[pt], atCds_[i],
                          atCds_[oatom], atomLines[i][j].angle1,
                          atomLines[i][j].angle2, bang, swapped);
      bondAngles.push_back(std::pair(bang, j));
      if (bang < minBondAngle) {
        minBondAngle = bang;
      }
      if (swapped) {
        std::swap(atomLines[i][j].line1, atomLines[i][j].line2);
      }
    }
    std::for_each(bondAngles.begin(), bondAngles.end(),
                  [&](std::pair<double, size_t> &ba) -> void {
                    ba.first -= minBondAngle;
                  });
    std::sort(bondAngles.begin(), bondAngles.end());
    std::vector<LinePair> newAtomLine;
    for (auto &ba : bondAngles) {
      newAtomLine.push_back(atomLines[i][ba.second]);
    }
    atomLines[i] = newAtomLine;
  }
}

namespace {
std::pair<Point2D, Point2D> getArcEnds(const DrawShapeArc &arc) {
  std::pair<Point2D, Point2D> retVal;
  // for these purposes, it's always a circle, so just use the x
  // radius
  retVal.first.x =
      arc.points_[0].x + arc.points_[1].x * cos(arc.ang1_ * M_PI / 180.0);
  retVal.first.y =
      arc.points_[0].y + arc.points_[1].x * sin(arc.ang1_ * M_PI / 180.0);
  retVal.second.x =
      arc.points_[0].x + arc.points_[1].x * cos(arc.ang2_ * M_PI / 180.0);
  retVal.second.y =
      arc.points_[0].y + arc.points_[1].x * sin(arc.ang2_ * M_PI / 180.0);
  return retVal;
}

}  // namespace
}  // namespace MolDraw2D_detail
}  // namespace RDKit
