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

#include <algorithm>
#include <list>
#include <vector>

#include <GraphMol/RWMol.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/DrawMolMCHLasso.h>

namespace RDKit {
namespace MolDraw2D_detail {

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
  extractAtomArcs(lassoNum, col, colAtoms, arcs);
  std::vector<std::unique_ptr<DrawShapeSimpleLine>> lines;
  extractBondLines(lassoNum, col, colAtoms, lines);
  fixArcsAndLines(arcs, lines);
  fixIntersectingLines(arcs, lines);
  fixIntersectingArcsAndLines(arcs, lines);
  fixProtrudingLines(lines);
  fixOrphanLines(arcs, lines);

  for (auto &it : arcs) {
    highlights_.push_back(std::move(it));
  }
  for (auto &it : lines) {
    highlights_.push_back(std::move(it));
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

// ****************************************************************************
void DrawMolMCHLasso::extractAtomArcs(
    size_t lassoNum, const RDKit::DrawColour &col,
    const std::vector<int> &colAtoms,
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs) const {
  // an empirically derived lineWidth.
  int lineWidth = 3;
  bool scaleLineWidth = true;
  for (auto ca : colAtoms) {
    if (ca < 0 || static_cast<unsigned int>(ca) >= drawMol_->getNumAtoms()) {
      // there's an error in the colour map
      continue;
    }
    double lassoWidth = getLassoWidth(this, ca, lassoNum);
    Point2D radii(lassoWidth, lassoWidth);
    std::vector<Point2D> pts{atCds_[ca], radii};
    DrawShapeArc *ell = new DrawShapeArc(pts, 0.0, 360.0, lineWidth,
                                         scaleLineWidth, col, false, ca);
    arcs.emplace_back(ell);
  }
}

// ****************************************************************************
void DrawMolMCHLasso::extractBondLines(
    size_t lassoNum, const RDKit::DrawColour &col,
    const std::vector<int> &colAtoms,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const {
  int lineWidth = 3;
  bool scaleLineWidth = true;
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
          for (auto m : {-1.0, +1.0}) {
            auto p1 = atCdsI + perp * dispI * m;
            auto p2 = atCdsJ + perp * dispJ * m;
            // if the mid point of the line is inside one of the circles
            // then it gets messy because it's then intersecting arcs to deal
            // with.  To duck the problem completely, push the line out to just
            // less than the radii of the circles (just less, so that they still
            // intersect rather than hitting on the tangent)
            auto mid = (p1 + p2) / 2.0;
            if ((atCdsI - mid).lengthSq() < lassoWidthI * lassoWidthI) {
              p1 = atCdsI + perp * lassoWidthI * 0.99 * m;
              p2 = atCdsJ + perp * lassoWidthJ * 0.99 * m;
            }
            DrawShapeSimpleLine *pl = new DrawShapeSimpleLine(
                {p1, p2}, lineWidth, scaleLineWidth, col, colAtoms[i],
                colAtoms[j], bond->getIdx(), noDash);
            lines.emplace_back(pl);
          }
        }
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
// Take the initial circles and lines, and:
// 1.trim the lines back to the circles they intersect
// 2.make a new set of arcs that go between the pairs of lines, cutting
// out the part that is inside each pair.
void DrawMolMCHLasso::fixArcsAndLines(
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const {
  std::vector<std::unique_ptr<DrawShapeArc>> newArcs;
  std::vector<std::pair<double, double>> lineEndAngles;
  std::list<double> singleLineEndAngles;
  for (auto &arc : arcs) {
    lineEndAngles.clear();
    singleLineEndAngles.clear();
    // The lines are arranged in pairs, and always in the same order within
    // the pair (-ve disp, then +ve disp).
    for (size_t i = 0; i < lines.size(); i += 2) {
      auto &line1 = lines[i];
      if (arc->atom1_ == line1->atom1_ || arc->atom1_ == line1->atom2_) {
        auto &line2 = lines[i + 1];
        auto adjEnd1 = adjustLineEnd(*arc, *line1);
        if (!adjEnd1) {
          continue;
        }
        auto adjEnd2 = adjustLineEnd(*arc, *line2);
        if (!adjEnd2) {
          continue;
        }
        double ang1, ang2;
        calcAnglesFromXAxis(arc->points_[0], *adjEnd1, *adjEnd2, ang1, ang2);
        lineEndAngles.push_back(std::pair{ang2, ang1});
      }
    }
    if (lineEndAngles.empty()) {
      newArcs.emplace_back(std::move(arc));
    } else {
      // check for the end of one going beyond the start of the next.
      // It's easiest if the first angle is 0, so there's no need to worry
      // about going from 360 to 0 on a full rotation.
      double transform = 0.0;
      if (lineEndAngles.size() > 1) {
        std::sort(lineEndAngles.begin(), lineEndAngles.end());
        transform = lineEndAngles.front().first;
        for (auto &ap : lineEndAngles) {
          ap.first -= transform;
          if (ap.first < 0) {
            ap.first += 360.0;
          }
          ap.second -= transform;
          if (ap.second < 0) {
            ap.second += 360.0;
          }
        }
        for (size_t i = 0; i < lineEndAngles.size() - 1; ++i) {
          if (lineEndAngles[i].first == std::numeric_limits<double>::max()) {
            continue;
          }
          // If an arc ends after the start of the next one, merge them to
          // have one from the start of the first to the end of the second.
          if (lineEndAngles[i].second > lineEndAngles[i + 1].first) {
            lineEndAngles[i].second = lineEndAngles[i + 1].second;
            lineEndAngles[i + 1].first = std::numeric_limits<double>::max();
            lineEndAngles[i + 1].second = std::numeric_limits<double>::max();
          }
        }
        // Also do the last and the first, where the overlap is only when
        // it's crossed over 360 again in transformed angles.
        if (lineEndAngles.back().first > lineEndAngles.back().second &&
            lineEndAngles.back().second > lineEndAngles.front().first) {
          lineEndAngles.back().second = lineEndAngles.front().second;
          lineEndAngles.front().first = std::numeric_limits<double>::max();
          lineEndAngles.front().second = std::numeric_limits<double>::max();
        }
      }
      // undo the transformation that started the angles at 0.
      for (const auto &lea : lineEndAngles) {
        if (lea.first != std::numeric_limits<double>::max()) {
          double ang = lea.first + transform;
          if (ang > 360.0) {
            ang -= 360.0;
          }
          singleLineEndAngles.push_back(ang);
          ang = lea.second + transform;
          if (ang > 360.0) {
            ang -= 360.0;
          }
          singleLineEndAngles.push_back(ang);
        }
      }
      // Move the front angle to the back, so the list is now the bits to draw,
      // not the bits to skip
      singleLineEndAngles.push_back(singleLineEndAngles.front());
      singleLineEndAngles.pop_front();
      for (auto ang = singleLineEndAngles.begin();
           ang != singleLineEndAngles.end();) {
        DrawShapeArc *newArc = new DrawShapeArc(
            arc->points_, *ang++, *ang++, arc->lineWidth_, arc->scaleLineWidth_,
            arc->lineColour_, arc->fill_, arc->atom1_);
        newArcs.emplace_back(newArc);
      }
    }
  }
  arcs.clear();
  for (auto &it : newArcs) {
    arcs.push_back(std::move(it));
  }
}

void DrawMolMCHLasso::fixIntersectingLines(
    const std::vector<std::unique_ptr<DrawShapeArc>> &arcs,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const {
  // Sometimes, particularly if the lines have been pushed out further
  // than normal because the circles were too large and overlapped a lot,
  // two lines can intersect, as they extend into a neighbouring circle.

  boost::dynamic_bitset<> intersectingEnds(drawMol_->getNumAtoms());
  for (auto &line1 : lines) {
    for (auto &line2 : lines) {
      if (line1 == line2 || line1->atom1_ > line2->atom1_) {
        continue;
      }
      Point2D ip;
      if (doLinesIntersect(line1->points_[0], line1->points_[1],
                           line2->points_[0], line2->points_[1], &ip)) {
        // for each line, the end that's inside one of the 3 arcs is the one
        // that has to move
        intersectingEnds.set(line1->atom1_);
        intersectingEnds.set(line1->atom2_);
        intersectingEnds.set(line2->atom1_);
        intersectingEnds.set(line2->atom2_);
        for (const auto &arc : arcs) {
          if (intersectingEnds[arc->atom1_]) {
            // The points are 0.99 of arc radius from the centre of the
            // arc.
            auto radsq = arc->points_[1].x * arc->points_[1].x * 0.98;
            if ((arc->points_[0] - line1->points_[0]).lengthSq() < radsq) {
              line1->points_[0] = ip;
            }
            if ((arc->points_[0] - line1->points_[1]).lengthSq() < radsq) {
              line1->points_[1] = ip;
            }
            if ((arc->points_[0] - line2->points_[0]).lengthSq() < radsq) {
              line2->points_[0] = ip;
            }
            if ((arc->points_[0] - line2->points_[1]).lengthSq() < radsq) {
              line2->points_[1] = ip;
            }
          }
        }
        intersectingEnds.reset(line1->atom1_);
        intersectingEnds.reset(line1->atom2_);
        intersectingEnds.reset(line2->atom1_);
        intersectingEnds.reset(line2->atom2_);
      }
    }
  }
}

void DrawMolMCHLasso::fixIntersectingArcsAndLines(
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const {
  // Occasionally, a line will intersect with an arc that it is not
  // directly associated with.  lasso_highlights_7.svg in the catch_test.cpp
  // being an example.
  for (auto &arc : arcs) {
    // The points are 0.99 of arc radius from the centre of the
    // arc.
    auto radsq = arc->points_[1].x * arc->points_[1].x * 0.98;
    for (auto &line : lines) {
      if ((arc->points_[0] - line->points_[0]).lengthSq() < radsq ||
          (arc->points_[0] - line->points_[1]).lengthSq() < radsq) {
        // make a temporary copy of the line because adjustLineEnd moves
        // the points_ in the line, which may be premature if line doesn't
        // intersect this line.
        DrawShapeSimpleLine tmpLine(line->points_, line->atom1_, line->atom2_);
        auto adjEnd1 = adjustLineEnd(*arc, tmpLine);
        if (!adjEnd1) {
          continue;
        }
        double ang1, ang2;
        // lazy reuse of existing function.
        calcAnglesFromXAxis(arc->points_[0], *adjEnd1, *adjEnd1, ang1, ang2);
        double tang1 = 0;
        double tang2 = arc->ang2_ - arc->ang1_;
        if (tang2 < 0.0) {
          tang2 += 360.0;
        }
        double tang3 = ang1 - arc->ang1_;
        if (tang3 < 0.0) {
          tang3 += 360.0;
        }
        if (tang3 > tang1 && tang3 < tang2) {
          // truncate the end that's nearest the intersection
          if (tang3 - tang1 > tang2 - tang3) {
            tang2 = tang3 + arc->ang1_;
            tang1 += arc->ang1_;
          } else {
            tang1 = tang3 + arc->ang1_;
            tang2 += arc->ang1_;
          }
          if (tang1 > 360.0) {
            tang1 -= 360.0;
          }
          if (tang2 > 360.0) {
            tang2 -= 360.0;
          }
          arc->ang1_ = tang1;
          arc->ang2_ = tang2;
          // and finally, move the end of the line
          adjustLineEnd(*arc, *line);
        }
      }
    }
  }
}

void DrawMolMCHLasso::fixProtrudingLines(
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const {
  // lasso_highlight_7.svg also had the problem where two lines didn't
  // intersect, but one protruded beyond the end of the other inside the
  // lasso.
  for (auto &line1 : lines) {
    for (auto &line2 : lines) {
      auto d1_0 = (line1->points_[0] - line2->points_[0]).length();
      auto d2_0 = (line2->points_[0] - line1->points_[1]).length();
      auto l_0 = (line1->points_[0] - line1->points_[1]).length();
      if (fabs(d1_0 + d2_0 - l_0) < 1.0e-4) {
        if (d1_0 < d2_0) {
          line1->points_[0] = line2->points_[0];
        } else {
          line1->points_[1] = line2->points_[0];
        }
      }
      auto d1_1 = (line1->points_[0] - line2->points_[1]).length();
      auto d2_1 = (line2->points_[1] - line1->points_[1]).length();
      auto l_1 = (line1->points_[0] - line1->points_[1]).length();
      if (fabs(d1_1 + d2_1 - l_1) < 1.0e-4) {
        if (d1_1 < d2_1) {
          line1->points_[0] = line2->points_[1];
        } else {
          line1->points_[1] = line2->points_[1];
        }
      }
    }
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
void DrawMolMCHLasso::fixOrphanLines(
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) {
  // lasso_highlights_7.svg had a line close to an arc at
  // one end, but not at the other.  Such lines are clearly
  // artifacts that need to be removed.  Takes out all lines
  // that aren't within a tolerance of something at both ends.
  std::vector<std::pair<Point2D, Point2D>> arcEnds;
  for (const auto &arc : arcs) {
    arcEnds.push_back(getArcEnds(*arc));
  }
  // This tolerance was arrived at empirically.  It's enough
  // to fix lasso_highlights_7.svg without breaking it and
  // lasso_highlights_6.svg.  A slightly tighter tolerance
  // (e.g. 0.003) causes green lines associated with bonds 9
  // and 12 to be removed incorrectly in both of these.
  static const double tol = 0.005;
  for (auto &line1 : lines) {
    bool attached0 = false, attached1 = false;
    for (const auto &arcEnd : arcEnds) {
      if ((line1->points_[0] - arcEnd.first).lengthSq() < tol) {
        attached0 = true;
      }
      if ((line1->points_[1] - arcEnd.first).lengthSq() < tol) {
        attached1 = true;
      }
      if ((line1->points_[0] - arcEnd.second).lengthSq() < tol) {
        attached0 = true;
      }
      if ((line1->points_[1] - arcEnd.second).lengthSq() < tol) {
        attached1 = true;
      }
    }
    if (!attached0 || !attached1) {
      for (auto &line2 : lines) {
        if (line1 == line2 || !line1 || !line2) {
          continue;
        }
        if ((line1->points_[0] - line2->points_[0]).lengthSq() < tol) {
          attached0 = true;
        }
        if ((line1->points_[0] - line2->points_[1]).lengthSq() < tol) {
          attached0 = true;
        }
        if ((line1->points_[1] - line2->points_[0]).lengthSq() < tol) {
          attached1 = true;
        }
        if ((line1->points_[1] - line2->points_[1]).lengthSq() < tol) {
          attached1 = true;
        }
        if (attached0 && attached1) {
          break;
        }
      }
    }
    if (!attached0 || !attached1) {
      line1.reset();
    }
  }
  lines.erase(std::remove_if(
                  lines.begin(), lines.end(),
                  [](const std::unique_ptr<DrawShapeSimpleLine> &line) -> bool {
                    return !line;
                  }),
              lines.end());
}
}  // namespace MolDraw2D_detail
}  // namespace RDKit
