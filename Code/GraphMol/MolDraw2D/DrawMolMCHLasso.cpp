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
void DrawMolMCHLasso::extractHighlights(double scale) { extractMCHighlights(); }

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
  for (const auto &cm : mcHighlightAtomMap_) {
    for (const auto &col : cm.second) {
      auto cpos = std::find(colours.begin(), colours.end(), col);
      if (cpos == colours.end()) {
        colours.push_back(col);
        colourAtoms.push_back(std::vector<int>(1, cm.first));
      } else {
        colourAtoms[std::distance(colours.begin(), cpos)].push_back(cm.first);
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
void DrawMolMCHLasso::drawLasso(size_t colNum, const RDKit::DrawColour &col,
                                const std::vector<int> &colAtoms) {
  if (colAtoms.empty()) {
    return;
  }
  std::vector<std::unique_ptr<DrawShapeArc>> arcs;
  extractAtomArcs(colNum, col, colAtoms, arcs);
  std::vector<std::unique_ptr<DrawShapeSimpleLine>> lines;
  extractBondLines(colNum, col, colAtoms, lines);
  fixArcsAndLines(arcs, lines);

  for (auto &it : arcs) {
    highlights_.push_back(std::move(it));
  }
  for (auto &it : lines) {
    highlights_.push_back(std::move(it));
  }
}

// ****************************************************************************
void DrawMolMCHLasso::extractAtomArcs(
    size_t colNum, const RDKit::DrawColour &col,
    const std::vector<int> &colAtoms,
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs) const {
  double xradius, yradius;
  getAtomRadius(colAtoms.front(), xradius, yradius);
  xradius += xradius * colNum * 0.5;
  Point2D radii(xradius, xradius);
  for (auto ca : colAtoms) {
    std::vector<Point2D> pts{atCds_[ca], radii};
    DrawShapeArc *ell = new DrawShapeArc(
        pts, 0.0, 360.0, drawOptions_.bondLineWidth, true, col, false, ca);
    arcs.emplace_back(ell);
  }
}

// ****************************************************************************
void DrawMolMCHLasso::extractBondLines(
    size_t colNum, const RDKit::DrawColour &col,
    const std::vector<int> &colAtoms,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const {
  double xradius, yradius;
  getAtomRadius(colAtoms.front(), xradius, yradius);
  xradius += xradius * colNum;
  if (colAtoms.size() > 1) {
    for (size_t i = 0U; i < colAtoms.size() - 1; ++i) {
      for (size_t j = i + 1; j < colAtoms.size(); ++j) {
        auto bond = drawMol_->getBondBetweenAtoms(colAtoms[i], colAtoms[j]);
        if (bond) {
          auto at1_cds = atCds_[colAtoms[i]];
          auto at2_cds = atCds_[colAtoms[j]];
          auto perp = calcPerpendicular(at1_cds, at2_cds);
          auto disp = xradius / 2;
          for (auto m : {1.0, -1.0}) {
            Point2D p1 = at1_cds + perp * disp * m;
            Point2D p2 = at2_cds + perp * disp * m;
            std::vector<Point2D> pts{p1, p2};
            DrawShapeSimpleLine *pl = new DrawShapeSimpleLine(
                pts, drawOptions_.bondLineWidth, drawOptions_.scaleBondWidth,
                col, colAtoms[i], colAtoms[j], bond->getIdx(), noDash);
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
}  // namespace

// ****************************************************************************
// Take the initial circles and lines, and:
// 1.trim the lines back to the circles they intersect
// 2.make a new set of arcs that go between the pairs of lines, cutting
// out the part that is inside each pair.
void DrawMolMCHLasso::fixArcsAndLines(
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const {
  static const Point2D index{1.0, 0.0};

  std::vector<std::unique_ptr<DrawShapeArc>> newArcs;
  std::vector<std::pair<double, double>> lineEndAngles;
  std::list<double> singleLineEndAngles;
  for (auto &arc : arcs) {
    lineEndAngles.clear();
    singleLineEndAngles.clear();
    // find 2 lines that intersect with this arc.  The arc should not
    // be drawn between them.
    for (auto &line1 : lines) {
      if (arc->atom1_ == line1->atom1_ || arc->atom1_ == line1->atom2_) {
        for (auto &line2 : lines) {
          if (line1 != line2 && line1->atom1_ == line2->atom1_ &&
              line1->atom2_ == line2->atom2_) {
            auto adjEnd1 = adjustLineEnd(*arc, *line1);
            double ang1, ang2;
            Point2D rad1, rad2;
            if (adjEnd1) {
              rad1 = arc->points_[0].directionVector(*adjEnd1);
              ang1 = 360.0 - rad1.signedAngleTo(index) * 180.0 / M_PI;
            } else {
              continue;
            }
            auto adjEnd2 = adjustLineEnd(*arc, *line2);
            if (adjEnd2) {
              rad2 = arc->points_[0].directionVector(*adjEnd2);
              ang2 = 360.0 - rad2.signedAngleTo(index) * 180.0 / M_PI;
            } else {
              continue;
            }
            // make sure they're going round in the same direction
            auto crossZ = rad1.x * rad2.y - rad1.y * rad2.x;
            if (crossZ > 0.0) {
              std::swap(ang1, ang2);
            }
            lineEndAngles.push_back(std::pair{ang2, ang1});
          }
        }
      }
    }
    // make sure they go round in ascending order.
    std::sort(lineEndAngles.begin(), lineEndAngles.end());
    // Move the front angle to the back, so the list is now the bits to draw,
    // not the bits to skip.
    for (const auto &lea : lineEndAngles) {
      singleLineEndAngles.push_back(lea.first);
      singleLineEndAngles.push_back(lea.second);
    }
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
  arcs.clear();
  for (auto &it : newArcs) {
    arcs.push_back(std::move(it));
  }
}
}  // namespace MolDraw2D_detail
}  // namespace RDKit
