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
// This is based on a suggestion and code from Christian Feldmann.
// It was discussion #4607.  His Python implementation (which I haven't
// followed to any great extent) is at
// https://github.com/c-feldmann/lassohighlight

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
  extractAtomColourLists(colours, colourAtoms);
  for (size_t i = 0U; i < colours.size(); ++i) {
    //    if (i != 2) {
    //      continue;
    //    }
    drawLasso(i, colours[i], colourAtoms[i]);
  }
}

// ****************************************************************************
void DrawMolMCHLasso::extractAtomColourLists(
    std::vector<DrawColour> &colours,
    std::vector<std::vector<int>> &colourAtoms) const {
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
  std::cout << "adjEnd distance : " << (arc.points_[0] - *adjEnd).lengthSq()
            << " vs " << arc.points_[1].x * arc.points_[1].x << std::endl;
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
void DrawMolMCHLasso::fixArcsAndLines(
    std::vector<std::unique_ptr<DrawShapeArc>> &arcs,
    std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const {
  static const Point2D index{1.0, 0.0};

  std::vector<std::unique_ptr<DrawShapeArc>> newArcs;
  std::list<double> lineEndAngles;
  for (auto &arc : arcs) {
    //    if (arc->atom1_ != 10 && arc->atom1_ != 12) {
    //      continue;
    //    }
    lineEndAngles.clear();
    for (auto &line1 : lines) {
      if (arc->atom1_ == line1->atom1_ || arc->atom1_ == line1->atom2_) {
        // find the other line
        for (auto &line2 : lines) {
          if (line1 != line2 && line1->atom1_ == line2->atom1_ &&
              line1->atom2_ == line2->atom2_) {
            std::cout << "arc at " << arc->atom1_ << " and lines between "
                      << line1->atom1_ << " and " << line1->atom2_ << " : "
                      << arc->ang1_ << " to " << arc->ang2_ << std::endl;
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
            auto crossZ = rad1.x * rad2.y - rad1.y * rad2.x;
            std::cout << "crossZ : " << crossZ << std::endl;
            auto mid = (*adjEnd1 + *adjEnd2) / 2.0;
            auto rad = arc->points_[0].directionVector(mid);
            auto ang3 = 360.0 - rad.signedAngleTo(index) * 180.0 / M_PI;
            std::cout << ang1 << " : " << ang3 << " : " << ang2 << std::endl;
            if (crossZ > 0.0) {
              std::cout << "swapping 1 and 2" << std::endl;
              std::swap(ang1, ang2);
            }
            //            auto insPos = lineEndAngles.begin();
            //            while (insPos != lineEndAngles.end() && *insPos <
            //            ang1) {
            //              ++insPos;
            //            }
            //            lineEndAngles.insert(insPos, ang1);
            //            lineEndAngles.insert(insPos, ang2);
            lineEndAngles.push_back(ang2);
            lineEndAngles.push_back(ang1);
#if 0
              std::cout << "centre : " << arc->points_[0]
                        << " to adjEnd1 : " << *adjEnd1 << " and : " << *adjEnd2
                        << std::endl;
              auto rad1 = arc->points_[0].directionVector(*adjEnd1);
              auto ang1 = 360.0 - rad1.signedAngleTo(index) * 180.0 / M_PI;
              auto rad2 = arc->points_[0].directionVector(*adjEnd2);
              auto ang2 = 360.0 - rad2.signedAngleTo(index) * 180.0 / M_PI;
              if (ang1 > ang2) {
                std::swap(ang1, ang2);
              }
              std::cout << "angles : " << ang1 << " and " << ang2 << std::endl;
              if (ang1 < arc->ang1_ || ang1 > arc->ang2_ || ang2 < arc->ang1_ ||
                  ang2 > arc->ang2_) {
                std::cout << "skipping" << std::endl;
                continue;
              }
              DrawShapeArc *newArc =
                  new DrawShapeArc(arc->points_, ang2, arc->ang2_,
                                   arc->lineWidth_, arc->scaleLineWidth_,
                                   arc->lineColour_, arc->fill_, arc->atom1_);
              newArcs.emplace_back(newArc);
              arc->ang2_ = ang1;
              std::cout << "angles now " << arc->ang1_ << " to " << arc->ang2_
                        << " and " << newArc->ang1_ << " to " << newArc->ang2_
                        << std::endl;
#endif
          }
        }
      }
    }
    std::cout << "Angles :";
    for (auto a : lineEndAngles) {
      std::cout << " " << a;
    }
    std::cout << std::endl;
    lineEndAngles.push_back(lineEndAngles.front());
    lineEndAngles.pop_front();
    std::cout << "Fiddled Angles :";
    for (auto a : lineEndAngles) {
      std::cout << " " << a;
    }
    std::cout << std::endl;
    for (auto ang = lineEndAngles.begin(); ang != lineEndAngles.end();) {
      DrawShapeArc *newArc = new DrawShapeArc(
          arc->points_, *ang++, *ang++, arc->lineWidth_, arc->scaleLineWidth_,
          arc->lineColour_, arc->fill_, arc->atom1_);
      newArcs.emplace_back(newArc);
      //      break;
    }
  }
  arcs.clear();
  for (auto &it : newArcs) {
    arcs.push_back(std::move(it));
  }
}
}  // namespace MolDraw2D_detail
}  // namespace RDKit
