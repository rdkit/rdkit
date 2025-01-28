//
//  Copyright (C) 2023 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <unordered_set>

#include <GraphMol/RWMol.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/DrawMolMCHCircleAndLine.h>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawMolMCHCircleAndLine::DrawMolMCHCircleAndLine(
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
void DrawMolMCHCircleAndLine::extractHighlights(double scale) {
  DrawMol::extractHighlights(scale);
  extractMCHighlights();
}

// ****************************************************************************
void DrawMolMCHCircleAndLine::extractMCHighlights() {
  std::vector<std::unique_ptr<DrawShape>> bondHighlights;
  makeBondHighlights(bondHighlights);
  std::vector<std::unique_ptr<DrawShape>> atomHighlights;
  makeAtomHighlights(atomHighlights);
  fixHighlightJoinProblems(atomHighlights, bondHighlights);

  for (auto &it : bondHighlights) {
    highlights_.push_back(std::move(it));
  }
  bondHighlights.clear();
  for (auto &it : atomHighlights) {
    highlights_.push_back(std::move(it));
  }
  atomHighlights.clear();
}

// ****************************************************************************
void DrawMolMCHCircleAndLine::makeBondHighlights(
    std::vector<std::unique_ptr<DrawShape>> &bondHighlights) {
  for (auto hb : mcHighlightBondMap_) {
    auto bond_idx = hb.first;
    auto lineWidth = drawOptions_.bondLineWidth;
    if (!drawOptions_.fillHighlights) {
      lineWidth = getHighlightBondWidth(drawOptions_, bond_idx,
                                        &highlightLinewidthMultipliers_);
    }
    auto bond = drawMol_->getBondWithIdx(bond_idx);
    auto at1_idx = bond->getBeginAtomIdx();
    auto at2_idx = bond->getEndAtomIdx();
    auto at1_cds = atCds_[at1_idx];
    auto at2_cds = atCds_[at2_idx];
    auto perp = calcPerpendicular(at1_cds, at2_cds);
    auto rad = 0.7 * drawOptions_.highlightRadius;

    auto make_adjusted_line = [&](Point2D p1, Point2D p2,
                                  const DrawColour &col) {
      adjustLineEndForHighlight(at1_idx, p2, p1);
      adjustLineEndForHighlight(at2_idx, p1, p2);
      std::vector<Point2D> pts{p1, p2};
      DrawShape *pl = new DrawShapeSimpleLine(
          pts, lineWidth, drawOptions_.scaleBondWidth, col, at1_idx, at2_idx,
          bond->getIdx(), noDash);
      bondHighlights.emplace_back(pl);
    };

    if (hb.second.size() < 2) {
      DrawColour col;
      if (hb.second.empty()) {
        col = drawOptions_.highlightColour;
      } else {
        col = hb.second.front();
      }
      if (drawOptions_.fillHighlights) {
        std::vector<Point2D> line_pts;
        line_pts.push_back(at1_cds + perp * rad);
        line_pts.push_back(at2_cds + perp * rad);
        line_pts.push_back(at2_cds - perp * rad);
        line_pts.push_back(at1_cds - perp * rad);
        DrawShape *pl = new DrawShapePolyLine(
            line_pts, lineWidth, drawOptions_.scaleBondWidth, col, true,
            at1_idx, at2_idx, bond->getIdx(), noDash);
        bondHighlights.emplace_back(pl);
      } else {
        make_adjusted_line(at1_cds + perp * rad, at2_cds + perp * rad, col);
        make_adjusted_line(at1_cds - perp * rad, at2_cds - perp * rad, col);
      }
    } else {
      auto col_rad = 2.0 * rad / hb.second.size();
      if (drawOptions_.fillHighlights) {
        auto p1 = at1_cds - perp * rad;
        auto p2 = at2_cds - perp * rad;
        std::vector<Point2D> line_pts;
        for (size_t i = 0; i < hb.second.size(); ++i) {
          line_pts.clear();
          line_pts.push_back(p1);
          line_pts.push_back(p1 + perp * col_rad);
          line_pts.push_back(p2 + perp * col_rad);
          line_pts.push_back(p2);
          DrawShape *pl = new DrawShapePolyLine(
              line_pts, lineWidth, drawOptions_.scaleBondWidth, hb.second[i],
              true, at1_idx, at2_idx, bond->getIdx(), noDash);
          bondHighlights.emplace_back(pl);
          p1 += perp * col_rad;
          p2 += perp * col_rad;
        }
      } else {
        std::vector<DrawColour> cols{hb.second};
        if (cols.size() % 2) {
          make_adjusted_line(at1_cds, at2_cds, cols[0]);
          cols.erase(cols.begin());
        }
        int step = 0;
        for (size_t i = 0; i < cols.size(); ++i) {
          // draw even numbers from the bottom, odd from the top
          auto offset = perp * (rad - step * col_rad);
          if (!(i % 2)) {
            make_adjusted_line(at1_cds - offset, at2_cds - offset, cols[i]);
          } else {
            make_adjusted_line(at1_cds + offset, at2_cds + offset, cols[i]);
            step++;
          }
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMolMCHCircleAndLine::makeAtomHighlights(
    std::vector<std::unique_ptr<DrawShape>> &atomHighlights) {
  for (auto &ha : mcHighlightAtomMap_) {
    double xradius, yradius;
    Point2D centre;
    int lineWidth = getHighlightBondWidth(drawOptions_, -1, nullptr);
    calcSymbolEllipse(ha.first, centre, xradius, yradius);
    if (ha.second.size() == 1) {
      Point2D radii(xradius, yradius);
      std::vector<Point2D> pts{centre, radii};
      DrawShape *ell =
          new DrawShapeEllipse(pts, lineWidth, true, ha.second.front(),
                               drawOptions_.fillHighlights, ha.first);
      atomHighlights.emplace_back(ell);
    } else {
      auto arc_size = 360.0 / double(ha.second.size());
      auto arc_start = 270.0;
      for (size_t i = 0; i < ha.second.size(); ++i) {
        auto arc_stop = arc_start + arc_size;
        if (arc_stop >= 360.0) {
          arc_stop -= 360.0;
        }
        std::vector<Point2D> pts{centre, Point2D(xradius, yradius)};
        DrawShape *arc = new DrawShapeArc(
            pts, arc_start, arc_stop, lineWidth, true, ha.second[i],
            drawOptions_.fillHighlights, ha.first);
        atomHighlights.emplace_back(arc);
        arc_start += arc_size;
        if (arc_start >= 360.0) {
          arc_start -= 360.0;
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMolMCHCircleAndLine::adjustLineEndForHighlight(int at_idx, Point2D p1,
                                                        Point2D &p2) const {
  // this code is transliterated from
  // http://csharphelper.com/blog/2017/08/calculate-where-a-line-segment-and-an-ellipse-intersect-in-c/
  // which has it in C#
  double xradius, yradius;
  Point2D centre;
  calcSymbolEllipse(at_idx, centre, xradius, yradius);
  if (xradius < 1.0e-6 || yradius < 1.0e-6) {
    return;
  }
  adjustLineEndForEllipse(centre, xradius, yradius, p1, p2);
}

// ****************************************************************************
void DrawMolMCHCircleAndLine::calcSymbolEllipse(unsigned int atomIdx,
                                                Point2D &centre,
                                                double &xradius,
                                                double &yradius) const {
  PRECONDITION(atomIdx < atCds_.size() && atomIdx < atomLabels_.size(), "")
  centre = atCds_[atomIdx];
  getAtomRadius(atomIdx, xradius, yradius);
  if (drawOptions_.atomHighlightsAreCircles || !atomLabels_[atomIdx] ||
      atomLabels_[atomIdx]->symbol_.empty()) {
    return;
  }

  double x_min, y_min, x_max, y_max;
  x_min = y_min = std::numeric_limits<double>::max();
  x_max = y_max = std::numeric_limits<double>::lowest();
  atomLabels_[atomIdx]->findExtremes(x_min, x_max, y_min, y_max);

  static const double root_2 = sqrt(2.0);
  xradius = std::max(xradius, root_2 * 0.5 * (x_max - x_min));
  yradius = std::max(yradius, root_2 * 0.5 * (y_max - y_min));
  centre.x = 0.5 * (x_max + x_min);
  centre.y = 0.5 * (y_max + y_min);
}

// ****************************************************************************
void DrawMolMCHCircleAndLine::fixHighlightJoinProblems(
    std::vector<std::unique_ptr<DrawShape>> &atomHighlights,
    std::vector<std::unique_ptr<DrawShape>> &bondHighlights) {
  // Sometimes, the bond highlights went outside the atom highlight ellipse
  // they were supposed to be intersecting.  This could occur particularly when
  // the atom highlights were circles rather than ellipses.  This code detects
  // these cases and makes the atom highlight radius larger.
  std::unordered_set<unsigned int> fettledAtoms;
  for (unsigned int i = 0; i < atomHighlights.size(); ++i) {
    auto &atomHL = atomHighlights[i];
    for (auto &bondHL : bondHighlights) {
      if (atomHL->atom1_ == bondHL->atom1_ ||
          atomHL->atom1_ == bondHL->atom2_) {
        // A multicoloured bond highlight is a polyline, with 4 points, so 2 is
        // opposite 0.  Normally a bond highlight is a single line, with
        // only 2 points.
        int p2 = bondHL->points_.size() == 2 ? 1 : 2;
        bool ins = doesLineIntersectArc(
            atomHL->points_[0], atomHL->points_[1].x, atomHL->points_[1].y, 0.0,
            360.0, 0.0, bondHL->points_[0], bondHL->points_[p2]);
        if (!ins) {
          fettledAtoms.insert(i);
          while (!ins) {
            double dist1 = (atomHL->points_[0] - bondHL->points_[0]).length();
            double dist2 = (atomHL->points_[0] - bondHL->points_[p2]).length();
            double maxrad =
                std::max(atomHL->points_[1].x, atomHL->points_[1].y);
            // If it's a filled highlight, the bond highlight won't have been
            // truncated at the outer radius of the circle, because there's no
            // need.  Set the dist2 to the radius now.
            if (dist2 < 1.0e-4) {
              dist2 = maxrad;
            }
            // This scales the radii of the ellipse by 1.25 times the
            // ratio of the shorter distance of the bond highlight ends
            // from the atom and the larger radius of the ellipse.  It's
            // to make sure the ellipse is definitely wide enough that the
            // bond highlight intersects it.
            double upscaler = 1.25 * std::min(dist1, dist2) / maxrad;
            // We certainly don't want to make the ellipse smaller.
            if (upscaler <= 1.0) {
              upscaler = 1.1;
            }
            // For radii of the ellipse are in points_[1].
            atomHL->points_[1] *= upscaler;
            ins = doesLineIntersectArc(atomHL->points_[0], atomHL->points_[1].x,
                                       atomHL->points_[1].y, 0.0, 360, 0.0,
                                       bondHL->points_[0], bondHL->points_[p2]);
          }
        }
      }
    }
  }
  for (auto fa : fettledAtoms) {
    // now adjust all the other bond highlights that end on this atom
    auto &atomHL = atomHighlights[fa];
    for (auto &bondHL : bondHighlights) {
      if (atomHL->atom1_ == bondHL->atom1_ ||
          atomHL->atom1_ == bondHL->atom2_) {
        if ((atomHL->points_[0] - bondHL->points_[0]).lengthSq() <
            (atomHL->points_[0] - bondHL->points_[1]).lengthSq()) {
          adjustLineEndForEllipse(atomHL->points_[0], atomHL->points_[1].x,
                                  atomHL->points_[1].y, bondHL->points_[1],
                                  bondHL->points_[0]);
        } else {
          adjustLineEndForEllipse(atomHL->points_[0], atomHL->points_[1].x,
                                  atomHL->points_[1].y, bondHL->points_[0],
                                  bondHL->points_[1]);
        }
      }
    }
  }
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit