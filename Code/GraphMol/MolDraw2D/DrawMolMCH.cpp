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

#include <GraphMol/RWMol.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/DrawMolMCH.h>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawMolMCH::DrawMolMCH(
    const ROMol &mol, const std::string &legend, int width, int height,
    MolDrawOptions &drawOptions, DrawText &textDrawer,
    const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
    const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
    const std::map<int, double> &highlight_radii,
    const std::map<int, int> &highlight_linewidth_multipliers, int confId)
    : DrawMol(mol, legend, width, height, drawOptions, textDrawer, nullptr,
              nullptr, nullptr, nullptr, nullptr, &highlight_radii, confId),
      mcHighlightAtomMap_(highlight_atom_map),
      mcHighlightBondMap_(highlight_bond_map),
      highlightLinewidthMultipliers_(highlight_linewidth_multipliers) {}

// ****************************************************************************
void DrawMolMCH::extractHighlights(double scale) {
  DrawMol::extractHighlights(scale);
  extractMCHighlights();
}

// ****************************************************************************
void DrawMolMCH::extractMCHighlights() {
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
void DrawMolMCH::makeBondHighlights(
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
void DrawMolMCH::makeAtomHighlights(
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
        arc_stop = arc_stop >= 360.0 ? arc_stop - 360.0 : arc_stop;
        std::vector<Point2D> pts{centre, Point2D(xradius, yradius)};
        DrawShape *arc = new DrawShapeArc(
            pts, arc_start, arc_stop, lineWidth, true, ha.second[i],
            drawOptions_.fillHighlights, ha.first);
        atomHighlights.emplace_back(arc);
        arc_start += arc_size;
        arc_start = arc_start >= 360.0 ? arc_start - 360.0 : arc_start;
      }
    }
  }
}

// ****************************************************************************
void DrawMolMCH::adjustLineEndForHighlight(int at_idx, Point2D p1,
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
void DrawMolMCH::calcSymbolEllipse(unsigned int atomIdx, Point2D &centre,
                                   double &xradius, double &yradius) const {
  centre = atCds_[atomIdx];
  xradius = drawOptions_.highlightRadius;
  yradius = xradius;
  if (highlightRadii_.find(atomIdx) != highlightRadii_.end()) {
    xradius = highlightRadii_.find(atomIdx)->second;
    yradius = xradius;
  }

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
void DrawMolMCH::fixHighlightJoinProblems(
    std::vector<std::unique_ptr<DrawShape>> &atomHighlights,
    std::vector<std::unique_ptr<DrawShape>> &bondHighlights) {
  std::vector<unsigned int> fettledAtoms;
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
          fettledAtoms.push_back(i);
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
            double upscaler = 1.25 * std::min(dist1, dist2) / maxrad;
            if (upscaler <= 1.0) {
              upscaler = 1.1;
            }
            atomHL->points_[1] *= upscaler;
            ins = doesLineIntersectArc(atomHL->points_[0], atomHL->points_[1].x,
                                       atomHL->points_[1].y, 0.0, 360, 0.0,
                                       bondHL->points_[0], bondHL->points_[p2]);
          }
        }
      }
    }
  }
  std::sort(fettledAtoms.begin(), fettledAtoms.end());
  fettledAtoms.erase(std::unique(fettledAtoms.begin(), fettledAtoms.end()),
                     fettledAtoms.end());
  for (unsigned int i = 0; i < fettledAtoms.size(); ++i) {
    // now adjust all the other bond highlights that end on this atom
    auto &atomHL = atomHighlights[fettledAtoms[i]];
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
    p1 += centre;
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
      p1 += centre;
      p2 += centre;
      return;
    }
    p2 = t_to_point(t);
  }
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit