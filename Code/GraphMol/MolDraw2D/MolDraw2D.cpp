//
//  Copyright (C) 2014-2022 David Cosgrove and Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (AstraZeneca)
// 27th May 2014
//
// Extensively modified by Greg Landrum
//

#include <GraphMol/QueryOps.h>
#include <GraphMol/MolDraw2D/AtomSymbol.h>
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawMolMCH.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/FileParsers/MolSGroupParsing.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <Geometry/point.h>
#include <Geometry/Transform2D.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/Matrix.h>

#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/MolEnumerator/LinkNode.h>

#include <Geometry/Transform3D.h>

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <memory>

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/format.hpp>

using namespace boost;
using namespace std;

namespace RDKit {

// ****************************************************************************
MolDraw2D::MolDraw2D(int width, int height, int panelWidth, int panelHeight)
    : width_(width),
      height_(height),
      panel_width_(panelWidth > 0 ? panelWidth : width),
      panel_height_(panelHeight > 0 ? panelHeight : height),
      legend_height_(0),
      scale_(1.0),
      fontScale_(1.0),
      x_offset_(0),
      y_offset_(0),
      fill_polys_(true),
      activeMolIdx_(-1),
      activeAtmIdx1_(-1),
      activeAtmIdx2_(-1),
      activeBndIdx_(-1) {}

// ****************************************************************************
MolDraw2D::~MolDraw2D() {}

// ****************************************************************************
void MolDraw2D::drawMolecule(const ROMol &mol, const std::string &legend,
                             const vector<int> *highlight_atoms,
                             const vector<int> *highlight_bonds,
                             const map<int, DrawColour> *highlight_atom_map,
                             const map<int, DrawColour> *highlight_bond_map,
                             const std::map<int, double> *highlight_radii,
                             int confId) {
  setupTextDrawer();
  drawMols_.emplace_back(new MolDraw2D_detail::DrawMol(
      mol, legend, panelWidth(), panelHeight(), drawOptions(), *text_drawer_,
      highlight_atoms, highlight_bonds, highlight_atom_map, highlight_bond_map,
      nullptr, highlight_radii, supportsAnnotations(), confId));
  drawMols_.back()->setOffsets(x_offset_, y_offset_);
  drawMols_.back()->createDrawObjects();
  fixVariableDimensions(*drawMols_.back());
  ++activeMolIdx_;
  startDrawing();
  drawTheMolecule(*drawMols_.back());
}

// ****************************************************************************
void MolDraw2D::drawMolecule(const ROMol &mol,
                             const vector<int> *highlight_atoms,
                             const map<int, DrawColour> *highlight_atom_map,
                             const std::map<int, double> *highlight_radii,
                             int confId) {
  drawMolecule(mol, "", highlight_atoms, highlight_atom_map, highlight_radii,
               confId);
}

// ****************************************************************************
void MolDraw2D::drawMolecule(const ROMol &mol, const std::string &legend,
                             const vector<int> *highlight_atoms,
                             const map<int, DrawColour> *highlight_atom_map,
                             const std::map<int, double> *highlight_radii,
                             int confId) {
  vector<int> highlight_bonds;
  if (highlight_atoms) {
    MolDraw2D_detail::getBondHighlightsForAtoms(mol, *highlight_atoms,
                                                highlight_bonds);
  }
  drawMolecule(mol, legend, highlight_atoms, &highlight_bonds,
               highlight_atom_map, nullptr, highlight_radii, confId);
}

// ****************************************************************************
// keeping this one to preserve the public API, although it is no longer
// very useful.
void MolDraw2D::drawMolecule(const ROMol &mol,
                             const vector<int> *highlight_atoms,
                             const vector<int> *highlight_bonds,
                             const map<int, DrawColour> *highlight_atom_map,
                             const map<int, DrawColour> *highlight_bond_map,
                             const std::map<int, double> *highlight_radii,
                             int confId) {
  drawMolecule(mol, "", highlight_atoms, highlight_bonds, highlight_atom_map,
               highlight_bond_map, highlight_radii, confId);
}

// ****************************************************************************
void MolDraw2D::drawMoleculeWithHighlights(
    const ROMol &mol, const string &legend,
    const map<int, vector<DrawColour>> &highlight_atom_map,
    const map<int, vector<DrawColour>> &highlight_bond_map,
    const map<int, double> &highlight_radii,
    const map<int, int> &highlight_linewidth_multipliers, int confId) {
  setupTextDrawer();
  MolDraw2D_detail::DrawMol *dm = new MolDraw2D_detail::DrawMolMCH(
      mol, legend, panelWidth(), panelHeight(), drawOptions(), *text_drawer_,
      highlight_atom_map, highlight_bond_map, highlight_radii,
      highlight_linewidth_multipliers, confId);
  drawMols_.emplace_back(dm);
  drawMols_.back()->createDrawObjects();
  fixVariableDimensions(*drawMols_.back());
  ++activeMolIdx_;
  startDrawing();
  drawTheMolecule((*drawMols_.back()));
  return;
}

// ****************************************************************************
void MolDraw2D::drawMolecules(
    const std::vector<ROMol *> &mols, const std::vector<std::string> *legends,
    const std::vector<std::vector<int>> *highlight_atoms,
    const std::vector<std::vector<int>> *highlight_bonds,
    const std::vector<std::map<int, DrawColour>> *highlight_atom_maps,
    const std::vector<std::map<int, DrawColour>> *highlight_bond_maps,
    const std::vector<std::map<int, double>> *highlight_radii,
    const std::vector<int> *confIds) {
  PRECONDITION(!legends || legends->size() == mols.size(), "bad size");
  PRECONDITION(!highlight_atoms || highlight_atoms->size() == mols.size(),
               "bad size");
  PRECONDITION(!highlight_bonds || highlight_bonds->size() == mols.size(),
               "bad size");
  PRECONDITION(
      !highlight_atom_maps || highlight_atom_maps->size() == mols.size(),
      "bad size");
  PRECONDITION(
      !highlight_bond_maps || highlight_bond_maps->size() == mols.size(),
      "bad size");
  PRECONDITION(!highlight_radii || highlight_radii->size() == mols.size(),
               "bad size");
  PRECONDITION(!confIds || confIds->size() == mols.size(), "bad size");
  PRECONDITION(panel_width_ != 0, "panel width cannot be zero");
  PRECONDITION(panel_height_ != 0, "panel height cannot be zero");
  PRECONDITION(width_ > 0 && height_ > 0,
               "drawMolecules() needs a fixed canvas size");
  if (mols.empty()) {
    return;
  }

  setupTextDrawer();
  int minScaleMol = 0;
  int minFontScaleMol = 0;
  int nCols = width() / panelWidth();
  int nRows = height() / panelHeight();

  for (size_t i = 0; i < mols.size(); ++i) {
    if (!mols[i]) {
      continue;
    }
    std::string legend = legends ? (*legends)[i] : "";
    const std::vector<int> *ha =
        highlight_atoms ? &(*highlight_atoms)[i] : nullptr;
    const std::map<int, DrawColour> *ham =
        highlight_atom_maps ? &(*highlight_atom_maps)[i] : nullptr;
    const std::map<int, DrawColour> *hbm =
        highlight_bond_maps ? &(*highlight_bond_maps)[i] : nullptr;
    int confId = confIds ? (*confIds)[i] : -1;
    const std::map<int, double> *hr =
        highlight_radii ? &(*highlight_radii)[i] : nullptr;
    unique_ptr<vector<int>> lhighlight_bonds;
    if (highlight_bonds) {
      lhighlight_bonds.reset(new std::vector<int>((*highlight_bonds)[i]));
    } else if (drawOptions().continuousHighlight && highlight_atoms) {
      lhighlight_bonds.reset(new vector<int>());
      MolDraw2D_detail::getBondHighlightsForAtoms(
          *mols[i], (*highlight_atoms)[i], *lhighlight_bonds);
    };
    auto prevSize = drawMols_.size();
    drawMols_.emplace_back(new MolDraw2D_detail::DrawMol(
        *mols[i], legend, panelWidth(), panelHeight(), drawOptions(),
        *text_drawer_, ha, lhighlight_bonds.get(), ham, hbm, nullptr, hr,
        supportsAnnotations(), confId));
    int row = 0;
    // note that this also works when no panel size is specified since
    // the panel dimensions defaults to -1
    if (nRows > 1) {
      row = i / nCols;
    }
    int col = 0;
    if (nCols > 1) {
      col = i % nCols;
    }
    drawMols_.back()->setOffsets(col * panelWidth(), row * panelHeight());
    drawMols_.back()->createDrawObjects();
    if (drawMols_.back()->getScale() < drawMols_[minScaleMol]->getScale()) {
      minScaleMol = prevSize;
    }
    if (drawMols_.back()->getFontScale() <
        drawMols_[minFontScaleMol]->getFontScale()) {
      minFontScaleMol = prevSize;
    }
  }

  for (auto &drawMol : drawMols_) {
    if (drawOptions().drawMolsSameScale) {
      drawMol->setScale(drawMols_[minScaleMol]->getScale(),
                        drawMols_[minFontScaleMol]->getFontScale(), true);
    }
    drawMol->tagAtomsWithCoords();
  }

  if (!drawMols_.empty()) {
    activeMolIdx_ = 0;
    startDrawing();
    for (size_t i = 0; i < drawMols_.size(); ++i) {
      activeMolIdx_ = i;
      drawTheMolecule((*drawMols_[i]));
    }
  } else {
    activeMolIdx_ = -1;
  }
}

// ****************************************************************************
void MolDraw2D::drawReaction(
    const ChemicalReaction &rxn, bool highlightByReactant,
    const std::vector<DrawColour> *highlightColorsReactants,
    const std::vector<int> *confIds) {
  std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> reagents;
  std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> products;
  std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> agents;

  // Copy the reaction because processing for drawing alters it.
  ChemicalReaction nrxn(rxn);
  // I think a larger minimum font size than the default works better for this
  int mfs = drawOptions().minFontSize > 12 ? drawOptions().minFontSize : 12;
  text_drawer_->setMinFontSize(mfs);
  int plusWidth;
  getReactionDrawMols(nrxn, highlightByReactant, highlightColorsReactants,
                      confIds, reagents, products, agents, plusWidth);

  if (height_ == -1) {
    for (auto &dm : drawMols_) {
      height_ = std::max(height_, dm->height_);
    }
  }
  std::vector<Point2D> offsets;
  Point2D arrowBeg, arrowEnd;
  calcReactionOffsets(reagents, products, agents, plusWidth, offsets, arrowBeg,
                      arrowEnd);
  activeMolIdx_ = -1;
  startDrawing();
  int xOffset = 0;
  xOffset = drawReactionPart(reagents, plusWidth, xOffset, offsets);
  // make sure the arrowhead is big enough
  double frac = 0.05;
  double delta = (arrowEnd - arrowBeg).length() * frac;
  // an arbitrary 5 pixel minimum
  if (delta < 5) {
    frac *= 5 / delta;
  }
  xOffset = drawReactionPart(agents, 0, xOffset, offsets);
  xOffset = drawReactionPart(products, plusWidth, xOffset, offsets);
  auto osbw = drawOptions().scaleBondWidth;
  if (reagents.empty() && products.empty() && agents.empty()) {
    // if it's an empty reaction, we won't have a DrawMol with a scale,
    // so we won't be able to do a scaled bond width
    drawOptions().scaleBondWidth = false;
  }
  drawArrow(arrowBeg, arrowEnd, false, frac, M_PI / 6,
            drawOptions().symbolColour, true);
  drawOptions().scaleBondWidth = osbw;
  if (drawOptions().includeMetadata) {
    this->updateMetadata(rxn);
  }
}

// ****************************************************************************
void MolDraw2D::drawLine(const Point2D &cds1, const Point2D &cds2,
                         const DrawColour &col1, const DrawColour &col2,
                         bool rawCoords) {
  if (drawOptions().comicMode) {
    // if rawCoords, we need a much bigger deviation.
    auto dev = rawCoords ? 0.5 : 0.03;
    auto scl = rawCoords ? 1.0 : scale_;
    setFillPolys(false);
    if (col1 == col2) {
      setColour(col1);
      auto pts =
          MolDraw2D_detail::handdrawnLine(cds1, cds2, scl, true, true, 4, dev);
      drawPolygon(pts, rawCoords);
    } else {
      auto mid = (cds1 + cds2) * 0.5;
      setColour(col1);
      auto pts =
          MolDraw2D_detail::handdrawnLine(cds1, mid, scl, true, false, 4, dev);
      drawPolygon(pts);
      setColour(col2);
      auto pts2 =
          MolDraw2D_detail::handdrawnLine(mid, cds2, scl, false, true, 4, dev);
      drawPolygon(pts2, rawCoords);
    }
  } else {
    if (col1 == col2) {
      setColour(col1);
      drawLine(cds1, cds2, rawCoords);
    } else {
      auto mid = (cds1 + cds2) * 0.5;
      setColour(col1);
      drawLine(cds1, mid, rawCoords);
      setColour(col2);
      drawLine(mid, cds2, rawCoords);
    }
  }
}

// ****************************************************************************
void MolDraw2D::drawTriangle(const Point2D &cds1, const Point2D &cds2,
                             const Point2D &cds3, bool rawCoords) {
  std::vector<Point2D> pts;
  if (!drawOptions().comicMode) {
    pts = {cds1, cds2, cds3};
  } else {
    double dev = rawCoords ? 0.5 : 0.03;
    double scl = rawCoords ? 1.0 : scale_;
    auto lpts =
        MolDraw2D_detail::handdrawnLine(cds1, cds2, scl, false, false, 4, dev);
    std::move(lpts.begin(), lpts.end(), std::back_inserter(pts));
    lpts = MolDraw2D_detail::handdrawnLine(cds2, cds3, scale_);
    std::move(lpts.begin(), lpts.end(), std::back_inserter(pts));
    lpts = MolDraw2D_detail::handdrawnLine(cds3, cds1, scale_);
    std::move(lpts.begin(), lpts.end(), std::back_inserter(pts));
  }
  drawPolygon(pts, rawCoords);
};

// ****************************************************************************
void MolDraw2D::drawEllipse(const Point2D &cds1, const Point2D &cds2,
                            bool rawCoords) {
  std::vector<Point2D> pts;
  MolDraw2D_detail::arcPoints(cds1, cds2, pts, 0, 360);
  drawPolygon(pts, rawCoords);
}

// ****************************************************************************
void MolDraw2D::drawArc(const Point2D &centre, double radius, double ang1,
                        double ang2, bool rawCoords) {
  drawArc(centre, radius, radius, ang1, ang2, rawCoords);
}

// ****************************************************************************
void MolDraw2D::drawArc(const Point2D &centre, double xradius, double yradius,
                        double ang1, double ang2, bool rawCoords) {
  std::vector<Point2D> pts;
  // 5 degree increments should be plenty, as the circles are probably
  // going to be small.
  int num_steps = 1 + int((ang2 - ang1) / 5.0);
  double ang_incr = double((ang2 - ang1) / num_steps) * M_PI / 180.0;
  double start_ang_rads = ang1 * M_PI / 180.0;
  for (int i = 0; i <= num_steps; ++i) {
    double ang = start_ang_rads + double(i) * ang_incr;
    double x = centre.x + xradius * cos(ang);
    double y = centre.y + yradius * sin(ang);
    pts.emplace_back(x, y);
  }

  if (fillPolys()) {
    // otherwise it draws an arc back to the pts.front() rather than filling
    // in the sector.
    pts.push_back(centre);
  }
  drawPolygon(pts, rawCoords);
}

// ****************************************************************************
void MolDraw2D::drawRect(const Point2D &cds1, const Point2D &cds2,
                         bool rawCoords) {
  std::vector<Point2D> pts(4);
  pts[0] = cds1;
  pts[1] = Point2D(cds1.x, cds2.y);
  pts[2] = cds2;
  pts[3] = Point2D(cds2.x, cds1.y);
  // if fillPolys() is false, it doesn't close the polygon because of
  // its use for drawing filled or open ellipse segments.
  if (!fillPolys()) {
    pts.push_back(cds1);
  }
  drawPolygon(pts, rawCoords);
}

// ****************************************************************************
//  we draw the line at cds2, perpendicular to the line cds1-cds2
void MolDraw2D::drawAttachmentLine(const Point2D &cds1, const Point2D &cds2,
                                   const DrawColour &col, double len,
                                   unsigned int nSegments, bool rawCoords) {
  auto perp = MolDraw2D_detail::calcPerpendicular(cds1, cds2);
  auto p1 = Point2D(cds2.x - perp.x * len / 2, cds2.y - perp.y * len / 2);
  auto p2 = Point2D(cds2.x + perp.x * len / 2, cds2.y + perp.y * len / 2);
  drawWavyLine(p1, p2, col, col, nSegments, rawCoords);
}

// ****************************************************************************
void MolDraw2D::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                             const DrawColour &col1, const DrawColour &col2,
                             unsigned int, double, bool rawCoords) {
  drawLine(cds1, cds2, col1, col2, rawCoords);
}

// ****************************************************************************
void MolDraw2D::drawArrow(const Point2D &arrowBegin, const Point2D &arrowEnd,
                          bool asPolygon, double frac, double angle,
                          const DrawColour &col, bool rawCoords) {
  Point2D ae(arrowEnd), p1, p2;
  MolDraw2D_detail::calcArrowHead(ae, p1, p2, arrowBegin, asPolygon, frac,
                                  angle);

  drawLine(arrowBegin, ae, col, col, rawCoords);
  if (!asPolygon) {
    drawLine(ae, p1, col, col, rawCoords);
    drawLine(ae, p2, col, col, rawCoords);
  } else {
    std::vector<Point2D> pts = {p1, ae, p2};
    bool fps = fillPolys();
    auto dc = colour();
    setFillPolys(true);
    setColour(col);
    drawPolygon(pts, rawCoords);
    setFillPolys(fps);
    setColour(dc);
  }
}

// ****************************************************************************
void MolDraw2D::drawPlus(const Point2D &plusPos, int plusWidth,
                         const DrawColour &col, bool rawCoords) {
  Point2D end1(plusPos.x, plusPos.y - plusWidth);
  Point2D end2(plusPos.x, plusPos.y + plusWidth);
  drawLine(end1, end2, col, col, rawCoords);
  end1 = Point2D(plusPos.x - plusWidth, plusPos.y);
  end2 = Point2D(plusPos.x + plusWidth, plusPos.y);
  drawLine(end1, end2, col, col, rawCoords);
}

// ****************************************************************************
// draws the string centred on cds
void MolDraw2D::drawString(const string &str, const Point2D &cds,
                           bool rawCoords) {
  auto draw_cds = rawCoords ? cds : getDrawCoords(cds);
  text_drawer_->drawString(str, draw_cds, MolDraw2D_detail::OrientType::N);
  //  int olw = lineWidth();
  //  setLineWidth(0);
  //  text_drawer_->drawStringRects(str, OrientType::N, TextAlignType::MIDDLE,
  //                                draw_cds, *this, rawCoords);
  //  setLineWidth(olw);
}

// ****************************************************************************
void MolDraw2D::drawString(const std::string &str, const Point2D &cds,
                           MolDraw2D_detail::TextAlignType talign,
                           bool rawCoords) {
  auto draw_cds = rawCoords ? cds : getDrawCoords(cds);
  text_drawer_->drawString(str, draw_cds, talign);
}

// ****************************************************************************
// transform a set of coords in the molecule's coordinate system
// to drawing system coordinates.  Prefers globalDrawTrans_ if it exists.
Point2D MolDraw2D::getDrawCoords(const Point2D &mol_cds) const {
  PRECONDITION(globalDrawTrans_ || !drawMols_.empty(), "no scaling info");
  if (globalDrawTrans_) {
    return globalDrawTrans_->getDrawCoords(mol_cds);
  } else {
    return drawMols_[activeMolIdx_]->getDrawCoords(mol_cds);
  }
}

// ****************************************************************************
Point2D MolDraw2D::getDrawCoords(int at_num) const {
  // this one can't use globalDrawTrans_, obviously.
  PRECONDITION(activeMolIdx_ >= 0 &&
                   static_cast<size_t>(activeMolIdx_) <= drawMols_.size(),
               "bad active mol index");
  PRECONDITION(!drawMols_.empty(), "no draw mols");
  PRECONDITION(
      static_cast<size_t>(at_num) < drawMols_[activeMolIdx_]->atCds_.size(),
      "bad atom number");
  return drawMols_[activeMolIdx_]->getDrawCoords(at_num);
}

// ****************************************************************************
Point2D MolDraw2D::getAtomCoords(const pair<int, int> &screen_cds) const {
  // Prefers globalDrawTrans_ if it exists.
  PRECONDITION(globalDrawTrans_ || !drawMols_.empty(), "no scaling info");
  return getAtomCoords(
      make_pair(double(screen_cds.first), double(screen_cds.second)));
}

// ****************************************************************************
// Prefers globalDrawTrans_ if it exists.
Point2D MolDraw2D::getAtomCoords(const pair<double, double> &screen_cds) const {
  // Prefers globalDrawTrans_ if it exists.
  PRECONDITION(globalDrawTrans_ || !drawMols_.empty(), "no scaling info");
  if (globalDrawTrans_) {
    return globalDrawTrans_->getAtomCoords(
        Point2D(screen_cds.first, screen_cds.second));
  } else {
    return drawMols_[activeMolIdx_]->getAtomCoords(
        Point2D(screen_cds.first, screen_cds.second));
  }
}

// ****************************************************************************
Point2D MolDraw2D::getAtomCoords(int at_num) const {
  PRECONDITION(!drawMols_.empty() && activeMolIdx_ >= 0 &&
                   static_cast<size_t>(activeMolIdx_) <= drawMols_.size(),
               "bad active mol index");
  return drawMols_[activeMolIdx_]->getAtomCoords(at_num);
}

// ****************************************************************************
const std::vector<Point2D> &MolDraw2D::atomCoords() const {
  PRECONDITION(activeMolIdx_ >= 0, "no index");
  return drawMols_[activeMolIdx_]->atCds_;
}

// ****************************************************************************
const std::vector<std::pair<std::string, MolDraw2D_detail::OrientType>>
    &MolDraw2D::atomSyms() const {
  PRECONDITION(activeMolIdx_ >= 0, "no index");
  return drawMols_[activeMolIdx_]->atomSyms_;
}

// ****************************************************************************
double MolDraw2D::getDrawLineWidth() const {
  double width = lineWidth();
  if (drawOptions().scaleBondWidth) {
    // lineWidthScaleFactor is defined in MolDraw2DHelpers.h
    width *= scale() * lineWidthScaleFactor;
    if (width < 0.0) {
      width = 0.0;
    }
  }
  return width;
}

// ****************************************************************************
Point2D MolDraw2D::minPt() const {
  // Prefers globalDrawTrans_ if it exists.
  PRECONDITION(globalDrawTrans_ || activeMolIdx_ >= 0, "bad active mol");
  // the ys are inverted in the DrawMol.
  if (globalDrawTrans_) {
    return Point2D(globalDrawTrans_->xMin_, -globalDrawTrans_->yMax_);
  } else {
    return Point2D(drawMols_[activeMolIdx_]->xMin_,
                   -drawMols_[activeMolIdx_]->yMax_);
  }
}
// ****************************************************************************
Point2D MolDraw2D::range() const {
  // Prefers globalDrawTrans_ if it exists.
  PRECONDITION(globalDrawTrans_ || activeMolIdx_ >= 0, "bad active mol");
  if (globalDrawTrans_) {
    return Point2D(globalDrawTrans_->xRange_, globalDrawTrans_->yRange_);
  } else {
    return Point2D(drawMols_[activeMolIdx_]->xRange_,
                   drawMols_[activeMolIdx_]->yRange_);
  }
}

// ****************************************************************************
double MolDraw2D::scale() const {
  PRECONDITION(activeMolIdx_ >= 0 &&
                   static_cast<size_t>(activeMolIdx_) <= drawMols_.size(),
               "bad active mol index");
  return drawMols_[activeMolIdx_]->getScale();
}

// ****************************************************************************
double MolDraw2D::fontSize() const { return text_drawer_->fontSize(); }

// ****************************************************************************
void MolDraw2D::setFontSize(double new_size) {
  text_drawer_->setFontSize(new_size);
}

// ****************************************************************************
void MolDraw2D::setScale(double newScale) { scale_ = newScale; }

// ****************************************************************************
void MolDraw2D::setScale(int width, int height, const Point2D &minv,
                         const Point2D &maxv, const ROMol *mol) {
  PRECONDITION(width > 0, "bad width");
  PRECONDITION(height > 0, "bad height");

  double x_min, x_max, x_range, y_min, y_max, y_range;
  bool setFontScale = false;
  if (mol) {
    setupTextDrawer();
    std::shared_ptr<MolDraw2D_detail::DrawMol> drawMol(
        new MolDraw2D_detail::DrawMol(*mol, "", panelWidth(), panelHeight(),
                                      drawOptions(), *text_drawer_));
    drawMol->createDrawObjects();
    // in the DrawMol, the ys are all inverted.
    x_min = min(minv.x, drawMol->xMin_);
    y_min = min(minv.y, -drawMol->yMax_);
    x_max = max(maxv.x, drawMol->xMax_);
    y_max = max(maxv.y, -drawMol->yMin_);
    fontScale_ = drawMol->getFontScale();
    setFontScale = true;
  } else {
    x_min = minv.x;
    y_min = minv.y;
    x_max = maxv.x;
    y_max = maxv.y;
  }
  x_range = x_max - x_min;
  y_range = y_max - y_min;

  if (x_range < 1.0e-4) {
    x_range = 1.0;
    x_min = -0.5;
  }
  if (y_range < 1.0e-4) {
    y_range = 1.0;
    y_min = -0.5;
  }

  // put a buffer round the drawing and calculate a final scale
  x_min -= drawOptions().padding * x_range;
  x_range *= 1 + 2 * drawOptions().padding;
  x_max = x_min + x_range;
  y_min -= drawOptions().padding * y_range;
  y_range *= 1 + 2 * drawOptions().padding;
  y_max = y_min + y_range;

  scale_ = std::min(double(width) / x_range, double(height) / y_range);
  // Absent any other information, we'll have to go with fontScale_ the
  // same as scale_.
  if (!setFontScale) {
    // get the text drawer to decide on a suitable fontScale.
    double ofs = text_drawer_->fontScale();
    text_drawer_->setFontScale(scale_);
    fontScale_ = text_drawer_->fontScale();
    text_drawer_->setFontScale(ofs, true);
  } else {
    fontScale_ *= scale_ / fontScale_;
  }
  forceScale_ = true;

  MolDraw2D_detail::DrawMol *drawMol = new MolDraw2D_detail::DrawMol(
      panelWidth(), panelHeight(), drawOptions(), *text_drawer_, x_min, x_max,
      y_min, y_max, scale_, fontScale_);
  Point2D trans, scale, toCentre;
  drawMol->getDrawTransformers(trans, scale, toCentre);
  globalDrawTrans_.reset(drawMol);
}

// ****************************************************************************
void MolDraw2D::getStringSize(const std::string &label, double &label_width,
                              double &label_height) const {
  text_drawer_->getStringSize(label, label_width, label_height);
  label_width /= scale();
  label_height /= scale();

  // cout << label << " : " << label_width << " by " << label_height
  //     << " : " << scale() << endl;
}

// ****************************************************************************
void MolDraw2D::getLabelSize(const string &label,
                             MolDraw2D_detail::OrientType orient,
                             double &label_width, double &label_height) const {
  if (orient == MolDraw2D_detail::OrientType::N ||
      orient == MolDraw2D_detail::OrientType::S) {
    label_height = 0.0;
    label_width = 0.0;
    vector<string> sym_bits =
        MolDraw2D_detail::atomLabelToPieces(label, orient);
    double height, width;
    for (auto bit : sym_bits) {
      getStringSize(bit, width, height);
      if (width > label_width) {
        label_width = width;
      }
      label_height += height;
    }
  } else {
    getStringSize(label, label_width, label_height);
  }
}

// ****************************************************************************
void MolDraw2D::getStringExtremes(const string &label,
                                  MolDraw2D_detail::OrientType orient,
                                  const Point2D &cds, double &x_min,
                                  double &y_min, double &x_max,
                                  double &y_max) const {
  text_drawer_->getStringExtremes(label, orient, x_min, y_min, x_max, y_max);
  auto draw_cds = getDrawCoords(cds);
  x_min += draw_cds.x;
  x_max += draw_cds.x;
  y_min += draw_cds.y;
  y_max += draw_cds.y;

  auto new_mins = getAtomCoords(make_pair(x_min, y_min));
  auto new_maxs = getAtomCoords(make_pair(x_max, y_max));
  x_min = new_mins.x;
  y_min = new_mins.y;
  x_max = new_maxs.x;
  y_max = new_maxs.y;

  // draw coords to atom coords reverses y
  if (y_min > y_max) {
    swap(y_min, y_max);
  }
}

// ****************************************************************************
void MolDraw2D::setActiveMolIdx(int newIdx) {
  PRECONDITION(newIdx >= -1 && newIdx < static_cast<int>(drawMols_.size()),
               "bad new activeMolIdx_");
  activeMolIdx_ = newIdx;
}

// ****************************************************************************
void MolDraw2D::setActiveAtmIdx(int at_idx1, int at_idx2) {
  at_idx1 = (at_idx1 < 0 ? -1 : at_idx1);
  at_idx2 = (at_idx2 < 0 ? -1 : at_idx2);
  if (at_idx2 >= 0 && at_idx1 < 0) {
    std::swap(at_idx1, at_idx2);
  }
  activeAtmIdx1_ = at_idx1;
  activeAtmIdx2_ = at_idx2;
}

// ****************************************************************************
void MolDraw2D::fixVariableDimensions(
    const MolDraw2D_detail::DrawMol &drawMol) {
  if (panel_width_ == -1) {
    width_ = panel_width_ = drawMol.width_;
  }
  if (panel_height_ == -1) {
    height_ = panel_height_ = drawMol.height_;
  }
}

// ****************************************************************************
void MolDraw2D::getReactionDrawMols(
    const ChemicalReaction &rxn, bool highlightByReactant,
    const std::vector<DrawColour> *highlightColorsReactants,
    const std::vector<int> *confIds,
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &reagents,
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &products,
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &agents,
    int &plusWidth) {
  ChemicalReaction nrxn(rxn);

  const double agentFrac = 0.4;
  double minScale = std::numeric_limits<double>::max(), minFontScale;

  std::map<int, DrawColour> atomColours;
  findReactionHighlights(rxn, highlightByReactant, highlightColorsReactants,
                         atomColours);
  // reactants & products
  makeReactionComponents(rxn.getReactants(), confIds, height(), atomColours,
                         reagents, minScale, minFontScale);
  makeReactionComponents(rxn.getProducts(), confIds, height(), atomColours,
                         products, minScale, minFontScale);
  for (auto &reagent : reagents) {
    reagent->setScale(minScale, minFontScale);
    reagent->shrinkToFit(false);
  }
  for (auto &product : products) {
    product->setScale(minScale, minFontScale);
    product->shrinkToFit(false);
  }
  // set the spacing for the plus signs to be 0.5 Angstrom.
  plusWidth = minScale;

  // agents
  int agentHeight = int(agentFrac * height_);
  minScale = std::numeric_limits<double>::max();
  makeReactionComponents(rxn.getAgents(), confIds, agentHeight, atomColours,
                         agents, minScale, minFontScale);
  for (auto &agent : agents) {
    agent->setScale(minScale, minFontScale);
    agent->shrinkToFit(false);
  }

  // set the active atom/bond indices so they run in series for all pieces.
  // DrawMols start each new set at 0 by default.  The original code had
  // reagents, products then agents, so do the same here.
  int atomIdxOffset = 0, bondIdxOffset = 0;
  for (auto &dm : drawMols_) {
    dm->activeAtmIdxOffset_ = atomIdxOffset;
    dm->activeBndIdxOffset_ = bondIdxOffset;
    atomIdxOffset += dm->drawMol_->getNumAtoms();
    bondIdxOffset += dm->drawMol_->getNumBonds();
  }
}

// ****************************************************************************
void MolDraw2D::makeReactionComponents(
    const std::vector<RDKit::ROMOL_SPTR> &bits, const std::vector<int> *confIds,
    int heightToUse, std::map<int, DrawColour> &atomColours,
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &dms,
    double &minScale, double &minFontScale) {
  for (size_t midx = 0; midx < bits.size(); ++midx) {
    std::vector<int> highlightAtoms, highlightBonds;
    std::map<int, DrawColour> highlightAtomMap, highlightBondMap;
    int cid = -1;
    if (confIds) {
      cid = (*confIds)[midx];
    }
    auto fragMol = bits[midx].get();
    for (auto atom : fragMol->atoms()) {
      auto ai = atomColours.find(atom->getAtomMapNum());
      for (auto bond : fragMol->bonds()) {
        auto beg = bond->getBeginAtom();
        auto begi = atomColours.find(beg->getAtomMapNum());
        if (begi != atomColours.end()) {
          auto end = bond->getEndAtom();
          auto endi = atomColours.find(end->getAtomMapNum());
          if (endi != atomColours.end() && begi->second == endi->second) {
            highlightBonds.push_back(bond->getIdx());
            highlightBondMap.insert(
                std::make_pair(bond->getIdx(), begi->second));
          }
        }
      }
      if (ai != atomColours.end()) {
        highlightAtoms.push_back(atom->getIdx());
        highlightAtomMap.insert(std::make_pair(atom->getIdx(), ai->second));
        atom->setAtomMapNum(0);
      }
    }
    makeReactionDrawMol(*(RWMol *)bits[midx].get(), cid, heightToUse,
                        highlightAtoms, highlightBonds, highlightAtomMap,
                        highlightBondMap, dms);
    if (dms.back()->getScale() < minScale) {
      minScale = dms.back()->getScale();
      minFontScale = dms.back()->getFontScale();
    }
  }
}

// ****************************************************************************
void MolDraw2D::makeReactionDrawMol(
    RWMol &mol, int confId, int molHeight,
    const std::vector<int> &highlightAtoms,
    const std::vector<int> &highlightBonds,
    const std::map<int, DrawColour> &highlightAtomMap,
    const std::map<int, DrawColour> &highlightBondMap,
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &mols) {
  mol.updatePropertyCache(false);
  if (drawOptions().prepareMolsBeforeDrawing) {
    RDLog::LogStateSetter blocker;
    MolOps::KekulizeIfPossible(mol, false);
    MolOps::setHybridization(mol);
  }
  if (!mol.getNumConformers()) {
    const bool canonOrient = true;
    RDDepict::compute2DCoords(mol, nullptr, canonOrient);
  } else {
    // we need to center the molecule
    MolDraw2D_detail::centerMolForDrawing(mol, confId);
  }
  // when preparing a reaction component to be drawn we should neither kekulize
  // (we did that above if required) nor add chiralHs
  const bool kekulize = false;
  const bool addChiralHs = false;
  MolDraw2DUtils::prepareMolForDrawing(mol, kekulize, addChiralHs);
  // the height is fixed, but the width is allowed to be as large as the
  // height and molecule dimensions dictate.
  mols.emplace_back(new MolDraw2D_detail::DrawMol(
      mol, "", -1, molHeight, drawOptions(), *text_drawer_, &highlightAtoms,
      &highlightBonds, &highlightAtomMap, &highlightBondMap, nullptr, nullptr,
      supportsAnnotations(), confId, true));
  mols.back()->createDrawObjects();
  drawMols_.push_back(mols.back());
  ++activeMolIdx_;
}

// ****************************************************************************
void MolDraw2D::calcReactionOffsets(
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &reagents,
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &products,
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &agents,
    int &plusWidth, std::vector<Point2D> &offsets, Point2D &arrowBeg,
    Point2D &arrowEnd) {
  // calculate the total width of the drawing - it may need re-scaling if
  // it's too wide for the panel.
  const int arrowMult = 2;  // number of plusWidths for an empty arrow.
  auto reactionWidth = [&](int gapWidth) -> int {
    int totWidth = 0;
    if (!reagents.empty()) {
      for (auto &dm : reagents) {
        totWidth += dm->width_;
      }
      totWidth += gapWidth * (reagents.size() - 1);
    }
    if (agents.empty()) {
      totWidth += arrowMult * gapWidth;
    } else {
      // the agent doesn't start at front of arrow
      for (auto &dm : agents) {
        totWidth += dm->width_ + gapWidth;
      }
      totWidth += gapWidth * (agents.size() - 1) / 2;
    }
    totWidth += gapWidth;  // either side of arrow
    if (!products.empty()) {
      for (auto &dm : products) {
        totWidth += dm->width_;
      }
      // we don't want a plus after the last product
      totWidth += gapWidth * (products.size() - 1);
    }
    return totWidth;
  };

  auto scaleDrawMols =
      [&](std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &dms,
          double stretch) {
        for (auto &dm : dms) {
          dm->setScale(stretch * dm->getScale(), stretch * dm->getFontScale(),
                       false);
          dm->shrinkToFit(false);
        }
      };

  if (width_ == -1) {
    width_ = reactionWidth(plusWidth);
  }

  // The DrawMols are sized/scaled according to the height() which is all
  // there is to go on initially, so they may be wider than the width() in
  // total. If so, shrink them to fit.  Because the shrinkng imposes min and max
  // font sizes, we may not get the smaller size we want first go, so iterate
  // until we do or we give up.
  int totWidth = reactionWidth(plusWidth);
  for (int i = 0; i < 5; ++i) {
    auto maxWidthIt =
        std::max_element(drawMols_.begin(), drawMols_.end(),
                         [&](std::shared_ptr<MolDraw2D_detail::DrawMol> &lhs,
                             std::shared_ptr<MolDraw2D_detail::DrawMol> &rhs) {
                           return lhs->width_ < rhs->width_;
                         });
    plusWidth = (*maxWidthIt)->width_ / 4;
    plusWidth = plusWidth > width() / 20 ? width() / 20 : plusWidth;
    auto oldTotWidth = totWidth;
    auto stretch = double(width_ * (1 - drawOptions().padding)) / totWidth;
    // If stretch < 1, we need to shrink the DrawMols to fit.  This isn't
    // necessary if we're just stretching them along the panel as they already
    // fit for height.
    if (stretch < 1.0) {
      scaleDrawMols(reagents, stretch);
      scaleDrawMols(agents, stretch);
      scaleDrawMols(products, stretch);
    } else {
      break;
    }
    totWidth = reactionWidth(plusWidth);
    if (fabs(totWidth - oldTotWidth) < 0.01 * width()) {
      break;
    }
  }

  // make sure plusWidth remains 1A, in pixels.
  if (!reagents.empty()) {
    plusWidth = reagents.front()->scale_;
  } else if (!products.empty()) {
    plusWidth = products.front()->scale_;
  }
  // if there's space, we can afford make the extras a bit bigger.
  int numGaps = reagents.empty() ? 0 : reagents.size() - 1;
  numGaps += products.empty() ? 0 : products.size() - 1;
  numGaps += 1;  // for either side of the arrow
  numGaps += agents.empty() ? 2 : agents.size() - 1;
  if (width() - totWidth > numGaps * 5) {
    plusWidth += 5;
  }
  totWidth = reactionWidth(plusWidth);

  // And finally work out where to put all the pieces, centring them.
  int xOffset = (width() - totWidth) / 2;
  for (size_t i = 0; i < reagents.size(); ++i) {
    offsets.emplace_back(xOffset, (height() - reagents[i]->height_) / 2);
    xOffset += reagents[i]->width_ + plusWidth;
  }
  if (reagents.empty()) {
    xOffset += plusWidth / 2;
  } else {
    // only half a plusWidth to the arrow
    xOffset -= plusWidth / 2;
  }
  arrowBeg.y = height() / 2;
  arrowBeg.x = xOffset;
  if (agents.empty()) {
    arrowEnd = Point2D(arrowBeg.x + arrowMult * plusWidth, height() / 2);
  } else {
    xOffset += plusWidth / 2;
    for (size_t i = 0; i < agents.size(); ++i) {
      offsets.emplace_back(xOffset, 0.45 * height() - agents[i]->height_);
      xOffset += agents[i]->width_ + plusWidth / 2;
    }
    // the overlap at the end of the arrow has already been added in the loop
    arrowEnd = Point2D(xOffset, height() / 2);
  }
  xOffset = arrowEnd.x + plusWidth / 2;
  for (size_t i = 0; i < products.size(); ++i) {
    offsets.emplace_back(xOffset, (height() - products[i]->height_) / 2);
    xOffset += products[i]->width_ + plusWidth;
  }
}

// ****************************************************************************
int MolDraw2D::drawReactionPart(
    std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &reactBit,
    int plusWidth, int initOffset, const std::vector<Point2D> &offsets) {
  if (reactBit.empty()) {
    return initOffset;
  }

  Point2D plusPos(0.0, height() / 2);
  for (size_t i = 0; i < reactBit.size(); ++i) {
    ++activeMolIdx_;
    reactBit[i]->setOffsets(offsets[initOffset].x, offsets[initOffset].y);
    reactBit[i]->draw(*this);
#if 0
// this is convenient for debugging
    setColour(DrawColour(0, 1.0, 1.0));
    drawLine(Point2D(offsets[initOffset].x + reactBit[i]->width_, 0),
             Point2D(offsets[initOffset].x + reactBit[i]->width_, height_),
             true);
    drawLine(Point2D(offsets[initOffset].x, 0),
             Point2D(offsets[initOffset].x, height_), true);
    setColour(DrawColour(1.0, 0, 1.0));
    drawLine(Point2D(offsets[initOffset].x, offsets[initOffset].y),
             Point2D(offsets[initOffset].x + reactBit[i]->width_,
                     offsets[initOffset].y),
             true);
    drawLine(Point2D(offsets[initOffset].x,
                     offsets[initOffset].y + reactBit[i]->height_),
             Point2D(offsets[initOffset].x + reactBit[i]->width_,
                     offsets[initOffset].y + reactBit[i]->height_),
             true);
#endif
    if (plusWidth && i < reactBit.size() - 1) {
      plusPos.x = (offsets[initOffset].x + reactBit[i]->width_ +
                   offsets[initOffset + 1].x) /
                  2;
      int plusStroke = plusWidth > 30 ? 10 : plusWidth / 3;
      drawPlus(plusPos, plusStroke, drawOptions().symbolColour, true);
    }
    ++initOffset;
  }
  return initOffset;
}

// ****************************************************************************
void MolDraw2D::findReactionHighlights(
    const ChemicalReaction &rxn, bool highlightByReactant,
    const std::vector<DrawColour> *highlightColorsReactants,
    std::map<int, DrawColour> &atomColours) const {
  std::unique_ptr<ROMol> tmol(ChemicalReactionToRxnMol(rxn));
  if (highlightByReactant) {
    const auto *colors = &drawOptions().highlightColourPalette;
    if (highlightColorsReactants) {
      colors = highlightColorsReactants;
    }
    for (size_t midx = 0; midx < rxn.getReactants().size(); ++midx) {
      auto fragMol = rxn.getReactants()[midx].get();
      for (auto &atom : fragMol->atoms()) {
        int atomRole = -1;
        if (atom->getPropIfPresent("molRxnRole", atomRole) && atomRole == 1 &&
            atom->getAtomMapNum()) {
          atomColours.insert(std::make_pair(atom->getAtomMapNum(),
                                            (*colors)[midx % colors->size()]));
        }
      }
    }
  }
}

// ****************************************************************************
void MolDraw2D::startDrawing() {
  if (needs_init_) {
    initDrawing();
    needs_init_ = false;
  }
  if (activeMolIdx_ <= 0 && drawOptions().clearBackground) {
    clearDrawing();
  }
}

// ****************************************************************************
void MolDraw2D::drawTheMolecule(MolDraw2D_detail::DrawMol &drawMol) {
  if (globalDrawTrans_) {
    drawMol.setTransformation(*globalDrawTrans_);
  } else if (forceScale_) {
    drawMol.setScale(scale_, fontScale_);
  }
  drawMol.draw(*this);
  if (drawOptions().includeMetadata) {
    this->updateMetadata(*drawMol.drawMol_, drawMol.confId_);
  }
}

// ****************************************************************************
void MolDraw2D::setupTextDrawer() {
  text_drawer_->setMaxFontSize(drawOptions().maxFontSize);
  text_drawer_->setMinFontSize(drawOptions().minFontSize);
  if (drawOptions().baseFontSize > 0.0) {
    text_drawer_->setBaseFontSize(drawOptions().baseFontSize);
  }
  try {
    text_drawer_->setFontFile(drawOptions().fontFile);
  } catch (std::runtime_error &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
    text_drawer_->setFontFile("");
    BOOST_LOG(rdWarningLog) << "Falling back to original font file "
                            << text_drawer_->getFontFile() << "." << std::endl;
  }
}

}  // namespace RDKit
