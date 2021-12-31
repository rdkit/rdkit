//
//  Copyright (C) 2014-2021 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <iostream>
#include <limits>

#include <Geometry/Transform2D.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolDraw2D/AtomSymbol.h>
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawShape.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

namespace RDKit {

// ****************************************************************************
DrawMol::DrawMol(const ROMol &mol, const std::string &legend,
                 int width, int height,
                 MolDrawOptions &drawOptions, DrawText &textDrawer,
                 const std::vector<int> *highlight_atoms,
                 const std::vector<int> *highlight_bonds,
                 const std::map<int, DrawColour> *highlight_atom_map,
                 const std::map<int, DrawColour> *highlight_bond_map,
                 const std::vector<std::pair<DrawColour, DrawColour>> *bond_colours,
                 const std::map<int, double> *highlight_radii,
                 int confId)
    : legend_(legend),
      drawOptions_(drawOptions),
      textDrawer_(textDrawer),
      highlightAtoms_(highlight_atoms),
      highlightBonds_(highlight_bonds),
      highlightAtomMap_(highlight_atom_map),
      highlightBondMap_(highlight_bond_map),
      bondColours_(bond_colours),
      highlightRadii_(highlight_radii),
      confId_(confId),
      width_(width),
      height_(height),
      scale_(1.0),
      fontScale_(1.0),
      xMin_(std::numeric_limits<double>::max() / 2.0),
      yMin_(std::numeric_limits<double>::max() / 2.0),
      xMax_(-std::numeric_limits<double>::max() / 2.0),
      yMax_(-std::numeric_limits<double>::max() / 2.0),
      xRange_(std::numeric_limits<double>::max()),
      yRange_(std::numeric_limits<double>::max()) {
  textDrawer_.setFontScale(fontScale_, true);
  std::cout << "Top of DrawMol c'tor" << std::endl;
  std::cout << "Width = " << width_ << "  height = " << height_ << std::endl;
  std::cout << "initial font scale : " << textDrawer_.fontScale() << std::endl;
  initDrawMolecule(mol);
}

// ****************************************************************************
void DrawMol::createDrawObjects() {
  extractAll();
  calculateScale();

  std::cout << "CCCCCCCCCCCCC" << std::endl;
  if (!textDrawer_.setFontScale(fontScale_, false)) {
    double nfs = textDrawer_.fontScale();
    std::cout << "font scale hit extreme.  We got " << nfs << " would have been "
              << fontScale_ << std::endl;
    textDrawer_.setFontScale(nfs / fontScale_, true);
    resetEverything();
    fontScale_ = textDrawer_.fontScale();
    extractAll();
    calculateScale();
    textDrawer_.setFontScale(fontScale_);
  } else {
    std::cout << "fontScale was fine" << std::endl;
  }
  std::cout << "XX font scale : " << textDrawer_.fontScale()
            << " : " << fontScale_ << std::endl;
  changeToDrawCoords();
  drawingInitialised_ = true;
}

// ****************************************************************************
void DrawMol::initDrawMolecule(const ROMol &mol) {
  drawMol_.reset(new RWMol(mol));
  if (drawOptions_.prepareMolsBeforeDrawing || !mol.getNumConformers()) {
    MolDraw2DUtils::prepareMolForDrawing(*drawMol_);
  }
  if (drawOptions_.centreMoleculesBeforeDrawing) {
    if (drawMol_->getNumConformers()) {
      centerMolForDrawing(*drawMol_, confId_);
    }
  }
  if (drawOptions_.simplifiedStereoGroupLabel &&
      !mol.hasProp(common_properties::molNote)) {
    prepareStereoGroups(*drawMol_);
  }
  if (drawOptions_.addStereoAnnotation) {
    MolDraw2D_detail::addStereoAnnotation(*drawMol_);
  }
  if (drawOptions_.addAtomIndices) {
    MolDraw2D_detail::addAtomIndices(*drawMol_);
  }
  if (drawOptions_.addBondIndices) {
    MolDraw2D_detail::addBondIndices(*drawMol_);
  }
}

// ****************************************************************************
void DrawMol::extractAll() {
  extractAtomCoords();
  extractAtomSymbols();
  extractBonds();
  extractRegions();
  extractHighlights();
  extractMolNotes();
  extractAtomNotes();
  extractBondNotes();
  extractRadicals();
  extractLegend();
}

// ****************************************************************************
void DrawMol::extractAtomCoords() {
  PRECONDITION(static_cast<int>(drawMol_->getNumConformers()) > 0, "no coords");

  const RDGeom::POINT3D_VECT &locs =
      drawMol_->getConformer(confId_).getPositions();

  // the transformation rotates anti-clockwise, as is conventional, but
  // probably not what our user expects.
  double rot = -drawOptions_.rotate * M_PI / 180.0;
  // assuming that if drawOptions_.rotate is set to 0.0, rot will be
  // exactly 0.0 without worrying about floating point number dust.  Does
  // anyone know if this is true?  It's not the end of the world if not,
  // as it's just an extra largely pointless rotation.
  // Floating point numbers are like piles of sand; every time you move
  // them around, you lose a little sand and pick up a little dirt.
  // — Brian Kernighan and P.J. Plauger
  // Nothing brings fear to my heart more than a floating point number.
  // — Gerald Jay Sussman
  // Some developers, when encountering a problem, say: “I know, I’ll
  // use floating-point numbers!”   Now, they have 1.9999999997 problems.
  // — unknown
  RDGeom::Transform2D trans;
  trans.SetTransform(Point2D(0.0, 0.0), rot);
  atCds_.clear();
  for (auto this_at : drawMol_->atoms()) {
    int thisIdx = this_at->getIdx();
    Point2D pt(locs[thisIdx].x, -locs[thisIdx].y);
    if (rot != 0.0) {
      trans.TransformPoint(pt);
    }
    atCds_.emplace_back(pt);
    // std::cout << "coords for " << thisIdx << " : " << pt << std::endl;
  }
}

// ****************************************************************************
void DrawMol::extractAtomSymbols() {
  PRECONDITION(atCds_.size() > 0, "no coords");

  atomicNums_.clear();
  for (auto at1 : drawMol_->atoms()) {
    std::pair<std::string, OrientType> atSym = getAtomSymbolAndOrientation(*at1);
    atomSyms_.emplace_back(atSym);
    if (!isComplexQuery(at1)) {
      atomicNums_.emplace_back(at1->getAtomicNum());
    } else {
      atomicNums_.push_back(0);
    }
    if (!atSym.first.empty()) {
      DrawColour atCol = getColour(at1->getIdx(), drawOptions_, atomicNums_,
                                   highlightAtoms_, highlightAtomMap_);
      AtomSymbol *al = new AtomSymbol(atSym.first, at1->getIdx(), atSym.second,
                                    atCds_[at1->getIdx()], atCol, textDrawer_);
      atomLabels_.emplace_back(std::unique_ptr<AtomSymbol>(al));
    } else {
      atomLabels_.emplace_back(std::unique_ptr<AtomSymbol>());
    }
  }
}

// ****************************************************************************
void DrawMol::extractBonds() {
  double doubleBondOffset = drawOptions_.multipleBondOffset;
  // mol files from, for example, Marvin use a bond length of 1 for just about
  // everything. When this is the case, the default multipleBondOffset is just
  // too much, so scale it back.
  calcMeanBondLengthSquare();
  if (meanBondLengthSquare_ < 1.4) {
    doubleBondOffset *= 0.6;
  }

  for (auto bond : drawMol_->bonds()) {
    bool isComplex = false;
    if (bond->hasQuery()) {
      std::string descr = bond->getQuery()->getDescription();
      if (bond->getQuery()->getNegation() || descr != "BondOrder") {
        isComplex = true;
        makeQueryBond(bond, doubleBondOffset);
      }
    }

    if (!isComplex) {
      makeStandardBond(bond, doubleBondOffset);
    }
  }

}

// ****************************************************************************
void DrawMol::extractHighlights() {
  std::cout << "DrawMol::extractHighlights" << std::endl;
  if (drawOptions_.continuousHighlight) {
    makeContinuousHighlights();
  } else {
    if (drawOptions_.circleAtoms && highlightAtoms_) {
      makeAtomCircleHighlights();
    }
  }
}

// ****************************************************************************
void DrawMol::extractRegions() {
  if (!drawOptions_.atomRegions.empty()) {
    for (auto &region : drawOptions_.atomRegions) {
      if (region.size() > 1) {
        Point2D minv = atCds_[region[0]];
        Point2D maxv = atCds_[region[0]];
        for (int idx : region) {
          const Point2D &pt = atCds_[idx];
          minv.x = std::min(minv.x, pt.x);
          minv.y = std::min(minv.y, pt.y);
          maxv.x = std::max(maxv.x, pt.x);
          maxv.y = std::max(maxv.y, pt.y);
        }
        Point2D center = (maxv + minv) / 2;
        Point2D size = (maxv - minv);
        size *= 0.2;
        minv -= size / 2;
        maxv += size / 2;
        std::vector<Point2D> pts(4);
        pts[0] = minv;
        pts[1] = Point2D(minv.x, maxv.y);
        pts[2] = maxv;
        pts[3] = Point2D(maxv.x, minv.y);
        DrawColour col(0.8, 0.8, 0.8);
        DrawShape *pl = new DrawShapePolyLine(pts, 1, false, col, true);
        highlights_.emplace_back(std::unique_ptr<DrawShape>(pl));
      }
    }
  }
}

// ****************************************************************************
void DrawMol::extractMolNotes() {
  std::string note;
  // the molNote property takes priority
  if (!drawMol_->getPropIfPresent(common_properties::molNote, note)) {
    unsigned int chiralFlag;
    if (drawOptions_.includeChiralFlagLabel &&
        drawMol_->getPropIfPresent(common_properties::_MolFileChiralFlag,
                             chiralFlag) &&
        chiralFlag) {
      note = "ABS";
    }
  }

  if (!note.empty()) {
    AnnotationType annot;
    annot.text_ = note;
    annot.align_ = TextAlignType::START;
    annot.scaleText_ = false;
    calcMolNotePosition(atCds_, textDrawer_, annot);
    if (annot.rect_.width_ < 0.0) {
      BOOST_LOG(rdWarningLog)
          << "Couldn't find good place for molecule note " << note << std::endl;
    } else {
      annotations_.push_back(annot);
    }
  }
}

// ****************************************************************************
void DrawMol::extractAtomNotes() {
  std::cout << "extractAtomNotes" << std::endl;
  for (auto atom : drawMol_->atoms()) {
    std::string note;
    if (atom->getPropIfPresent(common_properties::atomNote, note)) {
      if (!note.empty()) {
        std::cout << "atom " << atom->getIdx() << " : " << note << std::endl;
        AnnotationType annot;
        annot.text_ = note;
        calcAnnotationPosition(atom, annot);
        if (annot.rect_.width_ < 0.0) {
          BOOST_LOG(rdWarningLog)
              << "Couldn't find good place for note " << note << " for atom "
              << atom->getIdx() << std::endl;
        } else {
          annotations_.push_back(annot);
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMol::extractBondNotes() {
  std::cout << "extractBondNotes" << std::endl;
  for (auto bond : drawMol_->bonds()) {
    std::string note;
    if (bond->getPropIfPresent(common_properties::bondNote, note)) {
      if (!note.empty()) {
        AnnotationType annot;
        annot.text_ = note;
        calcAnnotationPosition(bond, annot);
        if (annot.rect_.width_ < 0.0) {
          BOOST_LOG(rdWarningLog)
              << "Couldn't find good place for note " << note << " for bond "
              << bond->getIdx() << std::endl;
        } else {
          annotations_.push_back(annot);
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMol::extractRadicals() {
  for (auto atom : drawMol_->atoms()) {
    if (!atom->getNumRadicalElectrons()) {
      continue;
    }
    StringRect rad_rect;
    OrientType orient = calcRadicalRect(atom, rad_rect);
    radicals_.emplace_back(std::make_pair(rad_rect, orient));
  }
}

// ****************************************************************************
void DrawMol::calculateScale() {
  findExtremes();

  // if either width_ or height_ is < 0 we are going to make a picture of
  // as yet unknown size, with fixed scale.
  bool setWidth = false;
  if (width_ < 0) {
    // FIX: technically we need to take the legend width into account too!
    width_ = drawOptions_.scalingFactor * xRange_;
    setWidth = true;
  }
  bool setHeight = false;
  double thisYRange = yRange_;
  if (height_ < 0) {
    // we need to adjust the range for the legend
    // if it's not present then legend_height_ will be zero and this will be a
    // no-op
    thisYRange += legendHeight_ / drawOptions_.scalingFactor;
    height_ = drawOptions_.scalingFactor * thisYRange;
    setHeight = true;
  }

  // put a 5% buffer round the drawing and calculate a final scale
  xMin_ -= drawOptions_.padding * xRange_;
  xRange_ *= 1 + 2 * drawOptions_.padding;
  xMax_ = xMin_ + xRange_;
  yMin_ -= drawOptions_.padding * yRange_;
  yRange_ *= 1 + 2 * drawOptions_.padding;
  yMax_ = yMin_ + yRange_;

  double newScale = 1.0;
  if (xRange_ > 1e-4 || yRange_ > 1e-4) {
    if (setWidth) {
      width_ = drawOptions_.scalingFactor * xRange_;
    }
    if (setHeight) {
      height_ = drawOptions_.scalingFactor * thisYRange;
    }

    newScale = std::min(double(width_) / xRange_,
                      double(height_ - legendHeight_) / yRange_);
    double fix_scale = newScale;
    // after all that, use the fixed scale unless it's too big, in which case
    // scale the drawing down to fit.
    // fixedScale takes precedence if both it and fixedBondLength are given.
    if (drawOptions_.fixedBondLength > 0.0) {
      fix_scale = drawOptions_.fixedBondLength;
    }
    if (drawOptions_.fixedScale > 0.0) {
      fix_scale = double(width_) * drawOptions_.fixedScale;
    }
    if (newScale > fix_scale) {
      newScale = fix_scale;
    }
  }
  double scale_mult = newScale / scale_;
  scale_ *= scale_mult;
  fontScale_ *= scale_mult;
  std::cout << "Final Scale : " << scale_ << " and fontScale_ : " << fontScale_ << std::endl;
  std::cout << "Padded mins : " << xMin_ << ", " << yMin_ << " with ranges : "
            << xRange_ << ", " << yRange_ << std::endl;
}

// ****************************************************************************
void DrawMol::findExtremes() {
  for (auto &bond : bonds_) {
    bond->findExtremes(xMin_, xMax_, yMin_, yMax_);
  }
  for (auto &atLab : atomLabels_) {
    if (atLab) {
      atLab->findExtremes(xMin_, xMax_, yMin_, yMax_);
    }
  }
  for (auto &hl : highlights_) {
    hl->findExtremes(xMin_, xMax_, yMin_, yMax_);
  }
  findAnnotationExtremes(annotations_, xMin_, xMax_, yMin_, yMax_);
  findRadicalExtremes(radicals_, xMin_, xMax_, yMin_, yMax_);

  // calculate the x and y spans
  xRange_ = xMax_ - xMin_;
  yRange_ = yMax_ - yMin_;
  if (xRange_ < 1e-4) {
    xRange_ = 2.0;
    xMin_ -= 1.0;
    xMax_ += 1.0;
  }
  if (yRange_ < 1e-4) {
    yRange_ = 2.0;
    yMin_ -= 1.0;
    yMax_ += 1.0;
  }
  std::cout << "Final mins : " << xMin_ << ", " << yMin_ << " with ranges : "
            << xRange_ << ", " << yRange_ << std::endl;
}

// ****************************************************************************
void DrawMol::changeToDrawCoords() {
  std::cout << "scaled mins and ranges : " << xMin_ * scale_ << ", "
            << yMin_ * scale_ << " :: " << xRange_ * scale_ << ", "
            << yRange_ * scale_ << std::endl;
  std::cout << "scales : " << scale_ << " and " << fontScale_ << std::endl;
  Point2D trans(-xMin_, -yMin_);
  Point2D scale(scale_, scale_);
  Point2D scaledRanges(scale_ * xRange_, scale_ * yRange_);
  int drawHeight = height_ - legendHeight_;
  Point2D toCentre((width_ - scaledRanges.x) / 2.0,
                   (drawHeight - scaledRanges.y) / 2.0);
  for (auto &bond : bonds_) {
    bond->move(trans);
    bond->scale(scale);
    bond->move(toCentre);
  }
  for (auto &label : atomLabels_) {
    if (label) {
      label->move(trans);
      label->scale(scale);
      label->move(toCentre);
    }
  }
  for (auto &hl : highlights_) {
    hl->move(trans);
    hl->scale(scale);
    hl->move(toCentre);
  }
  for (auto &annot : annotations_) {
    annot.rect_.trans_ += trans;
    annot.rect_.trans_.x *= scale.x;
    annot.rect_.trans_.y *= scale.y;
    annot.rect_.trans_ += toCentre;
    annot.rect_.width_ *= scale.x;
    annot.rect_.height_ *= scale.y;
  }
  for (auto &rad : radicals_) {
    rad.first.trans_ += trans;
    rad.first.trans_.x *= scale.x;
    rad.first.trans_.y *= scale.y;
    rad.first.trans_ += toCentre;
  }
}

// ****************************************************************************
void DrawMol::draw(MolDraw2D &drawer) const {
  PRECONDITION(drawingInitialised_,
               "you must call createDrawingObjects before calling draw")
  auto keepScale = drawer.scale();
  drawer.setScale(scale_);
  for (auto &hl : highlights_) {
    hl->draw(drawer);
  }
  for (auto &bond : bonds_) {
    bond->draw(drawer);
  }
  for (auto &label : atomLabels_) {
    if (label) {
      label->draw(drawer);
    }
  }
  drawAllAnnotations(drawer);
  drawRadicals(drawer);
  drawLegend(drawer);
  drawer.setScale(scale_);
}

// ****************************************************************************
void DrawMol::drawAllAnnotations(MolDraw2D &drawer) const {
  std::string currActClass = drawer.getActiveClass();
  if (currActClass.empty()) {
    drawer.setActiveClass("note");
  } else {
    drawer.setActiveClass(currActClass + " note");
  }
  for (auto &annot : annotations_) {
    drawAnnotation(annot);
  }
  drawer.setActiveClass(currActClass);
}

// ****************************************************************************
void DrawMol::drawAnnotation(const AnnotationType &annot) const {
  double full_font_scale = textDrawer_.fontScale();
  // turn off minFontSize for the annotation, as we do want it to be smaller
  // than the letters, even if that makes it tiny.  The annotation positions
  // have been calculated on the assumption that this is the case, and if
  // minFontSize is applied, they may well clash with the atom symbols.
  if (annot.scaleText_) {
    textDrawer_.setFontScale(
        drawOptions_.annotationFontScale * full_font_scale, true);
  }
  textDrawer_.setColour(annot.col_);
  textDrawer_.drawString(annot.text_, annot.rect_.trans_, annot.align_);
  if (annot.scaleText_) {
    textDrawer_.setFontScale(full_font_scale, true);
  }
}

// ****************************************************************************
void DrawMol::drawLegend(MolDraw2D &drawer) const {
  textDrawer_.setColour(drawOptions_.legendColour);
  double o_font_scale = textDrawer_.fontScale();
  double fsize = textDrawer_.fontSize();
  double new_font_scale = o_font_scale * drawOptions_.legendFontSize / fsize;
  textDrawer_.setFontScale(new_font_scale, true);
  std::string currActClass = drawer.getActiveClass();
  if (currActClass.empty()) {
    drawer.setActiveClass("legend");
  } else {
    drawer.setActiveClass(currActClass + " legend");
  }
  for (auto &leg : legends_) {
    drawAnnotation(leg);
  }
  drawer.setActiveClass(currActClass);
  textDrawer_.setFontScale(o_font_scale, true);
}

// ****************************************************************************
void DrawMol::drawRadicals(MolDraw2D &drawer) const {
  // take account of differing font scale and main scale if we've hit
  // max or min font size.
  double f_scale = textDrawer_.fontScale();
  double spot_rad = 0.2 * drawOptions_.multipleBondOffset * f_scale;
  drawer.setColour(DrawColour(0.0, 0.0, 0.0));
  auto draw_spot = [&](const Point2D &cds) {
    bool ofp = drawer.fillPolys();
    drawer.setFillPolys(true);
    int olw = drawer.lineWidth();
    drawer.setLineWidth(0);
    drawer.drawArc(cds, spot_rad, 0, 360);
    drawer.setLineWidth(olw);
    drawer.setFillPolys(ofp);
  };
  // cds in draw coords

  auto draw_spots = [&](const Point2D &cds, int num_spots, double width,
                        int dir = 0) {
    Point2D ncds = cds;
    switch (num_spots) {
      case 3:
        draw_spot(ncds);
        if (dir) {
          ncds.y = cds.y - 0.5 * width + spot_rad;
        } else {
          ncds.x = cds.x - 0.5 * width + spot_rad;
        }
        draw_spot(ncds);
        if (dir) {
          ncds.y = cds.y + 0.5 * width - spot_rad;
        } else {
          ncds.x = cds.x + 0.5 * width - spot_rad;
        }
        draw_spot(ncds);
        /* fallthrough */
      case 1:
        draw_spot(cds);
        break;
      case 4:
        if (dir) {
          ncds.y = cds.y + 6.0 * spot_rad;
        } else {
          ncds.x = cds.x + 6.0 * spot_rad;
        }
        draw_spot(ncds);
        if (dir) {
          ncds.y = cds.y - 6.0 * spot_rad;
        } else {
          ncds.x = cds.x - 6.0 * spot_rad;
        }
        draw_spot(ncds);
        /* fallthrough */
      case 2:
        if (dir) {
          ncds.y = cds.y + 2.0 * spot_rad;
        } else {
          ncds.x = cds.x + 2.0 * spot_rad;
        }
        draw_spot(ncds);
        if (dir) {
          ncds.y = cds.y - 2.0 * spot_rad;
        } else {
          ncds.x = cds.x - 2.0 * spot_rad;
        }
        draw_spot(ncds);
        break;
    }
  };

  size_t rad_num = 0;
  for (auto atom : drawMol_->atoms()) {
    int num_rade = atom->getNumRadicalElectrons();
    if (!num_rade) {
      continue;
    }
    auto rad_rect = radicals_[rad_num].first;
    OrientType draw_or = radicals_[rad_num].second;
    if (draw_or == OrientType::N || draw_or == OrientType::S ||
        draw_or == OrientType::C) {
      draw_spots(rad_rect.trans_, num_rade, rad_rect.width_, 0);
    } else {
      draw_spots(rad_rect.trans_, num_rade, rad_rect.height_, 1);
    }
    ++rad_num;
  }
}

// ****************************************************************************
void DrawMol::resetEverything() {
  scale_ = 1.0;
  fontScale_ = 1.0;
  xMin_ = std::numeric_limits<double>::max() / 2.0;
  yMin_ = std::numeric_limits<double>::max() / 2.0;
  xMax_ = -std::numeric_limits<double>::max() / 2.0;
  yMax_ = -std::numeric_limits<double>::max() / 2.0;
  xRange_ = std::numeric_limits<double>::max();
  yRange_ = std::numeric_limits<double>::max();
  meanBondLengthSquare_ = 0.0;
  legendHeight_ = 1.0;
  atCds_.clear();
  bonds_.clear();
  atomicNums_.clear();
  atomSyms_.clear();
  atomLabels_.clear();
  highlights_.clear();
  annotations_.clear();
  legends_.clear();
  radicals_.clear();
}

// ****************************************************************************
std::pair<std::string, OrientType> DrawMol::getAtomSymbolAndOrientation(
    const Atom &atom) const {
  OrientType orient = getAtomOrientation(atom);
  std::string symbol = getAtomSymbol(atom, orient);

  return std::make_pair(symbol, orient);
}

// ****************************************************************************
std::string getAtomListText(const Atom &atom) {
  PRECONDITION(atom.hasQuery(), "no query");
  PRECONDITION(atom.getQuery()->getDescription() == "AtomOr", "bad query type");

  std::string res = "";
  if (atom.getQuery()->getNegation()) {
    res += "!";
  }
  res += "[";
  std::vector<int> vals;
  getAtomListQueryVals(atom.getQuery(), vals);
  for (unsigned int i = 0; i < vals.size(); ++i) {
    if (i != 0) {
      res += ",";
    }
    res += PeriodicTable::getTable()->getElementSymbol(vals[i]);
  }

  return res + "]";
}

// ****************************************************************************
std::string DrawMol::getAtomSymbol(const RDKit::Atom &atom,
                                OrientType orientation) const {
  if (drawOptions_.noAtomLabels) {
    return "";
  }
  // adds XML-like annotation for super- and sub-script, in the same manner
  // as MolDrawing.py. My first thought was for a LaTeX-like system,
  // obviously...
  std::string symbol;
  bool literal_symbol = true;
  unsigned int iso = atom.getIsotope();
  if (drawOptions_.atomLabels.find(atom.getIdx()) !=
      drawOptions_.atomLabels.end()) {
    // specified labels are trump: no matter what else happens we will show
    // them.
    symbol = drawOptions_.atomLabels.find(atom.getIdx())->second;
  } else if (atom.hasProp(common_properties::_displayLabel) ||
             atom.hasProp(common_properties::_displayLabelW)) {
    // logic here: if either _displayLabel or _displayLabelW is set, we will
    // definitely use one of those. if only one is set, we'll use that one if
    // both are set and the orientation is W then we'll use _displayLabelW,
    // otherwise _displayLabel

    std::string lbl;
    std::string lblw;
    atom.getPropIfPresent(common_properties::_displayLabel, lbl);
    atom.getPropIfPresent(common_properties::_displayLabelW, lblw);
    if (lbl.empty()) {
      lbl = lblw;
    }
    if (orientation == OrientType::W && !lblw.empty()) {
      symbol = lblw;
    } else {
      symbol = lbl;
    }
  } else if (atom.hasProp(common_properties::atomLabel)) {
    symbol = atom.getProp<std::string>(common_properties::atomLabel);
  } else if (drawOptions_.dummiesAreAttachments && atom.getAtomicNum() == 0 &&
             atom.getDegree() == 1) {
    symbol = "";
    literal_symbol = false;
  } else if (isAtomListQuery(&atom)) {
    symbol = getAtomListText(atom);
  } else if (isComplexQuery(&atom)) {
    symbol = "?";
  } else if (drawOptions_.atomLabelDeuteriumTritium &&
             atom.getAtomicNum() == 1 && (iso == 2 || iso == 3)) {
    symbol = ((iso == 2) ? "D" : "T");
    iso = 0;
  } else {
    literal_symbol = false;
    std::vector<std::string> preText, postText;

    // first thing after the symbol is the atom map
    if (atom.hasProp("molAtomMapNumber")) {
      std::string map_num = "";
      atom.getProp("molAtomMapNumber", map_num);
      postText.push_back(std::string(":") + map_num);
    }

    if (0 != atom.getFormalCharge()) {
      // charge always comes post the symbol
      int ichg = atom.getFormalCharge();
      std::string sgn = ichg > 0 ? std::string("+") : std::string("-");
      ichg = abs(ichg);
      if (ichg > 1) {
        sgn = std::to_string(ichg) + sgn;
      }
      // put the charge as a superscript
      postText.push_back(std::string("<sup>") + sgn + std::string("</sup>"));
    }

    int num_h = (atom.getAtomicNum() == 6 && atom.getDegree() > 0)
                    ? 0
                    : atom.getTotalNumHs();  // FIX: still not quite right

    if (drawOptions_.explicitMethyl && atom.getAtomicNum() == 6 &&
        atom.getDegree() == 1) {
      symbol += atom.getSymbol();
      num_h = atom.getTotalNumHs();
    }

    if (num_h > 0 && !atom.hasQuery()) {
      // the H text comes after the atomic symbol
      std::string h = "H";
      if (num_h > 1) {
        // put the number as a subscript
        h += std::string("<sub>") + std::to_string(num_h) + std::string("</sub>");
      }
      postText.push_back(h);
    }

    if (0 != iso &&
        ((drawOptions_.isotopeLabels && atom.getAtomicNum() != 0) ||
         (drawOptions_.dummyIsotopeLabels && atom.getAtomicNum() == 0))) {
      // isotope always comes before the symbol
      preText.push_back(std::string("<sup>") + std::to_string(iso) +
                        std::string("</sup>"));
    }

    symbol = "";
    for (const std::string &se : preText) {
      symbol += se;
    }

    // allenes need a C, but extend to any atom with degree 2 and both
    // bonds in a line.
    if (isLinearAtom(atom, atCds_) ||
        (atom.getAtomicNum() != 6 || atom.getDegree() == 0 || preText.size() ||
         postText.size())) {
      symbol += atom.getSymbol();
    }
    for (const std::string &se : postText) {
      symbol += se;
    }
  }

  if (literal_symbol && !symbol.empty()) {
    symbol = "<lit>" + symbol + "</lit>";
  }
  // cout << "Atom symbol " << atom.getIdx() << " : " << symbol << endl;
  return symbol;
}  // namespace RDKit

// ****************************************************************************
OrientType DrawMol::getAtomOrientation(const RDKit::Atom &atom) const {
  // cout << "Atomic " << atom.getAtomicNum() << " degree : "
  //      << atom.getDegree() << " : " << atom.getTotalNumHs() << endl;
  // anything with a slope of more than 70 degrees is vertical. This way,
  // the NH in an indole is vertical as RDKit lays it out normally (72ish
  // degrees) but the 2 amino groups of c1ccccc1C1CCC(N)(N)CC1 are E and W
  // when they are drawn at the bottom of the molecule.
  // NB - this assumes that the atom coords have already been inverted
  // in Y to put them in the draw frame where N is down and S is up.
  static const double VERT_SLOPE = tan(70.0 * M_PI / 180.0);

  auto &mol = atom.getOwningMol();
  const Point2D &at1_cds = atCds_[atom.getIdx()];
  Point2D nbr_sum(0.0, 0.0);
  // cout << "Nbours for atom : " << at1->getIdx() << endl;
  for (const auto &nbri : make_iterator_range(mol.getAtomBonds(&atom))) {
    const Bond *bond = mol[nbri];
    const Point2D &at2_cds =
        atCds_[bond->getOtherAtomIdx(atom.getIdx())];
    nbr_sum += at2_cds - at1_cds;
  }

  OrientType orient = OrientType::C;
  if (atom.getDegree()) {
    double islope = 1000.0;
    if (fabs(nbr_sum.x) > 1.0e-4) {
      islope = nbr_sum.y / nbr_sum.x;
    }
    if (fabs(islope) <= VERT_SLOPE) {
      if (nbr_sum.x > 0.0) {
        orient = OrientType::W;
      } else {
        orient = OrientType::E;
      }
    } else {
      if (nbr_sum.y > 0.0) {
        orient = OrientType::S;
      } else {
        orient = OrientType::N;
      }
    }
    // atoms of single degree should always be either W or E, never N or S.  If
    // either of the latter, make it E if the slope is close to vertical,
    // otherwise have it either as required.
    if (orient == OrientType::N || orient == OrientType::S) {
      if (atom.getDegree() == 1) {
        if (fabs(islope) > VERT_SLOPE) {
          orient = OrientType::E;
        } else {
          if (nbr_sum.x > 0.0) {
            orient = OrientType::W;
          } else {
            orient = OrientType::E;
          }
        }
      } else if (atom.getDegree() == 3) {
        // Atoms of degree 3 can sometimes have a bond pointing down with S
        // orientation or up with N orientation, which puts the H on the bond.
        auto &mol = atom.getOwningMol();
        const Point2D &at1_cds = atCds_[atom.getIdx()];
        for (const auto &nbri : make_iterator_range(mol.getAtomBonds(&atom))) {
          const Bond *bond = mol[nbri];
          const Point2D &at2_cds =
              atCds_[bond->getOtherAtomIdx(atom.getIdx())];
          Point2D bond_vec = at2_cds - at1_cds;
          double ang = atan(bond_vec.y / bond_vec.x) * 180.0 / M_PI;
          if (ang > 80.0 && ang < 100.0 && orient == OrientType::S) {
            orient = OrientType::S;
            break;
          } else if (ang < -80.0 && ang > -100.0 && orient == OrientType::N) {
            orient = OrientType::N;
            break;
          }
        }
      }
    }
  } else {
    // last check: degree zero atoms from the last three periods should have
    // the Hs first
    static int HsListedFirstSrc[] = {8, 9, 16, 17, 34, 35, 52, 53, 84, 85};
    std::vector<int> HsListedFirst(
        HsListedFirstSrc,
        HsListedFirstSrc + sizeof(HsListedFirstSrc) / sizeof(int));
    if (std::find(HsListedFirst.begin(), HsListedFirst.end(),
                  atom.getAtomicNum()) != HsListedFirst.end()) {
      orient = OrientType::W;
    } else {
      orient = OrientType::E;
    }
  }

  return orient;
}

// ****************************************************************************
void DrawMol::calcMeanBondLengthSquare() {
  // meanBondLengthSquare_ initialised to 0.0 in class declaration
  if (meanBondLengthSquare_ == 0.0) {
    unsigned int nBonds = 0;
    for (const auto &bond : drawMol_->bonds()) {
      meanBondLengthSquare_ +=
          (atCds_[bond->getBeginAtomIdx()] - atCds_[bond->getEndAtomIdx()])
              .lengthSq();
      ++nBonds;
    }
    meanBondLengthSquare_ /= nBonds;
  }
}

// ****************************************************************************
void DrawMol::extractLegend() {
  if (legend_.empty()) {
    return;
  }
  legendHeight_ = int(0.05 * double(height_));

  auto calc_legend_height = [&](const std::vector<std::string> &legend_bits,
                                double &total_width, double &total_height) {
    total_width = total_height = 0;
    for (auto bit : legend_bits) {
      double xMin, yMin, xMax, yMax;
      textDrawer_.getStringExtremes(bit, OrientType::N, xMin, yMin, xMax,
                                    yMax, true);
      total_height += yMax - yMin;
      total_width = std::max(total_width, xMax - xMin);
    }
  };

  std::vector<std::string> legend_bits;
  // split any strings on newlines
  std::string next_piece;
  for (auto c : legend_) {
    if (c == '\n') {
      if (!next_piece.empty()) {
        legend_bits.push_back(next_piece);
      }
      next_piece = "";
    } else {
      next_piece += c;
    }
  }
  if (!next_piece.empty()) {
    legend_bits.push_back(next_piece);
  }

  // work out a font scale that allows the pieces to fit
  double o_font_scale = textDrawer_.fontScale();
  double fsize = textDrawer_.fontSize();
  double new_font_scale = o_font_scale * drawOptions_.legendFontSize / fsize;
  textDrawer_.setFontScale(new_font_scale, true);
  double total_width, total_height;
  calc_legend_height(legend_bits, total_width, total_height);
  if (total_height > legendHeight_) {
    new_font_scale *= double(legendHeight_) / total_height;
    textDrawer_.setFontScale(new_font_scale, true);
    calc_legend_height(legend_bits, total_width, total_height);
  }
  if (total_width > width_) {
    new_font_scale *= double(width_) / total_width;
    textDrawer_.setFontScale(new_font_scale, true);
    calc_legend_height(legend_bits, total_width, total_height);
  }
  Point2D loc(width_ / 2, height_ - total_height);
  for (auto bit : legend_bits) {
    AnnotationType annot;
    annot.text_ = bit;
    annot.orient_ = OrientType::N;
    annot.col_ = drawOptions_.legendColour;
    annot.fontScale_ = new_font_scale;
    annot.scaleText_ = false;
    annot.rect_.trans_ = loc;
    double xMin, yMin, xMax, yMax;
    textDrawer_.getStringExtremes(bit, OrientType::N, xMin, yMin, xMax,
                                  yMax, true);
    annot.rect_.width_ = xMax - xMin;
    annot.rect_.height_ = yMax - yMin;
    loc.y += annot.rect_.height_;
    legends_.emplace_back(annot);
  }
}

// ****************************************************************************
void DrawMol::makeStandardBond(Bond *bond, double doubleBondOffset) {
  int begAt = bond->getBeginAtomIdx();
  int endAt = bond->getEndAtomIdx();
  std::pair<DrawColour, DrawColour> cols = getBondColours(bond);

  auto bt = bond->getBondType();
  if (bt == Bond::DOUBLE || bt == Bond::AROMATIC) {
    makeDoubleBondLines(bond, doubleBondOffset, cols);
  } else if (bt == Bond::SINGLE && (bond->getBondDir() == Bond::BEGINWEDGE ||
                                    bond->getBondDir() == Bond::BEGINDASH)) {
    makeWedgedBond(bond, cols);
  } else if (bt == Bond::SINGLE && bond->getBondDir() == Bond::UNKNOWN) {
    makeWavyBond(bond, cols);
  } else if (bt == Bond::DATIVE || bt == Bond::DATIVEL || bt == Bond::DATIVER) {
    makeDativeBond(bond, cols);
  } else if (bt == Bond::ZERO) {
    makeZeroBond(bond, cols, shortDashes);
  } else if (bt == Bond::HYDROGEN) {
    makeZeroBond(bond, cols, dots);
  } else {
    // in all other cases, we will definitely want to draw a line between
    // the two atoms
    Point2D end1, end2;
    adjustBondEndsForLabels(begAt, endAt, end1, end2);
    newBondLine(end1, end2, cols.first, cols.second, begAt, endAt,
                bond->getIdx(), noDash);
    if (Bond::TRIPLE == bt) {
      makeTripleBondLines(bond, doubleBondOffset, cols);
    }
  }
}

// ****************************************************************************
void DrawMol::makeQueryBond(Bond *bond, double doubleBondOffset) {
  PRECONDITION(bond->hasQuery(), "no query");
  const auto qry = bond->getQuery();

  auto begAt = bond->getBeginAtom();
  auto endAt = bond->getEndAtom();
  Point2D end1, end2;
  adjustBondEndsForLabels(begAt->getIdx(), endAt->getIdx(), end1, end2);
  Point2D sat1 = atCds_[begAt->getIdx()];
  Point2D sat2 = atCds_[endAt->getIdx()];
  atCds_[begAt->getIdx()] = end1;
  atCds_[endAt->getIdx()] = end2;

  const Point2D &at1_cds = atCds_[begAt->getIdx()];
  const Point2D &at2_cds = atCds_[endAt->getIdx()];
  auto midp = (at2_cds + at1_cds) / 2.;
  auto tdash = shortDashes;
  DrawColour queryColour{0.5, 0.5, 0.5};
  int at1Idx, at2Idx;

  bool drawGenericQuery = false;
  if (qry->getDescription() == "SingleOrDoubleBond") {
    at1Idx = begAt->getIdx();
    at2Idx = drawOptions_.splitBonds ? -1 : endAt->getIdx();
    newBondLine(at1_cds, midp, queryColour, queryColour, at1Idx,
                at2Idx, bond->getIdx(), noDash);
    Point2D l1s, l1f, l2s, l2f;
    calcDoubleBondLines(*drawMol_, doubleBondOffset, *bond, atCds_, l1s, l1f,
                        l2s, l2f);
    at1Idx = drawOptions_.splitBonds ? endAt->getIdx() : begAt->getIdx();
    at2Idx = drawOptions_.splitBonds ? -1 : endAt->getIdx();
    midp = (l1s + l1f) / 2.0;
    newBondLine(midp, l1f, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), noDash);
    midp = (l2s + l2f) / 2.0;
    newBondLine(midp, l2f, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), noDash);
  } else if (qry->getDescription() == "SingleOrAromaticBond") {
    at1Idx = begAt->getIdx();
    at2Idx = drawOptions_.splitBonds ? -1 : endAt->getIdx();
    newBondLine(at1_cds, midp, queryColour, queryColour, at1Idx,
                at2Idx, bond->getIdx(), noDash);
    Point2D l1s, l1f, l2s, l2f;
    calcDoubleBondLines(*drawMol_, doubleBondOffset, *bond, atCds_, l1s, l1f,
                        l2s, l2f);
    at1Idx = drawOptions_.splitBonds ? endAt->getIdx() : begAt->getIdx();
    at2Idx = drawOptions_.splitBonds ? -1 : endAt->getIdx();
    midp = (l1s + l1f) / 2.0;
    newBondLine(midp, l1f, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), noDash);
    midp = (l2s + l2f) / 2.0;
    newBondLine(midp, l2f, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), tdash);
  } else if (qry->getDescription() == "DoubleOrAromaticBond") {
    at1Idx = begAt->getIdx();
    at2Idx = drawOptions_.splitBonds ? -1 : endAt->getIdx();
    Point2D l1s, l1f, l2s, l2f;
    calcDoubleBondLines(*drawMol_, doubleBondOffset, *bond, atCds_, l1s, l1f,
                        l2s, l2f);
    midp = (l1s + l1f) / 2.0;
    newBondLine(l1s, midp, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), noDash);
    newBondLine(midp, l1f, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), noDash);
    at1Idx = drawOptions_.splitBonds ? endAt->getIdx() : begAt->getIdx();
    at2Idx = drawOptions_.splitBonds ? -1 : endAt->getIdx();
    midp = (l2s + l2f) / 2.0;
    newBondLine(l2s, midp, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), noDash);
    newBondLine(midp, l2f, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), tdash);
  } else if (qry->getDescription() == "BondNull") {
    at1Idx = begAt->getIdx();
    at2Idx = endAt->getIdx();
    newBondLine(at1_cds, at2_cds, queryColour, queryColour, at1Idx,
                at2Idx, bond->getIdx(), tdash);
  } else if (qry->getDescription() == "BondAnd" &&
             qry->endChildren() - qry->beginChildren() == 2) {
    auto q1 = *(qry->beginChildren());
    auto q2 = *(qry->beginChildren() + 1);

    if (q2->getDescription() == "BondOrder") {
      std::swap(q1, q2);
    }
    if (q1->getDescription() == "BondOrder" &&
        q2->getDescription() == "BondInRing") {
      size_t currNumBonds = bonds_.size();
      makeStandardBond(bond, doubleBondOffset);
      for (size_t i = currNumBonds; i < bonds_.size(); ++i) {
        bonds_[i]->lineColour_ = queryColour;
      }

      Point2D segment = at2_cds - at1_cds;
      if (!q2->getNegation()) {
        segment /= segment.length() * 6;
        Point2D r1 = Point2D(0.5 * segment.x - 0.866 * segment.y,
                             0.866 * segment.x + 0.5 * segment.y);
        Point2D r2 =
            Point2D(0.5 * r1.x - 0.866 * r1.y, 0.866 * r1.x + 0.5 * r1.y);
        std::vector<Point2D> pts{midp + segment, midp + r1, midp + r2,
                                 midp - segment, midp - r1, midp - r2,
                                 midp + segment};
        DrawShapePolyLine *pl = new DrawShapePolyLine(
            pts, 1, false, queryColour, false, begAt->getIdx(), endAt->getIdx(),
            bond->getIdx(), noDash);
        bonds_.emplace_back(std::unique_ptr<DrawShape>(pl));
      } else {
        segment /= segment.length() * 10;
        auto l = segment.length();
        Point2D p1 = midp + segment + Point2D(l, l);
        Point2D p2 = midp + segment - Point2D(l, l);
        std::vector<Point2D> pts{p1, p2};
        DrawShapeEllipse *ell =
            new DrawShapeEllipse(pts, 1, false, queryColour, false);
        bonds_.emplace_back(std::unique_ptr<DrawShape>(ell));
        p1 = midp - segment + Point2D(l, l);
        p2 = midp - segment - Point2D(l, l);
        pts = std::vector<Point2D>{p1, p2};
        ell = new DrawShapeEllipse(pts, 1, false, queryColour, false);
        bonds_.emplace_back(std::unique_ptr<DrawShape>(ell));
      }
    } else {
      drawGenericQuery = true;
    }
  } else {
    drawGenericQuery = true;
  }
  if (drawGenericQuery) {
    newBondLine(at1_cds, at2_cds, queryColour, queryColour, at1Idx,
                at2Idx, bond->getIdx(), dots);
    bonds_.back()->lineWidth_ = 1;
    bonds_.back()->scaleLineWidth_ = false;
  }
  atCds_[begAt->getIdx()] = sat1;
  atCds_[endAt->getIdx()] = sat2;
}

// ****************************************************************************
void DrawMol::makeDoubleBondLines(Bond *bond, double doubleBondOffset,
                                  const std::pair<DrawColour, DrawColour> &cols) {
  Point2D end1, end2;
  int at1Idx = bond->getBeginAtomIdx();
  int at2Idx = bond->getEndAtomIdx();
  adjustBondEndsForLabels(at1Idx, at2Idx, end1, end2);

  Point2D l1s, l1f, l2s, l2f, sat1, sat2;
  sat1 = atCds_[at1Idx];
  atCds_[at1Idx] = end1;
  sat2 = atCds_[at2Idx];
  atCds_[at2Idx] = end2;
  calcDoubleBondLines(*drawMol_, doubleBondOffset, *bond, atCds_, l1s, l1f, l2s,
                      l2f);
  bool orig_slw = drawOptions_.scaleBondWidth;
  //  if (highlight_bond) {
  //    d2d.drawOptions().scaleBondWidth =
  //        d2d.drawOptions().scaleHighlightBondWidth;
  //  }
  int bondIdx = bond->getIdx();
  newBondLine(l1s, l1f, cols.first, cols.second, at1Idx, at2Idx,
              bondIdx, noDash);
  if (bond->getBondType() == Bond::AROMATIC) {
    newBondLine(l2s, l2f, cols.first, cols.second, at1Idx, at2Idx,
                bondIdx, dashes);
  } else {
    newBondLine(l2s, l2f, cols.first, cols.second, at1Idx, at2Idx,
                bondIdx, noDash);
  }
  drawOptions_.scaleBondWidth = orig_slw;
  atCds_[at1Idx] = sat1;
  atCds_[at2Idx] = sat2;
}

// ****************************************************************************
void DrawMol::makeTripleBondLines(
    Bond *bond, double doubleBondOffset,
    const std::pair<DrawColour, DrawColour> &cols) {
  Point2D end1, end2;
  int at1Idx = bond->getBeginAtomIdx();
  int at2Idx = bond->getEndAtomIdx();
  adjustBondEndsForLabels(at1Idx, at2Idx, end1, end2);

  Point2D l1s, l1f, l2s, l2f, sat1, sat2;
  sat1 = atCds_[at1Idx];
  atCds_[at1Idx] = end1;
  sat2 = atCds_[at2Idx];
  atCds_[at2Idx] = end2;
  int bondIdx = bond->getIdx();
  calcTripleBondLines(doubleBondOffset, *bond, atCds_, l1s, l1f, l2s, l2f);
  bool orig_slw = drawOptions_.scaleBondWidth;
  newBondLine(l1s, l1f, cols.first, cols.second, at1Idx, at2Idx,
              bondIdx, noDash);
  newBondLine(l2s, l2f, cols.first, cols.second, at1Idx, at2Idx,
              bondIdx, noDash);
  drawOptions_.scaleBondWidth = orig_slw;
  atCds_[at1Idx] = sat1;
  atCds_[at2Idx] = sat2;
}

// ****************************************************************************
void DrawMol::makeWedgedBond(Bond *bond,
                             const std::pair<DrawColour, DrawColour> &cols) {
  // swap the direction if at1 has does not have stereochem set
  // or if at2 does have stereochem set and the bond starts there
  auto at1 = bond->getBeginAtom();
  auto at2 = bond->getEndAtom();
  auto col1 = cols.first;
  auto col2 = cols.second;
  auto inverted = false;
  if ((at1->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
       at1->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW) ||
      (at1->getIdx() != bond->getBeginAtomIdx() &&
       (at2->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
        at2->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW))) {
    // std::cerr << "  swap" << std::endl;
    std::swap(at1, at2);
    std::swap(col1, col2);
    inverted = true;
  }
  if (drawOptions_.singleColourWedgeBonds) {
    col1 = drawOptions_.symbolColour;
    col2 = drawOptions_.symbolColour;
  }
  Point2D end1, end2;
  // If either of the atoms has a label, make the padding a bit bigger
  // so the end of the wedge doesn't run up to the atom symbol.
  // Obviously, if we ever change how the padding round the label is
  // calculated, this won't work.
  if (atomLabels_[at1->getIdx()] || atomLabels_[at2->getIdx()]) {
    meanBondLengthSquare_ *= 2.0;
  }
  adjustBondEndsForLabels(at1->getIdx(), at2->getIdx(), end1, end2);
  if (atomLabels_[at1->getIdx()] || atomLabels_[at2->getIdx()]) {
    meanBondLengthSquare_ /= 2.0;
  }
  const Point2D &at1_cds = atCds_[at1->getIdx()];
  const Point2D &at2_cds = atCds_[at2->getIdx()];

  Point2D perp = calcPerpendicular(at1_cds, at2_cds);
  Point2D disp = perp * 0.15;
  Point2D t1 = end2 + disp;
  Point2D t2 = end2 - disp;
  std::vector<Point2D> pts{end1, t1, t2};

  // deliberately not scaling highlighted bond width
  DrawShape *s;
  if (Bond::BEGINWEDGE == bond->getBondDir()) {
    s = new DrawShapeSolidWedge(pts, col1, col2, inverted,
                                drawOptions_.splitBonds, at1->getIdx(),
                                at2->getIdx(), bond->getIdx());
  } else {
    s = new DrawShapeDashedWedge(pts, col1, col2, inverted, at1->getIdx(),
                                 at2->getIdx(), bond->getIdx());
  }
  bonds_.emplace_back(std::unique_ptr<DrawShape>(s));
}

// ****************************************************************************
void DrawMol::makeWavyBond(Bond *bond,
                           const std::pair<DrawColour, DrawColour> &cols) {
  auto at1 = bond->getBeginAtom();
  auto at2 = bond->getEndAtom();
  Point2D end1, end2;
  adjustBondEndsForLabels(at1->getIdx(), at2->getIdx(), end1, end2);
  std::vector<Point2D> pts{end1, end2};
  DrawShapeWavyLine *s = new DrawShapeWavyLine(pts, drawOptions_.bondLineWidth,
                                               false, cols.first, cols.second);
  bonds_.emplace_back(std::unique_ptr<DrawShape>(s));
}

// ****************************************************************************
void DrawMol::makeDativeBond(Bond *bond,
                             const std::pair<DrawColour, DrawColour> &cols) {
  auto at1 = bond->getBeginAtom();
  auto at2 = bond->getEndAtom();
  Point2D end1, end2;
  adjustBondEndsForLabels(at1->getIdx(), at2->getIdx(), end1, end2);

  Point2D mid = (end1 + end2) * 0.5;
  newBondLine(end1, mid, cols.first, cols.first, at1->getIdx(),
              at2->getIdx(), bond->getIdx(), noDash);
  std::vector<Point2D> pts{mid, end2};
  DrawShapeArrow *a = new DrawShapeArrow(pts, drawOptions_.bondLineWidth, false,
                                         cols.second, true, 0.2, M_PI / 6);
  bonds_.emplace_back(std::unique_ptr<DrawShape>(a));
}

// ****************************************************************************
void DrawMol::makeZeroBond(Bond *bond,
                           const std::pair<DrawColour, DrawColour> &cols,
                           const DashPattern &dashPattern) {
  auto at1 = bond->getBeginAtom();
  auto at2 = bond->getEndAtom();
  Point2D end1, end2;
  adjustBondEndsForLabels(at1->getIdx(), at2->getIdx(), end1, end2);
  newBondLine(end1, end2, cols.first, cols.second, at1->getIdx(),
              at2->getIdx(), bond->getIdx(), dashPattern);
}

// ****************************************************************************
void DrawMol::adjustBondEndsForLabels(int begAtIdx, int endAtIdx,
                                      Point2D &begCds, Point2D &endCds) {
  double padding = 0.025 * meanBondLengthSquare_;
  begCds = atCds_[begAtIdx];
  endCds = atCds_[endAtIdx];
  if (atomLabels_[begAtIdx]) {
    adjustBondEndForString(atCds_[begAtIdx], atCds_[endAtIdx], padding,
                           atomLabels_[begAtIdx]->rects_, begCds);
  }
  if (atomLabels_[endAtIdx]) {
    adjustBondEndForString(atCds_[endAtIdx], atCds_[begAtIdx], padding,
                           atomLabels_[endAtIdx]->rects_, endCds);
  }
}

// ****************************************************************************
void DrawMol::newBondLine(const Point2D &pt1, const Point2D &pt2,
                          const DrawColour &col1, const DrawColour &col2,
                          int atom1Idx, int atom2Idx,
                          int bondIdx, const DashPattern &dashPattern) {
  if (col1 == col2 && !drawOptions_.splitBonds) {
    std::vector<Point2D> pts{pt1, pt2};
    DrawShape *b =
        new DrawShapeSimpleLine(pts, drawOptions_.bondLineWidth, false, col1,
                                atom1Idx, atom2Idx, bondIdx, dashPattern);
    bonds_.emplace_back(std::unique_ptr<DrawShape>(b));
  } else {
    Point2D mid = (pt1 + pt2) / 2.0;
    std::vector<Point2D> pts1{pt1, mid};
    int at1Idx = atom1Idx;
    int at2Idx = drawOptions_.splitBonds ? -1 : atom2Idx;
    DrawShape *b1 =
        new DrawShapeSimpleLine(pts1, drawOptions_.bondLineWidth, false, col1,
                                at1Idx, at2Idx, bondIdx, dashPattern);
    bonds_.emplace_back(std::unique_ptr<DrawShape>(b1));
    at1Idx = drawOptions_.splitBonds ? atom2Idx : atom1Idx;
    std::vector<Point2D> pts2{mid, pt2};
    DrawShape *b2 =
        new DrawShapeSimpleLine(pts2, drawOptions_.bondLineWidth, false, col2,
                                at1Idx, at2Idx, bondIdx, dashPattern);
    bonds_.emplace_back(std::unique_ptr<DrawShape>(b2));
  }
}

// ****************************************************************************
std::pair<DrawColour, DrawColour> DrawMol::getBondColours(Bond *bond) {
  DrawColour col1, col2;

  bool highlight_bond = false;
  if (highlightBonds_ &&
      std::find(highlightBonds_->begin(), highlightBonds_->end(),
                bond->getIdx()) != highlightBonds_->end()) {
    highlight_bond = true;
  }

  if (bondColours_) {
    col1 = (*bondColours_)[bond->getIdx()].first;
    col2 = (*bondColours_)[bond->getIdx()].second;
  } else {
    if (!highlight_bond || drawOptions_.continuousHighlight) {
      int at1_idx = bond->getBeginAtomIdx();
      col1 = getColour(at1_idx, drawOptions_, atomicNums_, highlightAtoms_,
                       highlightAtomMap_);
      int at2_idx = bond->getEndAtomIdx();
      col2 = getColour(at2_idx, drawOptions_, atomicNums_, highlightAtoms_,
                       highlightAtomMap_);
    } else {
      if (highlightBondMap_ && highlightBondMap_->find(bond->getIdx()) !=
                                    highlightBondMap_->end()) {
        col1 = col2 = highlightBondMap_->find(bond->getIdx())->second;
      } else {
        col1 = col2 = drawOptions_.highlightColour;
      }
    }
  }

  return std::make_pair(col1, col2);
}

// ****************************************************************************
void DrawMol::makeContinuousHighlights() {
  int tgt_lw = getHighlightBondWidth(drawOptions_, -1, nullptr);
  if (tgt_lw < 2) {
    tgt_lw = 2;
  }
  if (!drawOptions_.continuousHighlight) {
    tgt_lw /= 4;
  }

  if (highlightBonds_) {
    makeBondHighlightLines(tgt_lw);
  }
  std::cout << "number of highlights after bonds: " << highlights_.size() << std::endl;
  if (highlightAtoms_) {
    makeAtomEllipseHighlights(tgt_lw);
  }
  std::cout << "number of highlights after bonds + atoms : " << highlights_.size() << std::endl;
}

// ****************************************************************************
void DrawMol::makeAtomCircleHighlights() {
  PRECONDITION(highlightAtoms_ != nullptr, "no highlight atoms");
  DrawColour col;
  for (auto at : drawMol_->atoms()) {
    unsigned int thisIdx = at->getIdx();
    if (std::find(highlightAtoms_->begin(), highlightAtoms_->end(), thisIdx) !=
        highlightAtoms_->end()) {
      std::cout << "highlights for " << thisIdx << std::endl;
      if (highlightAtomMap_ &&
          highlightAtomMap_->find(thisIdx) != highlightAtomMap_->end()) {
        col = highlightAtomMap_->find(thisIdx)->second;
      } else {
        col = drawOptions_.highlightColour;
      }
      double radius = drawOptions_.highlightRadius;
      if (highlightRadii_ &&
          highlightRadii_->find(thisIdx) != highlightRadii_->end()) {
        radius = highlightRadii_->find(thisIdx)->second;
      }
      Point2D offset(radius, radius);
      Point2D p1 = atCds_[thisIdx] - offset;
      Point2D p2 = atCds_[thisIdx] + offset;
      std::vector<Point2D> pts{p1, p2};
      DrawShape *ell = new DrawShapeEllipse(pts, 2, false, col, true, thisIdx);
      highlights_.emplace_back(std::unique_ptr<DrawShape>(ell));
    }
  }
}

// ****************************************************************************
void DrawMol::makeAtomEllipseHighlights(int lineWidth) {
  PRECONDITION(highlightAtoms_ != nullptr, "no highlight atoms");
  if (!drawOptions_.fillHighlights) {
    // we need a narrower circle
    lineWidth /= 2;
  }
  for (auto atom : drawMol_->atoms()) {
    unsigned int thisIdx = atom->getIdx();
    if (std::find(highlightAtoms_->begin(), highlightAtoms_->end(),
                  thisIdx) != highlightAtoms_->end()) {
      DrawColour col = drawOptions_.highlightColour;
      if (highlightAtomMap_ &&
          highlightAtomMap_->find(thisIdx) != highlightAtomMap_->end()) {
        col = highlightAtomMap_->find(thisIdx)->second;
      }
      Point2D centre = atCds_[thisIdx];
      double xradius, yradius;
      if (highlightRadii_ && highlightRadii_->find(thisIdx) != highlightRadii_->end()) {
        xradius = highlightRadii_->find(thisIdx)->second;
      } else {
        xradius = drawOptions_.highlightRadius;
      }
      yradius = xradius;
      if (!drawOptions_.atomHighlightsAreCircles && atomLabels_[thisIdx]) {
        double xMin, yMin, xMax, yMax;
        xMin = yMin = std::numeric_limits<double>::max();
        xMax = yMax = -std::numeric_limits<double>::max();
        atomLabels_[thisIdx]->findExtremes(xMin, xMax, yMin, yMax);
        static const double root_2 = sqrt(2.0);
        xradius = std::max(xradius, root_2 * 0.5 * (xMax - xMin));
        yradius = std::max(yradius, root_2 * 0.5 * (yMax - yMin));
        centre.x = 0.5 * (xMax + xMin);
        centre.y = 0.5 * (yMax + yMin);
      }
      Point2D offset(xradius, yradius);
      Point2D p1 = centre - offset;
      Point2D p2 = centre + offset;
      std::vector<Point2D> pts{p1, p2};
      DrawShape *ell = new DrawShapeEllipse(pts, lineWidth, true, col, true,
                                            thisIdx);
      highlights_.emplace_back(std::unique_ptr<DrawShape>(ell));
    }
  }
  
}

// ****************************************************************************
void DrawMol::makeBondHighlightLines(int lineWidth) {
  PRECONDITION(highlightBonds_ != nullptr, "no highlight bonds");

  for (auto atom : drawMol_->atoms()) {
    unsigned int thisIdx = atom->getIdx();
    for (const auto &nbri : make_iterator_range(drawMol_->getAtomBonds(atom))) {
      const Bond *bond = (*drawMol_)[nbri];
      int nbrIdx = bond->getOtherAtomIdx(thisIdx);
      if (nbrIdx < static_cast<unsigned int>(atCds_.size()) &&
          nbrIdx > thisIdx) {
        if (std::find(highlightBonds_->begin(), highlightBonds_->end(),
                      bond->getIdx()) != highlightBonds_->end()) {
          DrawColour col = drawOptions_.highlightColour;
          if (highlightBondMap_ && highlightBondMap_->find(bond->getIdx()) !=
                                       highlightBondMap_->end()) {
            col = highlightBondMap_->find(bond->getIdx())->second;
          }
          std::vector<Point2D> pts{atCds_[thisIdx], atCds_[nbrIdx]};
          std::cout << "bond highlight for " << thisIdx << " to " << nbrIdx
                    << " width = " << lineWidth << std::endl;
          DrawShape *hb = new DrawShapeSimpleLine(
              pts, lineWidth, drawOptions_.scaleHighlightBondWidth, col,
              thisIdx, nbrIdx, bond->getIdx(), noDash);
          highlights_.emplace_back(std::unique_ptr<DrawShape>(hb));
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMol::calcAnnotationPosition(const Atom *atom,
                                     AnnotationType &annot) const {
  PRECONDITION(atom, "no atom");
  if (annot.text_.empty()) {
    annot.rect_.width_ = -1.0;  // so we know it's not valid.
    return;
  }

  std::cout << "calc annot pos for " << annot.text_ << std::endl;
  Point2D const &at_cds = atCds_[atom->getIdx()];
  annot.rect_.trans_.x = at_cds.x;
  annot.rect_.trans_.y = at_cds.y;
  calcAnnotationDims(annot);
  double start_ang = getNoteStartAngle(atom);
  calcAtomAnnotationPosition(atom, start_ang, annot);
}

// ****************************************************************************
void DrawMol::calcAnnotationPosition(const Bond *bond,
                                     AnnotationType &annot) const {
  PRECONDITION(bond, "no bond");
  if (annot.text_.empty()) {
    annot.rect_.width_ = -1.0;  // so we know it's not valid.
  }
  Point2D const &at1_cds = atCds_[bond->getBeginAtomIdx()];
  Point2D const &at2_cds = atCds_[bond->getEndAtomIdx()];
  Point2D perp = calcPerpendicular(at1_cds, at2_cds);
  Point2D bond_vec = at1_cds.directionVector(at2_cds);
  double bond_len = (at1_cds - at2_cds).length();
  std::vector<double> mid_offsets{0.5, 0.33, 0.66, 0.25, 0.75};
  double offset_step = drawOptions_.multipleBondOffset;
  StringRect least_worst_rect = StringRect();
  least_worst_rect.clash_score_ = 100;
  for (auto mo : mid_offsets) {
    Point2D mid = at1_cds + bond_vec * bond_len * mo;
    for (int j = 1; j < 6; ++j) {
      if (j == 1 && bond->getBondType() > 1) {
        continue;  // multiple bonds will need a bigger offset.
      }
      double offset = j * offset_step;
      annot.rect_.trans_ = mid + perp * offset;
      calcAnnotationDims(annot);
      int clash_score = doesNoteClash(annot);
      if (!clash_score) {
        return;
      }
      if (clash_score < least_worst_rect.clash_score_) {
        least_worst_rect = annot.rect_;
      }
      annot.rect_.trans_ = mid - perp * offset;
      clash_score = doesNoteClash(annot);
      if (!clash_score) {
        return;
      }
      if (clash_score < least_worst_rect.clash_score_) {
        least_worst_rect = annot.rect_;
      }
    }
  }
}

// ****************************************************************************
double DrawMol::getNoteStartAngle(const Atom *atom) const {
  if (atom->getDegree() == 0) {
    return M_PI / 2.0;
  }
  Point2D at_cds = atCds_[atom->getIdx()];
  std::vector<Point2D> bond_vecs;
  for (const auto &nbr : make_iterator_range(drawMol_->getAtomNeighbors(atom))) {
    Point2D bond_vec = at_cds.directionVector(atCds_[nbr]);
    bond_vec.normalize();
    bond_vecs.emplace_back(bond_vec);
  }

  Point2D ret_vec;
  if (bond_vecs.size() == 1) {
    if (!atomLabels_[atom->getIdx()]) {
      // go with perpendicular to bond.  This is mostly to avoid getting
      // a zero at the end of a bond to carbon, which looks like a black
      // oxygen atom in the default font in SVG and PNG.
      ret_vec.x = bond_vecs[0].y;
      ret_vec.y = -bond_vecs[0].x;
    } else {
      // go opposite end
      ret_vec = -bond_vecs[0];
    }
  } else if (bond_vecs.size() == 2) {
    ret_vec = bond_vecs[0] + bond_vecs[1];
    if (ret_vec.lengthSq() > 1.0e-6) {
      if (!atom->getNumImplicitHs() || atom->getAtomicNum() == 6) {
        // prefer outside the angle, unless there are Hs that will be in
        // the way, probably.
        ret_vec *= -1.0;
      }
    } else {
      // it must be a -# or == or some such.  Take perpendicular to
      // one of them
      ret_vec.x = -bond_vecs.front().y;
      ret_vec.y = bond_vecs.front().x;
      ret_vec.normalize();
    }
  } else {
    // just take 2 that are probably adjacent
    double discrim = 4.0 * M_PI / bond_vecs.size();
    for (size_t i = 0; i < bond_vecs.size() - 1; ++i) {
      for (size_t j = i + 1; j < bond_vecs.size(); ++j) {
        double ang = acos(bond_vecs[i].dotProduct(bond_vecs[j]));
        if (ang < discrim) {
          ret_vec = bond_vecs[i] + bond_vecs[j];
          ret_vec.normalize();
          discrim = -1.0;
          break;
        }
      }
    }
    if (discrim > 0.0) {
      ret_vec = bond_vecs[0] + bond_vecs[1];
      ret_vec *= -1.0;
    }
  }

  // start angle is the angle between ret_vec and the x axis
  return atan2(ret_vec.y, ret_vec.x);
}

// ****************************************************************************
void DrawMol::calcAtomAnnotationPosition(const Atom *atom,
                                         double start_ang,
                                         AnnotationType &annot) const {
  Point2D const &at_cds = atCds_[atom->getIdx()];

  double rad_step = 0.25;
  StringRect least_worst_rect = StringRect();
  least_worst_rect.clash_score_ = 100;
  for (int j = 1; j < 4; ++j) {
    double note_rad = j * rad_step;
    // experience suggests if there's an atom symbol, the close in
    // radius won't work.
    if (j == 1 && atomLabels_[atom->getIdx()]) {
      continue;
    }
    // scan at 30 degree intervals around the atom looking for somewhere
    // clear for the annotation.
    for (int i = 0; i < 12; ++i) {
      double ang = start_ang + i * 30.0 * M_PI / 180.0;
      annot.rect_.trans_.x = at_cds.x + cos(ang) * note_rad;
      annot.rect_.trans_.y = at_cds.y + sin(ang) * note_rad;
      int clash_score = doesNoteClash(annot);
      std::cout << "annot " << annot.text_ << " j = " << j << " i = " << i
                << " clash score = " << clash_score << std::endl;
      if (!clash_score) {
        return;
      } else {
        if (clash_score < least_worst_rect.clash_score_) {
          least_worst_rect = annot.rect_;
        }
      }
    }
  }
  std::cout << "XXXXZ " << annot.text_ << " gone with least worst rect"
            << std::endl << std::endl;
  annot.rect_ = least_worst_rect;
}

// ****************************************************************************
void DrawMol::calcAnnotationDims(AnnotationType &annot) const {
  std::vector<std::shared_ptr<StringRect>> rects;
  std::vector<TextDrawType> draw_modes;
  std::vector<char> draw_chars;

  double full_font_scale = textDrawer_.fontScale();
  textDrawer_.setFontScale(drawOptions_.annotationFontScale, true);
  textDrawer_.getStringRects(annot.text_, OrientType::C, rects, draw_modes,
                             draw_chars, false, annot.align_);
  textDrawer_.setFontScale(full_font_scale, true);
  annot.rect_.width_ = 0;
  annot.rect_.height_ = 0;
  for (auto &rect : rects) {
    annot.rect_.width_ += rect->width_;
    annot.rect_.height_ = std::max(rect->height_, annot.rect_.height_);
  }
}

// ****************************************************************************
int DrawMol::doesNoteClash(const AnnotationType &annot) const {
  double padding = scale_ * 0.02;
  return doesRectClash(annot.rect_, padding);
}

// ****************************************************************************
int DrawMol::doesRectClash(const StringRect &rect, double padding) const {
  for (auto &bond : bonds_) {
    if (bond->doesRectClash(rect, padding)) {
      return 1;
    }
  }
  for (auto &al : atomLabels_) {
    if (al && al->doesRectClash(rect, padding)) {
      return 2;
    }
  }
  for (auto &a : annotations_) {
    if (a.rect_.doesItIntersect(rect, padding)) {
      return 3;
    }
  }
  for (auto &hl : highlights_) {
    if (hl->doesRectClash(rect, padding)) {
      return 4;
    }
  }
  return 0;
}

// ****************************************************************************
OrientType DrawMol::calcRadicalRect(const Atom *atom, StringRect &rad_rect) const {
  int num_rade = atom->getNumRadicalElectrons();
  double spot_rad = 0.2 * drawOptions_.multipleBondOffset;
  Point2D const &at_cds = atCds_[atom->getIdx()];
  OrientType orient = atomSyms_[atom->getIdx()].second;
  double rad_size = (4 * num_rade - 2) * spot_rad;
  double x_min, y_min, x_max, y_max;
  if (atomLabels_[atom->getIdx()]) {
    x_min = y_min = std::numeric_limits<double>::max();
    x_max = y_max = -x_min;
    atomLabels_[atom->getIdx()]->findExtremes(x_min, x_max, y_min, y_max);
  } else {
    x_min = at_cds.x - 3 * spot_rad * textDrawer_.fontScale();
    x_max = at_cds.x + 3 * spot_rad * textDrawer_.fontScale();
    y_min = at_cds.y - 3 * spot_rad * textDrawer_.fontScale();
    y_max = at_cds.y + 3 * spot_rad * textDrawer_.fontScale();
  }

  auto try_all = [&](OrientType ornt) -> bool {
    if (!doesRectClash(rad_rect, 0.0)) {
      return true;
    } else {
      return false;
    }
  };

  auto try_north = [&]() -> bool {
    rad_rect.width_ = rad_size * textDrawer_.fontScale();
    rad_rect.height_ = spot_rad * 3.0 * textDrawer_.fontScale();
    rad_rect.trans_.x = at_cds.x;
    rad_rect.trans_.y = y_max + 0.5 * rad_rect.height_;
    return try_all(OrientType::N);
  };
  auto try_south = [&]() -> bool {
    rad_rect.width_ = rad_size * textDrawer_.fontScale();
    rad_rect.height_ = spot_rad * 3.0 * textDrawer_.fontScale();
    rad_rect.trans_.x = at_cds.x;
    rad_rect.trans_.y = y_min - 0.5 * rad_rect.height_;
    return try_all(OrientType::S);
  };
  auto try_east = [&]() -> bool {
    rad_rect.trans_.x = x_max + 3.0 * spot_rad * textDrawer_.fontScale();
    rad_rect.trans_.y = at_cds.y;
    rad_rect.width_ = spot_rad * 1.5 * textDrawer_.fontScale();
    rad_rect.height_ = rad_size * textDrawer_.fontScale();
    return try_all(OrientType::E);
  };
  auto try_west = [&]() -> bool {
    rad_rect.trans_.x = x_min - 3.0 * spot_rad * textDrawer_.fontScale();
    rad_rect.trans_.y = at_cds.y;
    rad_rect.width_ = spot_rad * 1.5 * textDrawer_.fontScale();
    rad_rect.height_ = rad_size * textDrawer_.fontScale();
    return try_all(OrientType::W);
  };

  auto try_rads = [&](OrientType ornt) -> bool {
    switch (ornt) {
      case OrientType::N:
      case OrientType::C:
        return try_north();
      case OrientType::E:
        return try_east();
      case OrientType::S:
        return try_south();
      case OrientType::W:
        return try_west();
    }
    return false;
  };
  if (try_rads(orient)) {
    return orient;
  }
  OrientType all_ors[4] = {OrientType::N, OrientType::E, OrientType::S,
                           OrientType::W};
  for (int io = 0; io < 4; ++io) {
    if (orient != all_ors[io]) {
      if (try_rads(all_ors[io])) {
        return all_ors[io];
      }
    }
  }
  // stick them N irrespective of a clash whilst muttering "sod it"
  // under our breath.
  try_north();
  return OrientType::N;
}

// ****************************************************************************
void centerMolForDrawing(RWMol &mol, int confId) {
  auto &conf = mol.getConformer(confId);
  RDGeom::Transform3D tf;
  auto centroid = MolTransforms::computeCentroid(conf);
  centroid *= -1;
  tf.SetTranslation(centroid);
  MolTransforms::transformConformer(conf, tf);
  MolTransforms::transformMolSubstanceGroups(mol, tf);
}

// ****************************************************************************
void prepareStereoGroups(RWMol &mol) {
  auto sgs = mol.getStereoGroups();
  if (sgs.size() == 1) {
    boost::dynamic_bitset<> chiralAts(mol.getNumAtoms());
    for (const auto atom : mol.atoms()) {
      if (atom->getChiralTag() > Atom::ChiralType::CHI_UNSPECIFIED &&
          atom->getChiralTag() < Atom::ChiralType::CHI_OTHER) {
        chiralAts.set(atom->getIdx(), 1);
      }
    }
    for (const auto atm : sgs[0].getAtoms()) {
      chiralAts.set(atm->getIdx(), 0);
    }
    if (chiralAts.none()) {
      // all specified chiral centers are accounted for by this StereoGroup.
      if (sgs[0].getGroupType() == StereoGroupType::STEREO_OR ||
          sgs[0].getGroupType() == StereoGroupType::STEREO_AND) {
        std::vector<StereoGroup> empty;
        mol.setStereoGroups(std::move(empty));
        std::string label = sgs[0].getGroupType() == StereoGroupType::STEREO_OR
                                ? "OR enantiomer"
                                : "AND enantiomer";
        mol.setProp(common_properties::molNote, label);
      }
      // clear the chiral codes on the atoms so that we don't
      // inadvertently draw them later
      for (const auto atm : sgs[0].getAtoms()) {
        mol.getAtomWithIdx(atm->getIdx())->clearProp(common_properties::_CIPCode);
      }
    }
  }
}

// ****************************************************************************
bool isLinearAtom(const Atom &atom, const std::vector<Point2D> &atCds) {
  if (atom.getDegree() == 2) {
    Point2D bond_vecs[2];
    Bond::BondType bts[2];
    Point2D const &at1_cds = atCds[atom.getIdx()];
    ROMol const &mol = atom.getOwningMol();
    int i = 0;
    for (const auto &nbr : make_iterator_range(mol.getAtomNeighbors(&atom))) {
      Point2D bond_vec = at1_cds.directionVector(atCds[nbr]);
      bond_vec.normalize();
      bond_vecs[i] = bond_vec;
      bts[i] = mol.getBondBetweenAtoms(atom.getIdx(), nbr)->getBondType();
      ++i;
    }
    return (bts[0] == bts[1] && bond_vecs[0].dotProduct(bond_vecs[1]) < -0.95);
  }
  return false;
}

// ****************************************************************************
DrawColour getColour(int atom_idx, const MolDrawOptions &drawOptions,
                     const std::vector<int> &atomicNums,
                     const std::vector<int> *highlightAtoms,
                     const std::map<int, DrawColour> *highlightMap) {
  PRECONDITION(atom_idx >= 0, "bad atom_idx");
  PRECONDITION(rdcast<int>(atomicNums.size()) > atom_idx, "bad atom_idx");

  DrawColour retval = getColourByAtomicNum(atomicNums[atom_idx], drawOptions);
  // set contents of highlight_atoms to red
  if (!drawOptions.circleAtoms && !drawOptions.continuousHighlight) {
    if (highlightAtoms &&
        highlightAtoms->end() !=
            find(highlightAtoms->begin(), highlightAtoms->end(), atom_idx)) {
      retval = drawOptions.highlightColour;
    }
    // over-ride with explicit colour from highlightMap if there is one
    if (highlightMap) {
      auto p = highlightMap->find(atom_idx);
      if (p != highlightMap->end()) {
        retval = p->second;
      }
    }
  }
  return retval;
}

// ****************************************************************************
DrawColour getColourByAtomicNum(int atomic_num,
                                const MolDrawOptions &drawOptions) {
  DrawColour res;
  if (drawOptions.atomColourPalette.find(atomic_num) !=
      drawOptions.atomColourPalette.end()) {
    res = drawOptions.atomColourPalette.find(atomic_num)->second;
  } else if (atomic_num != -1 &&
             drawOptions.atomColourPalette.find(-1) !=
                 drawOptions.atomColourPalette.end()) {
    // if -1 is in the palette, we use that for undefined colors
    res = drawOptions.atomColourPalette.find(-1)->second;
  } else {
    // if all else fails, default to black:
    res = DrawColour(0, 0, 0);
  }
  return res;
}

// ****************************************************************************
int getHighlightBondWidth(
    MolDrawOptions &drawOptions, int bond_idx,
    const std::map<int, int> *highlight_linewidth_multipliers) {
  int bwm = drawOptions.highlightBondWidthMultiplier;
  // if we're not doing filled highlights, the lines need to be narrower
  if (!drawOptions.fillHighlights) {
    bwm /= 2;
    if (bwm < 1) {
      bwm = 1;
    }
  }

  if (highlight_linewidth_multipliers &&
      !highlight_linewidth_multipliers->empty()) {
    auto it = highlight_linewidth_multipliers->find(bond_idx);
    if (it != highlight_linewidth_multipliers->end()) {
      bwm = it->second;
    }
  }
  int tgt_lw = drawOptions.bondLineWidth * bwm;

  return tgt_lw;
}

// ****************************************************************************
void calcDoubleBondLines(const ROMol &mol, double offset, const Bond &bond,
                         const std::vector<Point2D> &at_cds, Point2D &l1s,
                         Point2D &l1f, Point2D &l2s, Point2D &l2f) {
  // the percent shorter that the extra bonds in a double bond are
  const double multipleBondTruncation = 0.15;
  Atom *at1 = bond.getBeginAtom();
  Atom *at2 = bond.getEndAtom();
  const Point2D &at1_cds = at_cds[at1->getIdx()];
  const Point2D &at2_cds = at_cds[at2->getIdx()];

  Point2D perp;
  if (1 == at1->getDegree() || 1 == at2->getDegree() ||
      isLinearAtom(*at1, at_cds) || isLinearAtom(*at2, at_cds)) {
    perp = calcPerpendicular(at1_cds, at2_cds) * offset;
    l1s = at1_cds + perp;
    l1f = at2_cds + perp;
    l2s = at1_cds - perp;
    l2f = at2_cds - perp;
  } else if ((Bond::EITHERDOUBLE == bond.getBondDir()) ||
             (Bond::STEREOANY == bond.getStereo())) {
    // crossed bond
    perp = calcPerpendicular(at1_cds, at2_cds) * offset;
    l1s = at1_cds + perp;
    l1f = at2_cds - perp;
    l2s = at1_cds - perp;
    l2f = at2_cds + perp;
  } else {
    l1s = at1_cds;
    l1f = at2_cds;
    offset *= 2.0;
    if (mol.getRingInfo()->numBondRings(bond.getIdx())) {
      // in a ring, we need to draw the bond inside the ring
      perp = bondInsideRing(mol, bond, at1_cds, at2_cds, at_cds);
    } else {
      perp = bondInsideDoubleBond(mol, bond, at_cds);
    }
    Point2D bv = at1_cds - at2_cds;
    l2s = at1_cds - bv * multipleBondTruncation + perp * offset;
    l2f = at2_cds + bv * multipleBondTruncation + perp * offset;
  }
}

// ****************************************************************************
void calcTripleBondLines(double offset, const Bond &bond,
                         const std::vector<Point2D> &at_cds,
                         Point2D &l1s, Point2D &l1f, Point2D &l2s,
                         Point2D &l2f) {
  // the percent shorter that the extra bonds in a double bond are
  const double multipleBondTruncation = 0.15;

  Atom *at1 = bond.getBeginAtom();
  Atom *at2 = bond.getEndAtom();
  const Point2D &at1_cds = at_cds[at1->getIdx()];
  const Point2D &at2_cds = at_cds[at2->getIdx()];

  // 2 lines, a bit shorter and offset on the perpendicular
  double dbo = 2.0 * offset;
  Point2D perp = calcPerpendicular(at1_cds, at2_cds);
  double end1_trunc = 1 == at1->getDegree() ? 0.0 : multipleBondTruncation;
  double end2_trunc = 1 == at2->getDegree() ? 0.0 : multipleBondTruncation;
  Point2D bv = at1_cds - at2_cds;
  l1s = at1_cds - (bv * end1_trunc) + perp * dbo;
  l1f = at2_cds + (bv * end2_trunc) + perp * dbo;
  l2s = at1_cds - (bv * end1_trunc) - perp * dbo;
  l2f = at2_cds + (bv * end2_trunc) - perp * dbo;
}

// ****************************************************************************
// calculate normalised perpendicular to vector between two coords
Point2D calcPerpendicular(const Point2D &cds1, const Point2D &cds2) {
  double bv[2] = {cds1.x - cds2.x, cds1.y - cds2.y};
  double perp[2] = {-bv[1], bv[0]};
  double perp_len = sqrt(perp[0] * perp[0] + perp[1] * perp[1]);
  perp[0] /= perp_len;
  perp[1] /= perp_len;

  return Point2D(perp[0], perp[1]);
}

// ****************************************************************************
// calculate normalised perpendicular to vector between two coords, such that
// it's inside the angle made between (1 and 2) and (2 and 3).
Point2D calcInnerPerpendicular(const Point2D &cds1, const Point2D &cds2,
                               const Point2D &cds3) {
  Point2D perp = calcPerpendicular(cds1, cds2);
  double v1[2] = {cds1.x - cds2.x, cds1.y - cds2.y};
  double v2[2] = {cds2.x - cds3.x, cds2.y - cds3.y};
  double obv[2] = {v1[0] - v2[0], v1[1] - v2[1]};

  // if dot product of centre_dir and perp < 0.0, they're pointing in opposite
  // directions, so reverse perp
  if (obv[0] * perp.x + obv[1] * perp.y < 0.0) {
    perp.x *= -1.0;
    perp.y *= -1.0;
  }

  return perp;
}

// ****************************************************************************
// cds1 and cds2 are 2 atoms in a ring.  Returns the perpendicular pointing
// into the ring
Point2D bondInsideRing(const ROMol &mol, const Bond &bond, const Point2D &cds1,
                       const Point2D &cds2,
                       const std::vector<Point2D> &at_cds) {
  std::vector<size_t> bond_in_rings;
  auto bond_rings = mol.getRingInfo()->bondRings();
  for (size_t i = 0; i < bond_rings.size(); ++i) {
    if (find(bond_rings[i].begin(), bond_rings[i].end(), bond.getIdx()) !=
        bond_rings[i].end()) {
      bond_in_rings.push_back(i);
    }
  }

  // find another bond in the ring connected to bond, use the
  // other end of it as the 3rd atom.
  auto calc_perp = [&](const Bond *bond, const INT_VECT &ring) -> Point2D * {
    Atom *bgn_atom = bond->getBeginAtom();
    for (const auto &nbri2 : make_iterator_range(mol.getAtomBonds(bgn_atom))) {
      const Bond *bond2 = mol[nbri2];
      if (bond2 == bond) {
        continue;
      }
      if (find(ring.begin(), ring.end(), bond2->getIdx()) != ring.end()) {
        int atom3 = bond2->getOtherAtomIdx(bond->getBeginAtomIdx());
        Point2D *ret = new Point2D;
        *ret = calcInnerPerpendicular(cds1, cds2, at_cds[atom3]);
        return ret;
      }
    }
    return nullptr;
  };

  if (bond_in_rings.size() > 1) {
    // bond is in more than 1 ring.  Choose one that is the same aromaticity
    // as the bond, so that if bond is aromatic, the double bond is inside
    // the aromatic ring.  This is important for morphine, for example,
    // where there are fused aromatic and aliphatic rings.
    // morphine: CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5
    for (size_t i = 0; i < bond_in_rings.size(); ++i) {
      auto ring = bond_rings[bond_in_rings[i]];
      bool ring_ok = true;
      for (auto bond_idx : ring) {
        const Bond *bond2 = mol.getBondWithIdx(bond_idx);
        if (bond.getIsAromatic() != bond2->getIsAromatic()) {
          ring_ok = false;
          break;
        }
      }
      if (!ring_ok) {
        continue;
      }
      Point2D *ret = calc_perp(&bond, ring);
      if (ret) {
        Point2D real_ret(*ret);
        delete ret;
        return real_ret;
      }
    }
  }

  // either bond is in 1 ring, or we couldn't decide above, so just use the
  // first one
  auto ring = bond_rings[bond_in_rings.front()];
  Point2D *ret = calc_perp(&bond, ring);
  if (ret) {
    Point2D real_ret(*ret);
    delete ret;
    return real_ret;
  }

  // failsafe that it will hopefully never see.
  return calcPerpendicular(cds1, cds2);
}

// ****************************************************************************
// cds1 and cds2 are 2 atoms in a chain double bond.  Returns the
// perpendicular pointing into the inside of the bond
Point2D bondInsideDoubleBond(const ROMol &mol, const Bond &bond,
                             const std::vector<Point2D> &at_cds) {
  // a chain double bond, where it looks nicer IMO if the 2nd line is inside
  // the angle of outgoing bond. Unless it's an allene, where nothing
  // looks great.
  const Atom *at1 = bond.getBeginAtom();
  const Atom *at2 = bond.getEndAtom();
  const Atom *bondAt, *endAtom;
  if (at1->getDegree() > 1) {
    bondAt = at1;
    endAtom = at2;
  } else {
    bondAt = at2;
    endAtom = at1;
  }
  int at3 = -1;  // to stop the compiler whinging.
  for (const auto &nbri2 : make_iterator_range(mol.getAtomBonds(bondAt))) {
    const Bond *bond2 = mol[nbri2];
    if (&bond != bond2) {
      at3 = bond2->getOtherAtomIdx(bondAt->getIdx());
      break;
    }
  }

  return calcInnerPerpendicular(at_cds[endAtom->getIdx()],
                                at_cds[bondAt->getIdx()], at_cds[at3]);
}

// ****************************************************************************
void adjustBondEndForString(
    const Point2D &end1, const Point2D &end2, double padding,
    const std::vector<std::shared_ptr<StringRect>> &rects,
    Point2D &moveEnd) {
  Point2D labelPos = moveEnd;
  for (auto r : rects) {
    Point2D origTrans = r->trans_;
    r->trans_ += labelPos;

    Point2D tl, tr, bl, br;
    r->calcCorners(tl, tr, br, bl, padding);

    // if it's a wide label, such as C:7, the bond can intersect
    // more than 1 side of the rectangle, so check them all.
    std::unique_ptr<Point2D> ip(new Point2D);
    if (MolDraw2D_detail::doLinesIntersect(moveEnd, end2, tl, tr, ip.get())) {
      moveEnd = *ip;
    }
    if (MolDraw2D_detail::doLinesIntersect(moveEnd, end2, tr, br, ip.get())) {
      moveEnd = *ip;
    }
    if (MolDraw2D_detail::doLinesIntersect(moveEnd, end2, br, bl, ip.get())) {
      moveEnd = *ip;
    }
    if (MolDraw2D_detail::doLinesIntersect(moveEnd, end2, bl, tl, ip.get())) {
      moveEnd = *ip;
    }
    r->trans_ = origTrans;
  }
}

// ****************************************************************************
void calcMolNotePosition(const std::vector<Point2D> atCds, DrawText &textDrawer,
                         AnnotationType &annot) {
  if (annot.text_.empty()) {
    annot.rect_.width_ = -1.0;  // so we know it's not valid.
    return;
  }

  std::vector<std::shared_ptr<StringRect>> rects;
  std::vector<TextDrawType> draw_modes;
  std::vector<char> draw_chars;

  // at this point, the scale() should still be 1, so min and max font sizes
  // don't make sense, as we're effectively operating on atom coords rather
  // than draw.
  double full_font_scale = textDrawer.fontScale();
  textDrawer.setFontScale(1, true);
  textDrawer.getStringRects(annot.text_, OrientType::N, rects, draw_modes,
                               draw_chars);
  textDrawer.setFontScale(full_font_scale, true);
  // accumulate the widths of the rectangles so that we have the overall width
  for (const auto &rect : rects) {
    annot.rect_.width_ += rect->width_;
  }

  Point2D centroid{0., 0.};
  double minY = std::numeric_limits<double>::max();
  double maxX = -std::numeric_limits<double>::max();
  for (const auto &pt : atCds) {
    centroid += pt;
    maxX = std::max(pt.x, maxX);
    minY = std::min(pt.y, minY);
  }
  centroid /= atCds.size();

  // because we've inverted the Y coord, we need to use -minY, not +maxY.
  Point2D vect{maxX, -minY};
  vect.x -= centroid.x;
  vect.y += centroid.y;
  auto loc = centroid + vect * 0.9;
  loc = centroid;
  loc.x += vect.x * 0.9;
  loc.y -= vect.y * 0.9;
  annot.rect_.trans_ = loc;
  std::cout << "mol note pos : " << annot.rect_.trans_ << "   and dims "
            << annot.rect_.width_ << " by " << annot.rect_.height_ << std::endl;
}

// ****************************************************************************
void findAnnotationExtremes(const std::vector<AnnotationType> &annots,
                            double &xmin, double &xmax, double &ymin,
                            double &ymax) {
  for (auto const &pr : annots) {
    findRectExtremes(pr.rect_, pr.align_, xmin, xmax, ymin, ymax);
  }
}

// ****************************************************************************
void findRadicalExtremes(
    const std::vector<std::pair<StringRect, OrientType>> &radicals,
    double &xmin, double &xmax, double &ymin, double &ymax) {
  for (auto const &rad : radicals) {
    findRectExtremes(rad.first, TextAlignType::MIDDLE, xmin, xmax, ymin, ymax);
  }
}

// ****************************************************************************
void findRectExtremes(const StringRect &rect, const TextAlignType &align,
                      double &xmin, double &xmax, double &ymin, double &ymax){
  double this_xmax = rect.trans_.x;
  double this_xmin = rect.trans_.x;
  double this_ymax = rect.trans_.y;
  double this_ymin = rect.trans_.y;
  if (align == TextAlignType::START) {
    this_xmax += rect.width_;
  } else if (align == TextAlignType::END) {
    this_xmin -= rect.width_;
  } else {
    this_xmax += rect.width_ / 2.0;
    this_xmin -= rect.width_ / 2.0;
  }
  this_ymax += rect.height_ / 2.0;
  this_ymin -= rect.height_ / 2.0;

  xmax = std::max(xmax, this_xmax);
  xmin = std::min(xmin, this_xmin);
  ymax = std::max(ymax, this_ymax);
  ymin = std::min(ymin, this_ymin);
}

} // namespace RDKit