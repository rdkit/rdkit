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

#include <algorithm>
#include <iostream>
#include <limits>

#include <Geometry/Transform2D.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/FileParsers/MolSGroupParsing.h>
#include <GraphMol/MolDraw2D/AtomSymbol.h>
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawShape.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolEnumerator/LinkNode.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Atropisomers.h>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawMol::DrawMol(
    const ROMol &mol, const std::string &legend, int width, int height,
    const MolDrawOptions &drawOptions, DrawText &textDrawer,
    const std::vector<int> *highlight_atoms,
    const std::vector<int> *highlight_bonds,
    const std::map<int, DrawColour> *highlight_atom_map,
    const std::map<int, DrawColour> *highlight_bond_map,
    const std::vector<std::pair<DrawColour, DrawColour>> *bond_colours,
    const std::map<int, double> *highlight_radii, bool includeAnnotations,
    int confId, bool isReactionMol)
    : drawOptions_(drawOptions),
      textDrawer_(textDrawer),
      marginPadding_(drawOptions.padding),
      includeAnnotations_(includeAnnotations),
      isReactionMol_(isReactionMol),
      legend_(legend),
      confId_(confId),
      width_(width),
      height_(height),
      drawWidth_(width),
      drawHeight_(height),
      scale_(1.0),
      fontScale_(1.0),
      xMin_(std::numeric_limits<double>::max() / 2.0),
      yMin_(std::numeric_limits<double>::max() / 2.0),
      xMax_(std::numeric_limits<double>::lowest() / 2.0),
      yMax_(std::numeric_limits<double>::lowest() / 2.0),
      xRange_(std::numeric_limits<double>::max()),
      yRange_(std::numeric_limits<double>::max()),
      flexiCanvasX_(width < 0.0),
      flexiCanvasY_(height < 0.0) {
  if (highlight_atoms) {
    highlightAtoms_ = *highlight_atoms;
  }
  if (highlight_bonds) {
    highlightBonds_ = *highlight_bonds;
  }
  if (highlight_atom_map) {
    highlightAtomMap_ = *highlight_atom_map;
  }
  if (highlight_bond_map) {
    highlightBondMap_ = *highlight_bond_map;
  }
  if (bond_colours) {
    bondColours_ = *bond_colours;
  }
  if (highlight_radii) {
    highlightRadii_ = *highlight_radii;
  }
  textDrawer_.setFontScale(fontScale_, true);
  initDrawMolecule(mol);
}

// ****************************************************************************
DrawMol::DrawMol(int width, int height, const MolDrawOptions &drawOptions,
                 DrawText &textDrawer, double xmin, double xmax, double ymin,
                 double ymax, double scale, double fontscale)
    : drawOptions_(drawOptions),
      textDrawer_(textDrawer),
      marginPadding_(drawOptions.padding),
      isReactionMol_(false),
      confId_(-1),
      width_(width),
      height_(height),
      drawWidth_(width),
      drawHeight_(height),
      scale_(scale),
      fontScale_(fontscale),
      xMin_(xmin),
      yMin_(ymin),
      xMax_(xmax),
      yMax_(ymax),
      xRange_(xmax - xmin),
      yRange_(ymax - ymin),
      molHeight_(height) {
  textDrawer_.setFontScale(fontScale_, true);
  // we reverse the y coords of everything, so do that here, too
  yMin_ *= -1;
  yMax_ *= -1;
  std::swap(yMin_, yMax_);
}

// ****************************************************************************
void DrawMol::createDrawObjects() {
  textDrawer_.setFontScale(fontScale_, true);
  partitionForLegend();
  extractAll(scale_);
  calculateScale();

  bool ignoreFontLimits = drawOptions_.fixedFontSize != -1;
  if (!textDrawer_.setFontScale(fontScale_, ignoreFontLimits) ||
      ignoreFontLimits) {
    // in either of these cases, the relative font size isn't what we were
    // expecting, so we need to rebuild everything.

    // furthermore, if it's a fully flexible canvas and the font scale is
    // greater than the global scale, if there are characters at the edge
    // of the image, the canvas won't be big enough (Github6111). Rebuild
    // with an appropriate relative font size.
    if (flexiCanvasX_ && flexiCanvasY_ && (fontScale_ - scale_) > 1e-4) {
      width_ = -1;
      height_ = -1;
      auto currScale = textDrawer_.fontScale();
      auto relScale = fontScale_ / scale_;
      resetEverything();
      fontScale_ = relScale;
      textDrawer_.setFontScale(relScale, true);
      extractAll(scale_);
      calculateScale();
      textDrawer_.setFontScale(currScale, true);
    }
    setScale(scale_, textDrawer_.fontScale(), ignoreFontLimits);
  } else {
    finishCreateDrawObjects();
  }
}

// ****************************************************************************
void DrawMol::finishCreateDrawObjects() {
  // the legend and mol notes need the final scale to get the fonts the
  // correct size.
  extractLegend();
  changeToDrawCoords();
  // these need the draw coords.
  extractMolNotes();
  extractCloseContacts();
  drawingInitialised_ = true;
}

// ****************************************************************************
void DrawMol::initDrawMolecule(const ROMol &mol) {
  drawMol_.reset(new RWMol(mol));
  if (drawOptions_.centreMoleculesBeforeDrawing) {
    if (drawMol_->getNumConformers()) {
      centerMolForDrawing(*drawMol_, confId_);
    }
  }
  if (drawOptions_.unspecifiedStereoIsUnknown) {
    markUnspecifiedStereoAsUnknown(*drawMol_);
  }
  if (drawOptions_.useMolBlockWedging) {
    Chirality::reapplyMolBlockWedging(*drawMol_);
  }
  if (!isReactionMol_) {
    if (drawOptions_.prepareMolsBeforeDrawing) {
      MolDraw2DUtils::prepareMolForDrawing(*drawMol_);
    } else if (!mol.getNumConformers()) {
      const bool canonOrient = true;
      RDDepict::compute2DCoords(*drawMol_, nullptr, canonOrient);
    }
  }
  if (drawOptions_.simplifiedStereoGroupLabel &&
      !mol.hasProp(common_properties::molNote)) {
    bool removeAffectedStereoGroups = true;
    Chirality::simplifyEnhancedStereo(*drawMol_, removeAffectedStereoGroups);
  }
  if (drawOptions_.addStereoAnnotation) {
    Chirality::addStereoAnnotations(*drawMol_);
  }
  if (drawOptions_.addAtomIndices) {
    addAtomIndices(*drawMol_);
  }
  if (drawOptions_.addBondIndices) {
    addBondIndices(*drawMol_);
  }
}

// ****************************************************************************
void DrawMol::extractAll(double scale) {
  extractAtomCoords();
  extractAtomSymbols();
  // extractVariableBonds removes the * symbol from the end of dative bonds
  // that are showing a haptic bond, so this needs to be done before the
  // bonds are extracted, or it will shorten the bond so as not to clash with
  // the * which won't be in the final picture.
  extractVariableBonds();
  extractBonds();
  extractRegions();
  extractHighlights(scale);
  extractAttachments();
  extractAtomNotes();
  if (!drawOptions_.addStereoAnnotation) {
    extractStereoGroups();
  }
  extractBondNotes();
  extractRadicals();
  extractSGroupData();
  extractBrackets();
  extractLinkNodes();
}

// ****************************************************************************
void DrawMol::extractAtomCoords() {
  PRECONDITION(static_cast<int>(drawMol_->getNumConformers()) > 0, "no coords");

  const RDGeom::POINT3D_VECT &locs =
      drawMol_->getConformer(confId_).getPositions();

  // the transformation rotates anti-clockwise, as is conventional, but
  // probably not what our user expects.  But because we invert the y
  // coord, it needs to be applied in a positive direction.
  double rot = drawOptions_.rotate * M_PI / 180.0;
  // assuming that if drawOptions_.rotate is set to 0.0, rot will be
  // exactly 0.0 without worrying about floating point number dust.
  //
  // NB - the y coord is inverted, so that the molecule coords and the
  // draw coords always go down the page as y increases.  This is so
  // we can have all the draw entities in molecule coords or draw coords
  // and min and max y will be going in the same direction.
  RDGeom::Transform2D trans;
  trans.SetTransform(Point2D(0.0, 0.0), rot);
  atCds_.clear();
  for (auto pt3 : locs) {
    Point2D pt{pt3.x, -pt3.y};
    if (rot != 0.0) {
      trans.TransformPoint(pt);
    }
    atCds_.push_back(pt);
  }
}

// ****************************************************************************
void DrawMol::extractAtomSymbols() {
  atomicNums_.clear();
  for (auto at1 : drawMol_->atoms()) {
    if (!isComplexQuery(at1)) {
      atomicNums_.push_back(at1->getAtomicNum());
    } else {
      atomicNums_.push_back(0);
    }
    std::pair<std::string, OrientType> atSym =
        getAtomSymbolAndOrientation(*at1);
    atomSyms_.push_back(atSym);
    if (!atSym.first.empty()) {
      DrawColour atCol = getColour(at1->getIdx());
      AtomSymbol *al =
          new AtomSymbol(atSym.first, at1->getIdx(), atSym.second,
                         atCds_[at1->getIdx()], atCol, textDrawer_);
      atomLabels_.emplace_back(al);
    } else {
      atomLabels_.emplace_back(nullptr);
    }
  }
}

// ****************************************************************************
void DrawMol::extractBonds() {
  calcMeanBondLength();
  double doubleBondOffset = drawOptions_.multipleBondOffset * meanBondLength_;

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
  adjustBondsOnSolidWedgeEnds();
  smoothBondJoins();
}

// ****************************************************************************
void DrawMol::extractHighlights(double scale) {
  if (drawOptions_.continuousHighlight) {
    makeContinuousHighlights(scale);
  } else {
    if (drawOptions_.circleAtoms && !highlightAtoms_.empty()) {
      makeAtomCircleHighlights();
    }
  }
}

// ****************************************************************************
void DrawMol::extractRegions() {
  for (const auto &region : drawOptions_.atomRegions) {
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
      highlights_.emplace_back(pl);
    }
  }
}

// ****************************************************************************
void DrawMol::extractAttachments() {
  if (drawOptions_.dummiesAreAttachments) {
    for (const auto at1 : drawMol_->atoms()) {
      if (at1->hasProp(common_properties::atomLabel) ||
          drawOptions_.atomLabels.find(at1->getIdx()) !=
              drawOptions_.atomLabels.end()) {
        // skip dummies that explicitly have a label provided
        continue;
      }
      if (at1->getAtomicNum() == 0 && at1->getDegree() == 1) {
        Point2D &at1_cds = atCds_[at1->getIdx()];
        const auto &iter_pair = drawMol_->getAtomNeighbors(at1);
        const Atom *at2 = (*drawMol_)[*iter_pair.first];
        Point2D &at2_cds = atCds_[at2->getIdx()];
        Point2D perp = calcPerpendicular(at1_cds, at2_cds);
        Point2D p1 =
            Point2D(at1_cds.x - perp.x * 0.5, at1_cds.y - perp.y * 0.5);
        Point2D p2 =
            Point2D(at1_cds.x + perp.x * 0.5, at1_cds.y + perp.y * 0.5);
        DrawColour col(.5, .5, .5);
        std::vector<Point2D> points{p1, p2};
        double offset = drawOptions_.multipleBondOffset * meanBondLength_ / 2.0;
        DrawShapeWavyLine *wl = new DrawShapeWavyLine(
            points, drawOptions_.bondLineWidth, false, col, col, offset,
            at2->getIdx() + activeAtmIdxOffset_);
        bonds_.emplace_back(wl);
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
    // molecule annotations use a full-size font, hence the 1 below.
    DrawAnnotation tmp(note, TextAlignType::START, "note", 1, Point2D(0.0, 0.0),
                       drawOptions_.annotationColour, textDrawer_);
    double height, width;
    tmp.getDimensions(width, height);
    // Try all 4 corners until there's no clash with the underlying molecule.
    // Even though alignment is START, the DrawAnnotation puts the middle
    // of the first char at the location, so that needs to be adjusted for.
    std::vector<Point2D> locs = {
        {width_ - width, height},
        {0.0 + tmp.rects_[0]->width_ / 2.0, height},
        {0.0 + tmp.rects_[0]->width_ / 2.0, double(drawHeight_ - height)},
        {width_ - width, double(drawHeight_ - height)},
    };
    bool didIt = false;
    for (int i = 0; i < 3; ++i) {
      locs[i].x += xOffset_;
      locs[i].y += yOffset_;
      auto annot = std::make_unique<DrawAnnotation>(
          note, TextAlignType::START, "note", 1.0, locs[i],
          drawOptions_.annotationColour, textDrawer_);
      // Put it into the legends_, because it's already in draw coords, so
      // shouldn't be treated by changeToDrawCoords.
      if (!doesNoteClash(*annot)) {
        legends_.push_back(std::move(annot));
        didIt = true;
        break;
      }
    }
    if (!didIt) {
      // There was nowhere to put it that didn't clash, so live with it.
      legends_.emplace_back(
          new DrawAnnotation(note, TextAlignType::START, "note", 1.0, locs[0],
                             drawOptions_.annotationColour, textDrawer_));
    }
  }
}

// ****************************************************************************
void DrawMol::extractAtomNotes() {
  for (const auto atom : drawMol_->atoms()) {
    std::string note;
    if (atom->getPropIfPresent(common_properties::atomNote, note)) {
      if (!note.empty()) {
        DrawAnnotation *annot = new DrawAnnotation(
            note, TextAlignType::MIDDLE, "note",
            drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
            drawOptions_.annotationColour, textDrawer_);
        calcAnnotationPosition(atom, *annot);
        annotations_.emplace_back(annot);
      }
    }
  }
}

// ****************************************************************************
void DrawMol::extractStereoGroups() {
  int orCount(0), andCount(0);
  for (const StereoGroup &group : drawMol_->getStereoGroups()) {
    std::string stereoGroupType;

    switch (group.getGroupType()) {
      case RDKit::StereoGroupType::STEREO_ABSOLUTE:
        stereoGroupType = "abs";
        break;
      case RDKit::StereoGroupType::STEREO_OR:
        stereoGroupType = "or" + std::to_string(++orCount);
        break;
      case RDKit::StereoGroupType::STEREO_AND:
        stereoGroupType = "and" + std::to_string(++andCount);
        break;
      default:
        throw ValueErrorException("Unrecognized stereo group type");
    }

    std::vector<unsigned int> atomIds;
    std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        wedgeBonds;  // empty - all wedges should have been added to the mol, so
                     // this doesn't matter
    Atropisomers::getAllAtomIdsForStereoGroup(*drawMol_, group, atomIds,
                                              wedgeBonds);

    for (auto atomId : atomIds) {
      DrawAnnotation *annot = new DrawAnnotation(
          stereoGroupType, TextAlignType::MIDDLE, "stereoGroup",
          drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
          drawOptions_.annotationColour, textDrawer_);
      calcAnnotationPosition(drawMol_->getAtomWithIdx(atomId), *annot);
      annotations_.emplace_back(annot);
    }
  }
}

// ****************************************************************************
void DrawMol::extractBondNotes() {
  for (const auto bond : drawMol_->bonds()) {
    std::string note;
    if (bond->getPropIfPresent(common_properties::bondNote, note)) {
      if (!note.empty()) {
        DrawAnnotation *annot = new DrawAnnotation(
            note, TextAlignType::MIDDLE, "note",
            drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
            drawOptions_.annotationColour, textDrawer_);
        calcAnnotationPosition(bond, *annot);
        annotations_.emplace_back(annot);
      }
    }
  }
}

// ****************************************************************************
void DrawMol::extractRadicals() {
  if (!drawOptions_.includeRadicals) {
    return;
  }
  for (const auto atom : drawMol_->atoms()) {
    if (!atom->getNumRadicalElectrons()) {
      continue;
    }
    StringRect rad_rect;
    OrientType orient = calcRadicalRect(atom, rad_rect);
    radicals_.emplace_back(rad_rect, orient, atom->getIdx());
  }
}

// ****************************************************************************
void DrawMol::extractSGroupData() {
  if (!includeAnnotations_) {
    return;
  }
  const auto &sgs = getSubstanceGroups(*drawMol_);
  if (sgs.empty()) {
    return;
  }

  // details of this transformation are in extractAtomCoords
  double rot = drawOptions_.rotate * M_PI / 180.0;
  RDGeom::Transform2D tform;
  tform.SetTransform(Point2D(0.0, 0.0), rot);

  for (const auto &sg : sgs) {
    std::string typ;
    if (sg.getPropIfPresent("TYPE", typ) && typ == "DAT") {
      std::string text;
      // it seems like we should be rendering FIELDNAME, but
      // Marvin Sketch, Biovia Draw, and ChemDraw don't do it
      // if (sg.getPropIfPresent("FIELDNAME", text)) {
      //   text += "=";
      // };
      if (sg.hasProp("DATAFIELDS")) {
        STR_VECT dfs = sg.getProp<STR_VECT>("DATAFIELDS");
        for (const auto &df : dfs) {
          text += df + "|";
        }
        text.pop_back();
      }
      if (text.empty()) {
        continue;
      }
      int atomIdx = -1;
      if (!sg.getAtoms().empty()) {
        atomIdx = sg.getAtoms()[0];
      };
      bool located = false;
      std::string fieldDisp;
      Point2D origLoc(0.0, 0.0);
      if (sg.getPropIfPresent("FIELDDISP", fieldDisp)) {
        double xp = FileParserUtils::stripSpacesAndCast<double>(
            fieldDisp.substr(0, 10));
        double yp = FileParserUtils::stripSpacesAndCast<double>(
            fieldDisp.substr(10, 10));
        // we always invert y for the molecule coords
        origLoc = Point2D{xp, -yp};

        if (fieldDisp[25] == 'R') {
          if (atomIdx < 0) {
            // we will warn about this below
            text = "";
          } else if (fabs(xp) > 1e-3 || fabs(yp) > 1e-3) {
            // opposite sign for y
            origLoc.x += drawMol_->getConformer().getAtomPos(atomIdx).x;
            origLoc.y -= drawMol_->getConformer().getAtomPos(atomIdx).y;
            located = true;
          }
        } else {
          if (drawMol_->hasProp("_centroidx")) {
            Point2D centroid;
            drawMol_->getProp("_centroidx", centroid.x);
            drawMol_->getProp("_centroidy", centroid.y);
            // opposite sign for y
            origLoc.x += centroid.x;
            origLoc.y -= centroid.y;
          }
          located = true;
        }
        tform.TransformPoint(origLoc);
      }

      if (!text.empty()) {
        // looks like everybody renders these left justified
        DrawAnnotation *annot = new DrawAnnotation(
            text, TextAlignType::START, "note",
            drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
            drawOptions_.annotationColour, textDrawer_);
        if (!located) {
          if (atomIdx >= 0 && !text.empty()) {
            calcAnnotationPosition(drawMol_->getAtomWithIdx(atomIdx), *annot);
          }
        } else {
          annot->pos_ = origLoc;
        }
        annotations_.emplace_back(annot);
      } else {
        BOOST_LOG(rdWarningLog)
            << "FIELDDISP info not found for DAT SGroup which isn't "
               "associated with an atom. SGroup will not be rendered."
            << std::endl;
      }
    }
  }
}

// ****************************************************************************
void DrawMol::extractVariableBonds() {
  boost::dynamic_bitset<> atomsInvolved(drawMol_->getNumAtoms());
  for (const auto bond : drawMol_->bonds()) {
    std::string endpts;
    std::string attach;
    if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, endpts) &&
        bond->getPropIfPresent(common_properties::_MolFileBondAttach, attach)) {
      // FIX: maybe distinguish between "ANY" and "ALL" values of attach here?
      std::vector<unsigned int> oats =
          RDKit::SGroupParsing::ParseV3000Array<unsigned int>(endpts);
      atomsInvolved.reset();
      // decrement the indices and do error checking:
      for (auto &oat : oats) {
        if (oat == 0 || oat > drawMol_->getNumAtoms()) {
          throw ValueErrorException("Bad variation point index");
        }
        --oat;
        atomsInvolved.set(oat);
        auto center = atCds_[oat];
        Point2D offset{drawOptions_.variableAtomRadius,
                       drawOptions_.variableAtomRadius};
        std::vector<Point2D> points{center, offset};
        DrawShapeEllipse *ell = new DrawShapeEllipse(
            points, 1, true, drawOptions_.variableAttachmentColour, true, oat);
        preShapes_.emplace_back(ell);
      }

      for (const auto bond : drawMol_->bonds()) {
        if (atomsInvolved[bond->getBeginAtomIdx()] &&
            atomsInvolved[bond->getEndAtomIdx()]) {
          std::vector<Point2D> points{atCds_[bond->getBeginAtomIdx()],
                                      atCds_[bond->getEndAtomIdx()]};
          DrawShapeSimpleLine *sl = new DrawShapeSimpleLine(
              points, drawOptions_.variableBondWidthMultiplier, true,
              drawOptions_.variableAttachmentColour,
              bond->getBeginAtomIdx() + activeAtmIdxOffset_,
              bond->getEndAtomIdx() + activeAtmIdxOffset_,
              bond->getIdx() + activeBndIdxOffset_);
          preShapes_.emplace_back(sl);
        }
      }
      // correct the symbol of the end atom (remove the *):
      if (!bond->getBeginAtom()->getAtomicNum()) {
        atomSyms_[bond->getBeginAtomIdx()] = std::make_pair("", OrientType::C);
        atomLabels_[bond->getBeginAtomIdx()].reset();
      }
    }
  }
}

// ****************************************************************************
namespace {
// function to draw a label at the bottom of the bracket.
DrawAnnotation *drawBottomLabel(const std::string &label,
                                const DrawShape &brkShp,
                                const MolDrawOptions &drawOptions,
                                DrawText &textDrawer, bool horizontal) {
  // annotations go on the last bracket of an sgroup
  // LABEL goes at the bottom which is now the top
  auto topPt = brkShp.points_[1];
  auto brkPt = brkShp.points_[0];
  if ((!horizontal && brkShp.points_[2].y > topPt.y) ||
      (horizontal && brkShp.points_[2].x < topPt.x)) {
    topPt = brkShp.points_[2];
    brkPt = brkShp.points_[3];
  }
  DrawAnnotation *da = new DrawAnnotation(
      label, TextAlignType::MIDDLE, "connect", drawOptions.annotationFontScale,
      topPt + (topPt - brkPt), DrawColour(0.0, 0.0, 0.0), textDrawer);
  if (brkPt.x < topPt.x) {
    da->align_ = TextAlignType::START;
  }
  return da;
}
}  // namespace

void DrawMol::extractBrackets() {
  auto &sgs = getSubstanceGroups(*drawMol_);
  if (sgs.empty()) {
    return;
  }
  // details of this transformation are in extractAtomCoords
  double rot = drawOptions_.rotate * M_PI / 180.0;
  RDGeom::Transform2D trans;
  trans.SetTransform(Point2D(0.0, 0.0), rot);
  for (auto &sg : sgs) {
    if (sg.getBrackets().empty()) {
      continue;
    }
    // figure out the location of the reference point we'll use to figure out
    // which direction the bracket points
    // Thanks to John Mayfield for the thoughts on the best way to do this:
    //   http://efficientbits.blogspot.com/2015/11/bringing-molfile-sgroups-to-cdk.html
    Point2D refPt{0., 0.};

    if (!sg.getAtoms().empty()) {
      // use the average position of the atoms in the sgroup
      // Github5768 shows that this is a bit simplistic in some cases.  In
      // that molecule, there is a long chain that stretches outside the
      // bracket area that turns the last bracket the wrong way.
      // Just pick out the SGroup atoms that are inside brackets, rather
      // crudely.
      double xMin = std::numeric_limits<double>::max() / 2.0;
      double yMin = std::numeric_limits<double>::max() / 2.0;
      double xMax = std::numeric_limits<double>::lowest() / 2.0;
      double yMax = std::numeric_limits<double>::lowest() / 2.0;
      for (const auto &brk : sg.getBrackets()) {
        Point2D p1{brk[0].x, -brk[0].y};
        Point2D p2{brk[1].x, -brk[1].y};
        trans.TransformPoint(p1);
        trans.TransformPoint(p2);
        xMin = std::min({xMin, p1.x, p2.x});
        yMin = std::min({yMin, p1.y, p2.y});
        xMax = std::max({xMax, p1.x, p2.x});
        yMax = std::max({yMax, p1.y, p2.y});
      }

      int numIn = 0;
      for (auto aidx : sg.getAtoms()) {
        if (atCds_[aidx].x >= xMin && atCds_[aidx].x <= xMax &&
            atCds_[aidx].y >= yMin && atCds_[aidx].y <= yMax) {
          refPt += atCds_[aidx];
          ++numIn;
        }
      }
      if (numIn) {
        refPt /= numIn;
      } else {
        // we'll have to go with all of them, and live with the consequences
        for (auto aidx : sg.getAtoms()) {
          refPt += atCds_[aidx];
        }
        refPt /= sg.getAtoms().size();
      }
    }

    std::vector<std::pair<Point2D, Point2D>> sgBondSegments;
    for (auto bndIdx : sg.getBonds()) {
      const auto bnd = drawMol_->getBondWithIdx(bndIdx);
      if (std::find(sg.getAtoms().begin(), sg.getAtoms().end(),
                    bnd->getBeginAtomIdx()) != sg.getAtoms().end()) {
        sgBondSegments.push_back(std::make_pair(atCds_[bnd->getBeginAtomIdx()],
                                                atCds_[bnd->getEndAtomIdx()]));

      } else if (std::find(sg.getAtoms().begin(), sg.getAtoms().end(),
                           bnd->getEndAtomIdx()) != sg.getAtoms().end()) {
        sgBondSegments.push_back(std::make_pair(
            atCds_[bnd->getEndAtomIdx()], atCds_[bnd->getBeginAtomIdx()]));
      }
    }
    int numBrackets = 0;
    for (const auto &brk : sg.getBrackets()) {
      // the atom coords have been inverted in y, so the bracket coords
      // must be, too.
      ++numBrackets;
      Point2D p1{brk[0].x, -brk[0].y};
      Point2D p2{brk[1].x, -brk[1].y};
      trans.TransformPoint(p1);
      trans.TransformPoint(p2);
      auto points = getBracketPoints(p1, p2, refPt, sgBondSegments);
      DrawShapePolyLine *pl =
          new DrawShapePolyLine(points, drawOptions_.bondLineWidth, false,
                                DrawColour(0.0, 0.0, 0.0), false);
      postShapes_.emplace_back(pl);
    }
    if (includeAnnotations_) {
      // Find the bottom-most or right-most bracket.  First work out if the
      // bracket is largely horizontal or largely vertical.
      const auto &brkShp = *postShapes_.back();
      Point2D longline = brkShp.points_[1] - brkShp.points_[2];
      longline.normalize();
      static const double cos45 = 1.0 / sqrt(2.0);
      bool horizontal = fabs(longline.x) > cos45;
      size_t labelBrk = postShapes_.size() - 1;
      for (int i = 1; i < numBrackets; ++i) {
        const auto &brkShp = *postShapes_[postShapes_.size() - i - 1];
        if (horizontal) {
          if (brkShp.points_[2].y > postShapes_[labelBrk]->points_[2].y) {
            labelBrk = postShapes_.size() - i - 1;
          }
        } else {
          if (brkShp.points_[2].x > postShapes_[labelBrk]->points_[2].x) {
            labelBrk = postShapes_.size() - i - 1;
          }
        }
      }
      std::string connect;
      if (sg.getPropIfPresent("CONNECT", connect)) {
        // annotations go on the last bracket of an sgroup
        const auto &brkShp = *postShapes_[labelBrk];
        // CONNECT goes at the top, but that's now the bottom due to the y
        // inversion
        auto botPt = brkShp.points_[2];
        auto brkPt = brkShp.points_[3];
        if ((!horizontal && brkShp.points_[1].y < botPt.y) ||
            (horizontal && brkShp.points_[1].x > botPt.x)) {
          botPt = brkShp.points_[1];
          brkPt = brkShp.points_[0];
        }
        DrawAnnotation *da = new DrawAnnotation(
            connect, TextAlignType::MIDDLE, "connect",
            drawOptions_.annotationFontScale, botPt + (botPt - brkPt),
            DrawColour(0.0, 0.0, 0.0), textDrawer_);
        // if we're to the right of the bracket, we need to left justify,
        // otherwise things seem to work as is
        if (brkPt.x < botPt.x) {
          da->align_ = TextAlignType::START;
        }
        annotations_.emplace_back(da);
      }

      std::string label;
      if (sg.getPropIfPresent("LABEL", label)) {
        auto da = drawBottomLabel(label, *postShapes_[labelBrk], drawOptions_,
                                  textDrawer_, horizontal);
        annotations_.emplace_back(da);
      } else if (sg.getPropIfPresent("TYPE", label)) {
        if (label == "GEN") {
          // ChemDraw doesn't draw the GEN (type=generic) label.
          continue;
        }
        // draw the lowercase type if there's no label to go there.
        std::transform(label.begin(), label.end(), label.begin(), ::tolower);
        auto da = drawBottomLabel(label, *postShapes_[labelBrk], drawOptions_,
                                  textDrawer_, horizontal);
        annotations_.emplace_back(da);
      }
    }
  }
}

// ****************************************************************************
void DrawMol::extractLinkNodes() {
  if (!drawMol_->hasProp(common_properties::molFileLinkNodes)) {
    return;
  }

  bool strict = false;
  auto linkNodes = MolEnumerator::utils::getMolLinkNodes(*drawMol_, strict);
  for (const auto &node : linkNodes) {
    const double crossingFrac = 0.333;
    const double lengthFrac = 0.333;
    Point2D labelPt{-1000, -1000};
    Point2D labelPerp{0, 0};
    for (const auto &bAts : node.bondAtoms) {
      // unlike brackets, we know how these point
      Point2D startLoc = atCds_[bAts.first];
      Point2D endLoc = atCds_[bAts.second];
      auto vect = endLoc - startLoc;
      auto offset = vect * crossingFrac;
      auto crossingPt = startLoc + offset;
      Point2D perp{vect.y, -vect.x};
      perp *= lengthFrac;
      Point2D p1 = crossingPt + perp / 2.;
      Point2D p2 = crossingPt - perp / 2.;

      std::vector<std::pair<Point2D, Point2D>> bondSegments;  // not needed here
      std::vector<Point2D> points{
          getBracketPoints(p1, p2, startLoc, bondSegments)};
      DrawShapePolyLine *pl =
          new DrawShapePolyLine(points, drawOptions_.bondLineWidth, false,
                                DrawColour(0.0, 0.0, 0.0), false);
      postShapes_.emplace_back(pl);

      if (p1.x > labelPt.x) {
        labelPt = p1;
        labelPerp = crossingPt - startLoc;
      }
      if (p2.x > labelPt.x) {
        labelPt = p2;
        labelPerp = crossingPt - startLoc;
      }
    }

    // the label
    if (includeAnnotations_) {
      std::string label =
          (boost::format("(%d-%d)") % node.minRep % node.maxRep).str();
      Point2D perp = labelPerp;
      perp /= perp.length() * 5;
      DrawAnnotation *da =
          new DrawAnnotation(label, TextAlignType::START, "linknode",
                             drawOptions_.annotationFontScale, labelPt + perp,
                             DrawColour(0.0, 0.0, 0.0), textDrawer_);
      annotations_.emplace_back(da);
    }
  }
}

// ****************************************************************************
void DrawMol::extractCloseContacts() {
  if (drawOptions_.flagCloseContactsDist < 0) {
    return;
  }
  int tol =
      drawOptions_.flagCloseContactsDist * drawOptions_.flagCloseContactsDist;
  boost::dynamic_bitset<> flagged(atCds_.size());
  Point2D trans, scale, toCentre;
  getDrawTransformers(trans, scale, toCentre);
  for (unsigned int i = 0; i < atCds_.size(); ++i) {
    if (flagged[i]) {
      continue;
    }
    Point2D ci = transformPoint(atCds_[i], &trans, &scale, &toCentre);
    for (unsigned int j = i + 1; j < atCds_.size(); ++j) {
      if (flagged[j]) {
        continue;
      }
      Point2D cj = transformPoint(atCds_[j], &trans, &scale, &toCentre);
      double d = (cj - ci).lengthSq();
      if (d <= tol) {
        flagged.set(i);
        flagged.set(j);
        break;
      }
    }
    if (flagged[i]) {
      Point2D p1 = ci;
      Point2D p2 = p1;
      Point2D offset(0.1 * scale_, 0.1 * scale_);
      p1 -= offset;
      p2 += offset;
      std::vector<Point2D> points(5);
      points[0] = points[4] = p1;
      points[1] = Point2D{p1.x, p2.y};
      points[2] = Point2D{p2};
      points[3] = Point2D{p2.x, p1.y};
      DrawShapePolyLine *pl =
          new DrawShapePolyLine(points, drawOptions_.bondLineWidth, false,
                                DrawColour(1.0, 0.0, 0.0), false, i);
      postShapes_.emplace_back(pl);
    };
  }
}

// ****************************************************************************
void DrawMol::calculateScale() {
  findExtremes();

  // if width < 0, we'll take the scale off the yRange_, and likewise with
  // height and xRange_.  If both are negative, use drawOptions_scalingFactor.
  double newScale = 1.0;
  if (width_ < 0 && height_ < 0) {
    width_ = drawOptions_.scalingFactor * xRange_ * (1 + 2 * marginPadding_);
    molHeight_ =
        drawOptions_.scalingFactor * yRange_ * (1 + 2 * marginPadding_);
  } else if (width_ < 0 && yRange_ > 1.0e-4) {
    newScale = double(height_) / yRange_;
    // if the molecule is very wide and short (e.g. HO-NH2) don't let the
    // bonds get too long.
    double mbl = meanBondLength_ * newScale;
    if (mbl > molHeight_ / 2) {
      newScale *= (molHeight_ / 2) / mbl;
    }
    width_ = newScale * xRange_;
  } else if (height_ < 0 && xRange_ > 1.0e-4) {
    newScale = double(width_) / xRange_;
    double mbl = meanBondLength_ * newScale;
    if (mbl > width_ / 2) {
      newScale *= (width_ / 2) / mbl;
    }
    molHeight_ = newScale * yRange_;
  }
  if (height_ < 0) {
    height_ = molHeight_;
    if (legend_.empty()) {
      legendHeight_ = 0;
    }
  }
  drawWidth_ = width_ * (1 - 2 * marginPadding_);
  drawHeight_ = height_ * (1 - 2 * marginPadding_);
  partitionForLegend();

  if (xRange_ > 1e-4 || yRange_ > 1e-4) {
    newScale =
        std::min(double(drawWidth_) / xRange_, double(molHeight_) / yRange_);
    double fix_scale = newScale;
    // after all that, use the fixed scale unless it's too big, in which case
    // scale the drawing down to fit.
    // fixedScale takes precedence if both it and fixedBondLength are given.
    if (drawOptions_.fixedBondLength > 0.0) {
      fix_scale = drawOptions_.fixedBondLength;
    }
    if (drawOptions_.fixedScale > 0.0) {
      fix_scale = double(drawWidth_) * drawOptions_.fixedScale;
    }
    if (newScale > fix_scale) {
      newScale = fix_scale;
    }
  }
  double scale_mult = newScale / scale_;
  scale_ *= scale_mult;
  if (drawOptions_.fixedFontSize != -1) {
    fontScale_ = drawOptions_.fixedFontSize / textDrawer_.baseFontSize();
  } else {
    fontScale_ *= scale_mult;
  }
}

// ****************************************************************************
void DrawMol::findExtremes() {
  for (const auto &ps : preShapes_) {
    ps->findExtremes(xMin_, xMax_, yMin_, yMax_);
  }
  for (const auto &bond : bonds_) {
    bond->findExtremes(xMin_, xMax_, yMin_, yMax_);
  }
  for (const auto &atLab : atomLabels_) {
    if (atLab) {
      atLab->findExtremes(xMin_, xMax_, yMin_, yMax_);
    }
  }
  for (const auto &hl : highlights_) {
    hl->findExtremes(xMin_, xMax_, yMin_, yMax_);
  }
  if (includeAnnotations_) {
    for (const auto &a : annotations_) {
      a->findExtremes(xMin_, xMax_, yMin_, yMax_);
    }
  }
  findRadicalExtremes(radicals_, xMin_, xMax_, yMin_, yMax_);
  for (const auto &ps : postShapes_) {
    ps->findExtremes(xMin_, xMax_, yMin_, yMax_);
  }

  if (atCds_.empty()) {
    xMin_ = yMin_ = -1.0;
    xMax_ = yMax_ = 1.0;
  }
  // Calculate the x and y spans.  Don't include the padding, as that's
  // now taken into account with drawWidth_ and drawHeight_.
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
}

// ****************************************************************************
void DrawMol::changeToDrawCoords() {
  Point2D trans, scale, toCentre;
  getDrawTransformers(trans, scale, toCentre);
  transformAll(&trans, &scale, &toCentre);
}

// ****************************************************************************
void DrawMol::draw(MolDraw2D &drawer) const {
  PRECONDITION(drawingInitialised_,
               "you must call createDrawingObjects before calling draw")
  if (atCds_.empty()) {
    return;
  }
  auto keepScale = drawer.scale();
  drawer.setScale(scale_);
  auto keepFontScale = textDrawer_.fontScale();
  textDrawer_.setFontScale(fontScale_, true);

  for (auto &ps : preShapes_) {
    ps->draw(drawer);
  }
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
  if (includeAnnotations_) {
    for (auto &annot : annotations_) {
      annot->draw(drawer);
    }
  }
  if (drawOptions_.includeRadicals) {
    drawRadicals(drawer);
  }
  for (auto &ps : postShapes_) {
    ps->draw(drawer);
  }
  for (auto &leg : legends_) {
    leg->draw(drawer);
  }
  drawer.setScale(keepScale);
  textDrawer_.setFontScale(keepFontScale, true);
}

// ****************************************************************************
void DrawMol::drawRadicals(MolDraw2D &drawer) const {
  // take account of differing font scale and main scale if we've hit
  // max or min font size.
  double spot_rad = 0.2 * drawOptions_.multipleBondOffset * fontScale_;
  drawer.setColour(DrawColour(0.0, 0.0, 0.0));
  auto draw_spot = [&](const Point2D &cds) {
    bool ofp = drawer.fillPolys();
    drawer.setFillPolys(true);
    double olw = drawer.lineWidth();
    drawer.setLineWidth(0);
    drawer.drawArc(cds, spot_rad, 0, 360, true);
    drawer.setLineWidth(olw);
    drawer.setFillPolys(ofp);
  };

  // cds in draw coords
  auto draw_spots = [&](const Point2D &cds, int num_spots, double width,
                        int dir = 0) {
    Point2D ncds = cds;
    switch (num_spots) {
      case 3:
        if (dir) {
          ncds.y = cds.y - 0.6 * width + spot_rad;
        } else {
          ncds.x = cds.x - 0.6 * width + spot_rad;
        }
        draw_spot(ncds);
        if (dir) {
          ncds.y = cds.y + 0.6 * width - spot_rad;
        } else {
          ncds.x = cds.x + 0.6 * width - spot_rad;
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
  for (const auto &atom : drawMol_->atoms()) {
    int num_rade = atom->getNumRadicalElectrons();
    if (!num_rade) {
      continue;
    }
    auto rad_rect = std::get<0>(radicals_[rad_num]);
    OrientType draw_or = std::get<1>(radicals_[rad_num]);
    int atIdx = std::get<2>(radicals_[rad_num]);
    drawer.setActiveAtmIdx(atIdx);
    if (draw_or == OrientType::N || draw_or == OrientType::S ||
        draw_or == OrientType::C) {
      draw_spots(rad_rect.trans_, num_rade, rad_rect.width_, 0);
    } else {
      draw_spots(rad_rect.trans_, num_rade, rad_rect.height_, 1);
    }
    drawer.setActiveAtmIdx();
    ++rad_num;
  }
}

// ****************************************************************************
void DrawMol::resetEverything() {
  scale_ = 1.0;
  fontScale_ = 1.0;
  textDrawer_.setFontScale(1.0, true);
  xMin_ = std::numeric_limits<double>::max() / 2.0;
  yMin_ = std::numeric_limits<double>::max() / 2.0;
  xMax_ = std::numeric_limits<double>::lowest() / 2.0;
  yMax_ = std::numeric_limits<double>::lowest() / 2.0;
  xRange_ = std::numeric_limits<double>::max();
  yRange_ = std::numeric_limits<double>::max();
  meanBondLength_ = 0.0;
  atCds_.clear();
  bonds_.clear();
  preShapes_.clear();
  postShapes_.clear();
  atomicNums_.clear();
  atomSyms_.clear();
  atomLabels_.clear();
  highlights_.clear();
  annotations_.clear();
  legends_.clear();
  radicals_.clear();
  singleBondLines_.clear();
}

// ****************************************************************************
void DrawMol::shrinkToFit(bool withPadding) {
  double padding = withPadding ? marginPadding_ : 0;
  int newWidth = std::ceil((2 * padding + 1) * xRange_ * scale_);
  int newHeight = std::ceil((2 * padding + 1) * yRange_ * scale_);
  Point2D corr((newWidth - width_) / 2.0, (newHeight - height_) / 2.0);
  transformAll(&corr, nullptr, nullptr);
  width_ = newWidth;
  drawWidth_ = width_ * (1 - 2 * padding);
  height_ = newHeight;
  if (!legend_.empty()) {
    partitionForLegend();
    legends_.clear();
    extractLegend();
  } else {
    legendHeight_ = 0;
    molHeight_ = height_;
    drawHeight_ = height_ * (1 - 2 * padding);
  }
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
  PRECONDITION(atom.getQuery()->getNegation() ||
                   atom.getQuery()->getDescription() == "AtomOr",
               "bad query type");

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
const std::map<std::string, std::string> &getComplexQuerySymbolMap() {
  static const std::map<std::string, std::string> complexQuerySymbolMap{
      {"![H]", "A"},
      {"![C,H]", "Q"},
      {"![C]", "QH"},
      {"[F,Cl,Br,I,At]", "X"},
      {"[F,Cl,Br,I,At,H]", "XH"},
      {"![He,B,C,N,O,F,Ne,Si,P,S,Cl,Ar,As,Se,Br,Kr,Te,I,Xe,At,Rn,H]", "M"},
      {"![He,B,C,N,O,F,Ne,Si,P,S,Cl,Ar,As,Se,Br,Kr,Te,I,Xe,At,Rn]", "MH"},
  };
  return complexQuerySymbolMap;
}

std::set<std::string> createComplexQuerySymbolSet() {
  std::set<std::string> complexQuerySymbolSet;
  const auto &querySymbolMap = getComplexQuerySymbolMap();
  std::transform(
      querySymbolMap.begin(), querySymbolMap.end(),
      std::inserter(complexQuerySymbolSet, complexQuerySymbolSet.begin()),
      [](const auto &pair) { return pair.second; });
  return complexQuerySymbolSet;
}

const std::set<std::string> &getComplexQuerySymbolSet() {
  static const auto complexQuerySymbolSet = createComplexQuerySymbolSet();
  return complexQuerySymbolSet;
}

std::string getComplexQueryAtomEquivalent(const std::string &query) {
  const auto &complexQuerySymbolMap = getComplexQuerySymbolMap();
  auto it = complexQuerySymbolMap.find(query);
  return (it == complexQuerySymbolMap.end() ? query : it->second);
}

bool hasSymbolQueryType(const Atom &atom) {
  return getComplexQuerySymbolSet().count(atom.getQueryType()) > 0;
}

// ****************************************************************************
std::string DrawMol::getAtomSymbol(const Atom &atom,
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
  } else if (drawOptions_.useComplexQueryAtomSymbols &&
             hasSymbolQueryType(atom)) {
    symbol = atom.getQueryType();
  } else if (isAtomListQuery(&atom)) {
    symbol = getAtomListText(atom);
    if (drawOptions_.useComplexQueryAtomSymbols) {
      symbol = getComplexQueryAtomEquivalent(symbol);
    }
  } else if (isComplexQuery(&atom)) {
    symbol = "?";
    std::string mapNum;
    if (atom.getPropIfPresent("molAtomMapNumber", mapNum)) {
      symbol += ":" + mapNum;
    }
  } else if (drawOptions_.atomLabelDeuteriumTritium &&
             atom.getAtomicNum() == 1 && (iso == 2 || iso == 3)) {
    symbol = ((iso == 2) ? "D" : "T");
    iso = 0;
  } else {
    literal_symbol = false;
    std::vector<std::string> preText, postText;

    // first thing after the symbol is the atom map
    std::string mapNum;
    if (atom.getPropIfPresent("molAtomMapNumber", mapNum)) {
      postText.push_back(std::string(":") + mapNum);
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
        h += std::string("<sub>") + std::to_string(num_h) +
             std::string("</sub>");
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
  return symbol;
}

// ****************************************************************************
OrientType DrawMol::getAtomOrientation(const RDKit::Atom &atom) const {
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
  for (const auto bond : mol.atomBonds(&atom)) {
    const Point2D &at2_cds = atCds_[bond->getOtherAtomIdx(atom.getIdx())];
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
        for (const auto bond : mol.atomBonds(&atom)) {
          const Point2D &at2_cds = atCds_[bond->getOtherAtomIdx(atom.getIdx())];
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
void DrawMol::calcMeanBondLength() {
  // meanBondLength_ initialised to 0.0 in class declaration
  if (meanBondLength_ == 0.0) {
    meanBondLength_ = MolDraw2DUtils::meanBondLength(*drawMol_);
  }
}

// ****************************************************************************
void DrawMol::partitionForLegend() {
  if (legend_.empty()) {
    molHeight_ = drawHeight_;
    legendHeight_ = 0;
  } else {
    if (!flexiCanvasY_) {
      legendHeight_ = int(drawOptions_.legendFraction * float(drawHeight_));
      molHeight_ = drawHeight_ - legendHeight_;
    } else {
      molHeight_ = drawHeight_;
      // the legendHeight_ isn't needed for the flexiCanvas
    }
  }
}

// ****************************************************************************
// This must be called after calculateScale() because it needs the final
// font size to work out the legend font size which is given in
// drawOptions().legendFontSize in pixels, and then scaled down to fit
// the width_ and legendHeight_ if necessary.
void DrawMol::extractLegend() {
  if (legend_.empty()) {
    return;
  }
  auto calc_legend_height = [&](const std::vector<std::string> &legend_bits,
                                double relFontScale, double &total_width,
                                double &total_height) {
    total_width = total_height = 0;
    for (auto &bit : legend_bits) {
      double height, width;
      DrawAnnotation da(bit, TextAlignType::MIDDLE, "legend", relFontScale,
                        Point2D(0.0, 0.0), drawOptions_.legendColour,
                        textDrawer_);
      da.getDimensions(width, height);
      total_height += height;
      total_width = std::max(total_width, width);
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

  // work out a font scale that allows the pieces to fit, remembering there's
  // padding round the picture.
  double fsize = textDrawer_.fontSize();
  double relFontScale = drawOptions_.legendFontSize / fsize;
  double total_width, total_height;
  calc_legend_height(legend_bits, relFontScale, total_width, total_height);
  if (total_width >= drawWidth_) {
    if (!flexiCanvasX_) {
      relFontScale *= double(drawWidth_) / total_width;
      calc_legend_height(legend_bits, relFontScale, total_width, total_height);
    } else {
      width_ = total_width * (1 + 2 * marginPadding_);
      drawWidth_ = total_width;
    }
  }

  if (!flexiCanvasY_) {
    auto adjLegHt = drawHeight_ * drawOptions_.legendFraction;
    // subtract off space for the padding.
    if (total_height > adjLegHt) {
      relFontScale *= double(adjLegHt) / total_height;
      calc_legend_height(legend_bits, relFontScale, total_width, total_height);
    }
  } else {
    // a small gap between the legend and the picture looks better,
    // and make it at least 2 pixels.
    double extra_padding = total_height * marginPadding_;
    extra_padding = extra_padding < 2.0 ? 2.0 : extra_padding;
    legendHeight_ = total_height + extra_padding;
    drawHeight_ += legendHeight_;
    height_ += legendHeight_;
  }

  Point2D loc(drawWidth_ / 2 + xOffset_ + width_ * marginPadding_,
              marginPadding_ * height_ + drawHeight_ + yOffset_);
  for (auto bit : legend_bits) {
    DrawAnnotation *da =
        new DrawAnnotation(bit, TextAlignType::MIDDLE, "legend", relFontScale,
                           loc, drawOptions_.legendColour, textDrawer_);
    legends_.emplace_back(da);
  }
  // The letters have different amounts above and below the centre,
  // which matters when placing them vertically.
  // Draw them from the bottom up.
  double xmin, xmax, ymin, ymax;
  xmin = ymin = std::numeric_limits<double>::max();
  xmax = ymax = std::numeric_limits<double>::lowest();
  legends_.back()->findExtremes(xmin, xmax, ymin, ymax);
  double lastBelow = legends_.back()->pos_.y - ymax;
  double lastAbove = legends_.back()->pos_.y - ymin;
  legends_.back()->pos_.y += lastBelow;
  for (int i = legends_.size() - 2; i >= 0; --i) {
    xmin = ymin = std::numeric_limits<double>::max();
    xmax = ymax = std::numeric_limits<double>::lowest();
    legends_[i]->findExtremes(xmin, xmax, ymin, ymax);
    double thisBelow = legends_[i]->pos_.y - ymax;
    double thisAbove = legends_[i]->pos_.y - ymin;
    legends_[i]->pos_.y = legends_[i + 1]->pos_.y + thisBelow - lastAbove;
    lastAbove = thisAbove;
  }
}

// ****************************************************************************
void DrawMol::makeStandardBond(Bond *bond, double doubleBondOffset) {
  int begAt = bond->getBeginAtomIdx();
  int endAt = bond->getEndAtomIdx();
  // If the 2 atoms are on top of each other, don't do anything.  We can
  // end up with NaN for points in the shapes for things like chiral atoms
  // (issue 6569).
  const Point2D &at1_cds = atCds_[begAt];
  const Point2D &at2_cds = atCds_[endAt];
  if ((at1_cds - at2_cds).lengthSq() < 0.0001) {
    return;
  }
  std::pair<DrawColour, DrawColour> cols = getBondColours(bond);

  auto bt = bond->getBondType();
  if (bt == Bond::DOUBLE || bt == Bond::AROMATIC) {
    makeDoubleBondLines(bond, doubleBondOffset, cols);
  } else if (bt == Bond::SINGLE && (bond->getBondDir() == Bond::BEGINWEDGE ||
                                    bond->getBondDir() == Bond::BEGINDASH)) {
    makeWedgedBond(bond, cols);
  } else if (bt == Bond::SINGLE && bond->getBondDir() == Bond::UNKNOWN) {
    makeWavyBond(bond, doubleBondOffset, cols);
  } else if (bt == Bond::DATIVE || bt == Bond::DATIVEL || bt == Bond::DATIVER) {
    makeDativeBond(bond, doubleBondOffset, cols);
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
  const Point2D &at1_cds = atCds_[begAt->getIdx()];
  const Point2D &at2_cds = atCds_[endAt->getIdx()];
  // If the 2 atoms are on top of each other, don't do anything.  We can
  // end up with NaN for points in the shapes for things like chiral atoms
  // (issue 6569).
  if ((at1_cds - at2_cds).lengthSq() < 0.0001) {
    return;
  }

  Point2D end1, end2;
  adjustBondEndsForLabels(begAt->getIdx(), endAt->getIdx(), end1, end2);
  Point2D sat1 = atCds_[begAt->getIdx()];
  Point2D sat2 = atCds_[endAt->getIdx()];
  atCds_[begAt->getIdx()] = end1;
  atCds_[endAt->getIdx()] = end2;

  auto midp = (at2_cds + at1_cds) / 2.;
  auto tdash = shortDashes;
  const DrawColour &queryColour = drawOptions_.queryColour;

  bool drawGenericQuery = false;
  int at1Idx = begAt->getIdx();
  int at2Idx = endAt->getIdx();
  if (qry->getDescription() == "SingleOrDoubleBond") {
    at1Idx = begAt->getIdx();
    at2Idx = drawOptions_.splitBonds ? -1 : endAt->getIdx();
    newBondLine(at1_cds, midp, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), noDash);
    Point2D l1s, l1f, l2s, l2f;
    calcDoubleBondLines(doubleBondOffset, *bond, l1s, l1f, l2s, l2f);
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
    newBondLine(at1_cds, midp, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), noDash);
    Point2D l1s, l1f, l2s, l2f;
    calcDoubleBondLines(doubleBondOffset, *bond, l1s, l1f, l2s, l2f);
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
    calcDoubleBondLines(doubleBondOffset, *bond, l1s, l1f, l2s, l2f);
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
    newBondLine(at1_cds, at2_cds, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), tdash);
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
        DrawShapePolyLine *pl =
            new DrawShapePolyLine(pts, 1, false, queryColour, false,
                                  begAt->getIdx() + activeAtmIdxOffset_,
                                  endAt->getIdx() + activeAtmIdxOffset_,
                                  bond->getIdx() + activeBndIdxOffset_, noDash);
        bonds_.emplace_back(pl);
      } else {
        segment /= segment.length() * 10;
        auto l = segment.length();
        Point2D p1 = midp + segment;
        Point2D p2 = Point2D(l, l);
        std::vector<Point2D> pts{p1, p2};
        DrawShapeEllipse *ell =
            new DrawShapeEllipse(pts, 1, false, queryColour, false);
        bonds_.emplace_back(ell);
        p1 = midp - segment;
        p2 = Point2D(l, l);
        pts = std::vector<Point2D>{p1, p2};
        ell = new DrawShapeEllipse(pts, 1, false, queryColour, false);
        bonds_.emplace_back(ell);
      }
    } else {
      drawGenericQuery = true;
    }
  } else {
    drawGenericQuery = true;
  }
  if (drawGenericQuery) {
    newBondLine(at1_cds, at2_cds, queryColour, queryColour, at1Idx, at2Idx,
                bond->getIdx(), dots);
    bonds_.back()->lineWidth_ = drawOptions_.bondLineWidth;
    bonds_.back()->scaleLineWidth_ = false;
  }
  atCds_[begAt->getIdx()] = sat1;
  atCds_[endAt->getIdx()] = sat2;
}

// ****************************************************************************
void DrawMol::makeDoubleBondLines(
    Bond *bond, double doubleBondOffset,
    const std::pair<DrawColour, DrawColour> &cols) {
  Point2D end1, end2;
  int at1Idx = bond->getBeginAtomIdx();
  int at2Idx = bond->getEndAtomIdx();
  adjustBondEndsForLabels(at1Idx, at2Idx, end1, end2);

  bool skipOuterLine = false;
  if (bond->getBondDir() == Bond::BEGINWEDGE ||
      bond->getBondDir() == Bond::BEGINDASH) {
    makeWedgedBond(bond, cols);
    skipOuterLine = true;
  }

  Point2D l1s, l1f, l2s, l2f, sat1, sat2;
  sat1 = atCds_[at1Idx];
  atCds_[at1Idx] = end1;
  sat2 = atCds_[at2Idx];
  atCds_[at2Idx] = end2;
  calcDoubleBondLines(doubleBondOffset, *bond, l1s, l1f, l2s, l2f);
  int bondIdx = bond->getIdx();
  if (!skipOuterLine) {
    newBondLine(l1s, l1f, cols.first, cols.second, at1Idx, at2Idx, bondIdx,
                noDash);
  }
  if (bond->getBondType() == Bond::AROMATIC) {
    newBondLine(l2s, l2f, cols.first, cols.second, at1Idx, at2Idx, bondIdx,
                dashes);
  } else {
    // if it's a two colour line, then a simple line will be exactly half
    // one colour and half the other.  The second line to a terminal atom
    // in, for example, an aldehyde, such as in catch_tests.cpp's
    // testGithub_5269_2.svg, might be asymmetrically shorter, so we don't
    // want to change colour at halfway
    auto l1 = (l1s - l1f).lengthSq();
    auto l2 = (l2s - l2f).lengthSq();
    if ((bond->getBeginAtom()->getDegree() == 1 ||
         bond->getEndAtom()->getDegree() == 1) &&
        cols.first != cols.second && fabs(l1 - l2) > 0.01) {
      double midlen = sqrt(l1) / 2.0;
      Point2D notMid;
      if (bond->getBeginAtom()->getDegree() == 1) {
        Point2D lineDir = l2s.directionVector(l2f);
        notMid = l2s + lineDir * midlen;
      } else {
        Point2D lineDir = l2f.directionVector(l2s);
        notMid = l2f + lineDir * midlen;
      }
      newBondLine(l2s, notMid, cols.first, cols.first, at1Idx, at2Idx, bondIdx,
                  noDash);
      newBondLine(notMid, l2f, cols.second, cols.second, at1Idx, at2Idx,
                  bondIdx, noDash);
    } else {
      newBondLine(l2s, l2f, cols.first, cols.second, at1Idx, at2Idx, bondIdx,
                  noDash);
    }
  }
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
  calcTripleBondLines(doubleBondOffset, *bond, l1s, l1f, l2s, l2f);
  newBondLine(l1s, l1f, cols.first, cols.second, at1Idx, at2Idx, bondIdx,
              noDash);
  newBondLine(l2s, l2f, cols.first, cols.second, at1Idx, at2Idx, bondIdx,
              noDash);
  atCds_[at1Idx] = sat1;
  atCds_[at2Idx] = sat2;
}

// ****************************************************************************
void DrawMol::makeWedgedBond(Bond *bond,
                             const std::pair<DrawColour, DrawColour> &cols) {
  auto at1 = bond->getBeginAtom();
  auto at2 = bond->getEndAtom();
  auto col1 = cols.first;
  auto col2 = cols.second;
  if (drawOptions_.singleColourWedgeBonds) {
    col1 = drawOptions_.symbolColour;
    col2 = drawOptions_.symbolColour;
  }
  // If either of the atoms has a label, make the padding a bit bigger
  // so the end of the wedge doesn't run up to the atom symbol.
  // Obviously, if we ever change how the padding round the label is
  // calculated, currently a function of the mean bond length, this won't work.
  if (atomLabels_[at1->getIdx()] || atomLabels_[at2->getIdx()]) {
    meanBondLength_ *= 2.0;
  }
  Point2D end1, end2;
  adjustBondEndsForLabels(at1->getIdx(), at2->getIdx(), end1, end2);
  if (atomLabels_[at1->getIdx()] || atomLabels_[at2->getIdx()]) {
    meanBondLength_ /= 2.0;
  }
  const Point2D &at1_cds = atCds_[at1->getIdx()];
  const Point2D &at2_cds = atCds_[at2->getIdx()];

  Point2D perp = calcPerpendicular(at1_cds, at2_cds);
  // Set the 'fatness' of the wedge to be a fraction of the mean bond
  // length, so we should always see something.
  Point2D disp = perp * drawOptions_.multipleBondOffset * meanBondLength_ / 2.0;
  Point2D t1 = end2 + disp;
  Point2D t2 = end2 - disp;
  std::vector<Point2D> pts{end1, t1, t2};

  // deliberately not scaling highlighted bond width
  DrawShape *s;
  double lineWidth = drawOptions_.bondLineWidth < 1.0
                         ? drawOptions_.bondLineWidth
                         : drawOptions_.bondLineWidth / 2.0;
  if (Bond::BEGINWEDGE == bond->getBondDir()) {
    std::vector<Point2D> otherBondVecs;
    findOtherBondVecs(at2, at1, otherBondVecs);
    s = new DrawShapeSolidWedge(pts, col1, col2, drawOptions_.splitBonds,
                                otherBondVecs, lineWidth,
                                at1->getIdx() + activeAtmIdxOffset_,
                                at2->getIdx() + activeAtmIdxOffset_,
                                bond->getIdx() + activeBndIdxOffset_);
  } else {
    bool oneLessDash(at2->getDegree() > 1);
    s = new DrawShapeDashedWedge(pts, col1, col2, oneLessDash, lineWidth,
                                 at1->getIdx() + activeAtmIdxOffset_,
                                 at2->getIdx() + activeAtmIdxOffset_,
                                 bond->getIdx() + activeBndIdxOffset_);
  }
  bonds_.push_back(std::unique_ptr<DrawShape>(s));
}

// ****************************************************************************
void DrawMol::makeWavyBond(Bond *bond, double offset,
                           const std::pair<DrawColour, DrawColour> &cols) {
  auto at1 = bond->getBeginAtom();
  auto at2 = bond->getEndAtom();
  Point2D end1, end2;
  adjustBondEndsForLabels(at1->getIdx(), at2->getIdx(), end1, end2);
  std::vector<Point2D> pts{end1, end2};
  DrawShapeWavyLine *s = new DrawShapeWavyLine(
      pts, drawOptions_.bondLineWidth, false, cols.first, cols.second, offset,
      at1->getIdx() + activeAtmIdxOffset_, at2->getIdx() + activeAtmIdxOffset_,
      bond->getIdx() + activeBndIdxOffset_);
  bonds_.push_back(std::unique_ptr<DrawShape>(s));
}

// ****************************************************************************
void DrawMol::makeDativeBond(Bond *bond, double offset,
                             const std::pair<DrawColour, DrawColour> &cols) {
  auto at1 = bond->getBeginAtom();
  auto at2 = bond->getEndAtom();
  Point2D end1, end2;
  adjustBondEndsForLabels(at1->getIdx(), at2->getIdx(), end1, end2);

  Point2D mid = (end1 + end2) * 0.5;
  int atid2 = drawOptions_.splitBonds ? at1->getIdx() : at2->getIdx();
  newBondLine(end1, mid, cols.first, cols.first, at1->getIdx(), atid2,
              bond->getIdx(), noDash);
  std::vector<Point2D> pts{mid, end2};
  // Adjust the fraction of the line length that will be arrowhead so that
  // it is a consistent number of pixels.
  auto frac = 2.0 * offset / (end2 - end1).length();
  DrawShapeArrow *a = new DrawShapeArrow(
      pts, drawOptions_.bondLineWidth, false, cols.second, true,
      at1->getIdx() + activeAtmIdxOffset_, atid2 + activeAtmIdxOffset_,
      bond->getIdx() + activeBndIdxOffset_, frac, M_PI / 12);
  bonds_.push_back(std::unique_ptr<DrawShape>(a));
}

// ****************************************************************************
void DrawMol::makeZeroBond(Bond *bond,
                           const std::pair<DrawColour, DrawColour> &cols,
                           const DashPattern &dashPattern) {
  auto at1 = bond->getBeginAtom();
  auto at2 = bond->getEndAtom();
  Point2D end1, end2;
  adjustBondEndsForLabels(at1->getIdx(), at2->getIdx(), end1, end2);
  newBondLine(end1, end2, cols.first, cols.second, at1->getIdx(), at2->getIdx(),
              bond->getIdx(), dashPattern);
}

// ****************************************************************************
void DrawMol::adjustBondEndsForLabels(int begAtIdx, int endAtIdx,
                                      Point2D &begCds, Point2D &endCds) const {
  // The scale factor is empirical.
  double padding = 0.033 * meanBondLength_;
  if (drawOptions_.additionalAtomLabelPadding > 0.0) {
    padding += drawOptions_.additionalAtomLabelPadding;
  }
  begCds = atCds_[begAtIdx];
  endCds = atCds_[endAtIdx];
  if (atomLabels_[begAtIdx]) {
    adjustBondEndForString(endCds, padding, atomLabels_[begAtIdx]->rects_,
                           begCds);
  }
  if (atomLabels_[endAtIdx]) {
    adjustBondEndForString(begCds, padding, atomLabels_[endAtIdx]->rects_,
                           endCds);
  }
}

// ****************************************************************************
void DrawMol::newBondLine(const Point2D &pt1, const Point2D &pt2,
                          const DrawColour &col1, const DrawColour &col2,
                          int atom1Idx, int atom2Idx, int bondIdx,
                          const DashPattern &dashPattern) {
  bool scaleWidth = drawOptions_.scaleBondWidth;
  double lineWidth = drawOptions_.bondLineWidth;
  if (!drawOptions_.continuousHighlight &&
      std::find(highlightBonds_.begin(), highlightBonds_.end(), bondIdx) !=
          highlightBonds_.end()) {
    scaleWidth = drawOptions_.scaleHighlightBondWidth;
    lineWidth = getHighlightBondWidth(drawOptions_, bondIdx, nullptr) / 4;
  }

  if (col1 == col2 && !drawOptions_.splitBonds) {
    std::vector<Point2D> pts{pt1, pt2};
    DrawShape *b = new DrawShapeSimpleLine(
        pts, lineWidth, scaleWidth, col1, atom1Idx + activeAtmIdxOffset_,
        atom2Idx + activeAtmIdxOffset_, bondIdx + activeBndIdxOffset_,
        dashPattern);
    bonds_.push_back(std::unique_ptr<DrawShape>(b));
    if (dashPattern == noDash) {
      singleBondLines_.push_back(bonds_.size() - 1);
    }
  } else {
    Point2D mid = (pt1 + pt2) / 2.0;
    std::vector<Point2D> pts1{pt1, mid};
    int at1Idx = atom1Idx;
    int at2Idx = drawOptions_.splitBonds ? -1 : atom2Idx;
    DrawShape *b1 = new DrawShapeSimpleLine(
        pts1, lineWidth, scaleWidth, col1, at1Idx + activeAtmIdxOffset_,
        at2Idx + activeAtmIdxOffset_, bondIdx + activeBndIdxOffset_,
        dashPattern);
    bonds_.push_back(std::unique_ptr<DrawShape>(b1));
    if (dashPattern == noDash) {
      singleBondLines_.push_back(bonds_.size() - 1);
    }
    at1Idx = drawOptions_.splitBonds ? atom2Idx : atom1Idx;
    std::vector<Point2D> pts2{mid, pt2};
    DrawShape *b2 = new DrawShapeSimpleLine(
        pts2, lineWidth, scaleWidth, col2, at1Idx + activeAtmIdxOffset_,
        at2Idx + activeAtmIdxOffset_, bondIdx + activeBndIdxOffset_,
        dashPattern);
    bonds_.push_back(std::unique_ptr<DrawShape>(b2));
    if (dashPattern == noDash) {
      singleBondLines_.push_back(bonds_.size() - 1);
    }
  }
}

// ****************************************************************************
std::pair<DrawColour, DrawColour> DrawMol::getBondColours(Bond *bond) {
  DrawColour col1, col2;

  bool highlight_bond = false;
  if (std::find(highlightBonds_.begin(), highlightBonds_.end(),
                bond->getIdx()) != highlightBonds_.end()) {
    highlight_bond = true;
  }

  if (!bondColours_.empty()) {
    col1 = bondColours_[bond->getIdx()].first;
    col2 = bondColours_[bond->getIdx()].second;
  } else {
    if (!highlight_bond || drawOptions_.continuousHighlight) {
      col1 = getColour(bond->getBeginAtomIdx());
      col2 = getColour(bond->getEndAtomIdx());
    } else {
      if (highlightBondMap_.find(bond->getIdx()) != highlightBondMap_.end()) {
        col1 = col2 = highlightBondMap_.find(bond->getIdx())->second;
      } else {
        col1 = col2 = drawOptions_.highlightColour;
      }
    }
  }

  return std::make_pair(col1, col2);
}

// ****************************************************************************
void DrawMol::makeContinuousHighlights(double scale) {
  double tgt_lw = getHighlightBondWidth(drawOptions_, -1, nullptr);
  if (tgt_lw < 2.0) {
    tgt_lw = 2.0;
  }
  if (!drawOptions_.continuousHighlight) {
    tgt_lw /= 4.0;
  }

  if (!highlightBonds_.empty()) {
    makeBondHighlightLines(tgt_lw, scale);
  }
  if (!highlightAtoms_.empty()) {
    makeAtomEllipseHighlights(tgt_lw);
  }
}

// ****************************************************************************
void DrawMol::makeAtomCircleHighlights() {
  DrawColour col;
  for (const auto &at : drawMol_->atoms()) {
    unsigned int thisIdx = at->getIdx();
    if (std::find(highlightAtoms_.begin(), highlightAtoms_.end(), thisIdx) !=
        highlightAtoms_.end()) {
      if (highlightAtomMap_.find(thisIdx) != highlightAtomMap_.end()) {
        col = highlightAtomMap_.find(thisIdx)->second;
      } else {
        col = drawOptions_.highlightColour;
      }
      double radius = drawOptions_.highlightRadius;
      if (highlightRadii_.find(thisIdx) != highlightRadii_.end()) {
        radius = highlightRadii_.find(thisIdx)->second;
      }
      Point2D radii(radius, radius);
      std::vector<Point2D> pts{atCds_[thisIdx], radii};
      DrawShape *ell = new DrawShapeEllipse(pts, 2, false, col, true,
                                            thisIdx + activeAtmIdxOffset_);
      highlights_.push_back(std::unique_ptr<DrawShape>(ell));
    }
  }
}

// ****************************************************************************
void DrawMol::makeAtomEllipseHighlights(double lineWidth) {
  if (!drawOptions_.fillHighlights) {
    // we need a narrower circle
    lineWidth /= 2.0;
  }
  for (const auto &atom : drawMol_->atoms()) {
    unsigned int thisIdx = atom->getIdx();
    if (std::find(highlightAtoms_.begin(), highlightAtoms_.end(), thisIdx) !=
        highlightAtoms_.end()) {
      DrawColour col = drawOptions_.highlightColour;
      if (highlightAtomMap_.find(thisIdx) != highlightAtomMap_.end()) {
        col = highlightAtomMap_.find(thisIdx)->second;
      }
      Point2D centre = atCds_[thisIdx];
      double xradius, yradius;
      if (highlightRadii_.find(thisIdx) != highlightRadii_.end()) {
        xradius = highlightRadii_.find(thisIdx)->second;
      } else {
        xradius = drawOptions_.highlightRadius;
      }
      yradius = xradius;
      if (!drawOptions_.atomHighlightsAreCircles && atomLabels_[thisIdx]) {
        double xMin, yMin, xMax, yMax;
        xMin = yMin = std::numeric_limits<double>::max();
        xMax = yMax = std::numeric_limits<double>::lowest();
        atomLabels_[thisIdx]->findExtremes(xMin, xMax, yMin, yMax);
        static const double root_2 = sqrt(2.0);
        xradius = std::max(xradius, root_2 * 0.5 * (xMax - xMin));
        yradius = std::max(yradius, root_2 * 0.5 * (yMax - yMin));
        centre.x = 0.5 * (xMax + xMin);
        centre.y = 0.5 * (yMax + yMin);
      }
      Point2D radii(xradius, yradius);
      std::vector<Point2D> pts{centre, radii};
      DrawShape *ell = new DrawShapeEllipse(pts, lineWidth, true, col, true,
                                            thisIdx + activeAtmIdxOffset_);
      highlights_.push_back(std::unique_ptr<DrawShape>(ell));
    }
  }
}

// ****************************************************************************
void DrawMol::makeBondHighlightLines(double lineWidth, double scale) {
  // find the neighbours of atom that aren't otherAtom and that are
  // bonded to atom with a highlighted bond
  auto findHighBondNbrs = [&](const Atom *atom, const Atom *otherAtom,
                              std::vector<Atom *> &highNbrs) -> void {
    for (const auto bond : drawMol_->atomBonds(atom)) {
      auto nbr = bond->getOtherAtom(atom);
      if (nbr == otherAtom) {
        continue;
      }
      if (std::find(highlightBonds_.begin(), highlightBonds_.end(),
                    bond->getIdx()) != highlightBonds_.end()) {
        highNbrs.push_back(nbr);
      }
    }
  };

  // this is essentially the inverse of MolDraw2D::getDrawLineWidth
  // which ignores drawOptions_.scaleHighlightBondWidth and only
  // uses drawOptions_.scaleBondWidth.
  if (!drawOptions_.scaleHighlightBondWidth) {
    // so that when we scale it up again, it comes out the right size
    lineWidth /= scale;
  } else {
    // same conversion factor as in MolDraw2D::getDrawLineWidth()
    lineWidth *= lineWidthScaleFactor;
  }
  for (const auto atom : drawMol_->atoms()) {
    auto thisIdx = atom->getIdx();
    for (const auto bond : drawMol_->atomBonds(atom)) {
      unsigned int nbrIdx = bond->getOtherAtomIdx(thisIdx);
      if (nbrIdx < static_cast<unsigned int>(atCds_.size()) &&
          nbrIdx > thisIdx) {
        if (std::find(highlightBonds_.begin(), highlightBonds_.end(),
                      bond->getIdx()) != highlightBonds_.end()) {
          // This bond is to be highlighted by drawing a 4-6-sided
          // polygon underneath it.  If it is an isolated highlight, it
          // will be a rectangle underneath the bond.  If it joins
          // another highlighted bond, it will be mitred so the two join
          // without gaps.
          // These effects can be seen in bond_highlights_8.svg produced
          // by catch_tests.cpp.
          DrawColour col = getHighlightBondColour(
              bond, drawOptions_, highlightBonds_, highlightBondMap_,
              highlightAtoms_, highlightAtomMap_);
          std::vector<Atom *> thisHighNbrs;
          std::vector<Atom *> nbrHighNbrs;
          auto nbr = drawMol_->getAtomWithIdx(nbrIdx);
          findHighBondNbrs(atom, nbr, thisHighNbrs);
          findHighBondNbrs(nbr, atom, nbrHighNbrs);
          std::vector<Point2D> end1points;
          makeHighlightEnd(atom, nbr, lineWidth, thisHighNbrs, end1points);
          std::vector<Point2D> end2points;
          makeHighlightEnd(nbr, atom, lineWidth, nbrHighNbrs, end2points);
          std::vector<Point2D> points(end1points);
          points.insert(points.end(), end2points.begin(), end2points.end());
          // The end points are sometimes swapped round, such that a
          // butterfly-type shape is produced rather than a rectangle
          // (see Github5592).  Make a convex hull, using a simplified
          // form of Graham's scan algorithm - all the points
          // are in the convex hull so it's easier.  Graham's scan normally
          // has a second step that removes inner points, and this takes
          // care of any problems with floating point errors in the
          // comparisons below.  The shapes here are at most hexagons with
          // sharp angles so such issues have been deemed unlikely to
          // occur in practice.
          // Sort so the lowest y point is first, with lowest x as
          // tie-breaker.
          std::sort(points.begin(), points.end(),
                    [](Point2D &p1, Point2D &p2) -> bool {
                      if (p1.y < p2.y) {
                        return true;
                      } else if (p1.y == p2.y) {
                        return p1.x < p2.x;
                      }
                      return false;
                    });
          // Now sort points 1 -> end so they are all anti-clockwise
          // around points[0] by checking cross products.
          std::sort(points.begin() + 1, points.end(),
                    [&](Point2D &p1, Point2D &p2) -> bool {
                      auto &p0 = points.front();
                      auto val = (p1.y - p0.y) * (p2.x - p1.x) -
                                 (p1.x - p0.x) * (p2.y - p1.y);
                      if (val == 0.0) {
                        return (p0 - p2).lengthSq() < (p0 - p1).lengthSq();
                      } else if (val < 0.0) {
                        return true;
                      } else {
                        return false;
                      }
                    });
          DrawShape *hb = new DrawShapePolyLine(
              points, 0, false, col, true, thisIdx + activeAtmIdxOffset_,
              nbrIdx + activeAtmIdxOffset_,
              bond->getIdx() + activeBndIdxOffset_);
          highlights_.push_back(std::unique_ptr<DrawShape>(hb));
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMol::calcAnnotationPosition(const Atom *atom,
                                     DrawAnnotation &annot) const {
  PRECONDITION(atom, "no atom");
  double start_ang = getNoteStartAngle(atom);
  Point2D const &atCds = atCds_[atom->getIdx()];
  double radStep = 0.25;
  Point2D leastWorstPos = atCds;
  int leastWorstScore = 100;
  for (int j = 1; j < 4; ++j) {
    double note_rad = j * radStep;
    // experience suggests if there's an atom symbol, the close in
    // radius won't work.
    if (j == 1 && atomLabels_[atom->getIdx()]) {
      continue;
    }
    // scan at 30 degree intervals around the atom looking for somewhere
    // clear for the annotation.
    for (int i = 0; i < 12; ++i) {
      double ang = start_ang + i * 30.0 * M_PI / 180.0;
      annot.pos_.x = atCds.x + cos(ang) * note_rad;
      annot.pos_.y = atCds.y + sin(ang) * note_rad;
      int clashScore = doesNoteClash(annot);
      if (!clashScore) {
        return;
      } else {
        if (clashScore < leastWorstScore) {
          leastWorstScore = clashScore;
          leastWorstPos = annot.pos_;
        }
      }
    }
  }
  annot.pos_ = leastWorstPos;
}

// ****************************************************************************
void DrawMol::calcAnnotationPosition(const Bond *bond,
                                     DrawAnnotation &annot) const {
  PRECONDITION(bond, "no bond");
  Point2D const &at1_cds = atCds_[bond->getBeginAtomIdx()];
  Point2D at2_cds = atCds_[bond->getEndAtomIdx()];
  // If the atoms are on top of each other, perp comes out as NaN which
  // has very deleterious effects.  Issue 6569.  Move at2 by a small
  // amount in an arbitrary direction.
  if ((at1_cds - at2_cds).lengthSq() < 0.0001) {
    at2_cds.x += 0.1;
    at2_cds.y += 0.1;
  }
  Point2D perp = calcPerpendicular(at1_cds, at2_cds);
  Point2D bond_vec = at1_cds.directionVector(at2_cds);
  double bond_len = (at1_cds - at2_cds).length();
  std::vector<double> mid_offsets{0.5, 0.33, 0.66, 0.25, 0.75};
  double offset_step = drawOptions_.multipleBondOffset;
  Point2D leastWorstPos = (at1_cds + at2_cds) / 2.0;
  int leastWorstScore = 100;
  for (auto mo : mid_offsets) {
    Point2D mid = at1_cds + bond_vec * bond_len * mo;
    for (int j = 1; j < 6; ++j) {
      if (j == 1 && bond->getBondType() > 1) {
        continue;  // multiple bonds will need a bigger offset.
      }
      double offset = j * offset_step;
      annot.pos_ = mid + perp * offset;
      int clashScore = doesNoteClash(annot);
      if (!clashScore) {
        return;
      }
      if (clashScore < leastWorstScore) {
        leastWorstPos = annot.pos_;
        leastWorstScore = clashScore;
      }
      annot.pos_ = mid - perp * offset;
      clashScore = doesNoteClash(annot);
      if (!clashScore) {
        return;
      }
      if (clashScore < leastWorstScore) {
        leastWorstPos = annot.pos_;
        leastWorstScore = clashScore;
      }
    }
  }
  annot.pos_ = leastWorstPos;
}

// ****************************************************************************
double DrawMol::getNoteStartAngle(const Atom *atom) const {
  if (atom->getDegree() == 0) {
    return M_PI / 2.0;
  }
  const Point2D &at_cds = atCds_[atom->getIdx()];
  std::vector<Point2D> bond_vecs;
  for (auto nbr : make_iterator_range(drawMol_->getAtomNeighbors(atom))) {
    // If the nbr has the same coords as atom, bond_vec comes out as NaN, NaN
    // (issue 6559), so use a short arbitrary vector instead.
    Point2D bond_vec;
    if ((at_cds - atCds_[nbr]).lengthSq() < 0.0001) {
      bond_vec.x = 0.1;
      bond_vec.y = 0.1;
    } else {
      bond_vec = at_cds.directionVector(atCds_[nbr]);
    }
    bond_vec.normalize();
    bond_vecs.push_back(bond_vec);
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
int DrawMol::doesNoteClash(const DrawAnnotation &annot) const {
  // note that this will return a clash if annot is in annotations_.
  // It's intended only to be used when finding where to put the
  // annotation, so annot should only be added to annotations_ once
  // its position has been determined.
  for (auto &rect : annot.rects_) {
    Point2D otrans = rect->trans_;
    rect->trans_ += annot.pos_;
    // if padding is less than this, the letters can fit between the 2 lines
    // of a double bond, which can lead to sub-optimal placements.
    double padding = scale_ * 0.04;
    int clashScore = doesRectClash(*rect, padding);
    rect->trans_ = otrans;
    if (clashScore) {
      return clashScore;
    }
  }
  return 0;
}

// ****************************************************************************
int DrawMol::doesRectClash(const StringRect &rect, double padding) const {
  // No longer checks if it clashes with highlights.  This frequently
  // results in bad pictures and things look ok on top of highlights
  // (issues 5269 and 5195, PR 5272)

  // see if the rectangle clashes with any of the double bonds themselves,
  // as opposed to the draw shapes derived from them.  Github 5185 shows
  // that sometimes atom indices can just fit between the lines of a
  // double bond.
  // Also, no longer check if it clashes with highlights.  This frequently
  // results in bad pictures and things look ok on top of highlights.
  for (auto bond : drawMol_->bonds()) {
    if (bond->getBondType() == Bond::DOUBLE) {
      auto at1 = bond->getBeginAtomIdx();
      auto at2 = bond->getEndAtomIdx();
      if (doesLineIntersect(rect, atCds_[at1], atCds_[at2], 0.0)) {
        return 1;
      }
    }
  }
  for (const auto &bond : bonds_) {
    if (bond->doesRectClash(rect, padding)) {
      return 1;
    }
  }
  for (const auto &al : atomLabels_) {
    if (al && al->doesRectClash(rect, padding)) {
      return 2;
    }
  }
  for (const auto &a : annotations_) {
    if (a->doesRectClash(rect, padding)) {
      return 3;
    }
  }
  return 0;
}

// ****************************************************************************
OrientType DrawMol::calcRadicalRect(const Atom *atom,
                                    StringRect &rad_rect) const {
  int num_rade = atom->getNumRadicalElectrons();
  double spot_rad = 0.2 * drawOptions_.multipleBondOffset * fontScale_;
  Point2D atCds{atCds_[atom->getIdx()]};

  if (scale_ != 1.0) {
    atCds = getDrawCoords(atom->getIdx());
  }
  OrientType orient = atomSyms_[atom->getIdx()].second;
  double rad_size = (4 * num_rade - 2) * spot_rad / fontScale_;
  double x_min, y_min, x_max, y_max;
  if (atomLabels_[atom->getIdx()]) {
    x_min = y_min = std::numeric_limits<double>::max();
    x_max = y_max = std::numeric_limits<double>::lowest();
    atomLabels_[atom->getIdx()]->findExtremes(x_min, x_max, y_min, y_max);
  } else {
    x_min = atCds.x - 3 * spot_rad;

    x_max = atCds.x + 3 * spot_rad;
    y_min = atCds.y - 3 * spot_rad;
    y_max = atCds.y + 3 * spot_rad;
  }

  auto try_north = [&]() -> bool {
    rad_rect.width_ = rad_size * fontScale_;
    rad_rect.height_ = spot_rad * 3.0;
    rad_rect.trans_.x = atCds.x;
    rad_rect.trans_.y = y_max + 0.5 * rad_rect.height_;
    return !doesRectClash(rad_rect, 0.0);
  };
  auto try_south = [&]() -> bool {
    rad_rect.width_ = rad_size * fontScale_;
    rad_rect.height_ = spot_rad * 3.0;
    rad_rect.trans_.x = atCds.x;
    rad_rect.trans_.y = y_min - 0.5 * rad_rect.height_;
    return !doesRectClash(rad_rect, 0.0);
  };
  auto try_east = [&]() -> bool {
    rad_rect.trans_.x = x_max + 3.0 * spot_rad;
    rad_rect.trans_.y = atCds.y;
    rad_rect.width_ = spot_rad * 1.5;
    rad_rect.height_ = rad_size * fontScale_;
    return !doesRectClash(rad_rect, 0.0);
  };
  auto try_west = [&]() -> bool {
    rad_rect.trans_.x = x_min - 3.0 * spot_rad;
    rad_rect.trans_.y = atCds.y;
    rad_rect.width_ = spot_rad * 1.5;
    rad_rect.height_ = rad_size * fontScale_;
    return !doesRectClash(rad_rect, 0.0);
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
void DrawMol::getDrawTransformers(Point2D &trans, Point2D &scale,
                                  Point2D &toCentre) const {
  trans = Point2D(-xMin_, -yMin_);
  scale = Point2D(scale_, scale_);
  Point2D scaledRanges(scale_ * xRange_, scale_ * yRange_);
  toCentre = Point2D(
      (drawWidth_ - scaledRanges.x) / 2.0 + xOffset_ + width_ * marginPadding_,
      (molHeight_ - scaledRanges.y) / 2.0 + yOffset_ +
          height_ * marginPadding_);
}

// ****************************************************************************
Point2D DrawMol::getDrawCoords(const Point2D &atCds, const Point2D &trans,
                               const Point2D &scaleFactor,
                               const Point2D &toCentre) const {
  // we always invert y
  Point2D drawCoords{atCds.x, -atCds.y};
  drawCoords += trans;
  drawCoords.x *= scaleFactor.x;
  drawCoords.y *= scaleFactor.y;
  drawCoords += toCentre;
  return drawCoords;
}

// ****************************************************************************
Point2D DrawMol::getDrawCoords(const Point2D &atCds) const {
  // we always invert y
  Point2D drawCoords{atCds.x, -atCds.y};
  Point2D trans, scale, toCentre;
  getDrawTransformers(trans, scale, toCentre);
  drawCoords += trans;
  drawCoords.x *= scale.x;
  drawCoords.y *= scale.y;
  drawCoords += toCentre;
  return drawCoords;
}

// ****************************************************************************
Point2D DrawMol::getDrawCoords(int atnum) const {
  PRECONDITION(atnum >= 0 && atnum < static_cast<int>(atCds_.size()),
               "bad atom number");
  return getDrawCoords(Point2D(atCds_[atnum].x, -atCds_[atnum].y));
}

// ****************************************************************************
Point2D DrawMol::getAtomCoords(const Point2D &screenCds) const {
  Point2D trans, scale, toCentre;
  getDrawTransformers(trans, scale, toCentre);
  Point2D atCds{screenCds};
  atCds -= toCentre;
  atCds.x /= scale.x;
  atCds.y /= scale.y;
  atCds -= trans;
  // we always invert y
  return Point2D{atCds.x, -atCds.y};
}

// ****************************************************************************
Point2D DrawMol::getAtomCoords(int atnum) const {
  PRECONDITION(atnum >= 0 && atnum < static_cast<int>(atCds_.size()),
               "bad atom number");
  return atCds_[atnum];
}

// ****************************************************************************
void DrawMol::setScale(double newScale, double newFontScale,
                       bool ignoreFontLimits) {
  resetEverything();
  fontScale_ = newFontScale / newScale;
  textDrawer_.setFontScale(fontScale_, true);

  extractAll(newScale);
  findExtremes();

  textDrawer_.setFontScale(newFontScale, ignoreFontLimits);
  scale_ = newScale;
  fontScale_ = textDrawer_.fontScale();
  finishCreateDrawObjects();
}

// ****************************************************************************
void DrawMol::setTransformation(const DrawMol &sourceMol) {
  resetEverything();
  double relFontScale = sourceMol.fontScale_ / sourceMol.scale_;
  textDrawer_.setFontScale(relFontScale, true);
  xMin_ = sourceMol.xMin_;
  xMax_ = sourceMol.xMax_;
  yMin_ = sourceMol.yMin_;
  yMax_ = sourceMol.yMax_;
  xRange_ = sourceMol.xRange_;
  yRange_ = sourceMol.yRange_;

  extractAll(scale_);
  scale_ = sourceMol.scale_;
  fontScale_ = sourceMol.fontScale_;
  textDrawer_.setFontScale(fontScale_, true);
  finishCreateDrawObjects();
}

// ****************************************************************************
void DrawMol::setOffsets(double xOffset, double yOffset) {
  // Remove the existing offsets.  Presumably this will accumulate small
  // errors if it's done a lot.
  if (fabs(xOffset_) > 1e-4 || fabs(yOffset_) > 1e-4) {
    Point2D trans{-xOffset_, -yOffset_};
    transformAll(&trans, nullptr, nullptr);
  }
  xOffset_ = xOffset;
  yOffset_ = yOffset;
  Point2D trans{xOffset_, yOffset_};
  transformAll(&trans, nullptr, nullptr);
}

// ****************************************************************************
void DrawMol::tagAtomsWithCoords() {
  auto tag = boost::str(boost::format("_atomdrawpos_%d") % confId_);
  for (unsigned int j = 0; j < drawMol_->getNumAtoms(); ++j) {
    drawMol_->getAtomWithIdx(j)->setProp(tag, atCds_[j], true);
  }
}

// ****************************************************************************
void DrawMol::transformAll(const Point2D *trans, Point2D *scale,
                           const Point2D *toCentre) {
  for (auto &ps : preShapes_) {
    if (trans) {
      ps->move(*trans);
    }
    if (scale) {
      ps->scale(*scale);
    }
    if (toCentre) {
      ps->move(*toCentre);
    }
  }
  for (auto &bond : bonds_) {
    if (trans) {
      bond->move(*trans);
    }
    if (scale) {
      bond->scale(*scale);
    }
    if (toCentre) {
      bond->move(*toCentre);
    }
  }
  for (auto &hl : highlights_) {
    if (trans) {
      hl->move(*trans);
    }
    if (scale) {
      hl->scale(*scale);
    }
    if (toCentre) {
      hl->move(*toCentre);
    }
  }
  for (auto &annot : annotations_) {
    if (trans) {
      annot->move(*trans);
    }
    if (scale) {
      annot->scale(*scale);
    }
    if (toCentre) {
      annot->move(*toCentre);
    }
  }
  for (auto &label : atomLabels_) {
    if (label) {
      if (trans) {
        label->move(*trans);
      }
      if (scale) {
        label->scale(*scale);
      }
      if (toCentre) {
        label->move(*toCentre);
      }
    }
  }

  // radicals are based on StringRect so don't have their own class.
  // They need to be moved according to scale and scaled according to
  // fontscale.
  for (auto &rad : radicals_) {
    auto &r = get<0>(rad);
    r.trans_ = transformPoint(r.trans_, trans, scale, toCentre);
    r.width_ *= fontScale_;
    r.height_ *= fontScale_;
  }
  for (auto &ps : postShapes_) {
    if (trans) {
      ps->move(*trans);
    }
    if (scale) {
      ps->scale(*scale);
    }
    if (toCentre) {
      ps->move(*toCentre);
    }
  }
}

// ****************************************************************************
Point2D DrawMol::transformPoint(const Point2D &pt, const Point2D *trans,
                                Point2D *scale, const Point2D *toCentre) const {
  Point2D retPt{pt};
  if (trans) {
    retPt += *trans;
  }
  if (scale) {
    retPt.x *= scale->x;
    retPt.y *= scale->y;
  }
  if (toCentre) {
    retPt += *toCentre;
  }
  return retPt;
}

// ****************************************************************************
void DrawMol::calcDoubleBondLines(double offset, const Bond &bond, Point2D &l1s,
                                  Point2D &l1f, Point2D &l2s,
                                  Point2D &l2f) const {
  Atom *at1 = bond.getBeginAtom();
  Atom *at2 = bond.getEndAtom();
  Point2D perp;

  if (isLinearAtom(*at1, atCds_) || isLinearAtom(*at2, atCds_) ||
      (at1->getDegree() == 1 && at2->getDegree() == 1)) {
    const Point2D &at1_cds = atCds_[at1->getIdx()];
    const Point2D &at2_cds = atCds_[at2->getIdx()];
    perp = calcPerpendicular(at1_cds, at2_cds) * offset * 0.5;
    l1s = at1_cds + perp;
    l1f = at2_cds + perp;
    l2s = at1_cds - perp;
    l2f = at2_cds - perp;
  } else if ((at1->getDegree() == 1 || at2->getDegree() == 1)) {
    doubleBondTerminal(at1, at2, offset, l1s, l1f, l2s, l2f);
  } else {
    const Point2D &at1_cds = atCds_[at1->getIdx()];
    const Point2D &at2_cds = atCds_[at2->getIdx()];
    l1s = at1_cds;
    l1f = at2_cds;
    if (drawMol_->getRingInfo()->numBondRings(bond.getIdx())) {
      // in a ring, we need to draw the bond inside the ring
      bondInsideRing(bond, offset, l2s, l2f);
    } else {
      // if there are atom labels at both ends, straddle the atom-atom vector
      // rather than the normal 1 line on the vector, the other to the side.
      if (atomLabels_[at1->getIdx()] && atomLabels_[at2->getIdx()]) {
        doubleBondTerminal(at1, at2, offset, l1s, l1f, l2s, l2f);
        offset /= 2.0;
      } else {
        bondNonRing(bond, offset, l2s, l2f);
      }
    }

    // Occasionally, as seen in Github6170, a bad geometry about a bond can
    // result in the bonds being crossed as perpendiculars have become
    // confused.  Usually this is the result of when a bond to the
    // double bond is roughly linear with it.  This is a cheap test to see if
    // this has happened, uncrossing them if necessary.
    if (!areBondsParallel(l1s, l1f, l2f, l2s)) {
      std::swap(l1s, l2s);
    }
    if ((Bond::EITHERDOUBLE == bond.getBondDir()) ||
        (Bond::STEREOANY == bond.getStereo())) {
      // crossed bond
      std::swap(l1s, l2s);
    }
  }
}

// ****************************************************************************
// bond is in a ring, assumed to be double.
// Returns in l2s and l2f the start and finish points of the inner line
// of the double bond.
void DrawMol::bondInsideRing(const Bond &bond, double offset, Point2D &l2s,
                             Point2D &l2f) const {
  std::vector<size_t> bond_in_rings;
  const auto &bond_rings = drawMol_->getRingInfo()->bondRings();
  for (size_t i = 0; i < bond_rings.size(); ++i) {
    if (find(bond_rings[i].begin(), bond_rings[i].end(), bond.getIdx()) !=
        bond_rings[i].end()) {
      bond_in_rings.push_back(i);
    }
  }

  // given the bond and the atom at one end, find the ring atom connected to it
  // that isn't the other end of the bond.
  auto other_ring_atom = [&](unsigned int bondAtom, const Bond &bond,
                             const INT_VECT &ringBonds) -> int {
    auto atom = drawMol_->getAtomWithIdx(bondAtom);
    for (const auto bond2 : drawMol_->atomBonds(atom)) {
      if (bond2->getIdx() == bond.getIdx()) {
        continue;
      }
      if (find(ringBonds.begin(), ringBonds.end(), bond2->getIdx()) !=
          ringBonds.end()) {
        return bond2->getOtherAtomIdx(bondAtom);
      }
    }
    return -1;
  };

  const std::vector<int> *ringToUse = nullptr;
  if (bond_in_rings.size() > 1) {
    // bond is in more than 1 ring.  Choose one that is the same aromaticity
    // as the bond, so that if bond is aromatic, the double bond is inside
    // the aromatic ring.  This is important for morphine, for example,
    // where there are fused aromatic and aliphatic rings.
    // morphine: CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5
    for (size_t i = 0; i < bond_in_rings.size(); ++i) {
      ringToUse = &bond_rings[bond_in_rings[i]];
      bool ring_ok = true;
      for (auto bond_idx : *ringToUse) {
        const Bond *bond2 = drawMol_->getBondWithIdx(bond_idx);
        if (bond.getIsAromatic() != bond2->getIsAromatic()) {
          ring_ok = false;
          break;
        }
      }
      if (ring_ok) {
        break;
      }
    }
  } else {
    ringToUse = &bond_rings[bond_in_rings.front()];
  }

  // either bond is in 1 ring, or we couldn't decide above, so just use the
  // first one
  int thirdAtom = other_ring_atom(bond.getBeginAtomIdx(), bond, *ringToUse);
  int fourthAtom = other_ring_atom(bond.getEndAtomIdx(), bond, *ringToUse);
  // As seen in #5486, bonds in rings can be trans and the default code assumes
  // they are always cis.  If trans, treat as a non-ring bond.  It won't
  // necessarily come out on the inside of the ring, but that's quite
  // complicated to fix at this point.
  int begIdx = bond.getBeginAtomIdx();
  int endIdx = bond.getEndAtomIdx();
  bool isTrans = areBondsTrans(atCds_[thirdAtom], atCds_[begIdx],
                               atCds_[endIdx], atCds_[fourthAtom]);
  if (isTrans) {
    bondNonRing(bond, offset, l2s, l2f);
  } else {
    l2s = doubleBondEnd(thirdAtom, begIdx, endIdx, offset,
                        !bool(atomLabels_[bond.getBeginAtomIdx()]));
    l2f = doubleBondEnd(fourthAtom, endIdx, begIdx, offset,
                        !bool(atomLabels_[bond.getEndAtomIdx()]));
  }
}

// ****************************************************************************
// bond is not in a ring, assumed to be double.
// Returns in l2s and l2f the start and finish points of the inner line
void DrawMol::bondNonRing(const Bond &bond, double offset, Point2D &l2s,
                          Point2D &l2f) const {
  auto begAt = bond.getBeginAtom();
  auto endAt = bond.getEndAtom();
  const Atom *thirdAtom = nullptr;
  const Atom *fourthAtom = nullptr;

  bool begTrunc = !atomLabels_[begAt->getIdx()];
  bool endTrunc = !atomLabels_[endAt->getIdx()];

  // find a neighbour of at1 that isn't at2 and if possible isn't directly
  // opposite at1 to at2.
  auto nonColinearNbor = [&](Atom *at1, Atom *at2) -> const Atom * {
    const Atom *thirdAtom = nullptr;
    for (auto i = 1u; i < at1->getDegree(); ++i) {
      thirdAtom = otherNeighbor(at1, at2, i, *drawMol_);
      if (thirdAtom &&
          !areBondsParallel(atCds_[at1->getIdx()], atCds_[at2->getIdx()],
                            atCds_[at1->getIdx()],
                            atCds_[thirdAtom->getIdx()])) {
        return thirdAtom;
      }
    }
    if (thirdAtom == nullptr) {
      // we need something.
      thirdAtom = otherNeighbor(at1, at2, 1, *drawMol_);
    }
    return thirdAtom;
  };

  if (begAt->getDegree() == 2 && endAt->getDegree() == 2) {
    thirdAtom = otherNeighbor(begAt, endAt, 0, *drawMol_);
    fourthAtom = otherNeighbor(endAt, begAt, 0, *drawMol_);
    l2s = doubleBondEnd(thirdAtom->getIdx(), begAt->getIdx(), endAt->getIdx(),
                        offset, begTrunc);
    bool isTrans =
        areBondsTrans(atCds_[thirdAtom->getIdx()], atCds_[begAt->getIdx()],
                      atCds_[endAt->getIdx()], atCds_[fourthAtom->getIdx()]);
    if (isTrans) {
      Point2D perp = calcInnerPerpendicular(atCds_[endAt->getIdx()],
                                            atCds_[begAt->getIdx()],
                                            atCds_[thirdAtom->getIdx()]);
      l2f = atCds_[endAt->getIdx()] + perp * offset;
    } else {
      l2f = doubleBondEnd(fourthAtom->getIdx(), endAt->getIdx(),
                          begAt->getIdx(), offset, endTrunc);
    }
  } else if (begAt->getDegree() == 2 && endAt->getDegree() > 2) {
    thirdAtom = otherNeighbor(begAt, endAt, 0, *drawMol_);
    fourthAtom = otherNeighbor(endAt, begAt, 0, *drawMol_);
    l2s = doubleBondEnd(thirdAtom->getIdx(), begAt->getIdx(), endAt->getIdx(),
                        offset, begTrunc);
    bool isTrans =
        areBondsTrans(atCds_[thirdAtom->getIdx()], atCds_[begAt->getIdx()],
                      atCds_[endAt->getIdx()], atCds_[fourthAtom->getIdx()]);
    if (isTrans) {
      fourthAtom = nonColinearNbor(endAt, begAt);
    }
    l2f = doubleBondEnd(fourthAtom->getIdx(), endAt->getIdx(), begAt->getIdx(),
                        offset, endTrunc);
  } else if (begAt->getDegree() > 2 && endAt->getDegree() == 2) {
    thirdAtom = otherNeighbor(begAt, endAt, 0, *drawMol_);
    fourthAtom = otherNeighbor(endAt, begAt, 0, *drawMol_);
    l2s = doubleBondEnd(thirdAtom->getIdx(), begAt->getIdx(), endAt->getIdx(),
                        offset, begTrunc);
    bool isTrans =
        areBondsTrans(atCds_[thirdAtom->getIdx()], atCds_[begAt->getIdx()],
                      atCds_[endAt->getIdx()], atCds_[fourthAtom->getIdx()]);
    if (isTrans) {
      thirdAtom = nonColinearNbor(begAt, endAt);
      l2s = doubleBondEnd(thirdAtom->getIdx(), begAt->getIdx(), endAt->getIdx(),
                          offset, endTrunc);
    }
    l2f = doubleBondEnd(fourthAtom->getIdx(), endAt->getIdx(), begAt->getIdx(),
                        offset, endTrunc);
  } else if (begAt->getDegree() > 2 && endAt->getDegree() > 2) {
    thirdAtom = otherNeighbor(begAt, endAt, 0, *drawMol_);
    l2s = doubleBondEnd(thirdAtom->getIdx(), begAt->getIdx(), endAt->getIdx(),
                        offset, begTrunc);
    fourthAtom = otherNeighbor(endAt, begAt, 0, *drawMol_);
    bool isTrans =
        areBondsTrans(atCds_[thirdAtom->getIdx()], atCds_[begAt->getIdx()],
                      atCds_[endAt->getIdx()], atCds_[fourthAtom->getIdx()]);
    if (isTrans) {
      fourthAtom = nonColinearNbor(endAt, begAt);
    }
    l2f = doubleBondEnd(fourthAtom->getIdx(), endAt->getIdx(), begAt->getIdx(),
                        offset, endTrunc);
  }
}

// ****************************************************************************
void DrawMol::doubleBondTerminal(Atom *at1, Atom *at2, double offset,
                                 Point2D &l1s, Point2D &l1f, Point2D &l2s,
                                 Point2D &l2f) const {
  bool swapped = false;
  if (at1->getDegree() > 1 && at2->getDegree() == 1) {
    std::swap(at1, at2);
    swapped = true;
  }
  const Point2D &at1_cds = atCds_[at1->getIdx()];
  const Point2D &at2_cds = atCds_[at2->getIdx()];
  if (atomLabels_[at2->getIdx()]) {
    // either side of the bond line if going ot a label
    offset /= 2.0;
    Point2D perp = calcPerpendicular(at1_cds, at2_cds) * offset;
    l1s = at1_cds + perp;
    l1f = at2_cds + perp;
    l2s = at1_cds - perp;
    l2f = at2_cds - perp;
  } else if (at2->getDegree() > 2) {
    // lines either side of the bond line but at the at2 end,
    // the bonds extend to the intersection of the other bonds.
    // only need 1/2 the offset in this case.
    offset /= 2.0;
    Point2D perp = calcPerpendicular(at1_cds, at2_cds) * offset;
    l1s = at1_cds + perp;
    l1f = at2_cds + perp;
    l2s = at1_cds - perp;
    l2f = at2_cds - perp;
    // extend the two bond lines so they will intersect with the other bonds
    // from at2
    auto bl = std::max((l1s - l1f).length(), (l2s - l2f).length());
    Point2D l1 = l1s.directionVector(l1f);
    l1f = l1s + l1 * 2.0 * bl;
    Point2D l2 = l2s.directionVector(l2f);
    l2f = l2s + l2 * 2.0 * bl;
    Point2D ip;
    for (auto nbr : make_iterator_range(drawMol_->getAtomNeighbors(at2))) {
      auto nbr_cds = atCds_[nbr];
      if (doLinesIntersect(l1s, l1f, at2_cds, nbr_cds, &ip)) {
        l1f = ip;
      }
      if (doLinesIntersect(l2s, l2f, at2_cds, nbr_cds, &ip)) {
        l2f = ip;
      }
    }
  } else {
    // one line as normal, the 2nd truncates at the internal end only
    l1s = at1_cds;
    l1f = at2_cds;
    const Atom *thirdAtom = otherNeighbor(at2, at1, 0, *drawMol_);
    Point2D perp =
        calcInnerPerpendicular(at1_cds, at2_cds, atCds_[thirdAtom->getIdx()]);
    l2s = at1_cds + perp * offset;
    l2f = doubleBondEnd(at1->getIdx(), at2->getIdx(), thirdAtom->getIdx(),
                        offset, true);
    // If at1->at2->at3 is a straight line, l2f may have ended up on the
    // wrong side of the other bond from l2s because there is no inner
    // side of the bond.  Do it again with a negative offset if so.
    if (fabs(l1s.directionVector(l1f).dotProduct(l2s.directionVector(l2f))) <
        0.9999) {
      l2f = doubleBondEnd(at1->getIdx(), at2->getIdx(), thirdAtom->getIdx(),
                          -offset, true);
    }
    // if at1 has a label, need to move it so it's centred in between the
    // two lines (Github 5511).
    if (atomLabels_[at1->getIdx()]) {
      atomLabels_[at1->getIdx()]->cds_ += perp * offset * 0.5;
    }
  }
  if (swapped) {
    std::swap(l1s, l1f);
    std::swap(l2s, l2f);
  }
}

// ****************************************************************************
Point2D DrawMol::doubleBondEnd(unsigned int at1, unsigned int at2,
                               unsigned int at3, double offset,
                               bool trunc) const {
  Point2D v21 = atCds_[at2].directionVector(atCds_[at1]);
  Point2D v23 = atCds_[at2].directionVector(atCds_[at3]);
  Point2D v23perp(-v23.y, v23.x);
  v23perp.normalize();

  Point2D bis = v21 + v23;
  if (bis.lengthSq() < 1.0e-6) {
    // if the bonds are colinear, bis comes out as 0, and thus normalizes
    // to NaN which gives a very ugly result (Github #6027).  It's safe
    // to use v23perp in this case, so long as is on the right side of the
    // bond, which will be checked on return.
    return (atCds_[at2] - v23perp * offset);
  }

  bis.normalize();
  if (v23perp.dotProduct(bis) < 0.0) {
    v23perp = v23perp * -1.0;
  }
  Point2D ip;
  // if there's an atom label, we don't need to step the bond end back
  // because both ends are shortened to accommodate the letters.
  // likewise if the two lines don't intersect, it's already stepped
  // back enough (github 6025).
  bool ipAlreadySet = false;
  if (trunc) {
    ipAlreadySet = doLinesIntersect(atCds_[at2], atCds_[at2] + bis,
                                    atCds_[at2] + v23perp * offset,
                                    atCds_[at3] + v23perp * offset, &ip);
  }
  if (!ipAlreadySet) {
    ip = atCds_[at2] + v23perp * offset;
  }
  return ip;
}

// ****************************************************************************
void DrawMol::calcTripleBondLines(double offset, const Bond &bond, Point2D &l1s,
                                  Point2D &l1f, Point2D &l2s, Point2D &l2f) {
  Atom *at1 = bond.getBeginAtom();
  Atom *at2 = bond.getEndAtom();
  const Point2D &at1_cds = atCds_[at1->getIdx()];
  const Point2D &at2_cds = atCds_[at2->getIdx()];

  Point2D perp = calcPerpendicular(at1_cds, at2_cds);
  l1s = at1_cds + perp * offset;
  l1f = at2_cds + perp * offset;
  l2s = at1_cds - perp * offset;
  l2f = at2_cds - perp * offset;
}

// ****************************************************************************
void DrawMol::findOtherBondVecs(const Atom *atom, const Atom *otherAtom,
                                std::vector<Point2D> &otherBondVecs) const {
  if (atom->getDegree() == 1 || atomLabels_[atom->getIdx()]) {
    return;
  }
  for (unsigned int i = 1; i < atom->getDegree(); ++i) {
    auto thirdAtom = otherNeighbor(atom, otherAtom, i - 1, *drawMol_);
    Point2D const &at1_cds = atCds_[atom->getIdx()];
    Point2D const &at2_cds = atCds_[thirdAtom->getIdx()];
    otherBondVecs.push_back(at1_cds.directionVector(at2_cds));
  }
}

// ****************************************************************************
void DrawMol::adjustBondsOnSolidWedgeEnds() {
  for (auto &bond : drawMol_->bonds()) {
    if (bond->getBondDir() == Bond::BEGINWEDGE &&
        bond->getEndAtom()->getDegree() == 2 &&
        !atomLabels_[bond->getEndAtomIdx()]) {
      // find the bond at the end atom
      auto thirdAtom =
          otherNeighbor(bond->getEndAtom(), bond->getBeginAtom(), 0, *drawMol_);
      auto bond1 = drawMol_->getBondBetweenAtoms(bond->getEndAtomIdx(),
                                                 thirdAtom->getIdx());
      // If the bonds a co-linear, don't do anything (Github7036)
      auto b1 = atCds_[bond->getEndAtomIdx()].directionVector(
          atCds_[bond->getBeginAtomIdx()]);
      auto b2 = atCds_[bond1->getEndAtomIdx()].directionVector(
          atCds_[bond1->getBeginAtomIdx()]);
      if (fabs(1.0 - b1.dotProduct(b2)) < 0.001) {
        continue;
      }
      DrawShape *wedge = nullptr;
      DrawShape *bondLine = nullptr;
      double closestDist = 1.0;
      for (auto &shape : bonds_) {
        if (shape->bond_ == static_cast<int>(bond->getIdx())) {
          wedge = shape.get();
        }
        // there may be multiple lines for the bond, so we want one that
        // has an end as close as possible to the bond end atom.
        auto endCds = atCds_[bond->getEndAtomIdx()];
        if (shape->bond_ == static_cast<int>(bond1->getIdx())) {
          // only deal with simple lines, which I think should be the only
          // case, but...
          if (dynamic_cast<DrawShapeSimpleLine *>(shape.get()) == nullptr) {
            continue;
          }
          if ((shape->points_[0] - endCds).lengthSq() < closestDist) {
            closestDist = (shape->points_[0] - endCds).lengthSq();
            bondLine = shape.get();
          }
          if ((shape->points_[1] - endCds).lengthSq() < closestDist) {
            closestDist = (shape->points_[1] - endCds).lengthSq();
            bondLine = shape.get();
          }
        }
      }
      if (wedge != nullptr && bondLine != nullptr) {
        int p1 = -1, p2 = -1;
        // find the points that are the top of the wedge.  Clearly, this
        // assumes the order that the triangles are created in the
        // DrawShapeSolidWedge.
        if (wedge->points_.size() == 3) {
          p1 = 1;
          p2 = 2;
        } else if (wedge->points_.size() == 9) {
          p1 = 4;
          p2 = 5;
        }
        // want the p1 or p2 that is furthest from the 3rd atom - make it p1
        if (p1 != -1 && p2 != -1) {
          if ((atCds_[thirdAtom->getIdx()] - wedge->points_[p1]).lengthSq() <
              (atCds_[thirdAtom->getIdx()] - wedge->points_[p2]).lengthSq()) {
            std::swap(p1, p2);
          }
          // now make the coords of the end of the bondLine that matches p1 the
          // same as p1
          if (bondLine->atom1_ == wedge->atom2_) {
            bondLine->points_[0] = wedge->points_[p1];
          } else {
            bondLine->points_[1] = wedge->points_[p1];
          }
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMol::smoothBondJoins() {
  // Because the bonds are drawn as individual lines rather than as paths
  // through the molecule, they don't join up nicely.  Put a little path
  // round the join where it's needed to hide the problem.
  // The bonds aren't drawn as paths because in SVGs each line is given
  // classes for the atoms and bond it involves, and people use this to
  // identify the lines for other purposes.
  for (auto atom : drawMol_->atoms()) {
    bool doIt = false;
    if (atom->getDegree() == 2) {
      doIt = true;
    } else if (atom->getDegree() == 3) {
      for (const auto nbr : drawMol_->atomNeighbors(atom)) {
        auto bond =
            drawMol_->getBondBetweenAtoms(atom->getIdx(), nbr->getIdx());
        if ((nbr->getDegree() == 1 && bond->getBondType() == Bond::DOUBLE) ||
            bond->getBondDir() == Bond::BEGINDASH ||
            bond->getBondDir() == Bond::BEGINWEDGE) {
          doIt = true;
        }
      }
    }
    if (doIt) {
      bool done = false;
      for (unsigned int i = 0; i < singleBondLines_.size(); ++i) {
        auto &sbl1 = bonds_[singleBondLines_[i]];
        int p1 = -1;
        int p2 = -1;
        if (static_cast<int>(atom->getIdx()) == sbl1->atom1_) {
          p1 = 0;
        } else if (static_cast<int>(atom->getIdx()) == sbl1->atom2_) {
          p1 = 1;
        }
        if (p1 != -1) {
          for (unsigned int j = 0; j < singleBondLines_.size(); ++j) {
            if (i == j) {
              continue;
            }
            auto &sbl2 = bonds_[singleBondLines_[j]];
            if (static_cast<int>(atom->getIdx()) == sbl2->atom1_) {
              p2 = 0;
            } else if (static_cast<int>(atom->getIdx()) == sbl2->atom2_) {
              p2 = 1;
            }
            if (p2 != -1) {
              double dist = (sbl1->points_[p1] - sbl2->points_[p2]).lengthSq();
              if (dist < 1.0e-6) {
                // make a small polyline to paper over the cracks.
                int p12 = p1 == 1 ? 0 : 1;
                int p22 = p2 == 1 ? 0 : 1;
                // If the lines are different colours, make the line round
                // the corner shorter so that one colour doesn't extend
                // into the other one.  If they're the same colour, it's
                // better if they go round the corner a bit further to hide
                // the join better.  The numbers are empirical.
                double len =
                    sbl1->lineColour_ == sbl2->lineColour_ ? 0.05 : 0.025;
                Point2D dv1 = (sbl1->points_[p1] - sbl1->points_[p12]) * len;
                Point2D dv2 = (sbl1->points_[p1] - sbl2->points_[p22]) * len;
                std::vector<Point2D> pl_pts{sbl1->points_[p1] - dv1,
                                            sbl1->points_[p1],
                                            sbl1->points_[p1] - dv2};
                DrawShape *pl = new DrawShapePolyLine(pl_pts, sbl1->lineWidth_,
                                                      sbl1->scaleLineWidth_,
                                                      sbl1->lineColour_);
                bonds_.emplace_back(pl);
                done = true;
                break;
              }
            }
          }
        }
        if (done) {
          break;
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMol::makeHighlightEnd(const Atom *end1, const Atom *end2,
                               double lineWidth,
                               const std::vector<Atom *> &end1HighNbrs,
                               std::vector<Point2D> &points) {
  double halfLineWidth = lineWidth / 2.0;
  // find the intersection point of two lines parallel to lines from e2 to e1
  // and e3 to e1 and lineWidth from them.  If pm is 1, it's inside the
  // angle they make, if pm is -1, it's outside.  If the lines don't
  // intersect, it returns e1.
  auto innerPoint = [&](Point2D &e1, Point2D &e2, Point2D &e3,
                        double pm) -> Point2D {
    auto perp1 = calcInnerPerpendicular(e2, e1, e3);
    auto perp2 = calcInnerPerpendicular(e3, e1, e2);
    auto line12 = e2 + perp1 * pm * halfLineWidth;
    auto line11 = e1 + perp1 * pm * halfLineWidth;
    line11 = line12 + line12.directionVector(line11) * 2.0 * (e1 - e2).length();
    auto line22 = e3 + perp2 * pm * halfLineWidth;
    auto line21 = e1 + perp2 * pm * halfLineWidth;
    line21 = line22 + line22.directionVector(line21) * 2.0 * (e1 - e3).length();
    Point2D ins;
    if (doLinesIntersect(line12, line11, line22, line21, &ins)) {
      return ins;
    } else {
      return Point2D(e1);
    }
  };

  auto end1Cds = atCds_[end1->getIdx()];
  auto end2Cds = atCds_[end2->getIdx()];

  if (end1HighNbrs.empty()) {
    // If end1 doesn't have any highlighted neighbour bonds, then
    // it's a flat end.
    auto perp = calcPerpendicular(end1Cds, end2Cds);
    points.push_back(end1Cds + perp * halfLineWidth);
    points.push_back(end1Cds - perp * halfLineWidth);
  } else if (end1HighNbrs.size() == 1) {
    // There is only 1 intersection to deal with, which is easier - just
    // a slanted end.
    auto end3Cds = atCds_[end1HighNbrs[0]->getIdx()];
    auto b1 = end2Cds.directionVector(end1Cds);
    auto b2 = end2Cds.directionVector(end3Cds);
    if (1.0 - fabs(b1.dotProduct(b2)) < 1.0e-4) {
      // move end3 by a small amount to create an inner and outer
      auto d32 = end3Cds - end2Cds;
      Point2D d32transp(d32.y, -d32.x);
      d32transp *= 0.1;
      end3Cds += d32transp;
    }
    // The moved end is only used to construct ins1 and ins2 wrt
    // end1Cds and end2Cds so there's no need do anything more.
    auto ins1 = innerPoint(end1Cds, end2Cds, end3Cds, 1.0);
    points.push_back(ins1);
    auto ins2 = innerPoint(end1Cds, end2Cds, end3Cds, -1.0);
    points.push_back(ins2);
  } else if (end1HighNbrs.size() > 1) {
    // In this case, it needs a triangular end, as it's a junction
    // of at least 3 highlights.  The point of the triangle is
    // end1. The other points are defined by the first and last bond
    // vectors going round from the end1->end2 vector, so sort the
    // neighbours in order of increasing angle made with the end2->end1
    // vector.
    auto bvec = end1Cds.directionVector(end2Cds);
    std::vector<std::pair<int, double>> angs;
    for (unsigned i = 0; i < end1HighNbrs.size(); ++i) {
      auto ovec = end1Cds.directionVector(atCds_[end1HighNbrs[i]->getIdx()]);
      auto ang = bvec.signedAngleTo(ovec);
      angs.push_back(std::make_pair(i, ang));
    }
    std::sort(angs.begin(), angs.end(),
              [](const std::pair<int, double> &a1,
                 const std::pair<int, double> &a2) -> bool {
                return a1.second < a2.second;
              });
    // if both angles are on the same side as end1->end2, they need to
    // be the other way round.
    if (angs.front().second > M_PI && angs.back().second > M_PI) {
      std::reverse(angs.begin(), angs.end());
    }
    auto end3Cds = atCds_[end1HighNbrs[angs.front().first]->getIdx()];
    auto ins1 = innerPoint(end1Cds, end2Cds, end3Cds, 1.0);
    points.push_back(ins1);
    points.push_back(end1Cds);
    auto end4Cds = atCds_[end1HighNbrs[angs.back().first]->getIdx()];
    // if both angles are on the same side as end1->end2, they need to
    // be the other way round.
    double pm = 1.0;
    if ((angs.front().second > M_PI && angs.back().second > M_PI) ||
        (angs.front().second < M_PI && angs.back().second < M_PI)) {
      pm = -1.0;
    }
    auto ins2 = innerPoint(end1Cds, end2Cds, end4Cds, pm);
    points.push_back(ins2);
  }
}

// ****************************************************************************
DrawColour DrawMol::getColour(int atom_idx) const {
  PRECONDITION(atom_idx >= 0, "bad atom_idx");
  PRECONDITION(rdcast<int>(atomicNums_.size()) > atom_idx, "bad atom_idx");

  DrawColour retval = getColourByAtomicNum(atomicNums_[atom_idx], drawOptions_);
  bool highlightedAtom = false;
  if (!drawOptions_.circleAtoms && !drawOptions_.continuousHighlight) {
    if (highlightAtoms_.end() !=
        find(highlightAtoms_.begin(), highlightAtoms_.end(), atom_idx)) {
      highlightedAtom = true;
      retval = drawOptions_.highlightColour;
    }
    // over-ride with explicit colour from highlightMap if there is one
    auto p = highlightAtomMap_.find(atom_idx);
    if (p != highlightAtomMap_.end()) {
      highlightedAtom = true;
      retval = p->second;
    }
    // if it's not a highlighted atom itself, but all the bonds off it
    // are highlighted, I think it's better if the atom itself adopts
    // the same highlight colour as the bonds.  It doesn't look right
    // if only some of the bonds are highlighted, IMO.
    if (!highlightedAtom) {
      Atom *atomPtr = drawMol_->getAtomWithIdx(atom_idx);
      int numBonds = 0, numHighBonds = 0;
      std::unique_ptr<DrawColour> highCol;
      for (const auto &nbri :
           boost::make_iterator_range(drawMol_->getAtomBonds(atomPtr))) {
        ++numBonds;
        const auto &nbr = (*drawMol_)[nbri];
        if (std::find(highlightBonds_.begin(), highlightBonds_.end(),
                      nbr->getIdx()) != highlightBonds_.end() ||
            highlightBondMap_.find(nbr->getIdx()) != highlightBondMap_.end()) {
          DrawColour hc = getHighlightBondColour(
              nbr, drawOptions_, highlightBonds_, highlightBondMap_,
              highlightAtoms_, highlightAtomMap_);
          if (!highCol) {
            highCol.reset(new DrawColour(hc));
          } else {
            if (!(hc == *highCol)) {
              numHighBonds = 0;
              break;
            }
          }
          ++numHighBonds;
        }
      }
      if (numBonds == numHighBonds && highCol) {
        retval = *highCol;
      }
    }
  }
  return retval;
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
bool isLinearAtom(const Atom &atom, const std::vector<Point2D> &atCds) {
  if (atom.getDegree() == 2) {
    Point2D bond_vecs[2];
    Bond::BondType bts[2];
    Point2D const &at1_cds = atCds[atom.getIdx()];
    ROMol const &mol = atom.getOwningMol();
    int i = 0;
    for (auto nbr : make_iterator_range(mol.getAtomNeighbors(&atom))) {
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
DrawColour getColourByAtomicNum(int atomic_num,
                                const MolDrawOptions &drawOptions) {
  DrawColour res;
  if (atomic_num == 1 && drawOptions.noAtomLabels) {
    atomic_num = 201;
  }
  if (drawOptions.atomColourPalette.find(atomic_num) !=
      drawOptions.atomColourPalette.end()) {
    res = drawOptions.atomColourPalette.find(atomic_num)->second;
  } else if (atomic_num != -1 && drawOptions.atomColourPalette.find(-1) !=
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
DrawColour getHighlightBondColour(
    const Bond *bond, const MolDrawOptions &drawOptions,
    const std::vector<int> &highlightBonds,
    const std::map<int, DrawColour> &highlightBondMap,
    const std::vector<int> &highlightAtoms,
    const std::map<int, DrawColour> &highlightAtomMap) {
  PRECONDITION(bond, "no bond provided");
  RDUNUSED_PARAM(highlightAtoms);

  DrawColour col(0.0, 0.0, 0.0);
  if (std::find(highlightBonds.begin(), highlightBonds.end(), bond->getIdx()) !=
      highlightBonds.end()) {
    col = drawOptions.highlightColour;
    if (highlightBondMap.find(bond->getIdx()) != highlightBondMap.end()) {
      col = highlightBondMap.find(bond->getIdx())->second;
    } else {
      // the highlight color of the bond is not explicitly provided. What about
      // the highlight colors of the begin/end atoms? Ideally these will both be
      // the same, but we want to set the coloring even if that's not the
      // case, so we'll use:
      //  - begin atom color if that is set
      //  - end atom color if that is set
      //  - the default highlight color otherwise
      if (highlightAtomMap.find(bond->getBeginAtomIdx()) !=
          highlightAtomMap.end()) {
        col = highlightAtomMap.find(bond->getBeginAtomIdx())->second;

      } else if (highlightAtomMap.find(bond->getEndAtomIdx()) !=
                 highlightAtomMap.end()) {
        col = highlightAtomMap.find(bond->getEndAtomIdx())->second;
      }
    }
  }
  return col;
}

// ****************************************************************************
double getHighlightBondWidth(
    const MolDrawOptions &drawOptions, int bond_idx,
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
  double tgt_lw = drawOptions.bondLineWidth * bwm;

  return tgt_lw;
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
void adjustBondEndForString(
    const Point2D &end2, double padding,
    const std::vector<std::shared_ptr<StringRect>> &rects, Point2D &moveEnd) {
  Point2D labelPos = moveEnd;
  for (auto r : rects) {
    Point2D origTrans = r->trans_;
    r->trans_ += labelPos;

    Point2D tl, tr, bl, br;
    r->calcCorners(tl, tr, br, bl, padding);

    // if it's a wide label, such as C:7, the bond can intersect
    // more than 1 side of the rectangle, so check them all.
    std::unique_ptr<Point2D> ip(new Point2D);
    if (doLinesIntersect(moveEnd, end2, tl, tr, ip.get())) {
      moveEnd = *ip;
    }
    if (doLinesIntersect(moveEnd, end2, tr, br, ip.get())) {
      moveEnd = *ip;
    }
    if (doLinesIntersect(moveEnd, end2, br, bl, ip.get())) {
      moveEnd = *ip;
    }
    if (doLinesIntersect(moveEnd, end2, bl, tl, ip.get())) {
      moveEnd = *ip;
    }
    r->trans_ = origTrans;
  }
}

// ****************************************************************************
void findRadicalExtremes(
    const std::vector<std::tuple<StringRect, OrientType, int>> &radicals,
    double &xmin, double &xmax, double &ymin, double &ymax) {
  for (const auto &rad : radicals) {
    findRectExtremes(get<0>(rad), TextAlignType::MIDDLE, xmin, xmax, ymin,
                     ymax);
  }
}

// ****************************************************************************
void findRectExtremes(const StringRect &rect, const TextAlignType &align,
                      double &xmin, double &xmax, double &ymin, double &ymax) {
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

// ****************************************************************************
void getBondHighlightsForAtoms(const ROMol &mol,
                               const std::vector<int> &highlight_atoms,
                               std::vector<int> &highlight_bonds) {
  highlight_bonds.clear();
  for (auto ai = highlight_atoms.begin(); ai != highlight_atoms.end(); ++ai) {
    for (auto aj = ai + 1; aj != highlight_atoms.end(); ++aj) {
      const Bond *bnd = mol.getBondBetweenAtoms(*ai, *aj);
      if (bnd) {
        highlight_bonds.push_back(bnd->getIdx());
      }
    }
  }
}

// ****************************************************************************
bool areBondsTrans(const Point2D &at1, const Point2D &at2, const Point2D &at3,
                   const Point2D &at4) {
  Point2D v21 = at1 - at2;
  Point2D v34 = at4 - at3;
  return (v21.dotProduct(v34) < 0.0);
}

// ****************************************************************************
bool areBondsParallel(const Point2D &at1, const Point2D &at2,
                      const Point2D &at3, const Point2D &at4, double tol) {
  Point2D v21 = at1.directionVector(at2);
  Point2D v34 = at4.directionVector(at3);
  return (fabs(1.0 - fabs(v21.dotProduct(v34))) < tol);
}

// ****************************************************************************
const Atom *otherNeighbor(const Atom *firstAtom, const Atom *secondAtom,
                          int nborNum, const ROMol &mol) {
  int nbourCount = 0;
  for (const auto nbr : mol.atomNeighbors(firstAtom)) {
    if (nbr->getIdx() != secondAtom->getIdx()) {
      if (nbourCount == nborNum) {
        return nbr;
      } else {
        nbourCount++;
      }
    }
  }
  return nullptr;
};

}  // namespace MolDraw2D_detail
}  // namespace RDKit
