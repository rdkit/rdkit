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
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawShape.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

namespace RDKit {

// ****************************************************************************
DrawMol::DrawMol(const ROMol &mol, const std::string &legend,
                 int width, int height,
                 MolDrawOptions &drawOptions,
                 const std::vector<int> *highlight_atoms,
                 const std::vector<int> *highlight_bonds,
                 const std::map<int, DrawColour> *highlight_atom_map,
                 const std::map<int, DrawColour> *highlight_bond_map,
                 const std::vector<std::pair<DrawColour, DrawColour>> *bond_colours,
                 const std::map<int, double> *highlight_radii,
                 int confId)
    : drawOptions_(drawOptions),
      highlightAtoms_(highlight_atoms),
      highlightBonds_(highlight_bonds),
      highlightAtomMap_(highlight_atom_map),
      highlightBondMap_(highlight_bond_map),
      bondColours_(bond_colours),
      highlightRadii_(highlight_radii),
      width_(width),
      height_(height),
      scale_(1.0),
      xMin_(std::numeric_limits<double>::max() / 2.0),
      yMin_(std::numeric_limits<double>::max() / 2.0),
      xMax_(-std::numeric_limits<double>::max() / 2.0),
      yMax_(-std::numeric_limits<double>::max() / 2.0),
      xRange_(std::numeric_limits<double>::max()),
      yRange_(std::numeric_limits<double>::max()) {
  std::cout << "Top of DrawMol c'tor" << std::endl;
  std::cout << "Width = " << width_ << "  height = " << height_ << std::endl;
  initDrawMolecule(mol, confId);
  extractAll(confId);

  calculateScale();
  changeToDrawCoords();
}

// ****************************************************************************
void DrawMol::initDrawMolecule(const ROMol &mol, int confId) {
  drawMol_.reset(new RWMol(mol));
  if (drawOptions_.prepareMolsBeforeDrawing || !mol.getNumConformers()) {
    MolDraw2DUtils::prepareMolForDrawing(*drawMol_);
  }
  if (drawOptions_.centreMoleculesBeforeDrawing) {
    if (drawMol_->getNumConformers()) {
      centerMolForDrawing(*drawMol_, confId);
    }
  }
  if (drawOptions_.simplifiedStereoGroupLabel &&
      !mol.hasProp(common_properties::molNote)) {
    prepareStereoGroups(*drawMol_);
  }
}

// ****************************************************************************
void DrawMol::extractAll(int confId) {
  extractAtomCoords(confId);
  extractAtomSymbols();
  extractBonds();
}

// ****************************************************************************
void DrawMol::extractBonds() {

  double double_bond_offset = drawOptions_.multipleBondOffset;
  // mol files from, for example, Marvin use a bond length of 1 for just about
  // everything. When this is the case, the default multipleBondOffset is just
  // too much, so scale it back.
  calcMeanBondLengthSquare();
  if (meanBondLengthSquare_ < 1.4) {
    double_bond_offset *= 0.6;
  }

  for (auto bond : drawMol_->bonds()) {
    int beg_at = bond->getBeginAtomIdx();
    int end_at = bond->getEndAtomIdx();
    std::pair<DrawColour, DrawColour> cols = getBondColours(bond);
    int olw = drawOptions_.bondLineWidth;
    newBondLine(atCds_[beg_at], atCds_[end_at], cols.first, cols.second);
    drawOptions_.bondLineWidth = olw;
  }

}

// ****************************************************************************
void DrawMol::extractAtomCoords(int confId) {
  PRECONDITION(static_cast<int>(drawMol_->getNumConformers()) > 0, "no coords");

  const RDGeom::POINT3D_VECT &locs = drawMol_->getConformer(confId).getPositions();

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
    int this_idx = this_at->getIdx();
    Point2D pt(locs[this_idx].x, locs[this_idx].y);
    if (rot != 0.0) {
      trans.TransformPoint(pt);
    }
    atCds_.emplace_back(pt);
    // std::cout << "coords for " << this_idx << " : " << pt << std::endl;
  }
}

// ****************************************************************************
void DrawMol::extractAtomSymbols() {
  atomicNums_.clear();
  for (auto at1 : drawMol_->atoms()) {
    atomSyms_.emplace_back(getAtomSymbolAndOrientation(*at1));
    if (!isComplexQuery(at1)) {
      atomicNums_.emplace_back(at1->getAtomicNum());
    } else {
      atomicNums_.push_back(0);
    }
  }
}

// ****************************************************************************
void DrawMol::calculateScale() {
  findExtremes();

  bool setWidth = false;
  if (width_ < 0) {
    // FIX: technically we need to take the legend width into account too!
    width_ = drawOptions_.scalingFactor * xRange_;
    setWidth = true;
  }
  bool setHeight = false;
  if (height_ < 0) {
    // we need to adjust the range for the legend
    // if it's not present then legend_height_ will be zero and this will be a
    // no-op
    height_ = drawOptions_.scalingFactor * yRange_;
    setHeight = true;
  }

  // put a 5% buffer round the drawing and calculate a final scale
  xMin_ -= drawOptions_.padding * xRange_;
  xRange_ *= 1 + 2 * drawOptions_.padding;
  xMax_ = xMin_ + xRange_;
  yMin_ -= drawOptions_.padding * yRange_;
  yRange_ *= 1 + 2 * drawOptions_.padding;
  yMax_ = yMin_ + yRange_;

  if (xRange_ > 1e-4 || yRange_ > 1e-4) {
    if (setWidth) {
      width_ = drawOptions_.scalingFactor * xRange_;
    }
    if (setHeight) {
      height_ = drawOptions_.scalingFactor * yRange_;
    }

    scale_ = std::min(double(width_) / xRange_, double(height_) / yRange_);
    double fix_scale = scale_;
    // after all that, use the fixed scale unless it's too big, in which case
    // scale the drawing down to fit.
    // fixedScale takes precedence if both it and fixedBondLength are given.
    if (drawOptions_.fixedBondLength > 0.0) {
      fix_scale = drawOptions_.fixedBondLength;
    }
    if (drawOptions_.fixedScale > 0.0) {
      fix_scale = double(width_) * drawOptions_.fixedScale;
    }
    if (scale_ > fix_scale) {
      scale_ = fix_scale;
    }
  } else {
    scale_ = 1;
  }
  std::cout << "Final Scale : " << scale_ << std::endl;
  std::cout << "Padded mins : " << xMin_ << ", " << yMin_ << " with ranges : "
            << xRange_ << ", " << yRange_ << std::endl;
}

// ****************************************************************************
void DrawMol::findExtremes() {
  std::cout << "Number of bonds : " << bonds_.size() << std::endl;
  for(auto &bond: bonds_) {
    bond->findExtremes(xMin_, xMax_, yMin_, yMax_);
  }

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
  Point2D trans(-xMin_, -yMin_);
  Point2D scale(scale_, -scale_);
  Point2D scaledRanges(scale_ * xRange_, scale_ * yRange_);
  Point2D toCentre((width_ - scaledRanges.x) / 2.0,
                   height_ - (height_ - scaledRanges.y) / 2);
  for(auto &bond: bonds_) {
    bond->move(trans);
    bond->scale(scale);
    bond->move(toCentre);
  }
}

// ****************************************************************************
void DrawMol::draw(MolDraw2D &drawer) const {
  for(auto &bond: bonds_) {
    bond->draw(drawer);
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
        orient = OrientType::N;
      } else {
        orient = OrientType::S;
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
            orient = OrientType::N;
            break;
          } else if (ang < -80.0 && ang > -100.0 && orient == OrientType::N) {
            orient = OrientType::S;
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
void DrawMol::newBondLine(const Point2D &pt1, const Point2D &pt2,
                          const DrawColour &col1, const DrawColour &col2) {
  if (col1 == col2) {
    std::vector<Point2D> pts{pt1, pt2};
    DrawShapePolyline *b =
        new DrawShapePolyline(pts, drawOptions_.bondLineWidth, false, col1);
    bonds_.emplace_back(std::unique_ptr<DrawShape>(b));
  } else {
    Point2D mid = (pt1 + pt2) / 2.0;
    std::vector<Point2D> pts1{pt1, mid};
    DrawShapePolyline *b1 =
        new DrawShapePolyline(pts1, drawOptions_.bondLineWidth, false, col1);
    bonds_.emplace_back(std::unique_ptr<DrawShape>(b1));
    std::vector<Point2D> pts2{mid, pt2};
    DrawShapePolyline *b2 =
        new DrawShapePolyline(pts2, drawOptions_.bondLineWidth, false, col2);
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
    if (!highlight_bond) {
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
      if (drawOptions_.continuousHighlight) {
        drawOptions_.bondLineWidth =
            getHighlightBondWidth(drawOptions_, bond->getIdx(), nullptr);
      } else {
        drawOptions_.bondLineWidth =
            getHighlightBondWidth(drawOptions_, bond->getIdx(), nullptr) / 4;
      }
    }
  }

  return std::make_pair(col1, col2);
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
    // over-ride with explicit colour from highlight_map if there is one
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

} // namespace RDKit