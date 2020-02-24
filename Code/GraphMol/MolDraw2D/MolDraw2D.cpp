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
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <Geometry/point.h>

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

namespace {
void getBondHighlightsForAtoms(const ROMol &mol,
                               const vector<int> &highlight_atoms,
                               vector<int> &highlight_bonds) {
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
}  // namespace

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
    getBondHighlightsForAtoms(mol, *highlight_atoms, highlight_bonds);
  }
  drawMolecule(mol, legend, highlight_atoms, &highlight_bonds,
               highlight_atom_map, nullptr, highlight_radii, confId);
}

// ****************************************************************************
void MolDraw2D::doContinuousHighlighting(
    const ROMol &mol, const vector<int> *highlight_atoms,
    const vector<int> *highlight_bonds,
    const map<int, DrawColour> *highlight_atom_map,
    const map<int, DrawColour> *highlight_bond_map,
    const std::map<int, double> *highlight_radii) {
  PRECONDITION(activeMolIdx_ >= 0, "bad active mol");
  int orig_lw = lineWidth();
  int tgt_lw = lineWidth() * 8;
  // try to scale lw to reflect the overall scaling:
  tgt_lw = max(
      orig_lw * 2,
      min(tgt_lw,
          (int)(scale_ / 25. * tgt_lw)));  // the 25 here is extremely empirical
  bool orig_fp = fillPolys();
  if (highlight_bonds) {
    for(auto this_at: mol.atoms()) {
      int this_idx = this_at->getIdx();
      for(const auto &nbri: make_iterator_range(mol.getAtomBonds(this_at))) {
        const Bond *bond = mol[nbri];
        int nbr_idx = bond->getOtherAtomIdx(this_idx);
        if (nbr_idx < static_cast<int>(at_cds_[activeMolIdx_].size()) &&
            nbr_idx > this_idx) {
          if (std::find(highlight_bonds->begin(), highlight_bonds->end(),
                        bond->getIdx()) != highlight_bonds->end()) {
            DrawColour col = drawOptions().highlightColour;
            if (highlight_bond_map &&
                highlight_bond_map->find(bond->getIdx()) !=
                    highlight_bond_map->end()) {
              col = highlight_bond_map->find(bond->getIdx())->second;
            }
            setLineWidth(tgt_lw);
            Point2D at1_cds = at_cds_[activeMolIdx_][this_idx];
            Point2D at2_cds = at_cds_[activeMolIdx_][nbr_idx];
            drawLine(at1_cds, at2_cds, col, col);
          }
        }
      }
    }
  }
  if (highlight_atoms) {
    for(auto this_at: mol.atoms()) {
      int this_idx = this_at->getIdx();
      if (std::find(highlight_atoms->begin(), highlight_atoms->end(),
                    this_idx) != highlight_atoms->end()) {
        DrawColour col = drawOptions().highlightColour;
        if (highlight_atom_map &&
            highlight_atom_map->find(this_idx) != highlight_atom_map->end()) {
          col = highlight_atom_map->find(this_idx)->second;
        }
        vector<DrawColour> cols(1, col);
        drawHighlightedAtom(this_idx, cols, highlight_radii);
      }
    }
  }
  setLineWidth(orig_lw);
  setFillPolys(orig_fp);
}

// ****************************************************************************
void MolDraw2D::drawMolecule(const ROMol &mol,
                             const vector<int> *highlight_atoms,
                             const vector<int> *highlight_bonds,
                             const map<int, DrawColour> *highlight_atom_map,
                             const map<int, DrawColour> *highlight_bond_map,
                             const std::map<int, double> *highlight_radii,
                             int confId) {

  int origWidth = curr_width_;
  unique_ptr<RWMol> rwmol = setupMoleculeDraw(mol, highlight_atoms,
                                              highlight_radii, confId);
  ROMol const &draw_mol = rwmol ? *(rwmol) : mol;
  if(!draw_mol.getNumConformers()) {
    // clearly, the molecule is in a sorry state.
    return;
  }

  if (drawOptions().continuousHighlight) {
    // if we're doing continuous highlighting, start by drawing the highlights
    doContinuousHighlighting(draw_mol, highlight_atoms, highlight_bonds,
                             highlight_atom_map, highlight_bond_map,
                             highlight_radii);
    // at this point we shouldn't be doing any more highlighting, so blow out
    // those variables:
    highlight_bonds = nullptr;
    highlight_atoms = nullptr;
  } else if (drawOptions().circleAtoms && highlight_atoms) {
    setFillPolys(drawOptions().fillHighlights);
    for(auto this_at: draw_mol.atoms()) {
      int this_idx = this_at->getIdx();
      if (std::find(highlight_atoms->begin(), highlight_atoms->end(),
                    this_idx) != highlight_atoms->end()) {
        if (highlight_atom_map &&
            highlight_atom_map->find(this_idx) != highlight_atom_map->end()) {
          setColour(highlight_atom_map->find(this_idx)->second);
        } else {
          setColour(drawOptions().highlightColour);
        }
        Point2D p1 = at_cds_[activeMolIdx_][this_idx];
        Point2D p2 = at_cds_[activeMolIdx_][this_idx];
        double radius = drawOptions().highlightRadius;
        if (highlight_radii &&
            highlight_radii->find(this_idx) != highlight_radii->end()) {
          radius = highlight_radii->find(this_idx)->second;
        }
        Point2D offset(radius, radius);
        p1 -= offset;
        p2 += offset;
        drawEllipse(p1, p2);
      }
    }
    setFillPolys(true);
  }

  drawBonds(draw_mol, highlight_atoms, highlight_atom_map,
            highlight_bonds, highlight_bond_map);

  vector<DrawColour> atom_colours;
  for(auto this_at: draw_mol.atoms()) {
    atom_colours.emplace_back(getColour(this_at->getIdx(), highlight_atoms,
                                        highlight_atom_map));
  }
  finishMoleculeDraw(draw_mol, atom_colours);
  curr_width_ = origWidth;

  // {
  //   Point2D p1(x_min_, y_min_), p2(x_min_ + x_range_, y_min_ + y_range_);
  //   setColour(DrawColour(0, 0, 0));
  //   setFillPolys(false);
  //   drawRect(p1, p2);
  // }
}

// ****************************************************************************
void MolDraw2D::drawMolecule(const ROMol &mol, const std::string &legend,
                             const vector<int> *highlight_atoms,
                             const vector<int> *highlight_bonds,
                             const map<int, DrawColour> *highlight_atom_map,
                             const map<int, DrawColour> *highlight_bond_map,
                             const std::map<int, double> *highlight_radii,
                             int confId) {
  drawMolecule(mol, highlight_atoms, highlight_bonds, highlight_atom_map,
               highlight_bond_map, highlight_radii, confId);
  drawLegend(legend);

}

// ****************************************************************************
void MolDraw2D::drawMoleculeWithHighlights(const ROMol &mol, const string &legend,
                                           const map<int, vector<DrawColour> > &highlight_atom_map,
                                           const map<int, vector<DrawColour> > &highlight_bond_map,
                                           const map<int, double> &highlight_radii,
                                           const map<int, int> &highlight_linewidth_multipliers,
                                           int confId) {

  int origWidth = curr_width_;
  vector<int> highlight_atoms;
  for(auto ha: highlight_atom_map) {
    highlight_atoms.emplace_back(ha.first);
  }
  unique_ptr<RWMol> rwmol = setupMoleculeDraw(mol, &highlight_atoms,
                                              &highlight_radii, confId);
  ROMol const &draw_mol = rwmol ? *(rwmol) : mol;
  if(!draw_mol.getNumConformers()) {
    // clearly, the molecule is in a sorry state.
    return;
  }

  for(auto ha: highlight_atom_map) {
    cout << "highlighting atom " << ha.first << " with " << ha.second.size()
         << " colours" << endl;
    drawHighlightedAtom(ha.first, ha.second, &highlight_radii);
  }

  // draw plain bonds - bond highlights will go either side
  drawBonds(draw_mol);

  vector<DrawColour> atom_colours;
  for(auto this_at: draw_mol.atoms()) {
    // passing nullptr means that we'll get a colour based on atomic number
    // only.  These colours are used for the atom labels, and it seems odd
    // to make the labels the same colour as the highlights as they won't be
    // visible and in any case we can't do multi-coloured text.
    atom_colours.emplace_back(getColour(this_at->getIdx(), nullptr, nullptr));
  }

  // this puts on atom labels and such
  finishMoleculeDraw(draw_mol, atom_colours);
  curr_width_ = origWidth;

  drawLegend(legend);

}

// ****************************************************************************
namespace {
void get2DCoordsMol(RWMol &mol, double &offset, double spacing, double &maxY,
                    double &minY, int confId, bool shiftAgents,
                    double coordScale) {
  try {
    MolOps::sanitizeMol(mol);
  } catch (const MolSanitizeException &) {
    mol.updatePropertyCache(false);
    try {
      MolOps::Kekulize(mol, false);  // kekulize, but keep the aromatic flags!
    } catch (const MolSanitizeException &) {
      // don't need to do anything
    }
    MolOps::setHybridization(mol);
  }

  const bool canonOrient = true;
  RDDepict::compute2DCoords(mol, nullptr, canonOrient);
  MolDraw2DUtils::prepareMolForDrawing(
      mol, false);  // don't kekulize, we just did that
  double minX = 1e8;
  double maxX = -1e8;
  double vShift = 0;
  if (shiftAgents) {
    vShift = 1.1 * maxY / 2;
  }

  Conformer &conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    RDGeom::Point3D &p = conf.getAtomPos(i);
    p *= coordScale;
    minX = std::min(minX, conf.getAtomPos(i).x);
  }
  offset += fabs(minX);
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    RDGeom::Point3D &p = conf.getAtomPos(i);
    p.y += vShift;
    if (!shiftAgents) {
      maxY = std::max(p.y, maxY);
      minY = std::min(p.y, minY);
    }
    p.x += offset;
    maxX = std::max(p.x, maxX);
  }
  offset = maxX + spacing;
}

// ****************************************************************************
void get2DCoordsForReaction(ChemicalReaction &rxn, const MolDrawOptions &opts,
                            Point2D &arrowBegin, Point2D &arrowEnd,
                            std::vector<double> &plusLocs, double spacing,
                            const std::vector<int> *confIds) {
  RDUNUSED_PARAM(opts);
  plusLocs.resize(0);
  double maxY = -1e8, minY = 1e8;
  double offset = 0.0;

  // reactants
  for (unsigned int midx = 0; midx < rxn.getNumReactantTemplates(); ++midx) {
    // add space for the "+" if required
    if (midx > 0) {
      plusLocs.push_back(offset);
      offset += spacing;
    }
    ROMOL_SPTR reactant = rxn.getReactants()[midx];
    int cid = -1;
    if (confIds) {
      cid = (*confIds)[midx];
    }
    get2DCoordsMol(*(RWMol *)reactant.get(), offset, spacing, maxY, minY, cid,
                   false, 1.0);
  }
  arrowBegin.x = offset;

  offset += spacing;

  double begAgentOffset = offset;

  // we need to do the products now so that we know the full y range.
  // these will have the wrong X coordinates, but we'll fix that later.
  offset = 0;
  for (unsigned int midx = 0; midx < rxn.getNumProductTemplates(); ++midx) {
    // add space for the "+" if required
    if (midx > 0) {
      plusLocs.push_back(offset);
      offset += spacing;
    }
    ROMOL_SPTR product = rxn.getProducts()[midx];
    int cid = -1;
    if (confIds) {
      cid = (*confIds)[rxn.getNumReactantTemplates() +
                       rxn.getNumAgentTemplates() + midx];
    }
    get2DCoordsMol(*(RWMol *)product.get(), offset, spacing, maxY, minY, cid,
                   false, 1.0);
  }

  offset = begAgentOffset;
  // agents
  for (unsigned int midx = 0; midx < rxn.getNumAgentTemplates(); ++midx) {
    ROMOL_SPTR agent = rxn.getAgents()[midx];
    int cid = -1;
    if (confIds) {
      cid = (*confIds)[rxn.getNumReactantTemplates() + midx];
    }
    get2DCoordsMol(*(RWMol *)agent.get(), offset, spacing, maxY, minY, cid,
                   true, 0.45);
  }
  if (rxn.getNumAgentTemplates()) {
    arrowEnd.x = offset;  //- spacing;
  } else {
    arrowEnd.x = offset + 3 * spacing;
  }
  offset = arrowEnd.x + 1.5 * spacing;

  // now translate the products over
  for (unsigned int midx = 0; midx < rxn.getNumProductTemplates(); ++midx) {
    ROMOL_SPTR product = rxn.getProducts()[midx];
    int cid = -1;
    if (confIds) {
      cid = (*confIds)[rxn.getNumReactantTemplates() +
                       rxn.getNumAgentTemplates() + midx];
    }
    Conformer &conf = product->getConformer(cid);
    for (unsigned int aidx = 0; aidx < product->getNumAtoms(); ++aidx) {
      conf.getAtomPos(aidx).x += offset;
    }
  }

  // fix the plus signs too
  unsigned int startP = 0;
  if (rxn.getNumReactantTemplates() > 1) {
    startP = rxn.getNumReactantTemplates() - 1;
  }
  for (unsigned int pidx = startP; pidx < plusLocs.size(); ++pidx) {
    plusLocs[pidx] += offset;
  }

  arrowBegin.y = arrowEnd.y = minY + (maxY - minY) / 2;
}

}  // namespace

// ****************************************************************************
void MolDraw2D::drawReaction(
    const ChemicalReaction &rxn, bool highlightByReactant,
    const std::vector<DrawColour> *highlightColorsReactants,
    const std::vector<int> *confIds) {
  ChemicalReaction nrxn(rxn);
  double spacing = 1.0;
  Point2D arrowBegin, arrowEnd;
  std::vector<double> plusLocs;
  get2DCoordsForReaction(nrxn, drawOptions(), arrowBegin, arrowEnd, plusLocs,
                         spacing, confIds);

  ROMol *tmol = ChemicalReactionToRxnMol(nrxn);
  MolOps::findSSSR(*tmol);

  if (needs_scale_ &&
      (!nrxn.getNumReactantTemplates() || !nrxn.getNumProductTemplates())) {
    // drawMolecule() will figure out the scaling so that the molecule
    // fits the drawing pane. In order to ensure that we have space for the
    // arrow, we need to figure out the scaling on our own.
    RWMol tmol2;
    tmol2.addAtom(new Atom(0), true, true);
    tmol2.addAtom(new Atom(0), true, true);
    tmol2.addConformer(new Conformer(2), true);
    tmol2.getConformer().getAtomPos(0) =
        RDGeom::Point3D(arrowBegin.x, arrowBegin.y, 0);
    tmol2.getConformer().getAtomPos(1) =
        RDGeom::Point3D(arrowEnd.x, arrowEnd.y, 0);

    tmol2.insertMol(*tmol);
    at_cds_.push_back(std::vector<Point2D>());
    atomic_nums_.push_back(std::vector<int>());
    atom_syms_.push_back(std::vector<std::pair<std::string, OrientType>>());
    activeMolIdx_++;
    extractAtomCoords(tmol2, 0, true);
    calculateScale();
    needs_scale_ = false;
    activeMolIdx_--;
    at_cds_.pop_back();
    atomic_nums_.pop_back();
    atom_syms_.pop_back();
  }

  std::vector<int> *atom_highlights = nullptr;
  std::map<int, DrawColour> *atom_highlight_colors = nullptr;
  std::vector<int> *bond_highlights = nullptr;
  std::map<int, DrawColour> *bond_highlight_colors = nullptr;
  if (highlightByReactant) {
    const std::vector<DrawColour> *colors =
        &drawOptions().highlightColourPalette;
    if (highlightColorsReactants) {
      colors = highlightColorsReactants;
    }
    std::vector<int> atomfragmap;
    MolOps::getMolFrags(*tmol, atomfragmap);

    atom_highlights = new std::vector<int>();
    atom_highlight_colors = new std::map<int, DrawColour>();
    bond_highlights = new std::vector<int>();
    bond_highlight_colors = new std::map<int, DrawColour>();
    std::map<int, int> atommap_fragmap;
    for (unsigned int aidx = 0; aidx < tmol->getNumAtoms(); ++aidx) {
      int atomRole = -1;
      Atom *atom = tmol->getAtomWithIdx(aidx);
      if (atom->getPropIfPresent("molRxnRole", atomRole) && atomRole == 1 &&
          atom->getAtomMapNum()) {
        atommap_fragmap[atom->getAtomMapNum()] = atomfragmap[aidx];
        atom_highlights->push_back(aidx);
        (*atom_highlight_colors)[aidx] =
            (*colors)[atomfragmap[aidx] % colors->size()];

        atom->setAtomMapNum(0);
        // add highlighted bonds to lower-numbered
        // (and thus already covered) neighbors
        for(const auto &nbri: make_iterator_range(tmol->getAtomNeighbors(atom))) {
          const Atom *nbr = (*tmol)[nbri];
          if (nbr->getIdx() < aidx &&
              atomfragmap[nbr->getIdx()] == atomfragmap[aidx]) {
            int bondIdx =
                tmol->getBondBetweenAtoms(aidx, nbr->getIdx())->getIdx();
            bond_highlights->push_back(bondIdx);
            (*bond_highlight_colors)[bondIdx] = (*atom_highlight_colors)[aidx];
          }
        }
      }
    }
    for (unsigned int aidx = 0; aidx < tmol->getNumAtoms(); ++aidx) {
      int atomRole = -1;
      Atom *atom = tmol->getAtomWithIdx(aidx);
      if (atom->getPropIfPresent("molRxnRole", atomRole) && atomRole == 2 &&
          atom->getAtomMapNum() &&
          atommap_fragmap.find(atom->getAtomMapNum()) !=
              atommap_fragmap.end()) {
        atom_highlights->push_back(aidx);
        (*atom_highlight_colors)[aidx] =
            (*colors)[atommap_fragmap[atom->getAtomMapNum()] % colors->size()];

        atom->setAtomMapNum(0);
        // add highlighted bonds to lower-numbered
        // (and thus already covered) neighbors
        for(const auto &nbri: make_iterator_range(tmol->getAtomNeighbors(atom))) {
          const Atom *nbr = (*tmol)[nbri];
          if (nbr->getIdx() < aidx && (*atom_highlight_colors)[nbr->getIdx()] ==
                                          (*atom_highlight_colors)[aidx]) {
            int bondIdx =
                tmol->getBondBetweenAtoms(aidx, nbr->getIdx())->getIdx();
            bond_highlights->push_back(bondIdx);
            (*bond_highlight_colors)[bondIdx] = (*atom_highlight_colors)[aidx];
          }
        }
      }
    }
  }

  drawMolecule(*tmol, "", atom_highlights, bond_highlights,
               atom_highlight_colors, bond_highlight_colors);

  delete tmol;
  delete atom_highlights;
  delete atom_highlight_colors;
  delete bond_highlights;
  delete bond_highlight_colors;

  double o_font_size = fontSize();
  setFontSize(2 * options_.legendFontSize / scale_);
  DrawColour odc = colour();
  setColour(options_.symbolColour);

  // now add the symbols
  for(double plusLoc: plusLocs) {
    Point2D loc(plusLoc, arrowBegin.y);
    drawString("+", loc);
  }

  // The arrow:
  drawArrow(arrowBegin, arrowEnd);

  setColour(odc);
  setFontSize(o_font_size);
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
  if (!mols.size()) {
    return;
  }

  std::vector<RWMol> tmols;
  tmols.reserve(mols.size());
  Point2D minP, maxP;
  minP.x = minP.y = numeric_limits<double>::max();
  maxP.x = maxP.y = -numeric_limits<double>::max();
  for (unsigned int i = 0; i < mols.size(); ++i) {
    if (!mols[i]) {
      tmols.push_back(RWMol());
      continue;
    }
    tmols.push_back(*(mols[i]));
    if (drawOptions().prepareMolsBeforeDrawing) {
      MolDraw2DUtils::prepareMolForDrawing(tmols[i]);
    }
    Conformer &conf = tmols[i].getConformer(confIds ? (*confIds)[i] : -1);
    RDGeom::Point3D centroid = MolTransforms::computeCentroid(conf, false);
    for (unsigned int j = 0; j < conf.getNumAtoms(); ++j) {
      RDGeom::Point3D &pj = conf.getAtomPos(j);
      pj -= centroid;
      minP.x = std::min(minP.x, pj.x);
      minP.y = std::min(minP.y, pj.y);
      maxP.x = std::max(maxP.x, pj.x);
      maxP.y = std::max(maxP.y, pj.y);
    }
  }
  setScale(panelWidth(), panelHeight(), minP, maxP);
  int nCols = width() / panelWidth();
  int nRows = height() / panelHeight();
  for (unsigned int i = 0; i < mols.size(); ++i) {
    if (!mols[i]) {
      continue;
    }

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
    setOffset(col * panelWidth(), row * panelHeight());

    vector<int> *lhighlight_bonds = nullptr;
    if (highlight_bonds) {
      lhighlight_bonds = new std::vector<int>((*highlight_bonds)[i]);
    } else if (drawOptions().continuousHighlight && highlight_atoms) {
      lhighlight_bonds = new vector<int>();
      getBondHighlightsForAtoms(tmols[i], (*highlight_atoms)[i],
                                *lhighlight_bonds);
    };

    drawMolecule(tmols[i], legends ? (*legends)[i] : "",
                 highlight_atoms ? &(*highlight_atoms)[i] : nullptr,
                 lhighlight_bonds,
                 highlight_atom_maps ? &(*highlight_atom_maps)[i] : nullptr,
                 highlight_bond_maps ? &(*highlight_bond_maps)[i] : nullptr,
                 highlight_radii ? &(*highlight_radii)[i] : nullptr,
                 confIds ? (*confIds)[i] : -1);
    // save the drawn positions of the atoms on the molecule. This is the only
    // way that we can later add metadata
    auto tag = boost::str(boost::format("_atomdrawpos_%d") %
                          (confIds ? (*confIds)[i] : -1));
    for (unsigned int j = 0; j < mols[i]->getNumAtoms(); ++j) {
      auto pt = getDrawCoords(j);
      mols[i]->getAtomWithIdx(j)->setProp(tag, pt, true);
    }
    delete lhighlight_bonds;
  }
};

// ****************************************************************************
void MolDraw2D::highlightCloseContacts() {
  if (drawOptions().flagCloseContactsDist < 0) {
    return;
  }
  int tol =
      drawOptions().flagCloseContactsDist * drawOptions().flagCloseContactsDist;
  boost::dynamic_bitset<> flagged(at_cds_[activeMolIdx_].size());
  for (unsigned int i = 0; i < at_cds_[activeMolIdx_].size(); ++i) {
    if (flagged[i]) {
      continue;
    }
    Point2D ci = getDrawCoords(at_cds_[activeMolIdx_][i]);
    for (unsigned int j = i + 1; j < at_cds_[activeMolIdx_].size(); ++j) {
      if (flagged[j]) {
        continue;
      }
      Point2D cj = getDrawCoords(at_cds_[activeMolIdx_][j]);
      double d = (cj - ci).lengthSq();
      if (d <= tol) {
        flagged.set(i);
        flagged.set(j);
        break;
      }
    }
    if (flagged[i]) {
      Point2D p1 = at_cds_[activeMolIdx_][i];
      Point2D p2 = p1;
      Point2D offset(0.1, 0.1);
      p1 -= offset;
      p2 += offset;
      bool ofp = fillPolys();
      setFillPolys(false);
      DrawColour odc = colour();
      setColour(DrawColour(1, 0, 0));
      drawRect(p1, p2);
      setColour(odc);
      setFillPolys(ofp);
    }
  }
}

// ****************************************************************************
// transform a set of coords in the molecule's coordinate system
// to drawing system coordinates
Point2D MolDraw2D::getDrawCoords(const Point2D &mol_cds) const {
  double x = scale_ * (mol_cds.x - x_min_ + x_trans_);
  double y = scale_ * (mol_cds.y - y_min_ + y_trans_);
  // y is now the distance from the top of the image, we need to
  // invert that:
  x += x_offset_;
  y -= y_offset_;
  y = height() - y;
  return Point2D(x, y);
}

// ****************************************************************************
Point2D MolDraw2D::getDrawCoords(int at_num) const {
  PRECONDITION(activeMolIdx_ >= 0, "bad mol idx");
  return getDrawCoords(at_cds_[activeMolIdx_][at_num]);
}

// ****************************************************************************
Point2D MolDraw2D::getAtomCoords(const pair<int, int> &screen_cds) const {
  int x = int(double(screen_cds.first) / scale_ + x_min_ - x_trans_);
  int y =
      int(double(y_min_ - y_trans_ - (screen_cds.second - height()) / scale_));
  return Point2D(x, y);
}

Point2D MolDraw2D::getAtomCoords(const pair<double, double> &screen_cds) const {
  auto x = double(screen_cds.first / scale_ + x_min_ - x_trans_);
  auto y = double(y_min_ - y_trans_ - (screen_cds.second - height()) / scale_);
  return Point2D(x, y);
}

// ****************************************************************************
Point2D MolDraw2D::getAtomCoords(int at_num) const {
  PRECONDITION(activeMolIdx_ >= 0, "bad active mol");
  return at_cds_[activeMolIdx_][at_num];
}

// ****************************************************************************
void MolDraw2D::setFontSize(double new_size) {
  font_size_ = new_size;
}

// ****************************************************************************
void MolDraw2D::setScale(int width, int height, const Point2D &minv,
                         const Point2D &maxv) {
  PRECONDITION(width > 0, "bad width");
  PRECONDITION(height > 0, "bad height");
  needs_scale_ = false;

  x_min_ = minv.x;
  y_min_ = minv.y;
  double x_max = maxv.x;
  double y_max = maxv.y;

  x_range_ = x_max - x_min_;
  y_range_ = y_max - y_min_;

  if (x_range_ < 1.0e-4) {
    x_range_ = 1.0;
    x_min_ = -0.5;
  }
  if (y_range_ < 1.0e-4) {
    y_range_ = 1.0;
    y_min_ = -0.5;
  }

  // put a buffer round the drawing and calculate a final scale
  x_min_ -= drawOptions().padding * x_range_;
  x_range_ *= 1 + 2 * drawOptions().padding;
  y_min_ -= drawOptions().padding * y_range_;
  y_range_ *= 1 + 2 * drawOptions().padding;

  scale_ = std::min(double(width) / x_range_, double(height) / y_range_);
  double y_mid = y_min_ + 0.5 * y_range_;
  double x_mid = x_min_ + 0.5 * x_range_;
  Point2D mid = getDrawCoords(Point2D(x_mid, y_mid));
  // that used the offset, we need to remove that:
  mid.x -= x_offset_;
  mid.y += y_offset_;
  x_trans_ = (width / 2 - mid.x) / scale_;
  y_trans_ = (mid.y - height / 2) / scale_;
}

// ****************************************************************************
void MolDraw2D::calculateScale(int width, int height,
                               const std::vector<int> *highlight_atoms,
                               const std::map<int, double> *highlight_radii) {
  PRECONDITION(width > 0, "bad width");
  PRECONDITION(height > 0, "bad height");
  PRECONDITION(activeMolIdx_ >= 0, "bad active mol");

  // cout << "calculateScale" << endl;
  auto centre_picture = [&]() {
    double y_mid = y_min_ + 0.5 * y_range_;
    double x_mid = x_min_ + 0.5 * x_range_;
    Point2D mid = getDrawCoords(Point2D(x_mid, y_mid));
    // that used the offset, we need to remove that:
    mid.x -= x_offset_;
    mid.y += y_offset_;
    x_trans_ = (width / 2 - mid.x) / scale_;
    y_trans_ = (mid.y - height / 2) / scale_;
  };

  x_min_ = y_min_ = numeric_limits<double>::max();
  double x_max(-numeric_limits<double>::max()),
      y_max(-numeric_limits<double>::max());
  for (auto &pt : at_cds_[activeMolIdx_]) {
    x_min_ = std::min(pt.x, x_min_);
    y_min_ = std::min(pt.y, y_min_);
    x_max = std::max(pt.x, x_max);
    y_max = std::max(pt.y, y_max);
  }

  x_range_ = x_max - x_min_;
  y_range_ = y_max - y_min_;
  if (x_range_ < 1e-4) {
    x_range_ = 1.;
    x_min_ -= 0.5;
    x_max += 0.5;
  }
  if (y_range_ < 1e-4) {
    y_range_ = 1.;
    y_min_ -= 0.5;
    y_max += 0.5;
  }

  scale_ = std::min(double(width) / x_range_, double(height) / y_range_);

  // we may need to adjust the scale if there are atom symbols that go off
  // the edges, and we probably need to do it iteratively because
  // get_string_size uses the current value of scale_.
  // We also need to adjust for highlighted atoms if there are any.
  // And now we need to take account of strings with N/S orientation
  // as well.
  while (scale_ > 1e-4) {
    // cout << "scale iteration : " << scale_ << endl;
    // just in case the drawer uses its own font size that depends on scale.
    // MolDraw2DCairo is one such.
    setFontSize(font_size_);
    for (int i = 0, is = atom_syms_[activeMolIdx_].size(); i < is; ++i) {
      if (!atom_syms_[activeMolIdx_][i].first.empty()) {
        double atsym_width, atsym_height;
        getStringSize(atom_syms_[activeMolIdx_][i].first, atsym_width,
                      atsym_height);
        double this_x_min = at_cds_[activeMolIdx_][i].x;
        double this_x_max = at_cds_[activeMolIdx_][i].x;
        double this_y_min = at_cds_[activeMolIdx_][i].y - atsym_height / 2;
        double this_y_max = at_cds_[activeMolIdx_][i].y + atsym_height / 2;
        OrientType orient = atom_syms_[activeMolIdx_][i].second;
        if (orient == W) {
          this_x_min -= atsym_width;
        } else if (orient == E) {
          this_x_max += atsym_width;
        } else if(orient == N || orient == S) {
          vector<string> sym_bits = atomLabelToPieces(i);
          double height, width, cumm_height = 0.0;
          for (auto bit: sym_bits) {
            getStringSize(bit, width, height);
            cumm_height += height;
          }
          this_y_min = at_cds_[activeMolIdx_][i].y - cumm_height / 2;
          this_y_max = at_cds_[activeMolIdx_][i].y + cumm_height / 2;
        } else {
          this_x_max += atsym_width / 2;
          this_x_min -= atsym_width / 2;
        }
        x_max = std::max(x_max, this_x_max);
        x_min_ = std::min(x_min_, this_x_min);
        y_max = std::max(y_max, this_y_max);
        y_min_ = std::min(y_min_, this_y_min);
      }
      if (highlight_atoms) {
        if (highlight_atoms->end() !=
            find(highlight_atoms->begin(), highlight_atoms->end(), i)) {
          double radius = 0.4;
          if (highlight_radii &&
              highlight_radii->find(i) != highlight_radii->end()) {
            radius = highlight_radii->find(i)->second;
          }
          double this_x_min = at_cds_[activeMolIdx_][i].x - radius;
          double this_x_max = at_cds_[activeMolIdx_][i].x + radius;
          double this_y_min = at_cds_[activeMolIdx_][i].y - radius;
          double this_y_max = at_cds_[activeMolIdx_][i].y + radius;
          x_max = std::max(x_max, this_x_max);
          x_min_ = std::min(x_min_, this_x_min);
          y_max = std::max(y_max, this_y_max);
          y_min_ = std::min(y_min_, this_y_min);
        }
      }
    }
    double old_scale = scale_;
    x_range_ = x_max - x_min_;
    y_range_ = y_max - y_min_;
    if (x_range_ < 1e-4) {
      x_range_ = 1.;
    }
    if (y_range_ < 1e-4) {
      y_range_ = 1.;
    }
    scale_ = std::min(double(width) / x_range_, double(height) / y_range_);
    if (fabs(scale_ - old_scale) < 0.1) {
      break;
    }
  }

  // put a 5% buffer round the drawing and calculate a final scale
  x_min_ -= drawOptions().padding * x_range_;
  x_range_ *= 1 + 2 * drawOptions().padding;
  y_min_ -= drawOptions().padding * y_range_;
  y_range_ *= 1 + 2 * drawOptions().padding;

  if (x_range_ > 1e-4 || y_range_ > 1e-4) {
    scale_ = std::min(double(width) / x_range_, double(height) / y_range_);
    double fix_scale = scale_;
    // after all that, use the fixed scale unless it's too big, in which case
    // scale the drawing down to fit.
    // fixedScale takes precedence if both it and fixedBondLength are given.
    if(drawOptions().fixedBondLength > 0.0) {
      fix_scale = drawOptions().fixedBondLength;
    }
    if(drawOptions().fixedScale > 0.0) {
      fix_scale = double(width) * drawOptions().fixedScale;
    }
    if(scale_ > fix_scale) {
      scale_ = fix_scale;
    }
    centre_picture();
  } else {
    scale_ = 1;
    x_trans_ = 0.;
    y_trans_ = 0.;
  }

  // cout << "leaving calculateScale" << endl;
  // cout << "final scale : " << scale_ << endl;

}

// ****************************************************************************
// establishes whether to put string draw mode into super- or sub-script
// mode based on contents of instring from i onwards. Increments i
// appropriately and returns true or false depending on whether it did
// something or not.
bool MolDraw2D::setStringDrawMode(const string &instring,
                                  TextDrawType &draw_mode, int &i) const {
  string bit1 = instring.substr(i, 5);
  string bit2 = instring.substr(i, 6);

  // could be markup for super- or sub-script
  if (string("<sub>") == bit1) {
    draw_mode = TextDrawSubscript;
    i += 4;
    return true;
  } else if (string("<sup>") == bit1) {
    draw_mode = TextDrawSuperscript;
    i += 4;
    return true;
  } else if (string("</sub>") == bit2) {
    draw_mode = TextDrawNormal;
    i += 5;
    return true;
  } else if (string("</sup>") == bit2) {
    draw_mode = TextDrawNormal;
    i += 5;
    return true;
  }

  return false;
}

// ****************************************************************************
void MolDraw2D::drawLine(const Point2D &cds1, const Point2D &cds2,
                         const DrawColour &col1, const DrawColour &col2) {
  if (col1 == col2) {
    setColour(col1);
    drawLine(cds1, cds2);
  } else {
    Point2D mid = (cds1 + cds2);
    mid *= .5;

    setColour(col1);
    drawLine(cds1, mid);
    setColour(col2);
    drawLine(mid, cds2);
  }
}

// ****************************************************************************
// draws the string centred on cds
void MolDraw2D::drawString(const string &str, const Point2D &cds) {
  double string_width, string_height;
  getStringSize(str, string_width, string_height);

  // FIX: this shouldn't stay
  double M_width, M_height;
  getStringSize(std::string("M"), M_width, M_height);

  double draw_x = cds.x - string_width / 2.0;
  double draw_y = cds.y - string_height / 2.0;

  double full_font_size = fontSize();
  TextDrawType draw_mode = TextDrawNormal;
  string next_char(" ");

  for (int i = 0, is = str.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == str[i] && setStringDrawMode(str, draw_mode, i)) {
      continue;
    }

    char next_c = str[i];
    next_char[0] = next_c;
    double char_width, char_height;
    getStringSize(next_char, char_width, char_height);

    // these font sizes and positions work best for Qt, IMO. They may want
    // tweaking for a more general solution.
    if (TextDrawSubscript == draw_mode) {
      // y goes from top to bottom, so add for a subscript!
      setFontSize(0.75 * full_font_size);
      char_width *= 0.5;
      drawChar(next_c,
               //               getDrawCoords(Point2D(draw_x, draw_y + 0.5 *
               //               char_height)));
               getDrawCoords(Point2D(draw_x, draw_y - 0.5 * char_height)));
      setFontSize(full_font_size);
    } else if (TextDrawSuperscript == draw_mode) {
      setFontSize(0.75 * full_font_size);
      char_width *= 0.5;
      drawChar(next_c,
               //               getDrawCoords(Point2D(draw_x, draw_y - 0.25 *
               //               char_height)));
               getDrawCoords(Point2D(draw_x, draw_y + .5 * M_height)));
      setFontSize(full_font_size);
    } else {
      setFontSize(full_font_size);
      drawChar(next_c, getDrawCoords(Point2D(draw_x, draw_y)));
    }
    draw_x += char_width;
  }
}

// ****************************************************************************
void MolDraw2D::drawString(const std::string &str, const Point2D &cds,
                           AlignType align) {

  RDUNUSED_PARAM(align);
  drawString(str, cds);

}

// ****************************************************************************
void MolDraw2D::drawStrings(const std::vector<std::string> &labels,
                            const Point2D &cds, OrientType orient) {

  if(orient == W) {
    // stick the pieces together again backwards and draw as one so there
    // aren't ugly splits in the string.
    string new_lab;
    for(auto i = labels.rbegin(); i != labels.rend(); ++i) {
      new_lab += *i;
    }
    Point2D new_cds;
    alignString(new_lab, labels.front(), 1, cds, new_cds);
    drawString(new_lab, new_cds, END);
  } else if(orient == E) {
    // likewise, but forwards
    string new_lab;
    for(auto lab: labels) {
      new_lab += lab;
    }
    Point2D new_cds;
    alignString(new_lab, labels.front(), 0, cds, new_cds);
    drawString(new_lab, new_cds, START);
  } else {
    double x_scale = 0.0, y_scale = 0.0;
    if(orient == N) {
      x_scale = 0.0;
      y_scale = -1.0;
    } else if(orient == S) {
      x_scale = 0.0;
      y_scale = 1.0;
    }

    Point2D next_cds(cds);
    for (auto lab: labels) {
      drawString(lab, next_cds, MIDDLE);
      double width, height;
      getStringSize(lab, width, height);
      next_cds.x += x_scale * width;
      next_cds.y += y_scale * height;
    }
  }

}

// ****************************************************************************
void MolDraw2D::alignString(const string &str, const string &align_char,
                            int align, const Point2D &in_cds,
                            Point2D &out_cds) const {

  if(align != 0 && align != 1) {
    out_cds = in_cds;
    return;
  }

  double str_width, str_height;
  getStringSize(str, str_width, str_height);
  double ac_width, ac_height;
  getStringSize(align_char, ac_width, ac_height);
  // align == 0 is left align - first char to go at in_cds.
  double dir = align == 0 ? 1.0 : -1.0;
  out_cds.x = in_cds.x + dir * 0.5 * (str_width - ac_width);
  // assuming we centre the string one the draw coords.
  out_cds.y = in_cds.y;

}

// ****************************************************************************
DrawColour MolDraw2D::getColour(
    int atom_idx, const std::vector<int> *highlight_atoms,
    const std::map<int, DrawColour> *highlight_map) {
  PRECONDITION(activeMolIdx_ >= 0, "bad mol idx");
  PRECONDITION(atom_idx >= 0, "bad atom_idx");
  PRECONDITION(rdcast<int>(atomic_nums_[activeMolIdx_].size()) > atom_idx,
               "bad atom_idx");
  DrawColour retval =
      getColourByAtomicNum(atomic_nums_[activeMolIdx_][atom_idx]);

  // set contents of highlight_atoms to red
  if (!drawOptions().circleAtoms && !drawOptions().continuousHighlight) {
    if (highlight_atoms &&
        highlight_atoms->end() !=
            find(highlight_atoms->begin(), highlight_atoms->end(), atom_idx)) {
      retval = drawOptions().highlightColour;
    }
    // over-ride with explicit colour from highlight_map if there is one
    if (highlight_map) {
      auto p = highlight_map->find(atom_idx);
      if (p != highlight_map->end()) {
        retval = p->second;
      }
    }
  }
  return retval;
}

// ****************************************************************************
DrawColour MolDraw2D::getColourByAtomicNum(int atomic_num) {
  DrawColour res;
  if (drawOptions().atomColourPalette.find(atomic_num) !=
      drawOptions().atomColourPalette.end()) {
    res = drawOptions().atomColourPalette[atomic_num];
  } else if (atomic_num != -1 && drawOptions().atomColourPalette.find(-1) !=
                                     drawOptions().atomColourPalette.end()) {
    // if -1 is in the palette, we use that for undefined colors
    res = drawOptions().atomColourPalette[-1];
  } else {
    // if all else fails, default to black:
    res = DrawColour(0, 0, 0);
  }
  return res;
}

// ****************************************************************************
unique_ptr<RWMol> MolDraw2D::setupMoleculeDraw(const ROMol &mol,
                                               const vector<int> *highlight_atoms,
                                               const map<int, double> *highlight_radii,
                                               int confId) {

  at_cds_.push_back(std::vector<Point2D>());
  atomic_nums_.push_back(std::vector<int>());
  atom_syms_.push_back(std::vector<std::pair<std::string, OrientType>>());
  activeMolIdx_++;

  // prepareMolForDrawing needs a RWMol but don't copy the original mol
  // if we don't need to
  unique_ptr<RWMol> rwmol;
  if(drawOptions().prepareMolsBeforeDrawing || !mol.getNumConformers()) {
    rwmol.reset(new RWMol(mol));
    MolDraw2DUtils::prepareMolForDrawing(*rwmol, false);
  }
  ROMol const &draw_mol = rwmol ? *(rwmol) : mol;
  if(!draw_mol.getNumConformers()) {
    // clearly, the molecule is in a sorry state.
    return rwmol;
  }
  if (drawOptions().bondLineWidth >= 0) {
    curr_width_ = drawOptions().bondLineWidth;
  }
  if (!activeMolIdx_) {  // on the first pass we need to do some work
    if (drawOptions().clearBackground) {
      clearDrawing();
    }
    extractAtomCoords(draw_mol, confId, true);
    extractAtomSymbols(draw_mol);
    if (needs_scale_) {
      calculateScale(panel_width_, panel_height_, highlight_atoms,
                     highlight_radii);
      needs_scale_ = false;
    }
    // make sure the font doesn't end up too large (the constants are
    // empirical)
//    if (scale_ <= 40.) {
//      setFontSize(font_size_);
//    } else {
//      setFontSize(font_size_ * 30. / scale_);
//    }
//    setFontSize(font_size_ * 30. / scale_);
  } else {
    extractAtomCoords(draw_mol, confId, false);
    extractAtomSymbols(draw_mol);
  }

  if (drawOptions().includeAtomTags) {
    tagAtoms(draw_mol);
  }
  if (drawOptions().atomRegions.size()) {
    for(const std::vector<int> &region: drawOptions().atomRegions) {
      if (region.size() > 1) {
        Point2D minv = at_cds_[activeMolIdx_][region[0]];
        Point2D maxv = at_cds_[activeMolIdx_][region[0]];
        for(int idx: region) {
          const Point2D &pt = at_cds_[activeMolIdx_][idx];
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
        setColour(DrawColour(.8, .8, .8));
        // drawEllipse(minv,maxv);
        drawRect(minv, maxv);
      }
    }
  }

  return rwmol;

}

// ****************************************************************************
void MolDraw2D::drawBonds(const ROMol &draw_mol,
                          const vector<int> *highlight_atoms,
                          const map<int, DrawColour> *highlight_atom_map,
                          const vector<int> *highlight_bonds,
                          const map<int, DrawColour> *highlight_bond_map) {

  for(auto this_at: draw_mol.atoms()) {
    int this_idx = this_at->getIdx();
    for (const auto &nbri : make_iterator_range(draw_mol.getAtomBonds(this_at))) {
      const Bond *bond = draw_mol[nbri];
      int nbr_idx = bond->getOtherAtomIdx(this_idx);
      if (nbr_idx < static_cast<int>(at_cds_[activeMolIdx_].size()) &&
          nbr_idx > this_idx) {
        drawBond(draw_mol, bond, this_idx, nbr_idx, highlight_atoms,
                 highlight_atom_map, highlight_bonds, highlight_bond_map);
      }
    }
  }

}

// ****************************************************************************
void MolDraw2D::finishMoleculeDraw(const RDKit::ROMol &draw_mol,
                                   const vector<DrawColour> &atom_colours) {

  if (drawOptions().dummiesAreAttachments) {
    for(auto at1: draw_mol.atoms()) {
      if (at1->hasProp(common_properties::atomLabel) ||
          drawOptions().atomLabels.find(at1->getIdx()) !=
          drawOptions().atomLabels.end()) {
        // skip dummies that explicitly have a label provided
        continue;
      }
      if (at1->getAtomicNum() == 0 && at1->getDegree() == 1) {
        Point2D &at1_cds = at_cds_[activeMolIdx_][at1->getIdx()];
        const auto &iter_pair = draw_mol.getAtomNeighbors(at1);
        const Atom *at2 = draw_mol[*iter_pair.first];
        Point2D &at2_cds = at_cds_[activeMolIdx_][at2->getIdx()];
        drawAttachmentLine(at2_cds, at1_cds, DrawColour(.5, .5, .5));
      }
    }
  }

  for (int i = 0, is = atom_syms_[activeMolIdx_].size(); i < is; ++i) {
    if (!atom_syms_[activeMolIdx_][i].first.empty()) {
      drawAtomLabel(i, atom_colours[i]);
    }
  }

  if (drawOptions().flagCloseContactsDist >= 0) {
    highlightCloseContacts();
  }

}

// ****************************************************************************
void MolDraw2D::drawLegend(const string &legend) {

  if (!legend.empty()) {
    // the 0.94 is completely empirical and was brought over from Python
    Point2D loc =
        getAtomCoords(std::make_pair(panel_width_ / 2., 0.94 * panel_height_));
    double o_font_size = fontSize();
    setFontSize(options_.legendFontSize /
                scale_);  // set the font size to about 12 pixels high

    DrawColour odc = colour();
    setColour(options_.legendColour);
    drawString(legend, loc);
    setColour(odc);
    setFontSize(o_font_size);
  }

}

// ****************************************************************************
void MolDraw2D::drawHighlightedAtom(int atom_idx, const vector<DrawColour> &colours,
                                    const map<int, double> *highlight_radii) {

  Point2D p1 = at_cds_[activeMolIdx_][atom_idx];
  Point2D p2 = at_cds_[activeMolIdx_][atom_idx];
  double radius = drawOptions().highlightRadius;
  if (highlight_radii &&
      highlight_radii->find(atom_idx) != highlight_radii->end()) {
    radius = highlight_radii->find(atom_idx)->second;
  }

  setFillPolys(true);
  setLineWidth(1);
  if(colours.size() == 1) {
    setColour(colours.front());
    Point2D offset(radius, radius);
    p1 -= offset;
    p2 += offset;
    drawEllipse(p1, p2);
  } else {
    double arc_size = 360.0 / double(colours.size());
    double arc_start = 90.0;
    for(size_t i = 0; i < colours.size() ; ++i) {
      cout << "atom " << atom_idx << " arc " << i << " from " << arc_start
           << " to " << arc_start + arc_size << " colour "
           << colours[i].r << ", " << colours[i].g << ", "
           << colours[i].b << endl;
      setColour(colours[i]);
      drawArc(p1, radius, arc_start, arc_start + arc_size);
      arc_start += arc_size;
    }
  }

}

// ****************************************************************************
void MolDraw2D::extractAtomCoords(const ROMol &mol, int confId,
                                  bool updateBBox) {
  PRECONDITION(activeMolIdx_ >= 0, "no mol id");
  PRECONDITION(static_cast<int>(at_cds_.size()) > activeMolIdx_, "no space");
  PRECONDITION(static_cast<int>(atomic_nums_.size()) > activeMolIdx_,
               "no space");
  PRECONDITION(static_cast<int>(mol.getNumConformers()) > 0,
               "no coords");

  at_cds_[activeMolIdx_].clear();
  atomic_nums_[activeMolIdx_].clear();
  if (updateBBox) {
    bbox_[0].x = bbox_[0].y = numeric_limits<double>::max();
    bbox_[1].x = bbox_[1].y = -1 * numeric_limits<double>::max();
  }
  const RDGeom::POINT3D_VECT &locs = mol.getConformer(confId).getPositions();

  // the transformation rotates anti-clockwise, as is conventional, but
  // probably not what our user expects.
  double rot = - drawOptions().rotate * M_PI / 180.0;
  // assuming that if drawOptions().rotate is set to 0.0, rot will be
  // exactly 0.0 without worrying about floating point number dust.  Does
  // anyone know if this is true?  It's not the end of the world if not,
  // as it's just an extra largely pointless rotation.
  // Floating point numbers are like piles of sand; every time you move
  // them around, you lose a little sand and pick up a little dirt.
  // — Brian Kernighan and P.J. Plauger
  // Nothing brings fear to my heart more than a floating point number.
  // — Gerald Jay Sussman
  // Some developers, when encountering a problem, say: “I know, I’ll
  // use floating-point numbers !”   Now, they have 1.9999999997 problems.
  // — unknown
  double costh = rot == 0.0 ? 0.0 : cos(rot);
  double sinth = rot == 0.0 ? 0.0 : sin(rot);

  for(auto this_at: mol.atoms()) {
    int this_idx = this_at->getIdx();
    Point2D pt(locs[this_idx].x, locs[this_idx].y);
    if(rot != 0.0) {
      double new_x = pt.x * costh - pt.y * sinth;
      double new_y = pt.x * sinth + pt.y * costh;
      pt.x = new_x;
      pt.y = new_y;
    }
    at_cds_[activeMolIdx_].push_back(pt);

    if (updateBBox) {
      bbox_[0].x = std::min(bbox_[0].x, pt.x);
      bbox_[0].y = std::min(bbox_[0].y, pt.y);
      bbox_[1].x = std::max(bbox_[1].x, pt.x);
      bbox_[1].y = std::max(bbox_[1].y, pt.y);
    }
  }
}

// ****************************************************************************
void MolDraw2D::extractAtomSymbols(const ROMol &mol) {
  PRECONDITION(activeMolIdx_ >= 0, "no mol id");
  PRECONDITION(static_cast<int>(atom_syms_.size()) > activeMolIdx_, "no space");
  PRECONDITION(static_cast<int>(atomic_nums_.size()) > activeMolIdx_,
               "no space");

  for(auto at1: mol.atoms()) {
    Point2D &at1_cds = at_cds_[activeMolIdx_][at1->getIdx()];
    Point2D nbr_sum(0.0, 0.0);
    // cout << "Nbours for atom : " << at1->getIdx() << endl;
    for(const auto &nbri: make_iterator_range(mol.getAtomBonds(at1))) {
      const Bond *bond = mol[nbri];
      Point2D &at2_cds =
          at_cds_[activeMolIdx_][bond->getOtherAtomIdx(at1->getIdx())];
      nbr_sum += at2_cds - at1_cds;
    }
    atom_syms_[activeMolIdx_].push_back(
        getAtomSymbolAndOrientation(*at1, nbr_sum));
    atomic_nums_[activeMolIdx_].push_back(at1->getAtomicNum());
  }

}

// ****************************************************************************
void MolDraw2D::drawBond(const ROMol &mol, const Bond *bond, int at1_idx,
                         int at2_idx, const vector<int> *highlight_atoms,
                         const map<int, DrawColour> *highlight_atom_map,
                         const vector<int> *highlight_bonds,
                         const map<int, DrawColour> *highlight_bond_map) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(activeMolIdx_ >= 0, "bad mol idx");
  RDUNUSED_PARAM(highlight_atoms);
  RDUNUSED_PARAM(highlight_atom_map);
  static const DashPattern noDash;
  static const DashPattern dots = assign::list_of(2)(6);
  static const DashPattern dashes = assign::list_of(6)(6);
  static const DashPattern shortDashes = assign::list_of(2)(2);
  // the percent shorter that the extra bonds in a double bond are
  const double multipleBondTruncation = 0.15;

  const Atom *at1 = mol.getAtomWithIdx(at1_idx);
  const Atom *at2 = mol.getAtomWithIdx(at2_idx);
  Point2D at1_cds = at_cds_[activeMolIdx_][at1_idx];
  Point2D at2_cds = at_cds_[activeMolIdx_][at2_idx];

  double double_bond_offset = options_.multipleBondOffset;
  // mol files from, for example, Marvin use a bond length of 1 for just about
  // everything. When this is the case, the default multipleBondOffset is just
  // too much, so scale it back.
  if ((at1_cds - at2_cds).lengthSq() < 1.4) {
    double_bond_offset *= 0.6;
  }

  adjustBondEndForLabel(at1_idx, at2_cds, at1_cds);
  adjustBondEndForLabel(at2_idx, at1_cds, at2_cds);

  bool highlight_bond = false;
  if (highlight_bonds &&
      std::find(highlight_bonds->begin(), highlight_bonds->end(),
                bond->getIdx()) != highlight_bonds->end()) {
    highlight_bond = true;
  }

  DrawColour col1, col2;
  int orig_lw = lineWidth();
  if (!highlight_bond) {
    col1 = getColour(at1_idx);
    col2 = getColour(at2_idx);
  } else {
    if (highlight_bond_map &&
        highlight_bond_map->find(bond->getIdx()) != highlight_bond_map->end()) {
      col1 = col2 = highlight_bond_map->find(bond->getIdx())->second;
    } else {
      col1 = col2 = drawOptions().highlightColour;
    }
    if (drawOptions().continuousHighlight) {
      setLineWidth(orig_lw * 8);
    } else {
      setLineWidth(orig_lw * 2);
    }
  }

  Bond::BondType bt = bond->getBondType();
  bool isComplex = false;
  if (bond->hasQuery()) {
    std::string descr = bond->getQuery()->getDescription();
    if (bond->getQuery()->getNegation() || descr != "BondOrder") {
      isComplex = true;
    }
    if (isComplex) {
      setDash(dots);
      drawLine(at1_cds, at2_cds, col1, col2);
      setDash(noDash);
    } else {
      bt = static_cast<Bond::BondType>(
          static_cast<BOND_EQUALS_QUERY *>(bond->getQuery())->getVal());
    }
  }

  if (!isComplex) {
    // it's a double bond and one end is 1-connected, do two lines parallel
    // to the atom-atom line.
    if ((bt == Bond::DOUBLE) &&
        (1 == at1->getDegree() || 1 == at2->getDegree())) {
      Point2D perp = calcPerpendicular(at1_cds, at2_cds) * double_bond_offset;
      drawLine(at1_cds + perp, at2_cds + perp, col1, col2);
      drawLine(at1_cds - perp, at2_cds - perp, col1, col2);
    } else if (Bond::SINGLE == bt && (Bond::BEGINWEDGE == bond->getBondDir() ||
                                      Bond::BEGINDASH == bond->getBondDir())) {
      // std::cerr << "WEDGE: from " << at1->getIdx() << " | "
      //           << bond->getBeginAtomIdx() << "-" << bond->getEndAtomIdx()
      //           << std::endl;
      // swap the direction if at1 has does not have stereochem set
      // or if at2 does have stereochem set and the bond starts there
      if ((at1->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
           at1->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW) ||
          (at1->getIdx() != bond->getBeginAtomIdx() &&
           (at2->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
            at2->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW))) {
        // std::cerr << "  swap" << std::endl;
        swap(at1_cds, at2_cds);
        swap(col1, col2);
      }
      if (Bond::BEGINWEDGE == bond->getBondDir()) {
        drawWedgedBond(at1_cds, at2_cds, false, col1, col2);
      } else {
        drawWedgedBond(at1_cds, at2_cds, true, col1, col2);
      }
    } else if (Bond::SINGLE == bt && Bond::UNKNOWN == bond->getBondDir()) {
      // unspecified stereo
      drawWavyLine(at1_cds, at2_cds, col1, col2);
    } else if (Bond::DOUBLE == bt &&
               ((Bond::EITHERDOUBLE == bond->getBondDir()) ||
                (Bond::STEREOANY == bond->getStereo()))) {
      // crossed bond
      Point2D perp = calcPerpendicular(at1_cds, at2_cds);
      perp *= double_bond_offset;
      drawLine(at1_cds + perp, at2_cds - perp, col1, col2);
      drawLine(at1_cds - perp, at2_cds + perp, col1, col2);
    } else if (Bond::DATIVE == bt || Bond::DATIVEL == bt ||
               Bond::DATIVER == bt) {
      bool fps = fillPolys();
      setFillPolys(true);
      // draw an arrow for dative bonds
      bool asPolygon = true;
      double frac = 0.1;
      double angle = M_PI / 8;
      if (static_cast<unsigned int>(at1_idx) == bond->getBeginAtomIdx()) {
        drawArrow(at1_cds, at2_cds, asPolygon, frac, angle);
      } else {
        drawArrow(at2_cds, at1_cds, asPolygon, frac, angle);
      }
      setFillPolys(fps);
    } else if (Bond::ZERO == bt) {
      setDash(shortDashes);
      drawLine(at1_cds, at2_cds, col1, col2);
      setDash(noDash);
    } else {
      // in all other cases, we will definitely want to draw a line between
      // the two atoms
      drawLine(at1_cds, at2_cds, col1, col2);
      if (Bond::TRIPLE == bt) {
        // 2 further lines, a bit shorter and offset on the perpendicular
        double dbo = 2.0 * double_bond_offset;
        Point2D perp = calcPerpendicular(at1_cds, at2_cds);
        double end1_trunc =
            1 == at1->getDegree() ? 0.0 : multipleBondTruncation;
        double end2_trunc =
            1 == at2->getDegree() ? 0.0 : multipleBondTruncation;
        Point2D bv = at1_cds - at2_cds;
        Point2D p1 = at1_cds - (bv * end1_trunc) + perp * dbo;
        Point2D p2 = at2_cds + (bv * end2_trunc) + perp * dbo;
        drawLine(p1, p2, col1, col2);
        p1 = at1_cds - (bv * end1_trunc) - perp * dbo;
        p2 = at2_cds + (bv * end2_trunc) - perp * dbo;
        drawLine(p1, p2, col1, col2);
      } else if (Bond::DOUBLE == bt || Bond::AROMATIC == bt) {
        // all we have left now are double bonds in a ring or not in a ring
        // and multiply connected
        Point2D perp;
        if (mol.getRingInfo()->numBondRings(bond->getIdx())) {
          // in a ring, we need to draw the bond inside the ring
          perp = bondInsideRing(mol, bond, at1_cds, at2_cds);
        } else {
          perp = bondInsideDoubleBond(mol, bond);
        }
        double dbo = 2.0 * double_bond_offset;
        Point2D bv = at1_cds - at2_cds;
        Point2D p1 = at1_cds - bv * multipleBondTruncation + perp * dbo;
        Point2D p2 = at2_cds + bv * multipleBondTruncation + perp * dbo;
        if (bt == Bond::AROMATIC) {
          setDash(dashes);
        }
        drawLine(p1, p2, col1, col2);
        if (bt == Bond::AROMATIC) {
          setDash(noDash);
        }
      }
    }
  }
  if (highlight_bond) {
    setLineWidth(orig_lw);
  }
}

// ****************************************************************************
void MolDraw2D::drawWedgedBond(const Point2D &cds1, const Point2D &cds2,
                               bool draw_dashed, const DrawColour &col1,
                               const DrawColour &col2) {
  Point2D perp = calcPerpendicular(cds1, cds2);
  Point2D disp = perp * 0.15;
  // make sure the displacement isn't too large using the current scale factor
  // (part of github #985)
  // the constants are empirical to make sure that the wedge is visible, but
  // not absurdly large.
  if (scale_ > 40) {
    disp *= .6;
  }
  Point2D end1 = cds2 + disp;
  Point2D end2 = cds2 - disp;

  setColour(col1);
  if (draw_dashed) {
    unsigned int nDashes = 10;
    // empirical cutoff to make sure we don't have too many dashes in the
    // wedge:
    if ((cds1 - cds2).lengthSq() < 1.0) {
      nDashes /= 2;
    }

    int orig_lw = lineWidth();
    int tgt_lw = 1;  // use the minimum line width
    setLineWidth(tgt_lw);

    Point2D e1 = end1 - cds1;
    Point2D e2 = end2 - cds1;
    for (unsigned int i = 1; i < nDashes + 1; ++i) {
      if ((nDashes / 2 + 1) == i) {
        setColour(col2);
      }
      Point2D e11 = cds1 + e1 * (rdcast<double>(i) / nDashes);
      Point2D e22 = cds1 + e2 * (rdcast<double>(i) / nDashes);
      drawLine(e11, e22);
    }
    setLineWidth(orig_lw);

  } else {
    if (col1 == col2) {
      drawTriangle(cds1, end1, end2);
    } else {
      Point2D e1 = end1 - cds1;
      Point2D e2 = end2 - cds1;
      Point2D mid1 = cds1 + e1 * 0.5;
      Point2D mid2 = cds1 + e2 * 0.5;
      drawTriangle(cds1, mid1, mid2);
      setColour(col2);
      drawTriangle(mid1, end2, end1);
      drawTriangle(mid1, mid2, end2);
    }
  }
}

// ****************************************************************************
void MolDraw2D::drawAtomLabel(int atom_num,
                              const std::vector<int> *highlight_atoms,
                              const std::map<int, DrawColour> *highlight_map) {

  drawAtomLabel(atom_num, getColour(atom_num, highlight_atoms, highlight_map));

}

// ****************************************************************************
void MolDraw2D::drawAtomLabel(int atom_num, const DrawColour &draw_colour) {

  setColour(draw_colour);
  OrientType orient = atom_syms_[activeMolIdx_][atom_num].second;
  vector<string> label_pieces = atomLabelToPieces(atom_num);
  drawStrings(label_pieces, at_cds_[activeMolIdx_][atom_num], orient);

}

// ****************************************************************************
vector<string> MolDraw2D::atomLabelToPieces(int atom_num) const {

  vector<string> label_pieces;
  const string &atsym = atom_syms_[activeMolIdx_][atom_num].first;
  // cout << "Splitting string " << atsym << endl;

  // if we have the mark-up <lit>XX</lit> the symbol is to be used
  // without modification
  if(atsym.substr(0, 5) == "<lit>") {
    string lit_sym = atsym.substr(5);
    size_t idx = lit_sym.find("</lit>");
    if(idx != string::npos) {
      lit_sym = lit_sym.substr(0, idx);
    }
    label_pieces.emplace_back(lit_sym);
    return label_pieces;
  }

  string next_piece;
  size_t i = 0;
  while(true) {
    if(i == atsym.length()) {
      if(!next_piece.empty()) {
        label_pieces.emplace_back(next_piece);
        break;
      }
    }
    if(atsym.substr(i, 2) == "<s" || atsym[i] == ':' || isupper(atsym[i])) {
      // save the old piece, start a new one
      if(!next_piece.empty()) {
        label_pieces.emplace_back(next_piece);
        next_piece.clear();
      }
    }
    next_piece += atsym[i++];
  }

  if(label_pieces.size() < 2) {
    return label_pieces;
  }

  // now some re-arrangement to make things look nicer
  // if there's an atom map (:nn) or a charge it needs to go after the
  // first atomic symbol which, because of <sup> might not be the first
  // piece.
  for(size_t j = 0; j < label_pieces.size(); ++j) {
    if(label_pieces[j][0] == ':'
       || label_pieces[j].substr(0, 6) == "<sup>+"
       || label_pieces[j].substr(0, 6) == "<sup>-") {
      if (label_pieces[0].substr(0, 5) == "<sup>") {
        label_pieces[1] += label_pieces[j];
      } else {
        label_pieces[0] += label_pieces[j];
      }
      label_pieces[j].clear();
      break;
    }
  }
  label_pieces.erase(remove(label_pieces.begin(), label_pieces.end(), ""),
                     label_pieces.end());

  // if there's a <sub> piece, attach it to the one before.
  for(size_t j = 1; j < label_pieces.size(); ++j) {
    if(label_pieces[j].substr(0, 5) == "<sub>") {
      label_pieces[j-1] += label_pieces[j];
      label_pieces[j].clear();
      break;
    }
  }
  label_pieces.erase(remove(label_pieces.begin(), label_pieces.end(), ""),
                     label_pieces.end());

  // if there's a <sup> piece, attach it to the one after.
  if(label_pieces.size() > 1) {
    for (size_t j = 0; j < label_pieces.size() - 1; ++j) {
      if (label_pieces[j].substr(0, 5) == "<sup>") {
        label_pieces[j+1] = label_pieces[j] + label_pieces[j+1];
        label_pieces[j].clear();
        break;
      }
    }
  }
  label_pieces.erase(remove(label_pieces.begin(), label_pieces.end(), ""),
                     label_pieces.end());

  // cout << "Final pieces : ";
  // for(auto l: label_pieces) {
  //   cout << l << " X ";
  // }
  // cout << endl;

  return label_pieces;

}

// ****************************************************************************
// calculate normalised perpendicular to vector between two coords
Point2D MolDraw2D::calcPerpendicular(const Point2D &cds1, const Point2D &cds2) {
  double bv[2] = {cds1.x - cds2.x, cds1.y - cds2.y};
  double perp[2] = {-bv[1], bv[0]};
  double perp_len = sqrt(perp[0] * perp[0] + perp[1] * perp[1]);
  perp[0] /= perp_len;
  perp[1] /= perp_len;

  return Point2D(perp[0], perp[1]);
}

// ****************************************************************************
// cds1 and cds2 are 2 atoms in a ring.  Returns the perpendicular pointing
// into the ring
Point2D MolDraw2D::bondInsideRing(const ROMol &mol, const Bond *bond,
                                  const Point2D &cds1, const Point2D &cds2) {
  Atom *bgn_atom = bond->getBeginAtom();
  for(const auto &nbri2: make_iterator_range(mol.getAtomBonds(bgn_atom))) {
    const Bond *bond2 = mol[nbri2];
    // morphine (CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5)
    // showed a problem where one of the dashed bonds of the aromatic ring
    // was not in the aromatic ring with the rest, because the first bond2
    // we come to is in the aliphatic ring fused to the aromatic.  So check
    // the types are not 1 aromatic, the other not.
    if (bond2->getIdx() == bond->getIdx() ||
        !mol.getRingInfo()->numBondRings(bond2->getIdx()) ||
        (bond->getIsAromatic() && !bond2->getIsAromatic()) ||
        (!bond->getIsAromatic() && bond2->getIsAromatic())) {
      continue;
    }
    bool same_ring = false;
    for(const INT_VECT &ring: mol.getRingInfo()->bondRings()) {
      if (find(ring.begin(), ring.end(), bond->getIdx()) != ring.end() &&
          find(ring.begin(), ring.end(), bond2->getIdx()) != ring.end()) {
        same_ring = true;
        break;
      }
    }
    if (same_ring) {
      // bond and bond2 are in the same ring, so use their vectors to define
      // the sign of the perpendicular.
      int atom3 = bond2->getOtherAtomIdx(bond->getBeginAtomIdx());
      return calcInnerPerpendicular(cds1, cds2, at_cds_[activeMolIdx_][atom3]);
    }
  }

  return calcPerpendicular(cds1, cds2);
}

// ****************************************************************************
// cds1 and cds2 are 2 atoms in a chain double bond.  Returns the
// perpendicular pointing into the inside of the bond
Point2D MolDraw2D::bondInsideDoubleBond(const ROMol &mol, const Bond *bond) {
  // a chain double bond, were it looks nicer IMO if the 2nd line is inside
  // the angle of outgoing bond. Unless it's an allene, where nothing
  // looks great.
  const Atom *at1 = bond->getBeginAtom();
  const Atom *at2 = bond->getEndAtom();
  const Atom *bond_atom, *end_atom;
  if (at1->getDegree() > 1) {
    bond_atom = at1;
    end_atom = at2;
  } else {
    bond_atom = at2;
    end_atom = at1;
  }
  int at3 = -1;  // to stop the compiler whinging.
  for(const auto &nbri2: make_iterator_range(mol.getAtomBonds(bond_atom))) {
    const Bond *bond2 = mol[nbri2];
    if (bond != bond2) {
      at3 = bond2->getOtherAtomIdx(bond_atom->getIdx());
      break;
    }
  }

  return calcInnerPerpendicular(at_cds_[activeMolIdx_][end_atom->getIdx()],
                                at_cds_[activeMolIdx_][bond_atom->getIdx()],
                                at_cds_[activeMolIdx_][at3]);
}

// ****************************************************************************
// calculate normalised perpendicular to vector between two coords, such that
// it's inside the angle made between (1 and 2) and (2 and 3).
Point2D MolDraw2D::calcInnerPerpendicular(const Point2D &cds1,
                                          const Point2D &cds2,
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
// take the coords for atnum, with neighbour nbr_cds, and move cds out to
// accommodate the label associated with it.
void MolDraw2D::adjustBondEndForLabel(int atnum, const Point2D &nbr_cds,
                                      Point2D &cds) const {
  if (atom_syms_[activeMolIdx_][atnum].first.empty()) {
    return;
  }

  double label_width, label_height;
  // These days, labels only have the first character in the way.
  vector<string> label_bits = atomLabelToPieces(atnum);
  getStringSize(label_bits.front(), label_width, label_height);

  double additional_width = 0.0;
  double additional_height = 0.0;
  if (drawOptions().additionalAtomLabelPadding > 0.0) {
    double M_width, M_height;
    getStringSize("M", M_width, M_height);
    additional_width = M_width * drawOptions().additionalAtomLabelPadding;
    additional_height = M_height * drawOptions().additionalAtomLabelPadding;
  }
  double lw2 = label_width / 2.0 + additional_width;
  double lh2 = label_height / 2.0 + additional_height;

  double x_offset = 0.0, y_offset = 0.0;
  if (fabs(nbr_cds.y - cds.y) < 1.0e-5) {
    // if the bond is horizontal
    x_offset = lw2;
  } else {
    x_offset = fabs(lh2 * (nbr_cds.x - cds.x) / (nbr_cds.y - cds.y));
    if (x_offset >= lw2) {
      x_offset = lw2;
    }
  }
  if (nbr_cds.x < cds.x) {
    x_offset *= -1.0;
  }

  if (fabs(nbr_cds.x - cds.x) < 1.0e-5) {
    // if the bond is vertical
    y_offset = lh2;
  } else {
    y_offset = fabs(lw2 * (cds.y - nbr_cds.y) / (nbr_cds.x - cds.x));
    if (y_offset >= lh2) {
      y_offset = lh2;
    }
  }
  if (nbr_cds.y < cds.y) {
    y_offset *= -1.0;
  }

  cds.x += x_offset;
  cds.y += y_offset;
}

// ****************************************************************************
pair<string, MolDraw2D::OrientType> MolDraw2D::getAtomSymbolAndOrientation(
    const Atom &atom, const Point2D &nbr_sum) const {

  OrientType orient = getAtomOrientation(atom, nbr_sum);
  string symbol = getAtomSymbol(atom);

  // std::cout << "   res: " << symbol << " orient: " << orient
  //           << " nbr_sum:" << nbr_sum << std::endl << std::endl;

  return std::make_pair(symbol, orient);
}

// ****************************************************************************
string MolDraw2D::getAtomSymbol(const RDKit::Atom &atom) const {
  // adds XML-like annotation for super- and sub-script, in the same manner
  // as MolDrawing.py. My first thought was for a LaTeX-like system,
  // obviously...

  string symbol;
  bool literal_symbol = true;
  unsigned int iso = atom.getIsotope();
  if (drawOptions().atomLabels.find(atom.getIdx()) !=
      drawOptions().atomLabels.end()) {
    // specified labels are trump: no matter what else happens we will show
    // them.
    symbol = drawOptions().atomLabels.find(atom.getIdx())->second;
  } else if (atom.hasProp(common_properties::atomLabel)) {
    symbol = atom.getProp<std::string>(common_properties::atomLabel);
  } else if (drawOptions().dummiesAreAttachments && atom.getAtomicNum() == 0 &&
             atom.getDegree() == 1) {
    symbol = "";
  } else if (isComplexQuery(&atom)) {
    symbol = "?";
  } else if (drawOptions().atomLabelDeuteriumTritium &&
             atom.getAtomicNum() == 1 && (iso == 2 || iso == 3)) {
    symbol = ((iso == 2) ? "D" : "T");
    iso = 0;
  } else {
    literal_symbol = false;
    std::vector<std::string> preText, postText;
    int num_h = (atom.getAtomicNum() == 6 && atom.getDegree() > 0)
                ? 0
                : atom.getTotalNumHs();  // FIX: still not quite right
    if (num_h > 0 && !atom.hasQuery()) {
      // the H text comes after the atomic symbol
      std::string h = "H";
      if (num_h > 1) {
        // put the number as a subscript
        h += string("<sub>") + std::to_string(num_h) + string("</sub>");
      }
      postText.push_back(h);
    }

    if (0 != iso) {
      // isotope always comes before the symbol
      preText.push_back(std::string("<sup>") + std::to_string(iso) +
                        std::string("</sup>"));
    }

    if (0 != atom.getFormalCharge()) {
      // charge always comes post the symbol
      int ichg = atom.getFormalCharge();
      string sgn = ichg > 0 ? string("+") : string("-");
      ichg = abs(ichg);
      if (ichg > 1) {
        sgn += std::to_string(ichg);
      }
      // put the charge as a superscript
      postText.push_back(string("<sup>") + sgn + string("</sup>"));
    }

    if (atom.hasProp("molAtomMapNumber")) {
      // atom map always comes at the end
      string map_num = "";
      atom.getProp("molAtomMapNumber", map_num);
      postText.push_back(std::string(":") + map_num);
    }

    symbol = "";
    for(const std::string &se: preText) { symbol += se; }
    if (atom.getAtomicNum() != 6 || atom.getDegree() == 0 || preText.size() ||
        postText.size()) {
      symbol += atom.getSymbol();
    }
    for(const std::string &se: postText) { symbol += se; }
  }

  if(literal_symbol) {
    symbol = "<lit>" + symbol + "</lit>";
  }
  // cout << "Atom symbol " << atom.getIdx() << " : " << symbol << endl;
  return symbol;

}

// ****************************************************************************
MolDraw2D::OrientType MolDraw2D::getAtomOrientation(const RDKit::Atom &atom,
                                                    const Point2D &nbr_sum) const {

  // cout << "Atomic " << atom.getAtomicNum() << " degree : "
  //      << atom.getDegree() << " : " << atom.getTotalNumHs() << endl;
  // anything with a slope of more than 70 degrees is vertical. This way,
  // the NH in an indole is vertical as RDKit lays it out normally (72ish
  // degrees) but the 2 amino groups of c1ccccc1C1CCC(N)(N)CC1 are E and W
  // when they are drawn at the bottom of the molecule.
  static const double VERT_SLOPE = tan(70.0 * M_PI / 180.0);

  OrientType orient = C;
  if(atom.getDegree()) {
    double islope = 1000.0;
    if(fabs(nbr_sum.x) > 1.0e-4) {
      islope = nbr_sum.y / nbr_sum.x;
    }
    // cout << "islope : " << islope << " : " << atan(islope) * 180.0 / M_PI << endl;
    if(fabs(islope) <= VERT_SLOPE) {
      if(nbr_sum.x > 0.0) {
        orient = W;
      } else {
        orient = E;
      }
    } else {
      if(nbr_sum.y > 0.0) {
        orient = N;
      } else {
        orient = S;
      }
    }
    // cout << "interim orient : " << orient << endl;
    // atoms of single degree should always be either W or E, never N or S.  If
    // either of the latter, make it E if the slope is close to vertical,
    // otherwise have it either as required.
    if(atom.getDegree() == 1 && (orient == N || orient == S)) {
      if (fabs(islope) > VERT_SLOPE) {
        orient = E;
      } else {
        if (nbr_sum.x > 0.0) {
          orient = W;
        } else {
          orient = E;
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
      orient = W;
    } else {
      orient = E;
    }
  }

  // cout << "orient : " << orient << endl;
  return orient;

}

// ****************************************************************************
void MolDraw2D::drawTriangle(const Point2D &cds1, const Point2D &cds2,
                             const Point2D &cds3) {
  std::vector<Point2D> pts(3);
  pts[0] = cds1;
  pts[1] = cds2;
  pts[2] = cds3;
  drawPolygon(pts);
};

// ****************************************************************************
void MolDraw2D::drawArrow(const Point2D &arrowBegin, const Point2D &arrowEnd,
                          bool asPolygon, double frac, double angle) {
  Point2D delta = arrowBegin - arrowEnd;
  double cos_angle = std::cos(angle), sin_angle = std::sin(angle);

  Point2D p1 = arrowEnd;
  p1.x += frac * (delta.x * cos_angle + delta.y * sin_angle);
  p1.y += frac * (delta.y * cos_angle - delta.x * sin_angle);

  Point2D p2 = arrowEnd;
  p2.x += frac * (delta.x * cos_angle - delta.y * sin_angle);
  p2.y += frac * (delta.y * cos_angle + delta.x * sin_angle);

  drawLine(arrowBegin, arrowEnd);
  if (!asPolygon) {
    drawLine(arrowEnd, p1);
    drawLine(arrowEnd, p2);
  } else {
    std::vector<Point2D> pts = {p1, arrowEnd, p2};
    drawPolygon(pts);
  }
}

// ****************************************************************************
void MolDraw2D::drawEllipse(const Point2D &cds1, const Point2D &cds2) {

  std::vector<Point2D> pts;
  MolDraw2D_detail::arcPoints(cds1, cds2, pts, 0, 360);
  drawPolygon(pts);

}

// ****************************************************************************
void MolDraw2D::drawArc(const Point2D &cds1, double radius,
                        double ang1, double ang2) {

  std::vector<Point2D> pts;
  // 5 degree increments should be plenty, as the circles are probably
  // going to be small.
  int num_steps = (ang2 - ang1) / 5.0;
  double ang_incr = 5.0 * M_PI / 180.0;
  double start_ang_rads = ang2 * M_PI / 180.0;
  for(int i = 0; i <= num_steps; ++i) {
    double ang = start_ang_rads + double(i) * ang_incr;
    double x = cds1.x + radius * cos(ang);
    double y = cds1.y + radius * sin(ang);
    pts.emplace_back(Point2D(x, y));
  }

  if(drawOptions().fillHighlights) {
    // otherwise it draws an arc back to the pts.front() rather than filling
    // in the sector.
    pts.emplace_back(cds1);
  }
  drawPolygon(pts);

}

// ****************************************************************************
void MolDraw2D::drawRect(const Point2D &cds1, const Point2D &cds2) {
  std::vector<Point2D> pts(4);
  pts[0] = cds1;
  pts[1] = Point2D(cds1.x, cds2.y);
  pts[2] = cds2;
  pts[3] = Point2D(cds2.x, cds1.y);
  drawPolygon(pts);
}

void MolDraw2D::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                             const DrawColour &col1, const DrawColour &col2,
                             unsigned int nSegments, double vertOffset) {
  RDUNUSED_PARAM(nSegments);
  RDUNUSED_PARAM(vertOffset);
  drawLine(cds1, cds2, col1, col2);
}
// ****************************************************************************
//  we draw the line at cds2, perpendicular to the line cds1-cds2
void MolDraw2D::drawAttachmentLine(const Point2D &cds1, const Point2D &cds2,
                                   const DrawColour &col, double len,
                                   unsigned int nSegments) {
  Point2D perp = calcPerpendicular(cds1, cds2);
  Point2D p1 = Point2D(cds2.x - perp.x * len / 2, cds2.y - perp.y * len / 2);
  Point2D p2 = Point2D(cds2.x + perp.x * len / 2, cds2.y + perp.y * len / 2);
  drawWavyLine(p1, p2, col, col, nSegments);
}

}  // namespace RDKit
