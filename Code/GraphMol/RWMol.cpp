//
//  Copyright (C) 2003-2016 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/foreach.hpp>

// our stuff
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "RWMol.h"
#include "Atom.h"
#include "Bond.h"
#include "BondIterators.h"
#include "RingInfo.h"

namespace RDKit {
void RWMol::destroy() {
  ROMol::destroy();
  d_partialBonds.clear();
  d_partialBonds.resize(0);
};

RWMol &RWMol::operator=(const RWMol &other) {
  if (this != &other) {
    this->clear();
    d_partialBonds.clear();
    numBonds = 0;
    initFromOther(other, false, -1);
  }
  return *this;
}

void RWMol::insertMol(const ROMol &other) {
  std::vector<unsigned int> newAtomIds(other.getNumAtoms());
  VERTEX_ITER firstA, lastA;
  boost::tie(firstA, lastA) = boost::vertices(other.d_graph);
  while (firstA != lastA) {
    Atom *newAt = other.d_graph[*firstA]->copy();
    unsigned int idx = addAtom(newAt, false, true);
    newAtomIds[other.d_graph[*firstA]->getIdx()] = idx;
    ++firstA;
  }
  for (unsigned int ati = 0; ati < other.getNumAtoms(); ++ati) {
    Atom *newAt = getAtomWithIdx(newAtomIds[ati]);
    // take care of atom-numbering-dependent properties:
    INT_VECT nAtoms;
    if (newAt->getPropIfPresent(common_properties::_ringStereoAtoms, nAtoms)) {
      BOOST_FOREACH (int &val, nAtoms) {
        if (val < 0) {
          val = -1 * (newAtomIds[(-val - 1)] + 1);
        } else {
          val = newAtomIds[val - 1] + 1;
        }
      }
      newAt->setProp(common_properties::_ringStereoAtoms, nAtoms, true);
    }
  }

  EDGE_ITER firstB, lastB;
  boost::tie(firstB, lastB) = boost::edges(other.d_graph);
  while (firstB != lastB) {
    Bond *bond_p = other.d_graph[*firstB]->copy();
    unsigned int idx1, idx2;
    idx1 = newAtomIds[bond_p->getBeginAtomIdx()];
    idx2 = newAtomIds[bond_p->getEndAtomIdx()];
    bond_p->setOwningMol(this);
    bond_p->setBeginAtomIdx(idx1);
    bond_p->setEndAtomIdx(idx2);
    BOOST_FOREACH (int &v, bond_p->getStereoAtoms()) { v = newAtomIds[v]; }
    addBond(bond_p, true);
    ++firstB;
  }

  // add atom to any conformers as well, if we have any
  if (other.getNumConformers() && !getNumConformers()) {
    for (ConstConformerIterator cfi = other.beginConformers();
         cfi != other.endConformers(); ++cfi) {
      Conformer *nconf = new Conformer(getNumAtoms());
      nconf->set3D((*cfi)->is3D());
      nconf->setId((*cfi)->getId());
      for (unsigned int i = 0; i < newAtomIds.size(); ++i)
        nconf->setAtomPos(newAtomIds[i], (*cfi)->getAtomPos(i));
      addConformer(nconf, false);
    }
  } else if (getNumConformers()) {
    if (other.getNumConformers() == getNumConformers()) {
      ConformerIterator cfi;
      ConstConformerIterator ocfi;
      for (cfi = beginConformers(), ocfi = other.beginConformers();
           cfi != endConformers(); ++cfi, ++ocfi) {
        (*cfi)->resize(getNumAtoms());
        for (unsigned int i = 0; i < newAtomIds.size(); ++i)
          (*cfi)->setAtomPos(newAtomIds[i], (*ocfi)->getAtomPos(i));
      }
    } else {
      for (ConformerIterator cfi = this->beginConformers();
           cfi != this->endConformers(); ++cfi) {
        (*cfi)->resize(getNumAtoms());
        for (unsigned int i = 0; i < newAtomIds.size(); ++i)
          (*cfi)->setAtomPos(newAtomIds[i], RDGeom::Point3D(0.0, 0.0, 0.0));
      }
    }
  }
}

unsigned int RWMol::addAtom(bool updateLabel) {
  Atom *atom_p = new Atom();
  atom_p->setOwningMol(this);
  MolGraph::vertex_descriptor which = boost::add_vertex(d_graph);
  d_graph[which].reset(atom_p);
  atom_p->setIdx(which);
  if (updateLabel) {
    clearAtomBookmark(ci_RIGHTMOST_ATOM);
    setAtomBookmark(atom_p, ci_RIGHTMOST_ATOM);
  }
  // add atom to any conformers as well, if we have any
  for (ConformerIterator cfi = this->beginConformers();
       cfi != this->endConformers(); ++cfi) {
    (*cfi)->setAtomPos(which, RDGeom::Point3D(0.0, 0.0, 0.0));
  }
  return rdcast<unsigned int>(which);
}

void RWMol::replaceAtom(unsigned int idx, Atom *atom_pin, bool updateLabel,
                        bool preserveProps) {
  RDUNUSED_PARAM(updateLabel);
  PRECONDITION(atom_pin, "bad atom passed to replaceAtom");
  URANGE_CHECK(idx, getNumAtoms());
  Atom *atom_p = atom_pin->copy();
  atom_p->setOwningMol(this);
  atom_p->setIdx(idx);
  MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
  if (preserveProps) {
    const bool replaceExistingData = false;
    atom_p->updateProps(*d_graph[vd].get(), replaceExistingData);
  }
  d_graph[vd].reset(atom_p);
  // FIX: do something about bookmarks
};

void RWMol::replaceBond(unsigned int idx, Bond *bond_pin, bool preserveProps) {
  PRECONDITION(bond_pin, "bad bond passed to replaceBond");
  URANGE_CHECK(idx, getNumBonds());
  BOND_ITER_PAIR bIter = getEdges();
  for (unsigned int i = 0; i < idx; i++) ++bIter.first;
  BOND_SPTR obond = d_graph[*(bIter.first)];
  Bond *bond_p = bond_pin->copy();
  bond_p->setOwningMol(this);
  bond_p->setIdx(idx);
  bond_p->setBeginAtomIdx(obond->getBeginAtomIdx());
  bond_p->setEndAtomIdx(obond->getEndAtomIdx());
  if (preserveProps) {
    const bool replaceExistingData = false;
    bond_p->updateProps( *d_graph[*(bIter.first)].get(), replaceExistingData );
  }

  d_graph[*(bIter.first)].reset(bond_p);
  // FIX: do something about bookmarks
};

Atom *RWMol::getActiveAtom() {
  if (hasAtomBookmark(ci_RIGHTMOST_ATOM))
    return getAtomWithBookmark(ci_RIGHTMOST_ATOM);
  else
    return getLastAtom();
};

void RWMol::setActiveAtom(Atom *at) {
  clearAtomBookmark(ci_RIGHTMOST_ATOM);
  setAtomBookmark(at, ci_RIGHTMOST_ATOM);
};
void RWMol::setActiveAtom(unsigned int idx) {
  setActiveAtom(getAtomWithIdx(idx));
};

void RWMol::removeAtom(unsigned int idx) { removeAtom(getAtomWithIdx(idx)); }

void RWMol::removeAtom(Atom *atom) {
  PRECONDITION(atom, "NULL atom provided");
  PRECONDITION(static_cast<RWMol *>(&atom->getOwningMol()) == this,
               "atom not owned by this molecule");
  unsigned int idx = atom->getIdx();

  // remove any bookmarks which point to this atom:
  ATOM_BOOKMARK_MAP *marks = getAtomBookmarks();
  ATOM_BOOKMARK_MAP::iterator markI = marks->begin();
  while (markI != marks->end()) {
    const ATOM_PTR_LIST &atoms = markI->second;
    // we need to copy the iterator then increment it, because the
    // deletion we're going to do in clearAtomBookmark will invalidate
    // it.
    ATOM_BOOKMARK_MAP::iterator tmpI = markI;
    ++markI;
    if (std::find(atoms.begin(), atoms.end(), atom) != atoms.end()) {
      clearAtomBookmark(tmpI->first, atom);
    }
  }

  // remove bonds attached to the atom
  std::vector<std::pair<unsigned int, unsigned int> > nbrs;
  ADJ_ITER b1, b2;
  boost::tie(b1, b2) = getAtomNeighbors(atom);
  while (b1 != b2) {
    nbrs.push_back(std::make_pair(atom->getIdx(), rdcast<unsigned int>(*b1)));
    ++b1;
  }
  for (unsigned int i = 0; i < nbrs.size(); ++i) {
    removeBond(nbrs[i].first, nbrs[i].second);
  }

  // loop over all atoms with higher indices and update their indices
  for (unsigned int i = idx + 1; i < getNumAtoms(); i++) {
    Atom *atom = getAtomWithIdx(i);
    atom->setIdx(i - 1);
  }

  // do the same with the coordinates in the conformations
  BOOST_FOREACH (CONFORMER_SPTR conf, d_confs) {
    RDGeom::POINT3D_VECT &positions = conf->getPositions();
    RDGeom::POINT3D_VECT_I pi = positions.begin();
    for (unsigned int i = 0; i < getNumAtoms() - 1; i++) {
      ++pi;
      if (i >= idx) {
        positions[i] = positions[i + 1];
      }
    }
    positions.erase(pi);
  }
  // now deal with bonds:
  //   their end indices may need to be decremented and their
  //   indices will need to be handled
  unsigned int nBonds = 0;
  EDGE_ITER beg, end;
  boost::tie(beg, end) = getEdges();
  while (beg != end) {
    BOND_SPTR bond = d_graph[*beg++];
    unsigned int tmpIdx = bond->getBeginAtomIdx();
    if (tmpIdx > idx) bond->setBeginAtomIdx(tmpIdx - 1);
    tmpIdx = bond->getEndAtomIdx();
    if (tmpIdx > idx) bond->setEndAtomIdx(tmpIdx - 1);
    bond->setIdx(nBonds++);
    for (INT_VECT::iterator bsi = bond->getStereoAtoms().begin();
         bsi != bond->getStereoAtoms().end(); ++bsi) {
      if ((*bsi) == rdcast<int>(idx)) {
        bond->getStereoAtoms().clear();
        break;
      } else if ((*bsi) > rdcast<int>(idx)) {
        --(*bsi);
      }
    }
  }

  // clear computed properties and reset our ring info structure
  // they are pretty likely to be wrong now:
  clearComputedProps(true);

  atom->setOwningMol(NULL);

  // remove all connections to the atom:
  MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
  boost::clear_vertex(vd, d_graph);
  // finally remove the vertex itself
  boost::remove_vertex(vd, d_graph);
}

unsigned int RWMol::addBond(unsigned int atomIdx1, unsigned int atomIdx2,
                            Bond::BondType bondType) {
  URANGE_CHECK(atomIdx1, getNumAtoms());
  URANGE_CHECK(atomIdx2, getNumAtoms());
  PRECONDITION(atomIdx1 != atomIdx2, "attempt to add self-bond");
  PRECONDITION(!(boost::edge(atomIdx1, atomIdx2, d_graph).second),
               "bond already exists");

  Bond *b = new Bond(bondType);
  b->setOwningMol(this);
  if (bondType == Bond::AROMATIC) {
    b->setIsAromatic(1);
    //
    // assume that aromatic bonds connect aromatic atoms
    //   This is relevant for file formats like MOL, where there
    //   is no such thing as an aromatic atom, but bonds can be
    //   marked aromatic.
    //
    getAtomWithIdx(atomIdx1)->setIsAromatic(1);
    getAtomWithIdx(atomIdx2)->setIsAromatic(1);
  }
  bool ok;
  MolGraph::edge_descriptor which;
  boost::tie(which, ok) = boost::add_edge(atomIdx1, atomIdx2, d_graph);
  d_graph[which].reset(b);
  // unsigned int res = rdcast<unsigned int>(boost::num_edges(d_graph));
  ++numBonds;
  b->setIdx(numBonds - 1);
  b->setBeginAtomIdx(atomIdx1);
  b->setEndAtomIdx(atomIdx2);

  // if both atoms have a degree>1, reset our ring info structure,
  // because there's a non-trivial chance that it's now wrong.
  if (dp_ringInfo && dp_ringInfo->isInitialized() &&
      boost::out_degree(atomIdx1, d_graph) > 1 &&
      boost::out_degree(atomIdx2, d_graph) > 1) {
    dp_ringInfo->reset();
  }

  return numBonds;  // res;
}

unsigned int RWMol::addBond(Atom *atom1, Atom *atom2, Bond::BondType bondType) {
  PRECONDITION(atom1 && atom2, "NULL atom passed in");
  return addBond(atom1->getIdx(), atom2->getIdx(), bondType);
}

unsigned int RWMol::addBond(Atom::ATOM_SPTR atom1, Atom::ATOM_SPTR atom2,
                            Bond::BondType bondType) {
  return addBond(atom1->getIdx(), atom2->getIdx(), bondType);
}

void RWMol::removeBond(unsigned int aid1, unsigned int aid2) {
  URANGE_CHECK(aid1, getNumAtoms());
  URANGE_CHECK(aid2, getNumAtoms());
  Bond *bnd = getBondBetweenAtoms(aid1, aid2);
  if (!bnd) return;
  unsigned int idx = bnd->getIdx();

  // remove any bookmarks which point to this bond:
  BOND_BOOKMARK_MAP *marks = getBondBookmarks();
  BOND_BOOKMARK_MAP::iterator markI = marks->begin();
  while (markI != marks->end()) {
    BOND_PTR_LIST &bonds = markI->second;
    // we need to copy the iterator then increment it, because the
    // deletion we're going to do in clearBondBookmark will invalidate
    // it.
    BOND_BOOKMARK_MAP::iterator tmpI = markI;
    ++markI;
    if (std::find(bonds.begin(), bonds.end(), bnd) != bonds.end()) {
      clearBondBookmark(tmpI->first, bnd);
    }
  }

  // loop over neighboring double bonds and remove their stereo atom
  //  information. This is definitely now invalid (was github issue 8)
  ADJ_ITER a1, a2;
  boost::tie(a1, a2) = boost::adjacent_vertices(aid1, d_graph);
  while (a1 != a2) {
    unsigned int oIdx = rdcast<unsigned int>(*a1);
    ++a1;
    if (oIdx == aid2) continue;
    Bond *obnd = getBondBetweenAtoms(aid1, oIdx);
    if (!obnd) continue;
    obnd->getStereoAtoms().clear();
  }
  boost::tie(a1, a2) = boost::adjacent_vertices(aid2, d_graph);
  while (a1 != a2) {
    unsigned int oIdx = rdcast<unsigned int>(*a1);
    ++a1;
    if (oIdx == aid1) continue;
    Bond *obnd = getBondBetweenAtoms(aid2, oIdx);
    if (!obnd) continue;
    obnd->getStereoAtoms().clear();
  }

  // reset our ring info structure, because it is pretty likely
  // to be wrong now:
  dp_ringInfo->reset();

  // loop over all bonds with higher indices and update their indices
  ROMol::EDGE_ITER firstB, lastB;
  boost::tie(firstB, lastB) = this->getEdges();
  while (firstB != lastB) {
    BOND_SPTR bond = (*this)[*firstB];
    if (bond->getIdx() > idx) {
      bond->setIdx(bond->getIdx() - 1);
    }
    ++firstB;
  }
  bnd->setOwningMol(NULL);

  MolGraph::vertex_descriptor vd1 = boost::vertex(aid1, d_graph);
  MolGraph::vertex_descriptor vd2 = boost::vertex(aid2, d_graph);
  boost::remove_edge(vd1, vd2, d_graph);
  --numBonds;
}

Bond *RWMol::createPartialBond(unsigned int atomIdx1, Bond::BondType bondType) {
  URANGE_CHECK(atomIdx1, getNumAtoms());

  Bond *b = new Bond(bondType);
  b->setOwningMol(this);
  b->setBeginAtomIdx(atomIdx1);
  return b;
}

unsigned int RWMol::finishPartialBond(unsigned int atomIdx2, int bondBookmark,
                                      Bond::BondType bondType) {
  PRECONDITION(hasBondBookmark(bondBookmark), "no such partial bond");
  URANGE_CHECK(atomIdx2, getNumAtoms());

  Bond *bsp = getBondWithBookmark(bondBookmark);
  if (bondType == Bond::UNSPECIFIED) {
    bondType = bsp->getBondType();
  }
  return addBond(bsp->getBeginAtomIdx(), atomIdx2, bondType);
}

}  // end o' namespace
