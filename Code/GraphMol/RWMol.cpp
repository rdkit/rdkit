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
#include "SubstanceGroup.h"

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
    initFromOther(other, false, -1);
  }
  return *this;
}

void RWMol::insertMol(const ROMol &other) {
  //std::vector<unsigned int> newAtomIds(other.getNumAtoms());
    unsigned int start_idx = getNumAtoms();
    for(auto *atom : other.atoms()) {
        addAtom(atom->copy(), false, true);
    }
    
  for (unsigned int ati = 0; ati < other.getNumAtoms(); ++ati) {
    Atom *newAt = getAtomWithIdx(ati+start_idx);
    // take care of atom-numbering-dependent properties:
    INT_VECT nAtoms;
    if (newAt->getPropIfPresent(common_properties::_ringStereoAtoms, nAtoms)) {
      BOOST_FOREACH (int &val, nAtoms) {
        if (val < 0) {
            unsigned int old_idx = -val - 1;
            unsigned new_idx = old_idx + start_idx;
          val = -1 * (new_idx + 1);
        } else {
            unsigned int old_idx = val - 1;
            unsigned new_idx = old_idx + start_idx;
          val = new_idx + 1;
        }
      }
      newAt->setProp(common_properties::_ringStereoAtoms, nAtoms, true);
    }
  }
    for(auto *bond : other.bonds()) {
        Bond *bond_p = bond->copy();
        unsigned int idx1 = bond->getBeginAtomIdx() + start_idx;
        unsigned int idx2 = bond->getEndAtomIdx() + start_idx;
        bond_p->setOwningMol(this);
        bond_p->setBeginAtomIdx(idx1);
        bond_p->setEndAtomIdx(idx2);
        BOOST_FOREACH (int &v, bond_p->getStereoAtoms()) { v += start_idx; }
        addBond(bond_p, true);
    }

  // add atom to any conformers as well, if we have any
  if (other.getNumConformers() && !getNumConformers()) {
    for (auto cfi = other.beginConformers(); cfi != other.endConformers();
         ++cfi) {
      auto *nconf = new Conformer(getNumAtoms());
      nconf->set3D((*cfi)->is3D());
      nconf->setId((*cfi)->getId());
      for (unsigned int i = 0; i < other.getNumAtoms(); ++i) {
        nconf->setAtomPos(start_idx + i, (*cfi)->getAtomPos(i));
      }
      addConformer(nconf, false);
    }
  } else if (getNumConformers()) {
    if (other.getNumConformers() == getNumConformers()) {
      ConformerIterator cfi;
      ConstConformerIterator ocfi;
      for (cfi = beginConformers(), ocfi = other.beginConformers();
           cfi != endConformers(); ++cfi, ++ocfi) {
        (*cfi)->resize(getNumAtoms());
        for (unsigned int i = 0; i < other.getNumAtoms(); ++i) {
          (*cfi)->setAtomPos(start_idx + i, (*ocfi)->getAtomPos(i));
        }
      }
    } else {
      for (auto cfi = this->beginConformers(); cfi != this->endConformers();
           ++cfi) {
        (*cfi)->resize(getNumAtoms());
        for (unsigned int i=start_idx; i<getNumAtoms(); ++i) {
          (*cfi)->setAtomPos(i, RDGeom::Point3D(0.0, 0.0, 0.0));
        }
      }
    }
  }
}

unsigned int RWMol::addAtom(bool updateLabel) {
  auto *atom_p = new Atom();
  atom_p->setOwningMol(this);
    unsigned int which = static_cast<unsigned int>(_atoms.size());
    _atoms.push_back(atom_p);
    atom_p->setIdx(which);
  if (updateLabel) {
    clearAtomBookmark(ci_RIGHTMOST_ATOM);
    setAtomBookmark(atom_p, ci_RIGHTMOST_ATOM);
  }

  // add atom to any conformers as well, if we have any
  for (auto cfi = this->beginConformers(); cfi != this->endConformers();
       ++cfi) {
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
  // this function gets more expensive.
  //  Could we call init from other?
 Atom *orig = _atoms[idx];
 
 atom_p->_oatoms = orig->_oatoms;
 atom_p->_bonds = orig->_bonds;
 for(auto *atom : _atoms) {
     for(auto &oatom: atom->_oatoms) {
         if(oatom == orig) {
             oatom = atom_p;
         }
     }
 }
    
  //MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
  if (preserveProps) {
    const bool replaceExistingData = false;
    atom_p->updateProps(*orig, replaceExistingData);
  }
  removeSubstanceGroupsReferencingAtom(*this, idx);
  delete orig;
  _atoms[idx] = atom_p;

  // FIX: do something about bookmarks
};

void RWMol::replaceBond(unsigned int idx, Bond *bond_pin, bool preserveProps) {
  PRECONDITION(bond_pin, "bad bond passed to replaceBond");
  URANGE_CHECK(idx, getNumBonds());
  Bond *obond = _bonds[idx];
  Bond *bond_p = bond_pin->copy();
  _bonds[idx] = bond_p;
  bond_p->setOwningMol(this);
  bond_p->setIdx(idx);
  bond_p->setBeginAtomIdx(obond->getBeginAtomIdx());
  bond_p->setEndAtomIdx(obond->getEndAtomIdx());
    
  {
    auto &obonds = _atoms[obond->getBeginAtomIdx()]->_bonds;
    for(size_t idx = 0; idx < obonds.size(); ++idx) {
        if(obonds[idx] == obond) {
            obonds[idx] = bond_p;
            break;
        }
    }
  }
  {
    auto &obonds = _atoms[obond->getEndAtomIdx()]->_bonds;
     for(size_t idx = 0; idx < obonds.size(); ++idx) {
         if(obonds[idx] == obond) {
             obonds[idx] = bond_p;
             break;
         }
     }
   }
  if (preserveProps) {
    const bool replaceExistingData = false;
    bond_p->updateProps(*obond, replaceExistingData);
  }

    delete obond;
    // FIX: do something about bookmarks

  removeSubstanceGroupsReferencingBond(*this, idx);
};

Atom *RWMol::getActiveAtom() {
  if (hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
    return getAtomWithBookmark(ci_RIGHTMOST_ATOM);
  } else {
    return getLastAtom();
  }
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
  auto markI = marks->begin();
  while (markI != marks->end()) {
    const ATOM_PTR_LIST &atoms = markI->second;
    // we need to copy the iterator then increment it, because the
    // deletion we're going to do in clearAtomBookmark will invalidate
    // it.
    auto tmpI = markI;
    ++markI;
    if (std::find(atoms.begin(), atoms.end(), atom) != atoms.end()) {
      clearAtomBookmark(tmpI->first, atom);
    }
  }

  // remove bonds attached to the atom
  std::vector<std::pair<unsigned int, unsigned int>> nbrs;
  ADJ_ITER b1, b2;
  boost::tie(b1, b2) = getAtomNeighbors(atom);
  while (b1 != b2) {
    nbrs.emplace_back(atom->getIdx(), rdcast<unsigned int>((*b1)->getIdx()));
    ++b1;
  }

  for (auto &nbr : nbrs) {
    removeBond(nbr.first, nbr.second);
  }

  // loop over all atoms with higher indices and update their indices
  for (unsigned int i = idx + 1; i < getNumAtoms(); i++) {
    Atom *higher_index_atom = getAtomWithIdx(i);
    higher_index_atom->setIdx(i - 1);
  }

  // do the same with the coordinates in the conformations
  BOOST_FOREACH (CONFORMER_SPTR conf, d_confs) {
    RDGeom::POINT3D_VECT &positions = conf->getPositions();
    auto pi = positions.begin();
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
  for(auto *bond: _bonds) {
    unsigned int tmpIdx = bond->getBeginAtomIdx();
    if (tmpIdx > idx) {
      bond->setBeginAtomIdx(tmpIdx - 1);
    }
    tmpIdx = bond->getEndAtomIdx();
    if (tmpIdx > idx) {
      bond->setEndAtomIdx(tmpIdx - 1);
    }
    bond->setIdx(nBonds++);
    for (auto bsi = bond->getStereoAtoms().begin();
         bsi != bond->getStereoAtoms().end(); ++bsi) {
      if ((*bsi) == rdcast<int>(idx)) {
        bond->getStereoAtoms().clear();
        break;
      } else if ((*bsi) > rdcast<int>(idx)) {
        --(*bsi);
      }
    }
  }

  // Finally remove the atmo
  _atoms.erase(_atoms.begin() + idx);
    
  removeSubstanceGroupsReferencingAtom(*this, idx);

  // Remove any stereo group which includes the atom being deleted
  removeGroupsWithAtom(atom, d_stereo_groups);

  // clear computed properties and reset our ring info structure
  // they are pretty likely to be wrong now:
  clearComputedProps(true);

  delete atom;
}

unsigned int RWMol::addBond(unsigned int atomIdx1, unsigned int atomIdx2,
                            Bond::BondType bondType) {
  URANGE_CHECK(atomIdx1, getNumAtoms());
  URANGE_CHECK(atomIdx2, getNumAtoms());
  PRECONDITION(atomIdx1 != atomIdx2, "attempt to add self-bond");
  PRECONDITION(!_atoms[atomIdx1]->hasNbr(atomIdx2),
               "bond already exists");

  auto *b = new Bond(bondType);
  b->setOwningMol(this);
    Atom *a1 = getAtomWithIdx(atomIdx1);
    Atom *a2 = getAtomWithIdx(atomIdx2);
  if (bondType == Bond::AROMATIC) {
    b->setIsAromatic(1);
    //
    // assume that aromatic bonds connect aromatic atoms
    //   This is relevant for file formats like MOL, where there
    //   is no such thing as an aromatic atom, but bonds can be
    //   marked aromatic.
    //
    a1->setIsAromatic(1);
    a2->setIsAromatic(1);
  }
  //bool ok;
  //MolGraph::edge_descriptor which;
  //boost::tie(which, ok) = boost::add_edge(atomIdx1, atomIdx2, d_graph);
  //d_graph[which] = b;
  // unsigned int res = rdcast<unsigned int>(boost::num_edges(d_graph));
    auto idx = _bonds.size();
  b->setIdx(idx);
  b->setBeginAtomIdx(atomIdx1);
  b->setEndAtomIdx(atomIdx2);
    a1->_bonds.push_back(b);
    a1->_oatoms.push_back(a2);
    a2->_bonds.push_back(b);
    a2->_oatoms.push_back(a1);
  // if both atoms have a degree>1, reset our ring info structure,
  // because there's a non-trivial chance that it's now wrong.
  if (dp_ringInfo && dp_ringInfo->isInitialized() &&
      a1->nbrs().size() > 1 &&
      a2->nbrs().size() > 1 ) {
    dp_ringInfo->reset();
  }
    _bonds.push_back(b);
    return idx+1;  // doesn't return the bond idx, returns num bonds
}

unsigned int RWMol::addBond(Atom *atom1, Atom *atom2, Bond::BondType bondType) {
  PRECONDITION(atom1 && atom2, "NULL atom passed in");
  return addBond(atom1->getIdx(), atom2->getIdx(), bondType);
}

void RWMol::removeBond(unsigned int aid1, unsigned int aid2) {
  URANGE_CHECK(aid1, getNumAtoms());
  URANGE_CHECK(aid2, getNumAtoms());
  Bond *bnd = getBondBetweenAtoms(aid1, aid2);
  if (!bnd) {
    return;
  }
  unsigned int idx = bnd->getIdx();

  // remove any bookmarks which point to this bond:
  BOND_BOOKMARK_MAP *marks = getBondBookmarks();
  auto markI = marks->begin();
  while (markI != marks->end()) {
    BOND_PTR_LIST &bonds = markI->second;
    // we need to copy the iterator then increment it, because the
    // deletion we're going to do in clearBondBookmark will invalidate
    // it.
    auto tmpI = markI;
    ++markI;
    if (std::find(bonds.begin(), bonds.end(), bnd) != bonds.end()) {
      clearBondBookmark(tmpI->first, bnd);
    }
  }

  // loop over neighboring double bonds and remove their stereo atom
  //  information. This is definitely now invalid (was github issue 8)
    for(auto *obnd: _atoms[aid1]->bonds()) {
        if (std::find(obnd->getStereoAtoms().begin(), obnd->getStereoAtoms().end(),
                      aid2) != obnd->getStereoAtoms().end()) {
          obnd->getStereoAtoms().clear();
        }
    }
    
  for(auto *obnd: _atoms[aid2]->bonds()) {
      if (std::find(obnd->getStereoAtoms().begin(), obnd->getStereoAtoms().end(),
                    aid1) != obnd->getStereoAtoms().end()) {
        obnd->getStereoAtoms().clear();
      }
  }


  // reset our ring info structure, because it is pretty likely
  // to be wrong now:
  dp_ringInfo->reset();

  removeSubstanceGroupsReferencingBond(*this, idx);

  // loop over all bonds with higher indices and update their indices
  for(auto *bond : bonds()) {
    if (bond->getIdx() > idx) {
      bond->setIdx(bond->getIdx() - 1);
    }
  }
  bnd->setOwningMol(nullptr);
    _atoms[aid1]->removeNbr(_atoms[aid2]);
    _atoms[aid2]->removeNbr(_atoms[aid1]);
    _bonds.erase(_bonds.begin() + bnd->getIdx());
  delete bnd;
}

Bond *RWMol::createPartialBond(unsigned int atomIdx1, Bond::BondType bondType) {
  URANGE_CHECK(atomIdx1, getNumAtoms());

  auto *b = new Bond(bondType);
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

}  // namespace RDKit
