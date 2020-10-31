//
//  Copyright (C) 2003-2015 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <iostream>
#include <boost/foreach.hpp>

// our stuff
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "ROMol.h"
#include "Atom.h"
#include "QueryAtom.h"
#include "Bond.h"
#include "QueryBond.h"
#include "MolPickler.h"
#include "Conformer.h"
#include "SubstanceGroup.h"

namespace RDKit {
class QueryAtom;
class QueryBond;

const int ci_RIGHTMOST_ATOM = -0xBADBEEF;
const int ci_LEADING_BOND = -0xBADBEEF + 1;
const int ci_ATOM_HOLDER = -0xDEADD06;

void ROMol::destroy() {
  d_atomBookmarks.clear();
  d_bondBookmarks.clear();
    for(auto *atom : _atoms) {
        delete atom;
    }
  for(auto *bond : _bonds) {
      delete bond;
  }
    _atoms.clear();
    _bonds.clear();
  if (dp_ringInfo) {
    delete dp_ringInfo;
    dp_ringInfo = nullptr;
  }

  getSubstanceGroups(*this).clear();
  d_stereo_groups.clear();
};

ROMol::ROMol(const std::string &pickle) : RDProps() {
  initMol();
  MolPickler::molFromPickle(pickle, *this);
}

void ROMol::initFromOther(const ROMol &other, bool quickCopy, int confId) {
  if (this == &other) {
    return;
  }
  // std::cerr<<"    init from other: "<<this<<" "<<&other<<std::endl;
  // copy over the atoms
    
    // this is tricky?
    unsigned int idx = 0;
    for(auto *atom : other._atoms) {
        auto a = atom->copy();
        a->setOwningMol(this);
        a->setIdx(idx++);
        _atoms.push_back(a);
    }
    idx = 0;
    for(auto *bond : other._bonds) {
        auto b = bond->copy();
        b->setOwningMol(this);
        b->setIdx(idx++);
        _bonds.push_back(b);
    }
    
    // Copy over the internal graph structure
    {
        auto xi = _atoms.begin();
        auto yi = other._atoms.begin();
        while(xi != _atoms.end() && yi != other._atoms.end()) {
            for(auto *bond: (*yi)->bonds()) {
                (*xi)->bonds().push_back(_bonds[bond->getIdx()]);
            }
            for(auto *atom: (*yi)->nbrs()) {
                (*xi)->nbrs().push_back(_atoms[atom->getIdx()]);
            }
            ++xi;
            ++yi;
        }
        }

  // ring information
  if (dp_ringInfo) {
    delete dp_ringInfo;
  }
  if (other.dp_ringInfo) {
    dp_ringInfo = new RingInfo(*(other.dp_ringInfo));
  } else {
    dp_ringInfo = new RingInfo();
  }

  // enhanced stereochemical information
  d_stereo_groups.clear();
  for (auto &&otherGroup : other.d_stereo_groups) {
    std::vector<Atom *> atoms;
    for (auto &&otherAtom : otherGroup.getAtoms()) {
      atoms.push_back(getAtomWithIdx(otherAtom->getIdx()));
    }
    d_stereo_groups.emplace_back(otherGroup.getGroupType(), std::move(atoms));
  }

  if (!quickCopy) {
    // copy conformations
    for (auto ci = other.beginConformers(); ci != other.endConformers(); ++ci) {
      if (confId < 0 || rdcast<int>((*ci)->getId()) == confId) {
        auto *conf = new Conformer(*(*ci));
        this->addConformer(conf);
      }
    }

    // Copy sgroups
    for (const auto &sg : getSubstanceGroups(other)) {
      addSubstanceGroup(*this, sg);
    }

    d_props = other.d_props;

    // Bookmarks should be copied as well:
    BOOST_FOREACH (ATOM_BOOKMARK_MAP::value_type abmI, other.d_atomBookmarks) {
      BOOST_FOREACH (const Atom *aptr, abmI.second) {
        setAtomBookmark(getAtomWithIdx(aptr->getIdx()), abmI.first);
      }
    }
    BOOST_FOREACH (BOND_BOOKMARK_MAP::value_type bbmI, other.d_bondBookmarks) {
      BOOST_FOREACH (const Bond *bptr, bbmI.second) {
        setBondBookmark(getBondWithIdx(bptr->getIdx()), bbmI.first);
      }
    }
  } else {
    d_props.reset();
    STR_VECT computed;
    d_props.setVal(RDKit::detail::computedPropName, computed);
  }
  // std::cerr<<"---------    done init from other: "<<this<<"
  // "<<&other<<std::endl;
}

void ROMol::initMol() {
  d_props.reset();
  dp_ringInfo = new RingInfo();
  // ok every molecule contains a property entry called
  // RDKit::detail::computedPropName
  // which provides
  //  list of property keys that correspond to value that have been computed
  // this can used to blow out all computed properties while leaving the rest
  // along
  // initialize this list to an empty vector of strings
  STR_VECT computed;
  d_props.setVal(RDKit::detail::computedPropName, computed);
}

unsigned int ROMol::getAtomDegree(const Atom *at) const {
    return at->bonds().size();
};

unsigned int ROMol::getNumAtoms(bool onlyExplicit) const {
  int res = _atoms.size();
  if (!onlyExplicit) {
    // if we are interested in hydrogens as well add them up from
    // each
    for (ConstAtomIterator ai = beginAtoms(); ai != endAtoms(); ++ai) {
      res += (*ai)->getTotalNumHs();
    }
  }
  return res;
};
unsigned int ROMol::getNumHeavyAtoms() const {
  unsigned int res = 0;
  for (ConstAtomIterator ai = beginAtoms(); ai != endAtoms(); ++ai) {
    if ((*ai)->getAtomicNum() > 1) {
      ++res;
    }
  }
  return res;
};

// returns the first inserted atom with the given bookmark
Atom *ROMol::getAtomWithBookmark(int mark) {
  PRECONDITION(d_atomBookmarks.count(mark) != 0, "atom bookmark not found");
  PRECONDITION(d_atomBookmarks[mark].begin() != d_atomBookmarks[mark].end(),
               "atom bookmark not found");
  return *(d_atomBookmarks[mark].begin());
};

// returns all atoms with the given bookmark
ROMol::ATOM_PTR_LIST &ROMol::getAllAtomsWithBookmark(int mark) {
  PRECONDITION(d_atomBookmarks.count(mark) != 0, "atom bookmark not found");
  return d_atomBookmarks[mark];
};

// returns the unique atom with the given bookmark
Atom *ROMol::getUniqueAtomWithBookmark(int mark) {
  PRECONDITION(d_atomBookmarks.count(mark) == 1,
               "multiple atoms with same bookmark");
  return getAtomWithBookmark(mark);
}

// returns the first inserted bond with the given bookmark
Bond *ROMol::getBondWithBookmark(int mark) {
  PRECONDITION(d_bondBookmarks.count(mark) != 0, "bond bookmark not found");
  PRECONDITION(d_bondBookmarks[mark].begin() != d_bondBookmarks[mark].end(),
               "bond bookmark not found");
  return *(d_bondBookmarks[mark].begin());
};

// returns all bonds with the given bookmark
ROMol::BOND_PTR_LIST &ROMol::getAllBondsWithBookmark(int mark) {
  PRECONDITION(d_bondBookmarks.count(mark) != 0, "bond bookmark not found");
  return d_bondBookmarks[mark];
};

// returns the unique bond with the given bookmark
Bond *ROMol::getUniqueBondWithBookmark(int mark) {
  PRECONDITION(d_bondBookmarks.count(mark) == 1,
               "multiple bons with same bookmark");
  return getBondWithBookmark(mark);
}

void ROMol::clearAtomBookmark(const int mark) { d_atomBookmarks.erase(mark); }

void ROMol::clearAtomBookmark(int mark, const Atom *atom) {
  if (d_atomBookmarks.count(mark) != 0) {
    ATOM_PTR_LIST *entry = &d_atomBookmarks[mark];
    unsigned int tgtIdx = atom->getIdx();
    for (auto i = entry->begin(); i != entry->end(); ++i) {
      if ((*i)->getIdx() == tgtIdx) {
        entry->erase(i);
        break;
      }
    }
    if (entry->begin() == entry->end()) {
      d_atomBookmarks.erase(mark);
    }
  }
}

void ROMol::clearBondBookmark(int mark) { d_bondBookmarks.erase(mark); }
void ROMol::clearBondBookmark(int mark, const Bond *bond) {
  if (d_bondBookmarks.count(mark) != 0) {
    BOND_PTR_LIST *entry = &d_bondBookmarks[mark];
    unsigned int tgtIdx = bond->getIdx();
    for (auto i = entry->begin(); i != entry->end(); ++i) {
      if ((*i)->getIdx() == tgtIdx) {
        entry->erase(i);
        break;
      }
    }
    if (entry->begin() == entry->end()) {
      d_bondBookmarks.erase(mark);
    }
  }
}

unsigned int ROMol::getNumBonds(bool onlyHeavy) const {
  // By default return the bonds that connect only the heavy atoms
  // hydrogen connecting bonds are ignores
  int res = _bonds.size();
  if (!onlyHeavy) {
    // If we need hydrogen connecting bonds add them up
    for (ConstAtomIterator ai = beginAtoms(); ai != endAtoms(); ++ai) {
      res += (*ai)->getTotalNumHs();
    }
  }
  return res;
}

Bond *ROMol::getBondBetweenAtoms(int idx1, int idx2) {
  URANGE_CHECK(idx1, getNumAtoms());
  URANGE_CHECK(idx2, getNumAtoms());
  unsigned int uidx1 = rdcast<unsigned int>(idx1);
  unsigned int uidx2 = rdcast<unsigned int>(idx2);
  
    for (auto *bond : _bonds) {
        if ((bond->getBeginAtomIdx() == uidx1 && bond->getEndAtomIdx() == uidx2) ||
            (bond->getBeginAtomIdx() == uidx2 && bond->getEndAtomIdx() == uidx1)) {
            return bond;
        }
    }
    return nullptr;
}

const Bond *ROMol::getBondBetweenAtoms(int idx1,
                                       int idx2) const {
    
  URANGE_CHECK(idx1, getNumAtoms());
  URANGE_CHECK(idx2, getNumAtoms());
  unsigned int uidx1 = rdcast<unsigned int>(idx1);
  unsigned int uidx2 = rdcast<unsigned int>(idx2);
  
    for (auto *bond : _bonds) {
        if ((bond->getBeginAtomIdx() == uidx1 && bond->getEndAtomIdx() == uidx2) ||
            (bond->getBeginAtomIdx() == uidx2 && bond->getEndAtomIdx() == uidx1)) {
            return bond;
        }
    }
    return nullptr;
}

ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) {
  ADJ_ITER begin = at->nbrs().begin();
  ADJ_ITER end = at->nbrs().end();
  return std::make_pair(begin, end);
};
  
ROMol::CONST_ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
  CONST_ADJ_ITER begin = at->nbrs().begin();
  CONST_ADJ_ITER end = at->nbrs().end();
  return std::make_pair(begin, end);
};

ROMol::OBOND_ITER_PAIR ROMol::getAtomBonds(Atom const *at) const {
    return std::make_pair(at->bonds().begin(), at->bonds().end());
}


unsigned int ROMol::addAtom(Atom *atom_pin, bool updateLabel,
                            bool takeOwnership) {
  PRECONDITION(atom_pin, "null atom passed in");
  Atom *atom_p;
  if (!takeOwnership) {
    atom_p = atom_pin->copy();
  } else {
    atom_p = atom_pin;
  }

  atom_p->setOwningMol(this);
  auto which = _atoms.size();
  _atoms.push_back(atom_p);
  atom_p->setIdx(which);
  if (updateLabel) {
    replaceAtomBookmark(atom_p, ci_RIGHTMOST_ATOM);
  }
  for (auto cfi = this->beginConformers(); cfi != this->endConformers();
       ++cfi) {
    (*cfi)->setAtomPos(which, RDGeom::Point3D(0.0, 0.0, 0.0));
  }
  return rdcast<unsigned int>(which);
};

unsigned int ROMol::addBond(Bond *bond_pin, bool takeOwnership) {
  PRECONDITION(bond_pin, "null bond passed in");
  URANGE_CHECK(bond_pin->getBeginAtomIdx(), getNumAtoms());
  URANGE_CHECK(bond_pin->getEndAtomIdx(), getNumAtoms());
  PRECONDITION(bond_pin->getBeginAtomIdx() != bond_pin->getEndAtomIdx(),
               "attempt to add self-bond");

    PRECONDITION(!_atoms[bond_pin->getBeginAtomIdx()]->hasNbr(bond_pin->getEndAtomIdx()),
                 "bond already exists");

  Bond *bond_p;
  if (!takeOwnership) {
    bond_p = bond_pin->copy();
  } else {
    bond_p = bond_pin;
  }
    auto which = _bonds.size();
    bond_p->setIdx(which);
    bond_p->setOwningMol(this);
  // check to see if bond has been added???
    auto a1 = _atoms[bond_p->getBeginAtomIdx()];
    auto a2 = _atoms[bond_p->getEndAtomIdx()];
    a1->_bonds.push_back(bond_p);
    a2->_bonds.push_back(bond_p);
    a1->_oatoms.push_back(a2);
    a2->_oatoms.push_back(a1);
    _bonds.push_back(bond_p);
    
  return which+1;
}

void ROMol::setStereoGroups(std::vector<StereoGroup> stereo_groups) {
  d_stereo_groups = std::move(stereo_groups);
}

void ROMol::debugMol(std::ostream &str) const {
    str << "Atoms:" << std::endl;
    for(auto *atom : _atoms) {
       str << "\t" << *atom << " " << &atom->getOwningMol() << std::endl;
    }
    
    str << "Bonds:" << std::endl;
    for(auto *bond: _bonds) {
       str << "\t" << *bond << " " << &bond->getOwningMol() << std::endl;
    }
}

// --------------------------------------------
//
//  Iterators
//
// --------------------------------------------
ROMol::AtomIterator ROMol::beginAtoms() { return AtomIterator(this); }
ROMol::ConstAtomIterator ROMol::beginAtoms() const {
  return ConstAtomIterator(this);
}
ROMol::AtomIterator ROMol::endAtoms() {
  return AtomIterator(this, getNumAtoms());
}
ROMol::ConstAtomIterator ROMol::endAtoms() const {
  return ConstAtomIterator(this, getNumAtoms());
}

ROMol::AromaticAtomIterator ROMol::beginAromaticAtoms() {
  return AromaticAtomIterator(this);
}
ROMol::ConstAromaticAtomIterator ROMol::beginAromaticAtoms() const {
  return ConstAromaticAtomIterator(this);
}
ROMol::AromaticAtomIterator ROMol::endAromaticAtoms() {
  return AromaticAtomIterator(this, getNumAtoms());
}
ROMol::ConstAromaticAtomIterator ROMol::endAromaticAtoms() const {
  return ConstAromaticAtomIterator(this, getNumAtoms());
}

ROMol::HeteroatomIterator ROMol::beginHeteros() {
  return HeteroatomIterator(this);
}
ROMol::ConstHeteroatomIterator ROMol::beginHeteros() const {
  return ConstHeteroatomIterator(this);
}
ROMol::HeteroatomIterator ROMol::endHeteros() {
  return HeteroatomIterator(this, getNumAtoms());
}
ROMol::ConstHeteroatomIterator ROMol::endHeteros() const {
  return ConstHeteroatomIterator(this, getNumAtoms());
}

ROMol::QueryAtomIterator ROMol::beginQueryAtoms(QueryAtom const *what) {
  return QueryAtomIterator(this, what);
}
ROMol::ConstQueryAtomIterator ROMol::beginQueryAtoms(
    QueryAtom const *what) const {
  return ConstQueryAtomIterator(this, what);
}
ROMol::QueryAtomIterator ROMol::endQueryAtoms() {
  return QueryAtomIterator(this, getNumAtoms());
}
ROMol::ConstQueryAtomIterator ROMol::endQueryAtoms() const {
  return ConstQueryAtomIterator(this, getNumAtoms());
}
ROMol::MatchingAtomIterator ROMol::beginMatchingAtoms(bool (*what)(Atom *)) {
  return MatchingAtomIterator(this, what);
}
ROMol::ConstMatchingAtomIterator ROMol::beginMatchingAtoms(
    bool (*what)(const Atom *)) const {
  return ConstMatchingAtomIterator(this, what);
}
ROMol::MatchingAtomIterator ROMol::endMatchingAtoms() {
  return MatchingAtomIterator(this, getNumAtoms());
}
ROMol::ConstMatchingAtomIterator ROMol::endMatchingAtoms() const {
  return ConstMatchingAtomIterator(this, getNumAtoms());
}

ROMol::BondIterator ROMol::beginBonds() { return BondIterator(this); }
ROMol::ConstBondIterator ROMol::beginBonds() const {
  return ConstBondIterator(this);
}
ROMol::BondIterator ROMol::endBonds() {
  EDGE_ITER beg, end;
  boost::tie(beg, end) = getEdges();
  return BondIterator(this, end);
}
ROMol::ConstBondIterator ROMol::endBonds() const {
  CONST_EDGE_ITER beg, end;
  boost::tie(beg, end) = getEdges();
  return ConstBondIterator(this, end);
}
void ROMol::clearComputedProps(bool includeRings) const {
  // the SSSR information:
  if (includeRings) {
    this->dp_ringInfo->reset();
  }

  RDProps::clearComputedProps();

  for (auto atom : atoms()) {
    atom->clearComputedProps();
  }

  for (ConstBondIterator bondIt = this->beginBonds();
       bondIt != this->endBonds(); bondIt++) {
    (*bondIt)->clearComputedProps();
  }
}

void ROMol::updatePropertyCache(bool strict) {
    for(auto atom : _atoms) { atom->updatePropertyCache(strict); }
    for(auto bond: _bonds) { bond->updatePropertyCache(strict); }
}

bool ROMol::needsUpdatePropertyCache() const {
  for (ConstAtomIterator atomIt = this->beginAtoms();
       atomIt != this->endAtoms(); ++atomIt) {
    if ((*atomIt)->needsUpdatePropertyCache()) {
      return true;
    }
  }
  // there is no test for bonds yet since they do not obtain a valence property
  return false;
}

const Conformer &ROMol::getConformer(int id) const {
  // make sure we have more than one conformation
  if (d_confs.size() == 0) {
    throw ConformerException("No conformations available on the molecule");
  }

  if (id < 0) {
    return *(d_confs.front());
  }
  auto cid = (unsigned int)id;
  for (auto ci = this->beginConformers(); ci != this->endConformers(); ++ci) {
    if ((*ci)->getId() == cid) {
      return *(*ci);
    }
  }
  // we did not find a conformation with the specified ID
  std::string mesg = "Can't find conformation with ID: ";
  mesg += id;
  throw ConformerException(mesg);
}

Conformer &ROMol::getConformer(int id) {
  // make sure we have more than one conformation
  if (d_confs.size() == 0) {
    throw ConformerException("No conformations available on the molecule");
  }

  if (id < 0) {
    return *(d_confs.front());
  }
  auto cid = (unsigned int)id;
  for (auto ci = this->beginConformers(); ci != this->endConformers(); ++ci) {
    if ((*ci)->getId() == cid) {
      return *(*ci);
    }
  }
  // we did not find a conformation with the specified ID
  std::string mesg = "Can't find conformation with ID: ";
  mesg += id;
  throw ConformerException(mesg);
}

void ROMol::removeConformer(unsigned int id) {
  for (auto ci = d_confs.begin(); ci != d_confs.end(); ++ci) {
    if ((*ci)->getId() == id) {
      d_confs.erase(ci);
      return;
    }
  }
}

unsigned int ROMol::addConformer(Conformer *conf, bool assignId) {
  PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
               "Number of atom mismatch");
  if (assignId) {
    int maxId = -1;
    BOOST_FOREACH (CONFORMER_SPTR cptr, d_confs) {
      maxId = std::max((int)(cptr->getId()), maxId);
    }
    maxId++;
    conf->setId((unsigned int)maxId);
  }
  conf->setOwningMol(this);
  CONFORMER_SPTR nConf(conf);
  d_confs.push_back(nConf);
  return conf->getId();
}

}  // namespace RDKit
