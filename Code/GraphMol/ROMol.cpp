//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <iostream>

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

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit {
class QueryAtom;
class QueryBond;

const int ci_RIGHTMOST_ATOM = -0xBADBEEF;
const int ci_LEADING_BOND = -0xBADBEEF + 1;
const int ci_ATOM_HOLDER = -0xDEADD06;

void ROMol::destroy() {
  d_atomBookmarks.clear();
  d_bondBookmarks.clear();

  ATOM_ITER_PAIR atItP = boost::vertices(d_graph);
  while (atItP.first != atItP.second) {
    delete (d_graph)[*(atItP.first++)];
  }

  BOND_ITER_PAIR bondItP = boost::edges(d_graph);
  while (bondItP.first != bondItP.second) {
    delete (d_graph)[*(bondItP.first++)];
  }

  d_graph.clear();

  if (dp_ringInfo) {
    delete dp_ringInfo;
    dp_ringInfo = nullptr;
  }

  getSubstanceGroups(*this).clear();
  d_stereo_groups.clear();
};

ROMol::ROMol(const std::string &pickle) : RDProps() {
  initMol();
  numBonds = 0;
  MolPickler::molFromPickle(pickle, *this);
  numBonds = rdcast<unsigned int>(boost::num_edges(d_graph));
}

ROMol::ROMol(const std::string &pickle, unsigned int propertyFlags)
    : RDProps() {
  initMol();
  numBonds = 0;
  MolPickler::molFromPickle(pickle, *this, propertyFlags);
  numBonds = rdcast<unsigned int>(boost::num_edges(d_graph));
}

void ROMol::initFromOther(const ROMol &other, bool quickCopy, int confId) {
  if (this == &other) {
    return;
  }
  numBonds = 0;
  // std::cerr<<"    init from other: "<<this<<" "<<&other<<std::endl;
  // copy over the atoms
  const MolGraph &oGraph = other.d_graph;
  ROMol::ATOM_ITER_PAIR atItP = other.getVertices();
  while (atItP.first != atItP.second) {
    addAtom(oGraph[*atItP.first]->copy(), false, true);
    ++atItP.first;
  }

  // and the bonds:
  ROMol::BOND_ITER_PAIR bondItP = other.getEdges();
  while (bondItP.first != bondItP.second) {
    addBond(oGraph[*(bondItP.first++)]->copy(), true);
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

  if (other.dp_delAtoms) {
    dp_delAtoms.reset(new boost::dynamic_bitset<>(*other.dp_delAtoms));
  } else {
    dp_delAtoms.reset(nullptr);
  }
  if (other.dp_delBonds) {
    dp_delBonds.reset(new boost::dynamic_bitset<>(*other.dp_delBonds));
  } else {
    dp_delBonds.reset(nullptr);
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
    for (auto abmI : other.d_atomBookmarks) {
      for (const auto *aptr : abmI.second) {
        setAtomBookmark(getAtomWithIdx(aptr->getIdx()), abmI.first);
      }
    }
    for (auto bbmI : other.d_bondBookmarks) {
      for (const auto *bptr : bbmI.second) {
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
  PRECONDITION(at, "no atom");
  PRECONDITION(at->hasOwningMol(), "atom not associated with a molecule");
  PRECONDITION(&at->getOwningMol() == this,
               "atom not associated with this molecule");
  return rdcast<unsigned int>(boost::out_degree(at->getIdx(), d_graph));
};

unsigned int ROMol::getNumAtoms(bool onlyExplicit) const {
  int res = rdcast<int>(boost::num_vertices(d_graph));
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

Atom *ROMol::getAtomWithIdx(unsigned int idx) {
  PRECONDITION(getNumAtoms() > 0, "no atoms");
  URANGE_CHECK(idx, getNumAtoms());

  MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
  Atom *res = d_graph[vd];
  POSTCONDITION(res, "");
  return res;
}

const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
  PRECONDITION(getNumAtoms() > 0, "no atoms");
  URANGE_CHECK(idx, getNumAtoms());

  MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
  const Atom *res = d_graph[vd];

  POSTCONDITION(res, "");
  return res;
}

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
               "multiple bonds with same bookmark");
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
  int res = numBonds;
  if (!onlyHeavy) {
    // If we need hydrogen connecting bonds add them up
    for (ConstAtomIterator ai = beginAtoms(); ai != endAtoms(); ++ai) {
      res += (*ai)->getTotalNumHs();
    }
  }
  return res;
}

Bond *ROMol::getBondWithIdx(unsigned int idx) {
  PRECONDITION(getNumBonds() > 0, "no bonds");
  URANGE_CHECK(idx, getNumBonds());

  BOND_ITER_PAIR bIter = getEdges();
  for (unsigned int i = 0; i < idx; i++) {
    ++bIter.first;
  }
  Bond *res = d_graph[*(bIter.first)];

  POSTCONDITION(res != nullptr, "Invalid bond requested");
  return res;
}

const Bond *ROMol::getBondWithIdx(unsigned int idx) const {
  PRECONDITION(getNumBonds() > 0, "no bonds");
  URANGE_CHECK(idx, getNumBonds());

  BOND_ITER_PAIR bIter = getEdges();
  for (unsigned int i = 0; i < idx; i++) {
    ++bIter.first;
  }
  const Bond *res = d_graph[*(bIter.first)];

  POSTCONDITION(res != nullptr, "Invalid bond requested");
  return res;
}

Bond *ROMol::getBondBetweenAtoms(unsigned int idx1, unsigned int idx2) {
  URANGE_CHECK(idx1, getNumAtoms());
  URANGE_CHECK(idx2, getNumAtoms());
  Bond *res = nullptr;

  MolGraph::edge_descriptor edge;
  bool found;
  boost::tie(edge, found) = boost::edge(boost::vertex(idx1, d_graph),
                                        boost::vertex(idx2, d_graph), d_graph);
  if (found) {
    res = d_graph[edge];
  }
  return res;
}

const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
                                       unsigned int idx2) const {
  URANGE_CHECK(idx1, getNumAtoms());
  URANGE_CHECK(idx2, getNumAtoms());
  const Bond *res = nullptr;

  MolGraph::edge_descriptor edge;
  bool found;
  boost::tie(edge, found) = boost::edge(boost::vertex(idx1, d_graph),
                                        boost::vertex(idx2, d_graph), d_graph);
  if (found) {
    res = d_graph[edge];
  }
  return res;
}

ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
  PRECONDITION(at, "no atom");
  PRECONDITION(at->hasOwningMol(), "atom not associated with a molecule");
  PRECONDITION(&at->getOwningMol() == this,
               "atom not associated with this molecule");
  return boost::adjacent_vertices(at->getIdx(), d_graph);
};

ROMol::OBOND_ITER_PAIR ROMol::getAtomBonds(Atom const *at) const {
  PRECONDITION(at, "no atom");
  PRECONDITION(at->hasOwningMol(), "atom not associated with a molecule");
  PRECONDITION(&at->getOwningMol() == this,
               "atom not associated with this molecule");
  return boost::out_edges(at->getIdx(), d_graph);
}

ROMol::ATOM_ITER_PAIR ROMol::getVertices() { return boost::vertices(d_graph); }
ROMol::BOND_ITER_PAIR ROMol::getEdges() { return boost::edges(d_graph); }
ROMol::ATOM_ITER_PAIR ROMol::getVertices() const {
  return boost::vertices(d_graph);
}
ROMol::BOND_ITER_PAIR ROMol::getEdges() const { return boost::edges(d_graph); }

unsigned int ROMol::addAtom(Atom *atom_pin, bool updateLabel,
                            bool takeOwnership) {
  PRECONDITION(atom_pin, "null atom passed in");
  PRECONDITION(!takeOwnership || !atom_pin->hasOwningMol() ||
                   &atom_pin->getOwningMol() == this,
               "cannot take ownership of an atom which already has an owner");
  Atom *atom_p;
  if (!takeOwnership) {
    atom_p = atom_pin->copy();
  } else {
    atom_p = atom_pin;
  }

  atom_p->setOwningMol(this);
  MolGraph::vertex_descriptor which = boost::add_vertex(d_graph);
  d_graph[which] = atom_p;
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
  PRECONDITION(!takeOwnership || !bond_pin->hasOwningMol() ||
                   &bond_pin->getOwningMol() == this,
               "cannot take ownership of an bond which already has an owner");
  URANGE_CHECK(bond_pin->getBeginAtomIdx(), getNumAtoms());
  URANGE_CHECK(bond_pin->getEndAtomIdx(), getNumAtoms());
  PRECONDITION(bond_pin->getBeginAtomIdx() != bond_pin->getEndAtomIdx(),
               "attempt to add self-bond");
  PRECONDITION(!(boost::edge(bond_pin->getBeginAtomIdx(),
                             bond_pin->getEndAtomIdx(), d_graph)
                     .second),
               "bond already exists");

  Bond *bond_p;
  if (!takeOwnership) {
    bond_p = bond_pin->copy();
  } else {
    bond_p = bond_pin;
  }

  bond_p->setOwningMol(this);
  bool ok;
  MolGraph::edge_descriptor which;
  boost::tie(which, ok) = boost::add_edge(bond_p->getBeginAtomIdx(),
                                          bond_p->getEndAtomIdx(), d_graph);
  CHECK_INVARIANT(ok, "bond could not be added");
  d_graph[which] = bond_p;
  numBonds++;
  //  int res = rdcast<int>(boost::num_edges(d_graph));
  bond_p->setIdx(numBonds - 1);
  return numBonds;  // res;
}

void ROMol::setStereoGroups(std::vector<StereoGroup> stereo_groups) {
  d_stereo_groups = std::move(stereo_groups);
}

void ROMol::debugMol(std::ostream &str) const {
  ATOM_ITER_PAIR atItP = getVertices();
  BOND_ITER_PAIR bondItP = getEdges();

  str << "Atoms:" << std::endl;
  while (atItP.first != atItP.second) {
    str << "\t" << *d_graph[*(atItP.first++)] << std::endl;
  }

  str << "Bonds:" << std::endl;
  while (bondItP.first != bondItP.second) {
    str << "\t" << *d_graph[*(bondItP.first++)] << std::endl;
  }

  const auto &sgs = getSubstanceGroups(*this);
  if (!sgs.empty()) {
    str << "Substance Groups:" << std::endl;
    for (const auto &sg : sgs) {
      str << "\t" << sg << std::endl;
    }
  }

  const auto &stgs = getStereoGroups();
  if (!stgs.empty()) {
    unsigned idx = 0;
    str << "Stereo Groups:" << std::endl;
    for (const auto &stg : stgs) {
      str << "\t" << idx << ' ' << stg << std::endl;
      ++idx;
    }
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
  EDGE_ITER beg, end;
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
  for (AtomIterator atomIt = this->beginAtoms(); atomIt != this->endAtoms();
       ++atomIt) {
    (*atomIt)->updatePropertyCache(strict);
  }
  for (BondIterator bondIt = this->beginBonds(); bondIt != this->endBonds();
       ++bondIt) {
    (*bondIt)->updatePropertyCache(strict);
  }
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
  PRECONDITION(conf, "bad conformer");
  PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
               "Number of atom mismatch");
  if (assignId) {
    int maxId = -1;
    for (auto cptr : d_confs) {
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

#ifdef RDK_USE_BOOST_SERIALIZATION
template <class Archive>
void ROMol::save(Archive &ar, const unsigned int) const {
  std::string pkl;
  MolPickler::pickleMol(*this, pkl, PicklerOps::AllProps);
  ar << pkl;
}

template <class Archive>
void ROMol::load(Archive &ar, const unsigned int) {
  std::string pkl;
  ar >> pkl;

  delete dp_ringInfo;
  initMol();

  numBonds = 0;
  MolPickler::molFromPickle(pkl, *this, PicklerOps::AllProps);
  numBonds = rdcast<unsigned int>(boost::num_edges(d_graph));
}

template RDKIT_GRAPHMOL_EXPORT void ROMol::save<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, const unsigned int) const;
template RDKIT_GRAPHMOL_EXPORT void ROMol::load<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, const unsigned int);
#endif

}  // namespace RDKit
