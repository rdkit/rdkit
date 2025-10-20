//
//  Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

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
#include <GraphMol/rdmol_throw.h>
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

ROMol::ROMol() : d_ownedBySelf(true) {
  // d_ownedBySelf must be initialized before dp_mol, because dp_mol
  // compat data initialization depends on d_ownedBySelf.
  dp_mol = new RDKit::RDMol(this);
}
ROMol::ROMol(const ROMol &other, bool quickCopy , int confId) {
  // Self ownership must be set before dp_mol is created.
  d_ownedBySelf = true;
  dp_mol = new RDKit::RDMol(*other.dp_mol, quickCopy, confId, this);
}
ROMol::ROMol(const std::string &pickle) : ROMol() {
  MolPickler::molFromPickle(pickle, *this);
}
ROMol::ROMol(const std::string &pickle, unsigned int propertyFlags) : ROMol() {
  MolPickler::molFromPickle(pickle, *this, propertyFlags);
}
ROMol::ROMol(ROMol&& other) noexcept{
  dp_mol = other.dp_mol;
  d_ownedBySelf = true;
  if (other.d_ownedBySelf) {
    dp_mol->releaseCompatMolOwnership();
  }
  dp_mol->setROMolPointerCompat(this);
  dp_mol->setCompatPointersToSelf();
  // Moved from object must be usable per unit tests
  other.dp_mol = new RDKit::RDMol();
  other.d_ownedBySelf = true;
}

ROMol &ROMol::operator=(ROMol &&other) noexcept {
  if (this == &other) {
    return *this;
  }
  if (d_ownedBySelf && other.d_ownedBySelf) {
    delete dp_mol;
    dp_mol = other.dp_mol;
    // Moved from object must be usable per unit tests
    other.dp_mol = new RDKit::RDMol();
  } else {
    // In concept, something equivalent to this could happen:
    // RDMol mol;
    // mol.asROMol() = std::move(inputROMol);
    // in which case, this can't own the RDMol, so dp_mol must stay
    // the same and d_ownedBySelf must stay the same as before.
    // To avoid dp_mol deleting this upon assignment, temporarily release ownership of this.
    if (!d_ownedBySelf) {
      dp_mol->releaseCompatMolOwnership();
    }
    // Can't fully transfer ownership, so use copy assignment
    *dp_mol = *other.dp_mol;
  }
  dp_mol->setROMolPointerCompat(this);
  dp_mol->setCompatPointersToSelf();

  return *this;
}

ROMol::~ROMol() {
  if (d_ownedBySelf) {
    delete dp_mol;
  }
}

RingInfo* ROMol::getRingInfo() const {
  auto &ringInfo = dp_mol->getRingInfoCompat();
  return &ringInfo;
}

unsigned int ROMol::getNumAtoms() const { return dp_mol->getNumAtoms(); }
unsigned int ROMol::getNumAtoms(bool onlyExplicit) const {
  return dp_mol->getNumAtoms(onlyExplicit);
}
unsigned int ROMol::getNumHeavyAtoms() const {
  return dp_mol->getNumHeavyAtoms();
}
Atom* ROMol::getAtomWithIdx(unsigned int idx) {
  PRECONDITION(getNumAtoms() > 0, "no atoms");
  URANGE_CHECK(idx, getNumAtoms());
  return dp_mol->getAtomCompat(idx);
}
const Atom* ROMol::getAtomWithIdx(unsigned int idx) const {
  PRECONDITION(getNumAtoms() > 0, "no atoms");
  URANGE_CHECK(idx, getNumAtoms());
  return dp_mol->getAtomCompat(idx);
}
unsigned int ROMol::getAtomDegree(const Atom* at) const {
  PRECONDITION(at, "no atom");
  PRECONDITION(at->hasOwningMol(), "atom not associated with a molecule");
  PRECONDITION(&at->getOwningMol() == this,
               "atom not associated with this molecule");
  return dp_mol->getAtomDegree(at->getIdx());
}
unsigned int ROMol::getNumBonds(bool onlyHeavy) const {
  return dp_mol->getNumBonds(onlyHeavy);
}
Bond* ROMol::getBondWithIdx(unsigned int idx) {
  PRECONDITION(getNumBonds() > 0, "no bonds");
  URANGE_CHECK(idx, getNumBonds());
  return dp_mol->getBondCompat(idx);
}
const Bond* ROMol::getBondWithIdx(unsigned int idx) const {
  PRECONDITION(getNumBonds() > 0, "no bonds");
  URANGE_CHECK(idx, getNumBonds());
  return dp_mol->getBondCompat(idx);
}
Bond* ROMol::getBondBetweenAtoms(unsigned int idx1, unsigned int idx2) {
  URANGE_CHECK(idx1, getNumAtoms());
  URANGE_CHECK(idx2, getNumAtoms());
  const uint32_t bondIndex = dp_mol->getBondIndexBetweenAtoms(idx1, idx2);
  return (bondIndex == std::numeric_limits<std::uint32_t>::max())
             ? nullptr
             : getBondWithIdx(bondIndex);
}
const Bond* ROMol::getBondBetweenAtoms(unsigned int idx1, unsigned int idx2) const {
  URANGE_CHECK(idx1, getNumAtoms());
  URANGE_CHECK(idx2, getNumAtoms());
  const uint32_t bondIndex = dp_mol->getBondIndexBetweenAtoms(idx1, idx2);
  return (bondIndex == std::numeric_limits<std::uint32_t>::max())
             ? nullptr
             : getBondWithIdx(bondIndex);
}
ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const* at) const {
  PRECONDITION(at, "no atom");
  auto [beginNeighbors, endNeighbors] = dp_mol->getAtomNeighbors(at->getIdx());
  return ADJ_ITER_PAIR{ADJ_ITER{beginNeighbors}, ADJ_ITER{endNeighbors}};
}
ROMol::OBOND_ITER_PAIR ROMol::getAtomBonds(Atom const* at) const {
  PRECONDITION(at, "no atom");
  auto [beginBonds, endBonds] = dp_mol->getAtomBonds(at->getIdx());
  return OBOND_ITER_PAIR{OEDGE_ITER{beginBonds}, OEDGE_ITER{endBonds}};
}
ROMol::ATOM_ITER_PAIR ROMol::getVertices() const {
  return ATOM_ITER_PAIR{VERTEX_ITER{0}, VERTEX_ITER{getNumAtoms()}};
}
ROMol::BOND_ITER_PAIR ROMol::getEdges() const {
  return BOND_ITER_PAIR{EDGE_ITER{edge_descriptor{0}},
                        EDGE_ITER{edge_descriptor{getNumBonds()}}};
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

bool ROMol::hasQuery() const {
  for (auto atom : atoms()) {
    if (atom->hasQuery()) {
      return true;
    }
  }
  for (auto bond : bonds()) {
    if (bond->hasQuery()) {
      return true;
    }
  }
  return false;
}

ROMol::QueryAtomIterator ROMol::beginQueryAtoms(QueryAtom const* what) {
  return QueryAtomIterator(this, what);
}
ROMol::ConstQueryAtomIterator ROMol::beginQueryAtoms(
    QueryAtom const* what) const {
  return ConstQueryAtomIterator(this, what);
}
ROMol::QueryAtomIterator ROMol::endQueryAtoms() {
  return QueryAtomIterator(this, getNumAtoms());
}
ROMol::ConstQueryAtomIterator ROMol::endQueryAtoms() const {
  return ConstQueryAtomIterator(this, getNumAtoms());
}
ROMol::MatchingAtomIterator ROMol::beginMatchingAtoms(bool (*what)(Atom*)) {
  return MatchingAtomIterator(this, what);
}
ROMol::ConstMatchingAtomIterator ROMol::beginMatchingAtoms(
    bool (*what)(const Atom*)) const {
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
  auto [beg, end] = getEdges();
  return BondIterator(this, end);
}
ROMol::ConstBondIterator ROMol::endBonds() const {
  auto [beg, end] = getEdges();
  return ConstBondIterator(this, end);
}

unsigned int ROMol::getNumConformers() const {
  return dp_mol->getNumConformers();
}

void ROMol::clearComputedProps(bool includeRings) const {
  dp_mol->clearComputedProps(includeRings);
}

bool ROMol::needsUpdatePropertyCache() const {
  return dp_mol->needsUpdatePropertyCache();
}

void ROMol::updatePropertyCache(bool strict) {
  dp_mol->updatePropertyCache(strict);
}

void ROMol::clearPropertyCache() {
  dp_mol->clearPropertyCache();
}

const std::vector<StereoGroup> &ROMol::getStereoGroups() const {
  return dp_mol->getStereoGroupsCompat();
}

void ROMol::setStereoGroups(std::vector<StereoGroup> stereo_groups) {
  auto newGroups = std::make_unique<StereoGroups>();
  std::vector<StereoGroupType> &types = newGroups->stereoTypes;
  std::vector<uint32_t> &atomIndices = newGroups->atoms;
  std::vector<uint32_t> &bondIndices = newGroups->bonds;
  std::vector<uint32_t> &atomBegins = newGroups->atomBegins;
  std::vector<uint32_t> &bondBegins = newGroups->bondBegins;
  std::vector<uint32_t> &readIds = newGroups->readIds;
  std::vector<uint32_t> &writeIds = newGroups->writeIds;

  for (const auto &group : stereo_groups) {
    types.push_back(group.getGroupType());

    const auto &atoms = group.getAtoms();
    for (const auto *atom : atoms) {
      atomIndices.push_back(atom->getIdx());
    }
    const auto &bonds = group.getBonds();
    for (const auto *bond : bonds) {
      bondIndices.push_back(bond->getIdx());
    }

    atomBegins.push_back(atoms.size() + atomBegins.back());
    bondBegins.push_back(bonds.size() + bondBegins.back());
    readIds.push_back(group.getReadId());
    writeIds.push_back(group.getWriteId());
  }
  dp_mol->setStereoGroups(std::move(newGroups));
}

STR_VECT ROMol::getPropList(bool includePrivate, bool includeComputed) const {
  STR_VECT res = dp_mol->getPropList(includePrivate, includeComputed);
  if (includePrivate && includeComputed) {
    res.push_back(detail::computedPropName);
  }
  return res;
}

bool ROMol::hasProp(const std::string_view &key) const {
  PropToken token(key);
  if (token == detail::computedPropNameToken) {
    return true;
  }
  return dp_mol->hasProp(token);
}

void ROMol::clearProp(const std::string_view &key) const {
  dp_mol->clearMolPropIfPresent(PropToken(key));
}

void ROMol::updateProps(const RDKit::ROMol &source, bool preserveExisting) {
  if (&source == this) {
    return;
  }
  // Remove all properties if !preserveExisting
  if (!preserveExisting) {
    dp_mol->clearProps();
  }

  for (const RDMol::Property &prop : source.asRDMol().properties) {
    if (prop.scope() == RDMol::Scope::MOL) {
      dp_mol->copyProp(prop.name(), source.asRDMol(), prop.name());
    }
  }
}

void ROMol::updateProps(const RDProps &source, bool preserveExisting) {
  // Remove all properties if !preserveExisting
  if (!preserveExisting) {
    dp_mol->clearProps();
  }

  const STR_VECT keys = source.getPropList();
  std::vector<std::string> computedPropList;
  source.getPropIfPresent<std::vector<std::string>>(
      RDKit::detail::computedPropName, computedPropList);
  for (const auto &key : keys) {
    if (key == RDKit::detail::computedPropName) {
      continue;
    }
    const PropToken token(key);
    const RDValue &val = source.getPropRDValue(key);
    bool isComputed =
        std::find(computedPropList.begin(), computedPropList.end(), key) !=
        computedPropList.end();
    dp_mol->setMolProp(token, val, isComputed);
  }
}

void ROMol::clear(){dp_mol->clearProps();}

void ROMol::debugMol(std::ostream& str) const {
  str << "Atoms:" << std::endl;
  for (const auto atom : atoms()) {
    str << "\t" << *atom << std::endl;
  }

  str << "Bonds:" << std::endl;
  for (const auto bond : bonds()) {
    str << "\t" << *bond << std::endl;
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

// --------------
// Atom bookmarks
// --------------

void ROMol::setAtomBookmark(Atom *at, int mark) {
  // The atom may not be added to its owning molecule yet if it's an isolated
  // atom, but the ROMol still needs to be able to store a bookmark for it,
  // even though this isn't supported by the public RDMol interface.
  PRECONDITION(
      (&at->getOwningMol() == this &&
       ((&at->getOwningMol().asRDMol() != &at->getDataRDMol()) ||
        at->getOwningMol().getAtomWithIdx(at->getIdx()) == at)),
      "Bookmark atom must be for this ROMol, even if it's not added to the ROMol");

  auto &marks = (*dp_mol->getAtomBookmarksCompat())[mark];
  marks.push_back(at);
}

void ROMol::replaceAtomBookmark(Atom *at, int mark) {
  // The atom may not be added to its owning molecule yet if it's an isolated
  // atom, but the ROMol still needs to be able to store a bookmark for it,
  // even though this isn't supported by the public RDMol interface.
  PRECONDITION(
      (&at->getOwningMol() == this &&
       ((&at->getOwningMol().asRDMol() != &at->getDataRDMol()) ||
        at->getOwningMol().getAtomWithIdx(at->getIdx()) == at)),
      "Bookmark atom must be for this ROMol, even if it's not added to the ROMol");

  auto &marks = (*dp_mol->getAtomBookmarksCompat())[mark];
  marks.clear();
  marks.push_back(at);
}
Atom * ROMol::getAtomWithBookmark(int mark) {
  auto *allMarks = dp_mol->getAtomBookmarksCompat();
  PRECONDITION(allMarks != nullptr, "null atom bookmarks map");
  auto it = allMarks->find(mark);
  PRECONDITION(it != allMarks->end(), "atom bookmark not found");
  auto &marks = it->second;
  PRECONDITION(marks.size() != 0, "atom bookmark not found");
  return marks.front();
}

Atom * ROMol::getUniqueAtomWithBookmark(int mark) {
  return getAtomWithBookmark(mark);
}

ROMol::ATOM_PTR_LIST & ROMol::getAllAtomsWithBookmark(int mark) {
  return dp_mol->getAtomBookmarksCompat(mark);
}
void ROMol::clearAtomBookmark(int mark) {
  dp_mol->clearAtomBookmark(mark);
}
void ROMol::clearAtomBookmark(int mark, const Atom *at) {
  // The atom may not be added to its owning molecule yet if it's an isolated
  // atom, but the ROMol still needs to be able to store a bookmark for it,
  // even though this isn't supported by the public RDMol interface.
  PRECONDITION(
      (&at->getOwningMol() == this &&
       ((&at->getOwningMol().asRDMol() != &at->getDataRDMol()) ||
        at->getOwningMol().getAtomWithIdx(at->getIdx()) == at)),
      "Bookmark atom must be for this ROMol, even if it's not added to the ROMol");
  if (!dp_mol->hasAtomBookmark(mark)) {
    return;
  }

  auto &marks = dp_mol->getAtomBookmarksCompat(mark);
  auto removeIt = std::find(marks.begin(), marks.end(), at);
  if (removeIt != marks.end()) {
    marks.erase(removeIt);
    if (marks.size() == 0) {
      dp_mol->clearAtomBookmark(mark);
    }
  }
}

void ROMol::clearAllAtomBookmarks() {
  dp_mol->clearAllAtomBookmarks();
}

bool ROMol::hasAtomBookmark(int mark) const {
  return dp_mol->hasAtomBookmark(mark);
}
ROMol::ATOM_BOOKMARK_MAP *ROMol::getAtomBookmarks() {
  return dp_mol->getAtomBookmarksCompat();
}

// --------------
// Bond bookmarks
// --------------

void ROMol::setBondBookmark(Bond *bond, int mark) {
  // The bond may not be added to its owning molecule yet if it's an isolated
  // bond, but the ROMol still needs to be able to store a bookmark for it,
  // even though this isn't supported by the public RDMol interface.
  PRECONDITION(
      (&bond->getOwningMol() == this &&
       ((&bond->getOwningMol().asRDMol() != &bond->getDataRDMol()) ||
        bond->getOwningMol().getBondWithIdx(bond->getIdx()) == bond)),
      "Bookmark bond must be for this ROMol, even if it's not added to the ROMol");

  auto &marks = (*dp_mol->getBondBookmarksCompat())[mark];
  marks.push_back(bond);
}
Bond * ROMol::getBondWithBookmark(int mark) {
  auto *allMarks = dp_mol->getBondBookmarksCompat();
  PRECONDITION(allMarks != nullptr, "null bond bookmarks map");
  auto it = allMarks->find(mark);
  PRECONDITION(it != allMarks->end(), "bond bookmark not found");
  auto &marks = it->second;
  PRECONDITION(marks.size() != 0, "bond bookmark not found");
  return marks.front();
}

Bond * ROMol::getUniqueBondWithBookmark(int mark) {
  return getBondWithBookmark(mark);
}
ROMol::BOND_PTR_LIST & ROMol::getAllBondsWithBookmark(int mark) {
  return dp_mol->getBondBookmarksCompat(mark);
}
void ROMol::clearBondBookmark(int mark) {
  dp_mol->clearBondBookmark(mark);
}
void ROMol::clearBondBookmark(int mark, const Bond *bond) {
  // The bond may not be added to its owning molecule yet if it's an isolated
  // bond, but the ROMol still needs to be able to store a bookmark for it,
  // even though this isn't supported by the public RDMol interface.
  PRECONDITION(
      (&bond->getOwningMol() == this &&
       ((&bond->getOwningMol().asRDMol() != &bond->getDataRDMol()) ||
        bond->getOwningMol().getBondWithIdx(bond->getIdx()) == bond)),
      "Bookmark bond must be for this ROMol, even if it's not added to the ROMol");
  if (!dp_mol->hasBondBookmark(mark)) {
    return;
  }

  auto &marks = dp_mol->getBondBookmarksCompat(mark);
  auto removeIt = std::find(marks.begin(), marks.end(), bond);
  if (removeIt != marks.end()) {
    marks.erase(removeIt);
    if (marks.size() == 0) {
      dp_mol->clearBondBookmark(mark);
    }
  }
}
void ROMol::clearAllBondBookmarks() {
  dp_mol->clearAllBondBookmarks();
}
bool ROMol::hasBondBookmark(int mark) const {
  return dp_mol->hasBondBookmark(mark);
}
ROMol::BOND_BOOKMARK_MAP * ROMol::getBondBookmarks() {
  return dp_mol->getBondBookmarksCompat();
}


const Conformer& ROMol::getConformer(int id) const {
  const auto& confs = dp_mol->getConformersCompat();
  if (confs.size() == 0) {
    throw ConformerException("No conformations available on the molecule");
  }

  if (id < 0) {
    return *(confs.front());
  }
  auto cid = (unsigned int)id;
  for (auto& conf : confs) {
    if (conf->getId() == cid) {
      return *conf;
    }
  }
  // we did not find a conformation with the specified ID
  std::string mesg = "Can't find conformation with ID: ";
  mesg += id;
  throw ConformerException(mesg);
}

Conformer& ROMol::getConformer(int id) {
  auto& confs = dp_mol->getConformersCompat();
  if (confs.size() == 0) {
    throw ConformerException("No conformations available on the molecule");
  }

  if (id < 0) {
    return *(confs.front());
  }
  auto cid = (unsigned int)id;
  for (auto& conf : confs) {
    if (conf->getId() == cid) {
      return *conf;
    }
  }
  // we did not find a conformation with the specified ID
  std::string mesg = "Can't find conformation with ID: ";
  mesg += id;
  throw ConformerException(mesg);
}

void ROMol::removeConformer(unsigned int id) {
  dp_mol->removeConformer(id);
}

void ROMol::clearConformers() {
  dp_mol->clearConformers();
}

unsigned int ROMol::addConformer(Conformer* conf, bool assignId) {
  PRECONDITION(conf, "bad conformer");
  PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
               "Number of atom mismatch");
  auto& confs = dp_mol->getConformersCompat();
  int maxId = -1;
  if (assignId) {
    for (const auto& cptr : confs) {
      maxId = std::max((int)(cptr->getId()), maxId);
    }
    maxId++;
    conf->setId((unsigned int)maxId);
  } else {
    maxId = conf->getId();
  }
  conf->setOwningMol(this);
  confs.push_back(CONFORMER_SPTR(conf));

  auto [returnId, confData] = dp_mol->addConformer(maxId, conf->is3D());
  PRECONDITION(int(returnId) == maxId, "Id mismatch");
  PRECONDITION(getNumAtoms() == 0 || confData != nullptr, "Null conformer positions");

  // Special case of addConformer - non-const confData means that it marked the
  // conformer as being out of sync after synchronizing the positions, but
  // they're actually in sync. However, the caller may modify the positions from
  // the Conformer pointer it passed in without requesting it again, so mark
  // it as being modified.
  dp_mol->markConformersAsCompatModified();

  return conf->getId();
}

void ROMol::clearSubstanceGroups() {
  dp_mol->getSubstanceGroups().clear();
}

unsigned int ROMol::addAtom(Atom* atom, bool updateLabel, bool takeOwnership) {
  PRECONDITION(atom, "NULL atom provided");
  // A test in catch_graphmol.cpp requires that this precondition, copied from
  // replaceAtomPointerCompat, must fail before adding the atom.
  PRECONDITION((!takeOwnership || atom->dp_owningMol == nullptr ||
                atom->dp_owningMol == dp_mol),
               "Atom cannot have another owning mol");

  const uint32_t newAtomIdx = dp_mol->getNumAtoms();
  dp_mol->addAtom();
  getAtomWithIdx(newAtomIdx)->initFromOther(*atom);

  if (takeOwnership) {
    // Allow for stable pointer by resetting the compat entry.
    dp_mol->replaceAtomPointerCompat(atom, newAtomIdx);
  } else if (atom->hasQuery()) {
    dp_mol->setAtomQueryCompatFull(newAtomIdx, atom->getQuery()->copy());
  }

  if (updateLabel) {
    dp_mol->clearAtomBookmark(ci_RIGHTMOST_ATOM);
    dp_mol->setAtomBookmark(newAtomIdx, ci_RIGHTMOST_ATOM);
  }
  return newAtomIdx;
}

unsigned int ROMol::addBond(Bond* bond, bool takeOwnership) {
  PRECONDITION(bond, "NULL bond passed in");
  // A test in catch_graphmol.cpp requires that this precondition, copied from
  // replaceBondPointerCompat, must fail before adding the bond.
  PRECONDITION(!takeOwnership || bond->dp_owningMol == nullptr ||
                   bond->dp_owningMol == dp_mol,
               "Bond cannot have another owning mol");

  const uint32_t newBondIdx = getNumBonds();
  // Note that the ROMol version does not do aromaticity updates.
  dp_mol->addBond(bond->getBeginAtomIdx(),
                  bond->getEndAtomIdx(),
                  BondEnums::BondType(bond->getBondType()),
                  /*updateAromaticity=*/false);
  getBondWithIdx(newBondIdx)->initFromOther(*bond);

  if (takeOwnership) {
    // Allow for stable pointer by resetting the compat entry.
    dp_mol->replaceBondPointerCompat(bond, newBondIdx);
  } else if (bond->hasQuery()) {
    dp_mol->setBondQueryCompatFull(newBondIdx, bond->getQuery()->copy());
  }

  return getNumBonds();
}

void ROMol::initFromOther(
    [[maybe_unused]] const ROMol& other,
    [[maybe_unused]] bool quickCopy,
    [[maybe_unused]] int confId) {
  raiseNonImplementedFunction("initFromOther");
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

  MolPickler::molFromPickle(pkl, *this, PicklerOps::AllProps);
}

template RDKIT_GRAPHMOL_EXPORT void ROMol::save<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, const unsigned int) const;
template RDKIT_GRAPHMOL_EXPORT void ROMol::load<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, const unsigned int);
#endif // RDK_USE_BOOST_SERIALIZATION


uint32_t num_vertices(const ROMol& mol) {
  return mol.getNumAtoms();
}

uint32_t num_edges(const ROMol& mol) {
  return mol.getNumBonds();
}

}  // namespace RDKit
