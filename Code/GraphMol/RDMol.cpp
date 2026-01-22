//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit
//  contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDMol.h"
#include "Atom.h"
#include "Bond.h"
#include "Conformer.h"

#include "PeriodicTable.h"
#include "SanitException.h"
#include "RDGeneral/vector_utils.h"

#include "QueryAtom.h"
#include "QueryBond.h"
#include "QueryOps.h"

#include <boost/tokenizer.hpp>

#include <mutex>
#include <type_traits>

namespace Queries {
// Forward declaration
template <class NonCompatType>
class RDKIT_QUERY_EXPORT NonCompatQuery;

//! This class is for wrapping an existing Query that was created using the new
//! atom/bond types, for example, when constructing CompatibilityData, or when
//! setting a Query when CompatibilityData already exists.
//!
//! It's unlikely, but possible, for callers to call the modification functions,
//! so they need to work, but don't need to be efficient. This also needs to
//! reflect changes applied to the other Query. If callers hold onto pointers to
//! both, that's not feasible to support, though.
template <class CompatType>
class RDKIT_QUERY_EXPORT CompatQuery : public Query<int, CompatType, true> {
 public:
  using MatchFuncArgType = int;
  using DataFuncArgType = CompatType;
  static constexpr bool needsConversion = true;
  using BASE = Query<int, CompatType, true>;
  static constexpr bool isAtom =
      std::is_same_v<CompatType, const RDKit::Atom *>;
  using NonCompatType =
      std::conditional_t<isAtom, RDKit::ConstRDMolAtom, RDKit::ConstRDMolBond>;
  using INNER_TYPE = Query<int, NonCompatType, true>;

  INNER_TYPE *inner;
  bool ownsInner;

  CompatQuery(INNER_TYPE *innerIn, bool ownsInnerIn)
      : inner(innerIn), ownsInner(ownsInnerIn) {
    BASE::d_description = inner->getDescription();
    BASE::d_queryType = inner->getTypeLabel();
  }
  CompatQuery(CompatQuery &&) = default;
  ~CompatQuery() {
    if (ownsInner) {
      delete inner;
    }
  }

  bool Match(const CompatType what) const override {
    PRECONDITION(BASE::d_children.size() == 0,
                 "Children must be added to the inner query, not this");
    PRECONDITION(!BASE::getNegation(),
                 "Negation must be applied to the inner query, not this");
    return inner->Match(
        NonCompatType(&what->getOwningMol().asRDMol(), what->getIdx()));
  }
  BASE *copy() const override {
    PRECONDITION(BASE::d_children.size() == 0,
                 "Children must be added to the inner query, not this");
    // This copy will own the inner query
    CompatQuery<CompatType> *res =
        new CompatQuery<CompatType>(inner->copy(), true);
    res->d_description = this->d_description;
    res->d_queryType = this->d_queryType;
    return res;
  }

  void setMatchFunc([[maybe_unused]] bool (*what)(MatchFuncArgType)) override {
    raiseNonImplementedFunction("CompatQuery::setMatchFunc");
  }
  void setDataFunc(
      [[maybe_unused]] MatchFuncArgType (*what)(DataFuncArgType)) override {
    raiseNonImplementedFunction("CompatQuery::setDataFunc");
  }
  void setNegation(bool what) override { inner->setNegation(what); }
  void addChild(typename BASE::CHILD_TYPE child) override;

  const INNER_TYPE *getInnerQuery() const override { return inner; }
};

//! This class is for wrapping an existing Query that was created using the
//! compatibility atom/bond types, for example, when QueryAtom::setQuery is
//! called.
template <class NonCompatType>
class RDKIT_QUERY_EXPORT NonCompatQuery
    : public Query<int, NonCompatType, true> {
 public:
  using MatchFuncArgType = int;
  using DataFuncArgType = NonCompatType;
  static constexpr bool needsConversion = true;
  using BASE = Query<int, NonCompatType, true>;
  static constexpr bool isAtom =
      std::is_same_v<NonCompatType, RDKit::ConstRDMolAtom>;
  using CompatType =
      std::conditional_t<isAtom, const RDKit::Atom *, const RDKit::Bond *>;
  using INNER_TYPE = Query<int, CompatType, true>;

  INNER_TYPE *inner;
  bool ownsInner;

  NonCompatQuery(INNER_TYPE *innerIn, bool ownsInnerIn)
      : inner(innerIn), ownsInner(ownsInnerIn) {
    BASE::d_description = inner->getDescription();
    BASE::d_queryType = inner->getTypeLabel();
  }
  NonCompatQuery(NonCompatQuery &&) = default;
  ~NonCompatQuery() {
    if (ownsInner) {
      delete inner;
    }
  }

  bool Match(const NonCompatType what) const override {
    PRECONDITION(BASE::d_children.size() == 0,
                 "Children must be added to the inner query, not this");
    PRECONDITION(!BASE::getNegation(),
                 "Negation must be applied to the inner query, not this");
    if constexpr (isAtom) {
      return inner->Match(what.mol().asROMol().getAtomWithIdx(what.index()));
    } else {
      return inner->Match(what.mol().asROMol().getBondWithIdx(what.index()));
    }
  }
  BASE *copy() const override {
    PRECONDITION(BASE::d_children.size() == 0,
                 "Children must be added to the inner query, not this");
    // This copy will own the inner query
    NonCompatQuery<NonCompatType> *res =
        new NonCompatQuery<NonCompatType>(inner->copy(), true);
    res->d_description = this->d_description;
    res->d_queryType = this->d_queryType;
    return res;
  }

  void setMatchFunc([[maybe_unused]] bool (*what)(MatchFuncArgType)) override {
    raiseNonImplementedFunction("NonCompatQuery::setMatchFunc");
  }
  void setDataFunc(
      [[maybe_unused]] MatchFuncArgType (*what)(DataFuncArgType)) override {
    raiseNonImplementedFunction("NonCompatQuery::setDataFunc");
  }
  void setNegation(bool what) override { inner->setNegation(what); }
  void addChild(typename BASE::CHILD_TYPE child) override {
    inner->addChild(
        std::make_shared<CompatQuery<CompatType>>(child->copy(), true));
  }

  const INNER_TYPE *getInnerQuery() const override { return inner; }
};

// Defined after NonCompatQuery class due to dependence
template <class CompatType>
void CompatQuery<CompatType>::addChild(typename BASE::CHILD_TYPE child) {
  inner->addChild(
      std::make_shared<NonCompatQuery<NonCompatType>>(child->copy(), true));
}

template <class NonCompatType>
Query<int,
      std::conditional_t<std::is_same_v<NonCompatType, RDKit::ConstRDMolAtom>,
                         const RDKit::Atom *, const RDKit::Bond *>,
      true> *
makeNewCompatQuery(Query<int, NonCompatType, true> *inner) {
  static constexpr bool isAtom =
      std::is_same_v<NonCompatType, RDKit::ConstRDMolAtom>;
  using CompatType =
      std::conditional_t<isAtom, const RDKit::Atom *, const RDKit::Bond *>;
  using INNER_TYPE = NonCompatQuery<NonCompatType>;
  using OUTER_TYPE = CompatQuery<CompatType>;
  auto *nonCompatInner = dynamic_cast<INNER_TYPE *>(inner);
  if (nonCompatInner != nullptr) {
    // Don't re-wrap if already wrapped
    if (nonCompatInner->ownsInner) {
      // This was from a copy, so take ownership of inner, to avoid duplicating
      // it again
      nonCompatInner->ownsInner = false;
      return nonCompatInner->inner;
    }
    // Not already a copy, so make a new copy of the original type
    return nonCompatInner->inner->copy();
  }
  // Not already wrapped, so construct as normal
  return new OUTER_TYPE(inner, false);
}

template <class CompatType>
Query<int,
      std::conditional_t<std::is_same_v<CompatType, const RDKit::Atom *>,
                         RDKit::ConstRDMolAtom, RDKit::ConstRDMolBond>,
      true> *
makeNewNonCompatQuery(Query<int, CompatType, true> *inner) {
  if (inner == nullptr) {
    return nullptr;
  }
  static constexpr bool isAtom =
      std::is_same_v<CompatType, const RDKit::Atom *>;
  using NonCompatType =
      std::conditional_t<isAtom, RDKit::ConstRDMolAtom, RDKit::ConstRDMolBond>;
  using INNER_TYPE = CompatQuery<CompatType>;
  using OUTER_TYPE = NonCompatQuery<NonCompatType>;
  auto *compatInner = dynamic_cast<INNER_TYPE *>(inner);
  if (compatInner != nullptr) {
    // Don't re-wrap if already wrapped
    if (compatInner->ownsInner) {
      // This was from a copy, so take ownership of inner, to avoid duplicating
      // it again
      compatInner->ownsInner = false;
      return compatInner->inner;
    }
    // Not already a copy, so make a new copy of the original type
    return compatInner->inner->copy();
  }
  // Not already wrapped, so construct as normal
  return new OUTER_TYPE(inner, false);
}

}  // namespace Queries

namespace RDKit {

double AtomData::getMass() const {
  auto *table = PeriodicTable::getTable();
  if (isotope) {
    double res = table->getMassForIsotope(atomicNum, isotope);
    if (atomicNum != 0 && res == 0.0) {
      res = isotope;
    }
    return res;
  }
  return table->getAtomicWeight(atomicNum);
}

namespace {

enum class CompatSyncStatus {
  lastUpdatedRDMol = 0,
  lastUpdatedCompat,
  inSync
};

void copyRingInfoToCompatibilityData(const RingInfoCache &input,
                                     RingInfo &output,
                                     [[maybe_unused]] uint32_t numAtoms,
                                     [[maybe_unused]] uint32_t numBonds) {
  output.reset();

  // This handles the case where input is uninitialized, in case
  // the two structures must be synchronized to be equally uninitialized
  if (!input.isInitialized()) {
    return;
  }

  output.initialize(input.findRingType);
  INT_VECT atomIndices;
  INT_VECT bondIndices;
  const size_t numRings = input.numRings();

  for (size_t ringIndex = 0; ringIndex < numRings; ++ringIndex) {
    atomIndices.clear();
    bondIndices.clear();
    const uint32_t ringBeginIndex = input.ringBegins[ringIndex];
    const uint32_t ringEndIndex = input.ringBegins[ringIndex + 1];
    for (uint32_t atomBondIndex = ringBeginIndex; atomBondIndex < ringEndIndex;
         ++atomBondIndex) {
      atomIndices.push_back(input.atomsInRings[atomBondIndex]);
      bondIndices.push_back(input.bondsInRings[atomBondIndex]);
    }
    output.addRing(atomIndices, bondIndices);
  }
}
void copyRingInfoFromCompatibilityData(const RingInfo &input,
                                       RingInfoCache &output, uint32_t numAtoms,
                                       uint32_t numBonds) {
  output.reset();

  // This handles the case where input is uninitialized, in case
  // the two structures must be synchronized to be equally uninitialized
  if (!input.isInitialized()) {
    return;
  }

  output.isInit = true;
  output.findRingType = input.getRingType();

  const size_t numRings = input.numRings();
  output.ringBegins.reserve(numRings + 1);
  const VECT_INT_VECT &atomRings = input.atomRings();
  const VECT_INT_VECT &bondRings = input.bondRings();
  PRECONDITION(atomRings.size() == numRings, "atomRings size mismatch");
  PRECONDITION(bondRings.size() == numRings, "bondRings size mismatch");
  size_t offset = 0;
  for (size_t i = 0; i < numRings; ++i) {
    output.ringBegins.push_back(offset);
    offset += atomRings[i].size();
    PRECONDITION(atomRings[i].size() == bondRings[i].size(),
                 "atom vs. bond ring mismatch");
  }
  const size_t totalSize = offset;
  output.ringBegins.push_back(totalSize);

  output.atomsInRings.reserve(totalSize);
  output.bondsInRings.reserve(totalSize);
  for (size_t i = 0; i < numRings; ++i) {
    for (auto atom : atomRings[i]) {
      output.atomsInRings.push_back(atom);
    }
    for (auto bond : bondRings[i]) {
      output.bondsInRings.push_back(bond);
    }
  }

  output.atomMembershipBegins.clear();
  output.bondMembershipBegins.clear();
  output.atomMembershipBegins.resize(numAtoms + 1, 0);
  output.bondMembershipBegins.resize(numBonds + 1, 0);
  // Count the atoms and bonds
  for (size_t i = 0; i < numRings; ++i) {
    for (auto atom : atomRings[i]) {
      ++output.atomMembershipBegins[atom + 1];
    }
    for (auto bond : bondRings[i]) {
      ++output.bondMembershipBegins[bond + 1];
    }
  }
  // Partial sum to figure out where they go
  std::partial_sum(output.atomMembershipBegins.begin() + 1,
                   output.atomMembershipBegins.end(),
                   output.atomMembershipBegins.begin() + 1);
  std::partial_sum(output.bondMembershipBegins.begin() + 1,
                   output.bondMembershipBegins.end(),
                   output.bondMembershipBegins.begin() + 1);

  // Fill in the ring memberships for each atom and bond
  output.atomMemberships.resize(totalSize);
  output.bondMemberships.resize(totalSize);
  for (size_t i = 0; i < numRings; ++i) {
    for (auto atom : atomRings[i]) {
      auto &next = output.atomMembershipBegins[atom];
      output.atomMemberships[next] = i;
      ++next;
    }
    for (auto bond : bondRings[i]) {
      auto &next = output.bondMembershipBegins[bond];
      output.bondMemberships[next] = i;
      ++next;
    }
  }
  // Shift forward 1 to undo the increments
  uint32_t prev = 0;
  for (size_t i = 0; i < numAtoms; ++i) {
    std::swap(prev, output.atomMembershipBegins[i]);
  }
  prev = 0;
  for (size_t i = 0; i < numBonds; ++i) {
    std::swap(prev, output.bondMembershipBegins[i]);
  }

  output.initFusedInfoFromBondMemberships();
}
}  // namespace

struct RDMol::CompatibilityData {
  ROMol *compatMol = nullptr;
  std::unique_ptr<RWMol> mol;
  std::vector<std::unique_ptr<Atom>> atoms;
  std::vector<std::unique_ptr<Bond>> bonds;
  std::map<int, std::list<Atom *>> atomBookmarks;
  std::map<int, std::list<Bond *>> bondBookmarks;
  // Each INT_VECT is wrapped in a unique_ptr so that removing bonds doesn't
  // change the address of the INT_VECT. testGitHubIssue8() in molopstest.cpp
  // relies on this persistence.
  std::vector<std::unique_ptr<INT_VECT>> bondStereoAtoms;
  // stereoGroups is initialized on demand in getStereoGroupsCompat
  std::atomic<std::vector<StereoGroup> *> stereoGroups = nullptr;
  mutable std::atomic<CompatSyncStatus> stereoGroupsSyncStatus;

  ROMol::CONF_SPTR_LIST conformers;
  mutable std::atomic<CompatSyncStatus> conformerSyncStatus;

  std::unique_ptr<RingInfo> ringInfo;
  mutable std::atomic<CompatSyncStatus> ringInfoSyncStatus;

  CompatibilityData(RDMol &rdmol, ROMol *existingROMolPtr = nullptr)
      : atoms(rdmol.getNumAtoms()),
        bonds(rdmol.getNumBonds()),
        bondStereoAtoms(rdmol.getNumBonds()),
        stereoGroupsSyncStatus(CompatSyncStatus::inSync),
        conformerSyncStatus(CompatSyncStatus::inSync),
        ringInfoSyncStatus(CompatSyncStatus::inSync) {
    // For a self-owning compat mol, we need pointer compatibility.
    if (existingROMolPtr) {
      compatMol = existingROMolPtr;
      if (!existingROMolPtr->d_ownedBySelf) {
        auto rwmolPtr = dynamic_cast<RWMol *>(existingROMolPtr);
        if (!rwmolPtr) {
          throw ValueErrorException(
              "Existing ROMol ptr must be an RWMol if not owned by self");
        }
        mol.reset(rwmolPtr);
      }
    } else {
      mol.reset(new RWMol(&rdmol));
      mol->d_ownedBySelf = false;
      compatMol = mol.get();
    }

    for (uint32_t i = 0, n = rdmol.getNumBonds(); i < n; ++i) {
      bondStereoAtoms[i] = std::make_unique<INT_VECT>();
    }

    const uint32_t numAtoms = rdmol.getNumAtoms();
    for (uint32_t i = 0, n = numAtoms; i < n; ++i) {
      auto *query = rdmol.getAtomQuery(i);
      Atom *atom;
      if (query == nullptr) {
        atom = new Atom(&rdmol, i);
      } else {
        atom = new QueryAtom(&rdmol, i, makeNewCompatQuery(query));
      }
      atoms[i].reset(atom);
    }
    for (uint32_t i = 0, n = rdmol.getNumBonds(); i < n; ++i) {
      auto *query = rdmol.getBondQuery(i);
      Bond *bond;
      if (query == nullptr) {
        bond = new Bond(&rdmol, i);
      } else {
        bond = new QueryBond(&rdmol, i, makeNewCompatQuery(query));
      }
      bonds[i].reset(bond);
      const uint32_t *stereoAtoms = rdmol.getBond(i).stereoAtoms;
      if (*stereoAtoms != atomindex_t(-1)) {
        bondStereoAtoms[i]->reserve(2);
        bondStereoAtoms[i]->push_back(*stereoAtoms);
        bondStereoAtoms[i]->push_back(*(stereoAtoms + 1));
      }
    }

    const double *conformerPositions = rdmol.conformerPositionData.data();
    const size_t conformerOffsetStride = rdmol.conformerAtomCapacity * 3;
    for (uint32_t conformerIndex = 0, n = rdmol.getNumConformers();
         conformerIndex < n; ++conformerIndex) {
      const size_t conformerOffset = conformerIndex * conformerOffsetStride;
      auto *conformer = new Conformer();
      conformer->getPositions().reserve(numAtoms);
      for (uint32_t atomIndex = 0; atomIndex < numAtoms; ++atomIndex) {
        conformer->getPositions().push_back(RDGeom::Point3D(
            conformerPositions[conformerOffset + 3 * atomIndex + 0],
            conformerPositions[conformerOffset + 3 * atomIndex + 1],
            conformerPositions[conformerOffset + 3 * atomIndex + 2]));
      }
      if (rdmol.conformerIdsAndFlags.size() > 0) {
        conformer->setId(rdmol.conformerIdsAndFlags[conformerIndex].id);
        conformer->set3D(rdmol.conformerIdsAndFlags[conformerIndex].is3D);
      } else {
        conformer->setId(conformerIndex);
        conformer->set3D(true);
      }
      conformer->setOwningMol(compatMol);
      conformers.push_back(CONFORMER_SPTR(conformer));
    }

    for (const auto &[mark, atomIndices] : rdmol.atomBookmarks) {
      for (int atomIndex : atomIndices) {
        atomBookmarks[mark].push_back(atoms[atomIndex].get());
      }
    }
    for (const auto &[mark, bondIndices] : rdmol.bondBookmarks) {
      for (int bondIndex : bondIndices) {
        bondBookmarks[mark].push_back(bonds[bondIndex].get());
      }
    }

    ringInfo.reset(new RingInfo());
    if (rdmol.ringInfo.isInitialized()) {
      copyRingInfoToCompatibilityData(rdmol.ringInfo, *ringInfo,
                                      rdmol.getNumAtoms(), rdmol.getNumBonds());
    }
  }

  ~CompatibilityData() { delete stereoGroups.load(std::memory_order_relaxed); }

  // Copy conformers (including properties) from another CompatibilityData
  void copyConformersFrom(const CompatibilityData *source) {
    conformers.clear();
    for (const auto &otherConf : source->conformers) {
      auto *confCopy =
          new Conformer(*otherConf);  // Copy constructor preserves properties
      confCopy->setOwningMol(compatMol);
      conformers.push_back(CONFORMER_SPTR(confCopy));
    }
  }

  void setNewOwner(RDMol &rdmol) {
    compatMol->dp_mol = &rdmol;
    for (auto &atom : atoms) {
      atom->dp_dataMol = &rdmol;
      atom->dp_owningMol = &rdmol;
    }
    for (auto &bond : bonds) {
      bond->dp_dataMol = &rdmol;
      bond->dp_owningMol = &rdmol;
    }
    for (auto &conformer : conformers) {
      conformer->setOwningMol(compatMol);
    }
  }
};

void RDMol::releaseCompatMolOwnership() {
  std::ignore = ensureCompatInit()->mol.release();
};

void RDMol::setCompatPointersToSelf() {
  ensureCompatInit()->setNewOwner(*this);
}

void RDMol::setROMolPointerCompat(ROMol *ptr) {
  auto *compat = ensureCompatInit();
  compat->compatMol = ptr;
  if (!ptr->d_ownedBySelf) {
    auto *rwmol = dynamic_cast<RWMol *>(ptr);
    if (rwmol == nullptr) {
      throw ValueErrorException(
          "Existing ROMol ptr must be an RWMol if not owned by self");
    }
    compat->mol.reset(rwmol);
  }
  for (auto &group : substanceGroups) {
    group.setOwningMol(compat->compatMol);
  }
}

Atom *RDMol::getAtomCompat(uint32_t atomIndex) {
  return ensureCompatInit()->atoms[atomIndex].get();
}
const Atom *RDMol::getAtomCompat(uint32_t atomIndex) const {
  return ensureCompatInit()->atoms[atomIndex].get();
}

Bond *RDMol::getBondCompat(uint32_t bondIndex) {
  return ensureCompatInit()->bonds[bondIndex].get();
}
const Bond *RDMol::getBondCompat(uint32_t bondIndex) const {
  return ensureCompatInit()->bonds[bondIndex].get();
}

std::list<Atom *> &RDMol::getAtomBookmarksCompat(int mark) {
  RDMol::CompatibilityData *compat = ensureCompatInit();
  auto found = compat->atomBookmarks.find(mark);
  PRECONDITION(found != compat->atomBookmarks.end(),
               "No atom bookmark with mark " + std::to_string(mark));
  return found->second;
}

std::list<Bond *> &RDMol::getBondBookmarksCompat(int mark) {
  RDMol::CompatibilityData *compat = ensureCompatInit();
  auto found = compat->bondBookmarks.find(mark);
  PRECONDITION(found != compat->bondBookmarks.end(),
               "No bond bookmark with mark " + std::to_string(mark));
  return found->second;
}

std::map<int, std::list<Atom *>> *RDMol::getAtomBookmarksCompat() {
  RDMol::CompatibilityData *compat = ensureCompatInit();
  return &compat->atomBookmarks;
}

std::map<int, std::list<Bond *>> *RDMol::getBondBookmarksCompat() {
  RDMol::CompatibilityData *compat = ensureCompatInit();
  return &compat->bondBookmarks;
}

void RDMol::setAtomQueryCompat(atomindex_t atomIndex,
                               Atom::QUERYATOM_QUERY *query) {
  if (atomQueries.size() != getNumAtoms()) {
    atomQueries.resize(getNumAtoms());
  }
  atomQueries[atomIndex].reset(makeNewNonCompatQuery(query));
}

void RDMol::setBondQueryCompat(uint32_t bondIndex,
                               Bond::QUERYBOND_QUERY *query) {
  if (bondQueries.size() != getNumBonds()) {
    bondQueries.resize(getNumBonds());
  }
  bondQueries[bondIndex].reset(makeNewNonCompatQuery(query));
}

void RDMol::setAtomQueryCompatFull(
    atomindex_t atomIndex, Queries::Query<int, const Atom *, true> *query) {
  auto *compat = ensureCompatInit();
  Atom *atom = compat->atoms[atomIndex].get();
  QueryAtom *queryAtom = dynamic_cast<QueryAtom *>(atom);
  if (queryAtom == nullptr) {
    // Convert to a QueryAtom now that there's a query
    compat->atoms[atomIndex].reset(new QueryAtom(this, atomIndex, query));
  } else {
    queryAtom->setQueryPrivate(query);
  }

  setAtomQueryCompat(atomIndex, query);
}

void RDMol::setBondQueryCompatFull(
    uint32_t bondIndex, Queries::Query<int, const Bond *, true> *query) {
  auto *compat = ensureCompatInit();
  Bond *bond = compat->bonds[bondIndex].get();
  QueryBond *queryBond = dynamic_cast<QueryBond *>(bond);
  if (queryBond == nullptr) {
    // Convert to a QueryBond now that there's a query
    compat->bonds[bondIndex].reset(new QueryBond(this, bondIndex, query));
  } else {
    queryBond->setQueryPrivate(query);
  }

  setBondQueryCompat(bondIndex, query);
}

void RDMol::clearBondStereoAtomsCompat(uint32_t bondIndex) {
  ensureCompatInit()->bondStereoAtoms[bondIndex]->clear();
}
bool RDMol::hasBondStereoAtomsCompat(uint32_t bondIndex) const {
  return !ensureCompatInit()->bondStereoAtoms[bondIndex]->empty();
}
const INT_VECT *RDMol::getBondStereoAtomsCompat(uint32_t bondIndex) const {
  // Because RDMol::setBondStereoAtoms and RDMol::clearBondStereoAtoms update
  // these lists, no synchronization is needed here.
  return ensureCompatInit()->bondStereoAtoms[bondIndex].get();
}
INT_VECT *RDMol::getBondStereoAtomsCompat(uint32_t bondIndex) {
  // Because RDMol::setBondStereoAtoms and RDMol::clearBondStereoAtoms update
  // these lists, no synchronization is needed here.
  return ensureCompatInit()->bondStereoAtoms[bondIndex].get();
}

void RDMol::replaceAtomPointerCompat(Atom *atom, uint32_t index) {
  PRECONDITION(atom->dp_owningMol == nullptr || atom->dp_owningMol == this,
               "Atom cannot have another owning mol");
  PRECONDITION(atom->dp_dataMol != this, "About to delete this RDMol");

  // First clean up previous data. We're only taking the pointer, so no need to
  // copy the data.
  delete atom->dp_dataMol;

  auto *compat = ensureCompatInit();

  // Store original pointer for bookmark replacement.
  auto &atom_ptr = compat->atoms[index];
  Atom *origAtom = atom_ptr.get();

  // Now set the pointers for this RDMol.
  // The Atom pointed to by origAtom is invalid after calling reset,
  // but only the pointer value is used below.
  atom_ptr.reset(atom);
  atom->setIdx(index);
  atom->dp_dataMol = this;
  atom->dp_owningMol = this;

  // Replace pointers in bookmarks.
  for (auto &ab : compat->atomBookmarks) {
    for (auto &elem : ab.second) {
      if (elem == origAtom) {
        elem = atom;
      }
    }
  }

  // This could replace the atom pointer in compat->stereoGroups, if it's
  // not null, but because nothing should be reading while they're being
  // modified anyway, they can be cleared and recreated on read.
  // If this is ever a performance issue, some more extensive work may be
  // required, because StereoGroup doesn't have a function to replace atoms.
  compat->stereoGroups.store(nullptr, std::memory_order_relaxed);

  if (atom->hasQuery()) {
    auto *query = atom->getQuery();
    PRECONDITION(query != nullptr,
                 "Atom hasQuery returns true, but the query is null");
    setAtomQueryCompat(index, query);
  }
}

void RDMol::replaceBondPointerCompat(Bond *bond, uint32_t index) {
  PRECONDITION(bond->dp_owningMol == nullptr || bond->dp_owningMol == this,
               "Bond cannot have another owning mol");
  PRECONDITION(bond->dp_dataMol != this, "About to delete this RDMol");

  // First clean up previous data. We're only taking the pointer, so no need to
  // copy the data.
  delete bond->dp_dataMol;

  auto *compat = ensureCompatInit();

  // Store original pointer for bookmark replacement.
  auto &bond_ptr = compat->bonds[index];
  Bond *origBond = bond_ptr.get();

  // Now set the pointers for this RDMol.
  // The Bond pointed to by origBond is invalid after calling reset,
  // but only the pointer value is used below.
  bond_ptr.reset(bond);
  bond->setIdx(index);
  bond->dp_dataMol = this;
  bond->dp_owningMol = this;

  // Replace pointers in bookmarks.
  for (auto &ab : compat->bondBookmarks) {
    for (auto &elem : ab.second) {
      if (elem == origBond) {
        elem = bond;
      }
    }
  }

  // This could replace the bond pointer in compat->stereoGroups, if it's
  // not null, but because nothing should be reading while they're being
  // modified anyway, they can be cleared and recreated on read.
  // If this is ever a performance issue, some more extensive work may be
  // required, because StereoGroup doesn't have a function to replace bonds.
  compat->stereoGroups.store(nullptr, std::memory_order_relaxed);

  if (bond->hasQuery()) {
    auto *query = bond->getQuery();
    PRECONDITION(query != nullptr,
                 "Bond hasQuery returns true, but the query is null");
    setBondQueryCompat(index, query);
  }
}

template <bool toCompat, typename F>
void initCompatForWriteHelper(std::atomic<CompatSyncStatus> &syncStatus,
                              std::mutex &compatibilityMutex, F &&f) {
  const CompatSyncStatus fromStatus = toCompat
                                          ? CompatSyncStatus::lastUpdatedRDMol
                                          : CompatSyncStatus::lastUpdatedCompat;
  const CompatSyncStatus toStatus = toCompat
                                        ? CompatSyncStatus::lastUpdatedCompat
                                        : CompatSyncStatus::lastUpdatedRDMol;

  // Because the CPU treats reading the status an independent read from reading
  // the conformer data, memory_order_acquire must be used to prevent
  // speculative reads.
  auto status = syncStatus.load(std::memory_order_acquire);
  if (status != toStatus) {
    if (status == fromStatus) {
      std::lock_guard<std::mutex> lock_scope(compatibilityMutex);
      // Don't re-copy if another thread already copied while this thread was
      // waiting to acquire the lock.
      status = syncStatus.load(std::memory_order_acquire);
      if (status == fromStatus) {
        f();

        // memory_order_release must be used to prevent delayed writes of the
        // conformer data.
        syncStatus.store(toStatus, std::memory_order_release);
      } else if (status == CompatSyncStatus::inSync) {
        // Another thread must have requested read access at the same time.
        // Going from inSync to lastUpdatedCompat, conformer data has already
        // been written out to main memory, so can do a relaxed store.
        syncStatus.store(toStatus, std::memory_order_relaxed);
      }
    } else {
      // Already in sync. Conformer data has already been
      // written out to main memory, so can do a relaxed store.
      syncStatus.store(toStatus, std::memory_order_relaxed);
    }
  }
}

template <bool toCompat, typename F>
void initCompatForReadHelper(std::atomic<CompatSyncStatus> &syncStatus,
                             std::mutex &compatibilityMutex, F &&f) {
  const CompatSyncStatus fromStatus = toCompat
                                          ? CompatSyncStatus::lastUpdatedRDMol
                                          : CompatSyncStatus::lastUpdatedCompat;
  // Because the CPU treats reading the status an independent read from reading
  // the conformer data, memory_order_acquire must be used to prevent
  // speculative reads.
  auto status = syncStatus.load(std::memory_order_acquire);
  if (status == fromStatus) {
    std::lock_guard<std::mutex> lock_scope(compatibilityMutex);
    // Don't re-copy if another thread already copied while this thread was
    // waiting to acquire the lock.
    status = syncStatus.load(std::memory_order_acquire);
    if (status == fromStatus) {
      f();

      // Only requesting read access in this thread, so status can
      // be set to inSync. memory_order_release must be used to prevent
      // delayed writes of the conformer data.
      syncStatus.store(CompatSyncStatus::inSync, std::memory_order_release);
    }
  }
}

ROMol::CONF_SPTR_LIST &RDMol::getConformersCompat() {
  auto *const compat = ensureCompatInit();
  initCompatForWriteHelper<true>(
      compat->conformerSyncStatus, compatibilityMutex,
      [&]() { copyConformersToCompatibilityData(compat); });
  return compat->conformers;
}
const ROMol::CONF_SPTR_LIST &RDMol::getConformersCompat() const {
  auto *const compat = ensureCompatInit();
  initCompatForReadHelper<true>(
      compat->conformerSyncStatus, compatibilityMutex,
      [&]() { copyConformersToCompatibilityData(compat); });
  return compat->conformers;
}

RingInfo &RDMol::getRingInfoCompat() {
  auto *const compat = ensureCompatInit();
  initCompatForWriteHelper<true>(
      compat->ringInfoSyncStatus, compatibilityMutex, [&]() {
        copyRingInfoToCompatibilityData(ringInfo, *compat->ringInfo,
                                        getNumAtoms(), getNumBonds());
      });
  return *compat->ringInfo;
}
const RingInfo &RDMol::getRingInfoCompat() const {
  auto *const compat = ensureCompatInit();
  initCompatForReadHelper<true>(
      compat->ringInfoSyncStatus, compatibilityMutex, [&]() {
        copyRingInfoToCompatibilityData(ringInfo, *compat->ringInfo,
                                        getNumAtoms(), getNumBonds());
      });
  return *compat->ringInfo;
}

const std::vector<StereoGroup> &RDMol::getStereoGroupsCompat() {
  RDMol::CompatibilityData *compat = ensureCompatInit();

  // First, try without locking
  // We can do a relaxed load of the pointer, because accessing the other data
  // requires dereferencing the pointer, so would all be dependent reads.
  auto *stereoGroupsData = compat->stereoGroups.load(std::memory_order_relaxed);

  // Check if we need to repopulate (either nullptr or empty after clear)
  bool needsRepopulate =
      (stereoGroupsData == nullptr) || (stereoGroupsData->size() == 0);

  if (!needsRepopulate) {
    return *stereoGroupsData;
  }

  std::lock_guard<std::mutex> lock_scope(compatibilityMutex);
  // Check again after locking
  stereoGroupsData = compat->stereoGroups.load(std::memory_order_relaxed);
  needsRepopulate =
      (stereoGroupsData == nullptr) || (stereoGroupsData->size() == 0);

  if (!needsRepopulate) {
    return *stereoGroupsData;
  }

  // Create new vector if needed, otherwise reuse existing one
  if (stereoGroupsData == nullptr) {
    stereoGroupsData = new std::vector<StereoGroup>();
    compat->stereoGroups.store(stereoGroupsData, std::memory_order_relaxed);
  } else {
    // Clear existing vector to repopulate it
    stereoGroupsData->clear();
  }

  if (stereoGroups.get() != nullptr) {
    const std::vector<StereoGroupType> &types = stereoGroups->stereoTypes;
    const std::vector<uint32_t> &atomIndices = stereoGroups->atoms;
    const std::vector<uint32_t> &bondIndices = stereoGroups->bonds;
    const std::vector<uint32_t> &atomBegins = stereoGroups->atomBegins;
    const std::vector<uint32_t> &bondBegins = stereoGroups->bondBegins;

    const size_t numGroups = stereoGroups->getNumGroups();
    CHECK_INVARIANT(atomBegins.size() == numGroups + 1,
                    "Atom begins must be size 1 + number of groups");
    CHECK_INVARIANT(bondBegins.size() == numGroups + 1,
                    "Bond begins must be size 1 + number of groups");

    for (size_t i = 0; i < numGroups; i++) {
      std::vector<Atom *> atomPtrs;
      atomPtrs.reserve(atomBegins[i + 1] - atomBegins[i]);
      for (size_t j = atomBegins[i]; j < atomBegins[i + 1]; j++) {
        atomPtrs.push_back(compat->atoms[atomIndices[j]].get());
      }
      std::vector<Bond *> bondPtrs;
      bondPtrs.reserve(bondBegins[i + 1] - bondBegins[i]);
      for (size_t j = bondBegins[i]; j < bondBegins[i + 1]; j++) {
        bondPtrs.push_back(compat->bonds[bondIndices[j]].get());
      }
      stereoGroupsData->emplace_back(types[i], std::move(atomPtrs),
                                     std::move(bondPtrs),
                                     stereoGroups->readIds[i]);
      // Set writeId BEFORE setOwningMol to avoid marking as modified during
      // initialization
      stereoGroupsData->back().setWriteId(stereoGroups->writeIds[i]);
      // Now set the owner so future setWriteId calls will mark as modified
      stereoGroupsData->back().setOwningMol(compat->compatMol);
    }
  }

  return *stereoGroupsData;
}

void RDMol::clearStereoGroupsCompat() {
  // Nothing should be trying to read the stereo groups while it's being
  // written, so it's safe to clear the old compatibility data without locking.
  RDMol::CompatibilityData *compat =
      compatibilityData.load(std::memory_order_relaxed);
  if (compat != nullptr) {
    auto *compatStereoGroups =
        compat->stereoGroups.load(std::memory_order_relaxed);
    if (compatStereoGroups != nullptr) {
      // Clear the vector contents but keep the vector object alive to preserve
      // Python references. The vector will be repopulated on next access.
      compatStereoGroups->clear();
    }
  }
}

void RDMol::copyFromCompatibilityData(const CompatibilityData *source,
                                      bool quickCopy, int confId) {
  if (source == nullptr) {
    source = compatibilityData.load(std::memory_order_acquire);

    if (source == nullptr) {
      // Nothing to copy or clear
      return;
    }
  }

  // Copy back bookmarks. Bookmarks are all in CompatibilityData if it exists.
  atomBookmarks.clear();
  if (!quickCopy) {
    atomBookmarks.reserve(source->atomBookmarks.size());
    for (const auto &sourcePair : source->atomBookmarks) {
      std::vector<int> atomIndices;
      atomIndices.reserve(sourcePair.second.size());
      for (const auto *atom : sourcePair.second) {
        atomIndices.push_back(atom->getIdx());
      }
      atomBookmarks.insert(
          std::make_pair(sourcePair.first, std::move(atomIndices)));
    }
  }
  bondBookmarks.clear();
  if (!quickCopy) {
    bondBookmarks.reserve(source->bondBookmarks.size());
    for (const auto &sourcePair : source->bondBookmarks) {
      std::vector<int> bondIndices;
      bondIndices.reserve(sourcePair.second.size());
      for (const auto *bond : sourcePair.second) {
        bondIndices.push_back(bond->getIdx());
      }
      bondBookmarks.insert(
          std::make_pair(sourcePair.first, std::move(bondIndices)));
    }
  }

  // Copy back bond stereo atoms.  Values are only in CompatibilityData for
  // non-empty vectors. Any other values are expected to have been copied
  // separately, if source is not owned by this.
  PRECONDITION(
      source->bondStereoAtoms.size() == getNumBonds(),
      "CompatibilityData::bondStereoAtoms must have correct size for this RDMol");
  for (uint32_t bondIndex = 0, numBonds = getNumBonds(); bondIndex < numBonds;
       ++bondIndex) {
    const auto &atoms = *source->bondStereoAtoms[bondIndex];
    PRECONDITION(atoms.size() == 0 || atoms.size() == 2,
                 "bondStereoAtoms must have either no atoms or 2 atoms");
    if (atoms.size() == 2) {
      BondData &bond = getBond(bondIndex);
      bond.stereoAtoms[0] = atoms[0];
      bond.stereoAtoms[1] = atoms[1];
    }
  }

  // memory_order_acquire here is probably excessive, because no thread should
  // be synchronizing conformers at the same time, but play it safe.
  if (!quickCopy &&
      source->conformerSyncStatus.load(std::memory_order_acquire) ==
          CompatSyncStatus::lastUpdatedCompat) {
    copyConformersFromCompatibilityData(source, confId);
  }

  if (source->ringInfoSyncStatus.load(std::memory_order_acquire) ==
      CompatSyncStatus::lastUpdatedCompat) {
    copyRingInfoFromCompatibilityData(*source->ringInfo, ringInfo,
                                      getNumAtoms(), getNumBonds());
  }

  // Sync stereo groups from compat layer if they were modified
  if (source->stereoGroupsSyncStatus.load(std::memory_order_acquire) ==
      CompatSyncStatus::lastUpdatedCompat) {
    auto *compatStereoGroups =
        source->stereoGroups.load(std::memory_order_acquire);
    if (compatStereoGroups != nullptr && stereoGroups != nullptr &&
        compatStereoGroups->size() == stereoGroups->getNumGroups()) {
      // Sync writeIds from compat layer back to flat structure
      for (size_t i = 0; i < compatStereoGroups->size(); ++i) {
        stereoGroups->writeIds[i] = (*compatStereoGroups)[i].getWriteId();
      }
    }
  }

  delete compatibilityData.load(std::memory_order_relaxed);
  compatibilityData.store(nullptr, std::memory_order_relaxed);
}

RDMol::RDMol(ROMol *existingPtr) : RDMol() {
  auto *dataCopy = new CompatibilityData(*this, existingPtr);
  compatibilityData.store(dataCopy, std::memory_order_relaxed);
}

RDMol::RDMol(const RDMol &other, bool quickCopy, int confId,
             ROMol *existingPtr) {
  initFromOther(other, quickCopy, confId, existingPtr);
}

RDMol::RDMol(const RDMol &other, bool quickCopy, int confId)
    : RDMol(other, quickCopy, confId, nullptr) {}

RDMol::RDMol(const RDMol &other, ROMol *existingPtr)
    : RDMol(other, false, -1, existingPtr) {}

RDMol::RDMol(const RDMol &other) : RDMol(other, false, -1, nullptr) {}

RDMol &RDMol::operator=(const RDMol &other) {
  if (this == &other) {
    return *this;
  }
  atomData = other.atomData;
  atomBondStarts = other.atomBondStarts;
  otherAtomIndices = other.otherAtomIndices;
  bondDataIndices = other.bondDataIndices;
  bondData = other.bondData;
  ringInfo = other.ringInfo;
  numConformers = other.numConformers;
  conformerAtomCapacity = other.conformerAtomCapacity;
  conformerPositionData = other.conformerPositionData;
  conformerIdsAndFlags = other.conformerIdsAndFlags;
  substanceGroups = other.substanceGroups;
  isStereochemDone = other.isStereochemDone;
  properties = other.properties;
  atomBookmarks = other.atomBookmarks;
  bondBookmarks = other.bondBookmarks;

  if (other.dp_delAtoms != nullptr) {
    dp_delAtoms = std::make_unique<boost::dynamic_bitset<>>(*other.dp_delAtoms);
    dp_delBonds = std::make_unique<boost::dynamic_bitset<>>(*other.dp_delBonds);
  }

  for (const auto &[key, value] : other.monomerInfo) {
    monomerInfo[key] = std::unique_ptr<AtomMonomerInfo>(value->copy());
  }

  // Sync stereo groups from compat layer before copying if they were modified
  if (other.hasCompatibilityData()) {
    auto *otherCompat = other.getCompatibilityDataIfPresent();
    if (otherCompat != nullptr &&
        otherCompat->stereoGroupsSyncStatus.load(std::memory_order_acquire) ==
            CompatSyncStatus::lastUpdatedCompat) {
      auto *compatStereoGroups =
          otherCompat->stereoGroups.load(std::memory_order_acquire);
      if (compatStereoGroups != nullptr && other.stereoGroups != nullptr &&
          compatStereoGroups->size() == other.stereoGroups->getNumGroups()) {
        // Sync writeIds from compat layer to flat structure before copying
        for (size_t i = 0; i < compatStereoGroups->size(); ++i) {
          const_cast<RDMol &>(other).stereoGroups->writeIds[i] =
              (*compatStereoGroups)[i].getWriteId();
        }
      }
    }
  }

  if (other.stereoGroups != nullptr) {
    stereoGroups = std::make_unique<StereoGroups>(*other.stereoGroups);
  }

  delete compatibilityData.load(std::memory_order_relaxed);
  compatibilityData.store(nullptr, std::memory_order_relaxed);
  if (other.hasCompatibilityData()) {
    copyFromCompatibilityData(other.getCompatibilityDataIfPresent());
  }

  for (auto &sg : substanceGroups) {
    sg.setOwningMol(&asROMol());
  }

  return *this;
}

RDMol::RDMol(RDMol &&other) noexcept
    : atomData(std::move(other.atomData)),
      atomBondStarts(std::move(other.atomBondStarts)),
      otherAtomIndices(std::move(other.otherAtomIndices)),
      bondDataIndices(std::move(other.bondDataIndices)),
      bondData(std::move(other.bondData)),
      stereoGroups(std::move(other.stereoGroups)),
      ringInfo(std::move(other.ringInfo)),
      substanceGroups(std::move(other.substanceGroups)),
      numConformers(other.numConformers),
      conformerAtomCapacity(other.conformerAtomCapacity),
      conformerPositionData(std::move(other.conformerPositionData)),
      conformerIdsAndFlags(std::move(other.conformerIdsAndFlags)),
      isStereochemDone(std::move(other.isStereochemDone)),
      dp_delAtoms(std::move(other.dp_delAtoms)),
      dp_delBonds(std::move(other.dp_delBonds)),
      properties(std::move(other.properties)),
      atomBookmarks(std::move(other.atomBookmarks)),
      bondBookmarks(std::move(other.bondBookmarks)),
      monomerInfo(std::move(other.monomerInfo)) {
  CompatibilityData *data =
      other.compatibilityData.load(std::memory_order_relaxed);
  if (data != nullptr) {
    data->setNewOwner(*this);
  }
  compatibilityData.store(data, std::memory_order_relaxed);
  other.compatibilityData.store(nullptr, std::memory_order_relaxed);
}

RDMol &RDMol::operator=(RDMol &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  atomData = std::move(other.atomData);
  atomBondStarts = std::move(other.atomBondStarts);
  otherAtomIndices = std::move(other.otherAtomIndices);
  bondDataIndices = std::move(other.bondDataIndices);
  bondData = std::move(other.bondData);
  stereoGroups = std::move(other.stereoGroups);
  ringInfo = std::move(other.ringInfo);
  substanceGroups = std::move(other.substanceGroups);
  numConformers = other.numConformers;
  conformerAtomCapacity = other.conformerAtomCapacity;
  conformerPositionData = std::move(other.conformerPositionData);
  conformerIdsAndFlags = std::move(other.conformerIdsAndFlags);
  isStereochemDone = std::move(other.isStereochemDone);
  dp_delAtoms = std::move(other.dp_delAtoms);
  dp_delBonds = std::move(other.dp_delBonds);
  properties = std::move(other.properties);
  atomBookmarks = std::move(other.atomBookmarks);
  bondBookmarks = std::move(other.bondBookmarks);
  monomerInfo = std::move(other.monomerInfo);

  delete compatibilityData.load(std::memory_order_relaxed);
  CompatibilityData *data =
      other.compatibilityData.load(std::memory_order_relaxed);
  if (data != nullptr) {
    data->setNewOwner(*this);
  }
  compatibilityData.store(data, std::memory_order_relaxed);
  other.compatibilityData.store(nullptr, std::memory_order_relaxed);

  return *this;
}

void RDMol::initFromOther(const RDMol &other, bool quickCopy, int confId,
                          ROMol *existingPtr) {
  atomData = other.atomData;
  atomBondStarts = other.atomBondStarts;
  otherAtomIndices = other.otherAtomIndices;
  bondDataIndices = other.bondDataIndices;
  bondData = other.bondData;
  ringInfo = other.ringInfo;
  if (other.atomQueries.size() != 0) {
    atomQueries.reserve(other.atomQueries.size());
    for (auto &&query : other.atomQueries) {
      atomQueries.emplace_back((query != nullptr) ? query->copy() : nullptr);
    }
  }
  if (other.bondQueries.size() != 0) {
    bondQueries.reserve(other.bondQueries.size());
    for (auto &&query : other.bondQueries) {
      bondQueries.emplace_back((query != nullptr) ? query->copy() : nullptr);
    }
  }

  // Sync stereo groups from compat layer before copying if they were modified
  if (other.hasCompatibilityData()) {
    auto *otherCompat = other.getCompatibilityDataIfPresent();
    if (otherCompat != nullptr &&
        otherCompat->stereoGroupsSyncStatus.load(std::memory_order_acquire) ==
            CompatSyncStatus::lastUpdatedCompat) {
      auto *compatStereoGroups =
          otherCompat->stereoGroups.load(std::memory_order_acquire);
      if (compatStereoGroups != nullptr && other.stereoGroups != nullptr &&
          compatStereoGroups->size() == other.stereoGroups->getNumGroups()) {
        // Sync writeIds from compat layer to flat structure before copying
        for (size_t i = 0; i < compatStereoGroups->size(); ++i) {
          const_cast<RDMol &>(other).stereoGroups->writeIds[i] =
              (*compatStereoGroups)[i].getWriteId();
        }
      }
    }
  }

  if (other.stereoGroups != nullptr) {
    stereoGroups = std::make_unique<StereoGroups>(*other.stereoGroups);
  }

  if (other.dp_delAtoms != nullptr) {
    dp_delAtoms = std::make_unique<boost::dynamic_bitset<>>(*other.dp_delAtoms);
    dp_delBonds = std::make_unique<boost::dynamic_bitset<>>(*other.dp_delBonds);
  }

  for (const auto &[key, value] : other.monomerInfo) {
    monomerInfo[key] = std::unique_ptr<AtomMonomerInfo>(value->copy());
  }

  // This is currently only called in cases of construction or with quickCopy
  // being false, so it's unclear what should happen to existing bookmarks and
  // conformers if quickCopy is true, though the previous
  // ROMol::initFromOther code did explicitly clear props.
  PRECONDITION(
      !quickCopy || (atomBookmarks.size() == 0 && bondBookmarks.size() == 0 &&
                     numConformers == 0),
      "RDMol::initFromOther wasn't designed for quickCopy over existing data");
  if (quickCopy) {
    properties.clear();
  } else {
    properties = other.properties;
    atomBookmarks = other.atomBookmarks;
    bondBookmarks = other.bondBookmarks;

    if (confId < 0) {
      // Copy all conformers
      conformerIdsAndFlags = other.conformerIdsAndFlags;
      conformerPositionData = other.conformerPositionData;
      numConformers = other.numConformers;
      conformerAtomCapacity = other.conformerAtomCapacity;
    } else {
      uint32_t index = other.numConformers;
      try {
        index = other.findConformerIndex(confId);
      } catch (const ConformerException &) {
        // Noop, matching ROMol behavior.
      }
      if (index != other.numConformers) {
        allocateConformers(1, other.getNumAtoms());
        const double *otherData = other.conformerPositionData.data() +
                                  index * other.conformerAtomCapacity * 3;
        std::copy(otherData, otherData + other.getNumAtoms() * 3,
                  conformerPositionData.begin());
        if (other.conformerIdsAndFlags.size() > 0 &&
            (confId != 0 || !other.conformerIdsAndFlags[index].is3D)) {
          conformerIdsAndFlags.push_back(other.conformerIdsAndFlags[index]);
        }
        numConformers = 1;
      }
    }
  }

  const CompatibilityData *otherCompat = other.getCompatibilityDataIfPresent();
  const bool otherHasCompat = otherCompat != nullptr;
  if (otherHasCompat) {
    // The behaviour of quickCopy and confId above must match the behaviour
    // in copyFromCompatibilityData
    copyFromCompatibilityData(otherCompat, quickCopy, confId);
  }

  if (existingPtr) {
    PRECONDITION(
        !hasCompatibilityData(),
        "Cannot create RDMol with existing ROMol pointer and compatibility data");
    CompatibilityData *data = new CompatibilityData(*this, existingPtr);

    // If the other molecule had CompatibilityData with conformers containing
    // properties, copy those Conformer objects (including properties) instead
    // of using the newly created empty ones
    if (otherHasCompat && !quickCopy && confId < 0 && numConformers > 0 &&
        otherCompat->conformers.size() == numConformers) {
      data->copyConformersFrom(otherCompat);
    }

    existingPtr->dp_mol = this;
    // Ensure the compatibility data is written out to main memory before the
    // pointer to it (memory_order_release), in case another thread reads it
    // before the cache is flushed.
    compatibilityData.store(data, std::memory_order_release);
  }

  if (!quickCopy && other.substanceGroups.size() > 0) {
    // If substance groups exist, instantiate after we have an ROMol pointer.
    substanceGroups = other.substanceGroups;
    ROMol *thisROMol = &asROMol();
    for (auto &sg : substanceGroups) {
      sg.setOwningMol(thisROMol);
    }
  }
}

PropArray::PropArray() {
  size = 0;
  family = PropertyType::CHAR;
  data = nullptr;
}

PropArray::PropArray(uint32_t size, PropertyType family, bool isSet)
    : size(size), family(family) {
  construct(isSet);
}

namespace {

constexpr bool is8BitType(PropertyType family) {
  static_assert(sizeof(bool) == sizeof(char),
                "size assumption bool==char violated");
  return family == PropertyType::CHAR || family == PropertyType::BOOL;
}

constexpr bool is32BitType(PropertyType family) {
  static_assert(sizeof(float) == sizeof(int),
                "size assumption float==int violated");
  return family == PropertyType::INT32 || family == PropertyType::UINT32 ||
         family == PropertyType::FLOAT;
}

constexpr bool is64BitType(PropertyType family) {
  static_assert(sizeof(double) == sizeof(int64_t),
                "size assumption double==int64_t violated");
  return family == PropertyType::INT64 || family == PropertyType::UINT64 ||
         family == PropertyType::DOUBLE;
}

}  // namespace

PropArray::PropArray(const PropArray &other)
    : PropArray(other.size, other.family, false) {
  if (is8BitType(family)) {
    std::copy(static_cast<char *>(other.data),
              static_cast<char *>(other.data) + size,
              static_cast<char *>(data));
  } else if (is32BitType(family)) {
    std::copy(static_cast<int32_t *>(other.data),
              static_cast<int32_t *>(other.data) + size,
              static_cast<int32_t *>(data));
  } else if (is64BitType(family)) {
    std::copy(static_cast<int64_t *>(other.data),
              static_cast<int64_t *>(other.data) + size,
              static_cast<int64_t *>(data));
  } else {
    for (uint32_t i = 0; i < size; ++i) {
      copy_rdvalue(static_cast<RDValue *>(data)[i],
                   static_cast<RDValue *>(other.data)[i]);
    }
  }
  std::copy(other.isSetMask.get(), other.isSetMask.get() + size,
            isSetMask.get());
  numSet = other.numSet;
};

PropArray &PropArray::operator=(const PropArray &other) {
  if (this == &other) {
    return *this;
  }
  if (size > 0) {
    destroy();
  }
  size = other.size;
  family = other.family;
  construct(/*isSet=*/false);
  isSetMask = std::make_unique<bool[]>(size);
  std::copy(other.isSetMask.get(), other.isSetMask.get() + size,
            isSetMask.get());
  numSet = other.numSet;
  return *this;
}

PropArray::PropArray(PropArray &&other) {
  size = other.size;
  family = other.family;
  data = other.data;
  isSetMask = std::move(other.isSetMask);
  numSet = other.numSet;
  other.data = nullptr;
  other.size = 0;
}

PropArray &PropArray::operator=(PropArray &&other) {
  if (this == &other) {
    return *this;
  }
  if (size > 0) {
    destroy();
  }
  size = other.size;
  family = other.family;
  data = other.data;
  isSetMask = std::move(other.isSetMask);
  numSet = other.numSet;
  other.data = nullptr;
  other.size = 0;
  return *this;
}

PropArray::~PropArray() noexcept { destroy(); }

void PropArray::construct(bool isSet) {
  PRECONDITION(data == nullptr, "Constructing on existing data");
  isSetMask = std::make_unique<bool[]>(size);
  std::fill(isSetMask.get(), isSetMask.get() + size, isSet);
  numSet = isSet ? size : 0;
  if (is8BitType(family)) {
    data = new char[size];
  } else if (is32BitType(family)) {
    data = new int32_t[size];
  } else if (is64BitType(family)) {
    data = new int64_t[size];
  } else {
    // These are default constructed
    data = new RDValue[size];
  }
}

void PropArray::destroy() {
  if (is8BitType(family)) {
    delete[] static_cast<char *>(data);
  } else if (is32BitType(family)) {
    delete[] static_cast<int32_t *>(data);
  } else if (is64BitType(family)) {
    delete[] static_cast<double *>(data);
  } else {
    for (uint32_t i = 0; i < size; ++i) {
      RDValue::cleanup_rdvalue(static_cast<RDValue *>(data)[i]);
    }
    delete[] static_cast<RDValue *>(data);
  }
}

RDValue PropArray::toRDValue(uint32_t idx) const {
  PRECONDITION(data != nullptr, "Accessing null prop array");
  switch (family) {
    case PropertyType::CHAR:
      return RDValue(static_cast<char *>(data)[idx]);
    case PropertyType::BOOL:
      return RDValue(static_cast<bool *>(data)[idx]);
    case PropertyType::INT32:
      return RDValue(static_cast<int32_t *>(data)[idx]);
    case PropertyType::UINT32:
      return RDValue(static_cast<uint32_t *>(data)[idx]);
    case PropertyType::FLOAT:
      return RDValue(static_cast<float *>(data)[idx]);
    case PropertyType::INT64:
      return RDValue(static_cast<int64_t *>(data)[idx]);
    case PropertyType::UINT64:
      return RDValue(static_cast<uint64_t *>(data)[idx]);
    case PropertyType::DOUBLE:
      return RDValue(static_cast<double *>(data)[idx]);
    default:
      RDValue res;
      copy_rdvalue(res, static_cast<RDValue *>(data)[idx]);
      return res;
  }
}

void PropArray::appendElement() {
  ++size;
  std::unique_ptr<bool[]> newMask = std::make_unique<bool[]>(size);
  std::copy(isSetMask.get(), isSetMask.get() + size - 1, newMask.get());
  newMask[size - 1] = false;
  isSetMask = std::move(newMask);
  if (is8BitType(family)) {
    auto *newData = new char[size];
    std::copy(static_cast<char *>(data), static_cast<char *>(data) + size - 1,
              newData);
    delete[] static_cast<char *>(data);
    data = newData;
  } else if (is32BitType(family)) {
    auto *newData = new int32_t[size];
    std::copy(static_cast<int32_t *>(data),
              static_cast<int32_t *>(data) + size - 1, newData);
    delete[] static_cast<int32_t *>(data);
    data = newData;
  } else if (is64BitType(family)) {
    auto *newData = new int64_t[size];
    std::copy(static_cast<int64_t *>(data),
              static_cast<int64_t *>(data) + size - 1, newData);
    delete[] static_cast<int64_t *>(data);
    data = newData;
  } else {
    auto *newData = new RDValue[size];
    for (uint32_t i = 0; i < size - 1; ++i) {
      RDValue *orig = static_cast<RDValue *>(data) + i;
      // TODO: A move should be possible here.
      copy_rdvalue(newData[i], *orig);
      RDValue::cleanup_rdvalue(*orig);
    }
    newData[size - 1] = RDValue();
    delete[] static_cast<RDValue *>(data);
    data = newData;
  }
}

void PropArray::removeElement(uint32_t index) {
  URANGE_CHECK(index, size);

  // Shift down isSetMask. Note that isSetMask is now at least one element
  // larger than data, but this is never a problem. If elements are added
  // instead of removed, the mask is recreated.
  std::copy(isSetMask.get() + index + 1, isSetMask.get() + size,
            isSetMask.get() + index);

  if (is8BitType(family)) {
    std::copy(static_cast<char *>(data) + index + 1,
              static_cast<char *>(data) + size,
              static_cast<char *>(data) + index);
  } else if (is32BitType(family)) {
    std::copy(static_cast<int32_t *>(data) + index + 1,
              static_cast<int32_t *>(data) + size,
              static_cast<int32_t *>(data) + index);
  } else if (is64BitType(family)) {
    std::copy(static_cast<int64_t *>(data) + index + 1,
              static_cast<int64_t *>(data) + size,
              static_cast<int64_t *>(data) + index);
  } else {
    auto *dataPtr = static_cast<RDValue *>(data);
    for (uint32_t i = index + 1; i < size; ++i) {
      copy_rdvalue(dataPtr[i - 1], dataPtr[i]);
    }
    RDValue::cleanup_rdvalue(dataPtr[size - 1]);
  }
  --size;
}

namespace {
template <typename T>
RDValue *convertToRDValueHelper(const T *data, bool *isSetMask, uint32_t size) {
  RDValue *output = new RDValue[size];
  for (uint32_t i = 0; i < size; ++i) {
    if (isSetMask[i]) {
      output[i] = RDValue(data[i]);
    }
  }
  return output;
}
}  // namespace

void PropArray::convertToRDValue() {
  if (family == PropertyType::ANY) {
    return;
  }
  PRECONDITION(data != nullptr, "Accessing null prop array");
  // This null value will never be used. It's to avoid uninitialized value
  // warnings.
  RDValue *newData = nullptr;
  switch (family) {
    case PropertyType::CHAR:
      newData = convertToRDValueHelper(static_cast<char *>(data),
                                       isSetMask.get(), size);
      break;
    case PropertyType::BOOL:
      newData = convertToRDValueHelper(static_cast<bool *>(data),
                                       isSetMask.get(), size);
      break;
    case PropertyType::INT32:
      newData = convertToRDValueHelper(static_cast<int32_t *>(data),
                                       isSetMask.get(), size);
      break;
    case PropertyType::UINT32:
      newData = convertToRDValueHelper(static_cast<uint32_t *>(data),
                                       isSetMask.get(), size);
      break;
    case PropertyType::FLOAT:
      newData = convertToRDValueHelper(static_cast<float *>(data),
                                       isSetMask.get(), size);
      break;
    case PropertyType::INT64:
      newData = convertToRDValueHelper(static_cast<int64_t *>(data),
                                       isSetMask.get(), size);
      break;
    case PropertyType::UINT64:
      newData = convertToRDValueHelper(static_cast<uint64_t *>(data),
                                       isSetMask.get(), size);
      break;
    case PropertyType::DOUBLE:
      newData = convertToRDValueHelper(static_cast<double *>(data),
                                       isSetMask.get(), size);
      break;
    default:
      raiseNonImplementedDetail(
          "Unsupported type in PropArray::convertToRDValue");
  }
  destroy();
  data = newData;
  family = PropertyType::ANY;
}

void RingInfoCache::initFusedInfoFromBondMemberships() {
  // Call clear first to ensure all elements reinitialized to false or 0
  const size_t numRings = RingInfoCache::numRings();
  areRingsFused.clear();
  areRingsFused.resize(numRings * numRings, false);
  numFusedBonds.clear();
  numFusedBonds.resize(numRings, 0);
  if (numRings <= 1) {
    return;
  }

  // Create 2D bit vector of whether two rings share a bond
  const size_t numBonds = bondMembershipBegins.size() - 1;
  for (size_t bondIndex = 0; bondIndex < numBonds; ++bondIndex) {
    uint32_t begin = bondMembershipBegins[bondIndex];
    uint32_t end = bondMembershipBegins[bondIndex + 1];
    for (; begin + 1 < end; ++begin) {
      auto ring1 = bondMemberships[begin];
      for (uint32_t other = begin + 1; other < end; ++other) {
        auto ring2 = bondMemberships[other];
        areRingsFused[ring1 * numRings + ring2] = true;
        areRingsFused[ring2 * numRings + ring1] = true;
        ++numFusedBonds[ring1];
        ++numFusedBonds[ring2];
      }
    }
  }
}

void StereoGroups::addGroup(StereoGroupType type,
                            const std::vector<uint32_t> &atomIndices,
                            const std::vector<uint32_t> &bondIndices,
                            uint32_t readId) {
  const uint32_t currentNumGroups = stereoTypes.size();
  stereoTypes.push_back(type);
  readIds.push_back(readId);
  writeIds.push_back(undefinedGroupId);
  const int atomBegin = atomBegins[currentNumGroups];
  const int bondBegin = bondBegins[currentNumGroups];
  atomBegins.push_back(atomBegin + atomIndices.size());
  bondBegins.push_back(bondBegin + bondIndices.size());
  atoms.insert(atoms.end(), atomIndices.begin(), atomIndices.end());
  bonds.insert(bonds.end(), bondIndices.begin(), bondIndices.end());
}

bool StereoGroups::GroupsEq(uint32_t index1, uint32_t index2) {
  URANGE_CHECK(index1, stereoTypes.size());
  URANGE_CHECK(index2, stereoTypes.size());
  if (stereoTypes[index1] != stereoTypes[index2]) {
    return false;
  }
  // Note read and write IDs not part of comparison

  const uint32_t atomBegin1 = atomBegins[index1];
  const uint32_t atomEnd1 = atomBegins[index1 + 1];
  const uint32_t atomBegin2 = atomBegins[index2];
  const uint32_t atomEnd2 = atomBegins[index2 + 1];
  if (atomEnd1 - atomBegin1 != atomEnd2 - atomBegin2) {
    return false;
  }

  const uint32_t bondBegin1 = bondBegins[index1];
  const uint32_t bondEnd1 = bondBegins[index1 + 1];
  const uint32_t bondBegin2 = bondBegins[index2];
  const uint32_t bondEnd2 = bondBegins[index2 + 1];
  if (bondEnd1 - bondBegin1 != bondEnd2 - bondBegin2) {
    return false;
  }

  return std::equal(atoms.begin() + atomBegin1, atoms.begin() + atomEnd1,
                    atoms.begin() + atomBegin2) &&
         std::equal(bonds.begin() + bondBegin1, bonds.begin() + bondEnd1,
                    bonds.begin() + bondBegin2);
}

void StereoGroups::removeGroup(uint32_t groupIndex) {
  URANGE_CHECK(groupIndex, stereoTypes.size());
  const uint32_t atomBegin = atomBegins[groupIndex];
  const uint32_t atomEnd = atomBegins[groupIndex + 1];
  const uint32_t bondBegin = bondBegins[groupIndex];
  const uint32_t bondEnd = bondBegins[groupIndex + 1];

  const uint32_t numAtomsErased = atomEnd - atomBegin;
  const uint32_t numBondsErased = bondEnd - bondBegin;

  stereoTypes.erase(stereoTypes.begin() + groupIndex);
  readIds.erase(readIds.begin() + groupIndex);
  writeIds.erase(writeIds.begin() + groupIndex);
  atomBegins.erase(atomBegins.begin() + groupIndex);
  bondBegins.erase(bondBegins.begin() + groupIndex);

  atoms.erase(atoms.begin() + atomBegin, atoms.begin() + atomEnd);
  bonds.erase(bonds.begin() + bondBegin, bonds.begin() + bondEnd);

  for (uint32_t j = groupIndex; j < atomBegins.size(); ++j) {
    atomBegins[j] -= numAtomsErased;
  }
  for (uint32_t j = groupIndex; j < bondBegins.size(); ++j) {
    bondBegins[j] -= numBondsErased;
  }
}

namespace {
template <typename T>
void removeGroups(StereoGroups &groups, const T &keepPredicate) {
  // Do a single pass to remove any groups, to avoid n^2 time
  const uint32_t numGroups = groups.stereoTypes.size();
  uint32_t runningNumGroups = 0;
  uint32_t runningNumAtoms = 0;
  uint32_t runningNumBonds = 0;
  for (uint32_t i = 0; i < numGroups; i++) {
    if (!keepPredicate(groups, i)) {
      continue;
    }

    uint32_t beginNumAtoms = runningNumAtoms;
    uint32_t beginNumBonds = runningNumBonds;

    // Group is kept, so copy it to its new location
    if (runningNumAtoms != groups.atomBegins[i]) {
      for (uint32_t j = groups.atomBegins[i], end = groups.atomBegins[i + 1];
           j < end; ++j) {
        groups.atoms[runningNumAtoms] = groups.atoms[j];
        ++runningNumAtoms;
      }
    } else {
      runningNumAtoms = groups.atomBegins[i + 1];
    }

    if (runningNumBonds != groups.bondBegins[i]) {
      for (uint32_t j = groups.bondBegins[i], end = groups.bondBegins[i + 1];
           j < end; ++j) {
        groups.bonds[runningNumBonds] = groups.bonds[j];
        ++runningNumBonds;
      }
    } else {
      runningNumBonds = groups.bondBegins[i + 1];
    }

    groups.atomBegins[runningNumGroups] = beginNumAtoms;
    groups.bondBegins[runningNumGroups] = beginNumBonds;
    groups.stereoTypes[runningNumGroups] = groups.stereoTypes[i];
    groups.readIds[runningNumGroups] = groups.readIds[i];
    groups.writeIds[runningNumGroups] = groups.writeIds[i];

    ++runningNumGroups;
  }

  groups.atoms.resize(runningNumAtoms);
  groups.bonds.resize(runningNumBonds);

  groups.atomBegins[runningNumGroups] = runningNumAtoms;
  groups.bondBegins[runningNumGroups] = runningNumBonds;
  groups.atomBegins.resize(runningNumGroups + 1);
  groups.bondBegins.resize(runningNumGroups + 1);
  groups.stereoTypes.resize(runningNumGroups);
  groups.readIds.resize(runningNumGroups);
  groups.writeIds.resize(runningNumGroups);
}
}  // namespace

void StereoGroups::removeAtomFromGroups(uint32_t atomIndex,
                                        bool decrementIndices) {
  const uint32_t numGroups = stereoTypes.size();
  uint32_t runningNumAtoms = 0;
  for (uint32_t i = 0; i < numGroups; i++) {
    const uint32_t atomBegin = atomBegins[i];
    const uint32_t atomEnd = atomBegins[i + 1];

    // set atomBegins again now that we've stored the previous value.
    atomBegins[i] = runningNumAtoms;

    for (uint32_t j = atomBegin; j < atomEnd; ++j) {
      if (atoms[j] != atomIndex) {
        atoms[runningNumAtoms] = atoms[j];
        if (decrementIndices && atoms[j] > atomIndex) {
          atoms[runningNumAtoms]--;
        }
        ++runningNumAtoms;
      }
    }
  }
  atomBegins[numGroups] = runningNumAtoms;
  atoms.resize(runningNumAtoms);

  // Remove groups that have both no atoms and no bonds
  // Equivalently, keep groups that have at least one atom or bond
  removeGroups(*this, [](const StereoGroups &groups, uint32_t i) {
    return groups.atomBegins[i] != groups.atomBegins[i + 1] ||
           groups.bondBegins[i] != groups.bondBegins[i + 1];
  });
}

void StereoGroups::removeBondFromGroups(uint32_t bondIndex,
                                        bool decrementIndices) {
  const uint32_t numGroups = stereoTypes.size();
  uint32_t runningNumBonds = 0;
  for (uint32_t i = 0; i < numGroups; i++) {
    const uint32_t bondBegin = bondBegins[i];
    const uint32_t bondEnd = bondBegins[i + 1];

    // set atomBegins again now that we've stored the previous value.
    bondBegins[i] = runningNumBonds;

    for (uint32_t j = bondBegin; j < bondEnd; ++j) {
      if (bonds[j] != bondIndex) {
        bonds[runningNumBonds] = bonds[j];
        if (decrementIndices && bonds[j] > bondIndex) {
          bonds[runningNumBonds]--;
        }
        ++runningNumBonds;
      }
    }
  }
  bondBegins[numGroups] = runningNumBonds;
  bonds.resize(runningNumBonds);

  // Remove groups that have both no atoms and no bonds
  // Equivalently, keep groups that have at least one atom or bond
  removeGroups(*this, [](const StereoGroups &groups, uint32_t i) {
    return groups.atomBegins[i] != groups.atomBegins[i + 1] ||
           groups.bondBegins[i] != groups.bondBegins[i + 1];
  });
}

void StereoGroups::removeGroupsWithAtom(const uint32_t atomIndex) {
  removeGroups(*this, [atomIndex](const StereoGroups &groups, uint32_t i) {
    const uint32_t atomBegin = groups.atomBegins[i];
    const uint32_t atomEnd = groups.atomBegins[i + 1];
    for (uint32_t j = atomBegin; j < atomEnd; ++j) {
      if (groups.atoms[j] == atomIndex) {
        return false;
      }
    }
    return true;
  });
}

void StereoGroups::removeGroupsWithAtoms(
    const std::vector<uint32_t> &atomIndices) {
  removeGroups(*this, [&atomIndices](const StereoGroups &groups, uint32_t i) {
    const uint32_t atomBegin = groups.atomBegins[i];
    const uint32_t atomEnd = groups.atomBegins[i + 1];
    for (uint32_t j = atomBegin; j < atomEnd; ++j) {
      if (std::find(atomIndices.begin(), atomIndices.end(), groups.atoms[j]) !=
          atomIndices.end()) {
        return false;
      }
    }
    return true;
  });
}

void StereoGroups::removeGroupsWithBond(const uint32_t bondIndex) {
  removeGroups(*this, [bondIndex](const StereoGroups &groups, uint32_t i) {
    const uint32_t bondBegin = groups.bondBegins[i];
    const uint32_t bondEnd = groups.bondBegins[i + 1];
    for (uint32_t j = bondBegin; j < bondEnd; ++j) {
      if (groups.bonds[j] == bondIndex) {
        return false;
      }
    }
    return true;
  });
}

void StereoGroups::removeGroupsWithBonds(
    const std::vector<uint32_t> &bondIndices) {
  removeGroups(*this, [&bondIndices](const StereoGroups &groups, uint32_t i) {
    const uint32_t bondBegin = groups.bondBegins[i];
    const uint32_t bondEnd = groups.bondBegins[i + 1];
    for (uint32_t j = bondBegin; j < bondEnd; ++j) {
      if (std::find(bondIndices.begin(), bondIndices.end(), groups.bonds[j]) !=
          bondIndices.end()) {
        return false;
      }
    }
    return true;
  });
}

// Stero group ID assignment adapted directly from StereoGroup.cpp
namespace {
void storeIdsInUse(boost::dynamic_bitset<> &ids, unsigned int &groupId) {
  if (groupId == StereoGroups::undefinedGroupId) {
    return;
  } else if (groupId >= ids.size()) {
    ids.resize(groupId + 1);
  }
  if (ids[groupId]) {
    // This id is duplicate, let's reset it so we can reassign it later
    BOOST_LOG(rdWarningLog)
        << "StereoGroup ID " << groupId
        << " is used by more than one group, and will be reassigned"
        << std::endl;
    groupId = 0;
  } else {
    ids[groupId] = true;
  }
};

void assignMissingIds(const boost::dynamic_bitset<> &ids, unsigned &nextId,
                      unsigned &groupId) {
  if (groupId == StereoGroups::undefinedGroupId) {
    ++nextId;
    while (nextId < ids.size() && ids[nextId]) {
      ++nextId;
    }
    groupId = nextId;
  }
};
}  // namespace

void StereoGroups::assignStereoGroupIds() {
  if (getNumGroups() == 0) {
    return;
  }

  boost::dynamic_bitset<> andIds;
  boost::dynamic_bitset<> orIds;

  for (uint32_t i = 0; i < stereoTypes.size(); i++) {
    const StereoGroupType groupType = stereoTypes[i];
    if (groupType == StereoGroupType::STEREO_AND) {
      storeIdsInUse(andIds, writeIds[i]);
    } else if (groupType == StereoGroupType::STEREO_OR) {
      storeIdsInUse(orIds, writeIds[i]);
    }
  }

  unsigned andId = undefinedGroupId;
  unsigned orId = undefinedGroupId;
  for (uint32_t i = 0; i < stereoTypes.size(); i++) {
    const StereoGroupType groupType = stereoTypes[i];
    if (groupType == StereoGroupType::STEREO_AND) {
      assignMissingIds(andIds, andId, writeIds[i]);
    } else if (groupType == StereoGroupType::STEREO_OR) {
      assignMissingIds(orIds, orId, writeIds[i]);
    }
  }
}

void StereoGroups::forwardStereoGroupIds() {
  for (uint32_t i = 0; i < stereoTypes.size(); i++) {
    writeIds[i] = readIds[i];
  }
}

RDMol::~RDMol() {
  // This is needed because default destruction wouldn't delete
  // compatibilityData.
  clear();
}

void RDMol::clear() {
  atomData.resize(0);
  atomBondStarts.resize(0);
  otherAtomIndices.resize(0);
  bondDataIndices.resize(0);
  bondData.resize(0);
  stereoGroups.reset();
  ringInfo.reset();
  substanceGroups.clear();
  numConformers = 0;
  conformerAtomCapacity = 0;
  conformerPositionData.clear();
  conformerIdsAndFlags.clear();
  atomQueries.resize(0);
  bondQueries.resize(0);
  isStereochemDone = false;
  delete compatibilityData.load(std::memory_order_relaxed);
  compatibilityData.store(nullptr, std::memory_order_relaxed);
  properties.resize(0);
  atomBookmarks.clear();
  bondBookmarks.clear();
  monomerInfo.clear();
}

namespace {

unsigned int getEffectiveAtomicNum(const AtomData &atom, bool checkValue,
                                   const atomindex_t atomIndex) {
  auto effectiveAtomicNum =
      static_cast<int>(atom.getAtomicNum()) - atom.getFormalCharge();
  if (checkValue &&
      (effectiveAtomicNum < 0 ||
       effectiveAtomicNum >
           static_cast<int>(PeriodicTable::getTable()->getMaxAtomicNumber()))) {
    throw AtomValenceException("Effective atomic number out of range",
                               atomIndex);
  }
  effectiveAtomicNum = std::clamp(
      effectiveAtomicNum, 0,
      static_cast<int>(PeriodicTable::getTable()->getMaxAtomicNumber()));
  return static_cast<unsigned int>(effectiveAtomicNum);
}

bool canBeHypervalent(const AtomData &atom, unsigned int effectiveAtomicNum) {
  return (effectiveAtomicNum > 16 &&
          (atom.getAtomicNum() == 15 || atom.getAtomicNum() == 16)) ||
         (effectiveAtomicNum > 34 &&
          (atom.getAtomicNum() == 33 || atom.getAtomicNum() == 34));
}

}  // namespace

int RDMol::calculateAtomExplicitValence(atomindex_t atomIndex, bool strict,
                                        bool checkIt) const {
  const AtomData &atom = getAtom(atomIndex);

  // FIX: contributions of bonds to valence are being done at best
  // approximately
  uint32_t accumx2 = 0;
  accumx2 += 2 * atom.getNumExplicitHs();
  if (getNumBonds() > 0) {
    auto [atomBondBegin, atomBondEnd] = getAtomBonds(atomIndex);
    for (; atomBondBegin != atomBondEnd; ++atomBondBegin) {
      const BondData &bond = getBond(*atomBondBegin);
      accumx2 += bond.getTwiceValenceContrib(atomIndex);
    }
  }

  const auto &ovalens =
      PeriodicTable::getTable()->getValenceList(atom.getAtomicNum());
  // if we start with an atom that doesn't have specified valences, we stick
  // with that. otherwise we will use the effective valence
  unsigned int effectiveAtomicNum = atom.getAtomicNum();
  if (ovalens.size() > 1 || ovalens[0] != -1) {
    effectiveAtomicNum = getEffectiveAtomicNum(atom, checkIt, atomIndex);
  }
  unsigned int dv =
      PeriodicTable::getTable()->getDefaultValence(effectiveAtomicNum);
  const auto &valens =
      PeriodicTable::getTable()->getValenceList(effectiveAtomicNum);
  if (accumx2 > 2 * dv && isAromaticAtom(atomIndex)) {
    // this needs some explanation : if the atom is aromatic and
    // accum > dv we assume that no hydrogen can be added
    // to this atom.  We set x = (v + chr) such that x is the
    // closest possible integer to "accum" but less than
    // "accum".
    //
    // "v" here is one of the allowed valences. For example:
    //    sulfur here : O=c1ccs(=O)cc1
    //    nitrogen here : c1cccn1C

    int pval = dv;
    for (auto val : valens) {
      if (val == -1) {
        break;
      }
      if (2 * val > int(accumx2)) {
        break;
      } else {
        pval = val;
      }
    }
    // if we're within 1.5 of the allowed valence, go ahead and take it.
    // this reflects things like the N in c1cccn1C, which starts with
    // accum of 4, but which can be kekulized to C1=CC=CN1C, where
    // the valence is 3 or the bridging N in c1ccn2cncc2c1, which starts
    // with a valence of 4.5, but can be happily kekulized down to a valence
    // of 3
    if (accumx2 - 2 * pval <= 3) {
      accumx2 = 2 * pval;
    }
  }
  // despite promising to not to blame it on him - this a trick Greg
  // came up with: if we have a bond order sum of x.5 (i.e. 1.5, 2.5
  // etc) we would like it to round to the higher integer value --
  // 2.5 to 3 instead of 2 -- so we will add 1 to accumx2.
  // this plays a role in the number of hydrogen that are implicitly
  // added. This will only happen when the accum is a non-integer
  // value and less than the default valence (otherwise the above if
  // statement should have caught it). An example of where this can
  // happen is the following smiles:
  //    C1ccccC1
  // Daylight accepts this smiles and we should be able to Kekulize
  // correctly.
  uint32_t res = (accumx2 + 1) / 2;

  if (strict || checkIt) {
    int maxValence = valens.back();
    int offset = 0;
    // we have to include a special case here for negatively charged P, S, As,
    // and Se, which all support "hypervalent" forms, but which can be
    // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.
    if (canBeHypervalent(atom, effectiveAtomicNum)) {
      maxValence = ovalens.back();
      offset -= atom.getFormalCharge();
    }
    // we have historically accepted two-coordinate [H-] as a valid atom. This
    // is highly questionable, but changing it requires some thought. For now we
    // will just keep accepting it
    if (atom.getAtomicNum() == 1 && atom.getFormalCharge() == -1) {
      maxValence = 2;
    }
    // maxValence == -1 signifies that we'll take anything at the high end
    if (maxValence >= 0 && ovalens.back() >= 0 &&
        (int(res) + offset) > maxValence) {
      // the explicit valence is greater than any
      // allowed valence for the atoms

      if (strict) {
        // raise an error
        std::ostringstream errout;
        errout << "Explicit valence for atom # " << atomIndex << " "
               << PeriodicTable::getTable()->getElementSymbol(
                      atom.getAtomicNum())
               << ", " << res << ", is greater than permitted";
        std::string msg = errout.str();
        BOOST_LOG(rdErrorLog) << msg << std::endl;
        throw AtomValenceException(msg, atomIndex);
      } else {
        return -1;
      }
    }
  }
  return res;
}

// NOTE: this uses the explicitValence, so it will call
// calculateExplicitValence if it is not set on the given atom
int RDMol::calculateAtomImplicitValence(atomindex_t atomIndex, bool strict,
                                        bool checkIt) const {
  const AtomData &atom = getAtom(atomIndex);
  if (atom.noImplicit) {
    return 0;
  }
  auto explicitValence = atom.explicitValence;
  if (explicitValence == AtomData::unsetValenceVal) {
    explicitValence = calculateAtomExplicitValence(atomIndex, strict, checkIt);
  }
  // special cases
  auto atomicNum = atom.getAtomicNum();
  if (atomicNum == 0) {
    return 0;
  }
  for (auto [begin, end] = getAtomBonds(atomIndex); begin != end; ++begin) {
    if (QueryOps::hasComplexBondTypeQuery(ConstRDMolBond(this, *begin))) {
      return 0;
    }
  }

  auto formalCharge = atom.getFormalCharge();
  auto numRadicalElectrons = atom.getNumRadicalElectrons();
  if (explicitValence == 0 && numRadicalElectrons == 0 && atomicNum == 1) {
    if (formalCharge == 1 || formalCharge == -1) {
      return 0;
    } else if (formalCharge == 0) {
      return 1;
    } else {
      if (strict) {
        std::ostringstream errout;
        errout << "Unreasonable formal charge on atom # " << atomIndex << ".";
        std::string msg = errout.str();
        BOOST_LOG(rdErrorLog) << msg << std::endl;
        throw AtomValenceException(msg, atomIndex);
      } else if (checkIt) {
        return -1;
      } else {
        return 0;
      }
    }
  }
  int explicitPlusRadV = atom.explicitValence + atom.getNumRadicalElectrons();

  const auto &ovalens = PeriodicTable::getTable()->getValenceList(atomicNum);
  // if we start with an atom that doesn't have specified valences, we stick
  // with that. otherwise we will use the effective valence for the rest of
  // this.
  unsigned int effectiveAtomicNum = atomicNum;
  if (ovalens.size() > 1 || ovalens[0] != -1) {
    effectiveAtomicNum = getEffectiveAtomicNum(atom, checkIt, atomIndex);
  }
  if (effectiveAtomicNum == 0) {
    return 0;
  }

  // this is basically the difference between the allowed valence of
  // the atom and the explicit valence already specified - tells how
  // many Hs to add
  //

  // The d-block and f-block of the periodic table (i.e. transition metals,
  // lanthanoids and actinoids) have no default valence.
  int dv = PeriodicTable::getTable()->getDefaultValence(effectiveAtomicNum);
  if (dv == -1) {
    return 0;
  }

  // here is how we are going to deal with the possibility of
  // multiple valences
  // - check the explicit valence "ev"
  // - if it is already equal to one of the allowed valences for the
  //    atom return 0
  // - otherwise take return difference between next larger allowed
  //   valence and "ev"
  // if "ev" is greater than all allowed valences for the atom raise an
  // exception
  // finally aromatic cases are dealt with differently - these atoms are allowed
  // only default valences

  // we have to include a special case here for negatively charged P, S, As,
  // and Se, which all support "hypervalent" forms, but which can be
  // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.
  if (canBeHypervalent(atom, effectiveAtomicNum)) {
    effectiveAtomicNum = atomicNum;
    explicitPlusRadV -= formalCharge;
  }
  const auto &valens =
      PeriodicTable::getTable()->getValenceList(effectiveAtomicNum);

  int res = 0;
  // if we have an aromatic case treat it differently
  if (isAromaticAtom(atomIndex)) {
    if (explicitPlusRadV <= dv) {
      res = dv - explicitPlusRadV;
    } else {
      // As we assume when finding the explicitPlusRadValence if we are
      // aromatic we should not be adding any hydrogen and already
      // be at an accepted valence state,

      // FIX: this is just ERROR checking and probably moot - the
      // explicitPlusRadValence function called above should assure us that
      // we satisfy one of the accepted valence states for the
      // atom. The only diff I can think of is in the way we handle
      // formal charge here vs the explicit valence function.
      bool satis = false;
      for (auto vi = valens.begin(); vi != valens.end() && *vi > 0; ++vi) {
        if (explicitPlusRadV == *vi) {
          satis = true;
          break;
        }
      }
      if (!satis && (strict || checkIt)) {
        if (strict) {
          std::ostringstream errout;
          errout << "Explicit valence for aromatic atom # " << atomIndex
                 << " not equal to any accepted valence\n";
          std::string msg = errout.str();
          BOOST_LOG(rdErrorLog) << msg << std::endl;
          throw AtomValenceException(msg, atomIndex);
        } else {
          return -1;
        }
      }
      res = 0;
    }
  } else {
    // non-aromatic case we are allowed to have non default valences
    // and be able to add hydrogens
    res = -1;
    for (auto vi = valens.begin(); vi != valens.end() && *vi >= 0; ++vi) {
      int tot = *vi;
      if (explicitPlusRadV <= tot) {
        res = tot - explicitPlusRadV;
        break;
      }
    }
    if (res < 0) {
      if ((strict || checkIt) && valens.back() != -1 && ovalens.back() > 0) {
        // this means that the explicit valence is greater than any
        // allowed valence for the atoms
        if (strict) {
          // raise an error
          std::ostringstream errout;
          errout << "Explicit valence for atom # " << atomIndex << " "
                 << PeriodicTable::getTable()->getElementSymbol(atomicNum)
                 << " greater than permitted";
          std::string msg = errout.str();
          BOOST_LOG(rdErrorLog) << msg << std::endl;
          throw AtomValenceException(msg, atomIndex);
        } else {
          return -1;
        }
      } else {
        res = 0;
      }
    }
  }
  return res;
}

uint32_t RDMol::calcAtomExplicitValence(atomindex_t atomIndex, bool strict) {
  bool checkIt = false;
  AtomData &atom = getAtom(atomIndex);
  atom.explicitValence =
      calculateAtomExplicitValence(atomIndex, strict, checkIt);
  return atom.explicitValence;
}

uint32_t RDMol::calcAtomImplicitValence(atomindex_t atomIndex, bool strict) {
  AtomData &atom = getAtom(atomIndex);
  if (atom.explicitValence == AtomData::unsetValenceVal) {
    calcAtomExplicitValence(atomIndex, strict);
  }
  bool checkIt = false;
  atom.implicitValence =
      calculateAtomImplicitValence(atomIndex, strict, checkIt);
  return atom.implicitValence;
}

bool RDMol::hasValenceViolation(atomindex_t atomIndex) const {
  // Ignore dummy atoms, query atoms, or atoms attached to query bonds
  auto [begin, end] = getAtomBonds(atomIndex);
  auto is_query = [this](auto b) { return hasBondQuery(b); };
  const AtomData &atom = getAtom(atomIndex);
  auto atomicNum = atom.getAtomicNum();
  if (atomicNum == 0 || hasAtomQuery(atomIndex) ||
      std::any_of(begin, end, is_query)) {
    return false;
  }

  unsigned int effectiveAtomicNum;
  try {
    bool checkIt = true;
    effectiveAtomicNum = getEffectiveAtomicNum(atom, checkIt, atomIndex);
  } catch (const AtomValenceException &) {
    return true;
  }

  // special case for H:
  if (atomicNum == 1) {
    if (atom.getFormalCharge() > 1 || atom.getFormalCharge() < -1) {
      return true;
    }
  } else {
    // Non-H checks for absurd charge values:
    //   1. the formal charge is larger than the atomic number
    //   2. the formal charge moves us to a different row of the periodic table
    if (int(atom.getFormalCharge()) > int(atomicNum) ||
        PeriodicTable::getTable()->getRow(atomicNum) !=
            PeriodicTable::getTable()->getRow(effectiveAtomicNum)) {
      return true;
    }
  }

  bool strict = false;
  bool checkIt = true;
  if (calculateAtomExplicitValence(atomIndex, strict, checkIt) == -1 ||
      calculateAtomImplicitValence(atomIndex, strict, checkIt) == -1) {
    return true;
  }
  return false;
}

// These are copied from Atom.cpp
static const unsigned char octahedral_invert[31] = {
    0,   //  0 -> 0
    2,   //  1 -> 2
    1,   //  2 -> 1
    16,  //  3 -> 16
    14,  //  4 -> 14
    15,  //  5 -> 15
    18,  //  6 -> 18
    17,  //  7 -> 17
    10,  //  8 -> 10
    11,  //  9 -> 11
    8,   // 10 -> 8
    9,   // 11 -> 9
    13,  // 12 -> 13
    12,  // 13 -> 12
    4,   // 14 -> 4
    5,   // 15 -> 5
    3,   // 16 -> 3
    7,   // 17 -> 7
    6,   // 18 -> 6
    24,  // 19 -> 24
    23,  // 20 -> 23
    22,  // 21 -> 22
    21,  // 22 -> 21
    20,  // 23 -> 20
    19,  // 24 -> 19
    30,  // 25 -> 30
    29,  // 26 -> 29
    28,  // 27 -> 28
    27,  // 28 -> 27
    26,  // 29 -> 26
    25   // 30 -> 25
};

static const unsigned char trigonalbipyramidal_invert[21] = {
    0,   //  0 -> 0
    2,   //  1 -> 2
    1,   //  2 -> 1
    4,   //  3 -> 4
    3,   //  4 -> 3
    6,   //  5 -> 6
    5,   //  6 -> 5
    8,   //  7 -> 8
    7,   //  8 -> 7
    11,  //  9 -> 11
    12,  // 10 -> 12
    9,   // 11 -> 9
    10,  // 12 -> 10
    14,  // 13 -> 14
    13,  // 14 -> 13
    20,  // 15 -> 20
    19,  // 16 -> 19
    18,  // 17 -> 28
    17,  // 18 -> 17
    16,  // 19 -> 16
    15   // 20 -> 15
};

bool RDMol::invertAtomChirality(atomindex_t atomIndex) {
  AtomData &atom = getAtom(atomIndex);
  using AtomEnums::ChiralType;
  const ChiralType type = atom.getChiralTag();
  if (type == ChiralType::CHI_TETRAHEDRAL_CW) {
    atom.setChiralTag(ChiralType::CHI_TETRAHEDRAL_CCW);
    return true;
  }
  if (type == ChiralType::CHI_TETRAHEDRAL_CCW) {
    atom.setChiralTag(ChiralType::CHI_TETRAHEDRAL_CW);
    return true;
  }
  uint32_t *permProp = getAtomPropArrayIfPresent<uint32_t>(
      common_properties::_chiralPermutationToken);
  if (permProp == nullptr) {
    return false;
  }

  uint32_t perm = permProp[atomIndex];
  switch (type) {
    case ChiralType::CHI_TETRAHEDRAL:
      if (perm == 1) {
        perm = 2;
      } else if (perm == 2) {
        perm = 1;
      } else {
        perm = 0;
      }
      permProp[atomIndex] = perm;
      return perm != 0;
    case ChiralType::CHI_TRIGONALBIPYRAMIDAL:
      perm = (perm <= 20) ? trigonalbipyramidal_invert[perm] : 0;
      permProp[atomIndex] = perm;
      return perm != 0;
    case ChiralType::CHI_OCTAHEDRAL:
      perm = (perm <= 30) ? octahedral_invert[perm] : 0;
      permProp[atomIndex] = perm;
      return perm != 0;
    default:
      break;
  }
  return false;
}

void RDMol::clearComputedProps(bool includeRings) {
  // the SSSR information:
  if (includeRings) {
    ringInfo.reset();
    // Mark the compatibility layer's ring info as out of sync so it gets
    // updated
    if (hasCompatibilityData()) {
      auto *compat = getCompatibilityDataIfPresent();
      compat->ringInfoSyncStatus.store(CompatSyncStatus::lastUpdatedRDMol,
                                       std::memory_order_release);
    }
  }

  // Clear "computed" properties
  size_t destPropIndex = 0;
  for (size_t sourcePropIndex = 0, numProps = properties.size();
       sourcePropIndex < numProps; ++sourcePropIndex) {
    if (!properties[sourcePropIndex].isComputed()) {
      if (destPropIndex != sourcePropIndex) {
        properties[destPropIndex] = std::move(properties[sourcePropIndex]);
      }
      ++destPropIndex;
    }
  }
  if (destPropIndex != properties.size()) {
    properties.resize(destPropIndex);
  }
}

std::vector<std::string> RDMol::getPropList(bool includePrivate,
                                            bool includeComputed, Scope scope,
                                            uint32_t index) const {
  PRECONDITION(scope == Scope::MOL || index == PropIterator::anyIndexMarker ||
                   (scope == Scope::ATOM && index < getNumAtoms()) ||
                   (scope == Scope::BOND && index < getNumBonds()),
               "RDMol::getPropList index out of range");
  std::vector<std::string> res;
  for (const auto &prop : properties) {
    if (prop.scope() != scope) {
      continue;
    }
    if (scope != Scope::MOL && index != PropIterator::anyIndexMarker &&
        !prop.d_arrayData.isSetMask[index]) {
      continue;
    }
    if (!includeComputed && prop.isComputed()) {
      continue;
    }
    const std::string &name = prop.name().getString();
    if (!includePrivate && name[0] == '_') {
      continue;
    }
    res.push_back(name);
  }
  return res;
}

void RDMol::getComputedPropList(STR_VECT &res, Scope scope,
                                uint32_t index) const {
  res.clear();
  for (const auto &prop : properties) {
    if (prop.scope() != scope) {
      continue;
    }
    if (scope != Scope::MOL && index != PropIterator::anyIndexMarker &&
        !prop.d_arrayData.isSetMask[index]) {
      continue;
    }
    if (!prop.isComputed()) {
      continue;
    }
    res.push_back(prop.name().getString());
  }
}

void RDMol::clearProps() { properties.clear(); }

void RDMol::copyProp(const PropToken &destinationName, const RDMol &sourceMol,
                     const PropToken &sourceName, Scope scope) {
  PRECONDITION(
      scope == Scope::MOL ||
          (scope == Scope::ATOM && getNumAtoms() == sourceMol.getNumAtoms()) ||
          (scope == Scope::BOND && getNumBonds() == sourceMol.getNumBonds()),
      "Atom or bond counts must match in RDMol::copyProp");
  const auto *sourceProp = sourceMol.findProp(sourceName, scope);
  PRECONDITION(sourceProp != nullptr,
               "Source property missing in RDMol::copyProp");
  if (!sourceProp) {
    return;
  }
  auto *existingDestProp = findProp(destinationName, scope);
  if (existingDestProp != nullptr) {
    // No need to copy to self
    if (existingDestProp != sourceProp) {
      // Overwrite existing
      *existingDestProp = *sourceProp;
    }
    return;
  }
  auto &destProp = properties.emplace_back(*sourceProp);
  destProp.d_name = destinationName;
}

void RDMol::copySingleProp(const PropToken &destinationName,
                           uint32_t destinationIndex, const RDMol &sourceMol,
                           const PropToken &sourceName, uint32_t sourceIndex,
                           Scope scope) {
  PRECONDITION(scope != Scope::MOL,
               "RDMol::copyPropSingleIndex doesn't support molecule scope");
  PRECONDITION(
      (scope == Scope::ATOM && sourceIndex < sourceMol.getNumAtoms() &&
       destinationIndex < getNumAtoms()) ||
          (scope == Scope::BOND && sourceIndex < sourceMol.getNumBonds() &&
           destinationIndex < getNumBonds()),
      "atom or bond indices must be in bounds in RDMol::copyPropSingleIndex");
  const auto *sourceProp = sourceMol.findProp(sourceName, scope);
  PRECONDITION(sourceProp != nullptr,
               "Source property missing in RDMol::copyProp");

  auto *destProp = findProp(destinationName, scope);
  if (destProp == nullptr) {
    destProp = &properties.emplace_back();
    destProp->d_name = destinationName;
    destProp->d_isComputed = sourceProp->d_isComputed;
    destProp->d_scope = scope;
    uint32_t destSize = (scope == Scope::ATOM) ? getNumAtoms() : getNumBonds();
    destProp->d_arrayData =
        PropArray(destSize, sourceProp->d_arrayData.family, false);
  } else if (sourceProp->d_arrayData.family != destProp->d_arrayData.family) {
    // Check if types are compatible signed/unsigned integer pairs
    auto sourceFamily = sourceProp->d_arrayData.family;
    auto destFamily = destProp->d_arrayData.family;
    bool integerCompatible = (sourceFamily == PropertyType::INT32 &&
                              destFamily == PropertyType::UINT32) ||
                             (sourceFamily == PropertyType::UINT32 &&
                              destFamily == PropertyType::INT32) ||
                             (sourceFamily == PropertyType::INT64 &&
                              destFamily == PropertyType::UINT64) ||
                             (sourceFamily == PropertyType::UINT64 &&
                              destFamily == PropertyType::INT64);

    if (!integerCompatible) {
      // Convert to RDValue to support type mismatch
      if (destProp->d_arrayData.family != PropertyType::ANY) {
        destProp->d_arrayData.convertToRDValue();
      }
      PRECONDITION(destProp->d_arrayData.family == PropertyType::ANY,
                   "convertToRDValue should make family ANY");
      destProp->d_isComputed = sourceProp->d_isComputed;
      auto *destData = static_cast<RDValue *>(destProp->d_arrayData.data);
      RDValue::cleanup_rdvalue(destData[destinationIndex]);
      static_cast<RDValue *>(destData)[destinationIndex] =
          sourceProp->d_arrayData.getValueAs<RDValue>(sourceIndex);

      // Set the isSetMask when copying with type conversion
      bool &isSet = destProp->d_arrayData.isSetMask[destinationIndex];
      if (!isSet) {
        isSet = true;
        ++destProp->d_arrayData.numSet;
      }
      return;
    }
    // If integer-compatible, fall through to direct copy
  }

  // Copy directly (including integer-compatible types with different family
  // enums)
  destProp->d_isComputed = sourceProp->d_isComputed;
  auto sourceFamily = sourceProp->d_arrayData.family;
  auto destFamily = destProp->d_arrayData.family;
  const auto *sourceData = sourceProp->d_arrayData.data;
  auto *destData = destProp->d_arrayData.data;

  if ((sourceFamily == PropertyType::INT32 ||
       sourceFamily == PropertyType::UINT32) &&
      (destFamily == PropertyType::INT32 ||
       destFamily == PropertyType::UINT32)) {
    // Copy 32-bit integer (signed or unsigned)
    static_cast<int32_t *>(destData)[destinationIndex] =
        static_cast<const int32_t *>(sourceData)[sourceIndex];
  } else if ((sourceFamily == PropertyType::INT64 ||
              sourceFamily == PropertyType::UINT64) &&
             (destFamily == PropertyType::INT64 ||
              destFamily == PropertyType::UINT64)) {
    // Copy 64-bit integer (signed or unsigned)
    static_cast<int64_t *>(destData)[destinationIndex] =
        static_cast<const int64_t *>(sourceData)[sourceIndex];
  } else if (is8BitType(destFamily)) {
    static_cast<char *>(destData)[destinationIndex] =
        static_cast<const char *>(sourceData)[sourceIndex];
  } else if (is32BitType(destFamily)) {
    static_cast<int32_t *>(destData)[destinationIndex] =
        static_cast<const int32_t *>(sourceData)[sourceIndex];
  } else if (is64BitType(destFamily)) {
    static_cast<int64_t *>(destData)[destinationIndex] =
        static_cast<const int64_t *>(sourceData)[sourceIndex];
  } else {
    copy_rdvalue(static_cast<RDValue *>(destData)[destinationIndex],
                 static_cast<const RDValue *>(sourceData)[sourceIndex]);
  }
  bool &isSet = destProp->d_arrayData.isSetMask[destinationIndex];
  if (!isSet) {
    isSet = true;
    ++destProp->d_arrayData.numSet;
  }
}

bool RDMol::hasProp(const PropToken &name) const {
  return findProp(name, Scope::MOL) != nullptr;
}

bool RDMol::hasAtomProp(const PropToken &name,
                        const std::uint32_t index) const {
  const Property *prop = findProp(name, Scope::ATOM);
  if (prop == nullptr) {
    return false;
  }
  return prop->d_arrayData.isSetMask[index];
}

bool RDMol::hasBondProp(const PropToken &name,
                        const std::uint32_t index) const {
  const Property *prop = findProp(name, Scope::BOND);
  if (prop == nullptr) {
    return false;
  }
  return prop->d_arrayData.isSetMask[index];
}

void RDMol::updatePropertyCache(bool strict) {
  for (uint32_t atomIndex = 0, numAtoms = getNumAtoms(); atomIndex < numAtoms;
       ++atomIndex) {
    updateAtomPropertyCache(atomIndex, strict);
  }
  for (uint32_t bondIndex = 0, numBonds = getNumBonds(); bondIndex < numBonds;
       ++bondIndex) {
    updateBondPropertyCache(bondIndex, strict);
  }
}

bool RDMol::needsUpdatePropertyCache() const {
  for (uint32_t atomIndex = 0, numAtoms = getNumAtoms(); atomIndex < numAtoms;
       ++atomIndex) {
    if (getAtom(atomIndex).needsUpdatePropertyCache()) {
      return true;
    }
  }
  // there is no test for bonds yet since they do not obtain a valence property
  return false;
}

void RDMol::clearPropertyCache() {
  for (uint32_t atomIndex = 0, numAtoms = getNumAtoms(); atomIndex < numAtoms;
       ++atomIndex) {
    getAtom(atomIndex).clearPropertyCache();
  }
}

AtomData &RDMol::addAtom() {
  const atomindex_t newAtomIndex = atomData.size();
  auto &newAtom = atomData.emplace_back();
  atomBondStarts.push_back(bondDataIndices.size());

  // Handle properties
  for (Property &property : properties) {
    if (property.scope() != Scope::ATOM) {
      continue;
    }
    property.d_arrayData.appendElement();
  }

  // Resize atomQueries, but only if it's already in use
  if (atomQueries.size() != 0 && atomQueries.size() == newAtomIndex) {
    atomQueries.emplace_back(nullptr);
  }

  if (numConformers > 0) {
    allocateConformers(numConformers, getNumAtoms());
    for (size_t i = 0; i < numConformers; i++) {
      const size_t lastElementAtomIdx =
          i * 3 * conformerAtomCapacity + 3 * newAtomIndex;
      for (int j = 0; j < 3; j++) {
        conformerPositionData[lastElementAtomIdx + j] = 0.0;
      }
    }
  }

  // Handle compat
  auto *compat = getCompatibilityDataIfPresent();
  if (compat != nullptr) {
    // Assume that it's a regular Atom until a query is added
    compat->atoms.emplace_back(new Atom(this, newAtomIndex));

    for (auto &conf : compat->conformers) {
      conf->resize(getNumAtoms());
      conf->setAtomPos(newAtomIndex, RDGeom::Point3D(0.0, 0.0, 0.0));
    }
  }

  return newAtom;
}

BondData &RDMol::addBond(uint32_t beginAtomIdx, uint32_t endAtomIdx,
                         BondEnums::BondType bondType, bool updateAromaticity) {
  URANGE_CHECK(beginAtomIdx, getNumAtoms());
  URANGE_CHECK(endAtomIdx, getNumAtoms());
  PRECONDITION(beginAtomIdx != endAtomIdx,
               "Cannot create a bond between an atom and itself");
  PRECONDITION(getBondIndexBetweenAtoms(beginAtomIdx, endAtomIdx) ==
                   std::numeric_limits<uint32_t>::max(),
               "Bond already exists between atoms");
  const int firstAtomIdx = std::max(beginAtomIdx, endAtomIdx);
  const int secondAtomIdx = std::min(beginAtomIdx, endAtomIdx);

  // Insert the larger atom index first
  const uint32_t newBondIndex = bondData.size();
  auto &newBond = bondData.emplace_back();
  bondDataIndices.insert(
      bondDataIndices.begin() + atomBondStarts[firstAtomIdx + 1],
      bondData.size() - 1);
  bondDataIndices.insert(
      bondDataIndices.begin() + atomBondStarts[secondAtomIdx + 1],
      bondData.size() - 1);
  otherAtomIndices.insert(
      otherAtomIndices.begin() + atomBondStarts[firstAtomIdx + 1],
      secondAtomIdx);
  otherAtomIndices.insert(
      otherAtomIndices.begin() + atomBondStarts[secondAtomIdx + 1],
      firstAtomIdx);

  newBond.bondType = bondType;
  newBond.beginAtomIdx = beginAtomIdx;
  newBond.endAtomIdx = endAtomIdx;

  for (uint32_t i = 0; i < atomBondStarts.size(); i++) {
    auto &atomBondStart = atomBondStarts[i];
    if (i > endAtomIdx) {
      atomBondStart++;
    }
    if (i > beginAtomIdx) {
      atomBondStart++;
    }
  }

  if (bondType == BondEnums::BondType::AROMATIC) {
    newBond.setIsAromatic(true);
    if (updateAromaticity) {
      getAtom(beginAtomIdx).setIsAromatic(true);
      getAtom(endAtomIdx).setIsAromatic(true);
    }
  }

  // Handle properties
  for (Property &property : properties) {
    if (property.scope() != Scope::BOND) {
      continue;
    }
    property.d_arrayData.appendElement();
  }

  // Resize bondQueries, but only if it's already in use
  if (bondQueries.size() != 0 && bondQueries.size() == newBondIndex) {
    bondQueries.emplace_back(nullptr);
  }

  // we're in a batch edit, and at least one of the bond ends is scheduled
  // for deletion, so mark the new bond for deletion too:
  if (dp_delAtoms &&
      ((beginAtomIdx < dp_delAtoms->size() &&
        dp_delAtoms->test(beginAtomIdx)) ||
       (endAtomIdx < dp_delAtoms->size() && dp_delAtoms->test(endAtomIdx)))) {
    if (dp_delBonds->size() < newBondIndex + 1) {
      dp_delBonds->resize(newBondIndex + 1);
    }
    dp_delBonds->set(newBondIndex);
  }

  auto *compat = getCompatibilityDataIfPresent();

  // NOTE: This diverges from reference, but we decided this is the most
  // straightforward way to avoid other problems when unpickling atoms.
  // // reset property cache
  // getAtom(beginAtomIdx).clearPropertyCache();
  // getAtom(endAtomIdx).clearPropertyCache();

  // Handle compat
  if (compat != nullptr) {
    // Assume that it's a regular Bond until a query is added
    compat->bonds.emplace_back(new Bond(this, newBondIndex));
    compat->bondStereoAtoms.push_back(std::make_unique<INT_VECT>());
  }
  return newBond;
}

namespace {

//! Finds all instances of index and removes them from vectors. Decrements all
//! indices greater than index.
void removeIndexAndDecrementBookmarks(
    std::unordered_map<int, std::vector<int>> &bookmarks, int index) {
  for (auto &[idx, marks] : bookmarks) {
    auto removeIt = marks.end();
    for (auto it = marks.begin(); it != marks.end(); ++it) {
      if (*it == int(index)) {
        removeIt = it;
      }
      if (*it > int(index)) {
        --(*it);
      }
    }
    if (removeIt != marks.end()) {
      marks.erase(removeIt);
      // TODO: Should this remove the bookmark if marks is now empty?
    }
  }
}

}  // namespace

void RDMol::removeAtom(atomindex_t atomIndex, bool clearProps) {
  URANGE_CHECK(atomIndex, getNumAtoms());

  if (dp_delAtoms != nullptr) {
    if (dp_delAtoms->size() < getNumAtoms()) {
      dp_delAtoms->resize(getNumAtoms(), false);
    }
    dp_delAtoms->set(atomIndex);
    return;
  }

  // Remove bonds attached to the atom.  removeBond will change this array,
  // so we copy it, first.
  // TODO: To avoid having to copy the atom indices to a new array and delete
  // the bonds based on atom indices, and shift bond properties multiple times,
  // delete all of the bonds in bulk.
  const uint32_t numAtomsOrig = getNumAtoms();
  auto [atomBondsBegin, atomBondsEnd] = getAtomBonds(atomIndex);
  constexpr uint32_t localBufferSize = 8;
  uint32_t localBuffer[localBufferSize];
  uint32_t numBonds = atomBondsEnd - atomBondsBegin;
  std::unique_ptr<uint32_t[]> bondsDeleter;
  uint32_t *bondsToDelete;
  if (numBonds < localBufferSize) {
    bondsToDelete = localBuffer;
  } else {
    bondsToDelete = new uint32_t[numBonds];
    bondsDeleter.reset(bondsToDelete);
  }
  for (uint32_t i = 0; i < numBonds; i++) {
    bondsToDelete[i] = atomBondsBegin[i];
  }
  // Sort in descending order so that we can remove bonds without changing the
  // indexing
  std::sort(bondsToDelete, bondsToDelete + numBonds, std::greater<uint32_t>());

  for (uint32_t i = 0; i < numBonds; ++i) {
    removeBond(bondsToDelete[i]);
  }

  // Update atom indices in otherAtomIndices
  for (unsigned int &a : otherAtomIndices) {
    if (a > atomIndex) {
      --a;
    }
  }
  // Update bond data
  for (BondData &bond : bondData) {
    if (bond.beginAtomIdx > atomIndex) {
      --bond.beginAtomIdx;
    }
    if (bond.endAtomIdx > atomIndex) {
      --bond.endAtomIdx;
    }
  }

  // Update queries
  if (atomQueries.size() != 0) {
    atomQueries.erase(atomQueries.begin() + atomIndex);
  }

  // Handle conformers
  if (atomIndex != numAtomsOrig - 1) {
    const int conformerOffset = conformerAtomCapacity * 3;
    const int copySize =
        sizeof(conformerPositionData[0]) * 3 * (numAtomsOrig - atomIndex - 1);

    for (uint32_t i = 0; i < numConformers; i++) {
      std::memmove(
          &conformerPositionData[i * conformerOffset + atomIndex * 3],
          &conformerPositionData[i * conformerOffset + (atomIndex + 1) * 3],
          copySize);
    }
  }
  CompatibilityData *compat = getCompatibilityDataIfPresent();
  if (compat != nullptr) {
    for (auto &conf : compat->conformers) {
      RDGeom::POINT3D_VECT &positions = conf->getPositions();
      std::copy(positions.begin() + atomIndex + 1, positions.end(),
                positions.begin() + atomIndex);
      positions.pop_back();
    }
  }

  // Update bookmarks
  if (compat != nullptr) {
    const Atom *want = compat->atoms[atomIndex].get();
    std::vector<int> deleteList;
    for (auto &[key, val] : compat->atomBookmarks) {
      auto it = std::find(val.begin(), val.end(), want);
      if (it != val.end()) {
        val.erase(it);
      }
      if (val.empty()) {
        deleteList.push_back(key);
      }
    }
    for (int key : deleteList) {
      compat->atomBookmarks.erase(key);
    }
  } else {
    removeIndexAndDecrementBookmarks(atomBookmarks, atomIndex);
  }

  // Update atom indices in bond stereoAtoms
  for (size_t i = 0, n = getNumBonds(); i < n; ++i) {
    if (hasBondStereoAtoms(i)) {
      const atomindex_t *stereoAtoms = getBondStereoAtoms(i);
      if (stereoAtoms[0] == atomIndex || stereoAtoms[1] == atomIndex) {
        clearBondStereoAtoms(i);
      } else {
        atomindex_t newStereoAtom0 = stereoAtoms[0];
        atomindex_t newStereoAtom1 = stereoAtoms[1];
        if (newStereoAtom0 > atomIndex) {
          --newStereoAtom0;
        }
        if (newStereoAtom1 > atomIndex) {
          --newStereoAtom1;
        }
        setBondStereoAtoms(i, newStereoAtom0, newStereoAtom1);
      }
    }
  }

  // Update substance groups.
  removeSubstanceGroupsReferencingAtom(*this, atomIndex);

  // Update stereo groups.
  if (stereoGroups != nullptr) {
    stereoGroups->removeAtomFromGroups(atomIndex, /*decrement=*/true);
  }
  if (compat != nullptr) {
    std::vector<StereoGroup> *compatStereoGroups =
        compat->stereoGroups.load(std::memory_order_acquire);
    if (compatStereoGroups != nullptr) {
      removeAtomFromGroups(compat->atoms[atomIndex].get(), *compatStereoGroups);
    }
  }

  // Clear computed properties
  if (clearProps) {
    clearComputedProps(true);
  }

  // Shift atom property data back
  for (Property &property : properties) {
    if (property.scope() != Scope::ATOM) {
      continue;
    }
    property.d_arrayData.removeElement(atomIndex);
  }

  // Update special cases of properties referring to atom indices.
  // The _ringStereoAtomsAll, _ringStereoAtomsBegins, and
  // _ringStereoGroup properties are all "computed", so should be removed above.
  // TODO: Do they need to be removed in case !clearProps?

  const Property *molFileBondEndPts =
      findProp(RDKit::common_properties::_MolFileBondEndPtsToken, Scope::BOND);
  if (molFileBondEndPts != nullptr) {
    // Bond end indices may need to be decremented and their
    // indices will need to be handled and if they have an
    // ENDPTS prop that includes idx, it will need updating.
    std::string sprop;
    for (size_t i = 0, n = getNumBonds(); i < n; ++i) {
      if (!getBondPropIfPresent(
              RDKit::common_properties::_MolFileBondEndPtsToken, i, sprop)) {
        continue;
      }
      if (sprop.empty() || sprop.front() != '(' || sprop.back() != ')') {
        continue;
      }

      // This code could be optimized to avoid separate strings
      // or even avoid a vector of ints, if it becomes a bottleneck.
      sprop = sprop.substr(1, sprop.length() - 2);
      boost::char_separator<char> sep(" ");
      boost::tokenizer<boost::char_separator<char>> tokens(sprop, sep);
      unsigned int num_ats = std::stod(*tokens.begin());
      std::vector<unsigned int> oats;
      auto beg = tokens.begin();
      ++beg;
      std::transform(beg, tokens.end(), std::back_inserter(oats),
                     [](const std::string &a) { return std::stod(a); });

      auto idx_pos = std::find(oats.begin(), oats.end(), atomIndex + 1);
      if (idx_pos != oats.end()) {
        oats.erase(idx_pos);
        --num_ats;
      }
      if (!num_ats) {
        clearSingleBondProp(RDKit::common_properties::_MolFileBondEndPtsToken,
                            i);
        clearSingleBondProp(common_properties::_MolFileBondAttachToken, i);
      } else {
        sprop = "(" + std::to_string(num_ats) + " ";
        for (auto &i : oats) {
          if (i > atomIndex + 1) {
            --i;
          }
          sprop += std::to_string(i) + " ";
        }
        sprop[sprop.length() - 1] = ')';
        setSingleBondProp(RDKit::common_properties::_MolFileBondEndPtsToken, i,
                          sprop);
      }
    }
  }

  // There aren't any other special cases yet.
  // Shift back atomData and atomBondStarts
  atomData.erase(atomData.begin() + atomIndex);
  atomBondStarts.erase(atomBondStarts.begin() + atomIndex);

  // Reset our ring info structure, because it is pretty likely
  // to be wrong now. This might have already been done in removeBond or
  // clearComputedProps above, but removing an isolated atom can still
  // invalidate ringInfo.
  ringInfo.reset();

  if (compat != nullptr) {
    compat->ringInfo->reset();
    compat->ringInfoSyncStatus.store(CompatSyncStatus::inSync,
                                     std::memory_order_relaxed);

    for (uint32_t i = atomIndex + 1; i < numAtomsOrig; ++i) {
      compat->atoms[i]->d_index--;
    }
    compat->atoms.erase(compat->atoms.begin() + atomIndex);
  }
}

void RDMol::removeBond(uint32_t bondIndex) {
  // Match ROMol nullop for removing nonexistent bond.
  if (bondIndex >= getNumBonds()) {
    return;
  }

  if (dp_delBonds != nullptr) {
    if (dp_delBonds->size() < getNumBonds()) {
      dp_delBonds->resize(getNumBonds(), false);
    }
    dp_delBonds->set(bondIndex);
    return;
  }

  uint32_t startAtom = getBond(bondIndex).beginAtomIdx;
  uint32_t endAtom = getBond(bondIndex).endAtomIdx;
  const uint32_t numBondsOrig = getNumBonds();
  // Handle props
  for (Property &property : properties) {
    if (property.scope() != Scope::BOND) {
      continue;
    }
    property.d_arrayData.removeElement(bondIndex);
  }

  // Update queries
  if (bondQueries.size() != 0) {
    bondQueries.erase(bondQueries.begin() + bondIndex);
  }

  // remove any bookmarks which point to this bond:
  CompatibilityData *compat = getCompatibilityDataIfPresent();
  if (compat != nullptr) {
    const Bond *want = compat->bonds[bondIndex].get();
    std::vector<int> deleteList;
    for (auto &[key, val] : compat->bondBookmarks) {
      auto it = std::find(val.begin(), val.end(), want);
      if (it != val.end()) {
        val.erase(it);
      }
      if (val.empty()) {
        deleteList.push_back(key);
      }
    }
    for (int key : deleteList) {
      compat->bondBookmarks.erase(key);
    }
  } else {
    removeIndexAndDecrementBookmarks(bondBookmarks, bondIndex);
  }

  // Handle stereo. We loop over all bonds of both atoms, searching to see if
  // the second atom of the pair was listed as the stereo atom of one of these
  // bonds. See RWMol removeBond for details.
  std::array<std::pair<uint32_t, uint32_t>, 2> iterateAndLookupIndices = {
      std::pair<uint32_t, uint32_t>{startAtom, endAtom},
      std::pair<uint32_t, uint32_t>{endAtom, startAtom}};

  for (const auto &[iterateOverBondsIndex, lookForAtomAsStereoIndex] :
       iterateAndLookupIndices) {
    auto [begin, end] = getAtomBonds(iterateOverBondsIndex);
    for (auto bondIterator = begin; bondIterator != end; ++bondIterator) {
      if (*bondIterator == bondIndex) {
        // This is the self bond, it's going away anyway.
        continue;
      }
      auto stereoAtoms = getBondStereoAtoms(*bondIterator);
      if (stereoAtoms[0] == lookForAtomAsStereoIndex ||
          stereoAtoms[1] == lookForAtomAsStereoIndex) {
        clearBondStereoAtoms(*bondIterator);
        // github #6900 if we remove stereo atoms we need to remove
        //  the CIS and or TRANS since this requires stereo atoms
        BondData &oBondData = getBond(*bondIterator);
        const BondEnums::BondStereo stereo = oBondData.getStereo();
        if (stereo == BondEnums::BondStereo::STEREOCIS ||
            stereo == BondEnums::BondStereo::STEREOTRANS) {
          oBondData.setStereo(BondEnums::BondStereo::STEREONONE);
        }
      }
    }
  }

  // reset our ring info structure, because it is pretty likely
  // to be wrong now:
  ringInfo.reset();
  if (compat != nullptr) {
    compat->ringInfo->reset();
    compat->ringInfoSyncStatus.store(CompatSyncStatus::inSync,
                                     std::memory_order_relaxed);
  }
  // Handle substance groups
  removeSubstanceGroupsReferencingBond(*this, bondIndex);

  // Handle stereo groups
  if (stereoGroups != nullptr) {
    stereoGroups->removeBondFromGroups(bondIndex, /*decrement=*/true);
  }
  if (compat != nullptr) {
    std::vector<StereoGroup> *compatStereoGroups =
        compat->stereoGroups.load(std::memory_order_acquire);
    if (compatStereoGroups != nullptr) {
      removeBondFromGroups(compat->bonds[bondIndex].get(), *compatStereoGroups);
    }
  }

  // Finally remove the bond and update index vectors.
  bondData.erase(bondData.begin() + bondIndex);

  // Find the relevant entry in bondDataIndices. Note that these may be out of
  // order.
  uint32_t beginSearchIdx = atomBondStarts[endAtom];
  uint32_t endSearchIdx = atomBondStarts[endAtom + 1];
  while (beginSearchIdx < endSearchIdx) {
    if (bondDataIndices[beginSearchIdx] == bondIndex) {
      break;
    }
    ++beginSearchIdx;
  }
  if (beginSearchIdx == endSearchIdx) {
    throw ValueErrorException("Bond not found in atomBondStarts range");
  }
  uint32_t endToErase = beginSearchIdx;

  // Now do the same for the begin atom.
  beginSearchIdx = atomBondStarts[startAtom];
  endSearchIdx = atomBondStarts[startAtom + 1];
  while (beginSearchIdx < endSearchIdx) {
    if (bondDataIndices[beginSearchIdx] == bondIndex) {
      break;
    }
    ++beginSearchIdx;
  }
  if (beginSearchIdx == endSearchIdx) {
    throw ValueErrorException("Bond not found in atomBondStarts range");
  }
  const uint32_t beginToErase = beginSearchIdx;

  const uint32_t removeIndices[2]{std::min(beginToErase, endToErase),
                                  std::max(beginToErase, endToErase)};
  eraseMultipleIndices(bondDataIndices, removeIndices, 2);
  eraseMultipleIndices(otherAtomIndices, removeIndices, 2);

  for (uint32_t i = 0; i < atomBondStarts.size(); i++) {
    auto &atomBondStart = atomBondStarts[i];
    if (i > endAtom) {
      atomBondStart--;
    }
    if (i > startAtom) {
      atomBondStart--;
    }
  }

  for (uint32_t &dataIndex : bondDataIndices) {
    if (dataIndex > bondIndex) {
      --dataIndex;
    }
  }

  if (compat != nullptr) {
    for (uint32_t i = bondIndex + 1; i < numBondsOrig; ++i) {
      compat->bonds[i]->d_index--;
    }
    compat->bonds.erase(compat->bonds.begin() + bondIndex);
    compat->bondStereoAtoms.erase(compat->bondStereoAtoms.begin() + bondIndex);
  }
}

void RDMol::beginBatchEdit() {
  if (dp_delAtoms != nullptr) {
    throw ValueErrorException("Attempt to re-enter batchEdit mode");
  }
  dp_delAtoms = std::make_unique<boost::dynamic_bitset<>>(getNumAtoms(), false);
  dp_delBonds = std::make_unique<boost::dynamic_bitset<>>(getNumBonds(), false);
}

void RDMol::commitBatchEdit() {
  if (dp_delAtoms == nullptr) {
    return;
  }

  // Find all bonds that would be removed by their atom being removed.
  for (uint32_t i = 0; i < dp_delAtoms->size(); i++) {
    if (!dp_delAtoms->test(i)) {
      continue;
    }
    auto [begin, end] = getAtomBonds(i);
    for (auto bondIterator = begin; bondIterator != end; ++bondIterator) {
      dp_delBonds->set(*bondIterator);
    }
  }

  auto movedDelAtoms = std::move(dp_delAtoms);
  auto movedDelBonds = std::move(dp_delBonds);
  dp_delAtoms.reset();
  dp_delBonds.reset();

  // Now delete bonds in reverse order
  for (int i = movedDelBonds->size() - 1; i >= 0; i--) {
    if (movedDelBonds->test(i)) {
      removeBond(i);
    }
  }

  // Delete atoms in reverse order
  for (int i = movedDelAtoms->size() - 1; i >= 0; i--) {
    if (movedDelAtoms->test(i)) {
      removeAtom(i);
    }
  }

  dp_delAtoms.reset();
  dp_delBonds.reset();
}

void RDMol::rollbackBatchEdit() {
  dp_delAtoms.reset();
  dp_delBonds.reset();
}

// -----------------------------
// Atom bookmarks
// -----------------------------

void RDMol::setAtomBookmark(int atomIdx, int mark) {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    compatData->atomBookmarks[mark].push_back(compatData->atoms[atomIdx].get());
    return;
  }

  // Non legacy path.
  atomBookmarks[mark].push_back(atomIdx);
}
void RDMol::replaceAtomBookmark(int atomIdx, int mark) {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    auto &bookmark = compatData->atomBookmarks[mark];
    bookmark.clear();
    bookmark.push_back(compatData->atoms[atomIdx].get());
    return;
  }
  // Nonlegacy path
  std::vector<int> &marks = atomBookmarks[mark];
  marks.resize(1);
  marks[0] = atomIdx;
}
std::vector<int> RDMol::getAllAtomsWithBookmarks(int mark) {
  if (hasCompatibilityData()) {
    std::list<Atom *> &bookmarks = getAtomBookmarksCompat(mark);
    std::vector<int> res;
    res.reserve(bookmarks.size());
    for (Atom *atom : bookmarks) {
      res.push_back(atom->getIdx());
    }
    return res;
  }

  // Non-legacy path
  auto it = atomBookmarks.find(mark);
  PRECONDITION(it != atomBookmarks.end(), "Bookmark not found");
  return it->second;
}

int RDMol::getAtomWithBookmark(int mark) {
  const CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    auto found = compatData->atomBookmarks.find(mark);
    PRECONDITION(found != compatData->atomBookmarks.end(),
                 "Bookmark not found");
    return found->second.front()->getIdx();
  }

  // Non-legacy path
  auto it = atomBookmarks.find(mark);
  PRECONDITION(it != atomBookmarks.end(), "Bookmark not found");
  return *it->second.begin();
}
void RDMol::clearAtomBookmark(int mark) {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    compatData->atomBookmarks.erase(mark);
    return;
  }

  // non legacy path
  atomBookmarks.erase(mark);
}

void RDMol::clearAtomBookmark(int atomIdx, int mark) {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    auto foundMark = compatData->atomBookmarks.find(mark);
    if (foundMark == compatData->atomBookmarks.end()) {
      return;
    }
    auto &marks = foundMark->second;
    auto removeIt = std::find_if(
        marks.begin(), marks.end(),
        [atomIdx](const Atom *atom) { return int(atom->getIdx()) == atomIdx; });
    if (removeIt != marks.end()) {
      marks.erase(removeIt);
    }
    return;
  }

  auto it = atomBookmarks.find(mark);
  if (it == atomBookmarks.end()) {
    return;
  }
  std::vector<int> &marks = it->second;
  auto removeIt = std::find(marks.begin(), marks.end(), atomIdx);
  if (removeIt != marks.end()) {
    marks.erase(removeIt);
  }
}

void RDMol::clearAllAtomBookmarks() {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    compatData->atomBookmarks.clear();
  }
  atomBookmarks.clear();
}

bool RDMol::hasAtomBookmark(int mark) const {
  const CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    return compatData->atomBookmarks.count(mark) > 0;
  }

  if (auto it = atomBookmarks.find(mark); it != atomBookmarks.end()) {
    return !it->second.empty();
  }
  return false;
}

// -----------------------------
// Bond bookmarks
// -----------------------------

void RDMol::setBondBookmark(int bondIdx, int mark) {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    compatData->bondBookmarks[mark].push_back(compatData->bonds[bondIdx].get());
    return;
  }

  // Non legacy path.
  bondBookmarks[mark].push_back(bondIdx);
}

std::vector<int> RDMol::getAllBondsWithBookmarks(int mark) {
  if (hasCompatibilityData()) {
    std::list<Bond *> bookmarks = getBondBookmarksCompat(mark);
    std::vector<int> res;
    res.reserve(bookmarks.size());
    for (Bond *bond : bookmarks) {
      res.push_back(bond->getIdx());
    }
    return res;
  }

  // Non-legacy path
  auto it = bondBookmarks.find(mark);
  PRECONDITION(it != bondBookmarks.end(), "Bookmark not found");
  return it->second;
}

int RDMol::getBondWithBookmark(int mark) {
  const CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    auto found = compatData->bondBookmarks.find(mark);
    PRECONDITION(found != compatData->bondBookmarks.end(),
                 "Bookmark not found");
    return found->second.front()->getIdx();
  }

  // Non-legacy path
  auto it = bondBookmarks.find(mark);
  PRECONDITION(it != bondBookmarks.end(), "Bookmark not found");
  return *it->second.begin();
}

void RDMol::clearBondBookmark(int mark) {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    compatData->bondBookmarks.erase(mark);
    return;
  }

  // non legacy path
  bondBookmarks.erase(mark);
}
void RDMol::clearBondBookmark(int bondIdx, int mark) {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    auto foundMark = compatData->bondBookmarks.find(mark);
    if (foundMark == compatData->bondBookmarks.end()) {
      return;
    }
    auto &marks = foundMark->second;
    auto removeIt = std::find_if(
        marks.begin(), marks.end(),
        [bondIdx](const Bond *bond) { return int(bond->getIdx()) == bondIdx; });
    if (removeIt != marks.end()) {
      marks.erase(removeIt);
    }
    return;
  }

  auto it = bondBookmarks.find(mark);
  if (it == bondBookmarks.end()) {
    return;
  }
  std::vector<int> &marks = it->second;
  auto removeIt = std::find(marks.begin(), marks.end(), bondIdx);
  if (removeIt != marks.end()) {
    marks.erase(removeIt);
  }
}

void RDMol::clearAllBondBookmarks() {
  CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    compatData->bondBookmarks.clear();
  }
  bondBookmarks.clear();
}

bool RDMol::hasBondBookmark(int mark) const {
  const CompatibilityData *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    return compatData->bondBookmarks.count(mark) > 0;
  }

  if (auto it = bondBookmarks.find(mark); it != bondBookmarks.end()) {
    return !it->second.empty();
  }
  return false;
}

void RDMol::setStereoGroups(std::unique_ptr<StereoGroups> &&groups) {
  clearStereoGroupsCompat();

  stereoGroups = std::move(groups);

  // Repopulate the compatibility layer immediately to keep existing Python
  // references valid The vector object is kept alive by
  // clearStereoGroupsCompat, we just need to refill it
  if (hasCompatibilityData()) {
    getStereoGroupsCompat();
  }
}

void RDMol::markStereoGroupsAsCompatModified() const {
  if (const CompatibilityData *compatData = getCompatibilityDataIfPresent();
      compatData != nullptr) {
    compatData->stereoGroupsSyncStatus = CompatSyncStatus::lastUpdatedCompat;
  }
}

void RDMol::setAtomQuery(atomindex_t atomIndex,
                         std::unique_ptr<QUERYATOM_QUERY> query) {
  if (atomQueries.size() != getNumAtoms()) {
    atomQueries.resize(getNumAtoms());
  }
  atomQueries[atomIndex] = std::move(query);

  auto *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    Atom *atom = compatData->atoms[atomIndex].get();
    QueryAtom *queryAtom = dynamic_cast<QueryAtom *>(atom);
    auto *newQuery = makeNewCompatQuery(atomQueries[atomIndex].get());
    if (queryAtom == nullptr) {
      // Convert to a QueryAtom now that there's a query
      compatData->atoms[atomIndex].reset(
          new QueryAtom(this, atomIndex, newQuery));
    } else {
      queryAtom->setQueryPrivate(newQuery);
    }
  }
}
void RDMol::setBondQuery(uint32_t bondIndex,
                         std::unique_ptr<QUERYBOND_QUERY> query) {
  if (bondQueries.size() != getNumBonds()) {
    bondQueries.resize(getNumBonds());
  }
  bondQueries[bondIndex] = std::move(query);

  auto *compatData = getCompatibilityDataIfPresent();
  if (compatData != nullptr) {
    Bond *bond = compatData->bonds[bondIndex].get();
    QueryBond *queryBond = dynamic_cast<QueryBond *>(bond);
    auto *newQuery = makeNewCompatQuery(bondQueries[bondIndex].get());
    if (queryBond == nullptr) {
      // Convert to a QueryBond now that there's a query
      compatData->bonds[bondIndex].reset(
          new QueryBond(this, bondIndex, newQuery));
    } else {
      queryBond->setQueryPrivate(newQuery);
    }
  }
}

ROMol &RDMol::asROMol() { return *ensureCompatInit()->compatMol; }
const ROMol &RDMol::asROMol() const { return *ensureCompatInit()->compatMol; }
RWMol &RDMol::asRWMol() {
  auto *mol = dynamic_cast<RWMol *>(&asROMol());
  CHECK_INVARIANT(mol, "Failed to cast ROMol to RWMol");
  return *mol;
}

const RWMol &RDMol::asRWMol() const {
  const auto *mol = dynamic_cast<const RWMol *>(&asROMol());
  CHECK_INVARIANT(mol, "Failed to cast ROMol to RWMol");
  return *mol;
}

RDMol::CompatibilityData *RDMol::ensureCompatInit() const {
  // First, try without locking
  // We can do a relaxed load of the pointer, because accessing the other data
  // requires dereferencing the pointer, so would all be dependent reads.
  CompatibilityData *data = compatibilityData.load(std::memory_order_relaxed);
  if (data != nullptr) {
    return data;
  }

  std::lock_guard<std::mutex> lock_scope(compatibilityMutex);
  // Check again after locking
  data = compatibilityData.load(std::memory_order_relaxed);
  if (data != nullptr) {
    return data;
  }

  // NOTE: This casts away const because Atom and Bond structures may later
  // need write access.  This is in a private function for initializing
  // a mutable data member, so it should be okay.
  data = new CompatibilityData(const_cast<RDMol &>(*this));
  compatibilityData.store(data, std::memory_order_release);
  return data;
}

uint32_t RDMol::getNumConformers() const { return numConformers; }

namespace {

[[noreturn]] void throwIDNotFound(int id) {
  std::string exceptionMessage =
      "Can't find conformation with ID: " + std::to_string(id);
  throw ConformerException(std::move(exceptionMessage));
}

}  // namespace

uint32_t RDMol::findConformerIndex(int32_t id) const {
  // make sure we have more than one conformation
  if (numConformers == 0) {
    throw ConformerException("No conformations available on the molecule");
  }
  if (id < 0) {
    // Negative id means index 0
    return 0;
  }

  if (conformerIdsAndFlags.size() == 0) {
    // Fast path for ids being equal to indices
    if (size_t(id) >= numConformers) {
      throwIDNotFound(id);
    }
    return id;
  }

  CHECK_INVARIANT(conformerIdsAndFlags.size() == numConformers,
                  "conformerIdsAndFlags incorrect size in findConformerIndex");
  for (uint32_t i = 0; i < numConformers; ++i) {
    if (conformerIdsAndFlags[i].id == uint32_t(id)) {
      return i;
    }
  }
  throwIDNotFound(id);
}

const double *RDMol::getConformerPositions(int32_t id) const {
  auto *compat = getCompatibilityDataIfPresent();
  if (compat != nullptr) {
    initCompatForReadHelper<false>(
        compat->conformerSyncStatus, compatibilityMutex,
        [&]() { copyConformersFromCompatibilityData(compat); });
  }
  uint32_t index = findConformerIndex(id);
  CHECK_INVARIANT(conformerAtomCapacity >= getNumAtoms(),
                  "Conformer positions not allocated in getConformerPositions");
  return conformerPositionData.data() +
         size_t(3) * conformerAtomCapacity * index;
}

double *RDMol::getConformerPositions(int32_t id) {
  auto *compat = getCompatibilityDataIfPresent();
  if (compat != nullptr) {
    initCompatForWriteHelper<false>(
        compat->conformerSyncStatus, compatibilityMutex,
        [&]() { copyConformersFromCompatibilityData(compat); });
  }
  uint32_t index = findConformerIndex(id);
  CHECK_INVARIANT(conformerAtomCapacity >= getNumAtoms(),
                  "Conformer positions not allocated in getConformerPositions");
  return conformerPositionData.data() +
         size_t(3) * conformerAtomCapacity * index;
}

bool RDMol::is3DConformer(uint32_t id) const {
  auto *compat = getCompatibilityDataIfPresent();
  if (compat != nullptr) {
    initCompatForReadHelper<false>(
        compat->conformerSyncStatus, compatibilityMutex,
        [&]() { copyConformersFromCompatibilityData(compat); });
  }
  uint32_t index = findConformerIndex(id);
  if (index == numConformers) {
    throwIDNotFound(id);
  }
  return (conformerIdsAndFlags.size() == 0) || conformerIdsAndFlags[index].is3D;
}

void RDMol::clearConformers() {
  // Free the memory by swapping with empty vectors
  std::vector<double> tempDouble;
  conformerPositionData.swap(tempDouble);
  std::vector<ConformerIdAndFlags> tempIdsAndFlags;
  conformerIdsAndFlags.swap(tempIdsAndFlags);
  numConformers = 0;
  conformerAtomCapacity = 0;

  if (auto *compatData = getCompatibilityDataIfPresent();
      compatData != nullptr) {
    compatData->conformers.clear();
    // It might be excessive to use memory_order_release, but this will ensure
    // that the clearing of the conformers is completed before marking them
    // as in sync.
    compatData->conformerSyncStatus.store(CompatSyncStatus::inSync,
                                          std::memory_order_release);
  }
}

static void expandConformerIdsAndFlags(
    std::vector<ConformerIdAndFlags> &conformerIdsAndFlags,
    size_t numConformers) {
  PRECONDITION(conformerIdsAndFlags.size() == 0,
               "Expanding conformerIdsAndFlags that's already expanded");
  conformerIdsAndFlags.reserve(numConformers);
  for (size_t i = 0; i < numConformers; ++i) {
    auto &confData = conformerIdsAndFlags.emplace_back();
    confData.id = i;
    confData.is3D = true;
  }
}

void RDMol::copyConformersFromCompatibilityData(
    const CompatibilityData *compatData, int confId) const {
  PRECONDITION(compatData != nullptr, "No compatibility data to copy from");
  PRECONDITION(
      (confId >= 0 && numConformers <= 1) ||
          compatData->conformers.size() == numConformers,
      "Add space for conformers before copying from compatibility data");
  size_t confIndex = 0;
  const size_t numAtoms = getNumAtoms();
  for (const auto &conf : compatData->conformers) {
    if (confId < 0 || int(conf->getId()) == confId) {
      PRECONDITION(
          confIndex < numConformers,
          "Add space for conformers before copying from compatibility data");
      size_t atomPositionOffset = confIndex * conformerAtomCapacity * 3;
      const RDGeom::POINT3D_VECT &positions = conf->getPositions();
      PRECONDITION(positions.size() == numAtoms,
                   "Conformers must have one position per atom");
      for (size_t i = 0; i < numAtoms; ++i) {
        conformerPositionData[atomPositionOffset] = positions[i].x;
        conformerPositionData[atomPositionOffset + 1] = positions[i].y;
        conformerPositionData[atomPositionOffset + 2] = positions[i].z;
        atomPositionOffset += 3;
      }

      if (conformerIdsAndFlags.size() > 0) {
        conformerIdsAndFlags[confIndex].is3D = conf->is3D();
        conformerIdsAndFlags[confIndex].id = conf->getId();
      } else if (!conf->is3D() || conf->getId() != confIndex) {
        expandConformerIdsAndFlags(conformerIdsAndFlags, numConformers);
        conformerIdsAndFlags[confIndex].is3D = conf->is3D();
        conformerIdsAndFlags[confIndex].id = conf->getId();
      }
      confIndex++;
    }
  }
}

void RDMol::copyConformersToCompatibilityData(
    CompatibilityData *compatData) const {
  PRECONDITION(compatData != nullptr, "No compatibility data to copy to");
  PRECONDITION(compatData->conformers.size() == numConformers,
               "Add conformers before copying to compatibility data");
  size_t confIndex = 0;
  const size_t numAtoms = getNumAtoms();
  for (auto &conf : compatData->conformers) {
    size_t atomPositionOffset = confIndex * conformerAtomCapacity * 3;
    RDGeom::POINT3D_VECT &positions = conf->getPositions();
    PRECONDITION(positions.size() == numAtoms,
                 "Conformers must have one position per atom");
    for (size_t i = 0; i < numAtoms; ++i) {
      positions[i].x = conformerPositionData[atomPositionOffset];
      positions[i].y = conformerPositionData[atomPositionOffset + 1];
      positions[i].z = conformerPositionData[atomPositionOffset + 2];
      atomPositionOffset += 3;
    }
    conf->set3D(conformerIdsAndFlags.size() == 0 ||
                conformerIdsAndFlags[confIndex].is3D);
    confIndex++;
  }
}

void RDMol::markConformersAsCompatModified() const {
  if (const CompatibilityData *compatData = getCompatibilityDataIfPresent();
      compatData != nullptr) {
    compatData->conformerSyncStatus = CompatSyncStatus::lastUpdatedCompat;
  }
}

void RDMol::removeConformer(uint32_t id) {
  if (int32_t(id) < 0) {
    return;
  }
  uint32_t index;
  try {
    index = findConformerIndex(id);
  } catch (const ConformerException &) {
    return;
  }

  const uint32_t numAtoms = getNumAtoms();
  if (index != numConformers - 1) {
    // Overwrite the position at index with the position at numConformers - 1
    if (numAtoms > 0) {
      // Stride is 3*conformerAtomCapacity, but amount to copy is 3*numAtoms
      double *dest = conformerPositionData.data() +
                     size_t(3) * conformerAtomCapacity * index;
      const double *source =
          conformerPositionData.data() +
          size_t(3) * conformerAtomCapacity * (numConformers - 1);
      std::copy(source, source + size_t(3) * numAtoms, dest);
    }

    if (conformerIdsAndFlags.size() == 0) {
      // We can no longer have 0 to n based indexing, populate
      // conformerIdsAndFlags.
      expandConformerIdsAndFlags(conformerIdsAndFlags, numConformers - 1);
      conformerIdsAndFlags[index].id = numConformers - 1;
    } else {
      conformerIdsAndFlags[index] = conformerIdsAndFlags[numConformers - 1];
      conformerIdsAndFlags.pop_back();
    }
  } else if (conformerIdsAndFlags.size() != 0) {
    // Removing last conformer, so remove the last entry from
    // conformerIdsAndFlags
    conformerIdsAndFlags.pop_back();
  }

  if (auto *compatData = getCompatibilityDataIfPresent();
      compatData != nullptr) {
    for (auto iter = compatData->conformers.begin();
         iter != compatData->conformers.end(); ++iter) {
      if ((*iter)->getId() == id) {
        compatData->conformers.erase(iter);
        break;
      }
    }
  }

  numConformers--;
}

uint32_t RDMol::addConformer(const double *positions, int32_t id, bool is3D) {
  allocateConformers(numConformers + 1);
  const size_t numConformersPrev = numConformers;
  ++numConformers;

  if (getNumAtoms() > 0) {
    if (positions != nullptr) {
      std::copy(positions, positions + getNumAtoms() * 3,
                conformerPositionData.data() +
                    3 * conformerAtomCapacity * numConformersPrev);
    }
  }

  // Add to compat data. Only need if this is not being invoked as part of the
  // ROMol addConformer, so we check to make sure ROMol hasn't already appended
  // for itself.
  auto *compatData = getCompatibilityDataIfPresent();
  const bool needsCompatUpdate =
      compatData != nullptr && compatData->conformers.size() < numConformers;
  if (needsCompatUpdate) {
    auto *conf = new Conformer(getNumAtoms());
    if (positions != nullptr) {
      for (size_t i = 0; i < getNumAtoms(); ++i) {
        conf->setAtomPos(i,
                         RDGeom::Point3D(positions[i * 3], positions[i * 3 + 1],
                                         positions[i * 3 + 2]));
      }
    }
    conf->set3D(is3D);
    compatData->conformers.push_back(CONFORMER_SPTR(conf));
  }

  const bool hadLinearIds = conformerIdsAndFlags.size() == 0;
  // Case where we're allowed to set ID, and default flags.
  const bool allowsLinearIds =
      is3D && (id < 0 || size_t(id) == numConformersPrev);

  if (allowsLinearIds && hadLinearIds) {
    if (needsCompatUpdate) {
      compatData->conformers.back()->setId(numConformersPrev);
    }
    return numConformersPrev;
  }

  if (hadLinearIds) {
    // We now need to populate all IDs that were previously implicit.
    // We know they were 3D and the implict IDs were in order.
    expandConformerIdsAndFlags(conformerIdsAndFlags, numConformers);
    auto &newConfData = conformerIdsAndFlags.back();
    newConfData.is3D = is3D;
    newConfData.id = id < 0 ? numConformersPrev : id;
    if (needsCompatUpdate) {
      compatData->conformers.back()->setId(newConfData.id);
    }
    return newConfData.id;
  }

  // Previous IDs are populated.
  auto &newConfData = conformerIdsAndFlags.emplace_back();

  if (id < 0) {
    // We need to find an ID to give the new conformer.
    int64_t max = -1;
    for (const auto &confData : conformerIdsAndFlags) {
      max = std::max(max, int64_t(confData.id));
    }
    id = max + 1;
  }
  newConfData.id = id;
  newConfData.is3D = is3D;

  // Add to compat data
  if (needsCompatUpdate) {
    compatData->conformers.back()->setId(id);
  }
  return newConfData.id;
}

std::pair<uint32_t, double *> RDMol::addConformer(int32_t id, bool is3D) {
  uint32_t gotId = addConformer(nullptr, id, is3D);
  return {gotId, getConformerPositions(gotId)};
}

void RDMol::allocateConformers(size_t newNumConformers, size_t newNumAtoms) {
  size_t newConformerAtomCapacity = conformerAtomCapacity;
  if (conformerAtomCapacity < newNumAtoms) {
    newConformerAtomCapacity = std::max(2 * conformerAtomCapacity, newNumAtoms);
  }
  size_t numConformersCapacity =
      (conformerAtomCapacity != 0)
          ? (conformerPositionData.size() / (size_t(3) * conformerAtomCapacity))
          : 0;
  if (numConformersCapacity < newNumConformers) {
    numConformersCapacity =
        std::max(2 * numConformersCapacity, newNumConformers);
  }

  conformerPositionData.resize(size_t(3) * newConformerAtomCapacity *
                               numConformersCapacity);

  if (newConformerAtomCapacity != conformerAtomCapacity) {
    // This only works if we can't decrease the atom capacity
    assert(newConformerAtomCapacity > conformerAtomCapacity);

    // Copy the previous data to the new locations, with a backward loop and
    // memmove to avoid issues with overlap.
    size_t numConformersToCopy =
        (getNumAtoms() != 0) ? std::min(numConformers, newNumConformers) : 0;
    const size_t oldConformerStride = 3 * conformerAtomCapacity;
    const size_t newConformerStride = 3 * newConformerAtomCapacity;
    for (size_t i = numConformersToCopy; i > 0;) {
      --i;
      auto *dest = conformerPositionData.data() + i * newConformerStride;
      const auto *source =
          conformerPositionData.data() + i * oldConformerStride;
      // Move the existing data
      std::memmove(dest, source,
                   sizeof(conformerPositionData[0]) * 3 * getNumAtoms());
    }

    conformerAtomCapacity = newConformerAtomCapacity;
  }
}

const RingInfoCache &RDMol::getRingInfo() const {
  auto *compat = getCompatibilityDataIfPresent();
  if (compat != nullptr) {
    initCompatForReadHelper<false>(
        compat->ringInfoSyncStatus, compatibilityMutex, [&]() {
          copyRingInfoFromCompatibilityData(*compat->ringInfo, ringInfo,
                                            getNumAtoms(), getNumBonds());
        });
  }
  return ringInfo;
}
RingInfoCache &RDMol::getRingInfo() {
  auto *compat = getCompatibilityDataIfPresent();
  if (compat != nullptr) {
    initCompatForWriteHelper<false>(
        compat->ringInfoSyncStatus, compatibilityMutex, [&]() {
          copyRingInfoFromCompatibilityData(*compat->ringInfo, ringInfo,
                                            getNumAtoms(), getNumBonds());
        });
  }
  return ringInfo;
}

}  // namespace RDKit
