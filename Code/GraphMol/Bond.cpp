//
//  Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Bond.h"
#include "Atom.h"
#include "ROMol.h"
#include <RDGeneral/Invariant.h>
#include "Atropisomers.h"

namespace RDKit {

void Bond::initFromOther(const Bond &other, const bool preserveProps) {
  auto& newBondInfo = dp_dataMol->getBond(d_index);
  newBondInfo = other.dp_dataMol->getBond(other.d_index);

  if (!preserveProps) {
    dp_dataMol->clearSingleBondAllProps(d_index);
    for (auto it = other.dp_dataMol->beginProps(true, RDMol::Scope::BOND,
                                                other.d_index),
              end = other.dp_dataMol->endProps();
         it != end; ++it) {
      dp_dataMol->copySingleProp(*it, d_index, *other.dp_dataMol, *it,
                                 other.d_index, RDMol::Scope::BOND);
    }
  }

  // These may come from compat data
  const atomindex_t *stereoAtoms = other.dp_dataMol->getBondStereoAtoms(other.d_index);
  if (stereoAtoms != nullptr && stereoAtoms[0] != atomindex_t(-1) && stereoAtoms[1] != atomindex_t(-1)) {
    // Don't check the atom indices if this is a Bond outside a molecule
    dp_dataMol->setBondStereoAtoms(d_index, stereoAtoms[0], stereoAtoms[1],
                                   dp_dataMol->getNumAtoms() != 0);
  } else {
    dp_dataMol->clearBondStereoAtoms(d_index);
  }
}

Bond::Bond() {
  dp_dataMol = new RDMol();
  dp_dataMol->addUnconnectedBond();
  dp_owningMol = nullptr;
};

Bond::Bond(BondType bT): Bond() {
  setBondType(bT);
}

Bond::Bond(const Bond& other): Bond() {
  initFromOther(other);
}

Bond& Bond::operator=(const Bond &other) {
  if (this == &other) {
    return *this;
  }
  *this = Bond(other);
  return *this;
}

Bond::Bond(Bond &&o) noexcept {
  dp_dataMol = o.dp_dataMol;
  dp_owningMol = o.dp_owningMol;
  d_index = o.d_index;
  o.dp_dataMol = nullptr;
  o.dp_owningMol = nullptr;
}

Bond &Bond::operator=(Bond &&o) noexcept {
  if (this == &o) {
    return *this;
  }
  if (dp_dataMol != dp_owningMol) {
    delete dp_dataMol;
  }
  dp_dataMol = o.dp_dataMol;
  dp_owningMol = o.dp_owningMol;
  d_index = o.d_index;
  o.dp_dataMol = nullptr;
  o.dp_owningMol = nullptr;
  return *this;
}

Bond *Bond::copy() const {
  auto *bond = new Bond(*this);
  return bond;
}

ROMol &Bond::getOwningMol() const {
  PRECONDITION(dp_owningMol != nullptr, "no owner");
  return dp_owningMol->asROMol();
}

void Bond::setOwningMol(ROMol *other) {
  setOwningMol(&other->asRDMol());
}

void Bond::setOwningMol(RDMol* other) {
  PRECONDITION(dp_owningMol == nullptr || dp_owningMol == other,
               "setOwningMol called twice");
  dp_owningMol = other;
}

Bond::~Bond() {
  if (dp_dataMol != dp_owningMol) {
    delete dp_dataMol;
  }
}

Bond::BondType Bond::getBondType() const {
  return static_cast<BondType>(dp_dataMol->getBond(d_index).getBondType());
}
void Bond::setBondType(BondType bT) {
  dp_dataMol->getBond(d_index).setBondType(
      static_cast<RDKit::BondEnums::BondType>(bT));
}
double Bond::getBondTypeAsDouble() const {
  double res;
  switch (getBondType()) {
    case UNSPECIFIED:
    case IONIC:
    case ZERO:
      res = 0;
      break;
    case SINGLE:
      res = 1;
      break;
    case DOUBLE:
      res = 2;
      break;
    case TRIPLE:
      res = 3;
      break;
    case QUADRUPLE:
      res = 4;
      break;
    case QUINTUPLE:
      res = 5;
      break;
    case HEXTUPLE:
      res = 6;
      break;
    case ONEANDAHALF:
      res = 1.5;
      break;
    case TWOANDAHALF:
      res = 2.5;
      break;
    case THREEANDAHALF:
      res = 3.5;
      break;
    case FOURANDAHALF:
      res = 4.5;
      break;
    case FIVEANDAHALF:
      res = 5.5;
      break;
    case AROMATIC:
      res = 1.5;
      break;
    case DATIVEONE:
      res = 1.0;
      break;  // FIX: this should probably be different
    case DATIVE:
      res = 1.0;
      break;  // FIX: again probably wrong
    case HYDROGEN:
      res = 0.0;
      break;
    default:
      UNDER_CONSTRUCTION("Bad bond type");
  }
  return res;
}
double Bond::getValenceContrib(const Atom* atom) const {
  return 0.5 * dp_dataMol->getBond(d_index).getTwiceValenceContrib(atom->getIdx());
}
void Bond::setIsAromatic(bool what) {
  dp_dataMol->getBond(d_index).setIsAromatic(what);
}
bool Bond::getIsAromatic() const {
  return dp_dataMol->getBond(d_index).getIsAromatic();
}
void Bond::setIsConjugated(bool what) {
  dp_dataMol->getBond(d_index).setIsConjugated(what);
}
bool Bond::getIsConjugated() const {
  return dp_dataMol->getBond(d_index).getIsConjugated();
}
unsigned int Bond::getBeginAtomIdx() const {
  return dp_dataMol->getBond(d_index).getBeginAtomIdx();
}
unsigned int Bond::getEndAtomIdx() const {
  return dp_dataMol->getBond(d_index).getEndAtomIdx();
}
unsigned int Bond::getOtherAtomIdx(unsigned int thisIdx) const {
  return dp_dataMol->getBond(d_index).getOtherAtomIdx(thisIdx);
}
void Bond::setBeginAtomIdx(unsigned int what) {
  if (dp_owningMol) {
    URANGE_CHECK(what, getOwningMol().getNumAtoms());
  }
  dp_dataMol->getBond(d_index).setBeginAtomIdx(what);
}
void Bond::setEndAtomIdx(unsigned int what) {
  if (dp_owningMol) {
    URANGE_CHECK(what, getOwningMol().getNumAtoms());
  }
  dp_dataMol->getBond(d_index).setEndAtomIdx(what);
}
void Bond::setBeginAtom(Atom* at) {
  PRECONDITION(dp_owningMol != nullptr, "no owning molecule for bond");
  // TODO: Should this assert that at is from the same molecule?
  dp_dataMol->getBond(d_index).setBeginAtomIdx(at->getIdx());
}
void Bond::setEndAtom(Atom* at) {
  PRECONDITION(dp_owningMol != nullptr, "no owning molecule for bond");
  // TODO: Should this assert that at is from the same molecule?
  dp_dataMol->getBond(d_index).setEndAtomIdx(at->getIdx());
}
Atom* Bond::getBeginAtom() const {
  PRECONDITION(dp_owningMol, "no owning molecule for bond");
  return dp_owningMol->asROMol().getAtomWithIdx(
      dp_dataMol->getBond(d_index).getBeginAtomIdx());
}
Atom* Bond::getEndAtom() const {
  PRECONDITION(dp_owningMol, "no owning molecule for bond");
  return dp_owningMol->asROMol().getAtomWithIdx(
      dp_dataMol->getBond(d_index).getEndAtomIdx());
}
void Bond::setQuery([[maybe_unused]] QUERYBOND_QUERY* what) {
  //  Bonds don't have queries at the moment because I have not
  //  yet figured out what a good base query should be.
  //  It would be nice to be able to do substructure searches
  //  using molecules alone, so it'd be nice if we got this
  //  issue resolved ASAP.
  PRECONDITION(0, "plain bonds have no Query");
}
Bond::QUERYBOND_QUERY* Bond::getQuery() const {
  PRECONDITION(0, "plain bonds have no Query");
  return nullptr;
}
void Bond::expandQuery(
    [[maybe_unused]] QUERYBOND_QUERY* what,
    [[maybe_unused]] Queries::CompositeQueryType how,
    [[maybe_unused]] bool maintainOrder) {
  PRECONDITION(0, "plain bonds have no query");
}
bool Bond::Match(Bond const* what) const {
  bool res;
  if (getBondType() == Bond::UNSPECIFIED ||
      what->getBondType() == Bond::UNSPECIFIED) {
    res = true;
  } else {
    res = getBondType() == what->getBondType();
  }
  return res;
}
void Bond::setBondDir(BondDir what) {
  dp_dataMol->getBond(d_index).setBondDir(
      static_cast<BondEnums::BondDir>(what));
}
Bond::BondDir Bond::getBondDir() const {
  return static_cast<BondDir>(dp_dataMol->getBond(d_index).getBondDir());
}
void Bond::setStereo(BondStereo what) {
  dp_dataMol->getBond(d_index).setStereo(
      static_cast<BondEnums::BondStereo>(what));
}
Bond::BondStereo Bond::getStereo() const {
  return static_cast<BondStereo>(dp_dataMol->getBond(d_index).getStereo());
}
void Bond::setStereoAtoms(unsigned int bgnIdx, unsigned int endIdx) {
  PRECONDITION(
      getOwningMol().getBondBetweenAtoms(getBeginAtomIdx(), bgnIdx) != nullptr,
      "bgnIdx not connected to begin atom of bond");
  PRECONDITION(
      getOwningMol().getBondBetweenAtoms(getEndAtomIdx(), endIdx) != nullptr,
      "endIdx not connected to end atom of bond");

  BondData &bond = dp_dataMol->getBond(d_index);
  bond.stereoAtoms[0] = bgnIdx;
  bond.stereoAtoms[1] = endIdx;
  if (dp_dataMol->hasCompatibilityData()) {
    auto *compatStereo = dp_dataMol->getBondStereoAtomsCompat(d_index);
    CHECK_INVARIANT(compatStereo, "No stereo atoms in compatibility data");
    compatStereo->clear();
    compatStereo->push_back(bgnIdx);
    compatStereo->push_back(endIdx);
  }
}
const INT_VECT& Bond::getStereoAtoms() const {
  return *dp_dataMol->getBondStereoAtomsCompat(d_index);
}
INT_VECT& Bond::getStereoAtoms() {
  return *dp_dataMol->getBondStereoAtomsCompat(d_index);
}
void Bond::updatePropertyCache(bool strict) { (void)strict; }

void Bond::setFlags(std::uint64_t flags) {
  dp_dataMol->getBond(d_index).setFlags(flags);
}
std::uint64_t Bond::getFlags() const {
  return dp_dataMol->getBond(d_index).getFlags();
}
std::uint64_t &Bond::getFlags() {
  return dp_dataMol->getBond(d_index).getFlags();
}

uint8_t getTwiceBondType(Bond::BondType type) {
  return getTwiceBondType(static_cast<BondEnums::BondType>(type));
}

bool Bond::invertChirality() {
  switch (getStereo()) {
    case STEREOATROPCW:
      setStereo(STEREOATROPCCW);
      return true;
    case STEREOATROPCCW:
      setStereo(STEREOATROPCW);
      return true;

    default:
      break;
  }
  return false;
}

uint8_t getTwiceBondType(const Bond &b) {
  return getTwiceBondType(b.getBondType());
}

Atom * Bond::getOtherAtom(Atom const *what) const {
  PRECONDITION(what != nullptr, "null input Atom");
  PRECONDITION(dp_owningMol != nullptr, "no owning molecule for bond");
  return dp_owningMol->getAtomCompat(getOtherAtomIdx(what->getIdx()));
}

STR_VECT Bond::getPropList(bool includePrivate, bool includeComputed) const {
  STR_VECT res = dp_dataMol->getPropList(includePrivate, includeComputed,
                                         RDMol::Scope::BOND, d_index);
  if (includePrivate && includeComputed) {
    res.push_back(detail::computedPropName);
  }
  return res;
}

bool Bond::hasProp(const std::string &key) const{
  PropToken token(key);
  if (token == detail::computedPropNameToken) {
    return true;
  }
  return dp_dataMol->hasBondProp(PropToken(key), getIdx());
}

void Bond::clearProp(const std::string &key) const {
  dp_dataMol->clearSingleBondProp(PropToken(key), getIdx());
}
void Bond::clearComputedProps() const {
  dp_dataMol->clearSingleBondAllProps(getIdx(), true);
}

void Bond::updateProps(const RDProps &source, bool preserveExisting) {
  if (!preserveExisting) {
    clear();
  }
  std::vector<std::string> computedPropList;
  source.getPropIfPresent<std::vector<std::string>>(
      RDKit::detail::computedPropName, computedPropList);
  const STR_VECT keys = source.getPropList();
  for (const auto &key : keys) {
    if (key == RDKit::detail::computedPropName) {
      continue;
    }
    const RDValue &val = source.getPropRDValue(key);
    bool isComputed =
        std::find(computedPropList.begin(), computedPropList.end(), key) !=
        computedPropList.end();
    getDataRDMol().setSingleBondProp(PropToken(key), getIdx(), val, isComputed, true);
  }
}

void Bond::updateProps(const Bond &source, bool preserveExisting) {
  if (!preserveExisting) {
    clear();
  }
  for (auto it = source.dp_dataMol->beginProps(true, RDMol::Scope::BOND,
                                               source.d_index),
            end = source.dp_dataMol->endProps();
       it != end; ++it) {
    dp_dataMol->copySingleProp(*it, d_index, *source.dp_dataMol, *it,
                               source.d_index, RDMol::Scope::BOND);
  }
}

void Bond::clear() {
  dp_dataMol->clearSingleBondAllProps(getIdx(), false);
}

};  // namespace RDKit

namespace {
constexpr const char *bondTypeToString(RDKit::Bond::BondType d) {
  switch (d) {
    case RDKit::Bond::BondType::UNSPECIFIED:
      return "?";
    case RDKit::Bond::BondType::SINGLE:
      return "1";
    case RDKit::Bond::BondType::DOUBLE:
      return "2";
    case RDKit::Bond::BondType::TRIPLE:
      return "3";
    case RDKit::Bond::BondType::QUADRUPLE:
      return "4";
    case RDKit::Bond::BondType::QUINTUPLE:
      return "5";
    case RDKit::Bond::BondType::HEXTUPLE:
      return "6";
    case RDKit::Bond::BondType::ONEANDAHALF:
      return "1.5";
    case RDKit::Bond::BondType::TWOANDAHALF:
      return "2.5";
    case RDKit::Bond::BondType::THREEANDAHALF:
      return "3.5";
    case RDKit::Bond::BondType::FOURANDAHALF:
      return "4.5";
    case RDKit::Bond::BondType::FIVEANDAHALF:
      return "5.5";
    case RDKit::Bond::BondType::AROMATIC:
      return "a";
    case RDKit::Bond::BondType::IONIC:
      return "I";
    case RDKit::Bond::BondType::HYDROGEN:
      return "H";
    case RDKit::Bond::BondType::THREECENTER:
      return "3C";
    case RDKit::Bond::BondType::DATIVEONE:
      return "D1";
    case RDKit::Bond::BondType::DATIVE:
      return "D";
    case RDKit::Bond::BondType::DATIVEL:
      return "DL";
    case RDKit::Bond::BondType::DATIVER:
      return "DR";
    case RDKit::Bond::BondType::OTHER:
      return "Other";
    case RDKit::Bond::BondType::ZERO:
      return "0";
  }
  return ("");
}
constexpr const char *bondDirToString(RDKit::Bond::BondDir d) {
  switch (d) {
    case RDKit::Bond::BondDir::NONE:
      return "NONE";
    case RDKit::Bond::BondDir::BEGINWEDGE:
      return "wedge";
    case RDKit::Bond::BondDir::BEGINDASH:
      return "dash";
    case RDKit::Bond::BondDir::ENDDOWNRIGHT:
      return "\\";
    case RDKit::Bond::BondDir::ENDUPRIGHT:
      return "/";
    case RDKit::Bond::BondDir::EITHERDOUBLE:
      return "x";
    case RDKit::Bond::BondDir::UNKNOWN:
      return "?";
  }
  return ("");
}
constexpr const char *bondStereoToString(RDKit::Bond::BondStereo d) {
  switch (d) {
    case RDKit::Bond::BondStereo::STEREONONE:
      return "NONE";
    case RDKit::Bond::BondStereo::STEREOANY:
      return "ANY";
    case RDKit::Bond::BondStereo::STEREOZ:
      return "Z";
    case RDKit::Bond::BondStereo::STEREOE:
      return "E";
    case RDKit::Bond::BondStereo::STEREOCIS:
      return "CIS";
    case RDKit::Bond::BondStereo::STEREOTRANS:
      return "TRANS";
    case RDKit::Bond::BondStereo::STEREOATROPCW:
      return "CW";
    case RDKit::Bond::BondStereo::STEREOATROPCCW:
      return "CCW";
  }
  return ("");
}
}  // namespace

std::ostream &operator<<(std::ostream &target, const RDKit::Bond &bond) {
  target << bond.getIdx() << " ";
  target << bond.getBeginAtomIdx() << "->" << bond.getEndAtomIdx();
  target << " order: " << bondTypeToString(bond.getBondType());
  if (bond.getBondDir()) {
    target << " dir: " << bondDirToString(bond.getBondDir());
  }
  if (bond.getStereo()) {
    target << " stereo: " << bondStereoToString(bond.getStereo());
    if (bond.getStereoAtoms().size() == 2) {
      const auto &ats = bond.getStereoAtoms();
      target << " ats: (" << ats[0] << " " << ats[1] << ")";
    }
    if (bond.getStereo() == RDKit::Bond::BondStereo::STEREOATROPCCW ||
        bond.getStereo() == RDKit::Bond::BondStereo::STEREOATROPCW) {
      RDKit::Atropisomers::AtropAtomAndBondVec atomAndBonds[2];
      if (RDKit::Atropisomers::getAtropisomerAtomsAndBonds(
              &bond, atomAndBonds, bond.getOwningMol())) {
        target << " bonds: (";
        for (auto i = 0u; i < atomAndBonds[0].second.size(); ++i) {
          if (i) {
            target << " ";
          }
          target << atomAndBonds[0].second[i]->getIdx();
        }
        for (auto i = 0u; i < atomAndBonds[1].second.size(); ++i) {
          target << " " << atomAndBonds[1].second[i]->getIdx();
        }
        target << ")";
      }
    }
  }
  if (bond.getIsConjugated()) {
    target << " conj?: " << bond.getIsConjugated();
  }
  if (bond.getIsAromatic()) {
    target << " aromatic?: " << bond.getIsAromatic();
  }

  return target;
}
