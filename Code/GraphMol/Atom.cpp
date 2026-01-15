//
//  Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cmath>

#include "ROMol.h"
#include "RDMol.h"
#include "Atom.h"
#include "PeriodicTable.h"
#include "SanitException.h"
#include "QueryOps.h"
#include "MonomerInfo.h"

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Dict.h>
#include <GraphMol/rdmol_throw.h>

namespace RDKit {

bool isAromaticAtom(ConstRDMolAtom atom) {
  if (atom.data().getIsAromatic()) {
    return true;
  }
  auto [begin, end] = atom.mol().getAtomBonds(atom.index());
  for (; begin != end; ++begin) {
    const auto &bond = atom.mol().getBond(*begin);
    if (bond.getIsAromatic() ||
        bond.getBondType() == BondEnums::BondType::AROMATIC) {
      return true;
    }
  }
  return false;
}
bool isAromaticAtom(const Atom &atom) {
  if (atom.getIsAromatic()) {
    return true;
  }
  if (atom.hasOwningMol()) {
    for (const auto &bond : atom.getOwningMol().atomBonds(&atom)) {
      if (bond->getIsAromatic() ||
          bond->getBondType() == Bond::BondType::AROMATIC) {
        return true;
      }
    }
  }
  return false;
}

unsigned int getEffectiveAtomicNum(const Atom &atom, bool checkValue) {
  auto effectiveAtomicNum = atom.getAtomicNum() - atom.getFormalCharge();
  if (checkValue &&
      (effectiveAtomicNum < 0 ||
       effectiveAtomicNum >
           static_cast<int>(PeriodicTable::getTable()->getMaxAtomicNumber()))) {
    throw AtomValenceException("Effective atomic number out of range",
                               atom.getIdx());
  }
  effectiveAtomicNum = std::clamp(
      effectiveAtomicNum, 0,
      static_cast<int>(PeriodicTable::getTable()->getMaxAtomicNumber()));
  return static_cast<unsigned int>(effectiveAtomicNum);
}

// Determine whether or not an element is to the left of carbon.
bool isEarlyAtom(int atomicNum) {
  static const bool table[119] = {
      false,  // #0 *
      false,  // #1 H
      false,  // #2 He
      true,   // #3 Li
      true,   // #4 Be
      true,   // #5 B
      false,  // #6 C
      false,  // #7 N
      false,  // #8 O
      false,  // #9 F
      false,  // #10 Ne
      true,   // #11 Na
      true,   // #12 Mg
      true,   // #13 Al
      false,  // #14 Si
      false,  // #15 P
      false,  // #16 S
      false,  // #17 Cl
      false,  // #18 Ar
      true,   // #19 K
      true,   // #20 Ca
      true,   // #21 Sc
      true,   // #22 Ti
      false,  // #23 V
      false,  // #24 Cr
      false,  // #25 Mn
      false,  // #26 Fe
      false,  // #27 Co
      false,  // #28 Ni
      false,  // #29 Cu
      true,   // #30 Zn
      true,   // #31 Ga
      true,   // #32 Ge  see github #2606
      false,  // #33 As
      false,  // #34 Se
      false,  // #35 Br
      false,  // #36 Kr
      true,   // #37 Rb
      true,   // #38 Sr
      true,   // #39 Y
      true,   // #40 Zr
      true,   // #41 Nb
      false,  // #42 Mo
      false,  // #43 Tc
      false,  // #44 Ru
      false,  // #45 Rh
      false,  // #46 Pd
      false,  // #47 Ag
      true,   // #48 Cd
      true,   // #49 In
      true,   // #50 Sn  see github #2606
      true,   // #51 Sb  see github #2775
      false,  // #52 Te
      false,  // #53 I
      false,  // #54 Xe
      true,   // #55 Cs
      true,   // #56 Ba
      true,   // #57 La
      true,   // #58 Ce
      true,   // #59 Pr
      true,   // #60 Nd
      true,   // #61 Pm
      false,  // #62 Sm
      false,  // #63 Eu
      false,  // #64 Gd
      false,  // #65 Tb
      false,  // #66 Dy
      false,  // #67 Ho
      false,  // #68 Er
      false,  // #69 Tm
      false,  // #70 Yb
      false,  // #71 Lu
      true,   // #72 Hf
      true,   // #73 Ta
      false,  // #74 W
      false,  // #75 Re
      false,  // #76 Os
      false,  // #77 Ir
      false,  // #78 Pt
      false,  // #79 Au
      true,   // #80 Hg
      true,   // #81 Tl
      true,   // #82 Pb  see github #2606
      true,   // #83 Bi  see github #2775
      false,  // #84 Po
      false,  // #85 At
      false,  // #86 Rn
      true,   // #87 Fr
      true,   // #88 Ra
      true,   // #89 Ac
      true,   // #90 Th
      true,   // #91 Pa
      true,   // #92 U
      true,   // #93 Np
      false,  // #94 Pu
      false,  // #95 Am
      false,  // #96 Cm
      false,  // #97 Bk
      false,  // #98 Cf
      false,  // #99 Es
      false,  // #100 Fm
      false,  // #101 Md
      false,  // #102 No
      false,  // #103 Lr
      true,   // #104 Rf
      true,   // #105 Db
      true,   // #106 Sg
      true,   // #107 Bh
      true,   // #108 Hs
      true,   // #109 Mt
      true,   // #110 Ds
      true,   // #111 Rg
      true,   // #112 Cn
      true,   // #113 Nh
      true,   // #114 Fl
      true,   // #115 Mc
      true,   // #116 Lv
      true,   // #117 Ts
      true,   // #118 Og
  };
  return ((unsigned int)atomicNum < 119) && table[atomicNum];
}

Atom::Atom() {
  dp_dataMol = new RDMol();
  dp_dataMol->addAtom();
  dp_owningMol = nullptr;
}

Atom::Atom(unsigned int num): Atom() { setAtomicNum(num); }

Atom::Atom(const std::string &what) : Atom() {
  setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(what));
};

Atom::Atom(const Atom &other): Atom() {
  initFromOther(other);
}

Atom &Atom::operator=(const Atom &other) {
  if (this == &other) {
    return *this;
  }
  initFromOther(other);
  return *this;
}

Atom::Atom(Atom &&other) {
  dp_dataMol = other.dp_dataMol;
  dp_owningMol = other.dp_owningMol;
  d_index = other.d_index;
  other.dp_dataMol = nullptr;
  other.dp_owningMol = nullptr;
  other.d_index = 0;
}

Atom &Atom::operator=(Atom &&other) {
  if (this == &other) {
    return *this;
  }
  if (dp_dataMol != dp_owningMol) {
    delete dp_dataMol;
  }
  dp_dataMol = other.dp_dataMol;
  dp_owningMol = other.dp_owningMol;
  d_index = other.d_index;
  other.dp_dataMol = nullptr;
  other.dp_owningMol = nullptr;
  return *this;
}

void Atom::initFromOther(const Atom &other, const bool preserveProps) {
  dp_dataMol->getAtom(d_index) = other.dp_dataMol->getAtom(other.d_index);
  if (!preserveProps) {
    dp_dataMol->clearSingleAtomAllProps(d_index);
    for (auto it = other.dp_dataMol->beginProps(true, RDMol::Scope::ATOM,
                                                other.d_index),
              end = other.dp_dataMol->endProps();
         it != end; ++it) {
      dp_dataMol->copySingleProp(it->name(), d_index, *other.dp_dataMol, it->name(),
                                 other.d_index, RDMol::Scope::ATOM);
    }
  }
  const auto &otherMonomerInfo = other.dp_dataMol->monomerInfo;
  if (auto it = otherMonomerInfo.find(other.d_index); it != otherMonomerInfo.end()) {
    AtomMonomerInfo* monomerInfo = it->second.get();
    auto monomerInfoCopy = std::unique_ptr<AtomMonomerInfo>(monomerInfo->copy());
    dp_dataMol->monomerInfo[d_index] = std::move(monomerInfoCopy);
  } else {
    // If the atom doesn't have monomer info, remove it from the copy. Useful in
    // the case of replacement, where atomMonomerInfo should now not exist.
    dp_dataMol->monomerInfo.erase(d_index);
  }
}

Atom::~Atom() {
  if (dp_dataMol != dp_owningMol) {
    delete dp_dataMol;
  }
}

Atom *Atom::copy() const {
  auto *res = new Atom(*this);
  return res;
}

int Atom::getAtomicNum() const {
  return dp_dataMol->getAtom(d_index).getAtomicNum();
}
void Atom::setAtomicNum(int newNum) {
  dp_dataMol->getAtom(d_index).setAtomicNum(newNum);
}
std::string Atom::getSymbol() const {
  int atomicNum = getAtomicNum();
  std::string res;
  // handle dummies differently:
  if (atomicNum != 0 ||
      !dp_dataMol->getAtomPropIfPresent<std::string>(
          common_properties::dummyLabelToken, d_index, res)) {
    res = PeriodicTable::getTable()->getElementSymbol(atomicNum);
  }
  return res;
}
ROMol &Atom::getOwningMol() const {
  PRECONDITION(dp_owningMol != nullptr, "no owner");
  return dp_owningMol->asROMol();
}
unsigned int Atom::getDegree() const {
  return dp_owningMol ? dp_owningMol->getAtomDegree(d_index) : 0;
}
unsigned int Atom::getTotalDegree() const {
  return getTotalNumHs(false) + getDegree();
}
unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
  // Support for isolated Atoms by using data mol, instead of owning mol.
  return dp_dataMol->getTotalNumHs(d_index, includeNeighbors);
}
unsigned int Atom::getTotalValence() const {
  return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
}
unsigned int Atom::getNumImplicitHs() const {
  return dp_dataMol->getAtom(d_index).getNumImplicitHs();
}

unsigned int Atom::getValence(ValenceType which) const {
  if (!dp_owningMol) {
    return 0;
  }
  const uint32_t valence = dp_dataMol->getAtom(d_index).getValence(which);
  return valence == AtomData::unsetValenceVal ? (unsigned int)-1 : valence;
}

int Atom::getExplicitValence() const {
  const uint32_t explicitValence =
      dp_dataMol->getAtom(d_index).getExplicitValence();
  return explicitValence == AtomData::unsetValenceVal ? -1 : explicitValence;
}
int Atom::getImplicitValence() const {
  // RDMol returns a uint32_t, but we want to return -1 for the default value
  const uint32_t implicitValence =
      dp_dataMol->getAtom(d_index).getImplicitValence();
  return implicitValence == AtomData::unsetValenceVal ? -1 : implicitValence;
}
unsigned int Atom::getNumRadicalElectrons() const {
  return dp_dataMol->getAtom(d_index).getNumRadicalElectrons();
}
void Atom::setNumRadicalElectrons(unsigned int num) {
  dp_dataMol->getAtom(d_index).setNumRadicalElectrons(num);
}
int Atom::getFormalCharge() const {
  return dp_dataMol->getAtom(d_index).getFormalCharge();
}
void Atom::setFormalCharge(int what) {
  dp_dataMol->getAtom(d_index).setFormalCharge(what);
}
void Atom::setNoImplicit(bool what) {
  dp_dataMol->getAtom(d_index).setNoImplicit(what);
}
bool Atom::getNoImplicit() const {
  return dp_dataMol->getAtom(d_index).getNoImplicit();
}
void Atom::setNumExplicitHs(unsigned int what) {
  dp_dataMol->getAtom(d_index).setNumExplicitHs(what);
}
unsigned int Atom::getNumExplicitHs() const {
  return dp_dataMol->getAtom(d_index).getNumExplicitHs();
}
void Atom::setIsAromatic(bool what) {
  dp_dataMol->getAtom(d_index).setIsAromatic(what);
}
bool Atom::getIsAromatic() const {
  return dp_dataMol->getAtom(d_index).getIsAromatic();
}
double Atom::getMass() const {
  return dp_dataMol->getAtom(d_index).getMass();
}
void Atom::setIsotope(unsigned int what) {
  dp_dataMol->getAtom(d_index).setIsotope(what);
}
unsigned int Atom::getIsotope() const {
  return dp_dataMol->getAtom(d_index).getIsotope();
}
void Atom::setChiralTag(ChiralType what) {
  dp_dataMol->getAtom(d_index).setChiralTag(
      static_cast<RDKit::AtomEnums::ChiralType>(what));
}
bool Atom::invertChirality() {
  return dp_dataMol->invertAtomChirality(d_index);
}
Atom::ChiralType Atom::getChiralTag() const {
  return static_cast<ChiralType>(dp_dataMol->getAtom(d_index).getChiralTag());
}
void Atom::setHybridization(HybridizationType what) {
  dp_dataMol->getAtom(d_index).setHybridization(
      static_cast<RDKit::AtomEnums::HybridizationType>(what));
}
Atom::HybridizationType Atom::getHybridization() const {
  return static_cast<HybridizationType>(
      dp_dataMol->getAtom(d_index).getHybridization());
}

void Atom::setQuery([[maybe_unused]] QUERYATOM_QUERY* what) {
  //  Atoms don't have complex queries so this has to fail
  PRECONDITION(0, "plain atoms have no Query");
}
Atom::QUERYATOM_QUERY* Atom::getQuery() const { return nullptr; }
void Atom::expandQuery(
    [[maybe_unused]] QUERYATOM_QUERY* what,
    [[maybe_unused]] Queries::CompositeQueryType how,
    [[maybe_unused]] bool maintainOrder) {
  PRECONDITION(0, "plain atoms have no Query");
}
bool Atom::Match(Atom const* what) const {
  PRECONDITION(what, "bad query atom");
  bool res = getAtomicNum() == what->getAtomicNum();

  // special dummy--dummy match case:
  //   [*] matches [*],[1*],[2*],etc.
  //   [1*] only matches [*] and [1*]
  if (res) {
    // if (hasOwningMol() && what->hasOwningMol() &&
        // this->getOwningMol().getRingInfo()->isInitialized() &&
        // what->getOwningMol().getRingInfo()->isInitialized() &&
        // this->getOwningMol().getRingInfo()->numAtomRings(d_index) >
            // what->getOwningMol().getRingInfo()->numAtomRings(what->d_index)) {
      // res = false;
  // } else if (!this->getAtomicNum()) {
  if (!this->getAtomicNum()) {
      // this is the new behavior, based on the isotopes:
      int tgt = this->getIsotope();
      int test = what->getIsotope();
      if (tgt && test && tgt != test) {
        res = false;
      }
    } else {
      // standard atom-atom match: The general rule here is that if this atom
      // has a property that
      // deviates from the default, then the other atom should match that value.
      if ((this->getFormalCharge() &&
           this->getFormalCharge() != what->getFormalCharge()) ||
          (this->getIsotope() && this->getIsotope() != what->getIsotope()) ||
          (this->getNumRadicalElectrons() &&
           this->getNumRadicalElectrons() != what->getNumRadicalElectrons())) {
        res = false;
      }
    }
  }
  return res;
}
int Atom::getPerturbationOrder(const INT_LIST &probe) const {
  PRECONDITION(
      dp_owningMol != nullptr,
      "perturbation order not defined for atoms not associated with molecules");
  const size_t numBonds = dp_owningMol->getAtomDegree(d_index);
  PRECONDITION(numBonds == probe.size(), "size mismatch");

  auto [beginBonds, endBonds] = dp_owningMol->getAtomBonds(d_index);

  std::vector<int> copy(probe.begin(), probe.end());

  // The bond list on the atom should be sorted, but check just in case.
  bool isSorted = true;
  for (size_t i = 0; numBonds > 0 && i < numBonds - 1; ++i) {
    isSorted = isSorted && (beginBonds[i] < beginBonds[i + 1]);
  }

  if (isSorted) {
    // This is effectively selection sort, which minimizes the number of swaps,
    // so we can just count the number of swaps here.
    int nSwaps = 0;
    for (size_t desti = 0, n = copy.size(); desti < n; ++desti) {
      size_t besti = desti;
      for (size_t srci = desti + 1; srci < n; ++srci) {
        if (copy[srci] < copy[besti]) {
          besti = srci;
        }
      }
      if (besti != desti) {
        ++nSwaps;
        std::swap(copy[besti], copy[desti]);
      }
    }
    return nSwaps;
  }

  std::vector<int> bondsCopy(beginBonds, endBonds);
  int nSwaps = static_cast<int>(countSwapsToInterconvert(bondsCopy, copy));
  return nSwaps;
}
void Atom::updatePropertyCache(bool strict) {
  calcExplicitValence(strict);
  calcImplicitValence(strict);
}

bool Atom::needsUpdatePropertyCache() const {
  return dp_dataMol->getAtom(d_index).needsUpdatePropertyCache();
}
void Atom::clearPropertyCache() {
  dp_dataMol->getAtom(d_index).clearPropertyCache();
}
int Atom::calcExplicitValence(bool strict) {
  PRECONDITION(dp_owningMol,
               "calcExplicitValence requires an owning molecule");
  return dp_owningMol->calcExplicitValence(d_index, strict);
}
int Atom::calcImplicitValence(bool strict) {
  PRECONDITION(dp_owningMol,
               "calcImplicitValence requires an owning molecule");
  return dp_owningMol->calcImplicitValence(d_index, strict);
}

bool Atom::hasValenceViolation() const {
  PRECONDITION(dp_owningMol, "hasValenceViolation requires an owning molecule");
  return dp_owningMol->hasValenceViolation(d_index);
}

int Atom::getExplicitValencePrivate() const {
  PRECONDITION(dp_dataMol,
               "valence not defined for atoms not associated with molecules");
  return dp_dataMol->getAtom(d_index).explicitValence;
}
int Atom::getImplicitValencePrivate() const {
  PRECONDITION(dp_dataMol,
               "valence not defined for atoms not associated with molecules");
  return dp_dataMol->getAtom(d_index).implicitValence;
}

void setAtomRLabel(Atom *atom, int rlabel) {
  PRECONDITION(atom, "bad atom");
  // rlabel ==> n2 => 0..99
  PRECONDITION(rlabel >= 0 && rlabel < 100,
               "rlabel out of range for MDL files");
  // Default of zero indicates no RLabel
  if (rlabel == 0) {
    atom->getRDMol().clearSingleAtomProp(common_properties::_MolFileRLabelToken,
                                         atom->getIdx());
  } else {
    atom->getRDMol().setSingleAtomProp(common_properties::_MolFileRLabelToken,
                                       atom->getIdx(), uint32_t(rlabel));
  }
}
//! Gets the atom's RLabel
int getAtomRLabel(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  const uint32_t *prop =
      atom->getDataRDMol().getAtomPropArrayIfPresent<uint32_t>(
      common_properties::_MolFileRLabelToken);
  return (prop != nullptr) ? prop[atom->getIdx()] : 0;
}

void setAtomStringProp(Atom* atom, const PropToken& token, const std::string& value) {
  PRECONDITION(atom, "bad atom");
  if (value.empty()) {
    atom->getDataRDMol().clearSingleAtomProp(token, atom->getIdx());
  } else {
    atom->getDataRDMol().setSingleAtomProp<PropToken>(
        token, atom->getIdx(), PropToken(value), false, true);
  }
}

std::string getAtomStringProp(const Atom* atom, const PropToken& token) {
  PRECONDITION(atom, "bad atom");
  PropToken prop;
  bool found = atom->getDataRDMol().getAtomPropIfPresent<PropToken>(
      token, atom->getIdx(), prop);
  return found ? prop.getString() : std::string();
}

void setAtomAlias(Atom *atom, const std::string &alias) {
  setAtomStringProp(atom, common_properties::molFileAliasToken, alias);
}

std::string getAtomAlias(const Atom *atom) {
  return getAtomStringProp(atom, common_properties::molFileAliasToken);
}

void setAtomValue(Atom *atom, const std::string &value) {
  setAtomStringProp(atom, common_properties::molFileValueToken, value);
}

std::string getAtomValue(const Atom *atom) {
  return getAtomStringProp(atom, common_properties::molFileValueToken);
}

void setSupplementalSmilesLabel(Atom *atom, const std::string &label) {
  setAtomStringProp(atom, common_properties::_supplementalSmilesLabelToken, label);
}

std::string getSupplementalSmilesLabel(const Atom *atom) {
  return getAtomStringProp(atom, common_properties::_supplementalSmilesLabelToken);
}

unsigned int numPiElectrons(const Atom &atom) {
  unsigned int res = 0;
  if (atom.getIsAromatic()) {
    res = 1;
  } else if (atom.getHybridization() != Atom::SP3) {
    auto val =
        static_cast<unsigned int>(atom.getValence(Atom::ValenceType::EXPLICIT));
    unsigned int physical_bonds = atom.getNumExplicitHs();
    const auto &mol = atom.getOwningMol();
    for (const auto bond : mol.atomBonds(&atom)) {
      if (bond->getValenceContrib(&atom) != 0.0) {
        ++physical_bonds;
      }
    }
    CHECK_INVARIANT(val >= physical_bonds,
                    "explicit valence exceeds atom degree");
    res = val - physical_bonds;
  }
  return res;
}
}  // namespace RDKit

namespace {
constexpr const char *hybridizationToString(
    RDKit::Atom::HybridizationType type) {
  switch (type) {
    case RDKit::Atom::HybridizationType::UNSPECIFIED:
      return "";
    case RDKit::Atom::HybridizationType::S:
      return "S";
    case RDKit::Atom::HybridizationType::SP:
      return "SP";
    case RDKit::Atom::HybridizationType::SP2:
      return "SP2";
    case RDKit::Atom::HybridizationType::SP3:
      return "SP3";
    case RDKit::Atom::HybridizationType::SP3D:
      return "SP3D";
    case RDKit::Atom::HybridizationType::SP2D:
      return "SP2D";
    case RDKit::Atom::HybridizationType::SP3D2:
      return "SP3D2";
    case RDKit::Atom::HybridizationType::OTHER:
      return "OTHER";
  }
  return "";
}
constexpr const char *chiralityToString(RDKit::Atom::ChiralType type) {
  switch (type) {
    case RDKit::Atom::ChiralType::CHI_UNSPECIFIED:
      return "Unspecified";
    case RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CW:
      return "CW";
    case RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
      return "CCW";
    case RDKit::Atom::ChiralType::CHI_OTHER:
      return "Other";
    case RDKit::Atom::ChiralType::CHI_TETRAHEDRAL:
      return "Td";
    case RDKit::Atom::ChiralType::CHI_ALLENE:
      return "Allene";
    case RDKit::Atom::ChiralType::CHI_SQUAREPLANAR:
      return "SqP";
    case RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
      return "Tbp";
    case RDKit::Atom::ChiralType::CHI_OCTAHEDRAL:
      return "Oh";
  }
  return "";
}
}  // namespace
std::ostream &operator<<(std::ostream &target, const RDKit::Atom &at) {
  target << at.getIdx() << " " << at.getAtomicNum() << " " << at.getSymbol();
  target << " chg: " << at.getFormalCharge();
  target << "  deg: " << at.getDegree();
  target << " exp: ";
  target << (at.getExplicitValencePrivate() >= 0
                 ? std::to_string(at.getExplicitValencePrivate())
                                       : "N/A");

  target << " imp: ";
  if (at.getNoImplicit()) {
    target << "0";
  } else {
    target << (at.getImplicitValencePrivate() >= 0
                   ? std::to_string(at.getImplicitValencePrivate())
                                         : "N/A");
  }
  target << " hyb: " << hybridizationToString(at.getHybridization());
  if (at.getIsAromatic()) {
    target << " arom?: " << at.getIsAromatic();
  }
  if (at.getChiralTag() != RDKit::Atom::CHI_UNSPECIFIED) {
    target << " chi: " << chiralityToString(at.getChiralTag());
    int perm;
    if (at.getPropIfPresent(RDKit::common_properties::_chiralPermutation,
                            perm)) {
      target << "(" << perm << ")";
    }
    target << " nbrs:[";
    bool first = true;
    for (const auto nbr : at.getOwningMol().atomNeighbors(&at)) {
      if (!first) {
        target << " ";
      } else {
        first = false;
      }
      target << nbr->getIdx();
    }
    target << "]";
  }
  if (at.getNumRadicalElectrons()) {
    target << " rad: " << at.getNumRadicalElectrons();
  }
  if (at.getIsotope()) {
    target << " iso: " << at.getIsotope();
  }
  if (at.getAtomMapNum()) {
    target << " mapno: " << at.getAtomMapNum();
  }
  if (at.hasQuery()) {
    target << " query: " << at.getQuery()->getDescription();
  }
  return target;
};

namespace RDKit {

AtomMonomerInfo * Atom::getMonomerInfo() {
  auto found = dp_dataMol->monomerInfo.find(d_index);
  if (found != dp_dataMol->monomerInfo.end()) {
    return found->second.get();
  }
  return nullptr;
}

const AtomMonomerInfo * Atom::getMonomerInfo() const {
  const auto found = dp_dataMol->monomerInfo.find(d_index);
  if (found != dp_dataMol->monomerInfo.end()) {
    return found->second.get();
  }
  return nullptr;
}

void Atom::setMonomerInfo(AtomMonomerInfo *info) {
  dp_dataMol->monomerInfo[d_index].reset(info);
}

void Atom::setAtomMapNum(int mapno, bool strict) {
  PRECONDITION(
      !strict || (mapno >= 0 && mapno < 1000),
      "atom map number out of range [0..1000], use strict=false to override");
  const uint32_t idx = getIdx();
  if (mapno != 0) {
    dp_dataMol->setSingleAtomProp(common_properties::molAtomMapNumberToken, idx, mapno, false, true);
  } else if (dp_dataMol->hasAtomProp(common_properties::molAtomMapNumberToken, idx)) {
    dp_dataMol->clearSingleAtomProp(common_properties::molAtomMapNumberToken, idx);
  }
}

int Atom::getAtomMapNum() const {
  int mapno = 0;
  dp_dataMol->getAtomPropIfPresent(common_properties::molAtomMapNumberToken,
                                   getIdx(), mapno);
  return mapno;
}

void Atom::setFlags(std::uint64_t flags) {
  dp_dataMol->getAtom(d_index).setFlags(flags);
}
std::uint64_t Atom::getFlags() const {
  return dp_dataMol->getAtom(d_index).getFlags();
}
std::uint64_t &Atom::getFlags() {
  return dp_dataMol->getAtom(d_index).getFlags();
}

STR_VECT Atom::getPropList(bool includePrivate, bool includeComputed) const {
  STR_VECT res = dp_dataMol->getPropList(includePrivate, includeComputed,
                                         RDMol::Scope::ATOM, d_index);
  if (includePrivate && includeComputed) {
    // Only include __computedProps if there is a computed prop
    auto begin = dp_dataMol->beginProps(true, RDMol::Scope::ATOM, d_index);
    auto end = dp_dataMol->endProps();
    for (; begin != end; ++begin) {
      if (begin->isComputed()) {
        res.push_back(std::string(detail::computedPropName));
        break;
      }
    }
  }
  return res;
}


bool Atom::hasProp(const std::string_view &key) const {
  PropToken token(key);
  if (token == detail::computedPropNameToken) {
    return true;
  }
  return dp_dataMol->hasAtomProp(token, getIdx());
}

void Atom::clearProp(const std::string_view &key) const {
  dp_dataMol->clearSingleAtomProp(PropToken(key), getIdx());
}

void Atom::clearComputedProps() const {
  dp_dataMol->clearSingleAtomAllProps(getIdx(), true);
}

void Atom::updateProps(const RDProps &source, bool preserveExisting) {
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
    getDataRDMol().setSingleAtomProp(PropToken(key), getIdx(), val, isComputed, true);
  }
}

void Atom::updateProps(const Atom &source, bool preserveExisting) {
  if (!preserveExisting) {
    clear();
  }
  for (auto it = source.dp_dataMol->beginProps(true, RDMol::Scope::ATOM,
                                               source.d_index),
            end = source.dp_dataMol->endProps();
       it != end; ++it) {
    dp_dataMol->copySingleProp(it->name(), d_index, *source.dp_dataMol, it->name(),
                               source.d_index, RDMol::Scope::ATOM);
  }
}

void Atom::clear() { dp_dataMol->clearSingleAtomAllProps(getIdx(), false); }

void Atom::setOwningMol(ROMol *other) {
  setOwningMol(&other->asRDMol());
}

void Atom::setOwningMol(RDKit::RDMol *other) {
  PRECONDITION(dp_owningMol == nullptr || dp_owningMol == other,
               "setOwningMol called twice");
  dp_owningMol = other;
}

}  // namespace RDKit
