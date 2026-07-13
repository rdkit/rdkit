//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MacroMol.h"
#include "Atom.h"

#include <RDGeneral/Invariant.h>

#include <utility>

namespace RDKit {
namespace {
bool isMacroAtom(const Atom *atom) {
  return atom->getMacroAtomInfo() != nullptr;
}
}  // namespace

unsigned int MacroMol::addMacroAtom(MonomerClass monomerClass,
                                    std::string symbol) {
  auto atom = new Atom(0);

  atom->setMacroAtomInfo(new MacroAtomInfo(std::move(symbol), monomerClass));

  bool updateLabel = false;
  bool takeOwnership = true;
  return this->addAtom(atom, updateLabel, takeOwnership);
}

unsigned int MacroMol::addMacroBondHelper(unsigned int beginAtomIdx,
                                          unsigned int endAtomIdx,
                                          int beginAttachPt, int endAttachPt,
                                          Bond::BondType bondType) {
  auto bond = this->getBondBetweenAtoms(beginAtomIdx, endAtomIdx);
  if (!bond) {
    // The actual macro bond types live in MacroBondInfo.
    auto numBonds =
        RWMol::addBond(beginAtomIdx, endAtomIdx, Bond::BondType::UNSPECIFIED);
    auto bondIdx = numBonds - 1;
    bond = this->getBondWithIdx(bondIdx);
    bond->setMacroBondInfo(new MacroBondInfo(
        beginAttachPt, endAttachPt, static_cast<unsigned int>(bondType)));
    return numBonds;
  }
  bond->setBondType(Bond::BondType::UNSPECIFIED);
  auto *info = bond->getMacroBondInfo();
  if (!info) {
    bond->setMacroBondInfo(new MacroBondInfo());
    info = bond->getMacroBondInfo();
  }
  if (bond->getBeginAtomIdx() == beginAtomIdx) {
    info->addBond(beginAttachPt, endAttachPt,
                  static_cast<unsigned int>(bondType));
  } else {
    info->addBond(endAttachPt, beginAttachPt,
                  static_cast<unsigned int>(bondType));
  }
  auto numBonds = this->getNumBonds();
  return numBonds;
}

unsigned int MacroMol::addMacroBond(unsigned int beginAtomIdx,
                                    unsigned int endAtomIdx, int beginAttachPt,
                                    int endAttachPt, Bond::BondType bondType) {
  PRECONDITION(isMacroAtom(this->getAtomWithIdx(beginAtomIdx)),
               "begin atom is not a macro atom");
  PRECONDITION(isMacroAtom(this->getAtomWithIdx(endAtomIdx)),
               "end atom is not a macro atom");
  return addMacroBondHelper(beginAtomIdx, endAtomIdx, beginAttachPt,
                            endAttachPt, bondType);
}

unsigned int MacroMol::addAtomToMacroAtomBond(unsigned int beginAtomIdx,
                                              unsigned int endMacroAtomIdx,
                                              int endAttachPt,
                                              Bond::BondType bondType) {
  PRECONDITION(isMacroAtom(this->getAtomWithIdx(endMacroAtomIdx)),
               "end atom is not a macro atom");
  PRECONDITION(!isMacroAtom(this->getAtomWithIdx(beginAtomIdx)),
               "begin atom is a macro atom");
  return addMacroBondHelper(beginAtomIdx, endMacroAtomIdx, -1, endAttachPt,
                            bondType);
}

unsigned int MacroMol::addMacroAtomToAtomBond(unsigned int beginMacroAtomIdx,
                                              unsigned int endAtomIdx,
                                              int beginAttachPt,
                                              Bond::BondType bondType) {
  PRECONDITION(isMacroAtom(this->getAtomWithIdx(beginMacroAtomIdx)),
               "begin atom is not a macro atom");
  PRECONDITION(!isMacroAtom(this->getAtomWithIdx(endAtomIdx)),
               "end atom is a macro atom");
  return addMacroBondHelper(beginMacroAtomIdx, endAtomIdx, beginAttachPt, -1,
                            bondType);
}

unsigned int MacroMol::addBond(unsigned int beginAtomIdx,
                               unsigned int endAtomIdx,
                               Bond::BondType bondType) {
  PRECONDITION(!isMacroAtom(this->getAtomWithIdx(beginAtomIdx)),
               "begin atom is a macro atom");
  PRECONDITION(!isMacroAtom(this->getAtomWithIdx(endAtomIdx)),
               "end atom is a macro atom");

  return RWMol::addBond(beginAtomIdx, endAtomIdx, bondType);
}

unsigned int MacroMol::addBond(Atom *beginAtom, Atom *endAtom,
                               Bond::BondType bondType) {
  PRECONDITION(beginAtom && endAtom, "NULL atom passed in");
  return addBond(beginAtom->getIdx(), endAtom->getIdx(), bondType);
}

unsigned int MacroMol::addBond(Bond *bond, bool takeOwnership) {
  PRECONDITION(bond, "NULL bond passed in");
  PRECONDITION(!isMacroAtom(this->getAtomWithIdx(bond->getBeginAtomIdx())),
               "begin atom is a macro atom");
  PRECONDITION(!isMacroAtom(this->getAtomWithIdx(bond->getEndAtomIdx())),
               "end atom is a macro atom");
  return RWMol::addBond(bond, takeOwnership);
}

}  // namespace RDKit
