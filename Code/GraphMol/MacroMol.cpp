//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define USE_BETTER_ENUMS
#include "MacroMol.h"
#include "Atom.h"

namespace RDKit {
namespace {
bool isMacroAtom(const Atom *atom) {
  return atom->getMacroAtomInfo() != nullptr;
}
}  // namespace

unsigned int MacroMol::addMacroAtom(MonomerClass monomerClass,
                                    std::string symbol) {
  std::string className = monomerClass._to_string();
  auto atom = new Atom(0);

  atom->setMacroAtomInfo(new MacroAtomInfo(symbol, className));

  bool updateLabel = false;
  bool takeOwnership = true;
  return this->addAtom(atom, updateLabel, takeOwnership);
}

unsigned int MacroMol::addMacroBond(unsigned int beginAtomIdx,
                                    unsigned int endAtomIdx, int beginAttachPt,
                                    int endAttachPt, Bond::BondType bondType) {
  PRECONDITION((isMacroAtom(this->getAtomWithIdx(beginAtomIdx)) ||
                isMacroAtom(this->getAtomWithIdx(endAtomIdx))),
               "at least one atom must be a macro atom");
  auto bondIdx = RWMol::addBond(beginAtomIdx, endAtomIdx, bondType) - 1;
  auto bond = this->getBondWithIdx(bondIdx);
  bond->setProp(common_properties::_MacroMolBeginAttachPt, beginAttachPt);
  bond->setProp(common_properties::_MacroMolEndAttachPt, endAttachPt);
  return bondIdx;
}

unsigned int MacroMol::addAtomToMacroAtomBond(unsigned int beginAtomIdx,
                                              unsigned int endMacroAtomIdx,
                                              int endAttachPt,
                                              Bond::BondType bondType) {
  PRECONDITION(isMacroAtom(this->getAtomWithIdx(endMacroAtomIdx)),
               "end atom is not a macro atom");
  PRECONDITION(!isMacroAtom(this->getAtomWithIdx(beginAtomIdx)),
               "begin atom is a macro atom");
  return addMacroBond(beginAtomIdx, endMacroAtomIdx, -1, endAttachPt, bondType);
}

unsigned int MacroMol::addMacroAtomToAtomBond(unsigned int beginMacroAtomIdx,
                                              unsigned int endAtomIdx,
                                              int beginAttachPt,
                                              Bond::BondType bondType) {
  PRECONDITION(isMacroAtom(this->getAtomWithIdx(beginMacroAtomIdx)),
               "begin atom is not a macro atom");
  PRECONDITION(!isMacroAtom(this->getAtomWithIdx(endAtomIdx)),
               "end atom is a macro atom");
  return addMacroBond(beginMacroAtomIdx, endAtomIdx, beginAttachPt, -1,
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

}  // namespace RDKit
