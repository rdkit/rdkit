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

namespace RDKit {

std::string monomerClassToString(MonomerClass monomerClass) {
  switch (monomerClass) {
    case MonomerClass::AA:
      return "AA";
    case MonomerClass::NA:
      return "NA";
    case MonomerClass::CHEM:
      return "CHEM";
    case MonomerClass::OTHER:
      return "OTHER";
  }
  return "OTHER";
}

MonomerClass stringToMonomerClass(const std::string &className) {
  if (className == "AA") {
    return MonomerClass::AA;
  } else if (className == "NA") {
    return MonomerClass::NA;
  } else if (className == "CHEM") {
    return MonomerClass::CHEM;
  } else {
    return MonomerClass::OTHER;
  }
}

unsigned int MacroMol::addMacroAtom(MonomerClass monomerClass,
                                    std::string templateName) {
  auto className = monomerClassToString(monomerClass);
  auto atom = new Atom(0);

  atom->setProp(common_properties::isMacroAtom, true);
  atom->setProp(common_properties::dummyLabel, templateName);
  atom->setProp(common_properties::molAtomClass, className);

  return this->addAtom(atom, false, true);
}

unsigned int MacroMol::addMacroBond(unsigned int beginAtomIdx,
                                    unsigned int endAtomIdx, int beginAttachPt,
                                    int endAttachPt,
                                    Bond::BondType bondType) {
  PRECONDITION((this->getAtomWithIdx(beginAtomIdx)
                    ->hasProp(common_properties::isMacroAtom) ||
                this->getAtomWithIdx(endAtomIdx)
                    ->hasProp(common_properties::isMacroAtom)),
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
  PRECONDITION(this->getAtomWithIdx(endMacroAtomIdx)
                   ->hasProp(common_properties::isMacroAtom),
               "end atom is not a macro atom");
  PRECONDITION(!this->getAtomWithIdx(beginAtomIdx)
                    ->hasProp(common_properties::isMacroAtom),
               "begin atom is a macro atom");
  return addMacroBond(beginAtomIdx, endMacroAtomIdx, -1, endAttachPt,
                      bondType);
}

unsigned int MacroMol::addMacroAtomToAtomBond(unsigned int beginMacroAtomIdx,
                                              unsigned int endAtomIdx,
                                              int beginAttachPt,
                                              Bond::BondType bondType) {
  PRECONDITION(this->getAtomWithIdx(beginMacroAtomIdx)
                   ->hasProp(common_properties::isMacroAtom),
               "begin atom is not a macro atom");
  PRECONDITION(!this->getAtomWithIdx(endAtomIdx)
                    ->hasProp(common_properties::isMacroAtom),
               "end atom is a macro atom");
  return addMacroBond(beginMacroAtomIdx, endAtomIdx, beginAttachPt, -1,
                      bondType);
}

unsigned int MacroMol::addBond(unsigned int beginAtomIdx,
                               unsigned int endAtomIdx,
                               Bond::BondType bondType) {
  PRECONDITION(!this->getAtomWithIdx(beginAtomIdx)
                    ->hasProp(common_properties::isMacroAtom),
               "begin atom is a macro atom");
  PRECONDITION(!this->getAtomWithIdx(endAtomIdx)
                    ->hasProp(common_properties::isMacroAtom),
               "end atom is a macro atom");

  return RWMol::addBond(beginAtomIdx, endAtomIdx, bondType);
}

}  // namespace RDKit
