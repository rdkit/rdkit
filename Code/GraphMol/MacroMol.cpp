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

void MacroMol::addMacroBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                            int beginAttachPt, int endAttachPt,
                            bool isDirectional, Bond::BondType bondType) {
  auto bond = this->getBondBetweenAtoms(beginAtomIdx, endAtomIdx);
  if (!bond) {
    auto bondIdx =
        this->addBond(beginAtomIdx, endAtomIdx, Bond::BondType::ZERO) - 1;
    bond = this->getBondWithIdx(bondIdx);
  }
  MacroBondProps props;
  props.beginAttachPt = beginAttachPt;
  props.endAttachPt = endAttachPt;
  props.bondType = bondType;
  props.isDirectional = isDirectional;

  std::vector<MacroBondProps> propsList;
  bond->getPropIfPresent(common_properties::macroBondProps, propsList);
  propsList.push_back(props);
  bond->setProp(common_properties::macroBondProps, propsList);
}

void MacroMol::addAtomToMacroAtomBond(
    unsigned int beginAtomIdx, unsigned int endMacroAtomIdx, int endAttachPt,
    bool isDirectional, Bond::BondType bondType) {
  PRECONDITION(this->getAtomWithIdx(endMacroAtomIdx)
                   ->hasProp(common_properties::isMacroAtom),
               "end atom is not a macro atom");
  PRECONDITION(!this->getAtomWithIdx(beginAtomIdx)
                    ->hasProp(common_properties::isMacroAtom),
               "begin atom is a macro atom");
  addMacroBond(beginAtomIdx, endMacroAtomIdx, -1, endAttachPt, isDirectional,
               bondType);
}

void MacroMol::addMacroAtomToAtomBond(
    unsigned int beginMacroAtomIdx, unsigned int endAtomIdx, int beginAttachPt,
    bool isDirectional, Bond::BondType bondType) {
  PRECONDITION(this->getAtomWithIdx(beginMacroAtomIdx)
                   ->hasProp(common_properties::isMacroAtom),
               "begin atom is not a macro atom");
  PRECONDITION(!this->getAtomWithIdx(endAtomIdx)
                    ->hasProp(common_properties::isMacroAtom),
               "end atom is a macro atom");
  addMacroBond(beginMacroAtomIdx, endAtomIdx, beginAttachPt, -1, isDirectional,
               bondType);
}

void MacroMol::addBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                       int beginAttachPt, int endAttachPt,
                       std::optional<Bond::BondType> bondType) {
  PRECONDITION(!this->getAtomWithIdx(beginAtomIdx)
                    ->hasProp(common_properties::isMacroAtom),
               "begin atom is a macro atom");
  PRECONDITION(!this->getAtomWithIdx(endAtomIdx)
                    ->hasProp(common_properties::isMacroAtom),
               "end atom is a macro atom");

  auto bt = bondType.value_or(Bond::BondType::SINGLE);
  RWMol::addBond(beginAtomIdx, endAtomIdx, bt);
}

}  // namespace RDKit
