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
#include "MonomerInfo.h"

#include <optional>

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

  atom->setProp(common_properties::dummyLabel, templateName);
  atom->setProp(common_properties::molAtomClass, className);

  return this->addAtom(atom, false, true);
}

void MacroMol::addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                            int fromConnectionPoint, int toConnectionPoint,
                            std::optional<Bond::BondType> bondType) {
  const auto resolvedBondType = bondType.value_or(Bond::BondType::SINGLE);

  auto bondIdx = this->addBond(fromAtomIdx, toAtomIdx, resolvedBondType) - 1;

  auto bond = this->getBondWithIdx(bondIdx);

  this->setBondBookmark(bond, bondIdx);

  bond->setProp(common_properties::_MacroMolToAttachPt, toConnectionPoint);
  bond->setProp(common_properties::_MacroMolFromAttachPt, fromConnectionPoint);
}

}  // namespace RDKit
