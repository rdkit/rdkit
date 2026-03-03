//
// Copyright (C) 2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "DummyAtom.h"

#include <stdexcept>
#include <string>

#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>
#include <RDGeneral/types.h>

namespace RDKit {

std::unique_ptr<Atom> createDummyAtom() {
  auto atom = std::make_unique<QueryAtom>();
  atom->setAtomicNum(dummyAtomicNum);
  // This prevents the atom from being shown in SMARTS as [#0]
  atom->setQuery(makeAtomNullQuery());
  // This prevents the atom from having a VAL in MDL MOL file output
  atom->setNoImplicit(false);
  return atom;
}

bool isAttachmentPointDummy(const Atom &atom) {
  std::string label;
  return atom.getAtomicNum() == dummyAtomicNum && atom.getDegree() == 1 &&
         atom.getPropIfPresent(common_properties::atomLabel, label) &&
         label.find(attachmentPointLabelPrefix) == 0;
}

std::unique_ptr<Atom> makeNewRGroup(unsigned int rGroupNum) {
  if (rGroupNum == 0) {
    throw std::invalid_argument("R-groups cannot have an index of 0.");
  }
  auto atom = createDummyAtom();
  // Set atomLabel as _RX and dummyLabel as RX
  auto rLabel = std::string(rGroupLabelPrefix) + std::to_string(rGroupNum);
  auto dLabel = rLabel.substr(1);
  atom->setProp(common_properties::atomLabel, rLabel);
  atom->setProp(common_properties::dummyLabel, dLabel);
  atom->setProp(common_properties::_MolFileRLabel, rGroupNum);
  atom->setIsotope(rGroupNum);
  return atom;
}

std::optional<unsigned int> getRGroupNumber(const Atom *atom) {
  unsigned int rGroupNum = 0;
  if (!(atom->getAtomicNum() == dummyAtomicNum &&
        atom->getPropIfPresent(common_properties::_MolFileRLabel, rGroupNum))) {
    return std::nullopt;
  }
  if (rGroupNum == 0) {
    throw std::invalid_argument("R-groups cannot have an index of 0.");
  }
  return rGroupNum;
}

}  // namespace RDKit
