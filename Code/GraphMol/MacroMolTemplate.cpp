//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MacroMolTemplate.h"

#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>

#include <algorithm>
#include <set>
#include <unordered_set>

namespace RDKit {

namespace {
const std::string SUP_TYPE = "SUP";
const std::string LGRP_CLASS = "LGRP";

[[noreturn]] void invalidTemplate(const std::string &message) {
  // Give all builder validation failures a recognizable prefix.
  throw ValueErrorException("invalid MacroMolTemplate: " + message);
}

void validateUniqueInRange(const std::vector<unsigned int> &atomIdxs,
                           unsigned int numAtoms, const std::string &name) {
  // Require every index to be unique and present in the molecule.
  std::unordered_set<unsigned int> seen;
  for (const auto atomIdx : atomIdxs) {
    if (atomIdx >= numAtoms) {
      invalidTemplate(name + " contains an atom index out of range");
    }
    if (!seen.insert(atomIdx).second) {
      invalidTemplate(name + " contains duplicate atom indices");
    }
  }
}

void validateLeavingGroup(const ROMol &mol,
                          const MacroMolLeavingGroup &leavingGroup,
                          const std::unordered_set<unsigned int> &mainAtoms,
                          std::set<int> &attachmentPoints,
                          std::vector<unsigned int> &atomMembership) {
  // Validate one leaving group and record its atoms and attachment point.
  if (leavingGroup.atomIdxs.empty()) {
    invalidTemplate("leaving groups cannot be empty");
  }
  validateUniqueInRange(leavingGroup.atomIdxs, mol.getNumAtoms(),
                        "leaving group");
  if (leavingGroup.attachPoint <= 0) {
    invalidTemplate("attachment-point ids must be positive");
  }
  if (!attachmentPoints.insert(leavingGroup.attachPoint).second) {
    invalidTemplate("attachment-point ids must be unique");
  }
  if (mainAtoms.find(leavingGroup.attachAtomIdx) == mainAtoms.end()) {
    invalidTemplate("attachment atom must be part of the main group");
  }

  if (std::find(leavingGroup.atomIdxs.begin(), leavingGroup.atomIdxs.end(),
                leavingGroup.leavingAtomIdx) == leavingGroup.atomIdxs.end()) {
    invalidTemplate("leaving atom must be part of its leaving group");
  }
  if (!mol.getBondBetweenAtoms(leavingGroup.attachAtomIdx,
                               leavingGroup.leavingAtomIdx)) {
    invalidTemplate("attachment atom must be bonded to the leaving atom");
  }
  for (const auto atomIdx : leavingGroup.atomIdxs) {
    if (atomMembership[atomIdx] != 0) {
      invalidTemplate("main and leaving groups must not overlap");
    }
    atomMembership[atomIdx] = 2;
  }

  // Removing the declared bond must isolate exactly the leaving-group atoms.
  RWMol disconnected(mol);
  disconnected.removeBond(leavingGroup.attachAtomIdx,
                          leavingGroup.leavingAtomIdx);
  std::vector<int> fragmentIds;
  MolOps::getMolFrags(disconnected, fragmentIds);
  const auto leavingFragment = fragmentIds[leavingGroup.leavingAtomIdx];
  if (std::any_of(leavingGroup.atomIdxs.begin(), leavingGroup.atomIdxs.end(),
                  [&fragmentIds, leavingFragment](unsigned int atomIdx) {
                    return fragmentIds[atomIdx] != leavingFragment;
                  })) {
    invalidTemplate("each leaving group must be connected");
  }
  if (static_cast<size_t>(std::count(fragmentIds.begin(), fragmentIds.end(),
                                     leavingFragment)) !=
      leavingGroup.atomIdxs.size()) {
    invalidTemplate(
        "a leaving group must have exactly its declared attachment boundary");
  }
}

}  // namespace

const SubstanceGroup &MacroMolTemplate::getMainSgroup() const {
  return getSubstanceGroups(d_mol).at(d_mainSgroupIdx);
}

MacroMolTemplateBuilder &MacroMolTemplateBuilder::setMainGroup(
    std::vector<unsigned int> atomIdxs) {
  if (d_mainGroupSet) {
    invalidTemplate("main group has already been set");
  }
  d_mainAtomIdxs = std::move(atomIdxs);
  d_mainGroupSet = true;
  return *this;
}

MacroMolTemplateBuilder &MacroMolTemplateBuilder::addLeavingGroup(
    MacroMolLeavingGroup leavingGroup) {
  d_leavingGroups.push_back(std::move(leavingGroup));
  return *this;
}

std::unique_ptr<MacroMolTemplate> MacroMolTemplateBuilder::build() const {
  if (d_name.empty()) {
    invalidTemplate("name cannot be empty");
  }
  if (d_symbol.empty()) {
    invalidTemplate("template symbol cannot be empty");
  }
  if (!d_mainGroupSet || d_mainAtomIdxs.empty()) {
    invalidTemplate("main group must be nonempty");
  }

  // Validate the complete atom partition before deriving any SGroup metadata.
  validateUniqueInRange(d_mainAtomIdxs, d_mol.getNumAtoms(), "main group");
  std::unordered_set<unsigned int> mainAtoms(d_mainAtomIdxs.begin(),
                                             d_mainAtomIdxs.end());
  std::vector<unsigned int> atomMembership(d_mol.getNumAtoms(), 0);
  for (const auto atomIdx : d_mainAtomIdxs) {
    atomMembership[atomIdx] = 1;
  }

  std::set<int> attachmentPoints;
  for (const auto &leavingGroup : d_leavingGroups) {
    validateLeavingGroup(d_mol, leavingGroup, mainAtoms, attachmentPoints,
                         atomMembership);
  }
  if (std::find(atomMembership.begin(), atomMembership.end(), 0) !=
      atomMembership.end()) {
    invalidTemplate("main and leaving groups must partition all atoms");
  }

  // Mirror the canonical template definition as SUP SGroups for MOL/SDF I/O.
  RWMol mol(d_mol);
  SubstanceGroup mainSgroup(&mol, SUP_TYPE);
  mainSgroup.setProp("CLASS",
                     std::string(monomerClassToString(d_monomerClass)));
  mainSgroup.setProp("LABEL", d_name);
  mainSgroup.setAtoms(d_mainAtomIdxs);
  for (const auto &leavingGroup : d_leavingGroups) {
    mainSgroup.addAttachPoint(leavingGroup.attachAtomIdx,
                              static_cast<int>(leavingGroup.leavingAtomIdx),
                              std::to_string(leavingGroup.attachPoint));
  }
  const auto mainSgroupIdx = addSubstanceGroup(mol, std::move(mainSgroup));

  for (const auto &leavingGroup : d_leavingGroups) {
    SubstanceGroup sgroup(&mol, SUP_TYPE);
    sgroup.setProp("CLASS", LGRP_CLASS);
    sgroup.setAtoms(leavingGroup.atomIdxs);
    addSubstanceGroup(mol, std::move(sgroup));
  }

  return std::unique_ptr<MacroMolTemplate>(new MacroMolTemplate(
      std::move(mol), d_monomerClass, d_name, d_symbol, d_originalData,
      d_mainAtomIdxs, d_leavingGroups, mainSgroupIdx));
}

void MacroMolTemplateLibrary::addTemplate(
    std::unique_ptr<MacroMolTemplate> macroMolTemplate) {
  if (!macroMolTemplate) {
    throw ValueErrorException("cannot add a null MacroMolTemplate");
  }

  const MacroMolTemplateKey nameKey{
      macroMolTemplate->getMonomerClass(), macroMolTemplate->getName()};
  if (byName.find(nameKey) != byName.end()) {
    throw ValueErrorException(
        "MacroMolTemplateLibrary already contains an entry with the same "
        "monomer class and name");
  }

  const MacroMolTemplateKey symbolKey{macroMolTemplate->getMonomerClass(),
                                      macroMolTemplate->getSymbol()};
  if (bySymbol.find(symbolKey) != bySymbol.end()) {
    throw ValueErrorException(
        "MacroMolTemplateLibrary already contains an entry with the same "
        "monomer class and symbol");
  }

  const auto *templatePtr = macroMolTemplate.get();
  byName.emplace(nameKey, std::move(macroMolTemplate));
  bySymbol.emplace(symbolKey, templatePtr);
}

const MacroMolTemplate *MacroMolTemplateLibrary::getByName(
    MonomerClass monomerClass, const std::string &name) const {
  auto it = byName.find({monomerClass, name});
  if (it != byName.end()) {
    return it->second.get();
  }
  return nullptr;
}

const MacroMolTemplate *MacroMolTemplateLibrary::getBySymbol(
    MonomerClass monomerClass, const std::string &symbol) const {
  auto it = bySymbol.find({monomerClass, symbol});
  if (it != bySymbol.end()) {
    return it->second;
  }
  return nullptr;
}

}  // namespace RDKit
