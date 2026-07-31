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

#include <RDGeneral/Exceptions.h>

#include <algorithm>
#include <queue>
#include <set>
#include <unordered_set>

namespace RDKit {

namespace {
const std::string SUP_TYPE = "SUP";
const std::string LGRP_CLASS = "LGRP";

[[noreturn]] void invalidTemplate(const std::string &message) {
  throw ValueErrorException("invalid MacroMolTemplate: " + message);
}

void validateUniqueInRange(const std::vector<unsigned int> &atomIdxs,
                           unsigned int numAtoms, const std::string &name) {
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

bool isConnected(const ROMol &mol,
                 const std::vector<unsigned int> &atomIdxs) {
  std::unordered_set<unsigned int> groupAtoms(atomIdxs.begin(), atomIdxs.end());
  std::unordered_set<unsigned int> visited;
  std::queue<unsigned int> pending;
  pending.push(atomIdxs.front());
  visited.insert(atomIdxs.front());

  while (!pending.empty()) {
    const auto atomIdx = pending.front();
    pending.pop();
    for (const auto *bond : mol.atomBonds(mol.getAtomWithIdx(atomIdx))) {
      const auto otherAtomIdx = bond->getOtherAtomIdx(atomIdx);
      if (groupAtoms.find(otherAtomIdx) != groupAtoms.end() &&
          visited.insert(otherAtomIdx).second) {
        pending.push(otherAtomIdx);
      }
    }
  }
  return visited.size() == atomIdxs.size();
}

void validateLeavingGroup(const ROMol &mol,
                          const MacroMolLeavingGroup &leavingGroup,
                          const std::unordered_set<unsigned int> &mainAtoms,
                          std::set<int> &attachmentPoints,
                          std::vector<unsigned int> &atomMembership) {
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

  const std::unordered_set<unsigned int> leavingAtoms(
      leavingGroup.atomIdxs.begin(), leavingGroup.atomIdxs.end());
  if (leavingAtoms.find(leavingGroup.leavingAtomIdx) == leavingAtoms.end()) {
    invalidTemplate("leaving atom must be part of its leaving group");
  }
  if (!mol.getBondBetweenAtoms(leavingGroup.attachAtomIdx,
                               leavingGroup.leavingAtomIdx)) {
    invalidTemplate("attachment atom must be bonded to the leaving atom");
  }
  if (!isConnected(mol, leavingGroup.atomIdxs)) {
    invalidTemplate("each leaving group must be connected");
  }

  unsigned int boundaryBondCount = 0;
  bool foundAttachmentBond = false;
  for (const auto atomIdx : leavingGroup.atomIdxs) {
    if (atomMembership[atomIdx] != 0) {
      invalidTemplate("main and leaving groups must not overlap");
    }
    atomMembership[atomIdx] = 2;
    for (const auto *bond : mol.atomBonds(mol.getAtomWithIdx(atomIdx))) {
      const auto otherAtomIdx = bond->getOtherAtomIdx(atomIdx);
      if (leavingAtoms.find(otherAtomIdx) == leavingAtoms.end()) {
        ++boundaryBondCount;
        if (atomIdx == leavingGroup.leavingAtomIdx &&
            otherAtomIdx == leavingGroup.attachAtomIdx) {
          foundAttachmentBond = true;
        }
      }
    }
  }
  if (boundaryBondCount != 1 || !foundAttachmentBond) {
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

std::unique_ptr<MacroMolTemplate> MacroMolTemplateBuilder::build() && {
  if (d_name.empty()) {
    invalidTemplate("name cannot be empty");
  }
  if (d_symbol.empty()) {
    invalidTemplate("template symbol cannot be empty");
  }
  if (!d_mainGroupSet || d_mainAtomIdxs.empty()) {
    invalidTemplate("main group must be nonempty");
  }

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

  SubstanceGroup mainSgroup(&d_mol, SUP_TYPE);
  mainSgroup.setProp("CLASS",
                     std::string(monomerClassToString(d_monomerClass)));
  mainSgroup.setProp("LABEL", d_name);
  mainSgroup.setAtoms(d_mainAtomIdxs);
  for (const auto &leavingGroup : d_leavingGroups) {
    mainSgroup.addAttachPoint(leavingGroup.attachAtomIdx,
                              static_cast<int>(leavingGroup.leavingAtomIdx),
                              std::to_string(leavingGroup.attachPoint));
  }
  const auto mainSgroupIdx = addSubstanceGroup(d_mol, std::move(mainSgroup));

  for (const auto &leavingGroup : d_leavingGroups) {
    SubstanceGroup sgroup(&d_mol, SUP_TYPE);
    sgroup.setProp("CLASS", LGRP_CLASS);
    sgroup.setAtoms(leavingGroup.atomIdxs);
    addSubstanceGroup(d_mol, std::move(sgroup));
  }

  return std::unique_ptr<MacroMolTemplate>(new MacroMolTemplate(
      std::move(d_mol), d_monomerClass, std::move(d_name),
      std::move(d_symbol), std::move(d_originalData),
      std::move(d_mainAtomIdxs), std::move(d_leavingGroups), mainSgroupIdx));
}

namespace {
size_t getMainGroupSize(const MacroMolTemplate &templ) {
  return templ.getMainAtomIdxs().size();
}
}  // namespace

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

  const auto mainGroupSize = getMainGroupSize(*macroMolTemplate);
  const auto insertionPoint = std::upper_bound(
      orderedEntries.begin(), orderedEntries.end(), mainGroupSize,
      [](size_t size, const MacroMolTemplate *existingTemplate) {
        return size > getMainGroupSize(*existingTemplate);
      });

  const auto *templatePtr = macroMolTemplate.get();
  ownedTemplates.push_back(std::move(macroMolTemplate));
  byName.emplace(nameKey, templatePtr);
  bySymbol.emplace(symbolKey, templatePtr);
  orderedEntries.insert(insertionPoint, templatePtr);
}

const std::vector<const MacroMolTemplate *> &
MacroMolTemplateLibrary::entries() const {
  return orderedEntries;
}

const MacroMolTemplate *MacroMolTemplateLibrary::getByName(
    MonomerClass monomerClass, const std::string &name) const {
  auto it = byName.find({monomerClass, name});
  if (it != byName.end()) return it->second;
  return nullptr;
}

const MacroMolTemplate *MacroMolTemplateLibrary::getBySymbol(
    MonomerClass monomerClass, const std::string &symbol) const {
  auto it = bySymbol.find({monomerClass, symbol});
  if (it != bySymbol.end()) return it->second;
  return nullptr;
}

}  // namespace RDKit
