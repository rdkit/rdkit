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

namespace RDKit {

namespace {
const std::string SUP_TYPE = "SUP";
const std::string LGRP_CLASS = "LGRP";

size_t getMainGroupSize(const std::shared_ptr<MacroMolEntry> &entry) {
  if (!entry || !entry->molTemplate) {
    return 0;
  }
  const auto &templ = *entry->molTemplate;
  const auto *mainSgroup = templ.getMainSgroup();
  return mainSgroup ? mainSgroup->getAtoms().size() : 0;
}

bool isMainSgroup(const SubstanceGroup &sgroup) {
  std::string type;
  if (!sgroup.getPropIfPresent("TYPE", type) || type != SUP_TYPE) {
    return false;
  }
  std::string sgClass;
  return !sgroup.getPropIfPresent("CLASS", sgClass) || sgClass != LGRP_CLASS;
}

bool isLeavingGroup(const SubstanceGroup &sgroup) {
  std::string type;
  if (!sgroup.getPropIfPresent("TYPE", type) || type != SUP_TYPE) {
    return false;
  }
  std::string sgClass;
  return sgroup.getPropIfPresent("CLASS", sgClass) && sgClass == LGRP_CLASS;
}
}  // namespace

RWMol MacroMolTemplate::getMol() const { return d_mol; }

unsigned int MacroMolTemplate::setMainGroup(
    const std::vector<unsigned int> &atomIdxs, MonomerClass monomerClass) {
  PRECONDITION(getMainSgroup() == nullptr, "main group already set");
  SubstanceGroup sgroup(&d_mol, SUP_TYPE);
  sgroup.setProp("CLASS", monomerClassToString(monomerClass));
  sgroup.setAtoms(atomIdxs);
  return addSubstanceGroup(d_mol, sgroup);
}

unsigned int MacroMolTemplate::addLeavingGroup(
    const std::vector<unsigned int> &atomIdxs, unsigned int attachAtomIdx,
    unsigned int leavingAtomIdx, const int attachPt) {
  PRECONDITION(getMainSgroup() != nullptr, "main group must be set first");
  PRECONDITION(getMainSgroup()->includesAtom(attachAtomIdx),
               "attachment atom must be part of the main group");

  SubstanceGroup sgroup(&d_mol, SUP_TYPE);
  sgroup.setProp("CLASS", LGRP_CLASS);
  sgroup.setAtoms(atomIdxs);
  auto idx = addSubstanceGroup(d_mol, sgroup);

  getMainSgroup()->addAttachPoint(attachAtomIdx,
                                  static_cast<int>(leavingAtomIdx),
                                  std::to_string(attachPt));
  return idx;
}

const SubstanceGroup *MacroMolTemplate::getMainSgroup() const {
  for (const auto &sgroup : getSubstanceGroups(d_mol)) {
    if (isMainSgroup(sgroup)) {
      return &sgroup;
    }
  }
  return nullptr;
}

SubstanceGroup *MacroMolTemplate::getMainSgroup() {
  for (auto &sgroup : getSubstanceGroups(d_mol)) {
    if (isMainSgroup(sgroup)) {
      return &sgroup;
    }
  }
  return nullptr;
}

std::vector<const SubstanceGroup *> MacroMolTemplate::getLeavingGroups() const {
  std::vector<const SubstanceGroup *> leavingGroups;
  for (const auto &sgroup : getSubstanceGroups(d_mol)) {
    if (isLeavingGroup(sgroup)) {
      leavingGroups.push_back(&sgroup);
    }
  }
  return leavingGroups;
}

void MacroMolTemplateLibrary::addEntry(
    const std::shared_ptr<MacroMolEntry> &macroMolEntry) {
  const MacroMolTemplateKey templateKey{macroMolEntry->monomerClass,
                                        macroMolEntry->templateName};
  const auto templateIt = byTemplateName.find(templateKey);
  if (templateIt != byTemplateName.end() &&
      templateIt->second != macroMolEntry) {
    throw ValueErrorException(
        "MacroMolTemplateLibrary already contains an entry with the same "
        "monomer class and template name");
  }

  const MacroMolTemplateKey symbolKey{macroMolEntry->monomerClass,
                                      macroMolEntry->symbol};
  const auto symbolIt = bySymbol.find(symbolKey);
  if (symbolIt != bySymbol.end() && symbolIt->second != macroMolEntry) {
    throw ValueErrorException(
        "MacroMolTemplateLibrary already contains an entry with the same "
        "monomer class and symbol");
  }

  // If this exact entry is already registered, the maps and orderedEntries
  // already hold it; nothing to do. byTemplateName and bySymbol are always
  // written together below, so an entry is present in both maps or neither;
  // the conflict checks above guarantee any hit here is this same entry.
  const bool haveTemplate = templateIt != byTemplateName.end();
  const bool haveSymbol = symbolIt != bySymbol.end();
  if (haveTemplate || haveSymbol) {
    CHECK_INVARIANT(haveTemplate && haveSymbol,
                    "template/symbol maps out of sync");
    return;
  }

  byTemplateName[templateKey] = macroMolEntry;
  bySymbol[symbolKey] = macroMolEntry;

  orderedEntries.insert({getMainGroupSize(macroMolEntry), macroMolEntry});
}

std::vector<std::shared_ptr<MacroMolEntry>> MacroMolTemplateLibrary::entries()
    const {
  std::vector<std::shared_ptr<MacroMolEntry>> result;
  result.reserve(orderedEntries.size());
  for (const auto &entry : orderedEntries) {
    result.push_back(entry.second);
  }
  return result;
}

std::shared_ptr<MacroMolEntry> MacroMolTemplateLibrary::getByTemplateName(
    MonomerClass monomerClass, const std::string &templateName) const {
  auto it = byTemplateName.find({monomerClass, templateName});
  if (it != byTemplateName.end()) return it->second;
  return nullptr;
}

std::shared_ptr<MacroMolEntry> MacroMolTemplateLibrary::getBySymbol(
    MonomerClass monomerClass, const std::string &symbol) const {
  auto it = bySymbol.find({monomerClass, symbol});
  if (it != bySymbol.end()) return it->second;
  return nullptr;
}

}  // namespace RDKit
