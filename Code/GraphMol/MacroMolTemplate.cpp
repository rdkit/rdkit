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

#include <algorithm>

namespace RDKit {

namespace {
const std::string SUP_TYPE = "SUP";
const std::string LGRP_CLASS = "LGRP";

size_t getMainGroupSize(const MacroMolTemplate &templ) {
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

unsigned int MacroMolTemplate::setMainGroup(
    const std::vector<unsigned int> &atomIdxs) {
  PRECONDITION(getMainSgroup() == nullptr, "main group already set");
  SubstanceGroup sgroup(this, SUP_TYPE);
  sgroup.setProp("CLASS", monomerClassToString(d_monomerClass));
  sgroup.setAtoms(atomIdxs);
  return addSubstanceGroup(*this, sgroup);
}

unsigned int MacroMolTemplate::addLeavingGroup(
    const std::vector<unsigned int> &atomIdxs, unsigned int attachAtomIdx,
    unsigned int leavingAtomIdx, const int attachPt) {
  PRECONDITION(getMainSgroup() != nullptr, "main group must be set first");
  PRECONDITION(getMainSgroup()->includesAtom(attachAtomIdx),
               "attachment atom must be part of the main group");
  PRECONDITION(std::find(atomIdxs.begin(), atomIdxs.end(), leavingAtomIdx) !=
                   atomIdxs.end(),
               "leaving atom must be part of the leaving group");
  PRECONDITION(getBondBetweenAtoms(attachAtomIdx, leavingAtomIdx) != nullptr,
               "attachment atom must be bonded to the leaving atom");

  SubstanceGroup sgroup(this, SUP_TYPE);
  sgroup.setProp("CLASS", LGRP_CLASS);
  sgroup.setAtoms(atomIdxs);
  auto idx = addSubstanceGroup(*this, sgroup);

  getMainSgroup()->addAttachPoint(attachAtomIdx,
                                  static_cast<int>(leavingAtomIdx),
                                  std::to_string(attachPt));
  return idx;
}

const SubstanceGroup *MacroMolTemplate::getMainSgroup() const {
  for (const auto &sgroup : getSubstanceGroups(*this)) {
    if (isMainSgroup(sgroup)) {
      return &sgroup;
    }
  }
  return nullptr;
}

SubstanceGroup *MacroMolTemplate::getMainSgroup() {
  for (auto &sgroup : getSubstanceGroups(*this)) {
    if (isMainSgroup(sgroup)) {
      return &sgroup;
    }
  }
  return nullptr;
}

std::vector<const SubstanceGroup *> MacroMolTemplate::getLeavingGroups() const {
  std::vector<const SubstanceGroup *> leavingGroups;
  for (const auto &sgroup : getSubstanceGroups(*this)) {
    if (isLeavingGroup(sgroup)) {
      leavingGroups.push_back(&sgroup);
    }
  }
  return leavingGroups;
}

void MacroMolTemplateLibrary::addTemplate(
    std::unique_ptr<MacroMolTemplate> macroMolTemplate) {
  if (!macroMolTemplate) {
    throw ValueErrorException("cannot add a null MacroMolTemplate");
  }

  size_t mainGroupCount = 0;
  for (const auto &sgroup : getSubstanceGroups(*macroMolTemplate)) {
    if (isMainSgroup(sgroup)) {
      ++mainGroupCount;
    }
  }
  if (mainGroupCount != 1) {
    throw ValueErrorException(
        "MacroMolTemplate must contain exactly one main SUP SGroup");
  }

  const MacroMolTemplateKey templateKey{
      macroMolTemplate->getMonomerClass(),
      macroMolTemplate->getTemplateName()};
  if (byTemplateName.find(templateKey) != byTemplateName.end()) {
    throw ValueErrorException(
        "MacroMolTemplateLibrary already contains an entry with the same "
        "monomer class and template name");
  }

  const MacroMolTemplateKey symbolKey{macroMolTemplate->getMonomerClass(),
                                      macroMolTemplate->getSymbol()};
  if (bySymbol.find(symbolKey) != bySymbol.end()) {
    throw ValueErrorException(
        "MacroMolTemplateLibrary already contains an entry with the same "
        "monomer class and symbol");
  }

  const auto mainGroupSize = getMainGroupSize(*macroMolTemplate);
  // Keep larger main groups first while preserving insertion order among
  // templates with equal-sized main groups.
  const auto insertionPoint = std::upper_bound(
      orderedEntries.begin(), orderedEntries.end(), mainGroupSize,
      [](size_t size, const MacroMolTemplate *existingTemplate) {
        return size > getMainGroupSize(*existingTemplate);
      });

  const auto *templatePtr = macroMolTemplate.get();
  ownedTemplates.push_back(std::move(macroMolTemplate));
  byTemplateName.emplace(templateKey, templatePtr);
  bySymbol.emplace(symbolKey, templatePtr);
  orderedEntries.insert(insertionPoint, templatePtr);
}

const std::vector<const MacroMolTemplate *> &
MacroMolTemplateLibrary::entries() const {
  return orderedEntries;
}

const MacroMolTemplate *MacroMolTemplateLibrary::getByName(
    MonomerClass monomerClass, const std::string &templateName) const {
  auto it = byTemplateName.find({monomerClass, templateName});
  if (it != byTemplateName.end()) return it->second;
  return nullptr;
}

const MacroMolTemplate *MacroMolTemplateLibrary::getBySymbol(
    MonomerClass monomerClass, const std::string &symbol) const {
  auto it = bySymbol.find({monomerClass, symbol});
  if (it != bySymbol.end()) return it->second;
  return nullptr;
}

}  // namespace RDKit
