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
    const std::vector<unsigned int> &atomIdxs, const std::string &className) {
  PRECONDITION(getMainSgroup() == nullptr, "main group already set");
  SubstanceGroup sgroup(this, SUP_TYPE);
  sgroup.setProp("CLASS", className);
  sgroup.setAtoms(atomIdxs);
  return addSubstanceGroup(*this, sgroup);
}

unsigned int MacroMolTemplate::addLeavingGroup(
    const std::vector<unsigned int> &atomIdxs, unsigned int attachAtomIdx,
    unsigned int leavingAtomIdx, const int attachPt) {
  PRECONDITION(getMainSgroup() != nullptr, "main group must be set first");

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

void MacroMolTemplateLibrary::addEntry(
    const std::shared_ptr<MacroMolEntry> &macroMolEntry) {
  byTemplateName[{macroMolEntry->monomerClass, macroMolEntry->templateName}] =
      macroMolEntry;
  bySymbol[{macroMolEntry->monomerClass, macroMolEntry->symbol}] =
      macroMolEntry;
}

const std::shared_ptr<MacroMolEntry> &
MacroMolTemplateLibrary::getByTemplateName(
    const std::string &monomerClass, const std::string &templateName) const {
  static std::shared_ptr<MacroMolEntry> empty;
  auto it = byTemplateName.find({monomerClass, templateName});
  if (it != byTemplateName.end()) return it->second;
  return empty;
}

const std::shared_ptr<MacroMolEntry> &
MacroMolTemplateLibrary::getBySymbol(const std::string &monomerClass,
                                     const std::string &symbol) const {
  static std::shared_ptr<MacroMolEntry> empty;
  auto it = bySymbol.find({monomerClass, symbol});
  if (it != bySymbol.end()) return it->second;
  return empty;
}

}  // namespace RDKit
