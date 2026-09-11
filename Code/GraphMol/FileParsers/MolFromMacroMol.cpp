//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/FileParsers/MolFromMacroMol.h>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/FileParseException.h>

#include <map>
#include <memory>
#include <set>
#include <string>

namespace RDKit {
namespace {

struct MacroAtomDetails {
  MonomerClass monomerClass = MonomerClass::Other;
  std::string label;
};

struct ExpandedMacroAtom {
  std::map<unsigned int, unsigned int> templateToResultAtom;
  std::map<int, unsigned int> attachPointToTemplateAtom;
  std::set<int> usedAttachPoints;
};

struct BuildState {
  RWMol &result;
  std::map<unsigned int, unsigned int> atomToResultAtom;
  std::map<unsigned int, ExpandedMacroAtom> macroAtoms;
};

bool getMacroAtomDetails(const Atom &atom, MacroAtomDetails &info) {
  const auto *macroInfo = atom.getMacroAtomInfo();
  if (!macroInfo) {
    return false;
  }
  info.monomerClass = macroInfo->getMonomerClass();
  info.label = macroInfo->getSymbol();
  return !info.label.empty();
}

const MacroMolTemplate *findTemplate(
    const MacroMolTemplateLibrary &templates, const MacroAtomDetails &info) {
  const auto *templ = templates.getBySymbol(info.monomerClass, info.label);
  if (!templ) {
    templ = templates.getByName(info.monomerClass, info.label);
  }
  return templ;
}

void addBondCopy(RWMol &result, const Bond &oldBond,
                 unsigned int beginAtomIdx, unsigned int endAtomIdx,
                 Bond::BondType bondType) {
  auto *newBond = new Bond(bondType);
  newBond->setBeginAtomIdx(beginAtomIdx);
  newBond->setEndAtomIdx(endAtomIdx);
  newBond->setBondDir(oldBond.getBondDir());
  newBond->updateProps(oldBond, false);
  result.addBond(newBond, true);
}

void addBondCopy(RWMol &result, const Bond &oldBond,
                 unsigned int beginAtomIdx, unsigned int endAtomIdx) {
  addBondCopy(result, oldBond, beginAtomIdx, endAtomIdx,
              oldBond.getBondType());
}

void copyRegularAtom(BuildState &state, const Atom &atom) {
  auto *newAtom = new Atom(atom);
  state.atomToResultAtom[atom.getIdx()] =
      state.result.addAtom(newAtom, true, true);
}

std::set<int> findUsedAttachPoints(const MacroMol &macroMol,
                                   const Atom &macroAtom) {
  std::set<int> usedAttachPoints;
  for (const auto *macroBond : macroMol.bonds()) {
    const auto *macroBondInfo = macroBond->getMacroBondInfo();
    if (!macroBondInfo) {
      continue;
    }

    for (const auto &bondInfo : macroBondInfo->getBonds()) {
      int attachPoint = -1;
      if (macroBond->getBeginAtomIdx() == macroAtom.getIdx()) {
        attachPoint = bondInfo.beginAttachPt;
      } else if (macroBond->getEndAtomIdx() == macroAtom.getIdx()) {
        attachPoint = bondInfo.endAttachPt;
      }
      if (attachPoint >= 0) {
        usedAttachPoints.insert(attachPoint);
      }
    }
  }
  return usedAttachPoints;
}

std::set<unsigned int> findLeavingAtomsToSkip(
    const MacroMolTemplate &templ, const std::set<int> &usedAttachPoints) {
  std::set<unsigned int> atomsToSkip;
  for (const auto &leavingGroup : templ.getLeavingGroups()) {
    if (usedAttachPoints.find(leavingGroup.attachPoint) !=
        usedAttachPoints.end()) {
      atomsToSkip.insert(leavingGroup.atomIdxs.begin(),
                         leavingGroup.atomIdxs.end());
    }
  }
  return atomsToSkip;
}

void rememberTemplateAttachPoints(ExpandedMacroAtom &expanded,
                                  const MacroMolTemplate &templ) {
  for (const auto &leavingGroup : templ.getLeavingGroups()) {
    expanded.attachPointToTemplateAtom[leavingGroup.attachPoint] =
        leavingGroup.attachAtomIdx;
  }
}

void copyTemplateAtoms(BuildState &state, ExpandedMacroAtom &expanded,
                       const MacroMolTemplate &templ,
                       const std::set<unsigned int> &atomsToSkip) {
  for (const auto *oldAtom : templ.getMol().atoms()) {
    if (atomsToSkip.find(oldAtom->getIdx()) != atomsToSkip.end()) {
      continue;
    }

    auto *newAtom = new Atom(*oldAtom);
    expanded.templateToResultAtom[oldAtom->getIdx()] =
        state.result.addAtom(newAtom, true, true);
  }
}

void copyTemplateBonds(BuildState &state, const ExpandedMacroAtom &expanded,
                       const MacroMolTemplate &templ) {
  for (const auto *oldBond : templ.getMol().bonds()) {
    const auto beginAtom =
        expanded.templateToResultAtom.find(oldBond->getBeginAtomIdx());
    const auto endAtom =
        expanded.templateToResultAtom.find(oldBond->getEndAtomIdx());
    if (beginAtom == expanded.templateToResultAtom.end() ||
        endAtom == expanded.templateToResultAtom.end()) {
      continue;
    }
    addBondCopy(state.result, *oldBond, beginAtom->second, endAtom->second);
  }
}

void copyMainSgroup(BuildState &state, const ExpandedMacroAtom &expanded,
                    const MacroMolTemplate &templ) {
  SubstanceGroup mainSgroup(&state.result, "SUP");
  mainSgroup.setProp(
      "CLASS", std::string(monomerClassToString(templ.getMonomerClass())));
  mainSgroup.setProp("LABEL", templ.getName());

  std::vector<unsigned int> resultMainAtoms;
  resultMainAtoms.reserve(templ.getMainAtomIdxs().size());
  for (const auto templateAtomIdx : templ.getMainAtomIdxs()) {
    resultMainAtoms.push_back(
        expanded.templateToResultAtom.at(templateAtomIdx));
  }
  mainSgroup.setAtoms(std::move(resultMainAtoms));

  for (const auto &leavingGroup : templ.getLeavingGroups()) {
    const auto resultAttachAtom =
        expanded.templateToResultAtom.at(leavingGroup.attachAtomIdx);
    int resultLeavingAtom = -1;
    const auto leavingAtom =
        expanded.templateToResultAtom.find(leavingGroup.leavingAtomIdx);
    if (leavingAtom != expanded.templateToResultAtom.end()) {
      resultLeavingAtom = static_cast<int>(leavingAtom->second);
    }
    mainSgroup.addAttachPoint(resultAttachAtom, resultLeavingAtom,
                              std::to_string(leavingGroup.attachPoint));
  }
  addSubstanceGroup(state.result, std::move(mainSgroup));
}

void copyMacroAtom(BuildState &state, const MacroMol &macroMol,
                   const MacroMolTemplateLibrary &templates, const Atom &atom,
                   const MacroAtomDetails &info) {
  const auto *templ = findTemplate(templates, info);
  if (!templ) {
    throw FileParseException("No template found for macro atom " +
                             std::to_string(atom.getIdx()));
  }

  ExpandedMacroAtom expanded;
  expanded.usedAttachPoints = findUsedAttachPoints(macroMol, atom);
  const auto atomsToSkip =
      findLeavingAtomsToSkip(*templ, expanded.usedAttachPoints);

  rememberTemplateAttachPoints(expanded, *templ);
  copyTemplateAtoms(state, expanded, *templ, atomsToSkip);
  copyTemplateBonds(state, expanded, *templ);
  copyMainSgroup(state, expanded, *templ);

  state.macroAtoms[atom.getIdx()] = expanded;
}

void copyAtoms(BuildState &state, const MacroMol &macroMol,
               const MacroMolTemplateLibrary &templates) {
  for (const auto *atom : macroMol.atoms()) {
    MacroAtomDetails info;
    if (getMacroAtomDetails(*atom, info)) {
      copyMacroAtom(state, macroMol, templates, *atom, info);
    } else {
      copyRegularAtom(state, *atom);
    }
  }
}

unsigned int getResultAtomForBond(const BuildState &state, const Atom &atom,
                                  int attachPointId) {
  MacroAtomDetails info;
  if (!getMacroAtomDetails(atom, info)) {
    return state.atomToResultAtom.at(atom.getIdx());
  }

  const auto &expanded = state.macroAtoms.at(atom.getIdx());
  const auto attachPoint =
      expanded.attachPointToTemplateAtom.find(attachPointId);
  if (attachPoint == expanded.attachPointToTemplateAtom.end()) {
    throw FileParseException("Macro atom bond is missing an attachment point");
  }

  return expanded.templateToResultAtom.at(attachPoint->second);
}

void copyMacroMolBonds(BuildState &state, const MacroMol &macroMol) {
  for (const auto *oldBond : macroMol.bonds()) {
    const auto *macroBondInfo = oldBond->getMacroBondInfo();
    if (!macroBondInfo) {
      const auto beginAtomIdx =
          getResultAtomForBond(state, *oldBond->getBeginAtom(), -1);
      const auto endAtomIdx =
          getResultAtomForBond(state, *oldBond->getEndAtom(), -1);
      addBondCopy(state.result, *oldBond, beginAtomIdx, endAtomIdx);
      continue;
    }

    for (const auto &bondInfo : macroBondInfo->getBonds()) {
      const auto beginAtomIdx = getResultAtomForBond(
          state, *oldBond->getBeginAtom(), bondInfo.beginAttachPt);
      const auto endAtomIdx = getResultAtomForBond(
          state, *oldBond->getEndAtom(), bondInfo.endAttachPt);
      addBondCopy(state.result, *oldBond, beginAtomIdx, endAtomIdx,
                  static_cast<Bond::BondType>(bondInfo.bondType));
    }
  }
}

}  // namespace

std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol &macroMol, const MacroMolTemplateLibrary &templates) {
  auto result = std::make_unique<RWMol>();
  BuildState state{*result, {}, {}};

  copyAtoms(state, macroMol, templates);
  copyMacroMolBonds(state, macroMol);

  result->updatePropertyCache(false);
  return result;
}

}  // namespace RDKit
