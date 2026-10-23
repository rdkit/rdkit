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
#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/FileParseException.h>

#include <algorithm>
#include <functional>
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
  const MacroMolTemplate *templ = nullptr;
  std::map<unsigned int, unsigned int> templateToResultAtom;
  std::map<int, unsigned int> attachPointToTemplateAtom;
  std::map<int, int> attachPointToLeavingAtom;
  // result atom that took the place of each removed template leaving atom
  std::map<unsigned int, unsigned int> leavingAtomReplacements;
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
  newBond->setIsAromatic(oldBond.getIsAromatic());
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
    expanded.attachPointToLeavingAtom[leavingGroup.attachPoint] =
        leavingGroup.leavingAtomIdx;
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
  expanded.templ = templ;
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

void recordLeavingAtomReplacement(BuildState &state, const Atom &atom,
                                  int attachPointId,
                                  unsigned int replacementAtomIdx) {
  const auto expanded = state.macroAtoms.find(atom.getIdx());
  if (expanded == state.macroAtoms.end()) {
    return;
  }
  const auto leavingAtom =
      expanded->second.attachPointToLeavingAtom.find(attachPointId);
  if (leavingAtom == expanded->second.attachPointToLeavingAtom.end() ||
      leavingAtom->second < 0) {
    return;
  }
  expanded->second.leavingAtomReplacements[leavingAtom->second] =
      replacementAtomIdx;
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
      recordLeavingAtomReplacement(state, *oldBond->getBeginAtom(),
                                   bondInfo.beginAttachPt, endAtomIdx);
      recordLeavingAtomReplacement(state, *oldBond->getEndAtom(),
                                   bondInfo.endAttachPt, beginAtomIdx);
    }
  }
}

// Copies oldBond's double bond stereo to resultBond, mapping its stereo atoms
// with mapIdx. E/Z labels are converted to CIS/TRANS, which refer to the
// stereo atoms and so stay valid when a neighbor is replaced by an atom with a
// different CIP rank.
void copyBondStereo(const Bond &oldBond, Bond &resultBond,
                    const std::function<unsigned int(unsigned int)> &mapIdx) {
  auto stereo = oldBond.getStereo();
  if (stereo == Bond::STEREONONE) {
    return;
  }
  const auto &stereoAtoms = oldBond.getStereoAtoms();
  if (stereoAtoms.size() == 2) {
    resultBond.setStereoAtoms(mapIdx(stereoAtoms[0]), mapIdx(stereoAtoms[1]));
    stereo = Chirality::translateEZLabelToCisTrans(stereo);
  }
  resultBond.setStereo(stereo);
}

// Restores the stereo of template atoms and bonds. This must be done after all
// bonds have been added, since a removed leaving atom is replaced by the atom
// it is bonded to through an inter-monomer bond.
void fixTemplateStereo(BuildState &state) {
  for (const auto &[macroAtomIdx, expanded] : state.macroAtoms) {
    const auto mapIdx = [&expanded](unsigned int templAtomIdx) {
      const auto resultAtom = expanded.templateToResultAtom.find(templAtomIdx);
      if (resultAtom != expanded.templateToResultAtom.end()) {
        return resultAtom->second;
      }
      return expanded.leavingAtomReplacements.at(templAtomIdx);
    };
    const auto &templMol = expanded.templ->getMol();

    // Chiral tags are relative to the order of an atom's bonds. Replacing a
    // leaving atom's bond with an inter-monomer bond, which is added last,
    // can change that order, so restore the template's parity where needed.
    for (const auto *templAtom : templMol.atoms()) {
      const auto resultAtomIdx =
          expanded.templateToResultAtom.find(templAtom->getIdx());
      if (resultAtomIdx == expanded.templateToResultAtom.end()) {
        continue;  // a removed leaving atom
      }
      auto *resultAtom = state.result.getAtomWithIdx(resultAtomIdx->second);
      Chirality::restoreTetrahedralParity(*templAtom, *resultAtom, mapIdx);
    }

    for (const auto *templBond : templMol.bonds()) {
      const auto beginAtom =
          expanded.templateToResultAtom.find(templBond->getBeginAtomIdx());
      const auto endAtom =
          expanded.templateToResultAtom.find(templBond->getEndAtomIdx());
      if (beginAtom == expanded.templateToResultAtom.end() ||
          endAtom == expanded.templateToResultAtom.end()) {
        continue;  // a removed leaving group bond
      }
      auto *resultBond =
          state.result.getBondBetweenAtoms(beginAtom->second, endAtom->second);
      copyBondStereo(*templBond, *resultBond, mapIdx);
    }
  }
}

// Copies the double bond stereo of bonds between ordinary atoms. Stereo whose
// stereo atoms include a macro atom is not copied.
void fixMacroMolBondStereo(BuildState &state, const MacroMol &macroMol) {
  const auto isMapped = [&state](int atomIdx) {
    return state.atomToResultAtom.find(atomIdx) !=
           state.atomToResultAtom.end();
  };
  const auto mapIdx = [&state](unsigned int atomIdx) {
    return state.atomToResultAtom.at(atomIdx);
  };
  for (const auto *oldBond : macroMol.bonds()) {
    if (oldBond->getMacroBondInfo() ||
        oldBond->getStereo() == Bond::STEREONONE) {
      continue;
    }
    const auto &stereoAtoms = oldBond->getStereoAtoms();
    if (!std::all_of(stereoAtoms.begin(), stereoAtoms.end(), isMapped)) {
      continue;
    }
    auto *resultBond = state.result.getBondBetweenAtoms(
        mapIdx(oldBond->getBeginAtomIdx()), mapIdx(oldBond->getEndAtomIdx()));
    copyBondStereo(*oldBond, *resultBond, mapIdx);
  }
}

// Double bond stereo is fully described by the bond stereo and stereo atoms.
// The copied bond directions are incomplete where leaving atoms were replaced,
// and the legacy stereo perception would use them to override the stereo, so
// regenerate them.
void regenerateDoubleBondNeighborDirections(RWMol &mol) {
  for (auto *bond : mol.bonds()) {
    if (bond->getBondDir() == Bond::ENDUPRIGHT ||
        bond->getBondDir() == Bond::ENDDOWNRIGHT) {
      bond->setBondDir(Bond::NONE);
    }
  }
  MolOps::setDoubleBondNeighborDirections(mol);
}

}  // namespace

std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol &macroMol, const MacroMolTemplateLibrary &templates) {
  auto result = std::make_unique<RWMol>();
  BuildState state{*result, {}, {}};

  copyAtoms(state, macroMol, templates);
  copyMacroMolBonds(state, macroMol);
  fixTemplateStereo(state);
  fixMacroMolBondStereo(state, macroMol);

  result->updatePropertyCache(false);
  regenerateDoubleBondNeighborDirections(*result);
  return result;
}

}  // namespace RDKit
