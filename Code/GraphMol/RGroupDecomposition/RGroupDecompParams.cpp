//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RGroupDecomp.h"
#include "RGroupUtils.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FMCS/FMCS.h>

namespace RDKit
{
unsigned int RGroupDecompositionParameters::autoGetLabels(const RWMol &core) {
  unsigned int autoLabels = 0;
  if (!onlyMatchAtRGroups) {
    autoLabels = AtomIndexLabels;
  }
  bool hasMDLRGroup = false;
  bool hasAtomMapNum = false;
  bool hasIsotopes = false;
  bool hasDummies = false;
  for (auto atm : core.atoms()) {
    if (atm->getIsotope()) {
      hasIsotopes = true;
    }
    if (atm->getAtomMapNum()) {
      hasAtomMapNum = true;
    }
    if (atm->hasProp(common_properties::_MolFileRLabel)) {
      hasMDLRGroup = true;
    }
    if (atm->getAtomicNum() == 0) {
      hasDummies = true;
    }
  }

  if (hasMDLRGroup) {
    return autoLabels | MDLRGroupLabels;
  } else if (hasAtomMapNum) {
    return autoLabels | AtomMapLabels;
  } else if (hasIsotopes) {
    return autoLabels | IsotopeLabels;
  } else if (hasDummies) {
    return autoLabels | DummyAtomLabels;
  }

  return autoLabels;
}

bool rgdAtomCompare(const MCSAtomCompareParameters &p, const ROMol &mol1,
                    unsigned int atom1, const ROMol &mol2, unsigned int atom2,
                    void *userData) {
  if (!MCSAtomCompareElements(p, mol1, atom1, mol2, atom2, nullptr)) {
    return false;
  }
  unsigned int autoLabels = *reinterpret_cast<unsigned int *>(userData);
  bool atom1HasLabel = false;
  bool atom2HasLabel = false;
  const auto a1 = mol1.getAtomWithIdx(atom1);
  const auto a2 = mol2.getAtomWithIdx(atom2);
  if (autoLabels & MDLRGroupLabels) {
    atom1HasLabel |= a1->hasProp(common_properties::_MolFileRLabel);
    atom2HasLabel |= a2->hasProp(common_properties::_MolFileRLabel);
  }
  if (autoLabels & IsotopeLabels) {
    atom2HasLabel |= (a1->getIsotope() > 0);
    atom2HasLabel |= (a2->getIsotope() > 0);
  }
  if (autoLabels & AtomMapLabels) {
    atom1HasLabel |= (a1->getAtomMapNum() > 0);
    atom2HasLabel |= (a2->getAtomMapNum() > 0);
  }
  if (autoLabels & DummyAtomLabels) {
    atom1HasLabel |= (a1->getAtomicNum() == 0);
    atom2HasLabel |= (a2->getAtomicNum() == 0);
  }
  atom1HasLabel |= a1->hasProp(RLABEL);
  atom2HasLabel |= a2->hasProp(RLABEL);
  return !(atom1HasLabel ^ atom2HasLabel);
}

bool RGroupDecompositionParameters::prepareCore(RWMol &core,
                                                const RWMol *alignCore) {
  const bool relabel = labels & RelabelDuplicateLabels;
  unsigned int autoLabels = labels;
  if (labels == AutoDetect) {
    autoLabels = autoGetLabels(core);
    if (!autoLabels) {
      BOOST_LOG(rdWarningLog) << "RGroupDecomposition auto detect found no "
                                 "rgroups and onlyMatAtRgroups is set to true"
                              << std::endl;
      return false;
    }
  }

  int maxLabel = 1;
  if (alignCore && (alignment & MCS)) {
    std::vector<ROMOL_SPTR> mols;
    mols.push_back(ROMOL_SPTR(new ROMol(core)));
    mols.push_back(ROMOL_SPTR(new ROMol(*alignCore)));
    MCSParameters mcsParams;
    if (labels != AutoDetect) {
      mcsParams.AtomTyper = rgdAtomCompare;
      mcsParams.CompareFunctionsUserData = &autoLabels;
    }
    MCSResult res = findMCS(mols, &mcsParams);
    if (res.isCompleted()) {
      RWMol *m = SmartsToMol(res.SmartsString);
      if (m) {
        MatchVectType match1;
        MatchVectType match2;

        bool target_matched1 = SubstructMatch(core, *m, match1);
        bool target_matched2 = SubstructMatch(*alignCore, *m, match2);
        CHECK_INVARIANT(match1.size() == match2.size(),
                        "Matches should be the same size in prepareCore");

        if (target_matched1 && target_matched2) {
          for (size_t i = 0; i < match1.size(); ++i) {
            int queryAtomIdx1 = match1[i].first;
            int coreAtomIdx = match1[i].second;
            int queryAtomIdx2 = match2[i].first;
            int alignCoreAtomIdx = match2[i].second;
            CHECK_INVARIANT(queryAtomIdx1 == queryAtomIdx2,
                            "query atoms aren't the same");
            Atom *coreAtm = core.getAtomWithIdx(coreAtomIdx);
            const Atom *alignCoreAtm =
                alignCore->getAtomWithIdx(alignCoreAtomIdx);

            // clear up input rlabels
            coreAtm->setAtomMapNum(0);
            if (coreAtm->hasProp(common_properties::_MolFileRLabel)) {
              coreAtm->clearProp(common_properties::_MolFileRLabel);
              coreAtm->setIsotope(0);
            }
            if (alignCoreAtm->hasProp(RLABEL)) {
              int rlabel = alignCoreAtm->getProp<int>(RLABEL);
              maxLabel = (std::max)(maxLabel, rlabel + 1);
              coreAtm->setProp(RLABEL, rlabel);
            }
          }
        }
        delete m;
      }
    }
  }
  std::set<int> foundLabels;

  int nextOffset = 0;
  std::map<int, int> atomToLabel;

  for (auto atom : core.atoms()) {
    bool found = false;

    if (atom->hasProp(RLABEL)) {
      if (setLabel(atom, atom->getProp<int>(RLABEL), foundLabels, maxLabel,
                   relabel, Labelling::INTERNAL_LABELS)) {
        found = true;
      }
    }

    if (!found && (autoLabels & MDLRGroupLabels)) {
      unsigned int rgroup;
      if (atom->getPropIfPresent<unsigned int>(
              common_properties::_MolFileRLabel, rgroup)) {
        if (setLabel(atom, rdcast<int>(rgroup), foundLabels, maxLabel, relabel,
                     Labelling::RGROUP_LABELS)) {
          found = true;
        }
      }
    }

    if (!found && (autoLabels & IsotopeLabels) && atom->getIsotope() > 0) {
      if (setLabel(atom, rdcast<int>(atom->getIsotope()), foundLabels, maxLabel,
                   relabel, Labelling::ISOTOPE_LABELS)) {
        found = true;
      }
    }

    if (!found && (autoLabels & AtomMapLabels) && atom->getAtomMapNum() > 0) {
      if (setLabel(atom, rdcast<int>(atom->getAtomMapNum()), foundLabels,
                   maxLabel, relabel, Labelling::ATOMMAP_LABELS)) {
        found = true;
      }
    }

    if (!found && (autoLabels & DummyAtomLabels) && atom->getAtomicNum() == 0) {
      const bool forceRelabellingWithDummies = true;
      int defaultDummyStartLabel = maxLabel;
      if (setLabel(atom, defaultDummyStartLabel, foundLabels, maxLabel,
                   forceRelabellingWithDummies, Labelling::DUMMY_LABELS)) {
        found = true;
      }
    }

    // Unless there is an MCS match from above, we need to give different
    //  RLABELS to each core so keep track of which labels
    //  we have used (note that these are negative since they are
    //  potential rgroups and haven't been assigned yet)
    if (!found && (autoLabels & AtomIndexLabels)) {
      if (setLabel(atom, indexOffset - atom->getIdx(), foundLabels, maxLabel,
                   relabel, Labelling::INDEX_LABELS)) {
        nextOffset++;
      }
      found = true;
    }

    clearInputLabels(atom);

    int rlabel;
    if (atom->getPropIfPresent(RLABEL, rlabel)) {
      atomToLabel[atom->getIdx()] = rlabel;
    }
  }
  indexOffset -= nextOffset;

  MolOps::AdjustQueryParameters adjustParams;
  adjustParams.makeDummiesQueries = true;
  adjustParams.adjustDegree = false;
  adjustQueryProperties(core, &adjustParams);
  for (auto &it : atomToLabel) {
    core.getAtomWithIdx(it.first)->setProp(RLABEL, it.second);
  }
  return true;
}

}
