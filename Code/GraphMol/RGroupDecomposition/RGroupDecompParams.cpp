//
//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RGroupDecompParams.h"
#include "RGroupUtils.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/QueryBond.h>
#include <set>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

namespace {

bool hasLabel(const Atom *atom, unsigned int autoLabels) {
  bool atomHasLabel = false;
  if (autoLabels & MDLRGroupLabels) {
    atomHasLabel |= atom->hasProp(common_properties::_MolFileRLabel);
  }
  if (autoLabels & IsotopeLabels) {
    atomHasLabel |= (atom->getIsotope() > 0);
  }
  if (autoLabels & AtomMapLabels) {
    atomHasLabel |= (atom->getAtomMapNum() > 0);
  }
  if (autoLabels & DummyAtomLabels) {
    atomHasLabel |= (atom->getAtomicNum() == 0 && atom->getDegree() == 1);
  }
  // don't match negative rgroups as these are used by AtomIndexRLabels which
  // are set for the template core before the MCS and after the MCS for the
  // other core
  if (atom->hasProp(RLABEL)) {
    auto label = atom->getProp<int>(RLABEL);
    atomHasLabel |= label > 0;
  }
  return atomHasLabel;
}

/* When comparing atoms for the MCS overlay, we don't want to overlay a core
 * atom that has a labelled rgroup on top of one that hasn't a labelled rgroup.
 * The labelled rgroup may be on a terminal dummy atom bonded to the core atom.
 *
 * This function checks to see if there is a terminal dummy labelled rgroup
 * attached to the atom.
 */
bool hasAttachedLabels(const ROMol &mol, const Atom *atom,
                       unsigned int autoLabels) {
  RWMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
  while (nbrIdx != endNbrs) {
    const auto neighborAtom = mol.getAtomWithIdx(*nbrIdx);
    if (neighborAtom->getAtomicNum() == 0 && neighborAtom->getDegree() == 1 &&
        hasLabel(neighborAtom, autoLabels)) {
      return true;
    }
    ++nbrIdx;
  }
  return false;
}

}  // namespace

void updateRGroupDecompositionParametersFromJSON(
    RGroupDecompositionParameters &params, const std::string &details_json) {
  updateRGroupDecompositionParametersFromJSON(params, details_json.c_str());
}

void updateRGroupDecompositionParametersFromJSON(
    RGroupDecompositionParameters &params, const char *details_json) {
  static const std::map<std::string, RGroupLabels> rGroupLabelsMap{
      RGROUPLABELS_ENUM_ITEMS};
  static const std::map<std::string, RGroupMatching> rGroupMatchingMap{
      RGROUPMATCHING_ENUM_ITEMS};
  static const std::map<std::string, RGroupLabelling> rGroupLabellingMap{
      RGROUPLABELLING_ENUM_ITEMS};
  static const std::map<std::string, RGroupCoreAlignment>
      rGroupCoreAlignmentMap{RGROUPCOREALIGNMENT_ENUM_ITEMS};
  static const std::map<std::string, RGroupScore> rGroupScoreMap{
      RGROUPSCORE_ENUM_ITEMS};
  if (details_json && strlen(details_json)) {
    std::istringstream ss;
    boost::property_tree::ptree pt;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);

    std::string labels;
    labels = pt.get<std::string>("labels", labels);
    auto rGroupLabelsMapIt = rGroupLabelsMap.find(labels);
    if (rGroupLabelsMapIt != rGroupLabelsMap.end()) {
      params.labels = rGroupLabelsMapIt->second;
    }

    std::string matchingStrategy;
    matchingStrategy =
        pt.get<std::string>("matchingStrategy", matchingStrategy);
    auto rGroupMatchingMapIt = rGroupMatchingMap.find(matchingStrategy);
    if (rGroupMatchingMapIt != rGroupMatchingMap.end()) {
      params.matchingStrategy = rGroupMatchingMapIt->second;
    }

    std::string scoreMethod;
    scoreMethod = pt.get<std::string>("scoreMethod", scoreMethod);
    auto rGroupScoreMapIt = rGroupScoreMap.find(scoreMethod);
    if (rGroupScoreMapIt != rGroupScoreMap.end()) {
      params.scoreMethod = rGroupScoreMapIt->second;
    }

    std::string rgroupLabelling;
    rgroupLabelling = pt.get<std::string>("rgroupLabelling", rgroupLabelling);
    auto rGroupLabellingMapIt = rGroupLabellingMap.find(rgroupLabelling);
    if (rGroupLabellingMapIt != rGroupLabellingMap.end()) {
      params.rgroupLabelling = rGroupLabellingMapIt->second;
    }

    std::string alignment;
    alignment = pt.get<std::string>("alignment", alignment);
    auto rGroupCoreAlignmentMapIt = rGroupCoreAlignmentMap.find(alignment);
    if (rGroupCoreAlignmentMapIt != rGroupCoreAlignmentMap.end()) {
      params.alignment = rGroupCoreAlignmentMapIt->second;
    }

    params.chunkSize = pt.get<unsigned int>("chunkSize", params.chunkSize);
    params.onlyMatchAtRGroups =
        pt.get<bool>("onlyMatchAtRGroups", params.onlyMatchAtRGroups);
    params.removeAllHydrogenRGroups = pt.get<bool>(
        "removeAllHydrogenRGroups", params.removeAllHydrogenRGroups);
    params.removeAllHydrogenRGroupsAndLabels =
        pt.get<bool>("removeAllHydrogenRGroupsAndLabels",
                     params.removeAllHydrogenRGroupsAndLabels);
    params.removeHydrogensPostMatch = pt.get<bool>(
        "removeHydrogensPostMatch", params.removeHydrogensPostMatch);
    params.allowNonTerminalRGroups =
        pt.get<bool>("allowNonTerminalRGroups", params.allowNonTerminalRGroups);
    params.allowMultipleRGroupsOnUnlabelled =
        pt.get<bool>("allowMultipleRGroupsOnUnlabelled",
                     params.allowMultipleRGroupsOnUnlabelled);
    params.doTautomers = pt.get<bool>("doTautomers", params.doTautomers);
    params.doEnumeration = pt.get<bool>("doEnumeration", params.doEnumeration);
    params.includeTargetMolInResults = pt.get<bool>(
        "includeTargetMolInResults", params.includeTargetMolInResults);
    params.timeout = pt.get<double>("timeout", params.timeout);
  }
}

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
    if (atm->getAtomicNum() == 0 && atm->getDegree() == 1) {
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
  const auto a1 = mol1.getAtomWithIdx(atom1);
  const auto a2 = mol2.getAtomWithIdx(atom2);
  bool atom1HasLabel = hasLabel(a1, autoLabels);
  bool atom2HasLabel = hasLabel(a2, autoLabels);
  // check for the presence of rgroup labels on adjacent terminal dummy atoms
  atom1HasLabel |= hasAttachedLabels(mol1, a1, autoLabels);
  atom2HasLabel |= hasAttachedLabels(mol2, a2, autoLabels);
  return !(atom1HasLabel != atom2HasLabel);
}

bool RGroupDecompositionParameters::prepareCore(RWMol &core,
                                                const RWMol *alignCore) {
  const bool relabel = labels & RelabelDuplicateLabels;
  unsigned int autoLabels = labels;
  if (labels == AutoDetect) {
    autoLabels = autoGetLabels(core);
    if (!autoLabels) {
      BOOST_LOG(rdWarningLog) << "RGroupDecomposition auto detect found no "
                                 "rgroups and onlyMatchAtRgroups is set to true"
                              << std::endl;
      return false;
    }
  } else if (!onlyMatchAtRGroups) {
    autoLabels |= AtomIndexLabels;
  }

  // if we aren't doing stereochem matches, remove that info from the core
  if (!substructmatchParams.useChirality) {
    for (auto atom : core.atoms()) {
      atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    }
    for (auto bond : core.bonds()) {
      bond->setStereo(Bond::BondStereo::STEREONONE);
    }
  }
  // remove enhanced stereo info if not being used
  // or if chirality isn't being used
  if (!substructmatchParams.useChirality ||
      !substructmatchParams.useEnhancedStereo) {
    core.setStereoGroups(std::vector<StereoGroup>());
  }

  int maxLabel = 1;
  // makes no sense to do MCS alignment if we are only matching at user defined
  // R-Groups
  if (alignCore && !onlyMatchAtRGroups && (alignment & MCS)) {
    std::vector<ROMOL_SPTR> mols;
    mols.push_back(ROMOL_SPTR(new ROMol(core)));
    mols.push_back(ROMOL_SPTR(new ROMol(*alignCore)));
    MCSParameters mcsParams;
    if (autoLabels != AutoDetect) {
      mcsParams.AtomTyper = rgdAtomCompare;
      mcsParams.CompareFunctionsUserData = &autoLabels;
    }
    MCSResult res = findMCS(mols, &mcsParams);
    if (res.isCompleted()) {
      auto m = res.QueryMol;
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
      }
    }
  }
  std::set<int> foundLabels;

  std::map<int, int> atomToLabel;

  for (auto atom : core.atoms()) {
    bool found = false;

    if (atom->hasProp(RLABEL)) {
      // set from MCS match
      if (setLabel(atom, atom->getProp<int>(RLABEL), foundLabels, maxLabel,
                   relabel, Labelling::INTERNAL_LABELS)) {
        found = true;
      }
    }

    if (atom->getAtomicNum() == 0) {
      if (!found && (autoLabels & MDLRGroupLabels)) {
        unsigned int rgroup;
        if (atom->getPropIfPresent<unsigned int>(
                common_properties::_MolFileRLabel, rgroup)) {
          if (setLabel(atom, rdcast<int>(rgroup), foundLabels, maxLabel,
                       relabel, Labelling::RGROUP_LABELS)) {
            found = true;
            checkNonTerminal(*atom);
          }
        }
      }

      if (!found && (autoLabels & IsotopeLabels) && atom->getIsotope() > 0) {
        if (setLabel(atom, rdcast<int>(atom->getIsotope()), foundLabels,
                     maxLabel, relabel, Labelling::ISOTOPE_LABELS)) {
          checkNonTerminal(*atom);
          found = true;
        }
      }

      if (!found && (autoLabels & AtomMapLabels) && atom->getAtomMapNum() > 0) {
        if (setLabel(atom, rdcast<int>(atom->getAtomMapNum()), foundLabels,
                     maxLabel, relabel, Labelling::ATOMMAP_LABELS)) {
          checkNonTerminal(*atom);
          found = true;
        }
      }

      if (!found && (autoLabels & DummyAtomLabels) && atom->getDegree() == 1 &&
          !atom->hasProp(UNLABELED_CORE_ATTACHMENT)) {
        const bool forceRelabellingWithDummies = true;
        int defaultDummyStartLabel = maxLabel;
        if (setLabel(atom, defaultDummyStartLabel, foundLabels, maxLabel,
                     forceRelabellingWithDummies, Labelling::DUMMY_LABELS)) {
          found = true;
        }
      }
    }

    // Unless there is an MCS match from above, we need to give different
    //  RLABELS to each core so keep track of which labels
    //  we have used (note that these are negative since they are
    //  potential rgroups and haven't been assigned yet)
    if (!found && (autoLabels & AtomIndexLabels)) {
      // we should not need these on (non r group) core atoms when
      // allowMultipleRGroupsOnUnlabelled is set, but it is useful in case
      // insufficient dummy groups are added to the core
      if (setLabel(atom, indexOffset - atom->getIdx(), foundLabels, maxLabel,
                   relabel, Labelling::INDEX_LABELS)) {
      }
      found = true;
    }

    clearInputLabels(atom);

    int rlabel;
    if (atom->getPropIfPresent(RLABEL, rlabel)) {
      atomToLabel[atom->getIdx()] = rlabel;
    }
  }
  indexOffset -= core.getNumAtoms();

  MolOps::AdjustQueryParameters adjustParams;
  adjustParams.makeDummiesQueries = true;
  adjustParams.adjustDegree = false;
  adjustQueryProperties(core, &adjustParams);
  for (auto &it : atomToLabel) {
    core.getAtomWithIdx(it.first)->setProp(RLABEL, it.second);
  }

  return true;
}  // namespace RDKit

void RGroupDecompositionParameters::checkNonTerminal(const Atom &atom) const {
  if (allowNonTerminalRGroups || atom.getDegree() == 1) {
    return;
  }

  BOOST_LOG(rdWarningLog)
      << "Non terminal R group defined.  To allow set allowNonTerminalRGroups "
         "in RGroupDecompositionParameters"
      << std::endl;
  throw ValueErrorException("Non terminal R group defined.");
}

void RGroupDecompositionParameters::addDummyAtomsToUnlabelledCoreAtoms(
    RWMol &core) {
  if (!allowMultipleRGroupsOnUnlabelled) {
    return;
  }

  // add R group substitutions to fill atomic valence
  std::vector<Atom *> unlabeledCoreAtoms{};
  for (const auto atom : core.atoms()) {
    if (atom->getAtomicNum() == 1) {
      continue;
    }
    if (hasLabel(atom, labels)) {
      continue;
    }
    unlabeledCoreAtoms.push_back(atom);
  }

  for (const auto atom : unlabeledCoreAtoms) {
    atom->calcImplicitValence(false);
    const auto atomIndex = atom->getIdx();
    int dummiesToAdd;
    int maxNumDummies = 4 - static_cast<int>(atom->getDegree());

    // figure out the number of dummies to add
    if (atom->getAtomicNum() == 0) {
      dummiesToAdd = maxNumDummies;
    } else {
      double bondOrder = 0;
      for (const auto bond : core.atomBonds(atom)) {
        auto contrib = bond->getValenceContrib(atom);
        if (contrib == 0.0 && bond->hasQuery()) {
          contrib = 1.0;
        }
        bondOrder += contrib;
      }

      const auto &valances =
          PeriodicTable::getTable()->getValenceList(atom->getAtomicNum());
      auto valence = *std::max_element(valances.begin(), valances.end());
      // round up aromatic contributions
      dummiesToAdd = valence - (int)(bondOrder + .51);
      dummiesToAdd = std::min(dummiesToAdd, maxNumDummies);
    }

    std::vector<int> newIndices;
    for (int i = 0; i < dummiesToAdd; i++) {
      const auto newAtom = new Atom(0);
      newAtom->setProp<bool>(UNLABELED_CORE_ATTACHMENT, true);
      const auto newIdx = core.addAtom(newAtom, false, true);
      newIndices.push_back(newIdx);
      auto *qb = new QueryBond();
      qb->setQuery(makeBondNullQuery());
      qb->setBeginAtomIdx(atomIndex);
      qb->setEndAtomIdx(newIdx);
      core.addBond(qb, true);
      const auto dummy = core.getAtomWithIdx(newIdx);
      dummy->updatePropertyCache();
    }
    atom->updatePropertyCache(false);
    for (const auto newIdx : newIndices) {
      if (core.getNumConformers() > 0) {
        MolOps::setTerminalAtomCoords(core, newIdx, atomIndex);
      }
    }
  }
}

}  // namespace RDKit
