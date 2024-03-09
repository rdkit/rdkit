//
//  Copyright (c) 2017-2023, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RGroupDecompData.h"

#include "RGroupDecomp.h"
#include "RGroupMatch.h"
#include "RGroupGa.h"
#include "GraphMol/MolEnumerator/MolEnumerator.h"

// #define VERBOSE 1

namespace RDKit {

RGroupDecompData::RGroupDecompData(const RWMol &inputCore,
                                   RGroupDecompositionParameters inputParams)
    : params(std::move(inputParams)) {
  addInputCore(inputCore);
  prepareCores();
}

RGroupDecompData::RGroupDecompData(const std::vector<ROMOL_SPTR> &inputCores,
                                   RGroupDecompositionParameters inputParams)
    : params(std::move(inputParams)) {
  for (const auto &core : inputCores) {
    addInputCore(*core);
  }
  prepareCores();
}

void RGroupDecompData::addInputCore(const ROMol& inputCore) {
  if (params.doEnumeration) {
    if (const auto bundle = MolEnumerator::enumerate(inputCore);
        !bundle.empty()) {
      for (auto c : bundle.getMols()) {
        addCore(*c);
      }
    } else {
      addCore(inputCore);
    }
  } else  {
    addCore(inputCore);
  }
}


void RGroupDecompData::addCore(const ROMol &inputCore) {
  if (params.allowMultipleRGroupsOnUnlabelled && !params.onlyMatchAtRGroups) {
    RWMol core(inputCore);
    params.addDummyAtomsToUnlabelledCoreAtoms(core);
    cores[cores.size()] = RCore(core);
  } else {
    cores[cores.size()] = RCore(inputCore);
  }
}

void RGroupDecompData::prepareCores() {
  for (auto &core : cores) {
    RWMol *alignCore = core.first ? cores[0].core.get() : nullptr;
    CHECK_INVARIANT(params.prepareCore(*core.second.core, alignCore),
                    "Could not prepare at least one core");
    core.second.init();
    core.second.labelledCore.reset(new RWMol(*core.second.core));
  }
}

void RGroupDecompData::setRlabel(Atom *atom, int rlabel) {
  PRECONDITION(rlabel > 0, "RLabels must be >0");
  if (params.rgroupLabelling & AtomMap) {
    atom->setAtomMapNum(rlabel);
  }

  if (params.rgroupLabelling & MDLRGroup) {
    std::string dLabel = "R" + std::to_string(rlabel);
    atom->setProp(common_properties::dummyLabel, dLabel);
    setAtomRLabel(atom, rlabel);
  } else {
    atom->clearProp(common_properties::dummyLabel);
  }

  if (params.rgroupLabelling & Isotope) {
    atom->setIsotope(rlabel);
  }
}

int RGroupDecompData::getRlabel(Atom *atom) const {
  if (params.rgroupLabelling & AtomMap) {
    return atom->getAtomMapNum();
  }
  if (params.rgroupLabelling & Isotope) {
    return atom->getIsotope();
  }

  if (params.rgroupLabelling & MDLRGroup) {
    unsigned int label = 0;
    if (atom->getPropIfPresent(common_properties::_MolFileRLabel, label)) {
      return label;
    }
  }

  CHECK_INVARIANT(0, "no valid r label found");
}

double RGroupDecompData::scoreFromPrunedData(
    const std::vector<size_t> &permutation, bool reset) {
  PRECONDITION(
      static_cast<RGroupScore>(params.scoreMethod) == FingerprintVariance,
      "Scoring method is not fingerprint variance!");

  PRECONDITION(permutation.size() >= pruneLength,
               "Illegal permutation prune length");
  if (permutation.size() < pruneLength * 1.5) {
    for (unsigned int pos = pruneLength; pos < permutation.size(); ++pos) {
      prunedFingerprintVarianceScoreData.addVarianceData(pos, permutation[pos],
                                                         matches, labels);
    }
    double score =
        prunedFingerprintVarianceScoreData.fingerprintVarianceGroupScore();
    if (reset) {
      for (unsigned int pos = pruneLength; pos < permutation.size(); ++pos) {
        prunedFingerprintVarianceScoreData.removeVarianceData(
            pos, permutation[pos], matches, labels);
      }
    } else {
      pruneLength = permutation.size();
    }
    return score;
  } else {
    if (reset) {
      return fingerprintVarianceScore(permutation, matches, labels);
    } else {
      prunedFingerprintVarianceScoreData.clear();
      pruneLength = permutation.size();
      return fingerprintVarianceScore(permutation, matches, labels,
                                      &prunedFingerprintVarianceScoreData);
    }
  }
}

void RGroupDecompData::prune() {  // prune all but the current "best"
                                  // permutation of matches
  PRECONDITION(permutation.size() <= matches.size(),
               "permutation.size() should be <= matches.size()");
  size_t offset = matches.size() - permutation.size();
  for (size_t mol_idx = 0; mol_idx < permutation.size(); ++mol_idx) {
    std::vector<RGroupMatch> keepVector;
    size_t mi = mol_idx + offset;
    keepVector.push_back(matches[mi].at(permutation[mol_idx]));
    matches[mi] = keepVector;
  }

  permutation = std::vector<size_t>(permutation.size(), 0);
  if (params.scoreMethod == FingerprintVariance &&
      params.matchingStrategy != GA) {
    scoreFromPrunedData(permutation, false);
  }
}

// Return the RGroups with the current "best" permutation
//  of matches.
std::vector<RGroupMatch> RGroupDecompData::GetCurrentBestPermutation() const {
  const bool removeAllHydrogenRGroups =
      params.removeAllHydrogenRGroups ||
      params.removeAllHydrogenRGroupsAndLabels;

  std::vector<RGroupMatch> results;  // std::map<int, RGroup> > result;
  bool isPruned = (permutation.size() < matches.size());
  for (size_t i = 0; i < matches.size(); ++i) {
    size_t pi = (isPruned ? 0 : permutation.at(i));
    results.push_back(matches[i].at(pi));
  }

  // * if a dynamically-added RGroup (i.e., when onlyMatchAtRGroups=false)
  //   is all hydrogens, remove it
  // * if a user-defined RGroup is all hydrogens and either
  //   params.removeAllHydrogenRGroups==true or
  //   params.removeAllHydrogenRGroupsAndLabels==true, remove it

  // This logic is a bit tricky, find all labels that have common cores
  //  and analyze those sets independently.
  //  i.e. if core 1 doesn't have R1 then don't analyze it in when looking
  //  at label 1
  std::map<int, std::set<int>> labelCores;  // map from label->cores
  std::set<int> coresVisited;
  for (auto &position : results) {
    int core_idx = position.core_idx;
    if (coresVisited.find(core_idx) == coresVisited.end()) {
      coresVisited.insert(core_idx);
      auto core = cores.find(core_idx);
      if (core != cores.end()) {
        for (auto rlabels : getRlabels(*core->second.core)) {
          int rlabel = rlabels.first;
          labelCores[rlabel].insert(core_idx);
        }
      }
    }
  }

  std::set<int> labelsToErase;
  for (int label : labels) {
    if (label > 0 && !removeAllHydrogenRGroups) {
      continue;
    }
    bool allH = true;
    for (auto &position : results) {
      R_DECOMP::const_iterator rgroup = position.rgroups.find(label);
      bool labelHasCore =
          labelCores[label].find(position.core_idx) != labelCores[label].end();
      if (labelHasCore && rgroup != position.rgroups.end() &&
          !rgroup->second->is_hydrogen) {
        allH = false;
        break;
      }
    }

    if (allH) {
      labelsToErase.insert(label);
      for (auto &position : results) {
        position.rgroups.erase(label);
      }
    }
  }

  for (auto &position : results) {
    for (auto atom : position.matchedCore->atoms()) {
      if (int atomLabel; atom->getAtomicNum() == 0 &&
                         atom->getPropIfPresent(RLABEL, atomLabel)) {
        if (atomLabel > 0 && !params.removeAllHydrogenRGroupsAndLabels) {
          continue;
        }
        if (labelsToErase.find(atomLabel) != labelsToErase.end()) {
          atom->setAtomicNum(1);
          atom->clearProp(RLABEL);
          if (atom->hasProp(RLABEL_TYPE)) {
            atom->clearProp(RLABEL_TYPE);
          }
          if (atom->hasProp(UNLABELED_CORE_ATTACHMENT)) {
            atom->clearProp(UNLABELED_CORE_ATTACHMENT);
          }
          atom->updatePropertyCache(false);
        }
      }
    }
  }

  return results;
}

bool RGroupDecompData::UsedLabels::add(int rlabel) {
  if (labels_used.find(rlabel) != labels_used.end()) {
    return false;
  }
  labels_used.insert(rlabel);
  return true;
}

int RGroupDecompData::UsedLabels::next() {
  int i = 1;
  while (labels_used.find(i) != labels_used.end()) {
    ++i;
  }
  labels_used.insert(i);
  return i;
}

void RGroupDecompData::addCoreUserLabels(const RWMol &core,
                                         std::set<int> &userLabels) {
  auto atoms = getRlabels(core);
  for (const auto &p : atoms) {
    if (p.first > 0) {
      userLabels.insert(p.first);
    }
  }
}

void RGroupDecompData::addAtoms(
    RWMol &mol, const std::vector<std::pair<Atom *, Atom *>> &atomsToAdd) {
  for (const auto &i : atomsToAdd) {
    mol.addAtom(i.second, false, true);
    mol.addBond(i.first, i.second, Bond::SINGLE);
    if (mol.getNumConformers()) {
      MolOps::setTerminalAtomCoords(mol, i.second->getIdx(), i.first->getIdx());
    }
  }
}

bool RGroupDecompData::replaceHydrogenCoreDummy(const RGroupMatch &match,
                                                RWMol &core, const Atom &atom,
                                                const int currentLabel,
                                                const int rLabel) {
  // if the R group is just a hydrogen then the attachment point should
  // replace an existing hydrogen neighbor since all hydrogen neighbors
  // are copied from the input molecule to the extracted core.
  if (const auto group = match.rgroups.find(currentLabel);
      group != match.rgroups.end()) {
    if (group->second->is_hydrogen) {
      for (auto &neighbor : core.atomNeighbors(&atom)) {
        if (neighbor->getAtomicNum() == 1) {
          neighbor->setAtomicNum(0);
          setRlabel(neighbor, rLabel);
          return true;
        }
      }
    }
  }
  return false;
}

void RGroupDecompData::relabelCore(
    RWMol &core, std::map<int, int> &mappings, UsedLabels &used_labels,
    const std::set<int> &indexLabels,
    const std::map<int, std::vector<int>> &extraAtomRLabels,
    const RGroupMatch *const match) {
  // Now remap to proper rlabel ids
  //  if labels are positive, they come from User labels
  //  if they are negative, they come from indices and should be
  //  numbered *after* the user labels.
  //
  //  Some indices are attached to multiple bonds,
  //   these rlabels should be incrementally added last
  std::map<int, Atom *> atoms = getRlabels(core);
  // a core only has one labelled index
  //  a secondary structure extraAtomRLabels contains the number
  //  of bonds between this atom and the side chain

  // a sidechain atom has a vector of the attachments back to the
  //  core that takes the place of numBondsToRlabel

  std::vector<std::pair<Atom *, Atom *>> atomsToAdd;  // adds -R if necessary

  // Deal with user supplied labels
  for (const auto &rlabels : atoms) {
    int userLabel = rlabels.first;
    if (userLabel < 0) {
      continue;  // not a user specified label
    }
    Atom *atom = rlabels.second;
    mappings[userLabel] = userLabel;
    used_labels.add(userLabel);

    if (atom->getAtomicNum() == 0 &&
        atom->getDegree() == 1) {  // add to existing dummy/rlabel
      setRlabel(atom, userLabel);
    } else {
      // A non-terminal RGroup. Create a dummy by replacing an existing
      // hydrogen for hydrogen side chains, or a new dummy atom for heavy side
      // chains.
      bool addNew = true;
      if (match != nullptr) {
        addNew = !replaceHydrogenCoreDummy(*match, core, *atom, userLabel,
                                           userLabel);
        // If we can't replace a hydrogen only add the dummy if it exists in
        // the decomp This is unexpected.
        if (addNew && match->rgroups.find(userLabel) == match->rgroups.end()) {
          addNew = false;
        }
      }
      if (addNew) {
        auto *newAt = new Atom(0);
        setRlabel(newAt, userLabel);
        atomsToAdd.emplace_back(atom, newAt);
      }
    }
  }

  // Deal with non-user supplied labels
  for (auto newLabel : indexLabels) {
    auto atm = atoms.find(newLabel);
    if (atm == atoms.end()) {
      continue;
    }

    Atom *atom = atm->second;

    int rlabel;
    auto mapping = mappings.find(newLabel);
    if (mapping == mappings.end()) {
      rlabel = used_labels.next();
      mappings[newLabel] = rlabel;
    } else {
      rlabel = mapping->second;
    }

    if (atom->getAtomicNum() == 0 &&
        !isAnyAtomWithMultipleNeighborsOrNotUserRLabel(
            *atom)) {  // add to dummy
      setRlabel(atom, rlabel);
    } else {
      bool addNew = true;
      if (match != nullptr) {
        addNew =
            !replaceHydrogenCoreDummy(*match, core, *atom, newLabel, rlabel);
        // If we can't replace a hydrogen only add the dummy if it exists in
        // the decomp This is unexpected.
        if (addNew && match->rgroups.find(newLabel) == match->rgroups.end()) {
          addNew = false;
        }
      }
      if (addNew) {
        auto *newAt = new Atom(0);
        setRlabel(newAt, rlabel);
        atomsToAdd.emplace_back(atom, newAt);
      }
    }
  }

  // Deal with multiple bonds to the same label
  for (const auto &extraAtomRLabel : extraAtomRLabels) {
    auto atm = atoms.find(extraAtomRLabel.first);
    if (atm == atoms.end()) {
      continue;  // label not used in the rgroup
    }
    Atom *atom = atm->second;

    for (size_t i = 0; i < extraAtomRLabel.second.size(); ++i) {
      int rlabel = used_labels.next();
      // Is this necessary?
      CHECK_INVARIANT(
          atom->getAtomicNum() > 1,
          "Multiple attachments to a dummy (or hydrogen) is weird.");
      auto *newAt = new Atom(0);
      setRlabel(newAt, rlabel);
      atomsToAdd.emplace_back(atom, newAt);
    }
  }

  addAtoms(core, atomsToAdd);
  for (const auto &rlabels : atoms) {
    auto atom = rlabels.second;
    atom->clearProp(RLABEL);
    atom->clearProp(RLABEL_TYPE);
  }

  // Delay removing hydrogens from core until outputCoreMolecule is called,
  // If hydrogens are removed now and more dummies removed in
  // outputCoreMolecule then aromaticity perception in the core may be broken.

  core.updatePropertyCache(false);  // this was github #1550
}

void RGroupDecompData::relabelRGroup(RGroupData &rgroup,
                                     const std::map<int, int> &mappings) {
  PRECONDITION(rgroup.combinedMol.get(), "Unprocessed rgroup");

  RWMol &mol = *rgroup.combinedMol.get();

  if (rgroup.combinedMol->hasProp(done)) {
    rgroup.labelled = true;
    return;
  }

  mol.setProp(done, true);
  std::vector<std::pair<Atom *, Atom *>> atomsToAdd;  // adds -R if necessary
  std::map<int, int> rLabelCoreIndexToAtomicWt;

  for (RWMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    Atom *atom = *atIt;
    if (atom->hasProp(SIDECHAIN_RLABELS)) {
      atom->setIsotope(0);
      const std::vector<int> &rlabels =
          atom->getProp<std::vector<int>>(SIDECHAIN_RLABELS);
      // switch on atom mappings or rlabels....

      for (int rlabel : rlabels) {
        auto label = mappings.find(rlabel);
        CHECK_INVARIANT(label != mappings.end(), "Unprocessed mapping");

        if (atom->getAtomicNum() == 0) {
          if (!atom->hasProp(_rgroupInputDummy)) {
            setRlabel(atom, label->second);
          }
        } else if (atom->hasProp(RLABEL_CORE_INDEX)) {
          atom->setAtomicNum(0);
          setRlabel(atom, label->second);
        } else {
          auto *newAt = new Atom(0);
          setRlabel(newAt, label->second);
          atomsToAdd.emplace_back(atom, newAt);
        }
      }
    }
    if (atom->hasProp(RLABEL_CORE_INDEX)) {
      // convert to dummy as we don't want to collapse hydrogens onto the core
      // match
      auto rLabelCoreIndex = atom->getProp<int>(RLABEL_CORE_INDEX);
      rLabelCoreIndexToAtomicWt[rLabelCoreIndex] = atom->getAtomicNum();
      atom->setAtomicNum(0);
    }
  }

  addAtoms(mol, atomsToAdd);

  if (params.removeHydrogensPostMatch) {
    RDLog::LogStateSetter blocker;
    bool implicitOnly = false;
    bool updateExplicitCount = false;
    bool sanitize = false;
    MolOps::removeHs(mol, implicitOnly, updateExplicitCount, sanitize);
  }

  mol.updatePropertyCache(false);  // this was github #1550
  rgroup.labelled = true;

  // Restore any core matches that we have set to dummy
  for (RWMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    Atom *atom = *atIt;
    if (atom->hasProp(RLABEL_CORE_INDEX)) {
      // don't need to set IsAromatic on atom - that seems to have been saved
      atom->setAtomicNum(
          rLabelCoreIndexToAtomicWt[atom->getProp<int>(RLABEL_CORE_INDEX)]);
      atom->setNoImplicit(true);
      atom->clearProp(RLABEL_CORE_INDEX);
    }
    atom->clearProp(SIDECHAIN_RLABELS);
  }

#ifdef VERBOSE
  std::cerr << "Relabel Rgroup smiles " << MolToSmiles(mol) << std::endl;
#endif
}

// relabel the core and sidechains using the specified user labels
//  if matches exist for non labelled atoms, these are added as well
void RGroupDecompData::relabel() {
  std::vector<RGroupMatch> best = GetCurrentBestPermutation();

  // get the labels used
  std::set<int> userLabels;
  std::set<int> indexLabels;

  // Go through all the RGroups and find out which labels were
  //  actually used.

  // some atoms will have multiple attachment points, i.e. cycles
  //  split these up into new rlabels if necessary
  //  These are detected at match time
  //  This vector will hold the extra (new) labels required
  std::map<int, std::vector<int>> extraAtomRLabels;

  for (auto &it : best) {
    for (auto &rgroup : it.rgroups) {
      if (rgroup.first > 0) {
        userLabels.insert(rgroup.first);
      }
      if (rgroup.first < 0 && !params.onlyMatchAtRGroups) {
        indexLabels.insert(rgroup.first);
      }

      std::map<int, int> rlabelsUsedInRGroup =
          rgroup.second->getNumBondsToRlabels();
      for (auto &numBondsUsed : rlabelsUsedInRGroup) {
        // Make space for the extra labels
        if (numBondsUsed.second > 1) {  // multiple rgroup bonds to same atom
          extraAtomRLabels[numBondsUsed.first].resize(numBondsUsed.second - 1);
        }
      }
    }
  }

  // find user labels that are not present in the decomposition
  for (auto &core : cores) {
    core.second.labelledCore.reset(new RWMol(*core.second.core));
    addCoreUserLabels(*core.second.labelledCore, userLabels);
  }

  // Assign final RGroup labels to the cores and propagate these to
  //  the scaffold
  finalRlabelMapping.clear();

  UsedLabels used_labels;
  // Add all the user labels now to prevent an index label being assigned to a
  // user label when multiple cores are present (e.g. the user label is
  // present in the second core, but not the first).
  for (auto userLabel : userLabels) {
    used_labels.add(userLabel);
  }
  for (auto &core : cores) {
    relabelCore(*core.second.labelledCore, finalRlabelMapping, used_labels,
                indexLabels, extraAtomRLabels);
  }

  for (auto &it : best) {
    for (auto &rgroup : it.rgroups) {
      relabelRGroup(*rgroup.second, finalRlabelMapping);
    }
#ifdef VERBOSE
    std::cerr << "relabel core mol1 " << MolToSmiles(*it.matchedCore)
              << std::endl;
#endif
    relabelCore(*it.matchedCore, finalRlabelMapping, used_labels, indexLabels,
                extraAtomRLabels, &it);
#ifdef VERBOSE
    std::cerr << "relabel core mol2 " << MolToSmiles(*it.matchedCore)
              << std::endl;
#endif
  }

  std::set<int> uniqueMappedValues;
  std::transform(finalRlabelMapping.cbegin(), finalRlabelMapping.cend(),
                 std::inserter(uniqueMappedValues, uniqueMappedValues.end()),
                 [](const std::pair<int, int> &p) { return p.second; });
  CHECK_INVARIANT(finalRlabelMapping.size() == uniqueMappedValues.size(),
                  "Error in uniqueness of final RLabel mapping");
  CHECK_INVARIANT(
      uniqueMappedValues.size() == userLabels.size() + indexLabels.size(),
      "Error in final RMapping size");
}

double RGroupDecompData::score(
    const std::vector<size_t> &permutation,
    FingerprintVarianceScoreData *fingerprintVarianceScoreData) const {
  RGroupScore scoreMethod = static_cast<RGroupScore>(params.scoreMethod);
  switch (scoreMethod) {
    case Match:
      return rGroupScorer.matchScore(permutation, matches, labels);
      break;
    case FingerprintVariance:
      return fingerprintVarianceScore(permutation, matches, labels,
                                      fingerprintVarianceScoreData);
      break;
    default:;
  }
  return NAN;
}

RGroupDecompositionProcessResult RGroupDecompData::process(bool pruneMatches,
                                                           bool finalize) {
  if (matches.empty()) {
    return RGroupDecompositionProcessResult(false, -1);
  }
  auto t0 = std::chrono::steady_clock::now();
  std::unique_ptr<CartesianProduct> iterator;
  rGroupScorer.startProcessing();

  if (params.matchingStrategy == GA) {
    RGroupGa ga(*this, params.timeout >= 0 ? &t0 : nullptr);
    if (ga.numberPermutations() < 100 * ga.getPopsize()) {
      params.matchingStrategy = Exhaustive;
    } else {
      if (params.gaNumberRuns > 1) {
        auto results = ga.runBatch();
        auto best = max_element(results.begin(), results.end(),
                                [](const GaResult &a, const GaResult &b) {
                                  return a.rGroupScorer.getBestScore() <
                                         b.rGroupScorer.getBestScore();
                                });
        rGroupScorer = best->rGroupScorer;
      } else {
        auto result = ga.run();
        rGroupScorer = result.rGroupScorer;
      }
    }
  }
  size_t offset = 0;
  if (params.matchingStrategy != GA) {
    // Exhaustive search, get the MxN matrix
    // (M = matches.size(): number of molecules
    //  N = iterator.maxPermutations)
    std::vector<size_t> permutations;

    if (pruneMatches && params.scoreMethod != FingerprintVariance) {
      offset = previousMatchSize;
    }
    previousMatchSize = matches.size();
    permutations.reserve(matches.size() - offset);
    std::transform(matches.begin() + offset, matches.end(),
                   std::back_inserter(permutations),
                   [](const std::vector<RGroupMatch> &m) { return m.size(); });
    permutation = std::vector<size_t>(permutations.size(), 0);

    // run through all possible matches and score each set
    size_t count = 0;
#ifdef DEBUG
    std::cerr << "Processing" << std::endl;
#endif
    iterator.reset(new CartesianProduct(permutations));
    // Iterates through the permutation idx, i.e.
    //  [m1_permutation_idx,  m2_permutation_idx, m3_permutation_idx]

    while (iterator->next()) {
      if (count > iterator->maxPermutations) {
        throw ValueErrorException("next() did not finish");
      }
#ifdef DEBUG
      std::cerr << "**************************************************"
                << std::endl;
#endif
      double newscore = params.scoreMethod == FingerprintVariance
                            ? scoreFromPrunedData(iterator->permutation)
                            : score(iterator->permutation);

      if (fabs(newscore - rGroupScorer.getBestScore()) <
          1e-6) {  // heuristic to overcome floating point comparison issues
        rGroupScorer.pushTieToStore(iterator->permutation);
      } else if (newscore > rGroupScorer.getBestScore()) {
#ifdef DEBUG
        std::cerr << " ===> current best:" << newscore << ">"
                  << rGroupScorer.getBestScore() << std::endl;
#endif
        rGroupScorer.setBestPermutation(iterator->permutation, newscore);
        rGroupScorer.clearTieStore();
        rGroupScorer.pushTieToStore(iterator->permutation);
      }
      ++count;
    }
  }

  if (rGroupScorer.tieStoreSize() > 1) {
    rGroupScorer.breakTies(matches, labels, iterator, t0, params.timeout);
    rGroupScorer.clearTieStore();
  } else {
    checkForTimeout(t0, params.timeout);
  }
  permutation = rGroupScorer.getBestPermutation();
  if (pruneMatches || finalize) {
    prune();
  }

  if (finalize) {
    relabel();
  }

  return RGroupDecompositionProcessResult(true, rGroupScorer.getBestScore());
}

}  // namespace RDKit
