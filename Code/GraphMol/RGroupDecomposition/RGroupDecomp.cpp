//
//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "RGroupDecomp.h"
#include "RGroupDecompData.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <boost/dynamic_bitset.hpp>
#include <set>
#include <utility>
#include <vector>

#include "GraphMol/TautomerQuery/TautomerQuery.h"

// #define VERBOSE 1

namespace RDKit {

// Attachment Points
//  labeled cores => isotopes
//  atom mappings
//  atom indices => use -1 - atom index, range is [-1, ...., -num_atoms]
const std::string RLABEL = "tempRlabel";
const std::string RLABEL_TYPE = "tempRlabelType";
const std::string RLABEL_CORE_INDEX = "rLabelCoreIndex";
const std::string SIDECHAIN_RLABELS = "sideChainRlabels";
const std::string done = "RLABEL_PROCESSED";
const std::string CORE = "Core";
const std::string RPREFIX = "R";
const std::string _rgroupInputDummy = "_rgroupInputDummy";
const std::string UNLABELED_CORE_ATTACHMENT = "unlabeledCoreAttachment";

namespace {
void ADD_MATCH(R_DECOMP &match, int rlabel) {
  if (match.find(rlabel) == match.end()) {
    match[rlabel] = boost::make_shared<RGroupData>();
  }
}
}  // namespace

RGroupDecomposition::RGroupDecomposition(
    const ROMol &inputCore, const RGroupDecompositionParameters &params)
    : data(new RGroupDecompData(inputCore, params)) {}

RGroupDecomposition::RGroupDecomposition(
    const std::vector<ROMOL_SPTR> &cores,
    const RGroupDecompositionParameters &params)
    : data(new RGroupDecompData(cores, params)) {}

RGroupDecomposition::~RGroupDecomposition() { delete data; }

int RGroupDecomposition::getMatchingCoreIdx(
    const ROMol &mol, std::vector<MatchVectType> *matches) {
  RWMol rwmol(mol);
  std::vector<MatchVectType> matchesTmp;
  const RCore *rcore;
  auto coreIdx = getMatchingCoreInternal(rwmol, rcore, matchesTmp);
  if (matches) {
    std::set<MatchVectType> uniqueMatches;
    int numAtoms = mol.getNumAtoms();
    for (const auto &match : matchesTmp) {
      MatchVectType heavyMatch;
      heavyMatch.reserve(match.size());
      std::copy_if(std::make_move_iterator(match.begin()),
          std::make_move_iterator(match.end()), std::back_inserter(heavyMatch),
          [numAtoms](const auto &pair) { return pair.second < numAtoms; });
      std::sort(heavyMatch.begin(), heavyMatch.end());
      uniqueMatches.insert(heavyMatch);
    }
    *matches =
        std::vector<MatchVectType>(uniqueMatches.begin(), uniqueMatches.end());
  }
  return coreIdx;
}

int RGroupDecomposition::getMatchingCoreInternal(
    RWMol &mol, const RCore *&rcore, std::vector<MatchVectType> &matches) {
  rcore = nullptr;
  int core_idx = -1;
  const bool explicitOnly = false;
  const bool addCoords = true;
  MolOps::addHs(mol, explicitOnly, addCoords);
  std::vector<MatchVectType> tmatches;
  std::vector<MatchVectType> tmatches_filtered;

  // Find the first matching core (onlyMatchAtRGroups)
  // or the first core that requires the smallest number
  // of newly added labels and is a superstructure of
  // the first matching core
  int global_min_heavy_nbrs = -1;
  SubstructMatchParameters sssparams(params().substructmatchParams);
  sssparams.uniquify = false;
  sssparams.recursionPossible = true;
  for (auto &core : data->cores) {
    {
      // matching the core to the molecule is a two step process
      // First match to a reduced representation (the core minus terminal
      // R-groups). Next, match the R-groups. We do this as the core may not be
      // a substructure match for the molecule if a single molecule atom matches
      // 2 RGroup attachments (see https://github.com/rdkit/rdkit/pull/4002)

      // match the reduced representation:
      std::vector<MatchVectType> baseMatches;
      if (params().doTautomers) {
        // Here we are attempting to enumerate tautomers of the core
        if (auto tautomerQuery = core.second.getMatchingTautomerQuery();
            tautomerQuery != nullptr) {
          // query atom indices from the tautomer query are the same as the
          // template matching molecule
          baseMatches = tautomerQuery->substructOf(mol, sssparams);
        } else {
          // However, if it is not possible to Kekulize the core, we revert back
          // to the non-tautomer matching.
          baseMatches =
              SubstructMatch(mol, *core.second.matchingMol, sssparams);
        }
      } else {
        baseMatches = SubstructMatch(mol, *core.second.matchingMol, sssparams);
      }

      tmatches.clear();
      for (const auto &baseMatch : baseMatches) {
        // Match the R Groups
        // Important: there can be multiple core indices matching
        // the same target idx, because of #4002
        auto matchesIncludingRGroups =
            core.second.matchTerminalUserRGroups(mol, baseMatch, sssparams);
        /*
        std::cerr << "baseMatch ";
        for (const auto &pair : baseMatch) std::cerr << "(" << pair.first <<","
        << pair.second << "),"; std::cerr << std::endl; std::cerr <<
        "matchesIncludingRGroups "; for (const auto &matchWithDummy : matchesIncludingRGroups)
        { for (const auto &pair : matchWithDummy) std::cerr << "(" << pair.first
        <<"," << pair.second << "),"; std::cerr << " /// ";
        }
        std::cerr << std::endl;
        */
        tmatches.insert(tmatches.end(), std::make_move_iterator(matchesIncludingRGroups.cbegin()),
                        std::make_move_iterator(matchesIncludingRGroups.cend()));
      }
    }
    if (tmatches.empty()) {
      continue;
    }

    std::vector<int> tmatches_heavy_nbrs(tmatches.size(), 0);
    size_t i = 0;
    for (const auto &mv : tmatches) {
      bool passes_filter = data->params.onlyMatchAtRGroups;
      // targetToCoreIndices maps each atom idx in the molecule to a vector
      // of atom indices. This vector may be empty (if the atom in the molecule
      // has no match with core) or not. When not empty, it will most often
      // contain a single atom idx, corresponding to the matching index in the
      // core, as usually a core atom can only match a single molecule atom.
      // However, there is an important exception to this rule, i.e. when
      // the core bears a single R-group dummy at a certain position, while
      // the molecule has multiple substituents at the corresponding
      // position; in this case, the vector will contain the indices of the
      // root atom in all substituents which match a single R-group dummy on
      // the core.
      std::vector<std::vector<int>> targetToCoreIndices(mol.getNumAtoms());
      for (const auto &match : mv) {
        targetToCoreIndices[match.second].push_back(match.first);
      }

      for (const auto &match : mv) {
        const auto atm = mol.getAtomWithIdx(match.second);
        // is this a labelled rgroup or not?
        if (!core.second.isCoreAtomUserLabelled(match.first)) {
          // nope... if any neighbor is not part of the substructure
          // check if it is a hydrogen; otherwise, if onlyMatchAtRGroups
          // is true, skip the match
          for (const auto &nbri :
               boost::make_iterator_range(mol.getAtomNeighbors(atm))) {
            const auto &nbr = mol[nbri];
            if (nbr->getAtomicNum() != 1 &&
                targetToCoreIndices.at(nbr->getIdx()).empty()) {
              if (data->params.onlyMatchAtRGroups) {
                passes_filter = false;
                break;
              } else {
                // for each match, we keep track of the number of
                // R labels that need to be added to match all
                // non-user-labelled R groups in this molecule
                // if we use this core for RGD
                ++tmatches_heavy_nbrs[i];
              }
            }
          }
        } else if (core.second.isTerminalRGroupWithUserLabel(match.first)
              && data->params.onlyMatchAtRGroups && !core.second.checkAllBondsToRGroupPresent(
              mol, match.second, targetToCoreIndices)) {
          // labelled R-group
          passes_filter = false;
        }
        if (!passes_filter && data->params.onlyMatchAtRGroups) {
          break;
        }
      }

      if (passes_filter) {
        tmatches_filtered.push_back(std::move(mv));
      }
      ++i;
    }
    if (!data->params.onlyMatchAtRGroups) {
      // tmatches_heavy_nbrs.size() = tmatches.size(), and
      // tmatches.size() cannot be empty, otherwise we should not be here
      // but let's check it in case something changes upstream
      CHECK_INVARIANT(!tmatches_heavy_nbrs.empty(), "tmatches_heavy_nbrs must not be empty");
      int min_heavy_nbrs = *std::min_element(tmatches_heavy_nbrs.begin(),
                                             tmatches_heavy_nbrs.end());
      if (!rcore || (min_heavy_nbrs < global_min_heavy_nbrs &&
                     !SubstructMatch(*core.second.core, *rcore->core, sssparams)
                          .empty())) {
        i = 0;
        tmatches_filtered.clear();
        for (const auto heavy_nbrs : tmatches_heavy_nbrs) {
          if (heavy_nbrs <= min_heavy_nbrs) {
            tmatches_filtered.push_back(std::move(tmatches[i]));
          }
          ++i;
        }
        global_min_heavy_nbrs = min_heavy_nbrs;
        rcore = &core.second;
        core_idx = core.first;
        if (global_min_heavy_nbrs == 0) {
          break;
        }
      }
    } else if (!tmatches_filtered.empty()) {
      rcore = &core.second;
      core_idx = core.first;
      break;
    }
  }
  if (rcore) {
    matches = std::move(tmatches_filtered);
  }
  return core_idx;
}

int RGroupDecomposition::add(const ROMol &inmol) {
  // get the sidechains if possible
  //  Add hs for better symmetrization
  RWMol mol(inmol);
  const RCore *rcore;
  std::vector<MatchVectType> tmatches;
  auto core_idx = getMatchingCoreInternal(mol, rcore, tmatches);
  if (rcore == nullptr) {
    BOOST_LOG(rdDebugLog) << "No core matches" << std::endl;
    return -1;
  }

  if (tmatches.size() > 1) {
    if (data->params.matchingStrategy == NoSymmetrization) {
      tmatches.resize(1);
    } else if (data->matches.size() == 0) {
      // Greedy strategy just grabs the first match and
      //  takes the best matches from the rest
      if (data->params.matchingStrategy == Greedy) {
        tmatches.resize(1);
      }
    }
  }
  // mark any wildcards in input molecule:
  for (auto &atom : mol.atoms()) {
    if (atom->getAtomicNum() == 0) {
      atom->setProp(_rgroupInputDummy, true);
      // clean any existing R group numbers
      atom->setIsotope(0);
      atom->setAtomMapNum(0);
      atom->clearProp(common_properties::_MolFileRLabel);
      atom->setProp(common_properties::dummyLabel, "*");
    }
  }

  // strategies
  // ==========
  // Exhaustive - saves all matches and optimizes later exhaustive
  //               May never finish due to combinatorial complexity
  // Greedy - matches to *FIRST* available match
  // GreedyChunks - default - process every N chunks

  //  Should probably scan all mols first to find match with
  //  smallest number of matches...
  std::vector<RGroupMatch> potentialMatches;

  std::unique_ptr<ROMol> tMol;
  for (const auto &tmatche : tmatches) {
    const bool replaceDummies = false;
    const bool labelByIndex = true;
    const bool requireDummyMatch = false;
    // TODO see if we need replaceCoreAtomsWithMolMatches or can just use rcore->core
    auto coreCopy = rcore->replaceCoreAtomsWithMolMatches(mol, tmatche);
    tMol.reset(replaceCore(mol, *coreCopy, tmatche, replaceDummies,
                           labelByIndex, requireDummyMatch));
#ifdef VERBOSE
    std::cerr << "Core Match core_idx " << core_idx << " idx "
              << data->matches.size() << ": " << MolToSmarts(*coreCopy)
              << std::endl;
#endif
    if (tMol) {
#ifdef VERBOSE
      std::cerr << "All Fragments " << MolToSmiles(*tMol) << std::endl;
#endif
      R_DECOMP match;
      // rlabel rgroups
      MOL_SPTR_VECT fragments = MolOps::getMolFrags(*tMol, false);
      std::set<int> coreAtomAnyMatched;
      for (size_t i = 0; i < fragments.size(); ++i) {
        std::vector<int> rlabelsOnSideChain;
        const auto &newMol = fragments[i];
        newMol->setProp<int>("core", core_idx);
        newMol->setProp<int>("idx", data->matches.size());
        newMol->setProp<int>("frag_idx", i);
#ifdef VERBOSE
        std::cerr << "Fragment " << MolToSmiles(*newMol) << std::endl;
#endif
        for (auto sideChainAtom : newMol->atoms()) {
          if (sideChainAtom->getAtomicNum() != 0) {
            // we are only interested in sidechain R group atoms
            continue;
          }
          if (!sideChainAtom->hasProp(_rgroupInputDummy)) {
            // this is the index of the core atom that the R group
            // atom is attached to
            unsigned int coreAtomIndex = sideChainAtom->getIsotope();
            int rlabel;
            auto coreAtom = rcore->core->getAtomWithIdx(coreAtomIndex);
            coreAtomAnyMatched.insert(coreAtomIndex);
            if (coreAtom->getPropIfPresent(RLABEL, rlabel)) {
              std::vector<int> rlabelsOnSideChainAtom;
              sideChainAtom->getPropIfPresent(SIDECHAIN_RLABELS, rlabelsOnSideChainAtom);
              rlabelsOnSideChainAtom.push_back(rlabel);
              sideChainAtom->setProp(SIDECHAIN_RLABELS, rlabelsOnSideChainAtom);

              data->labels.insert(rlabel);  // keep track of all labels used
              rlabelsOnSideChain.push_back(rlabel);
              if (const auto [bondIdx, end] = newMol->getAtomBonds(sideChainAtom);
                  bondIdx != end) {
                auto connectingBond = (*newMol)[*bondIdx];
                if (connectingBond->getStereo() >
                    Bond::BondStereo::STEREOANY) {
                  // TODO: how to handle bond stereo on rgroups connected to
                  // core by stereo double bonds
                  connectingBond->setStereo(Bond::BondStereo::STEREOANY);
                }
              }
            }
          } else {
            // restore input wildcard
            sideChainAtom->clearProp(_rgroupInputDummy);
          }
        }

        if (!rlabelsOnSideChain.empty()) {
#ifdef VERBOSE
          std::string newCoreSmi = MolToSmiles(*newMol, true);
#endif

          for (auto rlabel : rlabelsOnSideChain) {
            ADD_MATCH(match, rlabel);
            match[rlabel]->add(newMol, rlabelsOnSideChain);
#ifdef VERBOSE
            std::cerr << "Fragment " << i << " R" << rlabel << " "
                      << MolToSmiles(*newMol) << std::endl;
#endif
          }
        } else {
          // special case, only one fragment
          if (fragments.size() == 1) {  // need to make a new core
            // remove the sidechains

            // GJ I think if we ever get here that it's really an error and I
            // believe that I've fixed the case where this code was called.
            // Still, I'm too scared to delete the block.
            RWMol newCore(mol);

            for (const auto &mvpair : tmatche) {
              const Atom *coreAtm = rcore->core->getAtomWithIdx(mvpair.first);
              Atom *newCoreAtm = newCore.getAtomWithIdx(mvpair.second);
              int rlabel;
              if (coreAtm->getPropIfPresent(RLABEL, rlabel)) {
                newCoreAtm->setProp<int>(RLABEL, rlabel);
              }
              newCoreAtm->setProp<bool>("keep", true);
            }

            newCore.beginBatchEdit();
            for (const auto atom : newCore.atoms()) {
              if (!atom->hasProp("keep")) {
                newCore.removeAtom(atom);
              }
            }
            newCore.commitBatchEdit();
            if (newCore.getNumAtoms()) {
              std::string newCoreSmi = MolToSmiles(newCore, true);
              // add a new core if possible
              auto newcore = data->newCores.find(newCoreSmi);
              int core_idx = 0;
              if (newcore == data->newCores.end()) {
                core_idx = data->newCores[newCoreSmi] = data->newCoreLabel--;
                data->cores[core_idx] = RCore(newCore);
                return add(inmol);
              }
            }
          }
        }
      }

      if (!match.empty()) {
        // this is the number of user-defined R labels associated with
        // non-hydrogen substituents
        auto numberUserGroupsInMatch = std::accumulate(
            match.begin(), match.end(), 0,
            [](int sum, const std::pair<int, boost::shared_ptr<RGroupData>> &p) {
              return p.first > 0 && !p.second->is_hydrogen ? ++sum : sum;
            });
        int numberMissingUserGroups =
            rcore->numberUserRGroups - numberUserGroupsInMatch;
        CHECK_INVARIANT(numberMissingUserGroups >= 0,
                        "Data error in missing user rgroup count");
        const auto extractedCore =
            rcore->extractCoreFromMolMatch(mol, tmatche, params());
        potentialMatches.emplace_back(core_idx, numberMissingUserGroups, match,
                                      extractedCore);
      }
    }
  }
  if (potentialMatches.empty()) {
    BOOST_LOG(rdDebugLog) << "No attachment points in side chains" << std::endl;
    return -2;
  }

  // in case the value ends up being changed in a future version of the code:
  if (data->prunePermutations) {
    data->permutationProduct = 1;
  }

  data->matches.push_back(std::move(potentialMatches));

  if (!data->matches.empty()) {
    if (data->params.matchingStrategy & Greedy ||
        (data->params.matchingStrategy & GreedyChunks &&
         data->matches.size() % data->params.chunkSize == 0)) {
      data->process(data->prunePermutations);
    }
  }
  return data->matches.size() - 1;
}

bool RGroupDecomposition::process() { return processAndScore().success; }

RGroupDecompositionProcessResult RGroupDecomposition::processAndScore() {
  try {
    const bool finalize = true;
    return data->process(data->prunePermutations, finalize);
  } catch (...) {
    return RGroupDecompositionProcessResult(false, -1);
  }
}

std::vector<std::string> RGroupDecomposition::getRGroupLabels() const {
  // this is a bit of a cheat
  RGroupColumns cols = getRGroupsAsColumns();
  std::vector<std::string> labels;
  for (auto it : cols) {
    labels.push_back(it.first);
  }
  std::sort(labels.begin(), labels.end());
  return labels;
}

RWMOL_SPTR RGroupDecomposition::outputCoreMolecule(
    const RGroupMatch &match, const UsedLabelMap &usedLabelMap) const {
  // this routine could probably be merged into RGroupDecompData::relabelCore

  const auto &core = data->cores[match.core_idx];
  if (!match.matchedCore) {
    return core.labelledCore;
  }
  auto coreWithMatches = match.matchedCore;
#ifdef VERBOSE
  std::cerr << "output core mol1 " << MolToSmarts(*coreWithMatches)
            << std::endl;
#endif
  std::map<Atom *, int> retainedRGroups;
  for (auto atomIdx = coreWithMatches->getNumAtoms(); atomIdx--;) {
    auto atom = coreWithMatches->getAtomWithIdx(atomIdx);
    if (atom->getAtomicNum()) {
      continue;
    }
    auto label = data->getRlabel(atom);
    // Always convert to hydrogen - then remove later if
    // removeHydrogensPostMatch is set
    Atom *nbrAtom = nullptr;
    for (const auto &nbri :
         boost::make_iterator_range(coreWithMatches->getAtomNeighbors(atom))) {
      nbrAtom = (*coreWithMatches)[nbri];
      break;
    }
    if (nbrAtom) {
      const bool isUserDefinedLabel =
          usedLabelMap.has(label) && usedLabelMap.isUserDefined(label);
      const bool isUsedLabel =
          usedLabelMap.has(label) && usedLabelMap.getIsUsed(label);
      if (!isUsedLabel && (!isUserDefinedLabel ||
                           data->params.removeAllHydrogenRGroupsAndLabels)) {
        // Always convert to hydrogen - then remove later if
        // removeHydrogensPostMatch is set
        atom->setAtomicNum(1);
        atom->updatePropertyCache(false);
      } else {
        retainedRGroups[atom] = label;
      }
    }
  }

#ifdef VERBOSE
  std::cerr << "output core mol2 " << MolToSmiles(*coreWithMatches)
            << std::endl;
#endif
  if (data->params.removeHydrogensPostMatch) {
    RDLog::LogStateSetter blocker;
    constexpr bool implicitOnly = false;
    constexpr bool updateExplicitCount = false;
    constexpr bool sanitize = false;
    MolOps::removeHs(*coreWithMatches, implicitOnly, updateExplicitCount,
                     sanitize);
    coreWithMatches->updatePropertyCache(false);
  }

  if (coreWithMatches->getNumConformers() > 0) {
    for (const auto &[atom, label] : retainedRGroups) {
      if (usedLabelMap.has(label) && usedLabelMap.isUserDefined(label)) {
        // coordinates of user defined R groups should already be copied over
        continue;
      }
      const auto neighbor = *coreWithMatches->atomNeighbors(atom).begin();
      const auto &mapping = data->finalRlabelMapping;
      if (const auto oldLabel = std::find_if(
              mapping.begin(), mapping.end(),
              [label = label](const auto &p) { return p.second == label; });
          oldLabel != mapping.end()) {
        if (auto iter = match.rgroups.find(oldLabel->first);
            iter != match.rgroups.end()) {
          MolOps::setTerminalAtomCoords(*coreWithMatches, atom->getIdx(),
                                        neighbor->getIdx());
        }
      }
    }
  }

  if (!coreWithMatches->getRingInfo()->isInitialized()) {
    MolOps::symmetrizeSSSR(*coreWithMatches);
  }
#ifdef VERBOSE
  std::cerr << "output core mol3 " << MolToSmiles(*coreWithMatches)
            << std::endl;
#endif
  return coreWithMatches;
}

RGroupRows RGroupDecomposition::getRGroupsAsRows() const {
  std::vector<RGroupMatch> permutation = data->GetCurrentBestPermutation();

  RGroupRows groups;

  auto usedLabelMap = UsedLabelMap(data->finalRlabelMapping);

  for (auto it = permutation.begin(); it != permutation.end(); ++it) {
    auto Rs_seen(usedLabelMap);
    // make a new rgroup entry
    groups.push_back(RGroupRow());
    RGroupRow &out_rgroups = groups.back();

    const R_DECOMP &in_rgroups = it->rgroups;

    for (const auto &rgroup : in_rgroups) {
      const auto realLabel = data->finalRlabelMapping.find(rgroup.first);
      CHECK_INVARIANT(realLabel != data->finalRlabelMapping.end(),
                      "unprocessed rlabel, please call process() first.");
      Rs_seen.setIsUsed(realLabel->second);
      out_rgroups[RPREFIX + std::to_string(realLabel->second)] =
          rgroup.second->combinedMol;
    }

    out_rgroups[CORE] = outputCoreMolecule(*it, Rs_seen);
  }
  return groups;
}

//! return rgroups in column order group[attachment_point][molidx] = ROMol
RGroupColumns RGroupDecomposition::getRGroupsAsColumns() const {
  std::vector<RGroupMatch> permutation = data->GetCurrentBestPermutation();

  RGroupColumns groups;
  std::unordered_set<std::string> rGroupWithRealMol{CORE};

  auto usedLabelMap = UsedLabelMap(data->finalRlabelMapping);

  unsigned int molidx = 0;
  for (auto it = permutation.begin(); it != permutation.end(); ++it, ++molidx) {
    auto Rs_seen(usedLabelMap);
    const R_DECOMP &in_rgroups = it->rgroups;

    for (const auto &rgroup : in_rgroups) {
      const auto realLabel = data->finalRlabelMapping.find(rgroup.first);
      CHECK_INVARIANT(realLabel != data->finalRlabelMapping.end(),
                      "unprocessed rlabel, please call process() first.");
      CHECK_INVARIANT(rgroup.second->combinedMol->hasProp(done),
                      "Not done! Call process()");

      CHECK_INVARIANT(!Rs_seen.getIsUsed(realLabel->second),
                      "R group label appears multiple times!");
      Rs_seen.setIsUsed(realLabel->second);
      std::string r = RPREFIX + std::to_string(realLabel->second);
      RGroupColumn &col = groups[r];
      if (molidx && col.size() < molidx - 1) {
        col.resize(molidx - 1);
      }
      col.push_back(rgroup.second->combinedMol);
      rGroupWithRealMol.insert(r);
    }
    groups[CORE].push_back(outputCoreMolecule(*it, Rs_seen));

    // add empty entries to columns where this molecule didn't appear
    for (const auto &realLabel : data->finalRlabelMapping) {
      if (!Rs_seen.getIsUsed(realLabel.second)) {
        std::string r = RPREFIX + std::to_string(realLabel.second);
        groups[r].push_back(boost::make_shared<RWMol>());
      }
    }
  }
  // purge R-group entries that have no mols
  for (auto it = groups.begin(); it != groups.end();) {
    auto itToErase = groups.end();
    if (!rGroupWithRealMol.count(it->first)) {
      itToErase = it;
    }
    ++it;
    if (itToErase != groups.end()) {
      groups.erase(itToErase);
    }
  }
  return groups;
}

const RGroupDecompositionParameters &RGroupDecomposition::params() const {
  return data->params;
}

namespace {
std::vector<unsigned int> Decomp(RGroupDecomposition &decomp,
                                 const std::vector<ROMOL_SPTR> &mols) {
  auto t0 = std::chrono::steady_clock::now();
  std::vector<unsigned int> unmatched;
  for (size_t i = 0; i < mols.size(); ++i) {
    int v = decomp.add(*mols[i].get());
    if (v == -1) {
      unmatched.push_back(i);
    }
    checkForTimeout(t0, decomp.params().timeout);
  }
  decomp.process();
  return unmatched;
}
}  // namespace
unsigned int RGroupDecompose(const std::vector<ROMOL_SPTR> &cores,
                             const std::vector<ROMOL_SPTR> &mols,
                             RGroupRows &rows,
                             std::vector<unsigned int> *unmatchedIndices,
                             const RGroupDecompositionParameters &options) {
  RGroupDecomposition decomp(cores, options);
  std::vector<unsigned int> unmatched = Decomp(decomp, mols);
  if (unmatchedIndices) {
    *unmatchedIndices = unmatched;
  }
  rows = decomp.getRGroupsAsRows();
  return mols.size() - unmatched.size();
}

unsigned int RGroupDecompose(const std::vector<ROMOL_SPTR> &cores,
                             const std::vector<ROMOL_SPTR> &mols,
                             RGroupColumns &columns,
                             std::vector<unsigned int> *unmatchedIndices,
                             const RGroupDecompositionParameters &options) {
  RGroupDecomposition decomp(cores, options);
  std::vector<unsigned int> unmatched = Decomp(decomp, mols);
  if (unmatchedIndices) {
    *unmatchedIndices = unmatched;
  }
  columns = decomp.getRGroupsAsColumns();
  return mols.size() - unmatched.size();
}
}  // namespace RDKit
