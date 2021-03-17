//  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FMCS/FMCS.h>
#include <boost/scoped_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <set>
#include <utility>
#include <vector>

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

int RGroupDecomposition::add(const ROMol &inmol) {
  // get the sidechains if possible
  //  Add hs for better symmetrization
  RWMol mol(inmol);
  bool explicitOnly = false;
  bool addCoords = true;
  MolOps::addHs(mol, explicitOnly, addCoords);

  int core_idx = 0;
  const RCore *rcore = nullptr;
  std::vector<MatchVectType> tmatches;

  // Find the first matching core.
  for (const auto &core : data->cores) {
    {
      const bool uniquify = false;
      const bool recursionPossible = true;
      const bool useChirality = true;
      SubstructMatch(mol, *core.second.core, tmatches, uniquify,
                     recursionPossible, useChirality);
    }

    if (data->params.onlyMatchAtRGroups) {
      std::vector<MatchVectType> tmatches_filtered;
      for (auto &mv : tmatches) {
        bool passes_filter = true;
        boost::dynamic_bitset<> target_match_indices(mol.getNumAtoms());
        for (auto &match : mv) {
          target_match_indices[match.second] = 1;
        }

        for (auto &match : mv) {
          const Atom *atm = mol.getAtomWithIdx(match.second);
          // is this a labelled rgroup or not?
          if (core.second.core_atoms_with_user_labels.find(match.first) ==
              core.second.core_atoms_with_user_labels.end()) {
            // nope... if any neighbor is not part of the substructure
            //  make sure we are a hydrogen, otherwise, skip the match
            for (const auto &nbri :
                 boost::make_iterator_range(mol.getAtomNeighbors(atm))) {
              const auto &nbr = mol[nbri];
              if (nbr->getAtomicNum() != 1 &&
                  !target_match_indices[nbr->getIdx()]) {
                passes_filter = false;
                break;
              }
            }
          }
          if (!passes_filter) {
            break;
          }
        }

        if (passes_filter) {
          tmatches_filtered.push_back(mv);
        }
      }
      tmatches = tmatches_filtered;
    }

    if (tmatches.empty()) {
      continue;
    } else {
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
      rcore = &core.second;
      core_idx = core.first;
      break;
    }
  }

  if (rcore == nullptr) {
    BOOST_LOG(rdWarningLog) << "No core matches" << std::endl;
    return -1;
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
  for (auto &tmatche : tmatches) {
    const bool replaceDummies = false;
    const bool labelByIndex = true;
    const bool requireDummyMatch = false;
    bool hasCoreDummies = false;
    auto coreCopy =
        rcore->replaceCoreAtomsWithMolMatches(hasCoreDummies, mol, tmatche);
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
        std::vector<int> attachments;
        boost::shared_ptr<ROMol> &newMol = fragments[i];
        newMol->setProp<int>("core", core_idx);
        newMol->setProp<int>("idx", data->matches.size());
        newMol->setProp<int>("frag_idx", i);

        for (auto at : newMol->atoms()) {
          unsigned int elno = at->getAtomicNum();
          if (elno == 0) {
            unsigned int index =
                at->getIsotope();  // this is the index into the core
            // it messes up when there are multiple ?
            int rlabel;
            auto coreAtom = rcore->core->getAtomWithIdx(index);
            coreAtomAnyMatched.insert(index);
            if (coreAtom->getPropIfPresent(RLABEL, rlabel)) {
              std::vector<int> rlabelsOnSideChain;
              at->getPropIfPresent(SIDECHAIN_RLABELS, rlabelsOnSideChain);
              rlabelsOnSideChain.push_back(rlabel);
              at->setProp(SIDECHAIN_RLABELS, rlabelsOnSideChain);

              data->labels.insert(rlabel);  // keep track of all labels used
              attachments.push_back(rlabel);
            }
          }
        }

        if (attachments.size() > 0) {
          // reject multiple attachments?
          // what to do with labelled cores ?
          std::string newCoreSmi = MolToSmiles(*newMol, true);

          for (size_t attach_idx = 0; attach_idx < attachments.size();
               ++attach_idx) {
            int rlabel = attachments[attach_idx];
            ADD_MATCH(match, rlabel);
            match[rlabel]->add(newMol, attachments);
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
                data->cores[core_idx] =
                    RCore(newCore, data->params.onlyMatchAtRGroups);
                return add(inmol);
              }
            }
          }
        }
      }

      if (match.size()) {
        auto numberUserGroupsInMatch = std::accumulate(
            match.begin(), match.end(), 0,
            [](int sum, std::pair<int, boost::shared_ptr<RGroupData>> p) {
              return p.first > 0 && !p.second->is_hydrogen ? ++sum : sum;
            });
        int numberMissingUserGroups =
            rcore->numberUserRGroups - numberUserGroupsInMatch;
        CHECK_INVARIANT(numberMissingUserGroups >= 0,
                        "Data error in missing user rgroup count");
        potentialMatches.emplace_back(
            core_idx, numberMissingUserGroups, match,
            hasCoreDummies || !data->params.onlyMatchAtRGroups ? coreCopy
                                                               : nullptr);
      }
    }
  }
  if (potentialMatches.size() == 0) {
    BOOST_LOG(rdWarningLog)
        << "No attachment points in side chains" << std::endl;
    return -1;
  }

  if (data->params.matchingStrategy != GA) {
    size_t N = 1;
    for (auto &matche : data->matches) {
      size_t sz = matche.size();
      N *= sz;
    }
    // oops, exponential is a pain
    if (N * potentialMatches.size() > 100000) {
      data->permutation = std::vector<size_t>(data->matches.size(), 0);
      data->process(true);
    }
  }

  data->matches.push_back(potentialMatches);
  data->permutation = std::vector<size_t>(data->matches.size(), 0);

  if (data->matches.size()) {
    if (data->params.matchingStrategy & Greedy ||
        (data->params.matchingStrategy & GreedyChunks &&
         data->matches.size() > 1 &&
         data->matches.size() % data->params.chunkSize == 0)) {
      data->process(true);
    }
  }
  return data->matches.size() - 1;
}

bool RGroupDecomposition::process() { return processAndScore().success; }

RGroupDecompositionProcessResult RGroupDecomposition::processAndScore() {
  try {
    const bool prune = true;
    const bool finalize = true;
    return data->process(prune, finalize);
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
  const auto &core = data->cores[match.core_idx];
  if (!match.matchedCore) {
    return core.labelledCore;
  }
  auto coreWithMatches = core.coreWithMatches(*match.matchedCore);
  for (auto atomIdx = coreWithMatches->getNumAtoms(); atomIdx--;) {
    auto atom = coreWithMatches->getAtomWithIdx(atomIdx);
    if (atom->getAtomicNum()) {
      continue;
    }
    auto label = atom->getAtomMapNum();
    Atom *nbrAtom = nullptr;
    for (const auto &nbri :
         boost::make_iterator_range(coreWithMatches->getAtomNeighbors(atom))) {
      nbrAtom = (*coreWithMatches)[nbri];
      break;
    }
    if (nbrAtom) {
      bool isUserDefinedLabel = usedLabelMap.isUserDefined(label);
      auto numExplicitHs = nbrAtom->getNumExplicitHs();
      if (usedLabelMap.getIsUsed(label)) {
        if (numExplicitHs) {
          nbrAtom->setNumExplicitHs(numExplicitHs - 1);
        }
      } else if (!isUserDefinedLabel ||
                 data->params.removeAllHydrogenRGroupsAndLabels) {
        coreWithMatches->removeAtom(atomIdx);
        // if we remove an unused label from an aromatic atom,
        // we need to check whether we need to adjust its explicit
        // H count, or it will fail to kekulize
        if (isUserDefinedLabel && nbrAtom->getIsAromatic()) {
          nbrAtom->updatePropertyCache(false);
          if (!numExplicitHs) {
            nbrAtom->setNumExplicitHs(nbrAtom->getExplicitValence() -
                                      nbrAtom->getDegree());
          }
        }
      }
      nbrAtom->updatePropertyCache(false);
    }
  }
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
