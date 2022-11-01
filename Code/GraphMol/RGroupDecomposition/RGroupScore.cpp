//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// #define DEBUG

#include "RGroupScore.h"
#include "RDGeneral/Invariant.h"
#include <vector>
#include <map>
#include <algorithm>
namespace RDKit {

// stupid total score
// This has to handle all permutations and doesn't do anything terribly smart
// For r-groups with large symmetries, this can take way too long.
double RGroupScorer::matchScore(
    const std::vector<size_t> &permutation,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels) {
  PRECONDITION(permutation.size() <= matches.size(),
               "permutation.size() should be <= matches.size()");
  double score = 0.;
  const std::string EMPTY_RGROUP = "";
  size_t offset = matches.size() - permutation.size();
#ifdef DEBUG
  std::cerr << "---------------------------------------------------"
            << std::endl;
  std::cerr << "Scoring permutation "
            << " num matches: " << matches.size() << std::endl;

  BOOST_LOG(rdDebugLog) << "Scoring" << std::endl;
  for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
    BOOST_LOG(rdDebugLog) << "Molecule " << m << " "
                          << matches[m + offset].at(permutation[m]).toString()
                          << std::endl;
  }
#endif

  // What is the largest rgroup count at any label
  restoreInitialState();
  std::map<int, int> num_rgroups;
  for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
    for (auto l : matches[m + offset].at(permutation[m]).rgroups) {
      d_current.N = std::max(d_current.N, ++num_rgroups[l.first]);
    }
  }
  // for each label (r-group)
  for (auto l : labels) {
#ifdef DEBUG
    std::cerr << "Label: " << l << std::endl;
#endif
    auto &labelData = d_current.labelDataMap[l];
    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule

      const auto &match = matches[m + offset][permutation[m]];
      auto rg = match.rgroups.find(l);
      if (rg == match.rgroups.end()) {
        continue;
      }
      ++labelData.numRGroups;
      if (rg->second->is_linker) {
        ++labelData.linkerMatchSet[rg->second->attachments];
#ifdef DEBUG
        std::cerr << "  combined: " << MolToSmiles(*rg->second->combinedMol)
                  << std::endl;
        std::cerr << " RGroup: " << rg->second->smiles << " "
                  << rg->second->is_hydrogen << std::endl;
        ;
#endif
      }
#ifdef DEBUG
      std::cerr << l << " rgroup count" << labelData.numRGroups << " num atoms"
                << rg->second->combinedMol->getNumAtoms(false)
                // looks like code has been edited round this define
                // << " score: " << count
                << std::endl;
#endif
      size_t i = 0;
      for (const auto &smiles : rg->second->smilesVect) {
        if (i == labelData.matchSetVect.size()) {
          labelData.matchSetVect.resize(i + 1);
        }
        unsigned int &count = labelData.matchSetVect[i][smiles];
        ++count;
#ifdef DEBUG
        std::cerr << i << " smiles:" << smiles << " " << count << std::endl;
        std::cerr << " Linker Score: "
                  << labelData.linkerMatchSet[rg->second->attachments]
                  << std::endl;
#endif
        ++i;
      }
    }

    double tempScore = 0.;
    for (auto &matchSet : labelData.matchSetVect) {
      // get the counts for each rgroup found and sort in reverse order
      // If we don't have as many rgroups as the largest set add a empty ones
      if (d_current.N - labelData.numRGroups > 0) {
        matchSet[EMPTY_RGROUP] = d_current.N - labelData.numRGroups;
      }
      std::vector<unsigned int> equivalentRGroupCount;

      std::transform(matchSet.begin(), matchSet.end(),
                     std::back_inserter(equivalentRGroupCount),
                     [](const std::pair<std::string, unsigned int> &p) {
                       return p.second;
                     });
      std::sort(equivalentRGroupCount.begin(), equivalentRGroupCount.end(),
                std::greater<unsigned int>());

      // score the sets from the largest to the smallest
      // each smaller set gets penalized (i+1) below
      for (size_t i = 0; i < equivalentRGroupCount.size(); ++i) {
        auto lscore = static_cast<double>(equivalentRGroupCount[i]) /
                      static_cast<double>(((i + 1) * matches.size()));
        tempScore += lscore * lscore;
#ifdef DEBUG
        std::cerr << "    lscore^2 " << i << ": " << lscore * lscore
                  << std::endl;
#endif
      }
      // make sure to rescale groups like [*:1].[*:1]C otherwise this will be
      // double counted
      // WE SHOULD PROBABLY REJECT THESE OUTRIGHT
      tempScore /= static_cast<double>(labelData.matchSetVect.size());
    }

    // overweight linkers with the same attachments points....
    //  because these belong to 2 (or more) rgroups we really want these to stay
    //  the size of the set is the number of labels that are being used
    //  ** this heuristic really should be taken care of above **
    unsigned int maxLinkerMatches = 0;
    for (const auto &it : labelData.linkerMatchSet) {
      if (it.first.size() > 1 || it.second > 1) {
        if (it.first.size() > maxLinkerMatches) {
          maxLinkerMatches = std::max(it.first.size(), it.second);
        }
      }
    }
#ifdef DEBUG
    std::cerr << "Max Linker Matches :" << maxLinkerMatches << std::endl;
#endif
    double increment = 1.0;        // no change in score
    double linkerIncrement = 1.0;  // no change in score
    if (maxLinkerMatches) {
      linkerIncrement = static_cast<double>(maxLinkerMatches) /
                        static_cast<double>(matches.size());
    } else {
      increment = tempScore;
    }
    score += increment * linkerIncrement;
#ifdef DEBUG
    std::cerr << "Increment: " << increment
              << " Linker_Increment: " << linkerIncrement << std::endl;
    std::cerr << "increment*linkerIncrement: " << increment * linkerIncrement
              << std::endl;
    std::cerr << "Score = " << score << std::endl;
#endif
  }  // end for each label

#ifdef DEBUG
  BOOST_LOG(rdDebugLog) << score << std::endl;
#endif

  return score;
}

RGroupScorer::RGroupScorer(const std::vector<std::vector<size_t>> &permutations,
                           double score) {
  d_bestScore = score;
  for (const auto &permutation : permutations) {
    pushTieToStore(permutation);
  }
  if (!d_store.empty()) {
    d_saved = d_store.front();
  }
}

void RGroupScorer::setBestPermutation(const std::vector<size_t> &permutation,
                                      double score) {
  d_bestScore = score;
  d_current.permutation = permutation;
  d_saved = d_current;
}

void RGroupScorer::startProcessing() {
  d_initial = d_saved;
  d_bestScore = -std::numeric_limits<double>::max();
  clearTieStore();
}

void RGroupScorer::breakTies(
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels,
    const std::unique_ptr<CartesianProduct> &iterator,
    const std::chrono::steady_clock::time_point &t0, double timeout) {
  size_t maxPermValue = 0;
  d_current = d_saved;
  d_current.numAddedRGroups = labels.size();
  std::vector<int> largestHeavyCounts;
  largestHeavyCounts.reserve(labels.size());
  std::vector<int> orderedLabels;
  orderedLabels.reserve(labels.size());
  std::copy_if(labels.begin(), labels.end(), std::back_inserter(orderedLabels),
               [](const int &i) { return !(i < 0); });
  std::copy_if(labels.rbegin(), labels.rend(),
               std::back_inserter(orderedLabels),
               [](const int &i) { return (i < 0); });
  // We only care about the sign of the ordered labels,
  // not about their value, so we convert the ordered map
  // into a vector for comparing with the tied permutations
  // If there is a change in sign, then it means a new
  // label was added compared to the cached version,
  // so we need to add a new counter initialized to 0
  auto it = d_current.heavyCountPerLabel.begin();
  for (auto label : orderedLabels) {
    int count = 0;
    if (it != d_current.heavyCountPerLabel.end()) {
      if (!((it->first > 0) ^ (label > 0))) {
        count = it->second;
      }
      ++it;
    }
    largestHeavyCounts.push_back(count);
  }
  std::vector<int> initialHeavyCounts(largestHeavyCounts);
  while (!d_store.empty()) {
    auto &state = d_store.front();
    std::vector<int> heavyCounts(initialHeavyCounts);
    state.computeTieBreakingCriteria(matches, orderedLabels, heavyCounts);
#ifdef DEBUG
    std::cerr << "tiedPermutation ";
    for (const auto &t : state.permutation) {
      std::cerr << t << ",";
    }
    std::cerr << " orderedLabels ";
    for (const auto &l : orderedLabels) {
      std::cerr << l << ",";
    }
    std::cerr << " heavyCounts ";
    for (auto hc : heavyCounts) {
      std::cerr << hc << ",";
    }
    std::cerr << " largestHeavyCounts ";
    for (auto hc : largestHeavyCounts) {
      std::cerr << hc << ",";
    }
    std::cerr << " state.numMatchedUserRGroups " << state.numMatchedUserRGroups
              << " d_current.numMatchedUserRGroups "
              << d_current.numMatchedUserRGroups << ", state.numAddedRGroups "
              << state.numAddedRGroups << ", d_current.numAddedRGroups "
              << d_current.numAddedRGroups << std::endl;
#endif
    size_t permValue =
        iterator ? iterator->value(state.permutation) : maxPermValue;
    if (state.numMatchedUserRGroups > d_current.numMatchedUserRGroups) {
      d_current = state;
      largestHeavyCounts = heavyCounts;
      maxPermValue = permValue;
    } else if (state.numMatchedUserRGroups == d_current.numMatchedUserRGroups) {
      if (state.numAddedRGroups < d_current.numAddedRGroups) {
        d_current = state;
        largestHeavyCounts = heavyCounts;
        maxPermValue = permValue;
      } else if (state.numAddedRGroups == d_current.numAddedRGroups) {
        if (heavyCounts > largestHeavyCounts) {
          d_current = state;
          largestHeavyCounts = heavyCounts;
          maxPermValue = permValue;
        } else if (heavyCounts == largestHeavyCounts) {
          if (permValue > maxPermValue) {
            d_current = state;
            largestHeavyCounts = heavyCounts;
            maxPermValue = permValue;
          }
        }
      }
    }
    checkForTimeout(t0, timeout);
    d_store.pop_front();
  }
  // convert back the heavy count vector into an ordered map
  // to store it in the saved cache
  d_current.heavyCountPerLabel.clear();
  auto count = largestHeavyCounts.begin();
  for (auto label : orderedLabels) {
    d_current.heavyCountPerLabel[label] = *count++;
  }
  d_saved = d_current;
}

void RGroupScorer::pushTieToStore(const std::vector<size_t> &permutation) {
  d_current.permutation = permutation;
  d_store.push_back(d_current);
}

void RGroupScorer::clearTieStore() { d_store.clear(); }

}  // namespace RDKit
