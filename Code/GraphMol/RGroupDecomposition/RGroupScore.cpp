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
#include <vector>
#include <map>
#include <algorithm>
namespace RDKit {

// stupid total score
// This has to handle all permutations and doesn't do anything terribly smart
//  For r-groups with large symmetries, this can take way too long.
double matchScore(const std::vector<size_t> &permutation,
                  const std::vector<std::vector<RGroupMatch>> &matches,
                  const std::set<int> &labels) {
  double score = 0.;
  const std::string EMPTY_RGROUP = "";
#ifdef DEBUG
  std::cerr << "---------------------------------------------------"
            << std::endl;
  std::cerr << "Scoring permutation "
            << " num matches: " << matches.size() << std::endl;


  BOOST_LOG(rdDebugLog) << "Scoring" << std::endl;
  for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
    BOOST_LOG(rdDebugLog)
      << "Molecule " << m << " " << matches[m][permutation[m]].toString()
      << std::endl;
  }
#endif

  // What is the largest rgroup count at any label
  int N = 0;
  std::map<int, int> num_rgroups;
  for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
    for (auto l : matches[m][permutation[m]].rgroups) {
       N = std::max(N, ++num_rgroups[l.first]);
    }
  }
  // for each label (r-group)
  for(auto l : labels ) {
#ifdef DEBUG
    std::cerr << "Label: " << l << std::endl;
#endif
    std::vector<std::map<std::string, unsigned int>> matchSetVect;
    std::map<std::set<int>, size_t> linkerMatchSet;

    int num_rgroups_for_label = 0;
    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule

      auto rg = matches[m][permutation[m]].rgroups.find(l);
      if (rg == matches[m][permutation[m]].rgroups.end()) {
        continue;
      }
      num_rgroups_for_label++;
      if (rg->second->is_linker) {
        ++linkerMatchSet[rg->second->attachments];
#ifdef DEBUG
        std::cerr << "  combined: " << MolToSmiles(*rg->second->combinedMol)
                  << std::endl;
        std::cerr << " RGroup: " << rg->second->smiles << " "
                  << rg->second->is_hydrogen << std::endl;
        ;
#endif
      }
#ifdef DEBUG
        std::cerr << l << " rgroup count" << num_rgroups_for_label << " num atoms" << rg->second->combinedMol->getNumAtoms(false)
                // looks like code has been edited round this define
                // << " score: " << count
                << std::endl;
#endif
      size_t i = 0;
      for (const auto &smiles : rg->second->smilesVect) {
        if (i == matchSetVect.size()) {
          matchSetVect.resize(i + 1);
        }
        unsigned int &count = matchSetVect[i][smiles];
        ++count;
#ifdef DEBUG
          std::cerr << i << " smiles:" << smiles << " " << count << std::endl;
        std::cerr << " Linker Score: "
                  << linkerMatchSet[rg->second->attachments] << std::endl;
#endif
        ++i;
      }
    }
    
    double tempScore = 0.;
    for (auto &matchSet : matchSetVect) {
      // get the counts for each rgroup found and sort in reverse order
      // If we don't have as many rgroups as the largest set add a empty ones
      if( N - num_rgroups_for_label > 0) {
          matchSet[EMPTY_RGROUP] = N - num_rgroups_for_label;
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
      tempScore /= matchSetVect.size();
    }

    // overweight linkers with the same attachments points....
    //  because these belong to 2 (or more) rgroups we really want these to stay
    //  the size of the set is the number of labels that are being used
    //  ** this heuristic really should be taken care of above **
    unsigned int maxLinkerMatches = 0;
    for (const auto &it : linkerMatchSet) {
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
  } // end for each label

#ifdef DEBUG
  BOOST_LOG(rdDebugLog) << score << std::endl;
#endif

  return score;
}

}  // namespace RDKit
