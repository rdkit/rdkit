//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RGroupScore.h"

namespace RDKit {
  
// stupid total score
// This has to handle all permutations and doesn't do anything terribly smart
//  For r-groups with large symmetries, this can take way too long.
double score(const std::vector<size_t> &permutation,
             const std::vector<std::vector<RGroupMatch>> &matches,
             const std::set<int> &labels) {
  double score = 0.;

#ifdef DEBUG
  std::cerr << "---------------------------------------------------"
            << std::endl;
  std::cerr << "Scoring permutation "
            << " num matches: " << matches.size() << std::endl;
#endif

  // For each label (group)
  for (int l : labels) {
#ifdef DEBUG
    std::cerr << "Label: " << l << std::endl;
#endif
    std::vector<std::map<std::string, unsigned int>> matchSetVect;
    std::map<std::set<int>, int> linkerMatchSet;

    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
      auto rg = matches[m][permutation[m]].rgroups.find(l);
      if (rg == matches[m][permutation[m]].rgroups.end()) {
        continue;
      }
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
      std::cerr << " " << rg->second->combinedMol->getNumAtoms(false)
                << " score: " << count << std::endl;
#endif
      size_t i = 0;
      for (const auto &smiles : rg->second->smilesVect) {
        if (i == matchSetVect.size()) {
          matchSetVect.resize(i + 1);
        }
        unsigned int &count = matchSetVect[i][smiles];
        ++count;
#ifdef DEBUG
        std::cerr << " Linker Score: "
                  << linkerMatchSet[rg->second->attachments] << std::endl;
#endif
        ++i;
      }
    }
    
    double tempScore = 0.;
    for (const auto &matchSet : matchSetVect) {
      // get the counts for each rgroup found and sort in reverse order
      std::vector<unsigned int> equivalentRGroupCount;

      std::transform(
          matchSet.begin(), matchSet.end(),
          std::back_inserter(equivalentRGroupCount),
          [](const std::pair<std::string, unsigned int> &p) { return p.second; });
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
    }

    // overweight linkers with the same attachments points....
    //  because these belong to 2 rgroups we really want these to stay
    //  ** this heuristic really should be taken care of above **
    int maxLinkerMatches = 0;
    for (const auto &it : linkerMatchSet ) {
      if (it.second > 1) {
        if (it.second > maxLinkerMatches) {
          maxLinkerMatches = it.second;
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
  }

  return score;
}  
}
