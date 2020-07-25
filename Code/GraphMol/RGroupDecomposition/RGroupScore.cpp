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
struct rgroup_column {
  unsigned int count;
  bool onlyH;
};
  
double score(const std::vector<size_t> &permutation,
             const std::vector<std::vector<RGroupMatch>> &matches,
             const std::set<int> &labels) {
  double score = 1.;

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
    std::map<std::string, rgroup_column> matchSet;
    std::map<std::set<int>, int> linkerMatchSet;

    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
      auto rg = matches[m][permutation[m]].rgroups.find(l);
      if (rg != matches[m][permutation[m]].rgroups.end()) {
#ifdef DEBUG
        std::cerr << "  combined: " << MolToSmiles(*rg->second->combinedMol)
                  << std::endl;
        std::cerr << " RGroup: " << rg->second->smiles << " "
                  << rg->second->is_hydrogen << std::endl;;
#endif
	rgroup_column &col = matchSet[rg->second->smiles];
	col.count += 1;
        // detect whether or not this is an H
        if (rg->second->is_hydrogen) {
	  col.onlyH = true;
        } else {
	  col.onlyH = false;
        }
#ifdef DEBUG
        std::cerr << " " << rg->second->combinedMol->getNumAtoms(false)
                  << " isH: " << onlyH[rg->second->smiles]
                  << " score: " << matchSet[rg->second->smiles] << std::endl;
#endif
        if (rg->second->is_linker) {
          linkerMatchSet[rg->second->attachments]++;
#ifdef DEBUG
          std::cerr << " Linker Score: "
                    << linkerMatchSet[rg->second->attachments]++ << std::endl;
#endif
        }
      }
    }
    
    // get the counts for each rgroup found and sort in reverse order
    std::vector<float> equivalentRGroupCount;

    for (const auto &kv :  matchSet) {
#ifdef DEBUG
      std::cerr << " equiv: " << it->first << " " << it->second << " "
                << permutation.size() << std::endl;
#endif

      // if the rgroup is hydrogens, only consider if the group is all
      //  hydrogen, otherwise score based on the non hydrogens
      if (kv.second.onlyH) {
        if (static_cast<size_t>(kv.second.count) == permutation.size()) {
          equivalentRGroupCount.push_back(static_cast<float>(kv.second.count));
        } else {
          // hydrogens in a mixed group don't contribute to the score
          equivalentRGroupCount.push_back(0.0);
        }
      } else {
        equivalentRGroupCount.push_back(static_cast<float>(kv.second.count));
      }
    }
    std::sort(equivalentRGroupCount.begin(), equivalentRGroupCount.end(),
              std::greater<float>());

    double tempScore = 1.;
    // score the sets from the largest to the smallest
    //  each smaller set gets penalized (i+1) below
    //  1.0 is the perfect score
    for (size_t i = 0; i < equivalentRGroupCount.size(); ++i) {
      auto lscore =
          equivalentRGroupCount[i] / ((i + 1) * (double)matches.size());
      if (lscore > 0) {
        tempScore *= lscore * lscore;
      }
#ifdef DEBUG
      std::cerr << "    lscore^2 " << i << ": " << lscore * lscore << std::endl;
#endif
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
      linkerIncrement = (double)(maxLinkerMatches) / (double)matches.size();
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
