#include "RGroupScore.h"

namespace RDKit {
// stupid total score
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

  for (int l : labels) {
#ifdef DEBUG
    std::cerr << "Label: " << l << std::endl;
#endif
    std::map<std::string, int> matchSet;
    std::map<std::set<int>, int> linkerMatchSet;
    std::map<std::string, int> onlyH;

    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
      auto rg = matches[m][permutation[m]].rgroups.find(l);
      if (rg != matches[m][permutation[m]].rgroups.end()) {
#ifdef DEBUG
        std::cerr << "  combined: " << MolToSmiles(*rg->second->combinedMol)
                  << std::endl;
        std::cerr << " RGroup: " << rg->second->smiles << " "
                  << rg->second->smiles.find_first_not_of("0123456789[]*H:.");
#endif
        matchSet[rg->second->smiles] += 1;
        // detect whether or not this is an H
        if (rg->second->smiles.find_first_not_of("0123456789[]*H:.") <
            rg->second->smiles.length()) {
          onlyH[rg->second->smiles] = 0;
        } else {
          onlyH[rg->second->smiles] = 1;
        }
#ifdef DEBUG
        std::cerr << " " << rg->second->combinedMol->getNumAtoms(false)
                  << " isH: " << onlyH[rg->second->smiles]
                  << " score: " << matchSet[rg->second->smiles] << std::endl;
#endif
        // XXX Use fragment counts to see if we are linking cycles?
        if (rg->second->smiles.find(".") == std::string::npos &&
            rg->second->attachments.size() > 1) {
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

    for (std::map<std::string, int>::const_iterator it = matchSet.begin();
         it != matchSet.end(); ++it) {
#ifdef DEBUG
      std::cerr << " equiv: " << it->first << " " << it->second << " "
                << permutation.size() << std::endl;
#endif

      // if the rgroup is hydrogens, only consider if the group is all
      //  hydrogen, otherwise score based on the non hydrogens
      if (onlyH[it->first]) {
        if (static_cast<size_t>(it->second) == permutation.size()) {
          equivalentRGroupCount.push_back(static_cast<float>(it->second));
        } else {
          // hydrogens in a mixed group don't contribute to the score
          equivalentRGroupCount.push_back(0.0);
        }
      } else {
        equivalentRGroupCount.push_back(static_cast<float>(it->second));
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
    for (std::map<std::set<int>, int>::const_iterator it =
             linkerMatchSet.begin();
         it != linkerMatchSet.end(); ++it) {
      if (it->second > 1) {
        if (it->second > maxLinkerMatches) {
          maxLinkerMatches = it->second;
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
