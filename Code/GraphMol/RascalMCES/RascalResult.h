//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

// A class to hold the results of a RASCAL MCES determination
// between 2 molecules.  Contains the bonds and atoms that
// correspond between the molecules, and also a SMARTS pattern
// defining the MCES.
//
#include <RDGeneral/export.h>

#ifndef RASCALRESULT_H
#define RASCALRESULT_H

#include <vector>

#include <GraphMol/ROMol.h>

namespace RDKit {

namespace RascalMCES {

class RDKIT_RASCALMCES_EXPORT RascalResult {
 public:
  RascalResult(const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
               const std::vector<std::vector<int>> &adjMatrix1,
               const std::vector<std::vector<int>> &adjMatrix2,
               const std::vector<unsigned int> &clique,
               const std::vector<std::pair<int, int>> &vtx_pairs, bool timedOut,
               bool swapped, double tier1Sim, double tier2Sim,
               bool ringMatchesRingOnly, bool singleLargestFrag,
               int minFragSep);
  // For when the tier[12]Sim didn't hit the threshold, but it
  // might be of interest what the estimates of similarity were.
  RascalResult(double tier1Sim, double tier2Sim);

  RascalResult(const RascalResult &other);

  RascalResult(RascalResult &&other) = default;

  ~RascalResult() = default;

  RascalResult &operator=(const RascalResult &other);

  RascalResult &operator=(RascalResult &&other) = default;

  // Cut the result down to the single largest fragment.  This is
  // irrecoverably destructive.
  void largestFragOnly();
  void largestFragsOnly(unsigned int numFrags = 2);
  void trimSmallFrags(unsigned int minFragSize = 3);

  std::vector<std::pair<int, int>> getBondMatches() const {
    return d_bondMatches;
  }

  std::vector<std::pair<int, int>> getAtomMatches() const {
    return d_atomMatches;
  }

  int getNumFrags() const;

  int getRingNonRingBondScore() const;

  int getAtomMatchScore() const;

  int getMaxDeltaAtomAtomDist() const;

  // returns number of atoms for largest fragment.
  int getLargestFragSize() const;

  std::string getSmarts() const;
  const std::shared_ptr<ROMol> getMcesMol() const;
  bool getTimedOut() const { return d_timedOut; };

  double getTier1Sim() const { return d_tier1Sim; }
  double getTier2Sim() const { return d_tier2Sim; }
  double getSimilarity() const;

 private:
  std::shared_ptr<ROMol> d_mol1;
  std::shared_ptr<ROMol> d_mol2;
  mutable std::shared_ptr<ROMol> d_mcesMol;
  std::vector<std::pair<int, int>> d_bondMatches;
  std::vector<std::pair<int, int>> d_atomMatches;

  mutable std::string d_smarts;
  bool d_timedOut{false};
  double d_tier1Sim;
  double d_tier2Sim;
  bool d_ringMatchesRingOnly{false};
  int d_maxFragSep{-1};

  // These are used for sorting the results.
  mutable int d_numFrags{-1};
  mutable int d_ringNonRingBondScore{-1};
  mutable int d_atomMatchScore{-1};
  mutable int d_maxDeltaAtomAtomDist{-1};
  mutable int d_largestFragSize{-1};

  // Assuming the frags are all part of the original MCES, just cut it
  // down to what's in the frags.
  void rebuildFromFrags(const std::vector<boost::shared_ptr<ROMol>> &frags);

  std::string createSmartsString() const;

  void matchCliqueAtoms(const std::vector<std::vector<int>> &mol1_adj_matrix);

  // If the clique involves a fragment that is more than d_maxFragSep from
  // any other frag in either molecule, discard the smaller frag.
  void applyMaxFragSep();

  // Make the fragments for either mol1 or mol2.  If molNum is not 1 or 2,
  // returns nullptr.
  RDKit::ROMol *makeMolFrags(int molNum) const;

  int calcRingNonRingScore() const;

  int calcAtomMatchScore() const;

  int calcLargestFragSize() const;

  // If there are multiple fragments, can be helpful as a tie-breaker.  It's the
  // maximum difference between through-bond distances between matching atoms in
  // the 2 molecules.
  int calcMaxDeltaAtomAtomDistScore() const;
};

RDKIT_RASCALMCES_EXPORT bool resultSort(const RascalResult &res1,
                                        const RascalResult &res2);

RDKIT_RASCALMCES_EXPORT void extractClique(
    const std::vector<unsigned int> &clique,
    const std::vector<std::pair<int, int>> &vtxPairs, bool swapped,
    std::vector<std::pair<int, int>> &bondMatches);

// do some simple cleaning of the SMARTS, to make it more user-friendly.
RDKIT_RASCALMCES_EXPORT void cleanSmarts(std::string &smarts);

// Primarily for debugging, these write out the corresponding bonds/atoms
// in Python list format, for ease of cut/paste into a highlighted image
// creation.
RDKIT_RASCALMCES_EXPORT void printBondMatches(const RascalResult &res,
                                              std::ostream &os);

RDKIT_RASCALMCES_EXPORT void printAtomMatches(const RascalResult &res,
                                              std::ostream &os);

// This prints out the scores in the order they are used in resultSort.
RDKIT_RASCALMCES_EXPORT void printScores(const RascalResult &res,
                                         std::ostream &os);

// Calculate the Johnson similarity between the two molecules using the given
// bondMatches.  It's the fraction of the 2 molecules that are in common,
// somewhat akin to the tanimoto - the square of the number of atoms plus
// number of bonds in the MCES divided by the product of the sums of the number
// of atoms and bonds in the 2 molecules.
// It has nothing to do with lying UK politicians.
RDKIT_RASCALMCES_EXPORT double johnsonSimilarity(
    const std::vector<std::pair<int, int>> &bondMatches,
    const std::vector<std::pair<int, int>> &atomMatches,
    const RDKit::ROMol &mol1, const RDKit::ROMol &mol2);

}  // namespace RascalMCES
}  // namespace RDKit

#endif  // RASCALRESULT_H
