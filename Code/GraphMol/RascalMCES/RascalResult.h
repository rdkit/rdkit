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
               bool ringMatchesRingOnly, bool singleLargestFrag, int minFragSep,
               bool exactConnectionsMatch = false);
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

  // The following 5 functions are used in resultCompare to rank
  // 2 MCES of the same size for the same pair of molecules.
  // returns the number of contiguous fragments in the MCES.
  int getNumFrags() const;

  // returns how many bonds in the clique don't match
  // cyclic/non-cyclic i.e. count as a matche in the MCES but
  // are ring bonds in one of the molecules and not in the other.
  int getRingNonRingBondScore() const;

  // returns a score for how well the atoms in the clique from mol1 match the
  // atoms for the clique in mol2.  Currently, the atom scores are the
  // difference in H count for matching atoms, and summed for the molecule. Its
  // so that, for example, an OH in mol1 that could match an OH or OMe matches
  // the OH for preference.
  int getAtomMatchScore() const;

  // returns a score for the maximum difference in through-bond distance for
  // pairs of matching atoms in the 2 molecules.  An MCES where 2 atoms
  // are far apart in one molecule and the corresponding atoms are close
  // together in the other will get a high score by this measure.
  int getMaxDeltaAtomAtomDist() const;

  // returns the number of atoms in the largest contiguous fragment
  // in the MCES.
  unsigned int getLargestFragSize() const;

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
  bool d_exactConnectionsMatch{false};

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

}  // namespace RascalMCES
}  // namespace RDKit

#endif  // RASCALRESULT_H
