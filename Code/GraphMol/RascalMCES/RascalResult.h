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

#ifndef RASCALRESULT_H
#define RASCALRESULT_H

#include <vector>

#include <GraphMol/ROMol.h>

namespace RDKit {

namespace RascalMCES {

class RascalResult {
 public:
  RascalResult(const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
               const std::vector<std::vector<int>> &adj_matrix1,
               const std::vector<std::vector<int>> &adj_matrix2,
               const std::vector<unsigned int> &clique,
               const std::vector<std::pair<int, int>> &vtx_pairs,
               bool timed_out, bool swapped, bool chiralSmarts, int minFragSep);

  RascalResult(const RascalResult &other);

  RascalResult(RascalResult &&other) = default;

  ~RascalResult() = default;

  RascalResult &operator=(const RascalResult &other);

  RascalResult &operator=(RascalResult &&other) = default;

  // Cut the result down to the single largest fragment.  This is
  // irrecoverably destructive.
  void largest_frag_only();

  std::shared_ptr<RDKit::ROMol> mol1() const { return d_mol1; };

  std::shared_ptr<RDKit::ROMol> mol2() const { return d_mol2; };

  std::vector<std::pair<int, int>> bond_matches() const {
    return d_bondMatches;
  }

  std::vector<std::pair<int, int>> atom_matches() const {
    return d_atomMatches;
  }

  int num_frags() const;

  int ring_non_ring_bond_score() const;

  int atom_match_score() const;

  int max_delta_atom_atom_dist() const;

  // returns number of atoms for largest fragment.
  int largest_frag_size() const;

  std::string smarts() const;

  bool timedout() const { return d_timedOut; };

  double similarity() const;

 private:
  std::shared_ptr<RDKit::ROMol> d_mol1;
  std::shared_ptr<RDKit::ROMol> d_mol2;
  std::vector<std::pair<int, int>> d_bondMatches;
  std::vector<std::pair<int, int>> d_atomMatches;

  mutable std::string d_smarts;
  bool d_timedOut{false};
  bool d_chiralSmarts{false};
  int d_maxFragSep{-1};

  // These are used for sorting the results.
  mutable int d_numFrags{-1};
  mutable int d_ringNonRingBondScore{-1};
  mutable int d_atomMatchScore{-1};
  mutable int d_maxDeltaAtomAtomDist{-1};
  mutable int d_largestFragSize{-1};

  std::string create_smarts_string() const;

  void match_clique_atoms(const std::vector<std::vector<int>> &mol1_adj_matrix);

  // If the clique involves a fragment that is more than d_maxFragSep from
  // any other frag in either molecule, discard the smaller frag.
  void apply_max_frag_sep();

  // Make the fragments for either mol1 or mol2.  If molNum is not 1 or 2,
  // returns nullptr.
  RDKit::ROMol *make_mol_frags(int molNum) const;

  int calc_ring_non_ring_score() const;

  int calc_atom_match_score() const;

  int calc_largest_frag_size() const;

  // If there are multiple fragments, can be helpful as a tie-breaker.  It's the
  // maximum difference between through-bond distances between matching atoms in
  // the 2 molecules.
  int calc_max_delta_atom_atom_dist_score() const;
};

bool resultSort(const RascalResult &res1, const RascalResult &res2);

void extract_clique(const std::vector<unsigned int> &clique,
                    const std::vector<std::pair<int, int>> &vtx_pairs,
                    bool swapped,
                    std::vector<std::pair<int, int>> &bond_matches);

// do some simple cleaning of the SMARTS, to make it more user-friendly.
void clean_smarts(std::string &smarts);

// Primarily for debugging, these write out the corresponding bonds/atoms
// in Python list format, for ease of cut/paste into a highlighted image
// creation.
void print_bond_matches(const RascalResult &res, std::ostream &os);

void print_atom_matches(const RascalResult &res, std::ostream &os);

double johnson_similarity(const std::vector<std::pair<int, int>> &bond_matches,
                          const RDKit::ROMol &mol1, const RDKit::ROMol &mol2);

}  // namespace RascalMCES
}  // namespace RDKit

#endif  // RASCALRESULT_H
