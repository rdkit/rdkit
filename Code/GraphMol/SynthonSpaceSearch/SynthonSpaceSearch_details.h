//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_SYNTHONSPACESEARCHDETAILS_H
#define RDKIT_SYNTHONSPACESEARCHDETAILS_H

#include <chrono>
#include <vector>

#include <RDGeneral/export.h>

using Clock = std::chrono::steady_clock;
using TimePoint = std::chrono::time_point<Clock>;

namespace RDKit {
class ROMol;
namespace SynthonSpaceSearch::details {

RDKIT_SYNTHONSPACESEARCH_EXPORT bool checkTimeOut(const TimePoint *endTime);

// Find all combinations of M things selected from N.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<std::vector<unsigned int>>
combMFromN(unsigned int m, unsigned int n);
// Find all permutations of M things selected from N.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<std::vector<unsigned int>>
permMFromN(unsigned int m, unsigned int n);
RDKIT_SYNTHONSPACESEARCH_EXPORT void fixAromaticRingSplits(
    std::vector<std::unique_ptr<ROMol>> &molFrags);

// Split the molecule into fragments.  maxBondSplits gives the maximum number
// of bonds to be used in each split.  There will a vector of vectors of
// molecules, 1 inner vector for each split i.e. maxBondSplits in total, the
// first with 1 split, the 2nd with 2 etc.  Each inner vector contains the
// fragments from a split molecule.  The maxBondSplits will be constrained to
// between 1 and 4 inclusive, so if it is supplied outside that range, it will
// be altered.  Also, you can't split a molecule on 3 bonds if it only contains
// 2.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<std::vector<std::unique_ptr<ROMol>>>
splitMolecule(const ROMol &query, unsigned int maxBondSplits,
              std::uint64_t maxNumFrags, TimePoint *endTime, bool &timedOut);
// Counts the number of [1*], [2*]...[4*] in the string.
RDKIT_SYNTHONSPACESEARCH_EXPORT int countConnections(const ROMol &frag);

// Return a bitset for each fragment giving the connector patterns
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<boost::dynamic_bitset<>>
getConnectorPatterns(const std::vector<std::unique_ptr<ROMol>> &fragSet);

// Return a bitset giving the different connector types in this
// molecule.
RDKIT_SYNTHONSPACESEARCH_EXPORT boost::dynamic_bitset<> getConnectorPattern(
    const std::vector<std::unique_ptr<ROMol>> &fragSet);

// Return copies of the mol fragments will all permutations of the connectors
// in the reaction onto the connectors in the fragments.
// E.g. if the reaction has 3 connectors, 1, 2 and 3 and the fragged mol has
// 2, return all permutations of 2 from 3.  It's ok if the fragged mol doesn't
// have all the connections in the reaction, although this may well result in
// a lot of hits.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<std::vector<std::unique_ptr<ROMol>>>
getConnectorPermutations(const std::vector<std::unique_ptr<ROMol>> &molFrags,
                         const boost::dynamic_bitset<> &fragConns,
                         const boost::dynamic_bitset<> &reactionConns);

// If all bits in one of the bitsets is unset, it means that nothing matched
// that synthon.  If at least one of the bitsets has a set bit, all products
// incorporating the synthon with no bits set must match the query so
// should be used because the query matches products that don't incorporate
// anything from 1 of the synthon lists.  For example, if the synthons are
// [1*]Nc1c([2*])cccc1 and [1*]=CC=C[2*] and the query is c1ccccc1.
RDKIT_SYNTHONSPACESEARCH_EXPORT void expandBitSet(
    std::vector<boost::dynamic_bitset<>> &bitSets);

RDKIT_SYNTHONSPACESEARCH_EXPORT void bitSetsToVectors(
    const std::vector<boost::dynamic_bitset<>> &bitSets,
    std::vector<std::vector<size_t>> &outVecs);

// class to step through all combinations of lists of different sizes.
// returns (0,0,0), (0,0,1), (0,1,0) etc.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT Stepper {
  explicit Stepper(const std::vector<size_t> &sizes) : d_sizes(sizes) {
    d_currState = std::vector<size_t>(sizes.size(), 0);
  }
  void step() {
    // Don't do anything if we're at the end, but expect an infinite
    // loop if the user isn't wise to this.
    if (d_currState[0] == d_sizes[0]) {
      return;
    }
    std::int64_t i = static_cast<std::int64_t>(d_currState.size()) - 1;
    while (i >= 0) {
      ++d_currState[i];
      if (d_currState[0] == d_sizes[0]) {
        return;
      }
      if (d_currState[i] == d_sizes[i]) {
        d_currState[i] = 0;
      } else {
        break;
      }
      --i;
    }
  }
  std::vector<size_t> d_currState;
  std::vector<size_t> d_sizes;
};

}  // namespace SynthonSpaceSearch::details
}  // namespace RDKit

#endif  // RDKIT_SYNTHONSPACESEARCHDETAILS_H
