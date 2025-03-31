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

#include <GraphMol/SynthonSpaceSearch/SynthonSpaceHitSet.h>
#include <RDGeneral/export.h>
#include <DataStructs/ExplicitBitVect.h>

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

// Split the molecule into fragments.  maxNumFrags gives the maximum number
// of fragments to be produced in each set.  There will a vector of vectors of
// molecules.  Each inner vector contains the fragments from a split molecule.
// The maxNumFrags will be constrained to the maximum number of synthons in
// the search space as there's no point making more fragments than that.
// Any complex query atoms will be stripped out of the fragments and replaced
// by a simple atom query.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<std::vector<std::unique_ptr<ROMol>>>
splitMolecule(const ROMol &query, unsigned int maxNumFrags,
              const std::uint64_t maxNumFragSets, const TimePoint *endTime,
              const int numThreads, bool &timedOut);
// Counts the number of [1*], [2*]...[4*] in the string.
RDKIT_SYNTHONSPACESEARCH_EXPORT int countConnections(const ROMol &mol);

// Return a bitset for each fragment giving the connector patterns
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<boost::dynamic_bitset<>>
getConnectorPatterns(const std::vector<std::unique_ptr<ROMol>> &fragSet);

// Return a bitset giving the different connector types in this
// molecule.
RDKIT_SYNTHONSPACESEARCH_EXPORT boost::dynamic_bitset<> getConnectorPattern(
    const std::vector<std::unique_ptr<ROMol>> &fragSet);

// Gets the permutations of connector numbers and the atoms they should
// be applied to in the molFrags.
// E.g. if the reaction has 3 connectors, 1, 2 and 3 and the fragged mol has
// 2, return all permutations of 2 from 3.  It's ok if the fragged mol doesn't
// have all the connections in the reaction, although this may well result in
// a lot of hits.
RDKIT_SYNTHONSPACESEARCH_EXPORT
std::vector<std::vector<std::vector<std::pair<Atom *, unsigned int>>>>
getConnectorPermutations(const std::vector<std::unique_ptr<ROMol>> &molFrags,
                         const boost::dynamic_bitset<> &fragConns,
                         const boost::dynamic_bitset<> &reactionConns);

// As above, but just returns the bitsets for the connector permutations,
// not the molecules.
RDKIT_SYNTHONSPACESEARCH_EXPORT
std::vector<std::vector<boost::dynamic_bitset<>>> getConnectorPermutations(
    const std::vector<boost::dynamic_bitset<>> &fragConnPatts,
    const boost::dynamic_bitset<> &reactionConns);

// If all bits in one of the bitsets is unset, it means that nothing matched
// that synthon.  If at least one of the bitsets has a set bit, all products
// incorporating the synthon with no bits set must match the query so
// should be used because the query matches products that don't incorporate
// anything from 1 of the synthon lists. Therefore those bits will all be
// set on exit.  For example, if the synthons are
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

// Return a molecule containing the portions of the molecule starting at
// each dummy atom and going out up to 3 bonds.  There may be more than
// 1 fragment if there are dummy atoms more than 3 bonds apart, and there
// may be fragments with more than 1 dummy atom if their dummy atoms fall
// within 3 bonds of each other.  E.g. the molecule [1*]CN(C[2*])Cc1ccccc1
// will give [1*]CN(C)C[1*].  The 2 dummy atoms are 4 bonds apart, but the
// fragments overlap.  All dummy atoms given isotope 1 whatever they had
// before.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::unique_ptr<ROMol> buildConnRegion(
    const ROMol &mol);

// Take any query atoms out of the molecule, replacing them with the
// nearest thing possible.  Probably this will just be the atomic
// number.  It doesn't change dummy atoms that have an isotope as
// these will be connectors, or anything with an AtomType query as
// they are uncontroversial. Returns true if it did something,
// false if the molecule was left unchanged.
RDKIT_SYNTHONSPACESEARCH_EXPORT bool removeQueryAtoms(RWMol &mol);

// Put together a product name in the Enamine style, which uses
// a semicolon as a separator and has the reagents names followed
// by the reaction name.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::string buildProductName(
    const std::string &reactionId, const std::vector<std::string> &fragIds);
RDKIT_SYNTHONSPACESEARCH_EXPORT std::string buildProductName(
    const RDKit::SynthonSpaceSearch::SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &fragNums);
// Zip the fragments together to make a molecule.  Assumes the connection
// points are marking by isotope numbers on dummy atoms.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::unique_ptr<ROMol> buildProduct(
    const std::vector<const ROMol *> &synthons);

// Make a map that has all the fragments with the same SMILES
// in a vector keyed by that SMILES.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::map<std::string, std::vector<ROMol *>>
mapFragsBySmiles(std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets,
                 bool &cancelled);

}  // namespace SynthonSpaceSearch::details
}  // namespace RDKit

#endif  // RDKIT_SYNTHONSPACESEARCHDETAILS_H
