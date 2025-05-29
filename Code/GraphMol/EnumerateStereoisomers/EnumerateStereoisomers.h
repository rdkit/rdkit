//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This is and the corresponding C++ file contain a transliteration
// of EnumerateStereoisomers.py.  It enumerates all tetrahedral,
// double bond and astropisomer stereochemistry.  The last being
// and enhancement over the Python original.

#ifndef RD_ENUMERATESTEREOISOMERS_H
#define RD_ENUMERATESTEREOISOMERS_H

#include <GraphMol/EnumerateStereoisomers/Flippers.h>

#include <random>
#include <unordered_set>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

namespace RDKit {
namespace EnumerateStereoisomers {
struct RDKIT_ENUMERATESTEREOISOMERS_EXPORT StereoEnumerationOptions {
  bool tryEmbedding{false};   // If true, the process attempts to generate
                              // a standard RDKit distance geometry
                              // conformation for the stereoisomer.  If this
                              // fails, we assume that the stereoisomer is
                              // non-physical and don't return it.  NOTE that
                              // this is computationally expensive and is just
                              // a heuristic that could result in
                              // stereoisomers being lost.
  bool onlyUnassigned{true};  // If true, stereocenters which have a
                              // specified stereochemistry will not be
                              // perturbed unless they are part of a relative
                              // stereo group.
  bool onlyStereoGroups{false};  // If true, only find stereoisomers that
                                 // differ at the StereoGroups associated
                                 // with the molecule.
  bool unique{true};  // If true, only stereoisomers that differ in canonical
                      // SMILES will be returned.
  std::uint64_t maxIsomers{0};  // The maximum number of isomers to yield.
                                // If the number of possible isomers is
                                // greater than maxIsomers, a random subset
                                // will be yielded.  If 0, all isomers are
                                // yielded.  Since every additional
                                // stereocenter doubles the number of results
                                // (and execution time) it's important to
                                // keep an eye on this.
  int randomSeed{-1};  // Seed for random number generator.  -1 means don't
                       // seed.
};

// Class that enumerates the stereoisomers of a molecule.  Acts like a
// Python generator so in principle has no limit on the number of stereoisomers
// it can produce.
class RDKIT_ENUMERATESTEREOISOMERS_EXPORT StereoisomerEnumerator {
 public:
  StereoisomerEnumerator() = delete;
  StereoisomerEnumerator(
      const ROMol &mol,
      const StereoEnumerationOptions &options = StereoEnumerationOptions(),
      bool verbose = false);
  StereoisomerEnumerator(const StereoisomerEnumerator &other) = delete;
  StereoisomerEnumerator(StereoisomerEnumerator &&other) = delete;
  ~StereoisomerEnumerator() = default;
  StereoisomerEnumerator &operator=(const StereoisomerEnumerator &other) =
      delete;
  StereoisomerEnumerator &operator=(StereoisomerEnumerator &&other) = delete;

  std::uint64_t getStereoisomerCount() const;

  // Return another stereoisomer, or an empty unique_ptr if we're done.
  std::unique_ptr<ROMol> next();

 private:
  RWMol d_mol;
  const StereoEnumerationOptions d_options;
  bool d_verbose;
  // The number of isomers successfully returned so far.  This may be
  // lower than d_seen.size() if d_options.tryEmbedding is true and there
  // have been embedding failures.
  std::uint64_t d_numReturned{0};
  // The number we need to return, the smaller of 2**N, where N is the
  // number of stereo centers, or d_options.maxIsomers.
  std::uint64_t d_numToReturn{1024};
  // 2**N.
  std::uint64_t d_totalPoss{0};
  std::unordered_set<std::string> d_generatedIsomers;

  // For the random bools
  std::unique_ptr<std::mt19937> d_randGen;
  std::bernoulli_distribution d_randDis{0.5};

  // Classes for setting the orientation at a particular stereocenter.
  std::vector<std::unique_ptr<details::Flipper>> d_flippers;

  // The stereo orientations we've already made
  std::unordered_set<boost::dynamic_bitset<>> d_seen;

  void buildFlippers();
  std::unique_ptr<ROMol> generateRandomIsomer();
  bool embeddable(ROMol &isomer);
};
}  // namespace EnumerateStereoisomers
}  // namespace RDKit

#endif  // RD_ENUMERATESTEREOISOMERS_H
