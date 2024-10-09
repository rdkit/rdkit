//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_HYPERSPACE_H
#define RDKIT_HYPERSPACE_H

#include <map>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/HyperspaceSearch/ReactionSet.h>

namespace RDKit {
class ROMol;

namespace HyperspaceSearch {

struct HyperspaceSearchParams {
  int maxBondSplits{3};  // The maximum number of bonds to break in the query.
                         // It should be no more than 1 less than the maximum
                         // number of Synthon sets in Hyperspace.  More than
                         // that doesn't matter, but will slow the search down
                         // to no good effect.
  int maxHits{1000};     // The maximum number of hits to return.  Use -1 for
                         // no maximum.
  bool buildHits{true};  // If false, reports the maximum number of hits that
                         // the search could produce, but doesn't return them.
};

// Holds the information about a set of hits.  The molecules can be built
// by making all combinations of reagents, one taken from each reagent set.
struct HyperspaceHitSet {
  std::string reactionId;
  std::vector<boost::dynamic_bitset<>> reagsToUse;
  size_t numHits{0};
};

class Hyperspace {
 public:
  // Create the hyperspace from a file in the correct format.
  explicit Hyperspace() = default;

  int numReactions() const { return d_reactions.size(); }
  const std::map<std::string, std::unique_ptr<ReactionSet>> &reactions() const {
    return d_reactions;
  }

  // Do a substructure search for query in the hyperspace.  Return vector of
  // molecules that match.
  std::vector<std::unique_ptr<ROMol>> substructureSearch(
      const ROMol &query,
      HyperspaceSearchParams params = HyperspaceSearchParams());

  // Search this particular fragmented molecule against the reactions.  The
  // fragments should be from 1 splitting, so between 1 and 4 members.
  // The fragments may be re-ordered in the process (largest fragment
  // heuristic).  This is in the public interface primarily for testing/
  // debugging purposes.
  std::vector<HyperspaceHitSet> searchFragSet(
      std::vector<std::unique_ptr<ROMol>> &fragSet);

  /*!
   *
   * @param inFile: name of the file containing the synthon-based library.
   *
   * The original format is:
   * all lines are tab-separated
   * first line:SMILES	synton_id	synton#	reaction_id
   * Note the spelling "synton" from the original paper/example file.
   * Subsequent lines have a single reagent e.g.
   * OCC([U])=NN=[Np]		1-1	0	triazole-1
   * C1CCCC1N([Pu])[U]		2-1	1	triazole-1
   * CC1CCN(C1)C(=[Np])[Pu]		3-1	2	triazole-1

   * Attachment points are U, Np, Pu and Am for up to 4 synthons per reaction.
   * A product is created by taking a synthon from each synton# value and
   * combining by replacing matching trans-uranic elements and replacing them
   * with a direct bond of the appropriate type.
   * A more (for RDKit) conventional connection flag of isotope labelled
   * dummy atoms is also accepted ([1*] etc.).
   * Throws a std::runtime_error if it doesn't think the format is correct,
   * which it does by checking that the first line is as above and subsequent
   * lines have 4 fields.
   * The formatting has been relaxed such that any whitespace may be used as
   * the field separator.
   */
  void readTextFile(const std::string &inFile);

  // Writes to/reads from a binary DB File in our format.
  void writeDBFile(const std::string &outFile) const;
  void readDBFile(const std::string &inFile);

  // Write a summary of the Hyperspace to given stream.
  void summarise(std::ostream &os) const;

 private:
  std::string d_fileName;
  std::map<std::string, std::unique_ptr<ReactionSet>> d_reactions;

  // Build the molecules from the reagents identified in reagentsToUse.
  // There should be bitset in reagentsToUse for each reagent set.
  // If not, it will fail.  Checks that all the results produced match the
  // query - it shouts if any don't, and doesn't include them.
  // If maxHits is -1, there's no limit on the number of hits returned.
  void buildHits(const std::vector<HyperspaceHitSet> &hitsets,
                 const ROMol &query,
                 std::vector<std::unique_ptr<ROMol>> &results, int maxHits);
  // get the subset of reagents for the given reaction to use for this
  // enumeration.
  std::vector<std::vector<ROMol *>> getReagentsToUse(
      const std::vector<boost::dynamic_bitset<>> &reagentsToUse,
      const std::string &reaction_id) const;
};

}  // namespace HyperspaceSearch
}  // namespace RDKit

#endif  // RDKIT_HYPERSPACE_H
