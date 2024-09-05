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

#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

namespace RDKit {
class ROMol;

namespace HyperspaceSSSearch {

struct Reagent {
  std::string d_smiles;
  std::string d_id;
  std::unique_ptr<ROMol> d_mol;
  explicit Reagent(const std::string &smi, const std::string &id)
      : d_smiles(smi), d_id(id) {}
};

struct ReactionSet {
  std::string d_id;
  std::vector<std::vector<std::unique_ptr<Reagent>>> d_reagents;
  boost::dynamic_bitset<> d_connectors;
};

class Hyperspace {
 public:
  // Create the hyperspace from a file in the correct format.
  /*!
   *
   * @param fileName: name of the file containing the synthon-based library.
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
   *
   * Throws a std::runtime_error if it doesn't think the format is correct,
   * which is does by checking that the first line is as above and subsequent
   * lines have 4 fields.
   * The formatting has been relaxed such that any whitespace may be used as
   * the field separator.
   */
  explicit Hyperspace() = default;
  explicit Hyperspace(const std::string &fileName);

  int numReactions() const { return d_reactions.size(); }
  const std::vector<std::unique_ptr<ReactionSet>> &reactions() const {
    return d_reactions;
  }

  std::vector<std::unique_ptr<ROMol>> search(const ROMol &query,
                                             unsigned int maxBondSplits);

 private:
  std::string d_fileName;
  std::vector<std::unique_ptr<ReactionSet>> d_reactions;

  void readFile();
  // scan through the connectors ([1*], [2*] etc.) in the reagents in reach
  // ReagentSet and set bits in d_connectors accordingly.
  void assignConnectorsUsed();

  // Search this particular fragmented molecule against the reactions.  The
  // fragments should be from 1 splitting, so between 2 and 4 members.
  std::vector<std::unique_ptr<ROMol>> searchFragSet(
      std::unique_ptr<ROMol> &fraggedMol);
};

}  // namespace HyperspaceSSSearch
}  // namespace RDKit

#endif  // RDKIT_HYPERSPACE_H
