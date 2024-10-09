//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_REACTIONSET_H
#define RDKIT_REACTIONSET_H

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/HyperspaceSearch/Reagent.h>

namespace RDKit {
class ROMol;

namespace HyperspaceSearch {
class Reagent;

class ReactionSet {
 public:
  ReactionSet() = default;
  ReactionSet(const std::string &id) : d_id(id) {}
  ReactionSet(const Reagent &rhs) = delete;
  ReactionSet(Reagent &&rhs) = delete;

  const std::string &id() const { return d_id; }
  const std::vector<std::vector<std::unique_ptr<Reagent>>> &reagents() const {
    return d_reagents;
  }
  const boost::dynamic_bitset<> &connectors() const { return d_connectors; }
  const std::vector<std::shared_ptr<ROMol>> &connectorRegions() const;

  const std::unique_ptr<ExplicitBitVect> &connRegFP() const;
  const std::vector<int> &numConnectors() const;

  // Writes to/reads from a binary stream.
  void writeToDBStream(std::ostream &os) const;
  void readFromDBStream(std::istream &is);

  void addReagent(int reagentSetNum, const std::string &smiles,
                  const std::string &reagentId);

  // Scan through the connectors ([1*], [2*] etc.) in the reagents
  // and set bits in d_connectors accordingly.  Also removes any empty
  // reagent sets, which might be because the synthon numbers start from
  // 1 rather than 0.
  void assignConnectorsUsed();

 private:
  std::string d_id;
  std::vector<std::vector<std::unique_ptr<Reagent>>> d_reagents;
  // 4 bits showing which connectors are present in reagents
  boost::dynamic_bitset<> d_connectors;

  // The connector regions of a molecule are the pieces of up to 3 bonds from
  // a connector atom into the molecule.  We keep a vector of all the ones
  // present in the reagents in the set, plus a fingerprint of all their
  // fingerprints folded into 1.  If a query fragment doesn't have a
  // connector region in common with any of the reagents it can be assumed that
  // the fragment won't have a match in this ReagentSet.
  mutable std::vector<std::shared_ptr<ROMol>> d_connectorRegions;
  mutable std::unique_ptr<ExplicitBitVect> d_connRegFP;
  // The number of connectors in the synthons in each reagent set.
  mutable std::vector<int> d_numConnectors;
};

}  // namespace HyperspaceSearch

}  // namespace RDKit

#endif  // RDKIT_REACTIONSET_H
