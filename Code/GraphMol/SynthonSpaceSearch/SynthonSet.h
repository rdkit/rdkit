//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_SYNTHONSET_H
#define RDKIT_SYNTHONSET_H

#include <iosfwd>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SynthonSpaceSearch/Synthon.h>

namespace RDKit {
class ROMol;

namespace SynthonSpaceSearch {
class Synthon;

// This class holds all the synthons for a particular reaction.
class RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSet {
 public:
  SynthonSet() = default;
  explicit SynthonSet(const std::string &id) : d_id(id) {}
  SynthonSet(const SynthonSet &rhs) = delete;
  SynthonSet(SynthonSet &&rhs) = delete;

  [[nodiscard]] const std::string &getId() const { return d_id; }
  [[nodiscard]] const std::vector<std::vector<std::unique_ptr<Synthon>>> &
  getSynthons() const {
    return d_synthons;
  }
  [[nodiscard]] const boost::dynamic_bitset<> &getConnectors() const {
    return d_connectors;
  }
  [[nodiscard]] const std::vector<std::shared_ptr<ROMol>> &getConnectorRegions()
      const;

  [[nodiscard]] const std::unique_ptr<ExplicitBitVect> &getConnRegFP() const;
  [[nodiscard]] const std::vector<int> &getNumConnectors() const;

  // Writes to/reads from a binary stream.
  void writeToDBStream(std::ostream &os) const;
  void readFromDBStream(std::istream &is, std::uint32_t version);

  // SynthonSet takes control of the newSynthon and manages it.
  void addSynthon(int synthonSetNum, std::unique_ptr<Synthon> newSynthon);

  // Build the connector regions and their fingerprints.  Only used when
  // creating a SynthonSpace from a text file.
  void buildConnectorRegions();

  // Scan through the connectors ([1*], [2*] etc.) in the synthons
  // and set bits in d_connectors accordingly.  Also removes any empty
  // reagent sets, which might be because the synthon numbers start from
  // 1 rather than 0.  Only used when creating a SynthonSpace from a text
  // file.
  void assignConnectorsUsed();

 private:
  std::string d_id;
  // The lists of synthons.  A product of the reaction is created by
  // combining 1 synthon from each of the outer vectors.
  std::vector<std::vector<std::unique_ptr<Synthon>>> d_synthons;
  // 4 bits showing which connectors are present in the synthons.
  boost::dynamic_bitset<> d_connectors;

  // The connector regions of a molecule are the pieces of up to 3 bonds from
  // a connector atom into the molecule.  We keep a vector of all the ones
  // present in the synthons in the set, plus a fingerprint of all their
  // fingerprints folded into 1.  If a query fragment doesn't have a
  // connector region in common with any of the synthons it can be assumed that
  // the fragment won't have a match in this SynthonSet.
  std::vector<std::shared_ptr<ROMol>> d_connectorRegions;
  // The fingerprint of the connector regions.  Fingerprints for all
  // connector regions are folded into the same fingerprint.
  std::unique_ptr<ExplicitBitVect> d_connRegFP;
  // The number of connectors in the synthons in each synthon set.
  std::vector<int> d_numConnectors;
};

}  // namespace SynthonSpaceSearch

}  // namespace RDKit

#endif  // RDKIT_SYNTHONSET_H
