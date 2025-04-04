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

#include <RDGeneral/export.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SynthonSpaceSearch/Synthon.h>

namespace RDKit {
class ROMol;

namespace SynthonSpaceSearch {
class Synthon;
class SynthonSpace;
struct SynthonSpaceSearchParams;

// This class holds pointers to all the synthons for a particular
// reaction.  The synthons themselves are in a pool in the
// SynthonSpace.
class RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSet {
 public:
  SynthonSet() = default;
  explicit SynthonSet(const std::string &id) : d_id(id) {}
  SynthonSet(const SynthonSet &rhs) = delete;
  SynthonSet(SynthonSet &&rhs) = delete;

  const std::string &getId() const { return d_id; }
  const std::vector<std::vector<std::pair<std::string, Synthon *>>> &
  getSynthons() const {
    return d_synthons;
  }
  const boost::dynamic_bitset<> &getConnectors() const { return d_connectors; }
  const std::vector<boost::dynamic_bitset<>> &getSynthonConnectorPatterns()
      const {
    return d_synthConnPatts;
  }
  const std::vector<std::shared_ptr<ROMol>> &getConnectorRegions() const;
  const std::vector<std::string> &getConnectorRegionSmiles() const;
  const std::vector<std::unique_ptr<ExplicitBitVect>> &getConnRegFPs() const;
  const std::unique_ptr<ExplicitBitVect> &getAddFP() const;
  const std::unique_ptr<ExplicitBitVect> &getSubtractFP() const;
  const std::vector<int> &getNumConnectors() const;
  std::uint64_t getNumProducts() const;
  bool hasFingerprints() const;
  bool hasAddAndSubtractFPs() const;

  // Writes to/reads from a binary stream.
  void writeToDBStream(std::ostream &os) const;
  void readFromDBStream(std::istream &is, const SynthonSpace &space,
                        std::uint32_t version);
  // write the enumerated molecules to the stream in SMILES format.
  void enumerateToStream(std::ostream &os) const;

  // This stores the pointer to the Synthon, but doesn't manage
  // it and should never delete it.
  void addSynthon(int synthonSetNum, Synthon *newSynthon,
                  const std::string &synthonId);

  // Sometimes the synthon sets are numbered from 1 in the text file,
  // in which case there'll be an empty set 0.
  void removeEmptySynthonSets();

  // The bonds in the synthons may not be the same as in the products, and
  // this is a problem for aromatic ring creation in particular.  Such as:
  // [1*]=CC=C[2*] and [1*]Nc1c([2*])cccc1 giving c1ccc2ncccc2c1.  So
  // make versions of the synthons that reflect this, storead as searchMol
  // in each synthon.
  void makeSynthonSearchMols();

  // Build the connector regions and their fingerprints.  Only used when
  // creating a SynthonSpace from a text file.
  void buildConnectorRegions();

  // Scan through the connectors ([1*], [2*] etc.) in the synthons
  // and set bits in d_connectors accordingly.  Also removes any empty
  // reagent sets, which might be because the synthon numbers start from
  // 1 rather than 0.  Only used when creating a SynthonSpace from a text
  // file.
  void assignConnectorsUsed();

  void buildSynthonFingerprints(
      const FingerprintGenerator<std::uint64_t> &fpGen);
  void buildAddAndSubtractFPs(const FingerprintGenerator<std::uint64_t> &fpGen,
                              unsigned int numBits);

  // Return the molecules for synthons for which the bits are true.
  // Obviously requires that reqSynths is the same dimensions as
  // d_synthons.
  std::vector<std::vector<ROMol *>> getSynthons(
      const std::vector<boost::dynamic_bitset<>> &reqSynths) const;

  std::string buildProductName(const std::vector<size_t> &synthNums) const;
  std::unique_ptr<ROMol> buildProduct(
      const std::vector<size_t> &synthNums) const;

 private:
  std::string d_id;
  // The lists of synthons.  A product of the reaction is created by
  // combining 1 synthon from each of the outer vectors.  The actual
  // Synthon objects are held in the SynthonSpace which manages all
  // the memory.  In different reactions/SynthonSets the same Synthon
  // can have different IDs, so we need to keep the ID here rather
  // than in the Synthon, whose primary key is its SMILES string.
  std::vector<std::vector<std::pair<std::string, Synthon *>>> d_synthons;
  // MAX_CONNECTOR_NUM+1 bits showing which connectors are present in all the
  // synthon sets.
  boost::dynamic_bitset<> d_connectors;
  // And the connector patterns for each synthon set. If synthon set 0
  // has connectors 1 and 3, then d_synthConnPatts[0] will have bits
  // 1 and 3 set.
  std::vector<boost::dynamic_bitset<>> d_synthConnPatts;

  // The connector regions of a molecule are the pieces of up to 3 bonds from
  // a connector atom into the molecule.  We keep a vector of all the ones
  // present in the synthons in the set, plus a fingerprint for each.
  // If a query fragment doesn't have a connector region in common with
  // any of the synthons it can be assumed that the fragment won't have
  // a match in this SynthonSet.
  std::vector<std::shared_ptr<ROMol>> d_connectorRegions;
  std::vector<std::string> d_connRegSmis;
  // The fingerprints of the connector regions.
  std::vector<std::unique_ptr<ExplicitBitVect>> d_connRegFPs;

  // When doing an approximate FP similarity by ORing together
  // the synthonFPs, adding d_addFP and subtracting d_subtractFP
  // accounts (a bit) for the joins and the dummy atoms
  // respectively.
  std::unique_ptr<ExplicitBitVect> d_addFP;
  std::unique_ptr<ExplicitBitVect> d_subtractFP;

  // The number of connectors in the synthons in each synthon set.
  std::vector<int> d_numConnectors;
};

}  // namespace SynthonSpaceSearch

}  // namespace RDKit

#endif  // RDKIT_SYNTHONSET_H
