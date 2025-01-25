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
struct SynthonSpaceSearchParams;

// This class holds all the synthons for a particular reaction.
class RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSet {
 public:
  SynthonSet() = default;
  explicit SynthonSet(const std::string &id) : d_id(id) {}
  SynthonSet(const SynthonSet &rhs) = delete;
  SynthonSet(SynthonSet &&rhs) = delete;

  const std::string &getId() const { return d_id; }
  const std::vector<std::vector<std::unique_ptr<Synthon>>> &getSynthons()
      const {
    return d_synthons;
  }
  const boost::dynamic_bitset<> &getConnectors() const { return d_connectors; }
  const std::vector<boost::dynamic_bitset<>> &getSynthonConnectorPatterns()
      const {
    return d_synthConnPatts;
  }
  const std::vector<std::shared_ptr<ROMol>> &getConnectorRegions() const;

  const std::unique_ptr<ExplicitBitVect> &getConnRegFP() const;
  const std::unique_ptr<ExplicitBitVect> &getAddFP() const;
  const std::unique_ptr<ExplicitBitVect> &getSubtractFP() const;
  const std::vector<int> &getNumConnectors() const;
  bool hasFingerprints() const;
  bool hasAddAndSubtractFPs() const;

  const std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> &
  getSynthonFPs() const;
  // Writes to/reads from a binary stream.
  void writeToDBStream(std::ostream &os) const;
  void readFromDBStream(std::istream &is, std::uint32_t version);
  // write the enumerated molecules to the stream in SMILES format.
  void enumerateToStream(std::ostream &os) const;

  // SynthonSet takes control of the newSynthon and manages it.
  void addSynthon(int synthonSetNum, std::unique_ptr<Synthon> newSynthon);

  // Sometimes the synthon sets are numbered from 1 in the text file,
  // in which case there'll be an empty set 0.
  void removeEmptySynthonSets();

  // The bonds in the synthons may not be the same as in the products, and
  // this is a problem for aromatic ring creation in particular.  Such as:
  // [1*]=CC=C[2*] and [1*]Nc1c([2*])cccc1 giving c1ccc2ncccc2c1.  So
  // transfer the types of bonds from the products to the synthons.
  void transferProductBondsToSynthons();

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
  void buildAddAndSubtractFPs(const FingerprintGenerator<std::uint64_t> &fpGen);

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
  // combining 1 synthon from each of the outer vectors.
  std::vector<std::vector<std::unique_ptr<Synthon>>> d_synthons;
  // 4 bits showing which connectors are present in all the
  // synthon sets.
  boost::dynamic_bitset<> d_connectors;
  // and the connector patterns for each synthon set.
  std::vector<boost::dynamic_bitset<>> d_synthConnPatts;

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

  // When doing an approximate FP similarity by ORing together
  // the synthonFPs, adding d_addFP and subtracting d_subtractFP
  // accounts (a bit) for the joins and the dummy atoms
  // respectively.
  std::unique_ptr<ExplicitBitVect> d_addFP;
  std::unique_ptr<ExplicitBitVect> d_subtractFP;

  // The number of connectors in the synthons in each synthon set.
  std::vector<int> d_numConnectors;

  // The fingerprints for the synthons for use with a fingerprint similarity
  // search. They are not properties of the Synthons because they are not
  // generated directly from them, as explained in buildSynthonFingerprints.
  std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> d_synthonFPs;

  // Tag each atom and bond in each synthon with its index and the synthon
  // set number it came from.
  void tagSynthonAtomsAndBonds() const;
};

}  // namespace SynthonSpaceSearch

}  // namespace RDKit

#endif  // RDKIT_SYNTHONSET_H
