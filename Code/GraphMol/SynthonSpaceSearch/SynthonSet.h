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
struct ShapeBuildParams;

// For holding a sample molecule from this set, based on the given
// synthon which is in the d_synthonSetNum vector of d_synthons
// of the SynthonSet.  Used for building shapes.  The d_mol isn't
// built immediately to save memory, so the information needed to
// produce it is all captured.
struct SampleMolRec {
  const class SynthonSet *d_synthonSet{nullptr};
  Synthon *d_synthon{nullptr};
  std::vector<size_t> d_synthonNums;
  size_t d_synthonSetNum;
  std::unique_ptr<ROMol> d_mol{nullptr};
  unsigned int d_numAtoms{0};
};

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
  // Return the synthon of the given number in the given set, based on
  // the sorted order rather than the input order.  Returns empty string
  // and nullptr if the numbers are invalid.
  const std::pair<std::string, Synthon *> getSubstructureOrderedSynthon(
      size_t setNum, size_t synthonNum) const;
  const std::pair<std::string, Synthon *> getFingerprintOrderedSynthon(
      size_t setNum, size_t synthonNum) const;
  const std::pair<std::string, Synthon *> getRascalOrderedSynthon(
      size_t setNum, size_t synthonNum) const;
  // Return the actual synthon number given its ordered number.  Returns
  // "-1" if the numbers are invalid.
  size_t getSubstructureOrderedSynthonNum(size_t setNum,
                                          size_t synthonNum) const;
  size_t getFingerprintOrderedSynthonNum(size_t setNum,
                                         size_t synthonNum) const;
  size_t getRascalOrderedSynthonNum(size_t setNum, size_t synthonNum) const;

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
  unsigned int getNumRingFormers() const { return d_numRingFormers; }

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
  // make versions of the synthons that reflect this, stored as searchMol
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

  void initializeSearchOrders();
  std::vector<std::vector<size_t>> orderSynthonsForSearch(
      const std::function<bool(const Synthon *synthon1,
                               const Synthon *synthon2)> &cmp);

  // make a SampleMolRec for this synthon, which is expected
  // to be in the SynthonSet.  Returns an empty object if it isn't,
  // but that really shouldn't happen.
  std::unique_ptr<SampleMolRec> makeSampleMolecule(Synthon *synthon) const;
  // Make a molecule from the given synthonNums, which index into d_synthons.
  std::unique_ptr<ROMol> buildMolecule(
      const std::vector<size_t> &synthonNums) const;

  void assessRingFormers();

 private:
  std::string d_id;
  // The lists of synthons.  A product of the reaction is created by
  // combining 1 synthon from each of the outer vectors.  The actual
  // Synthon objects are held in the SynthonSpace which manages all
  // the memory.  In different reactions/SynthonSets the same Synthon
  // can have different IDs, so we need to keep the ID here rather
  // than in the Synthon, whose primary key is its SMILES string.
  std::vector<std::vector<std::pair<std::string, Synthon *>>> d_synthons;

  // The order that the synthons should be searched in.  It will
  // be the same shape as d_synthons. The substructure
  // search will order in ascending size of number of heavy atoms,
  // the fingerprint search in ascending number of set bits and
  // the Rascal search in ascending number of heavy atoms and bonds.
  // This allows for a marginally more efficient search since either
  // or both ends of a synthon set can be ignored as not being able
  // to yield a match. There is no advantage in sorting the shapes.
  std::vector<std::vector<size_t>> d_substructureSearchOrders;
  std::vector<std::vector<size_t>> d_fingerprintSearchOrders;
  std::vector<std::vector<size_t>> d_rascalSearchOrders;

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
  // respectively.  Not used at present, but left in for
  // compatibility with old databases.
  std::unique_ptr<ExplicitBitVect> d_addFP;
  std::unique_ptr<ExplicitBitVect> d_subtractFP;

  // The number of connectors in the synthons in each synthon set.
  std::vector<int> d_numConnectors;

  // Synthons are shared, so sometimes we need to copy the molecules into a
  // new set that we can fiddle with without upsetting anything else.
  std::unique_ptr<RWMol> copySynthon(size_t synthonSetNum,
                                     size_t synthonIdx) const;
  std::vector<std::vector<std::unique_ptr<RWMol>>> copySynthons() const;

  // Take the synthons and build molecules from them.  longVecNum is the number
  // of the vector containing the synthon set of interest.  Each product is
  // a molecule made from the corresponding member of longVecNum and a small
  // element of the other vectors.  synthonMols are from copySynthons.
  std::vector<std::unique_ptr<ROMol>> buildSampleMolecules(
      const std::vector<std::vector<std::unique_ptr<RWMol>>> &synthonMols,
      size_t longVecNum) const;
  // The number of rings that may be formed by the synthons.  If there
  // are a pair of synthons A([1*])[2*] and B([1*])[2*] 1 ring can be
  // formed.
  unsigned int d_numRingFormers{0};
};

}  // namespace SynthonSpaceSearch

}  // namespace RDKit

#endif  // RDKIT_SYNTHONSET_H
