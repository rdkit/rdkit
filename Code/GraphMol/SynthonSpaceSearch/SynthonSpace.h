//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_SYNTHONSPACE_H
#define RDKIT_SYNTHONSPACE_H

/*! \file SynthonSpace.h

  \brief contains a class for searching combinatorial libraries in
         Synthon format such as Enamine REAL.

  \b Note that this functionality is experimental and the API may change
     in future releases.
*/

#include <map>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <RDGeneral/export.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>

namespace RDKit {
class ROMol;

namespace SynthonSpaceSearch {

// This the maximum number of connectors that we can deal with at the moment.
// In reality, there may be fewer than this.  However, the key limit is in
// The symbols used for the connectors in Enamine REAL etc.
const std::vector<std::string> CONNECTOR_SYMBOLS{"[U]", "[Np]", "[Pu]", "[Am]"};
constexpr unsigned int MAX_CONNECTOR_NUM{4};

struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceSearchParams {
  int maxBondSplits{MAX_CONNECTOR_NUM};  // The maximum number of bonds to break
                                         // in the query. It should be no more
                                         // than the maximum number of connector
                                         // types in the SynthonSpace.  At
                                         // present this is 4.  Specifying more
                                         // than that will not matter as it will
                                         // be reduced to 4.  Likewise, values
                                         // lower than 1 will be increased to 1.
  std::uint64_t maxNumFrags{
      100000};  // The maximum number of fragments the query can
                // be broken into.  Big molecules will create huge
                // numbers of fragments that may cause excessive
                // memory use.  If the number of fragments hits this number,
                // fragmentation stops and the search results will likely be
                // incomplete.
  std::int64_t maxHits{1000};  // The maximum number of hits to return.  Use
                               // -1 for no maximum.
  std::int64_t hitStart{0};    // Sequence number of hit to start from.  So that
                               // you can return the next N hits of a search
                               // having already obtained N-1.
  bool randomSample{false};    // If true, returns a random sample of the hit
                               // hits, up to maxHits in number.
  int randomSeed{-1};    // Seed for random-number generator.  -1 means use
                         // a random seed (std::random_device).
  bool buildHits{true};  // If false, reports the maximum number of hits that
                         // the search could produce, but doesn't return them.
  int numRandomSweeps{10};  // The random sampling doesn't always produce the
                            // required number of hits in 1 go.  This parameter
                            // controls how many loops it makes to try and get
                            // the hits before giving up.
  double similarityCutoff{0.5};  // Similarity cutoff for returning hits by
                                 // fingerprint similarity.  The default is
                                 // appropriate for a Morgan fingerprint of
                                 // radius=2, it may need changing for other
                                 // fingerprint types.
  double fragSimilarityAdjuster{
      0.1};  // Similarity values for fragments are generally low
             // due to low bit densities.  For the fragment
             // matching, reduce the similarity cutoff
             // by this amount.  A higher number will give slower search
             // times, a lower number will give faster searches at the
             // risk of missing some hits.  The value you use should have
             // a positive correlation with your FOMO.
  double approxSimilarityAdjuster{
      0.1};  // The fingerprint search uses an approximate similarity method
             // before building a product and doing a final check.  The
             // similarityCutoff is reduced by this value for the approximate
             // check.  A lower value will give faster run times at the
             // risk of missing some hits.  The value you use should have a
             // positive correlation with your FOMO.  The default is
             // appropriate for Morgan fingerprints.  With RDKit fingerprints,
             // 0.05 is adequate, and higher than that has been seen to
             // produce long run times.
  std::uint64_t timeOut{600};  // Maximum number of seconds to spend on a single
                               // search.  0 means no maximum.
};

// Holds the information about a set of hits.  The molecules can be built
// by making all combinations of synthons, one taken from each synthon set.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceHitSet {
  SynthonSpaceHitSet() = delete;
  SynthonSpaceHitSet(const std::string &id,
                     const std::vector<std::vector<size_t>> &stu)
      : reactionId(id), synthonsToUse(stu) {
    numHits = std::accumulate(
        synthonsToUse.begin(), synthonsToUse.end(), size_t(1),
        [](const int prevRes, const std::vector<size_t> &s2) -> size_t {
          return prevRes * s2.size();
        });
  }
  std::string reactionId;
  std::vector<std::vector<size_t>> synthonsToUse;
  size_t numHits{0};
};

class RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpace {
 public:
  // Create the synthonspace from a file in the correct format.
  explicit SynthonSpace() = default;
  SynthonSpace(const SynthonSpace &other) = delete;
  SynthonSpace &operator=(const SynthonSpace &other) = delete;
  // Get the number of different reactions in the SynthonSpace.
  /*!
   *
   * @return int
   */
  size_t getNumReactions() const { return d_reactions.size(); }
  const std::map<std::string, std::unique_ptr<SynthonSet>> &getReactions()
      const {
    return d_reactions;
  }

  // Get the total number of products that the SynthonSpace could produce.
  /*!
   *
   * @return std::int64_t
   */
  std::int64_t getNumProducts() const;

  std::string getSynthonFingerprintType() const { return d_fpType; }

  // Perform a substructure search with the given query molecule across
  // the synthonspace library.  Duplicate SMILES strings produced by different
  // reactions will be returned.
  /*!
   *
   * @param query : query molecule
   * @param params : (optional) settings for the search
   * @return : the hits as a SearchResults object.
   */
  SearchResults substructureSearch(
      const ROMol &query,
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  // Perform a fingerprint similarity search with the given query molecule
  // across the synthonspace library.  Duplicate SMILES strings produced by
  // different reactions will be returned.
  /*!
   *
   * @param query : query molecule
   * @param fpGen: a FingerprintGenerator object that will provide the
   *               fingerprints for the similarity calculation
   * @param params : (optional) settings for the search
   * @return : the hits as a SearchResults object.
   */
  SearchResults fingerprintSearch(
      const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  /*!
   *
   * @param inFilename: name of the file containing the synthon-based library.
   *
   * The original format is:
   * all lines are tab-separated
   * first line:SMILES	synton_id	synton#	reaction_id
   * Note the spelling "synton" from the original paper/example file.
   * Subsequent lines have a single reagent e.g.
   * OCC([U])=NN=[Np]		1-1	0	triazole-1
   * C1CCCC1N([Pu])[U]		2-1	1	triazole-1
   * CC1CCN(C1)C(=[Np])[Pu]		3-1	2	triazole-1
   *
   * Other acceptable formats are as above, but with a 5th column "release":
   * SMILES	synton_id	synton#	reaction_id release
   *
   * or a comma-separated equivalent of the first format:
   * SMILES,synton_id,synton_role,reaction_id
   * but with the 3rd column named differently but with the same meaning.
   * The formatting of the first 2 formats has been relaxed such that any
   * whitespace may be used as the field separator.
   *
   * Attachment points are U, Np, Pu and Am for up to 4 synthons per reaction.
   * A product is created by taking a synthon from each synton# value and
   * combining by replacing matching trans-uranic elements and replacing them
   * with a direct bond of the appropriate type.
   * A more (for RDKit) conventional connection flag of isotope labelled
   * dummy atoms is also accepted ([1*] etc.).
   * Throws a std::runtime_error if it doesn't think the format is correct,
   * which it does by checking that the first line is as above and subsequent
   * lines have appropriate number of fields.
   */
  void readTextFile(const std::string &inFilename);

  // Writes to/reads from a binary DB File in our format.
  /*!
   *
   * @param outFilename: the name of the file to write.
   */
  void writeDBFile(const std::string &outFilename) const;
  /*!
   *
   * @param inFilename: the name of the file to read.
   */
  void readDBFile(const std::string &inFilename);

  // Write a summary of the SynthonSpace to given stream.
  /*!
   *
   * @param os: stream
   */
  void summarise(std::ostream &os) const;

  // Writes the enumerated library to file in SMILES format (1 compound
  // per line, SMILES name
  /*!
  @param outFilename: name of the file to write
   */
  void writeEnumeratedFile(const std::string &outFilename) const;

  bool hasFingerprints() const;
  // Create the fingerprints for the synthons ready for fingerprint searches.
  // Will be done by the fingerprint search if not done ahead of time.
  void buildSynthonFingerprints(
      const FingerprintGenerator<std::uint64_t> &fpGen);

  bool hasAddAndSubstractFingerprints() const;
  // Create the add and substract fingerprints for the SynthonSets.
  // Will be done by the fingerprint search if not done ahead of time.
  void buildAddAndSubstractFingerprints(
      const FingerprintGenerator<std::uint64_t> &fpGen);

 private:
  std::string d_fileName;
  std::map<std::string, std::unique_ptr<SynthonSet>> d_reactions;

  // For the similarity search, this records the generator used for
  // creating synthon fingerprints that are read from a binary file.
  std::string d_fpType;
};

}  // namespace SynthonSpaceSearch
}  // namespace RDKit

#endif  // RDKIT_SYNTHONSPACE_H
