//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file and others here contain an implementation of
// synthonspace substructure search similar to that described in
// 'Fast Substructure Search in Combinatorial Library Spaces',
// Thomas Liphardt and Thomas Sander,
// J. Chem. Inf. Model. 2023, 63, 16, 5133–5141
// https://doi.org/10.1021/acs.jcim.3c00290

#ifndef RDKIT_SYNTHONSPACE_H
#define RDKIT_SYNTHONSPACE_H

/*! \file SynthonSpace.h

  \brief contains a class for searching combinatorial libraries in
         Synthon format such as Enamine REAL.

  \b Note that this functionality is experimental and the API may change
     in future releases.
*/

#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <RDGeneral/export.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearchHelpers.h>

namespace RDKit {
class ROMol;

namespace RascalMCES {
struct RascalOptions;
}

namespace SynthonSpaceSearch {

class Synthon;

class RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpace {
  friend class SynthonSet;
  friend class SynthonSpaceSearcher;
  friend class SynthonSpaceFingerprintSearcher;
  friend class SynthonSpaceRascalSearcher;
  friend class SynthonSpaceSubstructureSearcher;
  friend class SynthonSpaceShapeSearcher;

 public:
  explicit SynthonSpace() = default;
  ~SynthonSpace() = default;
  SynthonSpace(const SynthonSpace &other) = delete;
  SynthonSpace &operator=(const SynthonSpace &other) = delete;
  /*!
   * Get the number of different reactions in the SynthonSpace.
   *
   * @return int
   */
  size_t getNumReactions() const;

  /*!
   * Get the number of unique synthons in the SynthonSpace.  Synthons
   * can be used in more than 1 reaction.
   *
   * @return int
   */
  size_t getNumSynthons() const;
  /*!
   * Get a list of the names of all the reactions in the SynthonSpace.
   *
   * @return
   */
  std::vector<std::string> getReactionNames() const;
  std::shared_ptr<SynthonSet> getReaction(std::string reactionName);
  // The Synthons have a PatternFingerprint for screening in substructure
  // searches.  It's important that the screening process creates ones
  // of the same size, so this finds out what size that is.
  unsigned int getPatternFPSize() const;
  // Likewise for the fingerprints used for similarity searching
  unsigned int getFPSize() const;

  std::string getInputFileName() const;

  /*!
   * Get the total number of products that the SynthonSpace could produce.
   *
   * @return std::int64_t
   */
  std::uint64_t getNumProducts() const;

  /*!
   * Get the info string for the fingerprint generator used to
   * generate the stored fingerprints, so the user can query with
   * the same type.
   *
   * @return
   */
  std::string getSynthonFingerprintType() const { return d_fpType; }

  /*!
   * Perform a substructure search with the given query molecule across
   * the synthonspace library.  Duplicate SMILES strings produced by
   * different reactions will be returned.
   *
   * @param query : query molecule
   * @param params : (optional) settings for the search
   * @return : the hits as a SearchResults object.
   */
  SearchResults substructureSearch(
      const ROMol &query,
      const SubstructMatchParameters &matchParams = SubstructMatchParameters(),
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  /*!
   * Perform a substructure search with the given query molecule across
   * the synthonspace library.  Duplicate SMILES strings produced by
   * different reactions will be returned.  Search results are returned
   * incrementally through the provided callback, which will receive
   * at most `toTryChunkSize` sized lists of ROMols at a time, thereby
   * reducing the amount of memory required to hold search results.
   *
   * @param query : query molecule
   * @param callback: user-provided callback receiving chunks of ROMols.
   * @param params : (optional) settings for the search
   */
  void substructureSearch(
      const ROMol &query, const SearchResultCallback &callback,
      const SubstructMatchParameters &matchParams = SubstructMatchParameters(),
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  /*!
   * Perform a substructure search with the given generalized query
   * molecule across the synthonspace library.  Duplicate SMILES strings
   * produced by different reactions will be returned.
   *
   * @param query : query molecule
   * @param params : (optional) settings for the search
   * @return : the hits as a SearchResults object.
   */
  SearchResults substructureSearch(
      const GeneralizedSubstruct::ExtendedQueryMol &query,
      const SubstructMatchParameters &matchParams = SubstructMatchParameters(),
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  /*!
   * Perform a fingerprint similarity search with the given query molecule
   * across the synthonspace library.  Duplicate SMILES strings produced by
   * different reactions will be returned.
   * @param query : query molecule
   * @param fpGen: a FingerprintGenerator object that will provide the
   *               fingerprints for the similarity calculation
   * @param params : (optional) settings for the search
   * @return : the hits as a SearchResults object.
   */
  SearchResults fingerprintSearch(
      const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  // Perform a fingerprint similarity search with the given query molecule
  /*!
   * Perform a fingerprint similarity search with the given query molecule
   * across the synthonspace library.  Duplicate SMILES strings produced by
   * different reactions will be returned.  Search results are returned
   * incrementally through the provided callback, which will receive
   * at most `toTryChunkSize` sized lists of ROMols at a time, thereby
   * reducing the amount of memory required to hold search results.
   *
   * @param query : query molecule
   * @param fpGen: a FingerprintGenerator object that will provide the
   *               fingerprints for the similarity calculation
   * @param callback: user-provided callback receiving chunks of ROMols.
   * @param params : (optional) settings for the search
   */
  void fingerprintSearch(
      const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
      const SearchResultCallback &callback,
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  // Perform a RASCAL similarity search with the given query molecule
  // across the synthonspace library.  Duplicate SMILES strings produced by
  // different reactions will be returned.
  /*! Perform a RASCAL similarity search with the given query molecule
   * across the synthonspace library.  Duplicate SMILES strings produced by
   * different reactions will be returned.
   *
   * @param query : query molecule
   * @param rascalOptions: RASCAL options.  The similarityThreshold value
   *                       in the rascalOptions will be used rather than
   *                       params.similarityCutoff,
   *                       but params.fragSimilarityAdjuster will be used
   *                       to adjust the threshold for the fragment
   *                       comparisons.
   * @param params : (optional) settings for the search
   * @return : the hits as a SearchResults object.
   */
  SearchResults rascalSearch(
      const ROMol &query, const RascalMCES::RascalOptions &rascalOptions,
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  /*! Perform a shape similarity search with the given query molecule
   * across the synthonspace library.  Duplicate SMILES strings produced by
   * different reactions will be returned.
   *
   * @param query : query molecule
   * @param params : (optional) settings for the search
   * @return : the hits as a SearchResults object.
   */
  SearchResults shapeSearch(
      const ROMol &query,
      const SynthonSpaceSearchParams &params = SynthonSpaceSearchParams());

  // Perform a RASCAL similarity search with the given query molecule
  /* across the synthonspace library.  Duplicate SMILES strings produced by
   * different reactions will be returned.  Search results are returned
   * incrementally through the provided callback, which will receive
   * at most `toTryChunkSize` sized lists of ROMols at a time, thereby
   * reducing the amount of memory required to hold search results.
   *
   * @param query : query molecule
   * @param callback: user-provided callback receiving chunks of ROMols.
   * @param rascalOptions: RASCAL options.  The similarityThreshold value
   *                       in the rascalOptions will be used rather than
   *                       params.similarityCutoff,
   *                       but params.fragSimilarityAdjuster will be used
   *                       to adjust the threshold for the fragment
   *                       comparisons.
   * @param params : (optional) settings for the search
   */
  void rascalSearch(
      const ROMol &query, const RascalMCES::RascalOptions &rascalOptions,
      const SearchResultCallback &callback,
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
   * whitespace may be used as the field separator, but a tab is tried first
   * so that a tab-separated file may have spaces in the columns.
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
   * If it receives a SIGINT, returns cancelled=true.
   */
  void readTextFile(const std::string &inFilename, bool &cancelled);
  void readStream(std::istream &is, bool &cancelled);

  /*!
   * Writes to a binary DB File in our format.
   *
   * @param outFilename: the name of the file to write.
   */
  void writeDBFile(const std::string &outFilename) const;

  /*!
   * Reads from a binary DB File in our format.
   *
   * @param inFilename: the name of the file to read.
   * @param numThreads: number of threads to use in reading.  If negative,
   *                    adds the number to the number of hardware threads
   *                    available.
   */
  void readDBFile(const std::string &inFilename, int numThreads = 1);

  /*!
   * Write a summary of the SynthonSpace to given stream.
   *
   * @param os: stream
   */
  void summarise(std::ostream &os);

  /*!
   * Writes the enumerated library to file in SMILES format
   * (1 compound per line, SMILES name)
   *
   * @param outFilename: name of the file to write
   */
  void writeEnumeratedFile(const std::string &outFilename) const;
  void enumerateToStream(std::ostream &os) const;

  /*!
   * Create the fingerprints for the synthons ready for fingerprint searches.
   * Will be done by the fingerprint search if not done ahead of time.
   *
   * @param fpGen: a fingerprint generator of the appropriate type
   */
  void buildSynthonFingerprints(
      const FingerprintGenerator<std::uint64_t> &fpGen,
      unsigned int progressBarWidth = 0);

  /*!
   * Create conformers for the synthons ready for shape searching.  If a synthon
   * has unspecified stereochemistry, all possibilities will be enumerated and
   * shapes generated for each.
   *
   * @param shapeParams: controls the shape generation for each synthon
   */
  void buildSynthonShapes(bool &cancelled, const ShapeBuildParams &shapeParams =
                                               ShapeBuildParams());

  void reportSynthonUsage(std::ostream &os) const;
  std::uint64_t getNumSynthonsWithShapes() const;

 protected:
  unsigned int getMaxNumSynthons() const { return d_maxNumSynthons; }
  unsigned int getMaxNumConnectors() const;
  bool hasFingerprints() const;

  bool hasAddAndSubstractFingerprints() const;
  // Return whether the space contains a ring-forming reaction.
  bool getHasRingFormer() const { return d_hasRingFormer; }

  unsigned int getNumConformers() const;

  // Take the SMILES for a Synthon and if it's not in
  // d_synthonPool make it and add it.  If it is in the pool,
  // just look it up.  Either way, return a pointer to the
  // Synthon.
  Synthon *addSynthonToPool(const std::string &smiles);
  std::shared_ptr<SynthonSet> addReactionToPool(
      const std::string &reactionName);

  // Just do the lookup, and return nullptr if not found.
  Synthon *getSynthonFromPool(const std::string &smiles) const;

  void orderSynthonsForSearch();

 private:
  std::string d_fileName;
  // The reactions, keyed on their IDs as the first value
  // in the pair.
  std::vector<std::pair<std::string, std::shared_ptr<SynthonSet>>> d_reactions;
  // Keep the value of the maximum number of synthon sets used by
  // any of the reactions.  There's no point fragmenting any
  // query into more than this number of fragments.  Shouldn't
  // ever be higher than 4 at present.
  unsigned int d_maxNumSynthons{0};
  std::uint64_t d_numProducts{0};

  // This is actually 1000 * major version + 10 * minor
  // and hence the full version number.
  std::int32_t d_fileMajorVersion{-1};

  // The pool of all synthons, keyed on SMILES string.  Synthons
  // are frequently re-used in different reactions, so this means
  // they're only stored once.  They will be sorted and searched
  // for via first, which is its SMILES string.
  std::vector<std::pair<std::string, std::unique_ptr<Synthon>>> d_synthonPool;

  // This maps the Synthons (keyed by its SMILES) to the reactions/SynthonSets
  // they're in.
  std::unordered_map<std::string, std::vector<SynthonSet *>> d_synthonReactions;

  // For the fingerprint similarity search, this records the generator
  // used for creating synthon fingerprints that are read from a binary file.
  std::string d_fpType;
  // For the shape similarity search, this records the number of conformers
  // used.  It's not necessarily the same as the number of conformers each
  // synthon has.
  unsigned int d_numConformers{0};

  // Whether there is a ring-forming reaction in the space.
  bool d_hasRingFormer{false};

  SearchResults extendedSearch(const MolBundle &query,
                               const SubstructMatchParameters &matchParams,
                               const SynthonSpaceSearchParams &params);
  SearchResults extendedSearch(
      const GeneralizedSubstruct::ExtendedQueryMol::TautomerBundle_T &query,
      const SubstructMatchParameters &matchParams,
      const SynthonSpaceSearchParams &params);
  SearchResults extendedSearch(const TautomerQuery &query,
                               const SubstructMatchParameters &matchParams,
                               const SynthonSpaceSearchParams &params);

  // Fill up the d_synthonReactions object, listing all the reactions
  // each synthon is found in.
  void fillSynthonReactions();
  // For each synthon, build a sample molecule for each reaction it's in,
  // using d_synthonReactions.  Return them so the samples for a synthon
  // are in descending order of size so that when building shapes, we
  // do it on the smallest available example to save time, by popping
  // them off the back.  If maxSynthonAtoms > 0 and the synthon has
  // more heavy atoms than that (excluding dummy atoms) it is skipped.
  // The actual molecules in the SampleMolRecs aren't built until they're
  // needed, but all the information to do so is captured.
  void buildSynthonSampleMolecules(
      unsigned int maxSynthonAtoms,
      std::vector<std::vector<std::unique_ptr<SampleMolRec>>> &sampleMols)
      const;
};

/*!
 * Convert the text file into the binary DB file in our format.
 * Equivalent to readTextFile() followed by writeDBFile().
 * If a fingerprint generator is provided, fingerprints will
 * be created for all the synthons, which can be time-consuming.
 * @param inFilename name of the text file to read
 * @param outFilename name of the binary file to write
 * @param cancelled whether it received a SIGINT
 * @param fpGen optional fingerprint generator
 * @param shapeParams optional parameters for conformer generation
 */
RDKIT_SYNTHONSPACESEARCH_EXPORT void convertTextToDBFile(
    const std::string &inFilename, const std::string &outFilename,
    bool &cancelled, const FingerprintGenerator<std::uint64_t> *fpGen = nullptr,
    const ShapeBuildParams *shapeParams = nullptr);

/*!
 * Format an integer with spaces every 3 digits for ease
 * of reading.
 *
 * @return std::string
 */
RDKIT_SYNTHONSPACESEARCH_EXPORT std::string formattedIntegerString(
    std::int64_t value);

}  // namespace SynthonSpaceSearch
}  // namespace RDKit

#endif  // RDKIT_SYNTHONSPACE_H
