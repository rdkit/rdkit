//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef SYNTHONSPACESEARCHHELPERS_H
#define SYNTHONSPACESEARCHHELPERS_H

#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>
#include <GraphMol/SynthonSpaceSearch/SearchShapeInput.h>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit::SynthonSpaceSearch {

// This the maximum number of connectors that we can deal with at the moment.
// In reality, there may be fewer than this.  However, the key limit is in
// The symbols used for the connectors in Enamine REAL etc.
const std::vector<std::string> CONNECTOR_SYMBOLS{"[U]", "[Np]", "[Pu]", "[Am]"};
constexpr unsigned int MAX_CONNECTOR_NUM{4};

struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceSearchParams {
  std::int64_t maxHits{1000};  // The maximum number of hits to return.
                               // Use -1 for no maximum.
  std::uint64_t maxNumFragSets{
      100000};  // The maximum number of fragment sets the query can
                // be broken into.  Big molecules will create huge
                // numbers of fragment sets that may cause excessive
                // memory use.  If the number of fragment sets hits this
                // number, fragmentation stops and the search results
                // will likely be incomplete.
  std::int64_t toTryChunkSize{
      2500000};              // For similarity searching, especially
                             // fingerprint similarity, there can be a
                             // very large number of possible hits to
                             // screen which can use a lot of memory and
                             // crash the program.  It will also be very
                             // slow.  To alleviate the memory use, the
                             // possible hits are processed in chunks.
                             // This parameter sets the chunk size.
  std::int64_t hitStart{0};  // Sequence number of hit to start from.  So that
                             // you can return the next N hits of a search
                             // having already obtained N-1.
  bool randomSample{false};  // If true, returns a random sample of the hit
                             // hits, up to maxHits in number.
  int randomSeed{-1};        // Seed for random-number generator.  -1 means use
                             // a random seed (std::random_device).
  bool buildHits{true};  // If false, reports the maximum number of hits that
                         // the search could produce, but doesn't return them.
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
             // positive correlation with your FOMO.
  unsigned int minHitHeavyAtoms{0};  // Minimum number of heavy atoms in a hit.
  int maxHitHeavyAtoms{-1};          // Maximum number of heavy atoms in a hit.
  // -1 means no maximum.
  double minHitMolWt{0};  // Minimum molecular weight for a hit.
  double maxHitMolWt{0};  // Maximum molecular weight for a hit.  0.0 means
  // no maximum.
  unsigned int minHitChiralAtoms{
      0};                           // Minimum number of chiral atoms in a hit.
  int maxHitChiralAtoms{-1};        // Maximum number of chiral atoms in a hit.
                                    // -1 means no maximum.
  unsigned int numConformers{100};  // When doing a shape search, the number of
                                    // conformers to use for each molecule.
  double confRMSThreshold{0.5};  // When doing a shape search, the RMS threshold
                                 // to use when pruning conformers.  Passed
                                 // directly EmbedMultipleConfs.
  double shapePruneThreshold{
      0.95};  //! When doing shape searching, the shapes will be pruned so that
              //! no 2 shapes are more similar than this threshold.
  bool bestHit{false};  // If true, when doing a shape search it will return the
                        // hit conformer with the best shape match to a query
                        // conformer.  If false it will just return the first
                        // hit conformer that exceeds the similarity cutoff.
                        // The latter will be faster but the returned hit
                        // conformations likely to be less relevant.
  bool enumerateUnspecifiedStereo{
      false};  // When doing a shape search, if there is
               // unspecified stereochemistry in either
               // the query or potential hit, enumerate and
               // test all possibilities.
  EnumerateStereoisomers::StereoEnumerationOptions stereoEnumOpts{
      true, true, false, true,
      0,    0xdac};  // Options for stereoisomer enumeration.  Over-ride default
                     // tryEmbedding of false and use fixed randomSeed for
                     // consistency of results from run to run.
  std::uint64_t timeOut{600};  // Maximum number of seconds to spend on a single
                               // search.  0 means no maximum.
  int numThreads = 1;  // The number of threads to use.  If > 0, will use that
                       // number.  If <= 0, will use the number of hardware
                       // threads plus this number.  So if the number of
                       // hardware threads is 8, and numThreads is -1, it will
                       // use 7 threads.
  unsigned int useProgressBar{0};  // Makes a progress bar of given width.  The
                                   // number given is the number of '*'
                                   // characters in a full bar.  There will be
                                   // about another 35 characters or so
                                   // depending on the size of the job.  0
                                   // means no bar.
};

// Options to be passed to buildSynthonShapes.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT ShapeBuildParams {
  // The relevant ones are passed directly into EmbedMultipleConfs.
  unsigned int numConfs{100};  // Max number of conformations per synthon
  double rmsThreshold{0.5};    // RMS threshold used when pruning conformations
  double shapeSimThreshold{0.95};  // This is passed to pruneShapes(). For each
                                   // synthon, no 2 shapes will be more similar
                                   // to each other than the threshold.
  int numThreads{1};  // The number of threads to use.  If > 0, will use that
                      // number.  If <= 0, will use the number of hardware
                      // threads plus this number.
  int randomSeed{
      0xdac};  // Seed for random number generator.  Fixed by default
               // so the same synthon conformers are produced each time.
  EnumerateStereoisomers::StereoEnumerationOptions stereoEnumOpts{
      true, true, false, true,
      0,    -1};  // Options for stereoisomer enumeration.  Over-ride default
                  // tryEmbedding of false.
  unsigned int useProgressBar{0};   // Makes a progress bar of given width.  The
                                    // number given is the number of '*'
                                    // characters in a full bar.  There will be
                                    // about another 35 characters or so
                                    // depending on the size of the job.  0
                                    // means no bar.
  unsigned int maxSynthonAtoms{0};  // If > 0, sets a maximum number of heavy
                                    // atoms, excluding dummies, for a synthon
                                    // to have a shape made.
  unsigned int maxEmbedAttempts{10};  // Maximum attempts for an embedding.
  unsigned int timeOut{60};  // Maximum time in seconds to spend on each synthon
                             // when generating conformers.
  std::string
      interimFile;  // Interim file to write SynthonSpace to.  In the event of
                    // a failure, a restart from this file will be possible.
  std::uint64_t interimWrites{
      1000};  // If an interim file has been given, every this
              // many shapes write a new version of the file.  If 0, don't
              // do any writing.
};

using ShapeSet = std::vector<std::unique_ptr<GaussianShape::SearchShapeInput>>;

// Experiments have shown that it's quicker to parallelise the search
// at the SynthonSet level for substructure searches, but at the
// FragSet level for the slower searches such as fingerprints, Rascal
// and shape.  This enum is used to control which is done in a particular
// case.
enum class ThreadMode : unsigned char {
  ThreadReactions,
  ThreadFragments
};

}  // namespace RDKit::SynthonSpaceSearch
#endif  // SYNTHONSPACESEARCHHELPERS_H
