//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <DataStructs/BitOps.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceFingerprintSearcher.h>

namespace RDKit::SynthonSpaceSearch {

SynthonSpaceFingerprintSearcher::SynthonSpaceFingerprintSearcher(
    const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
    const SynthonSpaceSearchParams &params, SynthonSpace &space)
    : SynthonSpaceSearcher(query, params, space), d_fpGen(fpGen) {
  if (!getSpace().hasFingerprints() ||
      getSpace().getSynthonFingerprintType() != fpGen.infoString()) {
    getSpace().buildSynthonFingerprints(fpGen);
  }
  d_queryFP = std::unique_ptr<ExplicitBitVect>(d_fpGen.getFingerprint(query));
}

namespace {
// Take the fragged mol fps and flag all those synthons that have a fragment as
// a similarity match.
std::vector<boost::dynamic_bitset<>> getHitSynthons(
    const std::vector<std::unique_ptr<ExplicitBitVect>> &fragFPs,
    const double similarityCutoff, const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  synthonsToUse.reserve(reaction->getSynthons().size());
  for (const auto &synthonSet : reaction->getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  for (size_t i = 0; i < synthonOrder.size(); i++) {
    const auto &synthonFPs = reaction->getSynthonFPs()[synthonOrder[i]];
    bool fragMatched = false;
    for (size_t j = 0; j < synthonFPs.size(); j++) {
      if (const auto sim = TanimotoSimilarity(*fragFPs[i], *synthonFPs[j]);
          sim >= similarityCutoff) {
        synthonsToUse[synthonOrder[i]][j] = true;
        fragMatched = true;
      }
    }
    if (!fragMatched) {
      // No synthons matched this fragment, so the whole fragment set is a
      // bust.
      synthonsToUse.clear();
      return synthonsToUse;
    }
  }

  // Fill in any synthons where they all didn't match because there were
  // fewer fragments than synthons.
  details::expandBitSet(synthonsToUse);
  return synthonsToUse;
}
}  // namespace

std::vector<SynthonSpaceHitSet> SynthonSpaceFingerprintSearcher::searchFragSet(
    std::vector<std::unique_ptr<ROMol>> &fragSet) const {
  std::vector<SynthonSpaceHitSet> results;
  std::vector<std::unique_ptr<ExplicitBitVect>> fragFPs;
  fragFPs.reserve(fragSet.size());
  for (auto &frag : fragSet) {
    // For the fingerprints, ring info is required.
    unsigned int otf;
    sanitizeMol(*static_cast<RWMol *>(frag.get()), otf,
                MolOps::SANITIZE_SYMMRINGS);
    fragFPs.emplace_back(d_fpGen.getFingerprint(*frag));
  }

  const auto connPatterns = details::getConnectorPatterns(fragSet);
  boost::dynamic_bitset<> conns(MAX_CONNECTOR_NUM + 1);
  for (auto &connPattern : connPatterns) {
    conns |= connPattern;
  }

  for (const auto &[id, reaction] : getSpace().getReactions()) {
    // It can't be a hit if the number of fragments is more than the number
    // of synthon sets because some of the molecule won't be matched in any
    // of the potential products.  It can be less, in which case the unused
    // synthon set will be used completely, possibly resulting in a large
    // number of hits.
    if (fragSet.size() > reaction->getSynthons().size()) {
      continue;
    }
    auto synthConnPatts = reaction->getSynthonConnectorPatterns();

    // Need to try all combinations of synthon orders.
    auto synthonOrders =
        details::permMFromN(fragSet.size(), reaction->getSynthons().size());
    for (const auto &synthonOrder : synthonOrders) {
      // Get all the possible permutations of connector numbers compatible with
      // the number of synthon sets in this reaction.  So if the
      // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
      // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
      // that might be how they're labelled in the reaction database.
      auto connCombs = details::getConnectorPermutations(
          fragSet, conns, reaction->getConnectors());

      for (auto &connComb : connCombs) {
        // Make sure that for this connector combination, the synthons in this
        // order have something similar.  All query fragment connectors must
        // match something in the corresponding synthon.  The synthon can
        // have unused connectors.
        auto connCombConnPatterns = details::getConnectorPatterns(connComb);
        bool skip = false;
        for (size_t i = 0; i < connCombConnPatterns.size(); ++i) {
          if ((connCombConnPatterns[i] & synthConnPatts[synthonOrder[i]])
                  .count() < connCombConnPatterns[i].count()) {
            skip = true;
            break;
          }
        }
        if (skip) {
          continue;
        }

        // It appears that for fingerprints, the isotope numbers are
        // ignored so there's no need to worry about the connector numbers
        // in the fingerprints.
        auto theseSynthons = getHitSynthons(
            fragFPs,
            getParams().similarityCutoff - getParams().fragSimilarityAdjuster,
            reaction, synthonOrder);
        if (!theseSynthons.empty()) {
          const size_t numHits = std::accumulate(
              theseSynthons.begin(), theseSynthons.end(), 1,
              [](const int prevRes, const boost::dynamic_bitset<> &s2) {
                return prevRes * s2.count();
              });
          if (numHits) {
            results.push_back(
                SynthonSpaceHitSet{reaction->getId(), theseSynthons, numHits});
          }
        }
      }
    }
  }
  return results;
}

bool SynthonSpaceFingerprintSearcher::verifyHit(const ROMol &hit) const {
  const std::unique_ptr<ExplicitBitVect> fp(d_fpGen.getFingerprint(hit));
  if (const auto sim = TanimotoSimilarity(*fp, *d_queryFP);
      sim >= getParams().similarityCutoff) {
    hit.setProp<double>("Similarity", sim);
    return true;
  }
  return false;
}

}  // namespace RDKit::SynthonSpaceSearch