//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>

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
  if (!getSpace().hasAddAndSubstractFingerprints()) {
    getSpace().buildAddAndSubstractFingerprints(fpGen);
  }
  d_queryFP = std::unique_ptr<ExplicitBitVect>(d_fpGen.getFingerprint(query));
}

namespace {
// Take the fragged mol fps and flag all those synthons that have a fragment as
// a similarity match.
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<std::unique_ptr<ExplicitBitVect>> &fragFPs,
    const double similarityCutoff, const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;
  std::vector<std::vector<std::pair<size_t, double>>> fragSims(
      reaction->getSynthons().size());

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
        fragSims[synthonOrder[i]].emplace_back(j, sim);
        fragMatched = true;
      }
    }
    if (!fragMatched) {
      // No synthons matched this fragment, so the whole fragment set is a
      // bust.
      return retSynthons;
    }
  }

  // Fill in any synthons where they all didn't match because there were
  // fewer fragments than synthons.
  details::expandBitSet(synthonsToUse);
  details::bitSetsToVectors(synthonsToUse, retSynthons);

  // Now order the synthons in descending order of their similarity to
  // the corresponding fragFP.
  for (size_t i = 0; i < fragFPs.size(); i++) {
    if (fragSims[i].empty()) {
      // This one will have been filled in by expandBitSet so we need to use
      // all the synthons and a dummy similarity.
      fragSims[i].resize(synthonsToUse[i].size());
      for (size_t j = 0; j < fragSims[i].size(); j++) {
        fragSims[i][j] = std::make_pair(j, 0.0);
      }
    } else {
      std::sort(
          fragSims[i].begin(), fragSims[i].end(),
          [](const auto &a, const auto &b) { return a.second > b.second; });
    }
    retSynthons[i].clear();
    std::transform(
        fragSims[i].begin(), fragSims[i].end(),
        std::back_inserter(retSynthons[i]),
        [](const std::pair<size_t, double> &fs) { return fs.first; });
  }
  return retSynthons;
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
          SynthonSpaceHitSet hs{reaction->getId(), theseSynthons};
          if (hs.numHits) {
            results.push_back(hs);
          }
        }
      }
    }
  }
  return results;
}

bool SynthonSpaceFingerprintSearcher::quickVerify(
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<size_t> &synthNums) const {
  // Make an approximate fingerprint by combining the FPs for
  // these synthons, adding in the addFP and taking out the
  // subtractFP.
  const auto &synthFPs = reaction->getSynthonFPs();
  ExplicitBitVect fullFP(*synthFPs[0][synthNums[0]]);
  for (unsigned int i = 1; i < synthNums.size(); ++i) {
    fullFP |= *synthFPs[i][synthNums[i]];
  }
  fullFP |= *(reaction->getAddFP());
  // The subtract FP has already had its bits flipped, so just do a
  // straight AND.
  fullFP &= *(reaction->getSubtractFP());

  double approxSim = TanimotoSimilarity(fullFP, *d_queryFP);
  return approxSim >=
         getParams().similarityCutoff - getParams().approxSimilarityAdjuster;
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