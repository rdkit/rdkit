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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceFingerprintSearcher.h>
#include <RDGeneral/ControlCHandler.h>

namespace RDKit::SynthonSpaceSearch {

SynthonSpaceFingerprintSearcher::SynthonSpaceFingerprintSearcher(
    const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
    const SynthonSpaceSearchParams &params, SynthonSpace &space)
    : SynthonSpaceSearcher(query, params, space), d_fpGen(fpGen) {
  if (getSpace().hasFingerprints() &&
      d_fpGen.infoString() != getSpace().getSynthonFingerprintType()) {
    throw std::runtime_error(
        "The search fingerprints must match"
        " those in the database.  You are searching with " +
        d_fpGen.infoString() + " vs " + getSpace().getSynthonFingerprintType() +
        " in the database.");
  }
  d_queryFP = std::unique_ptr<ExplicitBitVect>(d_fpGen.getFingerprint(query));
}

namespace {
// Take the fragged mol fps and flag all those synthons that have a fragment as
// a similarity match.
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<ExplicitBitVect *> &fragFPs,
    const double similarityCutoff, const SynthonSet &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;
  std::vector<std::vector<std::pair<size_t, double>>> fragSims(
      reaction.getSynthons().size());

  synthonsToUse.reserve(reaction.getSynthons().size());
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  for (size_t i = 0; i < synthonOrder.size(); i++) {
    const auto &synthons = reaction.getSynthons()[synthonOrder[i]];
    bool fragMatched = false;
    for (size_t j = 0; j < synthons.size(); j++) {
      // There's a simple calculation for the maximum possible Tanimoto
      // Coefficient that these 2 fingerprints can achieve.
      const double maxSim =
          fragFPs[i]->getNumOnBits() <
                  synthons[j].second->getFP()->getNumOnBits()
              ? static_cast<double>(fragFPs[i]->getNumOnBits()) /
                    synthons[j].second->getFP()->getNumOnBits()
              : static_cast<double>(
                    synthons[j].second->getFP()->getNumOnBits()) /
                    fragFPs[i]->getNumOnBits();
      if (maxSim < similarityCutoff) {
        continue;
      }
      if (const auto sim =
              TanimotoSimilarity(*fragFPs[i], *synthons[j].second->getFP());
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
    std::transform(fragSims[i].begin(), fragSims[i].end(),
                   std::back_inserter(retSynthons[i]),
                   [](const auto &fs) { return fs.first; });
  }
  return retSynthons;
}
}  // namespace

void SynthonSpaceFingerprintSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) {
  if (!getSpace().hasFingerprints() ||
      getSpace().getSynthonFingerprintType() != d_fpGen.infoString()) {
    getSpace().buildSynthonFingerprints(d_fpGen);
  }
  if (ControlCHandler::getGotSignal()) {
    return;
  }

  // Slightly convoluted way of doing it to prepare for multi-threading.
  // Make a map of the unique SMILES strings for the fragments, keeping
  // track of them in the vector.
  bool cancelled = false;
  auto fragSmiToFrag = details::mapFragsBySmiles(fragSets, cancelled);
  if (cancelled) {
    return;
  }

  // Now generate the fingerprints for the fragments.  This is the
  // time-consuming bit that will be threaded.
  d_fragFPPool.resize(fragSmiToFrag.size());
  unsigned int fragNum = 0;
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    if (ControlCHandler::getGotSignal()) {
      return;
    }
    d_fragFPPool[fragNum++].reset(d_fpGen.getFingerprint(*frags.front()));
  }

  // Now use the pooled fps to populate the vectors for each fragSet.
  fragNum = 0;
  d_fragFPs.reserve(fragSmiToFrag.size());
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    for (auto &frag : frags) {
      d_fragFPs.emplace_back(frag, d_fragFPPool[fragNum].get());
    }
    ++fragNum;
  }
  std::sort(d_fragFPs.begin(), d_fragFPs.end(),
            [](const auto &p1, const auto &p2) -> bool {
              return p1.first > p2.first;
            });
}

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceFingerprintSearcher::searchFragSet(
    const std::vector<std::unique_ptr<ROMol>> &fragSet,
    const SynthonSet &reaction) const {
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> results;

  // It can't be a hit if the number of fragments is more than the number
  // of synthon sets because some of the molecule won't be matched in any
  // of the potential products.  It can be less, in which case the unused
  // synthon set will be used completely, possibly resulting in a large
  // number of hits.
  if (fragSet.size() > reaction.getSynthons().size()) {
    return results;
  }

  std::vector<ExplicitBitVect *> fragFPs;
  fragFPs.reserve(fragSet.size());
  for (auto &frag : fragSet) {
    std::pair<void *, ExplicitBitVect *> tmp{frag.get(), nullptr};
    const auto it =
        std::lower_bound(d_fragFPs.begin(), d_fragFPs.end(), tmp,
                         [](const auto &p1, const auto &p2) -> bool {
                           return p1.first > p2.first;
                         });
    fragFPs.push_back(it->second);
  }

  const auto connPatterns = details::getConnectorPatterns(fragSet);
  const auto synthConnPatts = reaction.getSynthonConnectorPatterns();

  // Get all the possible permutations of connector numbers compatible with
  // the number of synthon sets in this reaction.  So if the
  // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
  // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
  // that might be how they're labelled in the reaction database.
  const auto connCombConnPatterns =
      details::getConnectorPermutations(connPatterns, reaction.getConnectors());

  // Need to try all combinations of synthon orders.
  const auto synthonOrders =
      details::permMFromN(fragSet.size(), reaction.getSynthons().size());

  for (const auto &synthonOrder : synthonOrders) {
    for (auto &connCombPatt : connCombConnPatterns) {
      // Make sure that for this connector combination, the synthons in this
      // order have something similar.  All query fragment connectors must
      // match something in the corresponding synthon.  The synthon can
      // have unused connectors.
      bool skip = false;
      for (size_t i = 0; i < connCombPatt.size(); ++i) {
        if ((connCombPatt[i] & synthConnPatts[synthonOrder[i]]).count() <
            connCombPatt[i].count()) {
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
        std::unique_ptr<SynthonSpaceHitSet> hs(
            new SynthonSpaceFPHitSet(reaction, theseSynthons, fragSet));
        if (hs->numHits) {
          results.push_back(std::move(hs));
        }
      }
    }
  }
  return results;
}

bool SynthonSpaceFingerprintSearcher::quickVerify(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  // The hitsets produced by the fingerprint searcher are SynthonSpaceFPHitSets,
  // which have the synthon fps as well.
  const auto hs = dynamic_cast<const SynthonSpaceFPHitSet *>(hitset);
  // Make an approximate fingerprint by combining the FPs for
  // these synthons, adding in the addFP and taking out the
  // subtractFP.
  ExplicitBitVect fullFP(*hs->synthonFPs[0][synthNums[0]]);
  for (unsigned int i = 1; i < synthNums.size(); ++i) {
    fullFP |= *hs->synthonFPs[i][synthNums[i]];
  }
  fullFP |= *hs->addFP;
  // The subtract FP has already had its bits flipped, so just do a
  // straight AND.
  fullFP &= *hs->subtractFP;

  return TanimotoSimilarity(fullFP, *d_queryFP) >=
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