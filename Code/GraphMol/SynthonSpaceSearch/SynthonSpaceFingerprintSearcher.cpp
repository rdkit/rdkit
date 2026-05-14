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
#include <mutex>
#include <ranges>

#include <DataStructs/BitOps.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceFingerprintSearcher.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>

namespace RDKit::SynthonSpaceSearch {

SynthonSpaceFingerprintSearcher::SynthonSpaceFingerprintSearcher(
    const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
    const SynthonSpaceSearchParams &params, SynthonSpace *space)
    : SynthonSpaceSearcher(query, params, space), d_fpGen(fpGen) {
  if (getSpace() && getSpace()->hasFingerprints() &&
      d_fpGen.infoString() != getSpace()->getSynthonFingerprintType()) {
    throw std::runtime_error(
        "The search fingerprints must match"
        " those in the database.  You are searching with " +
        d_fpGen.infoString() + " vs " +
        getSpace()->getSynthonFingerprintType() + " in the database.");
  }
  d_queryFP = std::unique_ptr<ExplicitBitVect>(d_fpGen.getFingerprint(query));
}

namespace {

// Find the first synthon in the set that has no less than the required
// number of set bits, based on the threshold.
size_t findSynthonSearchStart(unsigned int numFragSetBits,
                              double similarityCutoff, size_t synthonSetNum,
                              const SynthonSet &reaction) {
  unsigned int minBits = similarityCutoff * numFragSetBits;
  auto s = reaction.getFingerprintOrderedSynthon(synthonSetNum, 0);

  size_t first = 0;
  // This is the procedure that https://en.wikipedia.org/wiki/Binary_search
  // calls binary_search_leftmost.
  size_t last = reaction.getSynthons()[synthonSetNum].size();
  while (first < last) {
    size_t mid = first + (last - first) / 2;
    if (reaction.getFingerprintOrderedSynthon(synthonSetNum, mid)
            .second->getFP()
            ->getNumOnBits() < minBits) {
      first = mid + 1;
    } else {
      last = mid;
    }
  }
  return first;
}

// Take the fragged mol fps and flag all those synthons that have a fragment as
// a similarity match.
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<ExplicitBitVect *> &fragFPs,
    const double similarityCutoff, const SynthonSet &reaction,
    const std::vector<unsigned int> &synthonSetOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;
  std::vector<std::vector<std::pair<size_t, double>>> fragSims(
      reaction.getSynthons().size());

  synthonsToUse.reserve(reaction.getSynthons().size());
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  unsigned int maxBits = 0;
  for (size_t i = 0; i < synthonSetOrder.size(); i++) {
    if (similarityCutoff <= 0.0) {
      maxBits = fragFPs[i]->getNumBits();
    } else {
      maxBits = 1 + fragFPs[i]->getNumOnBits() / similarityCutoff;
    }
    auto start =
        findSynthonSearchStart(fragFPs[i]->getNumOnBits(), similarityCutoff,
                               synthonSetOrder[i], reaction);
    const auto &synthons = reaction.getSynthons()[synthonSetOrder[i]];
    bool fragMatched = false;
    for (size_t j = start; j < synthons.size(); j++) {
      // Search them in the sorted order, stopping when the number of bits
      // goes above maxBits.
      auto synthonNum =
          reaction.getFingerprintOrderedSynthonNum(synthonSetOrder[i], j);

      if (synthons[synthonNum].second->getFP()->getNumOnBits() > maxBits) {
        break;
      }
      if (const auto sim = TanimotoSimilarity(
              *fragFPs[i], *synthons[synthonNum].second->getFP());
          sim >= similarityCutoff) {
        synthonsToUse[synthonSetOrder[i]][synthonNum] = true;
        fragSims[synthonSetOrder[i]].emplace_back(synthonNum, sim);
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
      std::ranges::sort(fragSims[i], [](const auto &a, const auto &b) {
        return a.second > b.second;
      });
    }
    retSynthons[i].clear();
    std::ranges::transform(fragSims[i], std::back_inserter(retSynthons[i]),
                           [](const auto &fs) { return fs.first; });
  }
  return retSynthons;
}

}  // namespace

bool SynthonSpaceFingerprintSearcher::extraSearchSetup(
    std::vector<std::vector<std::shared_ptr<ROMol>>> &fragSets,
    const TimePoint *endTime) {
  if (!getSpace()->hasFingerprints() ||
      getSpace()->getSynthonFingerprintType() != d_fpGen.infoString()) {
    getSpace()->buildSynthonFingerprints(d_fpGen);
  }
  if (ControlCHandler::getGotSignal() || details::checkTimeOut(endTime)) {
    return false;
  }

  // Slightly convoluted way of doing it to prepare for multi-threading.
  // Make a map of the unique SMILES strings for the fragments, keeping
  // track of them in the vector.
  bool cancelled = false;
  auto fragSmiToFrag = details::mapFragsBySmiles(fragSets, cancelled);
  if (cancelled) {
    return false;
  }

  // Generate the fingerprints for the fragments.
  d_fragFPPool.resize(fragSmiToFrag.size());
  unsigned int fragNum = 0;
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    if (ControlCHandler::getGotSignal()) {
      return false;
    }
    d_fragFPPool[fragNum++].reset(d_fpGen.getFingerprint(*frags.front()));
  }

  // Use the pooled fps to populate the vectors for each fragSet.
  fragNum = 0;
  d_fragFPs.reserve(fragSmiToFrag.size());
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    for (auto &frag : frags) {
      d_fragFPs.emplace_back(frag, d_fragFPPool[fragNum].get());
    }
    ++fragNum;
  }
  std::ranges::sort(d_fragFPs, [](const auto &p1, const auto &p2) -> bool {
    return p1.first > p2.first;
  });

  return true;
}

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceFingerprintSearcher::searchFragSet(
    const std::vector<std::shared_ptr<ROMol>> &fragSet,
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
    const auto it = std::ranges::lower_bound(
        d_fragFPs, tmp, [](const auto &p1, const auto &p2) -> bool {
          return p1.first > p2.first;
        });
    fragFPs.push_back(it->second);
  }

  // Need to try all combinations of synthon orders.
  const auto synthonOrders =
      details::permMFromN(fragSet.size(), reaction.getSynthons().size());

  for (const auto &synthonOrder : synthonOrders) {
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
  return results;
}

bool SynthonSpaceFingerprintSearcher::quickVerify(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  if (!SynthonSpaceSearcher::quickVerify(hitset, synthNums)) {
    return false;
  }
  auto approxSim = approxSimilarity(hitset, synthNums);
  return approxSim >=
         getParams().similarityCutoff - getParams().approxSimilarityAdjuster;
}

double SynthonSpaceFingerprintSearcher::approxSimilarity(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  // The hitsets produced by the fingerprint searcher are SynthonSpaceFPHitSets,
  // which have the synthon fps as well.
  const auto hs = dynamic_cast<const SynthonSpaceFPHitSet *>(hitset);
  // Make an approximate fingerprint by combining the FPs for
  // these synthons.
  ExplicitBitVect fullFP(*hs->synthonFPs[0][synthNums[0]]);
  for (unsigned int i = 1; i < synthNums.size(); ++i) {
    fullFP |= *hs->synthonFPs[i][synthNums[i]];
  }
  fullFP |= *hs->addFP;
  // The subtract FP has already had its bits flipped, so just do a
  // straight AND.
  fullFP &= *hs->subtractFP;
  return TanimotoSimilarity(fullFP, *d_queryFP);
}

bool SynthonSpaceFingerprintSearcher::verifyHit(
    ROMol &hit, const std::string &rxnId,
    const std::vector<const std::string *> &synthNames) {
  if (!SynthonSpaceSearcher::verifyHit(hit, rxnId, synthNames)) {
    return false;
  }
  const std::unique_ptr<ExplicitBitVect> fp(d_fpGen.getFingerprint(hit));
  const auto sim = TanimotoSimilarity(*fp, *d_queryFP);
  if (sim > getBestSimilaritySoFar()) {
    const auto prodName = details::buildProductName(rxnId, synthNames);
    hit.setProp<std::string>(common_properties::_Name, prodName);
    updateBestHitSoFar(hit, sim);
  }
  if (sim >= getParams().similarityCutoff) {
    hit.setProp<double>("Similarity", sim);
    const auto prodName = details::buildProductName(rxnId, synthNames);
    hit.setProp<std::string>(common_properties::_Name, prodName);
    return true;
  }
  return false;
}

}  // namespace RDKit::SynthonSpaceSearch