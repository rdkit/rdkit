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
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceRascalSearcher.h>

namespace RDKit::SynthonSpaceSearch {
SynthonSpaceRascalSearcher::SynthonSpaceRascalSearcher(
    const ROMol &query, const RascalMCES::RascalOptions &rascalOptions,
    const SynthonSpaceSearchParams &params, SynthonSpace *space)
    : SynthonSpaceSearcher(query, params, space),
      d_rascalOptions(rascalOptions),
      d_rascalFragOptions(rascalOptions) {
  d_rascalFragOptions.similarityThreshold -= params.fragSimilarityAdjuster;
}

namespace {

size_t findSynthonSearchStart(unsigned int fragAtomsAndBonds,
                              double similarityCutoff, size_t synthonSetNum,
                              const SynthonSet &reaction) {
  unsigned int minAtomsAndBonds = similarityCutoff * fragAtomsAndBonds;
  auto s = reaction.getRascalOrderedSynthon(synthonSetNum, 0);

  size_t first = 0;
  size_t last = reaction.getSynthons()[synthonSetNum].size();
  while (first < last) {
    size_t mid = first + (last - first) / 2;
    const auto &synthSM = reaction.getRascalOrderedSynthon(synthonSetNum, mid)
                              .second->getSearchMol();
    if (synthSM->getNumAtoms() + synthSM->getNumBonds() < minAtomsAndBonds) {
      first = mid + 1;
    } else {
      last = mid;
    }
  }
  return first;
}

std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<std::shared_ptr<ROMol>> &fragSet,
    const RascalMCES::RascalOptions &rascalOptions, const SynthonSet &reaction,
    const std::vector<unsigned int> &synthonSetOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;

  // It makes sense to match fragments against synthon sets in order of
  // smallest synthon set first because if a fragment doesn't have a match
  // in a synthon set, the whole thing's a bust.  So if fragShapes[0] is matched
  // against 1000 synthons and then fragShapes[1] is matched against 10 synthons
  // and doesn't match any of them, the first set of matches was wasted time.
  std::vector<std::pair<unsigned int, size_t>> fragOrders(
      synthonSetOrder.size());
  for (size_t i = 0; i < synthonSetOrder.size(); i++) {
    fragOrders[i].first = i;
    fragOrders[i].second = reaction.getSynthons()[synthonSetOrder[i]].size();
  }
  std::ranges::sort(fragOrders, [](const auto &a, const auto &b) {
    return a.second < b.second;
  });

  synthonsToUse.reserve(reaction.getSynthons().size());
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  for (size_t i = 0; i < synthonSetOrder.size(); i++) {
    const auto fragNum = fragOrders[i].first;
    unsigned int fragAtomsAndBonds =
        fragSet[fragNum]->getNumAtoms() + fragSet[fragNum]->getNumBonds();
    unsigned int maxAtomsAndBonds =
        1 + fragAtomsAndBonds / rascalOptions.similarityThreshold;
    auto start = findSynthonSearchStart(fragAtomsAndBonds,
                                        rascalOptions.similarityThreshold,
                                        synthonSetOrder[fragNum], reaction);
    const auto &synthons = reaction.getSynthons()[synthonSetOrder[fragNum]];
    bool fragMatched = false;
    for (size_t j = start; j < synthons.size(); j++) {
      // Search them in the sorted order, stopping when the number of atoms and
      // bonds goes above maxAtomsAndBonds.
      auto synthonNum =
          reaction.getRascalOrderedSynthonNum(synthonSetOrder[fragNum], j);
      const auto &synthonSM = synthons[synthonNum].second->getSearchMol();
      if (synthonSM->getNumAtoms() + synthonSM->getNumBonds() >
          maxAtomsAndBonds) {
        break;
      }
      const auto rascalResults =
          RascalMCES::rascalMCES(*fragSet[fragNum], *synthonSM, rascalOptions);
      if (!rascalResults.empty()) {
        synthonsToUse[synthonSetOrder[fragNum]][synthonNum] = true;
        fragMatched = true;
      }
    }
    if (!fragMatched) {
      return retSynthons;
    }
  }
  // Fill in any synthons where they all didn't match because there were
  // fewer fragments than synthons.
  details::expandBitSet(synthonsToUse);
  details::bitSetsToVectors(synthonsToUse, retSynthons);

  // Now sort the selected synthons into ascending order of number of
  // atoms, since smaller molecules are likely to be of more interest.
  for (size_t i = 0; i < retSynthons.size(); ++i) {
    const auto &synthonsi = reaction.getSynthons()[i];
    std::sort(retSynthons[i].begin(), retSynthons[i].end(),
              [&](const size_t a, const size_t b) {
                return synthonsi[a].second->getOrigMol()->getNumAtoms() <
                       synthonsi[b].second->getOrigMol()->getNumAtoms();
              });
  }
  return retSynthons;
}
}  // namespace

bool SynthonSpaceRascalSearcher::extraSearchSetup(
    std::vector<std::vector<std::shared_ptr<ROMol>>> &fragSets,
    const TimePoint *endTime) {
  int numDone = 100;
  for (const auto &fragSet : fragSets) {
    for (const auto &frag : fragSet) {
      unsigned int otf;
      sanitizeMol(*static_cast<RWMol *>(frag.get()), otf,
                  MolOps::SANITIZE_SYMMRINGS);
      --numDone;
      if (!numDone) {
        numDone = 100;
        if (details::checkTimeOut(endTime)) {
          return false;
        }
      }
    }
  }
  return true;
}

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceRascalSearcher::searchFragSet(
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

  const auto connPatterns = details::getConnectorPatterns(fragSet);
  boost::dynamic_bitset<> conns(MAX_CONNECTOR_NUM + 1);
  for (auto &connPattern : connPatterns) {
    conns |= connPattern;
  }

  auto synthConnPatts = reaction.getSynthonConnectorPatterns();

  // Get all the possible permutations of connector numbers compatible with
  // the number of synthon sets in this reaction.  So if the
  // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
  // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
  // that might be how they're labelled in the reaction database.
  const auto connCombConnPatterns =
      details::getConnectorPermutations(connPatterns, reaction.getConnectors());

  // Need to try all combinations of synthon orders.
  auto synthonOrders =
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

      // Rascal ignores isotope numbers, which makes things easier.
      auto theseSynthons =
          getHitSynthons(fragSet, d_rascalFragOptions, reaction, synthonOrder);
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

bool SynthonSpaceRascalSearcher::quickVerify(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  if (!SynthonSpaceSearcher::quickVerify(hitset, synthNums)) {
    return false;
  }
  auto approxSim = approxSimilarity(hitset, synthNums);
  return approxSim >= d_rascalOptions.similarityThreshold;
}

double SynthonSpaceRascalSearcher::approxSimilarity(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  // If the query is an exact substructure of the product, then that's an upper
  // bound on the Johnson similarity.  Check that that is not below the
  // threshold.
  int qbit = getQuery().getNumAtoms() + getQuery().getNumBonds();
  int numAtoms = 0, numBonds = 0, numConns = 0;
  for (size_t i = 0; i < synthNums.size(); i++) {
    const auto &synth = hitset->synthonsToUse[i][synthNums[i]].second;
    numAtoms += synth->getNumHeavyAtoms();
    numBonds += synth->getOrigMol()->getNumBonds() - synth->getNumDummies();
    numConns += synth->getNumDummies();
  }
  // Pairs of dummies will form a bond in the final result, which need to
  // be accounted for.  If there are hanging dummy atoms the result will
  // be an over-estimate of the similarity because we won't be adding enough
  // bonds to the denominator.  That's ok for this.
  numConns /= 2;
  // The Johnson similarity is
  // (commonNatoms + commonNbonds)**2 /
  // ((Natoms1 + Nbonds1) * (Natoms2 + Natoms2))
  // and in this case the common atoms are the whole query, so the square
  // cancels out.
  double bestSim = static_cast<double>(qbit) /
                   static_cast<double>(numAtoms + numBonds + numConns);
  return bestSim;
}

bool SynthonSpaceRascalSearcher::verifyHit(
    ROMol &hit, const std::string &rxnId,
    const std::vector<const std::string *> &synthNames) {
  if (!SynthonSpaceSearcher::verifyHit(hit, rxnId, synthNames)) {
    return false;
  }
  auto res = RascalMCES::rascalMCES(hit, getQuery(), d_rascalOptions);
  // Rascal reports all matches that proceed to full MCES elucidation,
  // even if the final similarity value ends up below the threshold.
  // We only want those over the threshold.
  if (!res.empty()) {
    if (res.front().getSimilarity() > getBestSimilaritySoFar()) {
      const auto prodName = details::buildProductName(rxnId, synthNames);
      hit.setProp<std::string>(common_properties::_Name, prodName);
      updateBestHitSoFar(hit, res.front().getSimilarity());
    }
    if (res.front().getSimilarity() >= d_rascalOptions.similarityThreshold) {
      hit.setProp<double>("Similarity", res.front().getSimilarity());
      const auto prodName = details::buildProductName(rxnId, synthNames);
      hit.setProp<std::string>(common_properties::_Name, prodName);
      return true;
    }
  }
  return false;
}

}  // namespace RDKit::SynthonSpaceSearch
