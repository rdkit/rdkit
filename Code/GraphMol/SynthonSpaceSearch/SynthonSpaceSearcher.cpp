//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <random>
#include <boost/random/discrete_distribution.hpp>

#include <RDGeneral/ControlCHandler.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>

namespace RDKit::SynthonSpaceSearch {

SynthonSpaceSearcher::SynthonSpaceSearcher(
    const ROMol &query, const SynthonSpaceSearchParams &params,
    SynthonSpace &space)
    : d_query(query), d_params(params), d_space(space) {}

SearchResults SynthonSpaceSearcher::search() {
  if (d_params.randomSample && d_params.maxHits == -1) {
    throw std::runtime_error(
        "Random sample is incompatible with maxHits of -1.");
  }

  if (d_params.randomSample) {
    if (!d_randGen) {
      // Use boost random number code because it is the same across
      // different platforms, which isn't crucial but makes the tests
      // work on all platforms.
      if (d_params.randomSeed == -1) {
        std::random_device rd;
        d_randGen = std::make_unique<boost::mt19937>(rd());
      } else {
        d_randGen = std::make_unique<boost::mt19937>(d_params.randomSeed);
      }
    }
  }
  std::vector<std::unique_ptr<ROMol>> results;

  TimePoint *endTime = nullptr;
  TimePoint endTimePt;
  if (d_params.timeOut > 0) {
    endTimePt = Clock::now() + std::chrono::seconds(d_params.timeOut);
    endTime = &endTimePt;
  }
  ControlCHandler::reset();
  bool timedOut = false;
  auto fragments = details::splitMolecule(
      d_query, d_params.maxBondSplits, d_params.maxNumFrags, endTime, timedOut);
  if (timedOut || ControlCHandler::getGotSignal()) {
    return SearchResults{std::move(results), 0UL, timedOut,
                         ControlCHandler::getGotSignal()};
  }

  std::vector<SynthonSpaceHitSet> allHits;
  size_t totHits = 0;
  for (auto &fragSet : fragments) {
    timedOut = details::checkTimeOut(endTime);
    if (timedOut) {
      break;
    }
    if (ControlCHandler::getGotSignal()) {
      break;
    }
    if (auto theseHits = searchFragSet(fragSet); !theseHits.empty()) {
      totHits += std::accumulate(
          theseHits.begin(), theseHits.end(), 0,
          [](const size_t prevVal, const SynthonSpaceHitSet &hs) -> size_t {
            return prevVal + hs.numHits;
          });
      allHits.insert(allHits.end(), theseHits.begin(), theseHits.end());
    }
  }

  if (!timedOut && !ControlCHandler::getGotSignal() && d_params.buildHits) {
    buildHits(allHits, totHits, endTime, timedOut, results);
  }

  return SearchResults{std::move(results), totHits, timedOut,
                       ControlCHandler::getGotSignal()};
}

std::unique_ptr<ROMol> SynthonSpaceSearcher::buildAndVerifyHit(
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<size_t> &synthNums,
    std::set<std::string> &resultsNames) const {
  const auto prodName = reaction->buildProductName(synthNums);

  std::unique_ptr<ROMol> prod;
  if (resultsNames.insert(prodName).second) {
    if (resultsNames.size() < static_cast<size_t>(d_params.hitStart)) {
      return prod;
    }
    if (!quickVerify(reaction, synthNums)) {
      return prod;
    }
    prod = reaction->buildProduct(synthNums);

    // Do a final check of the whole thing.  It can happen that the
    // fragments match synthons but the final product doesn't match.
    // A key example is when the 2 synthons come together to form an
    // aromatic ring.  For substructure searching, an aliphatic query
    // can match the aliphatic synthon so they are selected as a hit,
    // but the final aromatic ring isn't a match.
    // E.g. Cc1cccc(C(=O)N[1*])c1N=[2*] and c1ccoc1C(=[2*])[1*]
    // making Cc1cccc2c(=O)[nH]c(-c3ccco3)nc12.  The query c1ccc(CN)o1
    // when split is a match to the synthons (c1ccc(C[1*])o1 and [1*]N)
    // but the product the hydroxyquinazoline is aromatic, at least in
    // the RDKit model so the N in the query doesn't match.
    if (!verifyHit(*prod)) {
      prod.reset();
    }
  }
  if (prod) {
    prod->setProp<std::string>(common_properties::_Name, prodName);
  }
  return prod;
}

namespace {
void sortHits(std::vector<std::unique_ptr<ROMol>> &hits) {
  if (!hits.empty() && hits.front()->hasProp("Similarity")) {
    std::sort(hits.begin(), hits.end(),
              [](const std::unique_ptr<ROMol> &lhs,
                 const std::unique_ptr<ROMol> &rhs) {
                const auto lsim = lhs->getProp<double>("Similarity");
                const auto rsim = rhs->getProp<double>("Similarity");
                return lsim > rsim;
              });
  }
}
}  // namespace

void SynthonSpaceSearcher::buildHits(
    std::vector<SynthonSpaceHitSet> &hitsets, const size_t totHits,
    const TimePoint *endTime, bool &timedOut,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }
  std::sort(
      hitsets.begin(), hitsets.end(),
      [](const SynthonSpaceHitSet &hs1, const SynthonSpaceHitSet &hs2) -> bool {
        if (hs1.reactionId == hs2.reactionId) {
          return hs1.numHits < hs2.numHits;
        }
        return hs1.reactionId < hs2.reactionId;
      });
  // Keep track of the result names so we can weed out duplicates by
  // reaction and synthons.  Different splits may give rise to the same
  // synthon combination.  This will keep the same molecule produced via
  // different reactions which I think makes sense.  The resultsNames will
  // be accumulated even if the molecule itself doesn't make it into the
  // results set, for example if it isn't a random selection or it's
  // outside maxHits or hitStart.
  std::set<std::string> resultsNames;

  if (d_params.randomSample) {
    buildRandomHits(hitsets, totHits, resultsNames, endTime, timedOut, results);
  } else {
    buildAllHits(hitsets, resultsNames, endTime, timedOut, results);
  }
  sortHits(results);
}

void SynthonSpaceSearcher::buildAllHits(
    const std::vector<SynthonSpaceHitSet> &hitsets,
    std::set<std::string> &resultsNames, const TimePoint *endTime,
    bool &timedOut, std::vector<std::unique_ptr<ROMol>> &results) const {
  std::uint64_t numTries = 100;
  for (const auto &[reactionId, synthonsToUse, numHits] : hitsets) {
    std::vector<size_t> numSynthons;
    numSynthons.reserve(synthonsToUse.size());
    for (auto &stu : synthonsToUse) {
      numSynthons.push_back(stu.size());
    }
    const auto &reaction = getSpace().getReactions().find(reactionId)->second;
    details::Stepper stepper(numSynthons);
    std::vector<size_t> theseSynthNums(synthonsToUse.size(), 0);
    while (stepper.d_currState[0] != numSynthons[0]) {
      for (size_t i = 0; i < stepper.d_currState.size(); ++i) {
        theseSynthNums[i] = synthonsToUse[i][stepper.d_currState[i]];
      }
      if (auto prod =
              buildAndVerifyHit(reaction, theseSynthNums, resultsNames)) {
        results.push_back(std::move(prod));
      }
      if (results.size() == static_cast<size_t>(d_params.maxHits)) {
        return;
      }
      stepper.step();
      // Don't check the time every go, as it's quite expensive.
      --numTries;
      if (!numTries) {
        numTries = 100;
        timedOut = details::checkTimeOut(endTime);
        if (timedOut) {
          break;
        }
      }
      if (ControlCHandler::getGotSignal()) {
        break;
      }
    }
    if (timedOut || ControlCHandler::getGotSignal()) {
      break;
    }
  }
}

namespace {
struct RandomHitSelector {
  // Uses a weighted random number selector to give a random representation
  // of the hits from each reaction proportionate to the total number of hits
  // expected from each reaction.
  RandomHitSelector(const std::vector<SynthonSpaceHitSet> &hitsets,
                    const SynthonSpace &space)
      : d_hitsets(hitsets), d_synthSpace(space) {
    d_hitSetWeights.reserve(hitsets.size());
    std::transform(
        hitsets.begin(), hitsets.end(), std::back_inserter(d_hitSetWeights),
        [](const SynthonSpaceHitSet &hs) -> size_t { return hs.numHits; });
    d_hitSetSel = boost::random::discrete_distribution<size_t>(
        d_hitSetWeights.begin(), d_hitSetWeights.end());
    d_synthSels.resize(hitsets.size());
    for (size_t hi = 0; hi < hitsets.size(); ++hi) {
      d_synthSels[hi] =
          std::vector<boost::random::uniform_int_distribution<size_t>>(
              hitsets[hi].synthonsToUse.size());
      for (size_t i = 0; i < hitsets[hi].synthonsToUse.size(); ++i) {
        d_synthSels[hi][i] = boost::random::uniform_int_distribution<size_t>(
            0, hitsets[hi].synthonsToUse[i].size() - 1);
      }
    }
  }

  std::pair<std::string, std::vector<size_t>> selectSynthComb(
      boost::random::mt19937 &randGen) {
    std::vector<size_t> synths;
    const size_t hitSetNum = d_hitSetSel(randGen);
    for (size_t i = 0; i < d_hitsets[hitSetNum].synthonsToUse.size(); ++i) {
      const size_t synthNum = d_synthSels[hitSetNum][i](randGen);
      synths.push_back(d_hitsets[hitSetNum].synthonsToUse[i][synthNum]);
    }
    return std::make_pair(d_hitsets[hitSetNum].reactionId, synths);
  }

  const std::vector<SynthonSpaceHitSet> &d_hitsets;
  const SynthonSpace &d_synthSpace;

  std::vector<size_t> d_hitSetWeights;
  boost::random::discrete_distribution<size_t> d_hitSetSel;
  std::vector<std::vector<boost::random::uniform_int_distribution<size_t>>>
      d_synthSels;
};
}  // namespace

void SynthonSpaceSearcher::buildRandomHits(
    const std::vector<SynthonSpaceHitSet> &hitsets, const size_t totHits,
    std::set<std::string> &resultsNames, const TimePoint *endTime,
    bool &timedOut, std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }
  auto rhs = RandomHitSelector(hitsets, d_space);

  std::uint64_t numFails = 0;
  std::uint64_t numTries = 100;
  while (results.size() < std::min(static_cast<std::uint64_t>(d_params.maxHits),
                                   static_cast<std::uint64_t>(totHits)) &&
         numFails < totHits * d_params.numRandomSweeps) {
    const auto &[reactionId, synths] = rhs.selectSynthComb(*d_randGen);
    const auto &reaction = getSpace().getReactions().find(reactionId)->second;
    if (auto prod = buildAndVerifyHit(reaction, synths, resultsNames)) {
      results.push_back(std::move(prod));
    } else {
      numFails++;
    }
    // Don't check the time every go, as it's quite expensive.
    --numTries;
    if (!numTries) {
      numTries = 100;
      timedOut = details::checkTimeOut(endTime);
      if (timedOut) {
        break;
      }
    }
    if (ControlCHandler::getGotSignal()) {
      break;
    }
  }
}

}  // namespace RDKit::SynthonSpaceSearch
