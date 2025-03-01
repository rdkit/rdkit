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
#include <thread>
#include <boost/random/discrete_distribution.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>

namespace RDKit::SynthonSpaceSearch {

SynthonSpaceSearcher::SynthonSpaceSearcher(
    const ROMol &query, const SynthonSpaceSearchParams &params,
    SynthonSpace &space)
    : d_query(query), d_params(params), d_space(space) {
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
}

SearchResults SynthonSpaceSearcher::search() {
  std::vector<std::unique_ptr<ROMol>> results;

  const TimePoint *endTime = nullptr;
  TimePoint endTimePt;
  if (d_params.timeOut > 0) {
    endTimePt = Clock::now() + std::chrono::seconds(d_params.timeOut);
    endTime = &endTimePt;
  }
  bool timedOut = false;
  auto fragments =
      details::splitMolecule(d_query, getSpace().getMaxNumSynthons(),
                             d_params.maxNumFragSets, endTime, timedOut);
  if (timedOut || ControlCHandler::getGotSignal()) {
    return SearchResults{std::move(results), 0UL, timedOut,
                         ControlCHandler::getGotSignal()};
  }

  extraSearchSetup(fragments);
  if (ControlCHandler::getGotSignal()) {
    return SearchResults{std::move(results), 0UL, timedOut, true};
  }
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> allHits;
  size_t totHits = 0;
  int numTries = 100;
  for (const auto &id : getSpace().getReactionNames()) {
    timedOut = details::checkTimeOut(endTime);
    if (timedOut) {
      break;
    }
    const auto &reaction = getSpace().getReaction(id);
    for (auto &fragSet : fragments) {
      if (ControlCHandler::getGotSignal()) {
        break;
      }
      --numTries;
      if (!numTries) {
        timedOut = details::checkTimeOut(endTime);
        numTries = 100;
      }
      if (timedOut) {
        break;
      }
      if (auto theseHits = searchFragSet(fragSet, *reaction);
          !theseHits.empty()) {
        totHits += std::accumulate(
            theseHits.begin(), theseHits.end(), 0,
            [](const size_t prevVal,
               const std::unique_ptr<SynthonSpaceHitSet> &hs) -> size_t {
              return prevVal + hs->numHits;
            });
        allHits.insert(allHits.end(),
                       std::make_move_iterator(theseHits.begin()),
                       std::make_move_iterator(theseHits.end()));
      }
    }
  }
  if (!timedOut && !ControlCHandler::getGotSignal() && d_params.buildHits) {
    buildHits(allHits, endTime, timedOut, results);
  }

  return SearchResults{std::move(results), totHits, timedOut,
                       ControlCHandler::getGotSignal()};
}

std::unique_ptr<ROMol> SynthonSpaceSearcher::buildAndVerifyHit(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  std::vector<const ROMol *> synths(synthNums.size());
  std::vector<std::string> synthNames(synthNums.size());
  for (size_t i = 0; i < synthNums.size(); i++) {
    synths[i] = hitset->synthonsToUse[i][synthNums[i]].second;
    synthNames[i] = hitset->synthonsToUse[i][synthNums[i]].first;
  }

  std::unique_ptr<ROMol> prod;
  if (!quickVerify(hitset, synthNums)) {
    return prod;
  }
  prod = details::buildProduct(synths);

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
  if (prod) {
    const auto prodName =
        details::buildProductName(hitset->reactionId, synthNames);
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
                if (lsim == rsim) {
                  return lhs->getNumAtoms() < rhs->getNumAtoms();
                }
                return lsim > rsim;
              });
  }
}

void sortAndUniquifyToTry(
    std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
        &toTry) {
  std::vector<std::pair<size_t, std::string>> tmp;
  tmp.reserve(toTry.size());
  for (size_t i = 0; i < toTry.size(); i++) {
    tmp.emplace_back(
        i, details::buildProductName(toTry[i].first, toTry[i].second));
  }
  std::sort(tmp.begin(), tmp.end(),
            [](const auto &lhs, const auto &rhs) -> bool {
              return lhs.second < rhs.second;
            });
  tmp.erase(std::unique(tmp.begin(), tmp.end(),
                        [](const auto &lhs, const auto &rhs) -> bool {
                          return lhs.second == rhs.second;
                        }),
            tmp.end());
  std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
      newToTry;
  newToTry.reserve(tmp.size());
  std::transform(tmp.begin(), tmp.end(), back_inserter(newToTry),
                 [&](const auto &p) -> auto { return toTry[p.first]; });
  toTry = newToTry;
}

bool haveEnoughHits(const std::vector<std::unique_ptr<ROMol>> &results,
                    const std::int64_t maxHits, const std::int64_t hitStart) {
  const std::int64_t numHits = std::accumulate(
      results.begin(), results.end(), 0,
      [](const size_t prevVal, const std::unique_ptr<ROMol> &m) -> size_t {
        if (m) {
          return prevVal + 1;
        }
        return prevVal;
      });
  // If there's a limit on the number of hits, we still need to keep the
  // first hitStart hits and remove them later.  They had to be built
  // to see if they passed verifyHit.
  if (maxHits != -1 && numHits >= maxHits + hitStart) {
    return true;
  }
  return false;
}

}  // namespace

void SynthonSpaceSearcher::buildHits(
    std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
    const TimePoint *endTime, bool &timedOut,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }
  if (d_params.randomSample) {
    std::shuffle(hitsets.begin(), hitsets.end(), *d_randGen);
  } else {
    std::sort(hitsets.begin(), hitsets.end(),
              [](const auto &hs1, const auto &hs2) -> bool {
                if (hs1->reactionId == hs2->reactionId) {
                  return hs1->numHits < hs2->numHits;
                }
                return hs1->reactionId < hs2->reactionId;
              });
  }
  buildAllHits(hitsets, endTime, timedOut, results);
}

void SynthonSpaceSearcher::buildAllHits(
    const std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
    const TimePoint *endTime, bool &timedOut,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  // toTry is the synthon numbers for a particular set of synthons
  // in the SynthonSet that will be zipped together to form a possible
  // hit molecule.  It will possibly hold possibilities from multiple
  // SynthonSets, and when there are enough to be worth processing,
  // they will be built into molecules, verified and accepted or
  // rejected as hits.
  std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>> toTry;
  bool enoughHits = false;

  // Each hitset contains possible hits from a single SynthonSet.
  for (const auto &hitset : hitsets) {
    // Set up the stepper to move through the synthons.
    std::vector<size_t> numSynthons;
    numSynthons.reserve(hitset->synthonsToUse.size());
    for (auto &stu : hitset->synthonsToUse) {
      numSynthons.push_back(stu.size());
    }
    details::Stepper stepper(numSynthons);

    const auto &reaction = getSpace().getReaction(hitset->reactionId);
    std::vector<size_t> theseSynthNums(reaction->getSynthons().size(), 0);
    // process the synthons
    while (stepper.d_currState[0] != numSynthons[0]) {
      toTry.emplace_back(hitset.get(), stepper.d_currState);
      if (toTry.size() == static_cast<size_t>(d_params.toTryChunkSize)) {
        std::vector<std::unique_ptr<ROMol>> partResults;
        processToTrySet(toTry, endTime, partResults);
        results.insert(results.end(),
                       std::make_move_iterator(partResults.begin()),
                       std::make_move_iterator(partResults.end()));
        partResults.clear();
        enoughHits =
            haveEnoughHits(results, d_params.maxHits, d_params.hitStart);
        timedOut = details::checkTimeOut(endTime);
        toTry.clear();
        if (enoughHits || timedOut || ControlCHandler::getGotSignal()) {
          break;
        }
      }
      stepper.step();
    }
    if (enoughHits || timedOut || ControlCHandler::getGotSignal()) {
      break;
    }
  }

  // Do any remaining.
  if (!enoughHits && !timedOut && !toTry.empty()) {
    processToTrySet(toTry, endTime, results);
  }

  sortHits(results);
  // The multi-threaded versions might produce more hits than requested.
  if (d_params.maxHits != -1 && static_cast<std::int64_t>(results.size()) >
                                    d_params.maxHits + d_params.hitStart) {
    results.erase(results.begin() + d_params.maxHits + d_params.hitStart,
                  results.end());
  }

  // Now get rid of any hits before d_params.hitStart.  It seems wasteful to
  // do it like this, but until a hit has been through verifyHit there's
  // no way of knowing whether it should be kept plus the threading makes it
  // very complicated to do otherwise.
  if (d_params.hitStart) {
    if (d_params.hitStart < static_cast<std::int64_t>(results.size())) {
      std::for_each(results.begin(), results.begin() + d_params.hitStart,
                    [](std::unique_ptr<ROMol> &m) { m.reset(); });
      results.erase(std::remove_if(results.begin(), results.end(),
                                   [](const std::unique_ptr<ROMol> &r) {
                                     return !static_cast<bool>(r);
                                   }),
                    results.end());
    } else {
      results.clear();
    }
  }
}

namespace {
void processPartHitsFromDetails(
    const std::vector<
        std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>> &toTry,
    const TimePoint *endTime, std::vector<std::unique_ptr<ROMol>> &results,
    const SynthonSpaceSearcher *searcher, const size_t startNum,
    size_t finishNum) {
  std::uint64_t numTries = 100;
  finishNum = finishNum > toTry.size() ? toTry.size() : finishNum;
  for (size_t i = startNum; i < finishNum; ++i) {
    if (auto prod =
            searcher->buildAndVerifyHit(toTry[i].first, toTry[i].second)) {
      results[i] = std::move(prod);
      if (haveEnoughHits(results, searcher->getParams().maxHits,
                         searcher->getParams().hitStart)) {
        break;
      }
    }
    // Don't check the time every go, as it's quite expensive.
    --numTries;
    if (!numTries) {
      numTries = 100;
      if (details::checkTimeOut(endTime)) {
        break;
      }
    }
    if (ControlCHandler::getGotSignal()) {
      break;
    }
  }
}
}  // namespace

void SynthonSpaceSearcher::makeHitsFromToTry(
    const std::vector<
        std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>> &toTry,
    const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  results.resize(toTry.size());
  std::cout << "numThreads : " << d_params.numThreads << " goes to "
            << getNumThreadsToUse(d_params.numThreads) << std::endl;
#if RDK_BUILD_THREADSAFE_SSS
  if (const auto numThreads = getNumThreadsToUse(d_params.numThreads);
      numThreads > 1) {
    const size_t eachThread = 1 + toTry.size() / numThreads;
    std::cout << "building " << eachThread << " on each of " << numThreads
              << " threads\n";
    size_t start = 0;
    std::vector<std::thread> threads;
    for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
      threads.push_back(std::thread(processPartHitsFromDetails, std::ref(toTry),
                                    endTime, std::ref(results), this, start,
                                    start + eachThread));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    std::cout << " doing singlethread" << std::endl;
    processPartHitsFromDetails(toTry, endTime, results, this, 0, toTry.size());
  }
#else
  processPartHitsFromDetails(toTry, endTime, results, this, 0, toTry.size());
#endif

  // Take out any gaps in the results set, where products didn't make the grade.
  results.erase(std::remove_if(results.begin(), results.end(),
                               [](const std::unique_ptr<ROMol> &r) {
                                 return !static_cast<bool>(r);
                               }),
                results.end());
}

void SynthonSpaceSearcher::processToTrySet(
    std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
        &toTry,
    const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  // There are possibly duplicate entries in toTry, because 2
  // different fragmentations might produce overlapping synthon lists in
  // the same reaction. The duplicates need to be removed.  Although
  // when doing the job in batches this is less likely.
  sortAndUniquifyToTry(toTry);

  if (d_params.randomSample) {
    std::shuffle(toTry.begin(), toTry.end(), *d_randGen);
  }
  makeHitsFromToTry(toTry, endTime, results);
}

#if 0
void SynthonSpaceSearcher::buildHits(
    std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
    const size_t totHits, const TimePoint *endTime, bool &timedOut,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }
  std::sort(hitsets.begin(), hitsets.end(),
            [](const std::unique_ptr<SynthonSpaceHitSet> &hs1,
               const std::unique_ptr<SynthonSpaceHitSet> &hs2) -> bool {
              if (hs1->reactionId == hs2->reactionId) {
                return hs1->numHits < hs2->numHits;
              }
              return hs1->reactionId < hs2->reactionId;
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
    const std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
    std::set<std::string> &resultsNames, const TimePoint *endTime,
    bool &timedOut, std::vector<std::unique_ptr<ROMol>> &results) const {
  std::uint64_t numTries = 100;
  for (const auto &hitset : hitsets) {
    std::vector<size_t> numSynthons;
    numSynthons.reserve(hitset->synthonsToUse.size());
    for (auto &stu : hitset->synthonsToUse) {
      numSynthons.push_back(stu.size());
    }
    details::Stepper stepper(numSynthons);
    std::vector<ROMol *> theseSynths(hitset->synthonsToUse.size(), nullptr);
    while (stepper.d_currState[0] != numSynthons[0]) {
      if (auto prod = buildAndVerifyHit(hitset.get(), stepper.d_currState,
                                        resultsNames)) {
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
  RandomHitSelector(
      const std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
      const SynthonSpace &space)
      : d_hitsets(hitsets), d_synthSpace(space) {
    d_hitSetWeights.reserve(hitsets.size());
    std::transform(hitsets.begin(), hitsets.end(),
                   std::back_inserter(d_hitSetWeights),
                   [](const auto &hs) -> size_t { return hs->numHits; });
    d_hitSetSel = boost::random::discrete_distribution<size_t>(
        d_hitSetWeights.begin(), d_hitSetWeights.end());
    d_synthSels.resize(hitsets.size());
    for (size_t hi = 0; hi < hitsets.size(); ++hi) {
      d_synthSels[hi] =
          std::vector<boost::random::uniform_int_distribution<size_t>>(
              hitsets[hi]->synthonsToUse.size());
      for (size_t i = 0; i < hitsets[hi]->synthonsToUse.size(); ++i) {
        d_synthSels[hi][i] = boost::random::uniform_int_distribution<size_t>(
            0, hitsets[hi]->synthonsToUse[i].size() - 1);
      }
    }
  }

  std::pair<size_t, std::vector<size_t>> selectSynthComb(
      boost::random::mt19937 &randGen) const {
    const size_t hitSetNum = d_hitSetSel(randGen);
    std::vector<size_t> synths(d_hitsets[hitSetNum]->synthonsToUse.size());
    for (size_t i = 0; i < d_hitsets[hitSetNum]->synthonsToUse.size(); ++i) {
      const size_t synthNum = d_synthSels[hitSetNum][i](randGen);
      synths[i] = synthNum;
    }
    return std::make_pair(hitSetNum, synths);
  }

  const std::vector<std::unique_ptr<SynthonSpaceHitSet>> &d_hitsets;
  const SynthonSpace &d_synthSpace;

  std::vector<size_t> d_hitSetWeights;
  boost::random::discrete_distribution<size_t> d_hitSetSel;
  std::vector<std::vector<boost::random::uniform_int_distribution<size_t>>>
      d_synthSels;
};
}  // namespace

void SynthonSpaceSearcher::buildRandomHits(
    const std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
    const size_t totHits, std::set<std::string> &resultsNames,
    const TimePoint *endTime, bool &timedOut,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }
  const auto rhs = RandomHitSelector(hitsets, d_space);

  std::uint64_t numFails = 0;
  std::uint64_t numTries = 100;
  while (results.size() < std::min(static_cast<std::uint64_t>(d_params.maxHits),
                                   static_cast<std::uint64_t>(totHits)) &&
         numFails < totHits * d_params.numRandomSweeps) {
    const auto &[hitSetNum, synths] = rhs.selectSynthComb(*d_randGen);
    if (auto prod =
            buildAndVerifyHit(hitsets[hitSetNum].get(), synths, resultsNames)) {
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
#endif
}  // namespace RDKit::SynthonSpaceSearch
