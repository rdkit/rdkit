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
        d_randGen = std::make_unique<std::mt19937>(rd());
      } else {
        d_randGen = std::make_unique<std::mt19937>(d_params.randomSeed);
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
  auto fragments = details::splitMolecule(
      d_query, getSpace().getMaxNumSynthons(), d_params.maxNumFragSets, endTime,
      d_params.numThreads, timedOut);
  if (timedOut || ControlCHandler::getGotSignal()) {
    return SearchResults{std::move(results), 0ULL, timedOut,
                         ControlCHandler::getGotSignal()};
  }
  extraSearchSetup(fragments);
  if (ControlCHandler::getGotSignal()) {
    return SearchResults{std::move(results), 0ULL, timedOut, true};
  }

  std::uint64_t totHits = 0;
  auto allHits = doTheSearch(fragments, endTime, timedOut, totHits);
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
        details::buildProductName(hitset->d_reaction->getId(), synthNames);
    prod->setProp<std::string>(common_properties::_Name, prodName);
  }
  return prod;
}

namespace {
std::vector<std::unique_ptr<SynthonSpaceHitSet>> searchReaction(
    SynthonSpaceSearcher *searcher, const SynthonSet &reaction,
    const TimePoint *endTime,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragments) {
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> hits;

  int numTries = 100;
  bool timedOut = false;
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
    if (auto theseHits = searcher->searchFragSet(fragSet, reaction);
        !theseHits.empty()) {
      hits.insert(hits.end(), std::make_move_iterator(theseHits.begin()),
                  std::make_move_iterator(theseHits.end()));
    }
  }

  return hits;
}

void processReactions(
    SynthonSpaceSearcher *searcher,
    const std::vector<std::string> &reactionNames,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragments,
    const TimePoint *endTime, std::atomic<std::int64_t> &mostRecentReaction,
    std::int64_t lastReaction,
    std::vector<std::vector<std::unique_ptr<SynthonSpaceHitSet>>>
        &reactionHits) {
  bool timedOut = false;

  while (true) {
    std::int64_t thisR = ++mostRecentReaction;
    // std::cout << thisR << " vs " << lastReaction << std::endl;
    if (thisR > lastReaction) {
      break;
    }
    timedOut = details::checkTimeOut(endTime);
    if (timedOut) {
      break;
    }
    const auto &reaction =
        searcher->getSpace().getReaction(reactionNames[thisR]);
    auto theseHits = searchReaction(searcher, *reaction, endTime, fragments);
    reactionHits[thisR] = std::move(theseHits);
    timedOut = details::checkTimeOut(endTime);
    if (timedOut) {
      break;
    }
  }
}
}  // namespace

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceSearcher::doTheSearch(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets,
    const TimePoint *endTime, bool &timedOut, std::uint64_t &totHits) {
  auto reactionNames = getSpace().getReactionNames();
  std::vector<std::vector<std::unique_ptr<SynthonSpaceHitSet>>> reactionHits(
      reactionNames.size());

  std::int64_t lastReaction = reactionNames.size() - 1;
  std::atomic<std::int64_t> mostRecentReaction = -1;
#if RDK_BUILD_THREADSAFE_SSS
  if (const auto numThreads = getNumThreadsToUse(d_params.numThreads);
      numThreads > 1) {
    std::vector<std::thread> threads;
    for (unsigned int i = 0U;
         i < std::min(static_cast<size_t>(numThreads), reactionNames.size());
         ++i) {
      threads.push_back(std::thread(processReactions, this,
                                    std::ref(reactionNames), std::ref(fragSets),
                                    endTime, std::ref(mostRecentReaction),
                                    lastReaction, std::ref(reactionHits)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    processReactions(this, reactionNames, fragSets, endTime, mostRecentReaction,
                     lastReaction, reactionHits);
  }
#else
  processReactions(this, reactionNames, fragSets, endTime, mostRecentReaction,
                   lastReaction, reactionHits);
#endif

  std::vector<std::unique_ptr<SynthonSpaceHitSet>> allHits;
  totHits = 0;
  for (size_t i = 0; i < reactionHits.size(); ++i) {
    totHits += std::accumulate(
        reactionHits[i].begin(), reactionHits[i].end(), 0,
        [](const size_t prevVal, const std::unique_ptr<SynthonSpaceHitSet> &hs)
            -> size_t { return prevVal + hs->numHits; });
    allHits.insert(allHits.end(),
                   std::make_move_iterator(reactionHits[i].begin()),
                   std::make_move_iterator(reactionHits[i].end()));
  }
  timedOut = details::checkTimeOut(endTime);
  return allHits;
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
                if (hs1->d_reaction->getId() == hs2->d_reaction->getId()) {
                  return hs1->numHits < hs2->numHits;
                }
                return hs1->d_reaction->getId() < hs2->d_reaction->getId();
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

    const auto &reaction = getSpace().getReaction(hitset->d_reaction->getId());
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
    const SynthonSpaceSearcher *searcher,
    std::atomic<std::int64_t> &mostRecentTry, std::int64_t lastTry) {
  std::uint64_t numTries = 100;
  while (true) {
    std::int64_t thisTry = ++mostRecentTry;
    if (thisTry > lastTry) {
      break;
    }
    if (auto prod = searcher->buildAndVerifyHit(toTry[thisTry].first,
                                                toTry[thisTry].second)) {
      results[thisTry] = std::move(prod);
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
  std::int64_t lastTry = toTry.size() - 1;
  std::atomic<std::int64_t> mostRecentTry = -1;

#if RDK_BUILD_THREADSAFE_SSS
  // This assumes that each chunk of the toTry list will take roughly the
  // same amount of time to process.  To a first approximation, that's
  // probably reasonable.  Some entries in toTry will fail quickVerify
  // so the more time-consuming construction of the hit and final
  // check in verifyHit won't be needed, so the chunks won't take exactly
  // equal time. It can always be re-visited if the threads run for very
  // different lengths of time in an average search.
  if (const auto numThreads = getNumThreadsToUse(d_params.numThreads);
      numThreads > 1) {
    const size_t eachThread = 1 + toTry.size() / numThreads;
    size_t start = 0;
    std::vector<std::thread> threads;
    for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
      threads.push_back(std::thread(processPartHitsFromDetails, std::ref(toTry),
                                    endTime, std::ref(results), this,
                                    std::ref(mostRecentTry), lastTry));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    processPartHitsFromDetails(toTry, endTime, results, this, mostRecentTry,
                               lastTry);
  }
#else
  processPartHitsFromDetails(toTry, endTime, results, this, mostRecentTry,
                             lastTry);
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
}  // namespace RDKit::SynthonSpaceSearch
