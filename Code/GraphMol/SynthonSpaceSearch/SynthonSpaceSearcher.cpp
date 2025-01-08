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
#include <atomic>
#include <random>
#include <thread>
#include <boost/random/discrete_distribution.hpp>

#include <RDGeneral/RDThreads.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>
#include <boost/fusion/container/vector/vector.hpp>

namespace RDKit::SynthonSpaceSearch {

static std::atomic<bool> timedOut = false;

namespace {
void checkTimeOut(const TimePoint *endTime) {
  if (!timedOut && endTime != nullptr && Clock::now() > *endTime) {
    BOOST_LOG(rdWarningLog) << "Timed out.\n";
    timedOut = true;
  }
}
}  // namespace

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

  auto fragments = details::splitMolecule(d_query, d_params.maxBondSplits);
  TimePoint *endTime = nullptr;
  TimePoint endTimePt;
  timedOut = false;
  if (d_params.timeOut > 0) {
    endTimePt = Clock::now() + std::chrono::seconds(d_params.timeOut);
    endTime = &endTimePt;
  }
  std::vector<std::vector<SynthonSpaceHitSet>> allFragHits(fragments.size());
#if RDK_BUILD_THREADSAFE_SSS
  auto processPartFragSet =
      [](std::vector<std::vector<std::unique_ptr<ROMol>>> &fragments,
         TimePoint *endTime,
         std::vector<std::vector<SynthonSpaceHitSet>> &allFragHits,
         SynthonSpaceSearcher *self, size_t start, size_t finish) -> void {
    finish = finish > fragments.size() ? fragments.size() : finish;
    for (size_t i = start; i < finish; i++) {
      checkTimeOut(endTime);
      if (timedOut) {
        break;
      }
      allFragHits[i] = self->searchFragSet(fragments[i]);
    }
  };

  auto numThreads = getNumThreadsToUse(d_params.numThreads);
  std::cout << "numThreads: " << numThreads << std::endl;
  if (numThreads > 1) {
    size_t eachThread = 1 + (fragments.size() / numThreads);
    size_t start = 0;
    std::vector<std::thread> threads;
    for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
      threads.push_back(std::thread(processPartFragSet, std::ref(fragments),
                                    endTime, std::ref(allFragHits), this, start,
                                    start + eachThread));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    processPartFragSet(fragments, endTime, allFragHits, this, 0,
                       fragments.size());
  }
#else
  processPartFragSet(fragments, endTime, allFragHits, this, 0,
                     fragments.size());
#endif

  size_t totHits = 0;
  std::vector<SynthonSpaceHitSet> allHits;
  for (const auto &fh : allFragHits) {
    totHits += std::accumulate(
        fh.begin(), fh.end(), 0,
        [](const size_t prevVal, const SynthonSpaceHitSet &hs) -> size_t {
          return prevVal + hs.numHits;
        });
    allHits.insert(allHits.end(), fh.begin(), fh.end());
  }

  if (!timedOut && d_params.buildHits) {
    buildHits(allHits, totHits, endTime, results);
  }

  return SearchResults{std::move(results), totHits, timedOut};
}

std::unique_ptr<ROMol> SynthonSpaceSearcher::buildAndVerifyHit(
    const SynthonSet *reaction, const std::vector<size_t> &synthNums) const {
  std::unique_ptr<ROMol> prod;
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
  if (prod) {
    const auto prodName = reaction->buildProductName(synthNums);
    prod->setProp<std::string>(common_properties::_Name, prodName);
  }
  return prod;
}

void SynthonSpaceSearcher::buildHits(
    std::vector<SynthonSpaceHitSet> &hitsets, const size_t totHits,
    const TimePoint *endTime,
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
  buildAllHits(hitsets, endTime, results);
  // if (d_params.randomSample) {
  //   buildRandomHits(hitsets, totHits, endTime, results);
  // } else {
  //   buildAllHits(hitsets, endTime, results);
  // }
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

void SynthonSpaceSearcher::buildAllHits(
    const std::vector<SynthonSpaceHitSet> &hitsets, const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  std::vector<std::tuple<const SynthonSet *, std::vector<size_t>>> toTry;
  for (const auto &[reactionId, synthonsToUse, numHits] : hitsets) {
    std::vector<std::vector<size_t>> synthonNums;
    synthonNums.reserve(synthonsToUse.size());
    std::vector<size_t> numSynthons;
    numSynthons.reserve(synthonsToUse.size());
    for (auto &stu : synthonsToUse) {
      numSynthons.push_back(stu.count());
      synthonNums.emplace_back();
      synthonNums.back().reserve(stu.count());
      for (size_t j = 0; j < stu.size(); ++j) {
        if (stu[j]) {
          synthonNums.back().push_back(j);
        }
      }
    }
    const auto &reaction = getSpace().getReactions().find(reactionId)->second;
    details::Stepper stepper(numSynthons);
    std::vector<size_t> theseSynthNums(synthonNums.size(), 0);
    while (stepper.d_currState[0] != numSynthons[0]) {
      for (size_t i = 0; i < stepper.d_currState.size(); ++i) {
        theseSynthNums[i] = synthonNums[i][stepper.d_currState[i]];
      }
      toTry.emplace_back(reaction.get(), theseSynthNums);
      stepper.step();
    }
  }

  // There are more than likely duplicate entries in toTry, because 2 different
  // fragmentations might produce overlapping synthon lists in the same
  // reaction. The duplicates need to be removed.
  // For most cases, it's quicker and just as effective to sort using the
  // address of the SynthonSet.  If the randomSeed has been set, use the
  // SynthonSet id for consistency - the addresses of the SynthonSets aren't
  // guaranteed to be in the same order from run to run, but the ids must be.
  // By default, std::sort just uses the first element of the tuple, which
  // isn't enough for this.
  std::cout << "Number to try: " << toTry.size() << std::endl;
  if (d_params.randomSeed == -1) {
    std::sort(
        toTry.begin(), toTry.end(),
        [](const std::tuple<const SynthonSet *, std::vector<size_t>> &rec1,
           const std::tuple<const SynthonSet *, std::vector<size_t>> &rec2) {
          if (std::get<0>(rec1) == std::get<0>(rec2)) {
            return std::get<1>(rec1) < std::get<1>(rec2);
          }
          return std::get<0>(rec1) < std::get<0>(rec2);
        });
  } else {
    std::sort(
        toTry.begin(), toTry.end(),
        [](const std::tuple<const SynthonSet *, std::vector<size_t>> &rec1,
           const std::tuple<const SynthonSet *, std::vector<size_t>> &rec2) {
          if (std::get<0>(rec1)->getId() == std::get<0>(rec2)->getId()) {
            return std::get<1>(rec1) < std::get<1>(rec2);
          }
          return std::get<0>(rec1)->getId() < std::get<0>(rec2)->getId();
        });
  }
  toTry.erase(std::unique(toTry.begin(), toTry.end()), toTry.end());
  if (d_params.randomSample) {
    std::shuffle(toTry.begin(), toTry.end(), *d_randGen);
  }
  makeHitsFromDetails(toTry, endTime, results);

  // Now get rid of any hits before d_params.hitStart.  It seems wasteful to
  // do it like this, but until a hit has been through verifyHit there's
  // no way of knowing whether it should be kept plus the threading makes it
  // very complicated to do otherwise.
  std::cout << "Number of hits : " << results.size()
            << " hitStart : " << d_params.hitStart << std::endl;
  if (d_params.hitStart) {
    if (d_params.hitStart < static_cast<std::int64_t>(results.size())) {
      std::for_each(results.begin(), results.begin() + d_params.hitStart,
                    [](std::unique_ptr<ROMol> &m) { m.reset(); });
      results.erase(std::remove_if(results.begin(), results.end(),
                                   [](const std::unique_ptr<ROMol> &r) {
                                     return !bool(r);
                                   }),
                    results.end());
    } else {
      results.clear();
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
    d_synthons.resize(hitsets.size());
    for (size_t hi = 0; hi < hitsets.size(); ++hi) {
      const SynthonSpaceHitSet &hs = hitsets[hi];
      d_synthons[hi] =
          std::vector<std::vector<size_t>>(hs.synthonsToUse.size());
      d_synthSels[hi] =
          std::vector<boost::random::uniform_int_distribution<size_t>>(
              hs.synthonsToUse.size());
      d_synthons[hi].resize(hs.synthonsToUse.size());
      for (size_t i = 0; i < hs.synthonsToUse.size(); ++i) {
        d_synthons[hi][i].reserve(hs.synthonsToUse[i].count());
        d_synthSels[hi][i] = boost::random::uniform_int_distribution<size_t>(
            0, hs.synthonsToUse[i].count() - 1);
        for (size_t j = 0; j < hs.synthonsToUse[i].size(); ++j) {
          if (hs.synthonsToUse[i][j]) {
            d_synthons[hi][i].push_back(j);
          }
        }
      }
    }
  }

  std::pair<std::string, std::vector<size_t>> selectSynthComb(
      boost::random::mt19937 &randGen) {
    std::vector<size_t> synths;
    const size_t hitSetNum = d_hitSetSel(randGen);
    for (size_t i = 0; i < d_hitsets[hitSetNum].synthonsToUse.size(); ++i) {
      const size_t synthNum = d_synthSels[hitSetNum][i](randGen);
      synths.push_back(d_synthons[hitSetNum][i][synthNum]);
    }
    return std::make_pair(d_hitsets[hitSetNum].reactionId, synths);
  }

  const std::vector<SynthonSpaceHitSet> &d_hitsets;
  const SynthonSpace &d_synthSpace;

  std::vector<size_t> d_hitSetWeights;
  boost::random::discrete_distribution<size_t> d_hitSetSel;
  std::vector<std::vector<std::vector<size_t>>> d_synthons;
  std::vector<std::vector<boost::random::uniform_int_distribution<size_t>>>
      d_synthSels;
};
}  // namespace

void SynthonSpaceSearcher::buildRandomHits(
    const std::vector<SynthonSpaceHitSet> &hitsets, const size_t totHits,
    const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
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
    if (auto prod = buildAndVerifyHit(reaction.get(), synths)) {
      results.push_back(std::move(prod));
    } else {
      numFails++;
    }
    // Don't check the time every go, as it's quite expensive.
    --numTries;
    if (!numTries) {
      numTries = 100;
      checkTimeOut(endTime);
      if (timedOut) {
        break;
      }
    }
  }
  sortHits(results);
}

namespace {
void processPartHitsFromDetails(
    const std::vector<std::tuple<const SynthonSet *, std::vector<size_t>>>
        &toDo,
    const TimePoint *endTime, std::vector<std::unique_ptr<ROMol>> &results,
    const SynthonSpaceSearcher *searcher, size_t startNum, size_t finishNum) {
  std::uint64_t numTries = 100;
  finishNum = finishNum > toDo.size() ? toDo.size() : finishNum;
  for (size_t i = startNum; i < finishNum; ++i) {
    auto prod =
        searcher->buildAndVerifyHit(std::get<0>(toDo[i]), std::get<1>(toDo[i]));
    results[i] = std::move(prod);
    // Don't check the time every go, as it's quite expensive.
    --numTries;
    if (results[i]) {
      if (!numTries) {
        numTries = 100;
        checkTimeOut(endTime);
        if (timedOut) {
          break;
        }
      }
      std::int64_t numHits = std::accumulate(
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
      if (searcher->getParams().maxHits != -1 &&
          numHits ==
              searcher->getParams().maxHits + searcher->getParams().hitStart) {
        break;
      }
    }
  }
}
}  // namespace

void SynthonSpaceSearcher::makeHitsFromDetails(
    const std::vector<std::tuple<const SynthonSet *, std::vector<size_t>>>
        &toDo,
    const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  results.resize(toDo.size());
#if RDK_BUILD_THREADSAFE_SSS
  auto numThreads = getNumThreadsToUse(d_params.numThreads);
  std::cout << "numThreads: " << numThreads << std::endl;
  if (numThreads > 1) {
    size_t eachThread = 1 + (toDo.size() / numThreads);
    size_t start = 0;
    std::vector<std::thread> threads;
    for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
      threads.push_back(std::thread(processPartHitsFromDetails, toDo, endTime,
                                    std::ref(results), this, start,
                                    start + eachThread));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    processPartHitsFromDetails(toDo, endTime, results, this, 0, toDo.size());
  }
#else
  processPartHitsFromDetails(toDo, endTime, results, this, 0, toDo.size());
#endif

  results.erase(
      std::remove_if(results.begin(), results.end(),
                     [](const std::unique_ptr<ROMol> &r) { return !bool(r); }),
      results.end());
  sortHits(results);
}

}  // namespace RDKit::SynthonSpaceSearch
