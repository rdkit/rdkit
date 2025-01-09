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
#include <iostream>
#include <locale>
#include <random>
#include <sstream>
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

// When debugging in the synthon space the numbers get very large
// and it is useful to have them formatted like this.
template <class Char>
class MyFacet : public std::numpunct<Char> {
 public:
  std::string do_grouping() const { return "\3"; }
  Char do_thousands_sep() const { return ' '; }
};

std::string formatLargeInt(std::int64_t n) {
  std::ostringstream oss;
  oss.imbue(std::locale(oss.getloc(), new MyFacet<char>));
  oss << n;
  return oss.str();
}

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
      allFragHits[i] = self->searchFragSet(std::ref(fragments[i]));
    }
  };

  auto numThreads = getNumThreadsToUse(d_params.numThreads);
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

  std::int64_t totHits = 0;
  std::vector<SynthonSpaceHitSet> allHits;
  for (const auto &fh : allFragHits) {
    totHits += std::accumulate(
        fh.begin(), fh.end(), std::int64_t(0),
        [](const size_t prevVal, const SynthonSpaceHitSet &hs) -> std::int64_t {
          return prevVal + hs.numHits;
        });
    allHits.insert(allHits.end(), fh.begin(), fh.end());
  }

  std::vector<std::unique_ptr<ROMol>> results;
  if (!timedOut && d_params.buildHits) {
    buildHits(allHits, endTime, results);
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
    std::vector<SynthonSpaceHitSet> &hitsets, const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }
  if (d_params.randomSample) {
    std::shuffle(hitsets.begin(), hitsets.end(), *d_randGen);
  } else {
    std::sort(hitsets.begin(), hitsets.end(),
              [](const SynthonSpaceHitSet &hs1,
                 const SynthonSpaceHitSet &hs2) -> bool {
                if (hs1.reactionId == hs2.reactionId) {
                  return hs1.numHits < hs2.numHits;
                }
                return hs1.reactionId < hs2.reactionId;
              });
  }
  buildAllHits(hitsets, endTime, results);
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

void sortAndUniquifyToTry(
    std::vector<std::tuple<const SynthonSet *, std::vector<size_t>>> &toTry,
    int randomSeed) {
  // For most cases, it's quicker and just as effective to sort using the
  // address of the SynthonSet.  If the randomSeed has been set, use the
  // SynthonSet id for consistency - the addresses of the SynthonSets aren't
  // guaranteed to be in the same order from run to run, but the ids must be.
  // By default, std::sort just uses the first element of the tuple, which
  // isn't enough for this.
  if (randomSeed == -1) {
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
}

bool haveEnoughHits(const std::vector<std::unique_ptr<ROMol>> &results,
                    std::int64_t maxHits, std::int64_t hitStart) {
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
  if (maxHits != -1 && numHits >= maxHits + hitStart) {
    return true;
  }
  return false;
}
}  // namespace

void SynthonSpaceSearcher::buildAllHits(
    const std::vector<SynthonSpaceHitSet> &hitsets, const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  std::vector<std::tuple<const SynthonSet *, std::vector<size_t>>> toTry;
  bool enoughHits = false;
  for (const auto &[reactionId, synthonsToUse, numHits] : hitsets) {
    // Set up the stepper to move through the synthons.
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
    details::Stepper stepper(numSynthons);

    std::vector<size_t> theseSynthNums(synthonNums.size(), 0);
    const auto &reaction = getSpace().getReactions().find(reactionId)->second;
    // process the synthons
    while (stepper.d_currState[0] != numSynthons[0]) {
      for (size_t i = 0; i < stepper.d_currState.size(); ++i) {
        theseSynthNums[i] = synthonNums[i][stepper.d_currState[i]];
      }
      toTry.emplace_back(reaction.get(), theseSynthNums);
      if (toTry.size() == static_cast<size_t>(d_params.toTryChunkSize)) {
        std::vector<std::unique_ptr<ROMol>> partResults;
        processToTrySet(toTry, endTime, partResults);
        results.insert(results.end(),
                       std::make_move_iterator(partResults.begin()),
                       std::make_move_iterator(partResults.end()));
        partResults.clear();
        enoughHits =
            haveEnoughHits(results, d_params.maxHits, d_params.hitStart);
        toTry.clear();
        if (enoughHits || timedOut) {
          break;
        }
      }
      stepper.step();
    }
    if (enoughHits || timedOut) {
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
                                     return !bool(r);
                                   }),
                    results.end());
    } else {
      results.clear();
    }
  }
}

namespace {
void processPartHitsFromDetails(
    const std::vector<std::tuple<const SynthonSet *, std::vector<size_t>>>
        &toTry,
    const TimePoint *endTime, std::vector<std::unique_ptr<ROMol>> &results,
    const SynthonSpaceSearcher *searcher, size_t startNum, size_t finishNum) {
  std::uint64_t numTries = 100;
  finishNum = finishNum > toTry.size() ? toTry.size() : finishNum;
  for (size_t i = startNum; i < finishNum; ++i) {
    auto prod = searcher->buildAndVerifyHit(std::get<0>(toTry[i]),
                                            std::get<1>(toTry[i]));
    if (prod) {
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
      checkTimeOut(endTime);
      if (timedOut) {
        break;
      }
    }
  }
}
}  // namespace

void SynthonSpaceSearcher::makeHitsFromToTry(
    const std::vector<std::tuple<const SynthonSet *, std::vector<size_t>>>
        &toTry,
    const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  results.resize(toTry.size());
#if RDK_BUILD_THREADSAFE_SSS
  auto numThreads = getNumThreadsToUse(d_params.numThreads);
  if (numThreads > 1) {
    size_t eachThread = 1 + (toTry.size() / numThreads);
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
    processPartHitsFromDetails(toTry, endTime, results, this, 0, toTry.size());
  }
#else
  processPartHitsFromDetails(toTry, endTime, results, this, 0, toTry.size());
#endif

  // Take out any gaps in the results set, where products didn't make the grade.
  results.erase(
      std::remove_if(results.begin(), results.end(),
                     [](const std::unique_ptr<ROMol> &r) { return !bool(r); }),
      results.end());
}

void SynthonSpaceSearcher::processToTrySet(
    std::vector<std::tuple<const SynthonSet *, std::vector<size_t>>> &toTry,
    const TimePoint *endTime,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  // There are possibly duplicate entries in toTry, because 2
  // different fragmentations might produce overlapping synthon lists in
  // the same reaction. The duplicates need to be removed.  Although
  // when doing the job in batches this is less likely.
  sortAndUniquifyToTry(toTry, d_params.randomSeed);

  if (d_params.randomSample) {
    std::shuffle(toTry.begin(), toTry.end(), *d_randGen);
  }
  makeHitsFromToTry(toTry, endTime, results);
}

}  // namespace RDKit::SynthonSpaceSearch
