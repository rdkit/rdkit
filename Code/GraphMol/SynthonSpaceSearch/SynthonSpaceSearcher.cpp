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

#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SynthonSpaceSearch/ProgressBar.h>
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
      if (d_params.randomSeed == -1) {
        std::random_device rd;
        d_randGen = std::make_unique<std::mt19937>(rd());
      } else {
        d_randGen = std::make_unique<std::mt19937>(d_params.randomSeed);
      }
    }
  }
  // For the fragmentation, it is often useful to be able to keep track of the
  // original indices.
  for (auto atom : getQuery().atoms()) {
    atom->setProp<unsigned int>("ORIG_IDX", atom->getIdx());
  }
  for (auto bond : getQuery().bonds()) {
    bond->setProp<unsigned int>("ORIG_IDX", bond->getIdx());
  }
}

SearchResults SynthonSpaceSearcher::search(ThreadMode threadMode) {
  ControlCHandler::reset();
  std::vector<std::unique_ptr<ROMol>> results;
  const TimePoint *endTime = nullptr;
  TimePoint endTimePt;
  if (d_params.timeOut > 0) {
    endTimePt = Clock::now() + std::chrono::seconds(d_params.timeOut);
    endTime = &endTimePt;
  }
  bool timedOut = false;
  std::uint64_t totHits = 0;
  auto fragments = details::splitMolecule(
      d_query, getSpace().getMaxNumSynthons(), d_params.maxNumFragSets, endTime,
      d_params.numThreads, timedOut);
  std::cout << "Number of fragment sets : " << fragments.size() << std::endl;
  auto allHits = assembleHitSets(endTime, timedOut, totHits, threadMode);
  if (!timedOut && !ControlCHandler::getGotSignal() && d_params.buildHits) {
    buildHits(allHits, endTime, timedOut, results);
  }
  return SearchResults{std::move(results), std::move(d_bestHitFound), totHits,
                       timedOut, ControlCHandler::getGotSignal()};
}

void SynthonSpaceSearcher::search(const SearchResultCallback &cb,
                                  ThreadMode threadMode) {
  bool timedOut = false;
  const TimePoint *endTime = nullptr;
  std::uint64_t totHits = 0;
  auto allHits = assembleHitSets(endTime, timedOut, totHits, threadMode);
  std::sort(allHits.begin(), allHits.end(),
            [](const auto &hs1, const auto &hs2) -> bool {
              if (hs1->d_reaction->getId() == hs2->d_reaction->getId()) {
                return hs1->numHits < hs2->numHits;
              }
              return hs1->d_reaction->getId() < hs2->d_reaction->getId();
            });

  // from buildAllhits
  std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>> toTry;
  std::int64_t hitCount = 0;
  bool stop = false;  // set by callback

  // Each hitset contains possible hits from a single SynthonSet.
  for (const auto &hitset : allHits) {
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
        hitCount += partResults.size();
        stop = cb(partResults);
        toTry.clear();
        if (stop || (d_params.maxHits != -1 && hitCount >= d_params.maxHits)) {
          break;
        }
      }
      stepper.step();
    }
    if (stop || (d_params.maxHits != -1 && hitCount >= d_params.maxHits)) {
      break;
    }
  }

  // Do any remaining.
  if ((d_params.maxHits == -1 || hitCount < d_params.maxHits) && !stop &&
      !toTry.empty()) {
    std::vector<std::unique_ptr<ROMol>> partResults;
    processToTrySet(toTry, endTime, partResults);
    cb(partResults);
  }
}

std::unique_ptr<ROMol> SynthonSpaceSearcher::buildAndVerifyHit(
    const SynthonSpaceHitSet *hitset, const std::vector<size_t> &synthNums) {
  std::vector<const ROMol *> synths(synthNums.size());
  std::vector<const std::string *> synthNames(synthNums.size());

  std::unique_ptr<ROMol> prod;
  if (!quickVerify(hitset, synthNums)) {
    return prod;
  }

  for (size_t i = 0; i < synthNums.size(); i++) {
    synths[i] =
        hitset->synthonsToUse[i][synthNums[i]].second->getOrigMol().get();
    synthNames[i] = &(hitset->synthonsToUse[i][synthNums[i]].first);
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
  if (!verifyHit(*prod, hitset->d_reaction->getId(), synthNames)) {
    prod.reset();
  }
  return prod;
}

namespace {

void searchReactionPart(
    const std::vector<std::vector<std::unique_ptr<ROMol>>> &fragments,
    const TimePoint *endTime, std::atomic<std::int64_t> &mostRecentFrag,
    SynthonSpaceSearcher *searcher, const SynthonSet &reaction,
    std::vector<std::vector<std::unique_ptr<SynthonSpaceHitSet>>> &allSetHits,
    std::unique_ptr<ProgressBar> &pbar) {
  bool timedOut = false;
  int numTries = 1;

  std::int64_t lastFrag = fragments.size();
  while (true) {
    std::int64_t nextFrag = ++mostRecentFrag;
    if (nextFrag >= lastFrag) {
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
    allSetHits[nextFrag] =
        searcher->searchFragSet(fragments[nextFrag], reaction);
    if (pbar) {
      pbar->increment();
    }
  }
}

std::vector<std::unique_ptr<SynthonSpaceHitSet>> searchReaction(
    SynthonSpaceSearcher *searcher, const SynthonSet &reaction,
    const TimePoint *endTime, unsigned numThreads,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragments,
    std::unique_ptr<ProgressBar> &pbar) {
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> hits;
  std::vector<std::vector<std::unique_ptr<SynthonSpaceHitSet>>> allSetHits(
      fragments.size());
  std::atomic<std::int64_t> mostRecentFrag = -1;
  if (numThreads > 1) {
    std::vector<std::thread> threads;
    for (unsigned int i = 0u;
         i < std::min(static_cast<size_t>(numThreads), fragments.size()); ++i) {
      threads.push_back(std::thread(searchReactionPart, std::ref(fragments),
                                    endTime, std::ref(mostRecentFrag), searcher,
                                    std::ref(reaction), std::ref(allSetHits),
                                    std::ref(pbar)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    searchReactionPart(fragments, endTime, mostRecentFrag, searcher, reaction,
                       allSetHits, pbar);
  }
  for (auto &ash : allSetHits) {
    hits.insert(hits.end(), std::make_move_iterator(ash.begin()),
                std::make_move_iterator(ash.end()));
  }

  return hits;
}

void processReactions(
    SynthonSpaceSearcher *searcher,
    const std::vector<std::string> &reactionNames,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragments,
    const TimePoint *endTime, std::atomic<std::int64_t> &mostRecentReaction,
    std::int64_t lastReaction,
    std::vector<std::vector<std::unique_ptr<SynthonSpaceHitSet>>> &reactionHits,
    std::unique_ptr<ProgressBar> &pbar) {
  bool timedOut = false;

  while (true) {
    std::int64_t thisR = ++mostRecentReaction;
    if (thisR > lastReaction) {
      break;
    }
    timedOut = details::checkTimeOut(endTime);
    if (timedOut) {
      break;
    }
    const auto &reaction =
        searcher->getSpace().getReaction(reactionNames[thisR]);
    auto theseHits =
        searchReaction(searcher, *reaction, endTime, 1, fragments, pbar);
    reactionHits[thisR] = std::move(theseHits);
  }
}
}  // namespace

unsigned int SynthonSpaceSearcher::getNumQueryFragmentsRequired() {
  return getSpace().getMaxNumSynthons();
}

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceSearcher::assembleHitSets(const TimePoint *endTime, bool &timedOut,
                                      std::uint64_t &totHits,
                                      ThreadMode threadMode) {
  std::vector<std::unique_ptr<ROMol>> results;
  auto fragments = details::splitMolecule(
      d_query, getNumQueryFragmentsRequired(), d_params.maxNumFragSets, endTime,
      d_params.numThreads, timedOut);
  if (timedOut || ControlCHandler::getGotSignal()) {
    return {};
  }
  if (!extraSearchSetup(fragments, endTime) ||
      ControlCHandler::getGotSignal()) {
    return {};
  }
  return doTheSearch(fragments, endTime, timedOut, totHits, threadMode);
}

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceSearcher::doTheSearch(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets,
    const TimePoint *endTime, bool &timedOut, std::uint64_t &totHits,
    ThreadMode threadMode) {
  auto reactionNames = getSpace().getReactionNames();
  std::vector<std::vector<std::unique_ptr<SynthonSpaceHitSet>>> reactionHits;
  if (threadMode == ThreadMode::ThreadFragments) {
    const unsigned int numThreads = getNumThreadsToUse(d_params.numThreads);
    // Do the reactions one at a time, parallelising the fragSet searches.
    // For the slower searches, this minimises the amount of time that
    // threads are left idle.  ThreadReactions runs the risk of 1 thread
    // doing all the work if all the hits come from 1 reaction.
    std::unique_ptr<ProgressBar> pbar;
    if (getParams().useProgressBar) {
      std::cout << "\nSearching fragment sets" << std::endl;
      pbar.reset(new ProgressBar(getParams().useProgressBar,
                                 reactionNames.size() * fragSets.size()));
    }
    for (const auto &reactionName : reactionNames) {
      const auto &reaction = getSpace().getReaction(reactionName);
      reactionHits.push_back(
          searchReaction(this, *reaction, endTime, numThreads, fragSets, pbar));
    }
    if (pbar) {
      std::cout << std::endl;
    }
  } else {
    // Parallelise at the reaction level - each thread does all the fragSets
    // for a reaction in one go.  This is quicker for substructure searches.
    std::int64_t lastReaction = reactionNames.size() - 1;
    std::atomic<std::int64_t> mostRecentReaction = -1;
    reactionHits.resize(reactionNames.size());
    std::unique_ptr<ProgressBar> pbar;
    if (const auto numThreads = getNumThreadsToUse(d_params.numThreads);
        numThreads > 1) {
      std::vector<std::thread> threads;
      for (unsigned int i = 0u;
           i < std::min(static_cast<size_t>(numThreads), reactionNames.size());
           ++i) {
        threads.push_back(std::thread(
            processReactions, this, std::ref(reactionNames), std::ref(fragSets),
            endTime, std::ref(mostRecentReaction), lastReaction,
            std::ref(reactionHits), std::ref(pbar)));
      }
      for (auto &t : threads) {
        t.join();
      }
    } else {
      processReactions(this, reactionNames, fragSets, endTime,
                       mostRecentReaction, lastReaction, reactionHits, pbar);
    }
  }

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

bool SynthonSpaceSearcher::quickVerify(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  if (getParams().minHitHeavyAtoms || getParams().maxHitHeavyAtoms > -1) {
    unsigned int numHeavyAtoms = 0;
    for (unsigned int i = 0; i < synthNums.size(); ++i) {
      numHeavyAtoms +=
          hitset->synthonsToUse[i][synthNums[i]].second->getNumHeavyAtoms();
    }
    if (numHeavyAtoms < getParams().minHitHeavyAtoms ||
        (getParams().maxHitHeavyAtoms > -1 &&
         std::cmp_greater(numHeavyAtoms, getParams().maxHitHeavyAtoms))) {
      return false;
    }
  }
  if (getParams().minHitMolWt > 0.0 || getParams().maxHitMolWt > 0.0) {
    double molWt = 0.0;
    for (unsigned int i = 0; i < synthNums.size(); ++i) {
      molWt += hitset->synthonsToUse[i][synthNums[i]].second->getMolWt();
    }
    if (molWt < getParams().minHitMolWt ||
        (getParams().maxHitMolWt > 0.0 && molWt > getParams().maxHitMolWt)) {
      return false;
    }
  }
  if (getParams().minHitChiralAtoms || getParams().maxHitChiralAtoms > -1) {
    unsigned int numChiralAtoms = 0;
    unsigned int numExcDummies = 0;
    for (unsigned int i = 0; i < synthNums.size(); ++i) {
      numChiralAtoms +=
          hitset->synthonsToUse[i][synthNums[i]].second->getNumChiralAtoms();
      numExcDummies += hitset->synthonsToUse[i][synthNums[i]]
                           .second->getNumChiralAtomsExcDummies();
    }
    // numChiralAtoms is the upper bound on the number of chiral atoms in
    // the hit, numExcDummies the lower bound.
    if (numChiralAtoms < getParams().minHitChiralAtoms ||
        (getParams().maxHitChiralAtoms > -1 &&
         std::cmp_greater(numExcDummies, getParams().maxHitChiralAtoms))) {
      return false;
    }
  }
  return true;
}

void SynthonSpaceSearcher::updateBestHitSoFar(const ROMol &possBest,
                                              double sim) {
  if (sim > d_bestSimilarity) {
    d_bestSimilarity = sim;
    d_bestHitFound.reset(new ROMol(possBest));
    d_bestHitFound->setProp<double>("Similarity", sim);
  }
}

// It's conceivable that building the full molecule has changed the
// chiral atom count.
bool SynthonSpaceSearcher::verifyHit(ROMol &mol, const std::string &,
                                     const std::vector<const std::string *> &) {
  if (getParams().minHitChiralAtoms || getParams().maxHitChiralAtoms) {
    auto numChiralAtoms = details::countChiralAtoms(mol);
    if (numChiralAtoms < getParams().minHitChiralAtoms ||
        (getParams().maxHitChiralAtoms > -1 &&
         std::cmp_greater(numChiralAtoms, getParams().maxHitChiralAtoms))) {
      return false;
    }
  }
  return true;
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
    std::vector<std::unique_ptr<ROMol>> &results) {
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
    std::vector<std::unique_ptr<ROMol>> &results) {
  // toTry is the synthon numbers for a particular set of synthons
  // in the SynthonSet that will be zipped together to form a possible
  // hit molecule.  It will possibly hold possibilities from multiple
  // SynthonSets, and when there are enough to be worth processing,
  // they will be built into molecules, verified and accepted or
  // rejected as hits.
  std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>> toTry;
  bool enoughHits = false;

  if (hitsets.empty()) {
    BOOST_LOG(rdWarningLog) << "No possible synthon matches." << std::endl;
    return;
  }

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
    // process the synthons
    while (stepper.d_currState[0] != numSynthons[0]) {
      toTry.emplace_back(hitset.get(), stepper.d_currState);
      if (std::cmp_equal(toTry.size(), d_params.toTryChunkSize)) {
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
  if (d_params.maxHits != -1 &&
      std::cmp_greater(results.size(), d_params.maxHits + d_params.hitStart)) {
    results.erase(results.begin() + d_params.maxHits + d_params.hitStart,
                  results.end());
  }

  // Now get rid of any hits before d_params.hitStart.  It seems wasteful to
  // do it like this, but until a hit has been through verifyHit there's
  // no way of knowing whether it should be kept plus the threading makes it
  // very complicated to do otherwise.
  if (d_params.hitStart) {
    if (std::cmp_less(d_params.hitStart, results.size())) {
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
    SynthonSpaceSearcher *searcher, std::atomic<std::int64_t> &mostRecentTry,
    std::int64_t lastTry, std::unique_ptr<ProgressBar> &pbar) {
  std::uint64_t numTries = 1;
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
    if (pbar) {
      pbar->increment();
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
    const TimePoint *endTime, std::vector<std::unique_ptr<ROMol>> &results) {
  results.resize(toTry.size());
  std::int64_t lastTry = toTry.size() - 1;
  std::atomic<std::int64_t> mostRecentTry = -1;
  std::unique_ptr<ProgressBar> pbar;
  if (getParams().useProgressBar) {
    std::cout << "\nBuilding and checking hits." << std::endl;
    pbar.reset(new ProgressBar(getParams().useProgressBar, toTry.size()));
  }
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
                                    std::ref(mostRecentTry), lastTry,
                                    std::ref(pbar)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    processPartHitsFromDetails(toTry, endTime, results, this, mostRecentTry,
                               lastTry, pbar);
  }
  if (pbar) {
    std::cout << std::endl;
  }

  // Take out any gaps in the results set, where products didn't make the grade.
  results.erase(std::remove_if(results.begin(), results.end(),
                               [](const std::unique_ptr<ROMol> &r) {
                                 return !static_cast<bool>(r);
                               }),
                results.end());
}

void SynthonSpaceSearcher::sortToTryByApproxSimilarity(
    std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
        &toTry) const {
  std::vector<std::pair<size_t, double>> tmp;
  // toTry will always have something in it at this point.
  tmp.reserve(toTry.size());
  for (size_t i = 0; i < toTry.size(); i++) {
    tmp.emplace_back(i, approxSimilarity(toTry[i].first, toTry[i].second));
  }
  std::sort(tmp.begin(), tmp.end(),
            [](const auto &lhs, const auto &rhs) -> bool {
              return lhs.second > rhs.second;
            });
  BOOST_LOG(rdInfoLog) << "Best approx similarity : " << tmp.front().second
                       << " and worst : " << tmp.back().second << std::endl;
  std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
      newToTry;
  newToTry.reserve(tmp.size());
  std::transform(tmp.begin(), tmp.end(), back_inserter(newToTry),
                 [&](const auto &p) -> auto { return toTry[p.first]; });
  toTry = newToTry;
}

void SynthonSpaceSearcher::processToTrySet(
    std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
        &toTry,
    const TimePoint *endTime, std::vector<std::unique_ptr<ROMol>> &results) {
  // There are possibly duplicate entries in toTry, because 2
  // different fragmentations might produce overlapping synthon lists in
  // the same reaction. The duplicates need to be removed.
  details::sortAndUniquifyToTry(toTry);

  if (d_params.randomSample) {
    std::shuffle(toTry.begin(), toTry.end(), *d_randGen);
  } else {
    sortToTryByApproxSimilarity(toTry);
  }
  makeHitsFromToTry(toTry, endTime, results);
}
}  // namespace RDKit::SynthonSpaceSearch
