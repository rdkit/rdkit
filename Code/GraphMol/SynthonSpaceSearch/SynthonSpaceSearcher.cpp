//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/MolOps.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>

namespace RDKit::SynthonSpaceSearch {

SubstructureResults SynthonSpaceSearcher::search() {
  if (d_params.randomSample && d_params.maxHits == -1) {
    throw std::runtime_error(
        "Random sample is incompatible with maxHits of -1.");
  }

  if (d_params.randomSample) {
    if (!d_randGen) {
      std::random_device rd;
      d_randGen = std::make_unique<std::mt19937>(rd());
    }
    if (d_params.randomSeed != -1) {
      d_randGen->seed(d_params.randomSeed);
    }
  }
  std::vector<std::unique_ptr<ROMol>> results;

  auto fragments = details::splitMolecule(d_query, d_params.maxBondSplits);
  std::cout << "Number of fragments: " << fragments.size() << std::endl;
  std::vector<SynthonSpaceHitSet> allHits;
  size_t totHits = 0;
  for (auto &fragSet : fragments) {
    auto theseHits = searchFragSet(fragSet);
    if (!theseHits.empty()) {
      totHits += std::accumulate(
          theseHits.begin(), theseHits.end(), 0,
          [](const size_t prevVal, const SynthonSpaceHitSet &hs) -> size_t {
            return prevVal + hs.numHits;
          });
      allHits.insert(allHits.end(), theseHits.begin(), theseHits.end());
    }
  }

  if (d_params.buildHits) {
    buildHits(allHits, totHits, results);
  }
  return SubstructureResults{std::move(results), totHits};
}

namespace {
// class to step through all combinations of list of different sizes.
// returns (0,0,0), (0,0,1), (0,1,0) etc.
struct Stepper {
  explicit Stepper(std::vector<size_t> &sizes) : d_sizes(sizes) {
    d_currState = std::vector<size_t>(sizes.size(), 0);
  }
  void step() {
    // Don't do anything if we're at the end, but expect an infinite
    // loop if the user isn't wise to this.
    if (d_currState[0] == d_sizes[0]) {
      return;
    }
    std::int64_t i = static_cast<std::int64_t>(d_currState.size()) - 1;
    while (i >= 0) {
      ++d_currState[i];
      if (d_currState[0] == d_sizes[0]) {
        return;
      }
      if (d_currState[i] == d_sizes[i]) {
        d_currState[i] = 0;
      } else {
        break;
      }
      --i;
    }
  }
  std::vector<size_t> d_currState;
  std::vector<size_t> d_sizes;
};

}  // namespace

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

void SynthonSpaceSearcher::buildHits(
    std::vector<SynthonSpaceHitSet> &hitsets, const size_t totHits,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }
  std::sort(
      hitsets.begin(), hitsets.end(),
      [](const SynthonSpaceHitSet &hs1, const SynthonSpaceHitSet &hs2) -> bool {
        if (hs1.reactionId == hs2.reactionId) {
          return hs1.numHits < hs2.numHits;
        } else {
          return hs1.reactionId < hs2.reactionId;
        }
      });
  // Keep track of the result names so we can weed out duplicates by
  // reaction and synthons.  Different splits may give rise to the same
  // synthon combination.  This will keep the same molecule produced via
  // different reactions which I think makes sense.  The resultsNames will
  // be accumulated even if the molecule itself doesn't make it into the
  // results set, for example if it isn't a random selection or it's
  // outside maxHits or hitStart.
  std::set<std::string> resultsNames;

  std::cout << "Upper bound on hits : " << totHits << std::endl;
  if (d_params.randomSample) {
    buildRandomHits(hitsets, totHits, resultsNames, results);
  } else {
    buildAllHits(hitsets, resultsNames, results);
  }
}

namespace {
void sortHits(std::vector<std::unique_ptr<ROMol>> &hits) {
  if (!hits.empty() && hits.front()->hasProp("Similarity")) {
    std::sort(hits.begin(), hits.end(),
              [](const std::unique_ptr<ROMol> &lhs,
                 const std::unique_ptr<ROMol> &rhs) {
                double lsim = lhs->getProp<double>("Similarity");
                double rsim = rhs->getProp<double>("Similarity");
                return lsim > rsim;
              });
  }
}
}  // namespace

void SynthonSpaceSearcher::buildAllHits(
    const std::vector<SynthonSpaceHitSet> &hitsets,
    std::set<std::string> &resultsNames,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  for (const auto &hitset : hitsets) {
    const auto &synthonsToUse = hitset.synthonsToUse;
    std::vector<std::vector<size_t>> synthonNums;
    synthonNums.reserve(synthonsToUse.size());
    std::vector<size_t> numSynthons;
    numSynthons.reserve(synthonsToUse.size());
    for (auto &stu : synthonsToUse) {
      numSynthons.push_back(stu.count());
      synthonNums.push_back(std::vector<size_t>());
      synthonNums.back().reserve(stu.count());
      for (size_t j = 0; j < stu.size(); ++j) {
        if (stu[j]) {
          synthonNums.back().push_back(j);
        }
      }
    }
    const auto &reaction =
        getSpace().getReactions().find(hitset.reactionId)->second;
    Stepper stepper(numSynthons);
    std::vector<size_t> theseSynthNums(synthonNums.size(), 0);
    while (stepper.d_currState[0] != numSynthons[0]) {
      for (size_t i = 0; i < stepper.d_currState.size(); ++i) {
        theseSynthNums[i] = synthonNums[i][stepper.d_currState[i]];
      }
      auto prod = buildAndVerifyHit(reaction, theseSynthNums, resultsNames);
      if (prod) {
        results.push_back(std::move(prod));
      }
      if (results.size() == static_cast<size_t>(d_params.maxHits)) {
        return;
      }
      stepper.step();
    }
  }
  sortHits(results);
}

namespace {
struct RandomHitSelector {
  RandomHitSelector(const std::vector<SynthonSpaceHitSet> &hitsets,
                    const SynthonSpace &space)
      : d_hitsets(hitsets), d_synthSpace(space) {
    d_hitSetSel = std::uniform_int_distribution<size_t>(0, hitsets.size() - 1);
    d_synthons.resize(hitsets.size());
    d_synthSels.resize(hitsets.size());
    for (size_t hi = 0; hi < hitsets.size(); ++hi) {
      const SynthonSpaceHitSet &hs = hitsets[hi];
      d_synthons[hi] =
          std::vector<std::vector<size_t>>(hs.synthonsToUse.size());
      d_synthSels[hi] = std::vector<std::uniform_int_distribution<size_t>>(
          hs.synthonsToUse.size());
      d_synthons[hi].resize(hs.synthonsToUse.size());
      for (size_t i = 0; i < hs.synthonsToUse.size(); ++i) {
        d_synthons[hi][i].reserve(hs.synthonsToUse[i].count());
        d_synthSels[hi][i] = std::uniform_int_distribution<size_t>(
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
      std::mt19937 &randGen) {
    std::vector<size_t> synths;
    size_t hitSetNum = d_hitSetSel(randGen);
    for (size_t i = 0; i < d_hitsets[hitSetNum].synthonsToUse.size(); ++i) {
      size_t synthNum = d_synthSels[hitSetNum][i](randGen);
      synths.push_back(d_synthons[hitSetNum][i][synthNum]);
    }
    return std::make_pair(d_hitsets[hitSetNum].reactionId, synths);
  }

  const std::vector<SynthonSpaceHitSet> &d_hitsets;
  const SynthonSpace &d_synthSpace;

  std::uniform_int_distribution<size_t> d_hitSetSel;
  std::vector<std::vector<std::vector<size_t>>> d_synthons;
  std::vector<std::vector<std::uniform_int_distribution<size_t>>> d_synthSels;
};
}  // namespace

void SynthonSpaceSearcher::buildRandomHits(
    const std::vector<SynthonSpaceHitSet> &hitsets, const size_t totHits,
    std::set<std::string> &resultsNames,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }
  auto rhs = RandomHitSelector(hitsets, d_space);

  uint64_t numFails = 0;
  while (results.size() <
             std::min(static_cast<const std::uint64_t>(d_params.maxHits),
                      static_cast<std::uint64_t>(totHits)) &&
         numFails < totHits * d_params.numRandomSweeps) {
    const auto &[reactionId, synths] = rhs.selectSynthComb(*d_randGen);
    const auto &reaction = getSpace().getReactions().find(reactionId)->second;
    auto prod = buildAndVerifyHit(reaction, synths, resultsNames);
    if (prod) {
      results.push_back(std::move(prod));
    } else {
      numFails++;
    }
  }
  sortHits(results);
}

std::vector<std::vector<ROMol *>> SynthonSpaceSearcher::getSynthonsToUse(
    const std::vector<boost::dynamic_bitset<>> &synthonsToUse,
    const std::string &reaction_id) const {
  if (const auto &it = getSpace().getReactions().find(reaction_id);
      it == getSpace().getReactions().end()) {
    throw std::runtime_error("Reaction " + reaction_id +
                             "not in the reaction set.");
  }
  const auto &reaction = getSpace().getReactions().find(reaction_id)->second;
  return reaction->getSynthons(synthonsToUse);
}

}  // namespace RDKit::SynthonSpaceSearch