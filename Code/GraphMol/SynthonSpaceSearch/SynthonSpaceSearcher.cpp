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
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>

namespace RDKit::SynthonSpaceSearch {

SubstructureResults SynthonSpaceSearcher::search() {
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
  RDKit::MatchVectType dontCare;

  auto fragments = details::splitMolecule(d_query, d_params.maxBondSplits);
  std::vector<SynthonSpaceHitSet> allHits;
  size_t totHits = 0;
  for (auto &fragSet : fragments) {
    // std::cout << "Next frag set : ";
    // for (const auto &f : fragSet) {
    //   std::cout << MolToSmiles(*f) << "  ";
    // }
    // std::cout << std::endl;
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
    std::sort(allHits.begin(), allHits.end(),
              [](const SynthonSpaceHitSet &hs1,
                 const SynthonSpaceHitSet &hs2) -> bool {
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
    // The random sampling, if requested, doesn't always produce the required
    // number of hits in 1 sweep.  Also, totHits, being an upper bound, may
    // not be achievable.
    size_t tmpMaxHits(d_params.maxHits);
    int numSweeps = 0;
    while (results.size() < std::min(tmpMaxHits, totHits) &&
           numSweeps < d_params.numRandomSweeps) {
      buildHits(allHits, totHits, resultsNames, results);
      // totHits is an upper bound, so may not be reached.
      if (d_params.maxHits == -1 || !d_params.randomSample) {
        break;
      }
      numSweeps++;
    }
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

void SynthonSpaceSearcher::buildHits(
    const std::vector<SynthonSpaceHitSet> &hitsets, const size_t totHits,
    std::set<std::string> &resultsNames,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  const double randDiscrim =
      static_cast<double>(d_params.maxHits) / static_cast<double>(totHits);

  for (const auto &hitset : hitsets) {
    const auto &synthonsToUse = hitset.synthonsToUse;
    auto synthons = getSynthonsToUse(synthonsToUse, hitset.reactionId);
    if (synthons.empty()) {
      return;
    }

    std::vector<size_t> numSynthons;
    numSynthons.reserve(synthons.size());
    for (auto &s : synthons) {
      numSynthons.push_back(s.size());
    }
    Stepper stepper(numSynthons);
    const size_t numReactions = synthonsToUse.size();
    MolzipParams mzparams;
    mzparams.label = MolzipLabel::Isotope;
    while (stepper.d_currState[0] != numSynthons[0]) {
      std::string combName =
          hitset.reactionId + "_" +
          synthons[0][stepper.d_currState[0]]->getProp<std::string>(
              common_properties::_Name);
      for (size_t i = 1; i < numReactions; ++i) {
        combName +=
            "_" + synthons[i][stepper.d_currState[i]]->getProp<std::string>(
                      common_properties::_Name);
      }
      if (!d_params.randomSample || d_params.maxHits == -1 ||
          (d_params.randomSample && dist(*d_randGen) < randDiscrim)) {
        if (resultsNames.insert(combName).second) {
          if (resultsNames.size() < static_cast<size_t>(d_params.hitStart)) {
            continue;
          }
          auto combMol = std::make_unique<ROMol>(
              ROMol(*synthons[0][stepper.d_currState[0]]));
          for (size_t i = 1; i < numReactions; ++i) {
            combMol.reset(
                combineMols(*combMol, *synthons[i][stepper.d_currState[i]]));
          }
          auto prod = molzip(*combMol, mzparams);
          MolOps::sanitizeMol(*dynamic_cast<RWMol *>(prod.get()));
          // Do a final check of the whole thing.  It can happen that the
          // fragments match synthons but the final product doesn't match, and
          // a key example is when the 2 synthons come together to form an
          // aromatic ring.  For substructure searching, an aliphatic query
          // can match the aliphatic synthon so they are selected as a hit,
          // but the final aromatic ring isn't a match.
          // E.g. Cc1cccc(C(=O)N[1*])c1N=[2*] and c1ccoc1C(=[2*])[1*]
          // making Cc1cccc2c(=O)[nH]c(-c3ccco3)nc12.  The query c1ccc(CN)o1
          // when split is a match to the synthons (c1ccc(C[1*])o1 and [1*]N)
          // but the product the hydroxyquinazoline is aromatic, at least in
          // the RDKit model so the N in the query doesn't match.
          if (!verifyHit(*prod)) {
            continue;
          }
          prod->setProp<std::string>(common_properties::_Name, combName);
          results.push_back(std::move(prod));
        }
      }
      if (results.size() == static_cast<size_t>(d_params.maxHits)) {
        return;
      }
      stepper.step();
    }
  }
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
