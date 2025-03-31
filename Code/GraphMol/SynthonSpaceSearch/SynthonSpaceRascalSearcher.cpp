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
    const SynthonSpaceSearchParams &params, SynthonSpace &space)
    : SynthonSpaceSearcher(query, params, space),
      d_rascalOptions(rascalOptions),
      d_rascalFragOptions(rascalOptions) {
  d_rascalFragOptions.similarityThreshold -= params.fragSimilarityAdjuster;
}

namespace {
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<std::unique_ptr<ROMol>> &fragSet,
    const RascalMCES::RascalOptions &rascalOptions, const SynthonSet &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;

  synthonsToUse.reserve(reaction.getSynthons().size());
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  for (size_t i = 0; i < synthonOrder.size(); i++) {
    const auto &synthons = reaction.getSynthons()[synthonOrder[i]];
    bool fragMatched = false;
    for (size_t j = 0; j < synthons.size(); j++) {
      const auto rascalResults = RascalMCES::rascalMCES(
          *fragSet[i], *synthons[j].second->getSearchMol(), rascalOptions);
      if (!rascalResults.empty()) {
        synthonsToUse[synthonOrder[i]][j] = true;
        fragMatched = true;
      }
    }
    if (!fragMatched) {
      return retSynthons;
    }
  }
  // Fill in any synthons where they all didn't match.
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

void SynthonSpaceRascalSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) {
  for (const auto &fragSet : fragSets) {
    for (const auto &frag : fragSet) {
      unsigned int otf;
      sanitizeMol(*static_cast<RWMol *>(frag.get()), otf,
                  MolOps::SANITIZE_SYMMRINGS);
    }
  }
}

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceRascalSearcher::searchFragSet(
    const std::vector<std::unique_ptr<ROMol>> &fragSet,
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
  // If the query is an exact substructure of the product, then that's an upper
  // bound on the Johnson similarity.  Check that that is not below the
  // threshold.
  int qbit = getQuery().getNumAtoms() + getQuery().getNumBonds();
  int numAtoms = 0, numBonds = 0;
  for (size_t i = 0; i < synthNums.size(); i++) {
    const auto &smol = hitset->synthonsToUse[i][synthNums[i]].second;
    // Adjust for connector points that aren't in the final product.
    numAtoms += smol->getNumAtoms() -
                hitset->d_reaction->getSynthonConnectorPatterns()[i].count();
    numBonds += smol->getNumBonds() -
                hitset->d_reaction->getSynthonConnectorPatterns()[i].count();
  }
  // The Johnson similarity is
  // (commonNatoms + commonNbonds)**2 /
  // ((Natoms1 + Nbonds1) * (Natoms2 + Natoms2))
  // and in this case the common atoms are the whole query, so the square
  // cancels out.
  double bestSim =
      static_cast<double>(qbit) / static_cast<double>(numAtoms + numBonds);
  return bestSim >= d_rascalOptions.similarityThreshold;
}

bool SynthonSpaceRascalSearcher::verifyHit(const ROMol &hit) const {
  auto res = RascalMCES::rascalMCES(hit, getQuery(), d_rascalOptions);
  // Rascal reports all matches that proceed to full MCES elucidation,
  // even if the final similarity value ends up below the threshold.
  // We only want those over the threshold.
  if (!res.empty() &&
      res.front().getSimilarity() >= d_rascalOptions.similarityThreshold) {
    hit.setProp<double>("Similarity", res.front().getSimilarity());
    return true;
  }
  return false;
}

}  // namespace RDKit::SynthonSpaceSearch
