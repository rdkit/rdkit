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
std::vector<boost::dynamic_bitset<>> getHitSynthons(
    const std::vector<std::unique_ptr<ROMol>> &fragSet,
    const RascalMCES::RascalOptions &rascalOptions,
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  synthonsToUse.reserve(reaction->getSynthons().size());
  for (const auto &synthonSet : reaction->getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  for (size_t i = 0; i < synthonOrder.size(); i++) {
    const auto &synthons = reaction->getSynthons()[synthonOrder[i]];
    bool fragMatched = false;
    for (size_t j = 0; j < synthons.size(); j++) {
      const auto rascalResults = RascalMCES::rascalMCES(
          *fragSet[i], *synthons[j]->getSearchMol(), rascalOptions);
      if (!rascalResults.empty()) {
        synthonsToUse[synthonOrder[i]][j] = true;
        fragMatched = true;
      }
    }
    if (!fragMatched) {
      synthonsToUse.clear();
      return synthonsToUse;
    }
  }
  // Fill in any synthons where they all didn't match.
  details::expandBitSet(synthonsToUse);
  return synthonsToUse;
}
}  // namespace

std::vector<SynthonSpaceHitSet> SynthonSpaceRascalSearcher::searchFragSet(
    std::vector<std::unique_ptr<ROMol>> &fragSet) const {
  std::vector<SynthonSpaceHitSet> results;
  for (auto &frag : fragSet) {
    // Ring info is required.
    unsigned int otf;
    sanitizeMol(*static_cast<RWMol *>(frag.get()), otf,
                MolOps::SANITIZE_SYMMRINGS);
  }

  const auto connPatterns = details::getConnectorPatterns(fragSet);
  boost::dynamic_bitset<> conns(MAX_CONNECTOR_NUM + 1);
  for (auto &connPattern : connPatterns) {
    conns |= connPattern;
  }

  for (const auto &[id, reaction] : getSpace().getReactions()) {
    // It can't be a hit if the number of fragments is more than the number
    // of synthon sets because some of the molecule won't be matched in any
    // of the potential products.  It can be less, in which case the unused
    // synthon set will be used completely, possibly resulting in a large
    // number of hits.
    if (fragSet.size() > reaction->getSynthons().size()) {
      continue;
    }
    auto synthConnPatts = reaction->getSynthonConnectorPatterns();

    // Need to try all combinations of synthon orders.
    auto synthonOrders =
        details::permMFromN(fragSet.size(), reaction->getSynthons().size());
    for (const auto &synthonOrder : synthonOrders) {
      // Get all the possible permutations of connector numbers compatible with
      // the number of synthon sets in this reaction.  So if the
      // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
      // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
      // that might be how they're labelled in the reaction database.
      auto connCombs = details::getConnectorPermutations(
          fragSet, conns, reaction->getConnectors());

      for (auto &connComb : connCombs) {
        // Make sure that for this connector combination, the synthons in this
        // order have something similar.  All query fragment connectors must
        // match something in the corresponding synthon.  The synthon can
        // have unused connectors.
        auto connCombConnPatterns = details::getConnectorPatterns(connComb);
        bool skip = false;
        for (size_t i = 0; i < connCombConnPatterns.size(); ++i) {
          if ((connCombConnPatterns[i] & synthConnPatts[synthonOrder[i]])
                  .count() < connPatterns[i].count()) {
            skip = true;
            break;
          }
        }
        if (skip) {
          continue;
        }
        // Rascal ignores isotope numbers, which makes things easier.
        auto theseSynthons = getHitSynthons(fragSet, d_rascalFragOptions,
                                            reaction, synthonOrder);
        if (!theseSynthons.empty()) {
          const size_t numHits = std::accumulate(
              theseSynthons.begin(), theseSynthons.end(), 1,
              [](const int prevRes, const boost::dynamic_bitset<> &s2) {
                return prevRes * s2.count();
              });
          if (numHits) {
            results.push_back(
                SynthonSpaceHitSet{reaction->getId(), theseSynthons, numHits});
          }
        }
      }
    }
  }
  return results;
}

bool SynthonSpaceRascalSearcher::quickVerify(
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<size_t> &synthNums) const {
  // If the query is an exact substructure of the product, then that's an upper
  // bound on the Johnson similarity.  Check that that is not below the
  // threshold.
  int qbit = getQuery().getNumAtoms() + getQuery().getNumBonds();
  const auto &rsynths = reaction->getSynthons();
  int numAtoms = 0, numBonds = 0;
  for (size_t i = 0; i < synthNums.size(); i++) {
    const auto &smol = rsynths[i][synthNums[i]]->getSearchMol();
    // Adjust for connector points that aren't in the final product.
    numAtoms += smol->getNumAtoms() -
                reaction->getSynthonConnectorPatterns()[i].count();
    numBonds += smol->getNumBonds() -
                reaction->getSynthonConnectorPatterns()[i].count();
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
