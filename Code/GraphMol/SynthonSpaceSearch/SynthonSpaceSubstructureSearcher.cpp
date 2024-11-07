//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/QueryAtom.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSubstructureSearcher.h>

namespace RDKit::SynthonSpaceSearch {

namespace {

// Make pattern fps for the fragments and re-order fps and fragments into
// descending order of number of bits set in the fp (the "largest fragment
// heuristic").  The FP with the largest number of bits
// is the most likely to screen out a matching synthon set since smaller,
// less complex fragments are more likely to match something, so screen
// with that first.
std::vector<std::unique_ptr<ExplicitBitVect>> makePatternFPs(
    std::vector<std::unique_ptr<ROMol>> &molFrags) {
  std::vector<std::unique_ptr<ExplicitBitVect>> pattFPs;
  pattFPs.reserve(molFrags.size());
  for (const auto &frag : molFrags) {
    pattFPs.emplace_back(PatternFingerprintMol(*frag, 2048));
  }
  // Sort by descending number of bits set.
  std::vector<std::pair<size_t, ExplicitBitVect *>> fps(pattFPs.size());
  for (size_t i = 0; i < pattFPs.size(); ++i) {
    fps[i] = std::make_pair(i, pattFPs[i].get());
  }
  std::sort(fps.begin(), fps.end(),
            [](const std::pair<size_t, ExplicitBitVect *> &fp1,
               const std::pair<size_t, ExplicitBitVect *> &fp2) -> bool {
              return fp1.second->getNumOnBits() > fp2.second->getNumOnBits();
            });

  // Now put molFrags in the same order.
  std::vector<std::unique_ptr<ROMol>> newFrags(molFrags.size());
  std::vector<std::unique_ptr<ExplicitBitVect>> retFPs(molFrags.size());
  for (size_t i = 0; i < fps.size(); ++i) {
    newFrags[i] = std::move(molFrags[fps[i].first]);
    retFPs[i] = std::move(pattFPs[fps[i].first]);
  }
  molFrags = std::move(newFrags);
  return retFPs;
}

// Return a bitset giving the different connector types in this
// molecule.
boost::dynamic_bitset<> getConnectorPattern(
    const std::vector<std::unique_ptr<ROMol>> &fragSet) {
  boost::dynamic_bitset<> conns(MAX_CONNECTOR_NUM);
  for (const auto &frag : fragSet) {
    for (const auto &a : frag->atoms()) {
      if (!a->getAtomicNum() && a->getIsotope()) {
        conns.set(a->getIsotope() - 1);
      }
    }
  }
  return conns;
}

// Return copies of the mol fragments will all permutations of the connectors
// in the reaction onto the connectors in the fragments.
// E.g. if the reaction has 3 connectors, 1, 2 and 3 and the fragged mol has
// 2, return all permutations of 2 from 3.  It's ok if the fragged mol doesn't
// have all the connections in the reaction, although this may well result in
// a lot of hits.
std::vector<std::vector<std::unique_ptr<RWMol>>> getConnectorPermutations(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    const boost::dynamic_bitset<> &fragConns,
    const boost::dynamic_bitset<> &reactionConns) {
  std::vector<std::vector<std::unique_ptr<RWMol>>> connPerms;
  auto bitsToInts =
      [](const boost::dynamic_bitset<> &bits) -> std::vector<int> {
    std::vector<int> ints;
    for (size_t i = 0; i < bits.size(); ++i) {
      if (bits[i]) {
        ints.push_back(static_cast<int>(i));
      }
    }
    return ints;
  };
  auto numFragConns = fragConns.count();
  auto rConns = bitsToInts(reactionConns);
  auto perms = details::permMFromN(numFragConns, reactionConns.count());

  for (const auto &perm : perms) {
    connPerms.emplace_back();
    // Copy the fragments and set the isotope numbers according to this
    // permutation.
    for (const auto &f : molFrags) {
      connPerms.back().emplace_back(new RWMol(*f));
      boost::dynamic_bitset<> atomDone(f->getNumAtoms());
      for (auto atom : connPerms.back().back()->atoms()) {
        if (!atom->getAtomicNum()) {
          for (size_t i = 0; i < perm.size(); ++i) {
            if (!atomDone[atom->getIdx()] && atom->getIsotope() == i + 1) {
              atom->setIsotope(perm[i] + 1);
              if (atom->hasQuery()) {
                atom->setQuery(makeAtomTypeQuery(0, false));
                atom->expandQuery(makeAtomIsotopeQuery(perm[i] + 1));
              }
              atomDone[atom->getIdx()] = true;
            }
          }
        }
      }
    }
  }

  return connPerms;
}

void buildConnectorRegions(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &connRegs,
    std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> &connRegFPs) {
  for (const auto &frag : molFrags) {
    auto fragConnRegs = getConnRegion(*frag);
    if (!fragConnRegs) {
      // There were no connector atoms.
      continue;
    }
    std::vector<std::unique_ptr<ROMol>> splitConnRegs;
    MolOps::getMolFrags(*fragConnRegs, splitConnRegs, false);
    connRegFPs.emplace_back();
    for (auto &cr : splitConnRegs) {
      connRegFPs.back().emplace_back(PatternFingerprintMol(*cr));
    }
    connRegs.push_back(std::move(splitConnRegs));
  }
}

// Return true if all the fragments have a connector region that matches
// something in the reaction, false otherwise.
bool checkConnectorRegions(
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<std::vector<std::unique_ptr<ROMol>>> &connRegs,
    const std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>>
        &connRegFPs) {
  const auto &rxnConnRegs = reaction->getConnectorRegions();
  const auto &rxnConnRegsFP = reaction->getConnRegFP();
  MatchVectType dontCare;
  for (size_t i = 0; i < connRegFPs.size(); ++i) {
    bool connRegFound = false;
    for (size_t j = 0; j < connRegFPs[i].size(); ++j) {
      if (AllProbeBitsMatch(*connRegFPs[i][j], *rxnConnRegsFP)) {
        for (const auto &rxncr : rxnConnRegs) {
          if (SubstructMatch(*rxncr, *connRegs[i][j], dontCare)) {
            connRegFound = true;
            break;
          }
        }
      }
      if (connRegFound) {
        break;
      }
    }
    if (!connRegFound) {
      return false;
    }
  }
  return true;
}

// Take the molFrags and flag those synthons that have pattern fingerprints
// where all the bits match with the fragment.  The pattern fingerprints are
// insensitive to isotope numbers, so this can be done on the initial
// fragmentation, without generating the combinations of connector numbers.
// Matches the pattFPs with the synthon sets in the order synthonOrder, but
// returns the bitsets in the original order.
std::vector<boost::dynamic_bitset<>> screenSynthonsWithFPs(
    const std::vector<std::unique_ptr<ExplicitBitVect>> &pattFPs,
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> passedFPs;
  for (const auto &synthonSet : reaction->getSynthons()) {
    passedFPs.emplace_back(synthonSet.size());
  }

  boost::dynamic_bitset<> fragsMatched(synthonOrder.size());
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    const auto &synthonSet = reaction->getSynthons()[synthonOrder[i]];
    for (size_t j = 0; j < synthonSet.size(); ++j) {
      auto &synthon = synthonSet[j];
      if (AllProbeBitsMatch(*pattFPs[i], *synthon->getPattFP())) {
        passedFPs[synthonOrder[i]][j] = true;
        fragsMatched[i] = true;
      }
    }
    // If nothing matched this fragment, the whole thing's a bust.
    if (!fragsMatched[i]) {
      break;
    }
  }
  // If all the fragments had a match, these results are valid.
  if (fragsMatched.count() != fragsMatched.size()) {
    for (size_t i = 0; i < passedFPs.size(); ++i) {
      passedFPs[i].reset();
    }
  }

  return passedFPs;
}

// Take the fragged mol and flag all those synthons that have a fragment as
// a substructure match.  Only do this for those synthons that have already
// passed previous screening, and are flagged as such in passedScreens.
std::vector<boost::dynamic_bitset<>> getHitSynthons(
    std::vector<std::unique_ptr<RWMol>> &molFrags,
    const std::vector<boost::dynamic_bitset<>> &passedScreens,
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  RDKit::MatchVectType dontCare;
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  for (const auto &synthonSet : reaction->getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }

  // The tests must be applied for all permutations of synthon list against
  // fragment.
  auto synthonOrders =
      details::permMFromN(molFrags.size(), reaction->getSynthons().size());

  boost::dynamic_bitset<> fragsMatched(synthonOrder.size());
  // Match the fragment to the synthon set in this order.
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    const auto &synthonsSet = reaction->getSynthons()[synthonOrder[i]];
    const auto &passedScreensSet = passedScreens[synthonOrder[i]];
    for (size_t j = 0; j < synthonsSet.size(); ++j) {
      if (passedScreensSet[j]) {
        auto &synthon = synthonsSet[j];
        if (SubstructMatch(*synthon->getMol(), *molFrags[i], dontCare)) {
          synthonsToUse[synthonOrder[i]][j] = true;
          fragsMatched[i] = true;
        }
      }
    }
    // if the fragment didn't match anything, the whole thing's a bust.
    if (!fragsMatched[i]) {
      synthonsToUse.clear();
      return synthonsToUse;
    }
  }
  // If all bits in one of the bitsets is unset, it means that nothing matched
  // that synthon.  If at least one of the bitsets has a set bit, all products
  // incorporating the synthon with no bits set must match the query so
  // should be used because the query matches products that don't incorporate
  // anything from 1 of the synthon lists.  For example, if the synthons are
  // [1*]Nc1c([2*])cccc1 and [1*]=CC=C[2*] and the query is c1ccccc1.
  bool someSet = std::any_of(
      synthonsToUse.begin(), synthonsToUse.end(),
      [](const boost::dynamic_bitset<> &bs) -> bool { return bs.any(); });
  if (someSet) {
    for (auto &rtu : synthonsToUse) {
      if (!rtu.count()) {
        rtu.set();
      }
    }
  }
  return synthonsToUse;
}

}  // namespace

std::vector<SynthonSpaceHitSet> SynthonSpaceSubstructureSearcher::searchFragSet(
    std::vector<std::unique_ptr<ROMol>> &fragSet) const {
  std::vector<SynthonSpaceHitSet> results;

  details::fixAromaticRingSplits(fragSet);
  auto pattFPs = makePatternFPs(fragSet);
  std::vector<std::vector<std::unique_ptr<ROMol>>> connRegs;
  std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> connRegFPs;
  std::vector<int> numFragConns;
  numFragConns.reserve(fragSet.size());
  for (const auto &frag : fragSet) {
    numFragConns.push_back(details::countConnections(MolToSmiles(*frag)));
  }

  auto conns = getConnectorPattern(fragSet);
  for (auto &it : getSpace().getReactions()) {
    auto &reaction = it.second;
    // It can't be a hit if the number of fragments is more than the number
    // of synthon sets because some of the molecule won't be matched in any
    // of the potential products.  It can be less, in which case the unused
    // synthon set will be used completely, possibly resulting in a large
    // number of hits.
    if (fragSet.size() > reaction->getSynthons().size()) {
      continue;
    }

    // Check that all the frags have a connector region that matches something
    // in this reaction set.  Skip if not.
    if (connRegs.empty()) {
      buildConnectorRegions(fragSet, connRegs, connRegFPs);
    }
    if (!checkConnectorRegions(reaction, connRegs, connRegFPs)) {
      continue;
    }

    // Select only the synthons that have fingerprints that are a superset
    // of the fragment fingerprints.
    // Need to try all combinations of synthon orders.
    auto synthonOrders =
        details::permMFromN(pattFPs.size(), reaction->getSynthons().size());
    for (const auto &so : synthonOrders) {
      auto passedScreens = screenSynthonsWithFPs(pattFPs, reaction, so);
      // If none of the synthons passed the screens, move right along, nothing
      // to see.
      bool skip = std::all_of(
          passedScreens.begin(), passedScreens.end(),
          [](const boost::dynamic_bitset<> &s) -> bool { return s.none(); });
      if (skip) {
        continue;
      }

      // Get all the possible permutations of connector numbers compatible with
      // the number of synthon sets in this reaction.  So if the
      // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
      // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
      // that might be how they're labelled in the reaction database.
      auto connCombs =
          getConnectorPermutations(fragSet, conns, reaction->getConnectors());

      // Find all synthons that match the fragments with each connector
      // combination.
      for (auto &connComb : connCombs) {
        auto theseSynthons =
            getHitSynthons(connComb, passedScreens, reaction, so);
        if (!theseSynthons.empty()) {
          size_t numHits = std::accumulate(
              theseSynthons.begin(), theseSynthons.end(), 1,
              [](int prevRes, const boost::dynamic_bitset<> &s2) {
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

bool SynthonSpaceSubstructureSearcher::verifyHit(const ROMol &hit) const {
  MatchVectType dontCare;
  return SubstructMatch(hit, getQuery(), dontCare);
}

}  // namespace RDKit::SynthonSpaceSearch