//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/MolBundle.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
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
std::vector<ExplicitBitVect *> makePatternFPs(
    std::vector<std::unique_ptr<ROMol>> &molFrags,
    const std::map<void *, std::unique_ptr<ExplicitBitVect>> &allPattFPs) {
  std::vector<ExplicitBitVect *> pattFPs;
  pattFPs.reserve(molFrags.size());
  for (const auto &frag : molFrags) {
    const auto it = allPattFPs.find(frag.get());
    pattFPs.push_back(it->second.get());
  }
  // Sort by descending number of bits set.
  std::vector<std::pair<size_t, ExplicitBitVect *>> fps(pattFPs.size());
  for (size_t i = 0; i < pattFPs.size(); ++i) {
    fps[i] = std::make_pair(i, pattFPs[i]);
  }
  std::sort(fps.begin(), fps.end(),
            [](const std::pair<size_t, ExplicitBitVect *> &fp1,
               const std::pair<size_t, ExplicitBitVect *> &fp2) -> bool {
              return fp1.second->getNumOnBits() > fp2.second->getNumOnBits();
            });

  // Now put molFrags in the same order.
  std::vector<std::unique_ptr<ROMol>> newFrags(molFrags.size());
  std::vector<ExplicitBitVect *> retFPs(molFrags.size());
  for (size_t i = 0; i < fps.size(); ++i) {
    newFrags[i] = std::move(molFrags[fps[i].first]);
    retFPs[i] = pattFPs[fps[i].first];
  }
  molFrags = std::move(newFrags);
  return retFPs;
}

// Return true if all the fragments have a connector region that matches
// something in the reaction, false otherwise.
bool checkConnectorRegions(
    const SynthonSet &reaction,
    const std::vector<std::vector<ROMol *>> &connRegs,
    const std::vector<std::vector<const std::string *>> &connRegSmis,
    const std::vector<std::vector<ExplicitBitVect *>> &connRegFPs) {
  const auto &rxnConnRegs = reaction.getConnectorRegions();
  const auto &rxnConnRegSmis = reaction.getConnectorRegionSmiles();
  const auto &rxnConnRegFPs = reaction.getConnRegFPs();
  MatchVectType dontCare;
  for (size_t i = 0; i < connRegFPs.size(); ++i) {
    bool connRegFound = false;
    for (size_t j = 0; j < connRegFPs[i].size(); ++j) {
      for (size_t k = 0; k < rxnConnRegs.size(); ++k) {
        if (rxnConnRegSmis[k] == *connRegSmis[i][j]) {
          connRegFound = true;
          break;
        }
        if (connRegFPs[i][j]->getNumOnBits() <=
                rxnConnRegFPs[k]->getNumOnBits() &&
            AllProbeBitsMatch(*connRegFPs[i][j], *rxnConnRegFPs[k]) &&
            SubstructMatch(*rxnConnRegs[k], *connRegs[i][j], dontCare)) {
          connRegFound = true;
          break;
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
    const std::vector<ExplicitBitVect *> &pattFPs, const SynthonSet &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> passedFPs;
  for (const auto &synthonSet : reaction.getSynthons()) {
    passedFPs.emplace_back(synthonSet.size());
  }
  boost::dynamic_bitset<> fragsMatched(synthonOrder.size());
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    const auto &synthonSet = reaction.getSynthons()[synthonOrder[i]];

    for (size_t j = 0; j < synthonSet.size(); ++j) {
      if (auto &synthon = synthonSet[j].second;
          pattFPs[i]->getNumOnBits() <= synthon->getPattFP()->getNumOnBits() &&
          AllProbeBitsMatch(*pattFPs[i], *synthon->getPattFP())) {
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
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    const std::vector<boost::dynamic_bitset<>> &passedScreens,
    const SynthonSet &reaction, const std::vector<unsigned int> &synthonOrder) {
  MatchVectType dontCare;
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }

  // The tests must be applied for all permutations of synthon list against
  // fragment.
  auto synthonOrders =
      details::permMFromN(molFrags.size(), reaction.getSynthons().size());

  // Match the fragment to the synthon set in this order.
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    const auto &synthonsSet = reaction.getSynthons()[synthonOrder[i]];
    const auto &passedScreensSet = passedScreens[synthonOrder[i]];
    bool fragMatched = false;
    for (size_t j = 0; j < synthonsSet.size(); ++j) {
      if (passedScreensSet[j]) {
        if (const auto &[id, synthon] = synthonsSet[j];
            SubstructMatch(*synthon->getSearchMol(), *molFrags[i], dontCare)) {
          synthonsToUse[synthonOrder[i]][j] = true;
          fragMatched = true;
        }
      }
    }
    // if the fragment didn't match anything, the whole thing's a bust.
    if (!fragMatched) {
      synthonsToUse.clear();
      return retSynthons;
    }
  }

  // Fill in any synthons where they all didn't match.
  details::expandBitSet(synthonsToUse);
  details::bitSetsToVectors(synthonsToUse, retSynthons);

  // Now sort the selected synthons into ascending order of number of atoms,
  // since smaller molecules are likely to be of more interest.
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

void SynthonSpaceSubstructureSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) {
  bool saidSomething = false;

  const auto pattFPSize = getSpace().getPatternFPSize();
  for (auto &fragSet : fragSets) {
    for (auto &frag : fragSet) {
      // For the fingerprints, ring info is required.
      unsigned int otf;
      sanitizeMol(*static_cast<RWMol *>(frag.get()), otf,
                  MolOps::SANITIZE_SYMMRINGS);

      if (details::removeQueryAtoms(*static_cast<RWMol *>(frag.get())) &&
          !saidSomething) {
        saidSomething = true;
        BOOST_LOG(rdWarningLog) << "Complex queries can be slow." << std::endl;
      }

      d_pattFPs.insert(std::make_pair(
          frag.get(), std::unique_ptr<ExplicitBitVect>(
                          PatternFingerprintMol(*frag, pattFPSize))));

      if (auto fragConnRegs = details::buildConnRegion(*frag); fragConnRegs) {
        std::vector<std::unique_ptr<ROMol>> splitConnRegs;
        MolOps::getMolFrags(*fragConnRegs, splitConnRegs, false);
        std::vector<std::unique_ptr<ExplicitBitVect>> connRegFPs;
        connRegFPs.reserve(splitConnRegs.size());
        std::vector<std::string> connRegSmis;
        connRegSmis.reserve(splitConnRegs.size());
        connRegFPs.reserve(splitConnRegs.size());
        for (auto &cr : splitConnRegs) {
          connRegFPs.emplace_back(PatternFingerprintMol(*cr, PATT_FP_NUM_BITS));
          connRegSmis.emplace_back(MolToSmiles(*cr));
        }
        d_connRegs.insert(std::make_pair(frag.get(), std::move(splitConnRegs)));
        d_connRegSmis.insert(
            std::make_pair(frag.get(), std::move(connRegSmis)));
        d_connRegFPs.insert(std::make_pair(frag.get(), std::move(connRegFPs)));
      }
    }
  }
}

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceSubstructureSearcher::searchFragSet(
    std::vector<std::unique_ptr<ROMol>> &fragSet,
    const SynthonSet &reaction) const {
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> results;

  const auto pattFPs = makePatternFPs(fragSet, d_pattFPs);
  std::vector<int> numFragConns;
  numFragConns.reserve(fragSet.size());
  for (const auto &frag : fragSet) {
    numFragConns.push_back(details::countConnections(*frag));
  }

  const auto conns = details::getConnectorPattern(fragSet);
  // It can't be a hit if the number of fragments is more than the number
  // of synthon sets because some of the molecule won't be matched in any
  // of the potential products.  It can be less, in which case the unused
  // synthon set will be used completely, possibly resulting in a large
  // number of hits.
  if (fragSet.size() > reaction.getSynthons().size()) {
    return results;
  }

  // Check that all the frags have a connector region that matches something
  // in this reaction set.  Skip if not.
  std::vector<std::vector<ROMol *>> connRegs;
  std::vector<std::vector<const std::string *>> connRegSmis;
  std::vector<std::vector<ExplicitBitVect *>> connRegFPs;
  getConnectorRegions(fragSet, connRegs, connRegSmis, connRegFPs);
  if (!checkConnectorRegions(reaction, connRegs, connRegSmis, connRegFPs)) {
    return results;
  }

  // Get all the possible permutations of connector numbers compatible with
  // the number of synthon sets in this reaction.  So if the
  // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
  // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
  // that might be how they're labelled in the reaction database.

  // Need to copy fragSet into a working set, then apply the results to
  // the working copy below.
  std::vector<std::unique_ptr<ROMol>> fragSetCp(fragSet.size());
  for (unsigned int i = 0; i < fragSet.size(); ++i) {
    fragSetCp[i] = std::make_unique<ROMol>(*fragSet[i]);
  }
  const auto connCombs = details::getConnectorPermutations(
      fragSetCp, conns, reaction.getConnectors());

  // Select only the synthons that have fingerprints that are a superset
  // of the fragment fingerprints.
  // Need to try all combinations of synthon orders.
  const auto synthonOrders =
      details::permMFromN(pattFPs.size(), reaction.getSynthons().size());
  for (const auto &so : synthonOrders) {
    auto passedScreens = screenSynthonsWithFPs(pattFPs, reaction, so);
    // If none of the synthons passed the screens, move right along, nothing
    // to see.
    const bool skip = std::all_of(
        passedScreens.begin(), passedScreens.end(),
        [](const boost::dynamic_bitset<> &s) -> bool { return s.none(); });
    if (skip) {
      continue;
    }

    // Find all synthons that match the fragments with each connector
    // combination.
    for (const auto &connComb : connCombs) {
      for (size_t i = 0; i < connComb.size(); ++i) {
        for (const auto &[atom, isotopeNum] : connComb[i]) {
          atom->setIsotope(isotopeNum);
          if (atom->hasQuery()) {
            atom->setQuery(makeAtomTypeQuery(0, false));
            atom->expandQuery(makeAtomIsotopeQuery(isotopeNum));
          }
        }
      }
      auto theseSynthons =
          getHitSynthons(fragSetCp, passedScreens, reaction, so);
      if (!theseSynthons.empty()) {
        std::unique_ptr<SynthonSpaceHitSet> hs(
            new SynthonSpaceHitSet(reaction, theseSynthons, fragSet));
        if (hs->numHits) {
          results.push_back(std::move(hs));
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

void SynthonSpaceSubstructureSearcher::getConnectorRegions(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    std::vector<std::vector<ROMol *>> &connRegs,
    std::vector<std::vector<const std::string *>> &connRegSmis,
    std::vector<std::vector<ExplicitBitVect *>> &connRegFPs) const {
  for (const auto &frag : molFrags) {
    if (const auto it = d_connRegs.find(frag.get()); it != d_connRegs.end()) {
      connRegs.push_back(std::vector<ROMol *>());
      std::transform(
          it->second.begin(), it->second.end(),
          std::back_inserter(connRegs.back()),
          [](const std::unique_ptr<ROMol> &m) -> ROMol * { return m.get(); });
    }
    if (const auto &it = d_connRegSmis.find(frag.get());
        it != d_connRegSmis.end()) {
      connRegSmis.push_back(std::vector<const std::string *>());
      std::transform(
          it->second.begin(), it->second.end(),
          std::back_inserter(connRegSmis.back()),
          [](const std::string &s) -> const std::string * { return &s; });
    }
    if (const auto it = d_connRegFPs.find(frag.get());
        it != d_connRegFPs.end()) {
      connRegFPs.push_back(std::vector<ExplicitBitVect *>());
      std::transform(
          it->second.begin(), it->second.end(),
          std::back_inserter(connRegFPs.back()),
          [](const std::unique_ptr<ExplicitBitVect> &v) -> ExplicitBitVect * {
            return v.get();
          });
    }
  }
}

}  // namespace RDKit::SynthonSpaceSearch