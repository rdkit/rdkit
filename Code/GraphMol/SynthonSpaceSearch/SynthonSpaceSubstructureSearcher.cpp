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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSubstructureSearcher.h>
#include <RDGeneral/ControlCHandler.h>

namespace RDKit::SynthonSpaceSearch {

namespace {

// Re-order fragments into descending order of number of bits set in the
// pattern fp (the "largest fragment heuristic").  The FP with the
// largest number of bits is the most likely to screen out a matching
// synthon set since smaller, less complex fragments are more likely
// to match something, so screen with that first.
void reorderFragments(
    std::vector<std::unique_ptr<ROMol>> &molFrags,
    const std::vector<std::pair<void *, ExplicitBitVect *>> &allPattFPs) {
  std::vector<ExplicitBitVect *> pattFPs;
  pattFPs.reserve(molFrags.size());
  for (const auto &frag : molFrags) {
    std::pair<void *, ExplicitBitVect *> tmp(frag.get(), nullptr);
    const auto it =
        std::lower_bound(allPattFPs.begin(), allPattFPs.end(), tmp,
                         [](const auto &p1, const auto &p2) -> bool {
                           return p1.first > p2.first;
                         });
    pattFPs.push_back(it->second);
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

  // Now put orderedFrags in the same order.
  std::vector<std::unique_ptr<ROMol>> newFrags(molFrags.size());
  for (size_t i = 0; i < fps.size(); ++i) {
    newFrags[i] = std::move(molFrags[fps[i].first]);
  }
  molFrags = std::move(newFrags);
}

// Collect the pattern fps for the fragments.
std::vector<ExplicitBitVect *> gatherPatternFPs(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    const std::vector<std::pair<void *, ExplicitBitVect *>> &allPattFPs) {
  std::vector<ExplicitBitVect *> pattFPs;
  pattFPs.reserve(molFrags.size());
  for (const auto &frag : molFrags) {
    std::pair<void *, ExplicitBitVect *> tmp(frag.get(), nullptr);
    const auto it =
        std::lower_bound(allPattFPs.begin(), allPattFPs.end(), tmp,
                         [](const auto &p1, const auto &p2) -> bool {
                           return p1.first > p2.first;
                         });
    pattFPs.push_back(it->second);
  }
  return pattFPs;
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
            !SubstructMatch(*rxnConnRegs[k], *connRegs[i][j]).empty()) {
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
    const std::vector<std::unique_ptr<ExplicitBitVect>> &pattFPs,
    const SynthonSet &reaction, const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> passedFPs;
  for (const auto &synthonSet : reaction.getSynthons()) {
    passedFPs.emplace_back(synthonSet.size());
  }
  boost::dynamic_bitset<> fragsMatched(synthonOrder.size());
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    const auto &synthonSet = reaction.getSynthons()[synthonOrder[i]];

    for (size_t j = 0; j < synthonSet.size(); ++j) {
      std::cout << i << " vs " << j << " : " << synthonSet[j].first << " "
                << synthonSet[j].second->getSmiles() << " "
                << synthonSet[j].second->getPattFP()->getNumBits() << ":"
                << synthonSet[j].second->getPattFP()->getNumOnBits() << " vs "
                << pattFPs[i]->getNumBits() << ":" << pattFPs[i]->getNumOnBits()
                << std::endl;
      if (auto &synthon = synthonSet[j].second;
          pattFPs[i]->getNumOnBits() <= synthon->getPattFP()->getNumOnBits() &&
          AllProbeBitsMatch(*pattFPs[i], *synthon->getPattFP())) {
        passedFPs[synthonOrder[i]][j] = true;
        fragsMatched[i] = true;
      }
    }
    // If nothing matched this fragment, the whole thing's a bust.
    if (!fragsMatched[i]) {
      std::cout << "nothing matched for " << i << std::endl;
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
  std::cout << "\ngetHitSynthons with :: ";
  for (const auto &f : molFrags) {
    std::cout << MolToSmiles(*f) << " ";
  }
  std::cout << std::endl;
  std::cout << "synthon order : ";
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    std::cout << synthonOrder[i] << " ";
  }
  std::cout << std::endl;
#if 0
  std::unique_ptr<ROMol> joinedFrag;
  std::vector<ROMol *> molFragPtrs;
  if (molFrags.size() == synthonOrder.size() + 1) {
    // This is the case where there's a ring closure in the reaction, so
    // 2 of the synthons will be bidentate, but the query has only a part-ring
    // trying to match it. O=c1ncnc([c])c1[c] was ground zero for this.  It
    // would come out in the wash, but at the expense of selecting synthons
    // that can't match the query (only 1 of the connection points would be
    // checked).  So try and deal with it.
    // If there is more than 1 ring closure in a reaction it's really
    // complicated, it's likely to be very unusual, and it should also work
    // without special measures.
    // Find the query fragment with 2 dummy atoms, and combine
    // the corresponding 1-dummy-atom fragments into a single query.
    // So if there were 3 fragments [1*]c1ncnc(=O)c1[2*] [1*]c [2*]c
    // arrive at only 2, [1*]c1ncnc(=O)c1[2*] [1*]c.[2*]c
    std::cout << "Wahey, it must be a bidentate one" << std::endl;
    auto connPatts = details::getConnectorPatterns(molFrags);
    int twoConns = -1;
    for (size_t i = 0; i < connPatts.size(); ++i) {
      std::cout << i << " :: " << connPatts[i] << std::endl;
      if (connPatts[i].count() == 2) {
        twoConns = i;
      }
    }
    std::cout << twoConns << std::endl;
    if (twoConns != -1) {
      int firstFrag = -1;
      for (size_t i = 0; i < connPatts.size(); ++i) {
        if (i == twoConns) {
          molFragPtrs.push_back(molFrags[i].get());
          continue;
        }
        std::cout << "frag " << i << " :: " << connPatts[i] << " "
                  << MolToSmiles(*molFrags[i]) << std::endl;
        if ((connPatts[i] & connPatts[twoConns]).count()) {
          if (!joinedFrag) {
            joinedFrag.reset(new ROMol(*molFrags[i]));
            firstFrag = i;
            molFragPtrs.push_back(nullptr);
          } else {
            joinedFrag.reset(combineMols(*joinedFrag, *molFrags[i]));
          }
        } else {
          molFragPtrs.push_back(molFrags[i].get());
        }
      }
      if (joinedFrag) {
        molFragPtrs[firstFrag] = joinedFrag.get();
      }
    }
  } else {
    for (const auto &mol : molFrags) {
      molFragPtrs.push_back(mol.get());
    }
  }
  for (size_t i = 0; i < molFragPtrs.size(); ++i) {
    std::cout << i << " : " << MolToSmiles(*molFragPtrs[i]) << std::endl;
  }
#endif
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  for (size_t i = 0; i < reaction.getSynthons().size(); ++i) {
    std::cout << i << " : " << reaction.getSynthons()[i].size() << " : "
              << reaction.getSynthonConnectorPatterns()[i] << std::endl;
  }

  // Match the fragment to the synthon set in this order.
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    const auto &synthonsSet = reaction.getSynthons()[synthonOrder[i]];
    const auto &passedScreensSet = passedScreens[synthonOrder[i]];
    bool fragMatched = false;
    std::cout << "trying synthon set : ";
    for (size_t j = 0; j < synthonsSet.size(); ++j) {
      std::cout << synthonsSet[j].first << " ";
    }
    std::cout << std::endl;
    for (size_t j = 0; j < synthonsSet.size(); ++j) {
      if (passedScreensSet[j]) {
        if (const auto &[id, synthon] = synthonsSet[j];
            !SubstructMatch(*synthon->getSearchMol(), *molFrags[i]).empty()) {
          std::cout << id << " : " << synthon->getSmiles() << " matched frag "
                    << i << " : " << MolToSmiles(*molFrags[i]) << std::endl;
          synthonsToUse[synthonOrder[i]][j] = true;
          fragMatched = true;
        }
      }
    }
    // if the fragment didn't match anything, the whole thing's a bust.
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
  std::cout << "returning matching synthons " << std::endl;
  for (const auto &rs : retSynthons) {
    for (const auto i : rs) {
      std::cout << i << "  ";
    }
    std::cout << std::endl;
  }
  return retSynthons;
}

}  // namespace

unsigned int SynthonSpaceSubstructureSearcher::getNumQueryFragmentsRequired() {
  auto maxNumSynthons = getSpace().getMaxNumSynthons();
  auto maxNumConns = getSpace().getMaxNumConnectors();
  return std::max(maxNumSynthons, maxNumConns + 1);
}

void SynthonSpaceSubstructureSearcher::extraSearchSetup(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) {
  bool cancelled = false;
  auto fragSmiToFrag = details::mapFragsBySmiles(fragSets, cancelled);
  if (cancelled) {
    return;
  }
  // Now generate the pattern fingerprints for the fragments.
  const auto pattFPSize = getSpace().getPatternFPSize();
  d_pattFPsPool.resize(fragSmiToFrag.size());
  d_connRegsPool.resize(fragSmiToFrag.size());
  d_connRegSmisPool.resize(fragSmiToFrag.size());
  d_connRegFPsPool.resize(fragSmiToFrag.size());
  unsigned int fragNum = 0;
  bool saidSomething = false;
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    if (ControlCHandler::getGotSignal()) {
      return;
    }
    // For the fingerprints, ring info is required, and this needs to be
    // done for every copy of the fragment. We also need to fix any
    // query atoms.
    for (const auto &frag : frags) {
      unsigned int otf;
      sanitizeMol(*static_cast<RWMol *>(frag), otf, MolOps::SANITIZE_SYMMRINGS);

      // Query atoms may define the environment of the fragment (via recursive
      // SMARTS, for example) that a potentially matching synthon may not
      // have, so they need to be made more generic.  For example, if the
      // query is
      // [$(c1ccccc1),$(c1ccncc1),$(c1cnccc1)]C(=O)N1[C&!$(CC(=O))]CCC1
      // and it's split into [$(c1ccccc1),$(c1ccncc1),$(c1cnccc1)]C(=O)[1*]
      // that won't match the synthon [1*]c([*3])C(=O)[2*] even though the
      // synthon could be built into a ring that matches the query if it is
      // combined with a synthon such as [*1]ccccc[*3].  The downside is that
      // the fragment will become less discriminating leading to more false
      // positives in the initial screenout, hence the warning.
      if (details::removeQueryAtoms(*static_cast<RWMol *>(frag)) &&
          !saidSomething) {
        saidSomething = true;
        BOOST_LOG(rdWarningLog) << "Complex queries can be slow." << std::endl;
      }
    }
    d_pattFPsPool[fragNum] = std::unique_ptr<ExplicitBitVect>(
        PatternFingerprintMol(*frags.front(), pattFPSize));
    if (MolToSmiles(*frags.front()) == "[1*]n1cnc([2*])c([3*])c1=O") {
      for (const auto &f : frags) {
        auto fp = std::unique_ptr<ExplicitBitVect>(
            PatternFingerprintMol(*f, pattFPSize));
        std::cout << fragNum << " : " << MolToSmiles(*f) << " : "
                  << fp->getNumOnBits() << std::endl;
        f->debugMol(std::cout);
        std::cout << "Direct SMILES" << std::endl;
        auto m = "[1*]n1cnc([2*])c([3*])c1=O"_smiles;
        auto fp1 = std::unique_ptr<ExplicitBitVect>(
            PatternFingerprintMol(*m, pattFPSize));
        std::cout << fragNum << " : " << MolToSmiles(*m) << " : "
                  << fp1->getNumOnBits() << std::endl;
        m->debugMol(std::cout);
        std::cout << "MATCH = " << AllProbeBitsMatch(*fp, *fp1) << std::endl;
      }
    }
    // std::cout << fragNum << " : " << MolToSmiles(*frags.front()) << " : "
    //           << d_pattFPsPool[fragNum]->getNumOnBits() << std::endl;
    if (auto fragConnRegs = details::buildConnRegion(*frags.front());
        fragConnRegs) {
      MolOps::getMolFrags(*fragConnRegs, d_connRegsPool[fragNum], false);
      d_connRegSmisPool[fragNum].reserve(d_connRegsPool[fragNum].size());
      d_connRegFPsPool[fragNum].reserve(d_connRegsPool[fragNum].size());
      for (auto &cr : d_connRegsPool[fragNum]) {
        d_connRegFPsPool[fragNum].emplace_back(
            PatternFingerprintMol(*cr, PATT_FP_NUM_BITS));
        d_connRegSmisPool[fragNum].emplace_back(MolToSmiles(*cr));
      }
    }
    fragNum++;
  }
  // Now use the pooled info to populate the vectors for each fragSet.
  fragNum = 0;
  d_pattFPs.reserve(fragSmiToFrag.size());
  d_connRegs.reserve(fragSmiToFrag.size());
  d_connRegSmis.reserve(fragSmiToFrag.size());
  d_connRegFPs.reserve(fragSmiToFrag.size());
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    for (auto &frag : frags) {
      d_pattFPs.emplace_back(frag, d_pattFPsPool[fragNum].get());
      d_connRegs.emplace_back(frag, &d_connRegsPool[fragNum]);
      d_connRegSmis.emplace_back(frag, &d_connRegSmisPool[fragNum]);
      d_connRegFPs.emplace_back(frag, &d_connRegFPsPool[fragNum]);
    }
    ++fragNum;
  }
  std::sort(d_pattFPs.begin(), d_pattFPs.end(),
            [](const auto &p1, const auto &p2) -> bool {
              return p1.first > p2.first;
            });
  std::sort(d_connRegs.begin(), d_connRegs.end(),
            [](const auto &p1, const auto &p2) -> bool {
              return p1.first > p2.first;
            });
  std::sort(d_connRegSmis.begin(), d_connRegSmis.end(),
            [](const auto &p1, const auto &p2) -> bool {
              return p1.first > p2.first;
            });
  std::sort(d_connRegFPs.begin(), d_connRegFPs.end(),
            [](const auto &p1, const auto &p2) -> bool {
              return p1.first > p2.first;
            });

  // Now apply the largest fragment heuristic to the fragments
  for (auto &fs : fragSets) {
    reorderFragments(fs, d_pattFPs);
  }
}

namespace {
std::vector<std::unique_ptr<ROMol>> mergeFragments(
    size_t i, size_t one, size_t two, const std::vector<size_t> &others,
    const std::vector<std::unique_ptr<ROMol>> &fragSet) {
  std::vector<std::unique_ptr<ROMol>> retFrags;
  retFrags.emplace_back(std::make_unique<ROMol>(*fragSet[i]));
  retFrags.emplace_back(
      std::unique_ptr<ROMol>(combineMols(*fragSet[one], *fragSet[two])));
  for (auto o : others) {
    if (o != one && o != two) {
      retFrags.emplace_back(std::make_unique<ROMol>(*fragSet[o]));
    }
  }
  return retFrags;
}
std::vector<std::unique_ptr<ROMol>> mergeFragments(
    size_t i, size_t one, size_t two, size_t three, size_t four,
    const std::vector<std::unique_ptr<ROMol>> &fragSet) {
  std::vector<std::unique_ptr<ROMol>> retFrags;
  retFrags.emplace_back(std::make_unique<ROMol>(*fragSet[i]));
  retFrags.emplace_back(
      std::unique_ptr<ROMol>(combineMols(*fragSet[one], *fragSet[two])));
  retFrags.emplace_back(
      std::unique_ptr<ROMol>(combineMols(*fragSet[three], *fragSet[four])));
  return retFrags;
}

std::vector<std::unique_ptr<ExplicitBitVect>> mergeBitVects(
    size_t i, size_t one, size_t two, const std::vector<size_t> &others,
    const std::vector<ExplicitBitVect *> &bitVects) {
  std::vector<std::unique_ptr<ExplicitBitVect>> retVects;
  retVects.emplace_back(std::make_unique<ExplicitBitVect>(*bitVects[i]));
  retVects.emplace_back(
      std::make_unique<ExplicitBitVect>(*bitVects[one] | *bitVects[two]));
  for (auto o : others) {
    if (o != one && o != two) {
      retVects.emplace_back(std::make_unique<ExplicitBitVect>(*bitVects[o]));
    }
  }
  return retVects;
}
std::vector<std::unique_ptr<ExplicitBitVect>> mergeBitVects(
    size_t i, size_t one, size_t two, size_t three, size_t four,
    const std::vector<ExplicitBitVect *> &bitVects) {
  std::vector<std::unique_ptr<ExplicitBitVect>> retVects;
  retVects.emplace_back(std::make_unique<ExplicitBitVect>(*bitVects[i]));
  retVects.emplace_back(
      std::make_unique<ExplicitBitVect>(*bitVects[one] | *bitVects[two]));
  retVects.emplace_back(
      std::make_unique<ExplicitBitVect>(*bitVects[three] | *bitVects[four]));
  return retVects;
}

// If this is a ring-forming reaction, merge all possible pairs or fragments,
// creating multiple possible fragment sets
std::pair<std::vector<std::vector<std::unique_ptr<ROMol>>>,
          std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>>>
mergeRingFormingFrags(const std::vector<std::unique_ptr<ROMol>> &fragSet,
                      const std::vector<ExplicitBitVect *> &pattFPs,
                      const SynthonSet &reaction) {
  std::vector<std::vector<std::unique_ptr<ROMol>>> fragSetCps;
  std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> pattFPsMerged;
  // A copy of the input set will always be needed.
  fragSetCps.push_back(std::vector<std::unique_ptr<ROMol>>());
  pattFPsMerged.push_back(std::vector<std::unique_ptr<ExplicitBitVect>>());
  for (unsigned int i = 0; i < fragSet.size(); ++i) {
    fragSetCps.back().emplace_back(std::make_unique<ROMol>(*fragSet[i]));
    pattFPsMerged.back().emplace_back(
        std::make_unique<ExplicitBitVect>(*pattFPs[i]));
  }
  // If there are only 1 or 2 query fragments, no complications arise
  if (fragSet.size() < 3) {
    return std::make_pair(std::move(fragSetCps), std::move(pattFPsMerged));
  }

  if (reaction.getNumRingFormers() > 0) {
    // This is the case where there's 1 or 2 ring closures in the reaction, so
    // either 2 of the synthons will be bidentate (1 ring closure), or
    // 1 synthon will be 4-dentate, and 2 will be bidentate.
    // We need to allow for a query that only
    // has a part-ring trying to match it. O=c1ncnc([c])c1[c] was ground
    // zero for this.  For a 1 ring-former, take each bi-dentate fragment
    // and make all combinations of 1 pair from the rest and then the rest.
    // Within the current constraint of only 4 connectors being allowed,
    // 2 ring forming reactions must mean a single 4-connector fragment
    // and 4 1-connector fragments

    std::cout << "It's a ring former" << std::endl;
    const auto &connPatts = details::getConnectorPatterns(fragSet);
    for (size_t i = 0; i < connPatts.size(); ++i) {
      if (connPatts[i].count() > 1) {
        std::vector<size_t> others;
        for (size_t j = 0; j < connPatts.size(); ++j) {
          if (i != j) {
            others.push_back(j);
          }
        }
        std::cout << "doubler : " << i << "  " << connPatts[i]
                  << "  others :: ";
        for (auto o : others) {
          std::cout << o << "  " << connPatts[o] << "  ";
        }
        std::cout << std::endl;
        // If there are only 2 other fragments, it's easy
        if (others.size() == 2) {
          fragSetCps.push_back(
              mergeFragments(i, others[0], others[1], others, fragSet));
          pattFPsMerged.push_back(
              mergeBitVects(i, others[0], others[1], others, pattFPs));
        } else if (others.size() == 3) {
          // Merge each possible pair, leaving the other one
          fragSetCps.push_back(
              mergeFragments(i, others[0], others[1], others, fragSet));
          fragSetCps.push_back(
              mergeFragments(i, others[0], others[2], others, fragSet));
          fragSetCps.push_back(
              mergeFragments(i, others[1], others[2], others, fragSet));
          pattFPsMerged.push_back(
              mergeBitVects(i, others[0], others[1], others, pattFPs));
          pattFPsMerged.push_back(
              mergeBitVects(i, others[0], others[2], others, pattFPs));
          pattFPsMerged.push_back(
              mergeBitVects(i, others[1], others[2], others, pattFPs));
        } else if (others.size() == 4) {
          // merge all combinations of pairs
          fragSetCps.push_back(mergeFragments(i, others[0], others[1],
                                              others[2], others[3], fragSet));
          fragSetCps.push_back(mergeFragments(i, others[0], others[2],
                                              others[1], others[3], fragSet));
          fragSetCps.push_back(mergeFragments(i, others[1], others[2],
                                              others[0], others[3], fragSet));
          pattFPsMerged.push_back(mergeBitVects(i, others[0], others[1],
                                                others[2], others[3], pattFPs));
          pattFPsMerged.push_back(mergeBitVects(i, others[0], others[2],
                                                others[1], others[3], pattFPs));
          pattFPsMerged.push_back(mergeBitVects(i, others[1], others[2],
                                                others[0], others[3], pattFPs));
        }
      }
    }
  }
  std::cout << "fragSetCps" << std::endl;
  for (const auto &fragSetCp : fragSetCps) {
    for (const auto &frag : fragSetCp) {
      std::cout << MolToSmiles(*frag) << " ";
    }
    std::cout << std::endl;
  }
  for (const auto &pattFP : pattFPsMerged) {
    for (const auto &p : pattFP) {
      std::cout << p->getNumBits() << ", " << p->getNumOnBits() << " ";
    }
    std::cout << std::endl;
  }
  return std::make_pair(std::move(fragSetCps), std::move(pattFPsMerged));
}
}  // namespace

std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceSubstructureSearcher::searchFragSet(
    const std::vector<std::unique_ptr<ROMol>> &fragSet,
    const SynthonSet &reaction) const {
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> results;
  const auto pattFPs = gatherPatternFPs(fragSet, d_pattFPs);
  std::vector<int> numFragConns;
  numFragConns.reserve(fragSet.size());
  for (const auto &frag : fragSet) {
    numFragConns.push_back(details::countConnections(*frag));
  }

  std::cout << "searching fragset : " << fragSet.size() << " vs "
            << reaction.getConnectors().count() + 1 << " and "
            << reaction.getSynthons().size() << std::endl;
  for (const auto &f : fragSet) {
    std::cout << MolToSmiles(*f) << " ";
  }
  std::cout << std::endl;

  const auto conns = details::getConnectorPattern(fragSet);
  // It can't be a hit if the number of fragments is more than 1 plus the
  // number of bonds being formed in the reaction.  It can be less, in which
  // case the unused synthon set will be used completely, possibly resulting
  // in a large number of hits.
  if (fragSet.size() > reaction.getConnectors().count() + 1) {
    return results;
  }
  if (fragSet.size() > reaction.getSynthons().size()) {
    std::cout << "AWOOGA - here we go" << std::endl;
  }
  auto [fragSetCps, mergedPattFPs] =
      mergeRingFormingFrags(fragSet, pattFPs, reaction);

  for (size_t fragNum = 0; fragNum < fragSetCps.size(); ++fragNum) {
    const auto &fragSetCp = fragSetCps[fragNum];
    const auto &mergedPattFp = mergedPattFPs[fragNum];
    std::cout << "next fragSetCp :: ";
    for (const auto &frag : fragSetCp) {
      std::cout << MolToSmiles(*frag) << " ";
    }
    std::cout << std::endl;
    for (const auto &pattFP : mergedPattFp) {
      std::cout << pattFP->getNumBits() << "," << pattFP->getNumOnBits()
                << "  ";
    }
    std::cout << std::endl;
    // Check that all the frags have a connector region that matches something
    // in this reaction set.  Skip if not.
    std::vector<std::vector<ROMol *>> connRegs;
    std::vector<std::vector<const std::string *>> connRegSmis;
    std::vector<std::vector<ExplicitBitVect *>> connRegFPs;
    getConnectorRegions(fragSetCp, connRegs, connRegSmis, connRegFPs);
    for (const auto &cs : connRegSmis) {
      for (const auto &c : cs) {
        std::cout << *c << " ";
      }
      std::cout << std::endl;
    }
    if (!checkConnectorRegions(reaction, connRegs, connRegSmis, connRegFPs)) {
      std::cout << "return on connector regions" << std::endl;
      continue;
    }

    // Get all the possible permutations of connector numbers compatible with
    // the number of synthon sets in this reaction.  So if the
    // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
    // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
    // that might be how they're labelled in the reaction database.

    const auto connCombs = details::getConnectorPermutations(
        fragSetCp, conns, reaction.getConnectors());

    // Select only the synthons that have fingerprints that are a superset
    // of the fragment fingerprints.
    // Need to try all combinations of synthon orders.
    const auto synthonOrders =
        details::permMFromN(fragSetCp.size(), reaction.getSynthons().size());
    for (const auto &so : synthonOrders) {
      std::cout << "Doing Synthon order ";
      for (auto s : so) {
        std::cout << s << " ";
      }
      std::cout << std::endl;
      auto passedScreens = screenSynthonsWithFPs(mergedPattFp, reaction, so);
      // If none of the synthons passed the screens, move right along, nothing
      // to see.
      const bool skip = std::all_of(
          passedScreens.begin(), passedScreens.end(),
          [](const boost::dynamic_bitset<> &s) -> bool { return s.none(); });
      if (skip) {
        std::cout << "Skipping because fps" << std::endl;
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
          std::cout << "number of hit synthon sets : " << theseSynthons.size()
                    << std::endl;
          std::unique_ptr<SynthonSpaceHitSet> hs(
              new SynthonSpaceHitSet(reaction, theseSynthons, fragSet));
          for (const auto &stu : hs->synthonsToUse) {
            for (const auto &sp : stu) {
              std::cout << sp.first << " " << sp.second->getSmiles() << " :: ";
            }
            std::cout << std::endl;
          }
          if (hs->numHits) {
            results.push_back(std::move(hs));
          }
        }
      }
    }
  }
  return results;
}

bool SynthonSpaceSubstructureSearcher::verifyHit(ROMol &hit) const {
  if (!SynthonSpaceSearcher::verifyHit(hit)) {
    return false;
  }
  return !SubstructMatch(hit, getQuery(), d_matchParams).empty();
}

void SynthonSpaceSubstructureSearcher::getConnectorRegions(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    std::vector<std::vector<ROMol *>> &connRegs,
    std::vector<std::vector<const std::string *>> &connRegSmis,
    std::vector<std::vector<ExplicitBitVect *>> &connRegFPs) const {
  for (const auto &frag : molFrags) {
    std::pair<void *, void *> tmp(frag.get(), nullptr);
    const auto it1 =
        std::lower_bound(d_connRegs.begin(), d_connRegs.end(), tmp,
                         [](const auto &p1, const auto &p2) -> bool {
                           return p1.first > p2.first;
                         });
    if (it1->first == tmp.first && !it1->second->empty()) {
      connRegs.push_back(std::vector<ROMol *>());
      std::transform(it1->second->begin(), it1->second->end(),
                     std::back_inserter(connRegs.back()),
                     [](const auto &m) -> ROMol * { return m.get(); });
    }

    const auto it2 =
        std::lower_bound(d_connRegSmis.begin(), d_connRegSmis.end(), tmp,
                         [](const auto &p1, const auto &p2) -> bool {
                           return p1.first > p2.first;
                         });
    if (it2->first == tmp.first && !it2->second->empty()) {
      connRegSmis.push_back(std::vector<const std::string *>());
      std::transform(it2->second->begin(), it2->second->end(),
                     std::back_inserter(connRegSmis.back()),
                     [](const auto &s) -> const std::string * { return &s; });
    }

    const auto it3 =
        std::lower_bound(d_connRegFPs.begin(), d_connRegFPs.end(), tmp,
                         [](const auto &p1, const auto &p2) -> bool {
                           return p1.first > p2.first;
                         });
    if (it3->first == tmp.first && !it3->second->empty()) {
      connRegFPs.push_back(std::vector<ExplicitBitVect *>());
      std::transform(
          it3->second->begin(), it3->second->end(),
          std::back_inserter(connRegFPs.back()),
          [](const auto &v) -> ExplicitBitVect * { return v.get(); });
    }
  }
}

}  // namespace RDKit::SynthonSpaceSearch