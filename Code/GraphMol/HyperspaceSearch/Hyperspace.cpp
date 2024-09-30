//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <fstream>
#include <regex>
#include <string>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "Hyperspace.h"
#include "HyperspaceSubstructureSearch.h"

namespace RDKit::HyperspaceSSSearch {
Hyperspace::Hyperspace(const std::string &fileName) : d_fileName(fileName) {
  readFile();
  assignConnectorsUsed();
}

std::vector<std::unique_ptr<ROMol>> Hyperspace::search(
    const RDKit::ROMol &query, unsigned int maxBondSplits) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");

  std::vector<std::unique_ptr<ROMol>> results;
  RDKit::MatchVectType dontCare;

  auto fragments = details::splitMolecule(query, maxBondSplits);
  //  for (const auto &frags : fragments) {
  //    std::cout << "****************************" << std::endl;
  //    for (const auto &f : frags) {
  //      std::cout << MolToSmiles(*f) << std::endl;
  //    }
  //  }
  for (auto &fragSet : fragments) {
    auto theseResults = searchFragSet(fragSet);
    // Just for safety, do a final check
    for (auto &hitMol : theseResults) {
      if (SubstructMatch(*hitMol, query, dontCare)) {
        results.emplace_back(std::unique_ptr<ROMol>(hitMol.release()));
      }
    }
  }

  return results;
}

namespace {
inline std::vector<std::string> splitLine(const std::string &str,
                                          const std::regex &regexz) {
  return {std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
          std::sregex_token_iterator()};
}

}  // namespace

void Hyperspace::readFile() {
  std::ifstream ifs(d_fileName);
  // the first line is headers, check they're in the right format,
  // which is tab-separated:
  // SMILES	synton_id	synton#	reaction_id
  // Note that we keep the spelling "synton" from the original paper
  // and example input file, and allow any whitespace rather than just tab.
  std::string firstLine;
  getline(ifs, firstLine);
  //  std::cout << "parsing : " << firstLine << std::endl;
  std::regex regexz("\\s+");

  auto lineParts = splitLine(firstLine, regexz);
  if (lineParts != std::vector<std::string>{"SMILES", "synton_id", "synton#",
                                            "reaction_id"}) {
    throw std::runtime_error("Bad format for hyperspace file " + d_fileName);
  }
  std::string nextLine;
  int lineNum = 1;

  while (getline(ifs, nextLine)) {
    ++lineNum;
    auto nextReag = splitLine(nextLine, regexz);
    if (nextReag.size() != 4) {
      throw std::runtime_error("Bad format for hyperspace file " + d_fileName +
                               " on line " + std::to_string(lineNum));
    }
    //    std::cout << nextReag[0] << " : " << nextReag[1] << std::endl;
    if (auto it = d_reactions.find(nextReag[3]); it == d_reactions.end()) {
      d_reactions.insert(std::make_pair(
          nextReag[3], std::make_unique<ReagentSet>(nextReag[3])));
    }
    auto &currReaction = d_reactions[nextReag[3]];
    size_t synthonId = std::stoi(nextReag[2]);
    if (synthonId >= currReaction->d_reagents.size()) {
      for (size_t i = currReaction->d_reagents.size(); i < synthonId + 1; ++i) {
        currReaction->d_reagents.push_back(
            std::vector<std::unique_ptr<Reagent>>());
      }
    }
    currReaction->d_reagents[synthonId].emplace_back(
        new Reagent(nextReag[0], nextReag[1]));
  }
}

void Hyperspace::assignConnectorsUsed() {
  const static std::vector<std::regex> connRegexs{
      std::regex(R"(\[1\*\])"), std::regex(R"(\[2\*\])"),
      std::regex(R"(\[3\*\])"), std::regex(R"(\[4\*\])")};
  for (auto &it : d_reactions) {
    auto &reaction = it.second;
    reaction->d_connectors.resize(4, false);
    for (auto &reagSet : reaction->d_reagents) {
      for (auto &reag : reagSet) {
        for (size_t i = 0; i < 4; ++i) {
          if (std::regex_search(reag->d_smiles, connRegexs[i])) {
            reaction->d_connectors.set(i);
          }
        }
      }
    }
    //    std::cout << "Connectors : " << reaction->d_id << " : "
    //              << reaction->d_connectors << std::endl;
  }
}

namespace {
// Return a bitset giving the different connector types in this
// molecule.
boost::dynamic_bitset<> getConnectorPattern(
    const std::vector<std::shared_ptr<ROMol>> &fragSet) {
  boost::dynamic_bitset<> conns(4);
  for (const auto &frag : fragSet) {
    for (const auto &a : frag->atoms()) {
      if (!a->getAtomicNum()) {
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
std::vector<std::vector<std::shared_ptr<RWMol>>> getConnectorPermutations(
    const std::vector<std::shared_ptr<ROMol>> &molFrags,
    const boost::dynamic_bitset<> &fragConns,
    const boost::dynamic_bitset<> &reactionConns) {
  std::vector<std::vector<std::shared_ptr<RWMol>>> connPerms;
  //  for (const auto &mf : molFrags) {
  //    std::cout << MolToSmiles(*mf) << " ";
  //  }
  //  std::cout << " :: " << fragConns << " : " << reactionConns << std::endl;
  auto bitsToInts =
      [](const boost::dynamic_bitset<> &bits) -> std::vector<int> {
    std::vector<int> ints;
    for (size_t i = 0; i < bits.size(); ++i) {
      if (bits[i]) {
        ints.push_back(i);
      }
    }
    return ints;
  };
  auto numFragConns = fragConns.count();
  auto rConns = bitsToInts(reactionConns);
  //  std::cout << "Num frag conns : " << numFragConns
  //            << " num reaction conns : " << reactionConns.count() <<
  //            std::endl;
  auto perms = details::permMFromN(numFragConns, reactionConns.count());
  //  std::cout << "PERMS" << std::endl;

  for (const auto &perm : perms) {
    //    std::cout << "Next perm : ";
    //    for (auto &i : perm) {
    //      std::cout << i << " ";
    //    }
    //    std::cout << std::endl;
    connPerms.push_back(std::vector<std::shared_ptr<RWMol>>());
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
    //    std::cout << "Connector permutations" << std::endl;
    //    for (const auto &f : connPerms.back()) {
    //      std::cout << MolToSmiles(*f) << " " << MolToSmarts(*f) << " :: ";
    //    }
    //    std::cout << std::endl;
  }

  return connPerms;
}

// Take the molFrags and flag those reagents that have pattern fingerprints
// where all the bits match with the fragment.  The pattern fingerprints are
// insensitive to isotope numbers, so this can be done on the initial
// fragmentation, without generating the combinations of connector numbers.
std::vector<boost::dynamic_bitset<>> screenReagentsWithFPs(
    const std::vector<std::shared_ptr<ROMol>> &molFrags,
    std::unique_ptr<ReagentSet> &reaction) {
  std::vector<boost::dynamic_bitset<>> reagsToUse;
  std::vector<std::unique_ptr<ExplicitBitVect>> pattFPs;
  for (const auto &frag : molFrags) {
    pattFPs.emplace_back(
        std::unique_ptr<ExplicitBitVect>(PatternFingerprintMol(*frag, 2048)));
  }

  std::vector<boost::dynamic_bitset<>> passedFPs;
  for (const auto &reagSet : reaction->d_reagents) {
    passedFPs.push_back(boost::dynamic_bitset<>(reagSet.size()));
  }

  // Still need to try all combinations of reagent orders.
  auto reagentOrders =
      details::permMFromN(molFrags.size(), reaction->d_reagents.size());
  for (const auto &ro : reagentOrders) {
    //    std::cout << "Reagent orders : ";
    //    for (auto &i : ro) {
    //      std::cout << i << " ";
    //    }
    //    std::cout << std::endl;
    std::vector<boost::dynamic_bitset<>> thisPass;
    for (const auto &reagSet : reaction->d_reagents) {
      thisPass.push_back(boost::dynamic_bitset<>(reagSet.size()));
    }
    boost::dynamic_bitset<> fragsMatched(ro.size());
    for (size_t i = 0; i < ro.size(); ++i) {
      const auto &reagSet = reaction->d_reagents[ro[i]];
      for (size_t j = 0; j < reagSet.size(); ++j) {
        auto &reag = reagSet[j];
        if (AllProbeBitsMatch(*pattFPs[i], *reag->pattFP())) {
          //          std::cout << i << " vs " << j << " :: " <<
          //          MolToSmiles(*molFrags[i])
          //                    << " IS fp match to " <<
          //                    MolToSmiles(*reag->mol())
          //                    << std::endl;
          thisPass[ro[i]][j] = true;
          fragsMatched[i] = true;
        }
      }
    }
    // If all the fragments had a match, these results are valid.
    if (fragsMatched.count() == fragsMatched.size()) {
      for (size_t i = 0; i < passedFPs.size(); ++i) {
        passedFPs[i] |= thisPass[i];
      }
    }
  }

  //  for (const auto &pf : passedFPs) {
  //    std::cout << "passFPs : " << pf << std::endl;
  //  }
  return passedFPs;
}

// Take the fragged mol and flag all those reagents that have a fragment as
// a substructure match.  Only do this for those reagents that have already
// passed previous screening, and are flagged as such in passedScreens.
std::vector<boost::dynamic_bitset<>> getHitReagents(
    std::vector<std::shared_ptr<RWMol>> &molFrags,
    const std::vector<boost::dynamic_bitset<>> &passedScreens,
    std::unique_ptr<ReagentSet> &reaction) {
  RDKit::MatchVectType dontCare;
  std::vector<boost::dynamic_bitset<>> reagsToUse;
  for (const auto &reagSet : reaction->d_reagents) {
    reagsToUse.push_back(boost::dynamic_bitset<>(reagSet.size()));
  }

  // The tests must be applied for all permutations of reagent list against
  // fragment.
  auto reagentOrders =
      details::permMFromN(molFrags.size(), reaction->d_reagents.size());

  for (const auto &ro : reagentOrders) {
    //    std::cout << "Reagent orders : ";
    //    for (auto &i : ro) {
    //      std::cout << i << " ";
    //    }
    //    std::cout << std::endl;
    // Match the fragment to the reagent set in this order.
    std::vector<boost::dynamic_bitset<>> thisPass;
    for (const auto &reagSet : reaction->d_reagents) {
      thisPass.push_back(boost::dynamic_bitset<>(reagSet.size()));
    }
    boost::dynamic_bitset<> fragsMatched(ro.size());
    for (size_t i = 0; i < ro.size(); ++i) {
      const auto &reagSet = reaction->d_reagents[ro[i]];
      const auto &passedScreensSet = passedScreens[ro[i]];
      for (size_t j = 0; j < reagSet.size(); ++j) {
        if (passedScreensSet[j]) {
          auto &reag = reagSet[j];
          if (SubstructMatch(*reag->mol(), *molFrags[i], dontCare)) {
            thisPass[ro[i]][j] = true;
            fragsMatched[i] = true;
          }
        }
      }
    }
    // If all the fragments had a match, these results are valid.
    if (fragsMatched.count() == fragsMatched.size()) {
      for (size_t i = 0; i < reagsToUse.size(); ++i) {
        reagsToUse[i] |= thisPass[i];
      }
    }
  }
  // If all bits in one of the bitsets is unset, it means that nothing matched
  // that reagent.  If at least one of the bitsets has a set bit, all products
  // incorporating the reagent with no bits sets must match the query so
  // should be used.
  bool someSet = false;
  for (auto &rtu : reagsToUse) {
    if (rtu.count()) {
      someSet = true;
      break;
    }
  }
  if (someSet) {
    for (auto &rtu : reagsToUse) {
      if (!rtu.count()) {
        rtu.set();
      }
    }
  }
  return reagsToUse;
}

// Return true if all the fragments have a connector region that matches
// something in the reaction, false otherwise.
bool checkConnectorRegions(const std::vector<std::shared_ptr<ROMol>> &molFrags,
                           std::unique_ptr<ReagentSet> &reaction) {
  const auto &rxnConnRegs = reaction->connectorRegions();
  const auto &rxnConnRegsFP = reaction->connRegFP();
  RDKit::MatchVectType dontCare;
  for (const auto &frag : molFrags) {
    auto connRegs = getConnRegion(*frag);
    if (!connRegs) {
      // There were no connector atoms.
      continue;
    }
    std::vector<std::unique_ptr<ROMol>> splitConnRegs;
    MolOps::getMolFrags(*connRegs, splitConnRegs, false);
    bool connRegFound = false;
    for (const auto &cr : splitConnRegs) {
      std::unique_ptr<ExplicitBitVect> crfp(PatternFingerprintMol(*cr));
      //      std::cout << "Conn region : " << MolToSmiles(*cr) << std::endl;
      if (AllProbeBitsMatch(*crfp, *rxnConnRegsFP)) {
        //        std::cout << MolToSmiles(*cr) << " might be in " <<
        //        reaction->d_id
        //                  << std::endl;
        for (const auto &rxncr : rxnConnRegs) {
          if (SubstructMatch(*rxncr, *cr, dontCare)) {
            connRegFound = true;
            break;
          }
        }
      }
      if (connRegFound) {
        //        std::cout << "   It was" << std::endl;
        break;
      } else {
        //        std::cout << "   It wasn't" << std::endl;
      }
    }
    if (!connRegFound) {
      return false;
    }
  }
  //  std::cout << "CONN Reg OK" << std::endl;
  return true;
}
}  // namespace

std::vector<std::unique_ptr<ROMol>> Hyperspace::searchFragSet(
    const std::vector<std::shared_ptr<ROMol>> &fragSet) {
  std::vector<std::unique_ptr<ROMol>> results;

  auto conns = getConnectorPattern(fragSet);
  for (auto &it : d_reactions) {
    auto &reaction = it.second;
    //    std::cout << "Searching for " << fragSet.size() << " ::: ";
    //    for (const auto &f : fragSet) {
    //      std::cout << f->getNumAtoms() << " : " << MolToSmiles(*f) << " : "
    //                << MolToSmarts(*f) << " :: ";
    //    }
    //    std::cout << " in " << reaction->d_id << " : "
    //              << reaction->d_reagents.size() << std::endl;
    // It can't be a hit if the number of fragments is more than the number
    // of reagent sets because some of the molecule won't be matched in any
    // of the potential products.  It can be less, in which case the unused
    // reagent set will be used completely, possibly resulting in a large
    // number of hits.
    if (fragSet.size() > reaction->d_reagents.size()) {
      continue;
    }

    // Check that all the frags have a connector region that matches something
    // in this reaction set.  Skip if not.
    if (!checkConnectorRegions(fragSet, reaction)) {
      continue;
    }
    auto passedScreens = screenReagentsWithFPs(fragSet, reaction);

    // If none of the reagents passed the screens, move right along, nothing
    // to see.
    bool skip = true;
    for (const auto &ps : passedScreens) {
      if (ps.count()) {
        skip = false;
      }
    }
    if (skip) {
      continue;
    }

    // Get all the possible combinations of connector numbers.  So if the
    // fragmented molecule is C[1*].N[2*] we also try C[2*].N[1*] because
    // that might be how they're labelled in the reaction database.
    auto connCombs =
        getConnectorPermutations(fragSet, conns, reaction->d_connectors);

    // Find all reagents that match the fragments
    for (auto &connComb : connCombs) {
      auto theseReagents = getHitReagents(connComb, passedScreens, reaction);
      if (!theseReagents.empty()) {
        //        for (const auto &tr : theseReagents) {
        //          std::cout << "reagents : " << tr << std::endl;
        //        }
        int totReags =
            std::reduce(theseReagents.begin(), theseReagents.end(), 0,
                        [&](int prevRes, const boost::dynamic_bitset<> &s2) {
                          return prevRes + s2.count();
                        });
        //        std::cout << "totReags : " << totReags << std::endl;
        if (totReags) {
          buildHits(theseReagents, it.first, results);
        }
      }
    }
  }
  return results;
}

namespace {
struct Stepper {
  Stepper(std::vector<int> &sizes) : d_sizes(sizes) {
    d_currState = std::vector<int>(sizes.size(), 0);
  }
  void step() {
    // Don't do anything if we're at the end, but expect an infinite
    // loop if the user isn't wise to this.
    if (d_currState[0] == d_sizes[0]) {
      return;
    }
    int i = d_currState.size() - 1;
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
  std::vector<int> d_currState;
  std::vector<int> d_sizes;
};
}  // namespace
void Hyperspace::buildHits(
    const std::vector<boost::dynamic_bitset<>> &reagentsToUse,
    const std::string &reaction_id,
    std::vector<std::unique_ptr<ROMol>> &results) {
  if (reagentsToUse.empty()) {
    return;
  }
  //  std::cout << "Build hits with " << reaction_id << std::endl;
  auto reags = getReagentsToUse(reagentsToUse, reaction_id);
  if (reags.empty()) {
    //    std::cout << "Nothing to do" << std::endl;
    return;
  }

  std::vector<int> numReags;
  for (auto &r : reags) {
    numReags.push_back(r.size());
  }
  Stepper stepper(numReags);
  const int numReactions = reagentsToUse.size();
  MolzipParams params;
  params.label = MolzipLabel::Isotope;
  while (stepper.d_currState[0] != numReags[0]) {
    //    std::cout << "Step : ";
    //    for (auto &i : stepper.d_currState) {
    //      std::cout << i << " ";
    //    }
    //    std::cout << std::endl;
    std::unique_ptr<ROMol> combMol(
        new ROMol(*reags[0][stepper.d_currState[0]]));
    //    for (int i = 0; i < numReactions; ++i) {
    //      std::cout << MolToSmiles(*reags[i][stepper.d_currState[i]]) << "
    //      ";
    //    }
    //    std::cout << std::endl;
    for (int i = 1; i < numReactions; ++i) {
      combMol.reset(combineMols(*combMol, *reags[i][stepper.d_currState[i]]));
    }
    //    std::cout << "Comb mol : " << MolToSmiles(*combMol) << std::endl;
    auto prod = molzip(*combMol, params);
    //    std::cout << "Zipped : " << MolToSmiles(*prod) << std::endl;
    MolOps::sanitizeMol(*static_cast<RWMol *>(prod.get()));
    results.push_back(std::move(prod));
    stepper.step();
  }
}

std::vector<std::vector<ROMol *>> Hyperspace::getReagentsToUse(
    const std::vector<boost::dynamic_bitset<>> &reagentsToUse,
    const std::string &reaction_id) const {
  if (const auto &it = d_reactions.find(reaction_id); it == d_reactions.end()) {
    throw std::runtime_error("Reaction " + reaction_id +
                             "not in the reaction set.");
  }
  const auto &reaction = d_reactions.find(reaction_id)->second;

  std::vector<std::vector<ROMol *>> reags(reaction->d_reagents.size(),
                                          std::vector<ROMol *>());
  for (size_t i = 0; i < reagentsToUse.size(); ++i) {
    for (size_t j = 0; j < reagentsToUse[i].size(); ++j) {
      if (reagentsToUse[i][j]) {
        reags[i].push_back(reaction->d_reagents[i][j]->mol().get());
      }
    }
  }
  return reags;
}
}  // namespace RDKit::HyperspaceSSSearch
