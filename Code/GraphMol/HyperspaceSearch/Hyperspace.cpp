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
#include <GraphMol/QueryOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
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
  for (auto &fragSet : fragments) {
    for (auto &fraggedMol : fragSet) {
      auto theseResults = searchFragSet(*fraggedMol);
      // Just for safety, do a final check
      for (auto &hitMol : theseResults) {
        if (SubstructMatch(*hitMol, query, dontCare)) {
          results.emplace_back(std::unique_ptr<ROMol>(hitMol.release()));
        }
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
          nextReag[3], std::make_unique<ReactionSet>(nextReag[3])));
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
    std::cout << "Connectors : " << reaction->d_id << " : "
              << reaction->d_connectors << std::endl;
  }
}

namespace {
// Return a bitset giving the different connector types in this
// molecule.
boost::dynamic_bitset<> getConnectorPattern(const ROMol &mol) {
  boost::dynamic_bitset<> conns(4);
  for (const auto &a : mol.atoms()) {
    if (!a->getAtomicNum()) {
      conns.set(a->getIsotope() - 1);
    }
  }
  return conns;
}

// Return copies of the mol fragments will all permutations of the connectors
// in the reaction onto the connectors in the fragments.
// E.g. if the reaction has 3 connectors, 1, 2 and 3 and the fragged mol has
// 2, return all permutations of 2 from 3.  It's ok if the fragged mol doesn't
// have all the connections in the reaction, although this may well result in a
// lot of hits.
std::vector<std::vector<std::unique_ptr<RWMol>>> getConnectorPermutations(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    const boost::dynamic_bitset<> &fragConns,
    const boost::dynamic_bitset<> &reactionConns) {
  std::vector<std::vector<std::unique_ptr<RWMol>>> connPerms;
  for (const auto &mf : molFrags) {
    std::cout << MolToSmiles(*mf) << " ";
  }
  std::cout << " :: " << fragConns << " : " << reactionConns << std::endl;
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
  std::cout << "Num frag conns : " << numFragConns
            << " num reaction conns : " << reactionConns.count() << std::endl;
  auto perms = details::permMFromN(numFragConns, reactionConns.count());
  std::cout << "PERMS" << std::endl;

  for (const auto &perm : perms) {
    std::cout << "Next perm : ";
    for (auto &i : perm) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
    connPerms.push_back(std::vector<std::unique_ptr<RWMol>>());
    // Copy the fragments and set the isotope numbers according to this
    // permutation.
    for (const auto &f : molFrags) {
      connPerms.back().emplace_back(new RWMol(*f));
      for (auto atom : connPerms.back().back()->atoms()) {
        if (!atom->getAtomicNum()) {
          for (size_t i = 0; i < perm.size(); ++i) {
            if (atom->getIsotope() == i + 1) {
              atom->setIsotope(perm[i] + 1);
            }
          }
        }
      }
    }
    for (const auto &f : connPerms.back()) {
      std::cout << MolToSmiles(*f) << "::";
    }
    std::cout << std::endl;
  }

  return connPerms;
}

std::vector<boost::dynamic_bitset<>> getHitReagents(
    std::vector<std::unique_ptr<RWMol>> &fraggedMol,
    std::unique_ptr<ReactionSet> &reaction) {
  RDKit::MatchVectType dontCare;
  std::vector<boost::dynamic_bitset<>> reagsToUse;
  for (const auto &reagSet : reaction->d_reagents) {
    reagsToUse.push_back(boost::dynamic_bitset<>(reagSet.size()));
  }
  auto reagentOrders =
      details::permMFromN(fraggedMol.size(), reaction->d_reagents.size());
  boost::dynamic_bitset<> fragsMatched(fraggedMol.size());
  for (const auto &ro : reagentOrders) {
    //    std::cout << "Reagent orders : ";
    //    for (auto &i : ro) {
    //      std::cout << i << " ";
    //    }
    //    std::cout << std::endl;
    // Match the fragment to the reagent set in this order.
    for (size_t i = 0; i < ro.size(); ++i) {
      const auto &reagSet = reaction->d_reagents[ro[i]];
      for (size_t j = 0; j < reagSet.size(); ++j) {
        auto &reag = reagSet[j];
        if (!reag->d_mol) {
          //          std::cout << "making reagent " << reag->d_id << " from "
          //                    << reag->d_smiles << std::endl;
          reag->d_mol.reset(
              v2::SmilesParse::MolFromSmiles(reag->d_smiles).release());
        }
        //        std::cout << "seeing if reagent " << reag->d_smiles << " of "
        //        << ro[i]
        //                  << " has " << MolToSmiles(*fraggedMol[i]) <<
        //                  std::endl;
        if (SubstructMatch(*reag->d_mol, *fraggedMol[i], dontCare)) {
          //          std::cout << "match " << reag->d_smiles << " has "
          //                    << MolToSmiles(*fraggedMol[i]) << std::endl;
          reagsToUse[ro[i]][j] = true;
          fragsMatched[i] = true;
        }
      }
    }
  }
  std::cout << "Frags matched : " << fragsMatched << " : "
            << fragsMatched.count() << std::endl;
  // If all the fragments aren't matched, it's not a good solution, so
  // clear reagsToUse.  If all bits in one of the bitsets is unset, it
  // means that nothing matched that reagent.  Therefore all products
  // incorporating that reagent match the query so should be used.
  if (fragsMatched.count() != fragsMatched.size()) {
    reagsToUse.clear();
  } else {
    for (auto &rtu : reagsToUse) {
      if (!rtu.count()) {
        rtu.set();
      }
    }
  }
  return reagsToUse;
}
}  // namespace

std::vector<std::unique_ptr<ROMol>> Hyperspace::searchFragSet(
    const ROMol &fraggedMol) {
  PRECONDITION(fraggedMol.getNumAtoms() != 0,
               "Search query must contain atoms.");

  std::vector<std::unique_ptr<ROMol>> results;

  std::vector<std::unique_ptr<ROMol>> molFrags;
  auto numFrags = MolOps::getMolFrags(fraggedMol, molFrags, false);
  std::cout << "Num frags : " << numFrags << std::endl;
  auto conns = getConnectorPattern(fraggedMol);
  for (auto &it : d_reactions) {
    auto &reaction = it.second;
    std::cout << "Searching for " << MolToSmiles(fraggedMol) << " in "
              << reaction->d_id << " : " << reaction->d_reagents.size()
              << std::endl;
    // It can't be a hit if the number of fragments is more than the number
    // of reagent sets because some of the molecule won't be matched in any
    // of the potential products.  It can be less, in which case the unused
    // reagent set will be used completely, possibly resulting in a large
    // number of hits.
    if (numFrags > reaction->d_reagents.size()) {
      continue;
    }
    // Get all the possible combinations of connector numbers.  So if the
    // fragmented molecule is C[1*].N[2*] we also try C[2*].N[1*] because
    // that might be how they're labelled in the reaction database.
    auto connCombs =
        getConnectorPermutations(molFrags, conns, reaction->d_connectors);

    // Find all reagents that match the fragments
    for (auto &connComb : connCombs) {
      auto theseReagents = getHitReagents(connComb, reaction);
      if (!theseReagents.empty()) {
        buildHits(theseReagents, it.first, results);
        for (const auto &tr : theseReagents) {
          std::cout << "reagents : " << tr << std::endl;
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
  std::cout << "Build hits with " << reaction_id << std::endl;
  auto reags = getReagentsToUse(reagentsToUse, reaction_id);
  if (reags.empty()) {
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
    std::cout << "Step : ";
    for (auto &i : stepper.d_currState) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
    std::unique_ptr<ROMol> combMol(
        new ROMol(*reags[0][stepper.d_currState[0]]));
    for (size_t i = 0; i < numReactions; ++i) {
      std::cout << MolToSmiles(*reags[i][stepper.d_currState[i]]) << " ";
    }
    std::cout << std::endl;
    for (size_t i = 1; i < numReactions; ++i) {
      combMol.reset(combineMols(*combMol, *reags[i][stepper.d_currState[i]]));
    }
    std::cout << "Comb mol : " << MolToSmiles(*combMol) << std::endl;
    auto prod = molzip(*combMol, params);
    std::cout << "Zipped : " << MolToSmiles(*prod) << std::endl;
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
        reags[i].push_back(reaction->d_reagents[i][j]->d_mol.get());
      }
    }
  }
  return reags;
}
}  // namespace RDKit::HyperspaceSSSearch
