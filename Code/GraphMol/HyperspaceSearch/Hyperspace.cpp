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

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "Hyperspace.h"
#include "HyperspaceSubstructureSearch.h"

namespace RDKit {
namespace HyperspaceSSSearch {
Hyperspace::Hyperspace(const std::string &fileName) : d_fileName(fileName) {
  readFile();
  assignConnectorsUsed();
}

std::vector<std::unique_ptr<ROMol>> Hyperspace::search(
    const RDKit::ROMol &query, unsigned int maxBondSplits) {
  auto fragments = details::splitMolecule(query, maxBondSplits);
  for (auto &fragSet : fragments) {
    for (auto &fraggedMol : fragSet) {
      auto theseResults = searchFragSet(fraggedMol);
    }
  }

  std::vector<std::unique_ptr<ROMol>> results;
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
  std::cout << "parsing : " << firstLine << std::endl;
  std::regex regexz("\\s+");

  auto lineParts = splitLine(firstLine, regexz);
  if (lineParts != std::vector<std::string>{"SMILES", "synton_id", "synton#",
                                            "reaction_id"}) {
    throw std::runtime_error("Bad format for hyperspace file " + d_fileName);
  }
  std::string nextLine;
  int lineNum = 1;
  auto currReaction = std::make_unique<ReactionSet>();

  while (getline(ifs, nextLine)) {
    ++lineNum;
    auto nextReag = splitLine(nextLine, regexz);
    if (nextReag.size() != 4) {
      throw std::runtime_error("Bad format for hyperspace file " + d_fileName +
                               " on line " + std::to_string(lineNum));
    }
    auto reagent = std::make_unique<Reagent>(nextReag[0], nextReag[1]);
    if (currReaction->d_id != nextReag[3]) {
      if (!currReaction->d_reagents.empty()) {
        d_reactions.push_back(std::move(currReaction));
        currReaction = std::make_unique<ReactionSet>();
      }
      currReaction->d_id = nextReag[3];
    }
    size_t synthonId = std::stoi(nextReag[2]);
    if (synthonId >= currReaction->d_reagents.size()) {
      for (size_t i = currReaction->d_reagents.size(); i < synthonId + 1; ++i) {
        currReaction->d_reagents.push_back(
            std::vector<std::unique_ptr<Reagent>>());
      }
    }
    currReaction->d_reagents[synthonId].push_back(std::move(reagent));
  }
  if (!currReaction->d_reagents.empty()) {
    d_reactions.push_back(std::move(currReaction));
  }
}

void Hyperspace::assignConnectorsUsed() {
  const static std::vector<std::regex> connRegexs{
      std::regex(R"([U])"), std::regex(R"([Np])"), std::regex(R"([Pu])"),
      std::regex(R"([Am])")};
  for (auto &reaction : d_reactions) {
    reaction->d_connectors.resize(4, 0);
    for (auto &reagSet : reaction->d_reagents) {
      for (auto &reag : reagSet) {
        for (size_t i = 0; i < 4; ++i) {
          if (std::regex_search(reag->d_smiles, connRegexs[i])) {
            reaction->d_connectors.set(i);
          }
        }
      }
    }
    std::cout << reaction->d_id << " : ";
    for (int i = 0; i < 4; ++i) {
      std::cout << reaction->d_connectors[i] << " ";
    }
    std::cout << std::endl;
  }
}

std::vector<std::unique_ptr<ROMol>> Hyperspace::searchFragSet(
    std::unique_ptr<ROMol> &fraggedMol) {
  std::vector<std::unique_ptr<ROMol>> results;
  for (auto &reaction : d_reactions) {
    std::cout << "Searching for " << MolToSmiles(*fraggedMol) << " in "
              << reaction->d_id << std::endl;
    std::vector<boost::dynamic_bitset<>> reagentsToUse;
    for (const auto &reags : reaction->d_reagents) {
      reagentsToUse.push_back(boost::dynamic_bitset<>(reags.size()));
      reagentsToUse.back().set();
    }
  }
  // Next, split the fraggedMol into components, make copies of each
  // one with appropriate connectors based on the d_connectors and do
  // the search.  About the top of column 2, page 'D'.

  return results;
}

}  // namespace HyperspaceSSSearch
}  // namespace RDKit
