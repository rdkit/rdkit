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
#include <iostream>
#include <random>
#include <regex>
#include <set>
#include <string>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSubstructureSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace RDKit::SynthonSpaceSearch {

std::int64_t SynthonSpace::getNumProducts() const {
  std::int64_t totSize = 0;
  for (const auto &reaction : d_reactions) {
    const auto &rxn = reaction.second;
    size_t thisSize = 1;
    for (const auto &r : rxn->getSynthons()) {
      thisSize *= r.size();
    }
    totSize += thisSize;
  }
  return totSize;
}

SubstructureResults SynthonSpace::substructureSearch(
    const ROMol &query, SynthonSpaceSearchParams params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceSubstructureSearcher ssss(query, params, *this);
  return ssss.search();
}

namespace {
inline std::vector<std::string> splitLine(const std::string &str,
                                          const std::regex &regexz) {
  return {std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
          std::sregex_token_iterator()};
}

// The Enamine/Chemspace files come with the connection points on the
// synthons marked with [U], [Np], [Pu], [Am].  These need to be converted
// to dummy atoms with isotope labels (1, 2, 3, 4 respectively) which can
// be done safely on the SMILES string.
void fixConnectors(std::string &smiles) {
  for (int i = 0; i < MAX_CONNECTOR_NUM; ++i) {
    std::string regex =
        std::regex_replace(CONNECTOR_SYMBOLS[i], std::regex(R"(\[)"), R"(\[)");
    regex = std::regex_replace(regex, std::regex(R"(\])"), R"(\])");
    std::string repl = "[" + std::to_string(i + 1) + "*]";
    smiles = std::regex_replace(smiles, std::regex(regex), repl);
  }
}

}  // namespace

void SynthonSpace::readTextFile(const std::string &inFilename) {
  d_fileName = inFilename;
  std::ifstream ifs(d_fileName);
  if (!ifs.is_open() || ifs.bad()) {
    throw std::runtime_error("Couldn't open file " + d_fileName);
  }
  std::regex regexz("[\\s,]+");

  std::vector<std::vector<std::string>> firstLineOpts{
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id"},
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id",
                               "release"},
      std::vector<std::string>{"SMILES", "synton_id", "synton_role",
                               "reaction_id"}};
  int format = -1;
  std::string nextLine;
  int lineNum = 1;

  while (getline(ifs, nextLine)) {
    ++lineNum;
    if (nextLine.empty() || nextLine[0] == '#') {
      continue;
    }
    if (format == -1) {
      auto lineParts = splitLine(nextLine, regexz);
      for (size_t i = 0; i < firstLineOpts.size(); ++i) {
        if (lineParts == firstLineOpts[i]) {
          format = static_cast<int>(i);
        }
      }
      if (format == -1) {
        throw std::runtime_error("Bad format for SynthonSpace file " +
                                 d_fileName);
      }
      continue;
    }
    auto nextSynthon = splitLine(nextLine, regexz);
    if (nextSynthon.size() < 4) {
      throw std::runtime_error("Bad format for SynthonSpace file " +
                               d_fileName + " on line " +
                               std::to_string(lineNum));
    }
    if (auto it = d_reactions.find(nextSynthon[3]); it == d_reactions.end()) {
      d_reactions.insert(std::make_pair(
          nextSynthon[3], std::make_unique<SynthonSet>(nextSynthon[3])));
    }
    fixConnectors(nextSynthon[0]);
    auto &currReaction = d_reactions[nextSynthon[3]];
    int synthonNum{std::numeric_limits<int>::max()};
    if (format == 0 || format == 1) {
      synthonNum = std::stoi(nextSynthon[2]);
    } else if (format == 2) {
      // in this case it's a string "synton_2" etc.
      synthonNum = std::stoi(nextSynthon[2].substr(7));
    }
    currReaction->addSynthon(synthonNum, std::make_unique<Synthon>(Synthon(
                                             nextSynthon[0], nextSynthon[1])));
  }

  // Do some final processing.
  for (auto &[fst, snd] : d_reactions) {
    snd->buildConnectorRegions();
    snd->assignConnectorsUsed();
  }
}

void SynthonSpace::writeDBFile(const std::string &outFilename) const {
  std::ofstream os(outFilename, std::fstream::binary | std::fstream::trunc);
  streamWrite(os, d_reactions.size());
  for (const auto &[fst, snd] : d_reactions) {
    snd->writeToDBStream(os);
  }
  os.close();
}

void SynthonSpace::readDBFile(const std::string &inFilename) {
  d_fileName = inFilename;
  try {
    std::ifstream is(d_fileName, std::fstream::binary);
    size_t numRS;
    streamRead(is, numRS);
    for (size_t i = 0; i < numRS; ++i) {
      std::unique_ptr<SynthonSet> rs = std::make_unique<SynthonSet>();
      rs->readFromDBStream(is);
      d_reactions.insert(std::make_pair(rs->getId(), std::move(rs)));
    }
  } catch (std::exception &e) {
    std::cerr << "Error : " << e.what() << " for file " << d_fileName << "\n";
    exit(1);
  }
}

void SynthonSpace::summarise(std::ostream &os) const {
  os << "Read from file " << d_fileName << "\n"
     << "Number of reactions : " << d_reactions.size() << "\n";
  size_t totSize = 0;
  for (const auto &reaction : d_reactions) {
    const auto &rxn = reaction.second;
    os << "Reaction name " << rxn->getId() << "\n";
    size_t thisSize = 1;
    for (size_t i = 0; i < rxn->getSynthons().size(); ++i) {
      os << "  Synthon set " << i << " has " << rxn->getSynthons()[i].size()
         << " synthons"
         << "\n";
      thisSize *= rxn->getSynthons()[i].size();
    }
    totSize += thisSize;
  }
  os << "Approximate number of molecules in SynthonSpace : " << totSize
     << std::endl;
}

void SynthonSpace::buildSynthonFingerprints() {
  if (!dp_fpGenerator) {
    dp_fpGenerator.reset(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
  }
  for (const auto &it : d_reactions) {
    it.second->buildSynthonFingerprints(dp_fpGenerator);
  }
}

}  // namespace RDKit::SynthonSpaceSearch
