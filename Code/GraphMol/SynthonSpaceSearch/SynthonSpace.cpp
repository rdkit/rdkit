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
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceFingerprintSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSubstructureSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace RDKit::SynthonSpaceSearch {

// used for serialization
constexpr int32_t versionMajor = 1;
constexpr int32_t versionMinor = 0;
constexpr int32_t endianId = 0xa100f;

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

const std::unique_ptr<FingerprintGenerator<std::uint64_t>> &
SynthonSpace::getFPGenerator() {
  if (!dp_fpGenerator) {
    dp_fpGenerator.reset(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
  }
  return dp_fpGenerator;
}

SubstructureResults SynthonSpace::substructureSearch(
    const ROMol &query, SynthonSpaceSearchParams params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceSubstructureSearcher ssss(query, params, *this);
  return ssss.search();
}

SubstructureResults SynthonSpace::fingerprintSearch(
    const ROMol &query, SynthonSpaceSearchParams params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceFingerprintSearcher ssss(query, params, *this);
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
  for (unsigned int i = 0; i < MAX_CONNECTOR_NUM; ++i) {
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

  streamWrite(os, endianId);
  streamWrite(os, versionMajor);
  streamWrite(os, versionMinor);

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

    int32_t endianTest;
    streamRead(is, endianTest);
    if (endianTest != endianId) {
      throw std::runtime_error("Endianness mismatch in SynthonSpace file " +
                               d_fileName);
    }
    int32_t majorVersion;
    streamRead(is, majorVersion);
    int32_t minorVersion;
    streamRead(is, minorVersion);
    if (majorVersion > versionMajor ||
        (majorVersion == versionMajor && minorVersion > versionMinor)) {
      BOOST_LOG(rdWarningLog)
          << "Deserializing from a version number (" << majorVersion << "."
          << minorVersion << ")"
          << "that is higher than our version (" << versionMajor << "."
          << versionMinor << ").\nThis probably won't work." << std::endl;
    }
    // version sanity checking
    if (majorVersion > 1000 || minorVersion > 100) {
      throw std::runtime_error("unreasonable version numbers");
    }
    majorVersion = 1000 * majorVersion + minorVersion * 10;
    size_t numRS;
    streamRead(is, numRS);
    for (size_t i = 0; i < numRS; ++i) {
      std::unique_ptr<SynthonSet> rs = std::make_unique<SynthonSet>();
      rs->readFromDBStream(is, majorVersion);
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

bool SynthonSpace::hasFingerprints() const {
  if (d_reactions.empty()) {
    return false;
  }
  return d_reactions.begin()->second->hasFingerprints();
}

void SynthonSpace::buildSynthonFingerprints(const std::string &fpType) {
  if (fpType != d_fpType) {
    dp_fpGenerator.reset();
  }
  if (!dp_fpGenerator) {
    if (fpType == "Morgan_2") {
      dp_fpGenerator.reset(
          MorganFingerprint::getMorganGenerator<std::uint64_t>(2u));
    } else if (fpType == "Morgan_3") {
      dp_fpGenerator.reset(
          MorganFingerprint::getMorganGenerator<std::uint64_t>(3u));
    } else if (fpType == "RDKit_7") {
      dp_fpGenerator.reset(RDKitFP::getRDKitFPGenerator<std::uint64_t>(1u, 7u));
    } else {
      throw ValueErrorException("Unknown fingerprint type '" + fpType + "'.");
    }
  }
  for (const auto &it : d_reactions) {
    it.second->buildSynthonFingerprints(dp_fpGenerator);
  }
}

}  // namespace RDKit::SynthonSpaceSearch