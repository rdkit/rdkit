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
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceFingerprintSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSubstructureSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/StreamOps.h>

namespace RDKit::SynthonSpaceSearch {

// used for serialization
constexpr int32_t versionMajor = 2;
constexpr int32_t versionMinor = 1;
constexpr int32_t endianId = 0xa100f;

SynthonSpace::~SynthonSpace() {
  if (d_dbis) {
    d_dbis->close();
  }
}

size_t SynthonSpace::getNumReactions() const {
  if (d_lowMem) {
    return d_reactionPos.size();
  }
  return d_reactions.size();
}

std::vector<std::string> SynthonSpace::getReactionNames() const {
  std::vector<std::string> reactionNames;
  if (d_lowMem) {
    reactionNames.reserve(d_reactionPos.size());
    for (const auto &reaction : d_reactionPos) {
      reactionNames.push_back(reaction.first);
    }
  } else {
    reactionNames.reserve(d_reactions.size());
    for (const auto &reaction : d_reactions) {
      reactionNames.push_back(reaction.first);
    }
  }
  return reactionNames;
}

const std::shared_ptr<SynthonSet> SynthonSpace::getReaction(
    std::string reactionName) {
  // Even in lowMem mode, the reaction might be cached here.
  if (const auto &it = d_reactions.find(reactionName);
      it != d_reactions.end()) {
    return it->second;
  }
  if (d_lowMem) {
    if (!d_dbis && !d_fileName.empty()) {
      openAndCheckDBFile();
    }
    d_reactions.clear();
    if (const auto &it = d_reactionPos.find(reactionName);
        it != d_reactionPos.end()) {
      d_dbis->seekg(it->second);
      auto rs = std::make_shared<SynthonSet>();
      rs->readFromDBStream(*d_dbis, d_fileMajorVersion);
      d_reactions.insert(std::make_pair(rs->getId(), rs));
      return rs;
    }
  }
  throw std::runtime_error("Could not find synthon set for reaction " +
                           reactionName);
}

std::int64_t SynthonSpace::getNumProducts() const {
  std::int64_t totSize = 0;
  for (const auto &[id, rxn] : d_reactions) {
    size_t thisSize = 1;
    for (const auto &r : rxn->getSynthons()) {
      thisSize *= r.size();
    }
    totSize += thisSize;
  }
  return totSize;
}

SearchResults SynthonSpace::substructureSearch(
    const ROMol &query, const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceSubstructureSearcher ssss(query, params, *this);
  return ssss.search();
}

SearchResults SynthonSpace::fingerprintSearch(
    const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
    const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceFingerprintSearcher ssss(query, fpGen, params, *this);
  return ssss.search();
}

namespace {
std::vector<std::string> splitLine(const std::string &str,
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

int deduceFormat(const std::string &line) {
  // formats are based on the headers which should be one of the forms here.
  // If the columns are white-space separated the return value is 0-2, if
  // comma separated 3-5.
  static const std::vector<std::vector<std::string>> firstLineOpts{
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id"},
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id",
                               "release"},
      std::vector<std::string>{"SMILES", "synton_id", "synton_role",
                               "reaction_id"}};
  static const std::regex regexws("\\s+");
  auto lineParts = splitLine(line, regexws);
  for (size_t i = 0; i < firstLineOpts.size(); ++i) {
    if (lineParts == firstLineOpts[i]) {
      return static_cast<int>(i);
    }
  }
  static const std::regex regexc(",+");
  lineParts = splitLine(line, regexc);
  for (size_t i = 0; i < firstLineOpts.size(); ++i) {
    if (lineParts == firstLineOpts[i]) {
      return static_cast<int>(i + firstLineOpts.size());
    }
  }
  return -1;
}

std::vector<std::string> readSynthonLine(std::istream &is, int &lineNum,
                                         int &format,
                                         const std::string &fileName) {
  static const std::regex regexws("\\s+");
  static const std::regex regexc(",+");

  std::vector<std::string> nextSynthon;
  auto nextLine = getLine(is);
  ++lineNum;
  if (nextLine.empty() || nextLine[0] == '#') {
    return nextSynthon;
  }
  if (format == -1) {
    format = deduceFormat(nextLine);
    if (format == -1) {
      throw std::runtime_error("Bad format for SynthonSpace file " + fileName);
    }
    return nextSynthon;
  }
  if (format < 3) {
    nextSynthon = splitLine(nextLine, regexws);
  } else {
    nextSynthon = splitLine(nextLine, regexc);
  }
  if (nextSynthon.size() < 4) {
    throw std::runtime_error("Bad format for SynthonSpace file " + fileName +
                             " on line " + std::to_string(lineNum));
  }
  return nextSynthon;
}
int getSynthonNum(int format, const std::string &synthon) {
  int synthonNum{std::numeric_limits<int>::max()};
  if (format == 0 || format == 1 || format == 3 || format == 4) {
    synthonNum = std::stoi(synthon);
  } else if (format == 2 || format == 5) {
    // in this case it's a string "synton_2" etc.
    synthonNum = std::stoi(synthon.substr(7));
  }
  return synthonNum;
}
}  // namespace

void SynthonSpace::readTextFile(const std::string &inFilename) {
  d_fileName = inFilename;
  d_lowMem = false;
  std::ifstream ifs(d_fileName);
  if (!ifs.is_open() || ifs.bad()) {
    throw std::runtime_error("Couldn't open file " + d_fileName);
  }

  int format = -1;
  std::string nextLine;
  int lineNum = 1;
  while (!ifs.eof()) {
    auto nextSynthon = readSynthonLine(ifs, lineNum, format, d_fileName);
    if (nextSynthon.empty()) {
      continue;
    }
    if (auto it = d_reactions.find(nextSynthon[3]); it == d_reactions.end()) {
      d_reactions.insert(std::make_pair(
          nextSynthon[3], std::make_unique<SynthonSet>(nextSynthon[3])));
    }
    fixConnectors(nextSynthon[0]);
    auto &currReaction = d_reactions[nextSynthon[3]];
    auto synthonNum = getSynthonNum(format, nextSynthon[2]);
    currReaction->addSynthon(synthonNum, std::make_unique<Synthon>(Synthon(
                                             nextSynthon[0], nextSynthon[1])));
  }

  // Do some final processing.
  for (auto &[id, reaction] : d_reactions) {
    reaction->removeEmptySynthonSets();
    reaction->transferProductBondsToSynthons();
    reaction->buildConnectorRegions();
    reaction->assignConnectorsUsed();
  }
}

void SynthonSpace::writeDBFile(const std::string &outFilename) const {
  std::ofstream os(outFilename, std::fstream::binary | std::fstream::trunc);

  streamWrite(os, endianId);
  streamWrite(os, versionMajor);
  streamWrite(os, versionMinor);
  streamWrite(os, hasFingerprints());
  if (hasFingerprints()) {
    streamWrite(os, d_fpType);
  }
  streamWrite(os, static_cast<std::uint64_t>(d_reactions.size()));
  for (const auto &[reactionId, reaction] : d_reactions) {
    reaction->writeToDBStream(os);
  }
  os.close();
}

void SynthonSpace::readDBFile(const std::string &inFilename, bool lowMem) {
  d_fileName = inFilename;
  d_lowMem = lowMem;
  openAndCheckDBFile();
  bool hasFPs;
  streamRead(*d_dbis, hasFPs);
  if (hasFPs) {
    streamRead(*d_dbis, d_fpType, 0);
  }
  std::uint64_t numRS;
  streamRead(*d_dbis, numRS);
  for (size_t i = 0; i < numRS; ++i) {
    auto reactionPos = d_dbis->tellg();
    auto rs = std::make_shared<SynthonSet>();
    rs->readFromDBStream(*d_dbis, d_fileMajorVersion);
    if (!d_lowMem) {
      d_reactions.insert(std::make_pair(rs->getId(), rs));
    }
    d_reactionPos.insert(std::make_pair(rs->getId(), reactionPos));
  }
}

void SynthonSpace::summarise(std::ostream &os) const {
  os << "Read from file " << d_fileName << "\n"
     << "Number of reactions : " << d_reactions.size() << "\n";
  size_t totSize = 0;
  for (const auto &[id, rxn] : d_reactions) {
    os << "Reaction name " << id << "\n";
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

void SynthonSpace::writeEnumeratedFile(const std::string &outFilename) const {
  BOOST_LOG(rdWarningLog) << "Writing the enumerated file may take quite"
                          << " some time and result in a large file."
                          << std::endl;
  std::ofstream os(outFilename);
  for (const auto &[fst, snd] : d_reactions) {
    snd->enumerateToStream(os);
  }
  os.close();
}

bool SynthonSpace::hasFingerprints() const { return !d_fpType.empty(); }

void SynthonSpace::buildSynthonFingerprints(
    const FingerprintGenerator<std::uint64_t> &fpGen) {
  if (const auto fpType = fpGen.infoString();
      fpType != d_fpType || !hasFingerprints()) {
    BOOST_LOG(rdWarningLog)
        << "Building the fingerprints may take some time." << std::endl;
    d_fpType = fpType;
    for (const auto &[id, synthSet] : d_reactions) {
      synthSet->buildSynthonFingerprints(fpGen);
    }
  }
}

bool SynthonSpace::hasAddAndSubstractFingerprints() const {
  if (d_reactions.empty()) {
    return false;
  }
  return d_reactions.begin()->second->hasAddAndSubtractFPs();
}

void SynthonSpace::buildAddAndSubstractFingerprints(
    const FingerprintGenerator<std::uint64_t> &fpGen) {
  for (const auto &[id, synthSet] : d_reactions) {
    synthSet->buildAddAndSubtractFPs(fpGen);
  }
}

namespace {
void finalizeReaction(SynthonSet *reaction,
                      const FingerprintGenerator<std::uint64_t> *fpGen) {
  if (reaction) {
    reaction->removeEmptySynthonSets();
    reaction->transferProductBondsToSynthons();
    reaction->buildConnectorRegions();
    reaction->assignConnectorsUsed();
    if (fpGen) {
      reaction->buildSynthonFingerprints(*fpGen);
      reaction->buildAddAndSubtractFPs(*fpGen);
    }
  }
}
}  // namespace

void convertTextToDBFile(const std::string &inFilename,
                         const std::string &outFilename,
                         const FingerprintGenerator<std::uint64_t> *fpGen) {
  std::ifstream ifs(inFilename);
  if (!ifs.is_open() || ifs.bad()) {
    throw std::runtime_error("Couldn't open file " + inFilename);
  }
  std::ofstream os(outFilename, std::fstream::binary | std::fstream::trunc);
  streamWrite(os, endianId);
  streamWrite(os, versionMajor);
  streamWrite(os, versionMinor);
  streamWrite(os, fpGen != nullptr);
  if (fpGen != nullptr) {
    streamWrite(os, fpGen->infoString());
  }
  std::uint64_t numReacts = 0;
  auto numReactPos = os.tellp();
  streamWrite(os, numReacts);

  int format = -1;
  int lineNum = 1;
  std::unique_ptr<SynthonSet> currReaction;
  std::string oldReactionId;
  while (!ifs.eof()) {
    auto nextSynthon = readSynthonLine(ifs, lineNum, format, inFilename);
    if (nextSynthon.empty()) {
      continue;
    }
    fixConnectors(nextSynthon[0]);
    auto synthonNum = getSynthonNum(format, nextSynthon[2]);
    if (!currReaction) {
      ++numReacts;
      currReaction.reset(new SynthonSet(nextSynthon[3]));
      oldReactionId = nextSynthon[3];
    }
    if (oldReactionId != nextSynthon[3]) {
      finalizeReaction(currReaction.get(), fpGen);
      currReaction->writeToDBStream(os);
      ++numReacts;
      currReaction.reset(new SynthonSet(nextSynthon[3]));
      oldReactionId = nextSynthon[3];
    }
    currReaction->addSynthon(synthonNum, std::make_unique<Synthon>(Synthon(
                                             nextSynthon[0], nextSynthon[1])));
  }
  finalizeReaction(currReaction.get(), fpGen);
  currReaction->writeToDBStream(os);
  os.seekp(numReactPos, std::ios::beg);
  streamWrite(os, numReacts);
  os.close();
}

void SynthonSpace::openAndCheckDBFile() {
  if (d_dbis) {
    d_dbis->close();
  }
  d_dbis.reset(new std::ifstream);
  try {
    d_dbis->open(d_fileName, std::fstream::binary);

    int32_t endianTest;
    streamRead(*d_dbis, endianTest);
    if (endianTest != endianId) {
      throw std::runtime_error("Endianness mismatch in SynthonSpace file " +
                               d_fileName);
    }
    streamRead(*d_dbis, d_fileMajorVersion);
    int32_t minorVersion;
    streamRead(*d_dbis, minorVersion);
    if (d_fileMajorVersion > versionMajor ||
        (d_fileMajorVersion == versionMajor && minorVersion > versionMinor)) {
      BOOST_LOG(rdWarningLog)
          << "Deserializing from a version number (" << d_fileMajorVersion
          << "." << minorVersion << ")"
          << "that is higher than our version (" << versionMajor << "."
          << versionMinor << ").\nThis probably won't work." << std::endl;
    }
    // version sanity checking
    if (d_fileMajorVersion > 1000 || minorVersion > 100) {
      throw std::runtime_error("unreasonable version numbers");
    }
    d_fileMajorVersion = 1000 * d_fileMajorVersion + minorVersion * 10;
  } catch (std::exception &e) {
    std::cerr << "Error : " << e.what() << " for file " << d_fileName << "\n";
    exit(1);
  }
}

}  // namespace RDKit::SynthonSpaceSearch
