//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/SynthonSpaceSearch/MemoryMappedFileReader.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceFingerprintSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceRascalSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSubstructureSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>
#include <RDGeneral/StreamOps.h>

namespace RDKit::SynthonSpaceSearch {

// used for serialization
constexpr int32_t versionMajor = 3;
constexpr int32_t versionMinor = 0;
constexpr int32_t endianId = 0xa100f;

size_t SynthonSpace::getNumReactions() const { return d_reactions.size(); }

std::vector<std::string> SynthonSpace::getReactionNames() const {
  std::vector<std::string> reactionNames;
  reactionNames.reserve(d_reactions.size());
  for (const auto &reaction : d_reactions) {
    reactionNames.push_back(reaction.first);
  }
  return reactionNames;
}

const std::shared_ptr<SynthonSet> SynthonSpace::getReaction(
    std::string reactionName) {
  std::pair<std::string, std::shared_ptr<SynthonSet>> tmp =
      std::make_pair(reactionName, std::shared_ptr<SynthonSet>());
  if (const auto &it = std::lower_bound(
          d_reactions.begin(), d_reactions.end(), tmp,
          [](const std::pair<std::string, std::shared_ptr<SynthonSet>> &p1,
             const std::pair<std::string, std::shared_ptr<SynthonSet>> &p2)
              -> bool { return p1.first < p2.first; });
      it != d_reactions.end()) {
    return it->second;
  }
  throw std::runtime_error("Could not find synthon set for reaction " +
                           reactionName);
}

unsigned int SynthonSpace::getPatternFPSize() const {
  PRECONDITION(d_reactions.size(), "No synthon sets available.");
  for (const auto &[id, reaction] : d_reactions) {
    if (!reaction->getSynthons().empty()) {
      return reaction->getSynthons()
          .front()
          .front()
          .second->getPattFP()
          ->getNumBits();
    }
  }
  throw std::runtime_error(
      "Could not find pattern fingerprint for any synthon.");
}

unsigned int SynthonSpace::getFPSize() const {
  PRECONDITION(d_reactions.size(), "No synthon sets available.");
  for (const auto &[id, reaction] : d_reactions) {
    if (!reaction->getSynthons().empty()) {
      return reaction->getSynthons()
          .front()
          .front()
          .second->getFP()
          ->getNumBits();
    }
  }
  throw std::runtime_error("Could not find fingerprint for any synthon.");
}

std::uint64_t SynthonSpace::getNumProducts() const { return d_numProducts; }

std::string SynthonSpace::getInputFileName() const { return d_fileName; }

SearchResults SynthonSpace::substructureSearch(
    const ROMol &query, const SubstructMatchParameters &matchParams,
    const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  ControlCHandler::reset();
  SynthonSpaceSubstructureSearcher ssss(query, matchParams, params, *this);
  return ssss.search();
}

SearchResults SynthonSpace::substructureSearch(
    const GeneralizedSubstruct::ExtendedQueryMol &query,
    const SubstructMatchParameters &matchParams,
    const SynthonSpaceSearchParams &params) {
  if (std::holds_alternative<GeneralizedSubstruct::ExtendedQueryMol::RWMol_T>(
          query.xqmol)) {
    return substructureSearch(
        *std::get<GeneralizedSubstruct::ExtendedQueryMol::RWMol_T>(query.xqmol),
        matchParams, params);
#ifdef RDK_USE_BOOST_SERIALIZATION
  } else if (std::holds_alternative<
                 GeneralizedSubstruct::ExtendedQueryMol::MolBundle_T>(
                 query.xqmol)) {
    return extendedSearch(
        *std::get<GeneralizedSubstruct::ExtendedQueryMol::MolBundle_T>(
            query.xqmol),
        matchParams, params);
  } else if (std::holds_alternative<
                 GeneralizedSubstruct::ExtendedQueryMol::TautomerQuery_T>(
                 query.xqmol)) {
    return extendedSearch(
        *std::get<GeneralizedSubstruct::ExtendedQueryMol::TautomerQuery_T>(
            query.xqmol),
        matchParams, params);
  } else if (std::holds_alternative<
                 GeneralizedSubstruct::ExtendedQueryMol::TautomerBundle_T>(
                 query.xqmol)) {
    return extendedSearch(
        std::get<GeneralizedSubstruct::ExtendedQueryMol::TautomerBundle_T>(
            query.xqmol),
        matchParams, params);
  }
#endif
  else {
    UNDER_CONSTRUCTION("unrecognized type in ExtendedQueryMol");
  }
  return SearchResults();
}

SearchResults SynthonSpace::fingerprintSearch(
    const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
    const SynthonSpaceSearchParams &params) {
  ControlCHandler::reset();
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceFingerprintSearcher ssss(query, fpGen, params, *this);
  return ssss.search();
}

SearchResults SynthonSpace::rascalSearch(
    const ROMol &query, const RascalMCES::RascalOptions &rascalOptions,
    const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceRascalSearcher ssss(query, rascalOptions, params, *this);
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

void SynthonSpace::readTextFile(const std::string &inFilename,
                                bool &cancelled) {
  d_fileName = inFilename;
  std::ifstream ifs(d_fileName);
  if (!ifs.is_open() || ifs.bad()) {
    throw std::runtime_error("Couldn't open file " + d_fileName);
  }

  int format = -1;
  std::string nextLine;
  int lineNum = 1;
  ControlCHandler::reset();

  while (!ifs.eof()) {
    if (ControlCHandler::getGotSignal()) {
      cancelled = true;
      return;
    }
    auto nextSynthon = readSynthonLine(ifs, lineNum, format, d_fileName);
    if (nextSynthon.empty()) {
      continue;
    }
    auto currReaction = addReactionToPool(nextSynthon[3]);
    fixConnectors(nextSynthon[0]);
    auto synthonNum = getSynthonNum(format, nextSynthon[2]);
    auto newSynth = addSynthonToPool(nextSynthon[0]);
    currReaction->addSynthon(synthonNum, newSynth, nextSynthon[1]);
  }
  // Do some final processing.
  for (auto &[id, reaction] : d_reactions) {
    reaction->removeEmptySynthonSets();
    reaction->makeSynthonSearchMols();
    reaction->buildConnectorRegions();
    reaction->assignConnectorsUsed();
    d_numProducts += reaction->getNumProducts();
    if (reaction->getSynthons().size() > d_maxNumSynthons) {
      d_maxNumSynthons = reaction->getSynthons().size();
    }
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
  streamWrite(os, static_cast<std::uint64_t>(d_synthonPool.size()));
  streamWrite(os, static_cast<std::uint64_t>(d_reactions.size()));
  streamWrite(os, d_numProducts);

  // In case we ever want to do random access file reading, save the
  // positions of the synthons and reactions.  Put dummy positions
  // in for now then overwrite at the end.
  std::vector<std::uint64_t> synthonPos(d_synthonPool.size(), std::uint64_t(0));
  std::vector<std::uint64_t> reactionPos(d_reactions.size(), std::uint64_t(0));
  std::streampos synthonPoolStart = os.tellp();
  for (const auto p : synthonPos) {
    streamWrite(os, p);
  }
  for (const auto p : reactionPos) {
    streamWrite(os, p);
  }
  size_t synthonNum = 0;
  for (auto &[smiles, synthon] : d_synthonPool) {
    synthonPos[synthonNum++] = os.tellp();
    synthon->writeToDBStream(os);
  }
  size_t reactionNum = 0;
  for (auto &[id, reaction] : d_reactions) {
    reactionPos[reactionNum++] = os.tellp();
    reaction->writeToDBStream(os);
  }
  // Now write the positions properly
  os.seekp(synthonPoolStart);
  for (const auto p : synthonPos) {
    streamWrite(os, p);
  }
  for (const auto p : reactionPos) {
    streamWrite(os, p);
  }

  os.close();
}

namespace {
// Read the given synthons into array.  synthons is expected to be
// large enough to accept everything.
void readSynthons(
    const size_t startNum, size_t endNum, const char *fileMap,
    const std::vector<std::uint64_t> &synthonPos,
    std::vector<std::pair<std::string, std::unique_ptr<Synthon>>> &synthons) {
  if (endNum > synthons.size()) {
    endNum = synthons.size();
  }
  for (size_t i = startNum; i < endNum; i++) {
    std::string view(fileMap + synthonPos[i],
                     synthonPos[i + 1] - synthonPos[i]);
    std::istringstream is(view, std::ios::binary);
    auto tmp =
        std::make_pair(std::string(), std::unique_ptr<Synthon>(new Synthon));
    tmp.second->readFromDBStream(is);
    tmp.first = tmp.second->getSmiles();
    synthons[i] = std::move(tmp);
  }
}

void threadedReadSynthons(
    const char *fileMap, const std::vector<std::uint64_t> &synthonPos,
    unsigned int numThreads,
    std::vector<std::pair<std::string, std::unique_ptr<Synthon>>> &synthons) {
  size_t eachThread = 1 + (synthonPos.size() / numThreads);
  size_t start = 0;
  std::vector<std::thread> threads;
  for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
    threads.push_back(std::thread(readSynthons, start, start + eachThread,
                                  fileMap, std::ref(synthonPos),
                                  std::ref(synthons)));
  }
  for (auto &t : threads) {
    t.join();
  }
}

void readReactions(
    const size_t startNum, size_t endNum, const char *fileMap,
    const std::vector<std::uint64_t> &reactionPos, const SynthonSpace &space,
    std::uint32_t version,
    std::vector<std::pair<std::string, std::shared_ptr<SynthonSet>>>
        &reactions) {
  if (endNum > reactions.size()) {
    endNum = reactions.size();
  }
  for (size_t i = startNum; i < endNum; i++) {
    std::string view(fileMap + reactionPos[i],
                     reactionPos[i + 1] - reactionPos[i]);
    std::istringstream is(view, std::ios::binary);
    reactions[i] =
        std::make_pair(std::string(), std::make_shared<SynthonSet>());
    reactions[i].second->readFromDBStream(is, space, version);
    reactions[i].first = reactions[i].second->getId();
  }
}

void threadedReadReactions(
    const char *fileMap, const std::vector<std::uint64_t> &reactionPos,
    unsigned int numThreads, const SynthonSpace &space, std::uint32_t version,
    std::vector<std::pair<std::string, std::shared_ptr<SynthonSet>>>
        &reactions) {
  size_t eachThread = 1 + (reactionPos.size() / numThreads);
  size_t start = 0;
  std::vector<std::thread> threads;
  for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
    threads.push_back(std::thread(
        readReactions, start, start + eachThread, fileMap,
        std::ref(reactionPos), std::ref(space), version, std::ref(reactions)));
  }
  for (auto &t : threads) {
    t.join();
  }
}

}  // namespace

void SynthonSpace::readDBFile(const std::string &inFilename,
                              [[maybe_unused]] int numThreads) {
  d_fileName = inFilename;
  std::ifstream is(inFilename, std::fstream::binary);
  if (!is.is_open() || is.bad()) {
    throw std::runtime_error("Couldn't open file " + d_fileName);
  }
  int32_t endianTest;
  streamRead(is, endianTest);
  if (endianTest != endianId) {
    throw std::runtime_error("Endianness mismatch in SynthonSpace file " +
                             d_fileName);
  }
  streamRead(is, d_fileMajorVersion);
  int32_t minorVersion;
  streamRead(is, minorVersion);
  if (d_fileMajorVersion > versionMajor ||
      (d_fileMajorVersion == versionMajor && minorVersion > versionMinor)) {
    BOOST_LOG(rdWarningLog)
        << "Deserializing from a version number (" << d_fileMajorVersion << "."
        << minorVersion << ")"
        << "that is higher than our version (" << versionMajor << "."
        << versionMinor << ").\nThis probably won't work." << std::endl;
  }
  // version sanity checking
  if (d_fileMajorVersion > 1000 || minorVersion > 100) {
    throw std::runtime_error("unreasonable version numbers");
  }
  d_fileMajorVersion = 1000 * d_fileMajorVersion + minorVersion * 10;
  if (d_fileMajorVersion < 3000) {
    throw std::runtime_error(
        "This binary file version is no longer supported."
        "  Please re-build with a recent version of the RDKit.");
  }

  bool hasFPs;
  streamRead(is, hasFPs);
  if (hasFPs) {
    streamRead(is, d_fpType, 0);
  }
  std::uint64_t numSynthons;
  std::uint64_t numReactions;
  streamRead(is, numSynthons);
  streamRead(is, numReactions);
  streamRead(is, d_numProducts);

  std::vector<std::uint64_t> synthonPos(numSynthons, std::uint64_t(0));
  for (std::uint64_t i = 0; i < numSynthons; i++) {
    streamRead(is, synthonPos[i]);
  }
  std::vector<std::uint64_t> reactionPos(numReactions, std::uint64_t(0));
  for (std::uint64_t i = 0; i < numReactions; i++) {
    streamRead(is, reactionPos[i]);
  }
  is.close();

  details::MemoryMappedFileReader fileMap(d_fileName);
  // put the end of the last synthon and reaction into their respective arrays,
  synthonPos.push_back(reactionPos[0]);
  reactionPos.push_back(fileMap.d_size);
  d_synthonPool.resize(numSynthons);
#if RDK_BUILD_THREADSAFE_SSS
  unsigned int numThreadsToUse = getNumThreadsToUse(numThreads);
  if (numThreadsToUse > 1) {
    threadedReadSynthons(fileMap.d_mappedMemory, synthonPos, numThreadsToUse,
                         d_synthonPool);
  } else {
    readSynthons(0, numSynthons, fileMap.d_mappedMemory, synthonPos,
                 d_synthonPool);
  }
#else
    readSynthons(0, numSynthons, fileMap.d_mappedMemory, synthonPos,
                 d_synthonPool);
#endif
  if (!std::is_sorted(
          d_synthonPool.begin(), d_synthonPool.end(),
          [](const std::pair<std::string, std::unique_ptr<Synthon>> &p1,
             const std::pair<std::string, std::unique_ptr<Synthon>> &p2)
              -> bool { return p1.first < p2.first; })) {
    std::sort(
        d_synthonPool.begin(), d_synthonPool.end(),
        [](const std::pair<std::string, std::unique_ptr<Synthon>> &p1,
           const std::pair<std::string, std::unique_ptr<Synthon>> &p2) -> bool {
          return p1.first < p2.first;
        });
  }
  d_reactions.resize(numReactions);
#if RDK_BUILD_THREADSAFE_SSS
  if (numThreadsToUse > 1) {
    threadedReadReactions(fileMap.d_mappedMemory, reactionPos, numThreadsToUse,
                          *this, d_fileMajorVersion, d_reactions);
  } else {
    readReactions(0, numReactions, fileMap.d_mappedMemory, reactionPos, *this,
                  d_fileMajorVersion, d_reactions);
  }
#else
    readReactions(0, numReactions, fileMap.d_mappedMemory, reactionPos, *this,
                  d_fileMajorVersion, d_reactions);
#endif
  if (!std::is_sorted(
          d_reactions.begin(), d_reactions.end(),
          [](const std::pair<std::string, std::shared_ptr<SynthonSet>> &p1,
             const std::pair<std::string, std::shared_ptr<SynthonSet>> &p2)
              -> bool { return p1.first < p2.first; })) {
    std::sort(d_reactions.begin(), d_reactions.end(),
              [](const std::pair<std::string, std::shared_ptr<SynthonSet>> &p1,
                 const std::pair<std::string, std::shared_ptr<SynthonSet>> &p2)
                  -> bool { return p1.first < p2.first; });
  }
  for (const auto &[id, reaction] : d_reactions) {
    if (reaction->getSynthons().size() > d_maxNumSynthons) {
      d_maxNumSynthons = reaction->getSynthons().size();
    }
  }
}

void SynthonSpace::summarise(std::ostream &os) {
  os << "Read from file " << d_fileName << "\n"
     << "Number of reactions : " << d_reactions.size() << "\n";
  int synthCounts[MAX_CONNECTOR_NUM + 1]{0, 0, 0, 0, 0};
  for (const auto &id : getReactionNames()) {
    auto rxn = getReaction(id);
    os << "Reaction name " << id << "\n";
    for (size_t i = 0; i < rxn->getSynthons().size(); ++i) {
      os << "  Synthon set " << i << " has " << rxn->getSynthons()[i].size()
         << " synthons"
         << "\n";
    }
    synthCounts[rxn->getSynthons().size()]++;
  }
  os << "Approximate number of molecules in SynthonSpace : "
     << formattedIntegerString(getNumProducts()) << std::endl;
  os << "Number of unique synthons : " << d_synthonPool.size() << std::endl;
  for (unsigned int i = 0; i < MAX_CONNECTOR_NUM; ++i) {
    if (synthCounts[i] > 0) {
      os << "Number of " << i << " molecule reactions : " << synthCounts[i]
         << std::endl;
    }
  }
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
    unsigned int numBits = 0;
    for (const auto &[id, synthSet] : d_reactions) {
      if (ControlCHandler::getGotSignal()) {
        return;
      }
      synthSet->buildSynthonFingerprints(fpGen);
      if (!numBits) {
        numBits = synthSet->getSynthons()
                      .front()
                      .front()
                      .second->getFP()
                      ->getNumBits();
      }
      synthSet->buildAddAndSubtractFPs(fpGen, numBits);
    }
  }
}

bool SynthonSpace::hasAddAndSubstractFingerprints() const {
  if (d_reactions.empty()) {
    return false;
  }
  return d_reactions.begin()->second->hasAddAndSubtractFPs();
}

Synthon *SynthonSpace::addSynthonToPool(const std::string &smiles) {
  // Clearly this is going to be inefficient, but it won't be done
  // except when reading a Text file which should only be done
  // occasionally, when converting a new database to binary.
  auto tmp = std::make_pair(smiles, std::unique_ptr<Synthon>());
  if (auto it = std::lower_bound(
          d_synthonPool.begin(), d_synthonPool.end(), tmp,
          [](const std::pair<std::string, std::unique_ptr<Synthon>> &p1,
             const std::pair<std::string, std::unique_ptr<Synthon>> &p2)
              -> bool { return p1.first < p2.first; });
      it != d_synthonPool.end() && it->first == smiles) {
    return it->second.get();
  } else {
    tmp.second.reset(new Synthon(smiles));
    auto retVal = tmp.second.get();
    d_synthonPool.insert(it, std::move(tmp));
    return retVal;
  }
}

std::shared_ptr<SynthonSet> SynthonSpace::addReactionToPool(
    const std::string &reactionName) {
  std::pair<std::string, std::shared_ptr<SynthonSet>> tmp =
      std::make_pair(reactionName, std::shared_ptr<SynthonSet>());
  if (const auto &it = std::lower_bound(
          d_reactions.begin(), d_reactions.end(), tmp,
          [](const std::pair<std::string, std::shared_ptr<SynthonSet>> &p1,
             const std::pair<std::string, std::shared_ptr<SynthonSet>> &p2)
              -> bool { return p1.first < p2.first; });
      it != d_reactions.end() && it->first == reactionName) {
    return it->second;
  } else {
    tmp.second.reset(new SynthonSet(reactionName));
    d_reactions.insert(it, tmp);
    return tmp.second;
  }
}

Synthon *SynthonSpace::getSynthonFromPool(const std::string &smiles) const {
  auto tmp = std::make_pair(smiles, std::unique_ptr<Synthon>());
  if (auto it = std::lower_bound(
          d_synthonPool.begin(), d_synthonPool.end(), tmp,
          [](const std::pair<std::string, std::unique_ptr<Synthon>> &p1,
             const std::pair<std::string, std::unique_ptr<Synthon>> &p2)
              -> bool { return p1.first < p2.first; });
      it != d_synthonPool.end()) {
    return it->second.get();
  }
  return nullptr;
}

SearchResults SynthonSpace::extendedSearch(
    const MolBundle &query, const SubstructMatchParameters &matchParams,
    const SynthonSpaceSearchParams &params) {
  SearchResults results;
  SynthonSpaceSearchParams tmpParams(params);
  for (unsigned int i = 0; i < query.size(); ++i) {
    auto theseResults = substructureSearch(*query[i], matchParams, tmpParams);
    if (tmpParams.maxHits != -1) {
      tmpParams.maxHits -= theseResults.getHitMolecules().size();
    }
    results.mergeResults(theseResults);
  }
  return results;
}

SearchResults SynthonSpace::extendedSearch(
    const GeneralizedSubstruct::ExtendedQueryMol::TautomerBundle_T &query,
    const SubstructMatchParameters &matchParams,
    const SynthonSpaceSearchParams &params) {
  SearchResults results;
  SynthonSpaceSearchParams tmpParams(params);
  for (const auto &tq : *query) {
    auto theseResults = extendedSearch(*tq, matchParams, params);
    if (tmpParams.maxHits != -1) {
      tmpParams.maxHits -= theseResults.getHitMolecules().size();
    }
    results.mergeResults(theseResults);
  }
  return results;
}

SearchResults SynthonSpace::extendedSearch(
    const TautomerQuery &query, const SubstructMatchParameters &matchParams,
    const SynthonSpaceSearchParams &params) {
  SearchResults results;
  SynthonSpaceSearchParams tmpParams(params);
  for (const auto &tq : query.getTautomers()) {
    auto theseResults = substructureSearch(*tq, matchParams, params);
    if (tmpParams.maxHits != -1) {
      tmpParams.maxHits -= theseResults.getHitMolecules().size();
    }
    results.mergeResults(theseResults);
  }
  return results;
}

void convertTextToDBFile(const std::string &inFilename,
                         const std::string &outFilename, bool &cancelled,
                         const FingerprintGenerator<std::uint64_t> *fpGen) {
  SynthonSpace synthSpace;
  cancelled = false;
  synthSpace.readTextFile(inFilename, cancelled);
  if (ControlCHandler::getGotSignal()) {
    cancelled = true;
    return;
  }
  if (fpGen) {
    synthSpace.buildSynthonFingerprints(*fpGen);
    if (ControlCHandler::getGotSignal()) {
      cancelled = true;
      return;
    }
  }
  synthSpace.writeDBFile(outFilename);
}

namespace {
// Stuff for formatting integers with spaces every 3 digits.
template <class Char>
class MyFacet : public std::numpunct<Char> {
 public:
  std::string do_grouping() const { return "\3"; }
  Char do_thousands_sep() const { return ' '; }
};
}  // namespace

std::string formattedIntegerString(std::int64_t value) {
  std::ostringstream oss;
  oss.imbue(std::locale(oss.getloc(), new MyFacet<char>));
  oss << value;
  return oss.str();
}

}  // namespace RDKit::SynthonSpaceSearch
