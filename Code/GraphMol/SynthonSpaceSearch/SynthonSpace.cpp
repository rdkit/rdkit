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
#include <filesystem>
#include <fstream>
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
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceShapeSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSubstructureSearcher.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SynthonSpaceSearch/ProgressBar.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>
#include <RDGeneral/StreamOps.h>

namespace RDKit {
namespace SynthonSpaceSearch {
// used for serialization
constexpr int32_t versionMajor = 3;
constexpr int32_t versionMinor = 2;
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

std::shared_ptr<SynthonSet> SynthonSpace::getReaction(
    std::string reactionName) {
  std::pair<std::string, std::shared_ptr<SynthonSet>> tmp =
      std::make_pair(reactionName, std::shared_ptr<SynthonSet>());
  if (const auto &it = std::lower_bound(
          d_reactions.begin(), d_reactions.end(), tmp,
          [](const std::pair<std::string, std::shared_ptr<SynthonSet>> &p1,
             const std::pair<std::string, std::shared_ptr<SynthonSet>> &p2)
              -> bool { return p1.first < p2.first; });
      it != d_reactions.end() && it->first == reactionName) {
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
size_t SynthonSpace::getNumSynthons() const { return d_synthonPool.size(); }

std::string SynthonSpace::getInputFileName() const { return d_fileName; }

SearchResults SynthonSpace::substructureSearch(
    const ROMol &query, const SubstructMatchParameters &matchParams,
    const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  ControlCHandler::reset();

  SynthonSpaceSubstructureSearcher ssss(query, matchParams, params, *this);
  return ssss.search(ThreadMode::ThreadReactions);
}

void SynthonSpace::substructureSearch(
    const ROMol &query, const SearchResultCallback &cb,
    const SubstructMatchParameters &matchParams,
    const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  ControlCHandler::reset();
  SynthonSpaceSubstructureSearcher ssss(query, matchParams, params, *this);
  ssss.search(cb, ThreadMode::ThreadReactions);
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
  }
#ifdef RDK_USE_BOOST_SERIALIZATION
  if (std::holds_alternative<
          GeneralizedSubstruct::ExtendedQueryMol::MolBundle_T>(query.xqmol)) {
    return extendedSearch(
        *std::get<GeneralizedSubstruct::ExtendedQueryMol::MolBundle_T>(
            query.xqmol),
        matchParams, params);
  }
  if (std::holds_alternative<
          GeneralizedSubstruct::ExtendedQueryMol::TautomerQuery_T>(
          query.xqmol)) {
    return extendedSearch(
        *std::get<GeneralizedSubstruct::ExtendedQueryMol::TautomerQuery_T>(
            query.xqmol),
        matchParams, params);
  }
  if (std::holds_alternative<
          GeneralizedSubstruct::ExtendedQueryMol::TautomerBundle_T>(
          query.xqmol)) {
    return extendedSearch(
        std::get<GeneralizedSubstruct::ExtendedQueryMol::TautomerBundle_T>(
            query.xqmol),
        matchParams, params);
  }
#endif
  UNDER_CONSTRUCTION("unrecognized type in ExtendedQueryMol");
}

SearchResults SynthonSpace::fingerprintSearch(
    const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
    const SynthonSpaceSearchParams &params) {
  ControlCHandler::reset();
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");

  SynthonSpaceFingerprintSearcher ssss(query, fpGen, params, *this);
  return ssss.search(ThreadMode::ThreadFragments);
}

void SynthonSpace::fingerprintSearch(
    const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
    const SearchResultCallback &cb, const SynthonSpaceSearchParams &params) {
  ControlCHandler::reset();
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceFingerprintSearcher ssss(query, fpGen, params, *this);
  ssss.search(cb, ThreadMode::ThreadFragments);
}

SearchResults SynthonSpace::rascalSearch(
    const ROMol &query, const RascalMCES::RascalOptions &rascalOptions,
    const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceRascalSearcher ssrs(query, rascalOptions, params, *this);
  return ssrs.search(ThreadMode::ThreadFragments);
}

void SynthonSpace::rascalSearch(const ROMol &query,
                                const RascalMCES::RascalOptions &rascalOptions,
                                const SearchResultCallback &cb,
                                const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  SynthonSpaceRascalSearcher ssss(query, rascalOptions, params, *this);
  ssss.search(cb, ThreadMode::ThreadFragments);
}

SearchResults SynthonSpace::shapeSearch(
    const ROMol &query, const SynthonSpaceSearchParams &params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");
  if (!params.enumerateUnspecifiedStereo) {
    SynthonSpaceShapeSearcher ssss(query, params, *this);
    return ssss.search(ThreadMode::ThreadFragments);
  }
  // Otherwise, we need to enumerate any isomers, search each one
  // and accumulate the results.
  auto paramsCp = params;
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.numThreads = params.numThreads;
  dgParams.pruneRmsThresh = params.confRMSThreshold;
  dgParams.randomSeed = params.randomSeed;
  auto isomers = details::generateIsomerConformers(
      query, params.numConformers, params.enumerateUnspecifiedStereo,
      params.stereoEnumOpts, dgParams);
  SearchResults allResults;
  for (auto &isomer : isomers) {
    SynthonSpaceShapeSearcher ssss(*isomer, paramsCp, *this);
    auto results = ssss.search(ThreadMode::ThreadFragments);
    allResults.mergeResults(results);
    if (params.maxHits > 0) {
      if (allResults.getHitMolecules().size() >
          static_cast<std::uint64_t>(params.maxHits)) {
        return allResults;
      }
      paramsCp.maxHits = params.maxHits - allResults.getHitMolecules().size();
    }
  }
  return allResults;
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
  // comma separated 3-5, if tab-separated, 6-8.  Tab-separated takes
  // precedence over space-separated to allow for CXSMILES for the synthons
  // which might have a tab in them.
  static const std::vector<std::vector<std::string>> firstLineOpts{
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id"},
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id",
                               "release"},
      std::vector<std::string>{"SMILES", "synton_id", "synton_role",
                               "reaction_id"}};
  static const std::regex regext("\\t");
  auto lineParts = splitLine(line, regext);
  for (size_t i = 0; i < firstLineOpts.size(); ++i) {
    if (lineParts == firstLineOpts[i]) {
      return static_cast<int>(i + 6);
    }
  }

  // This includes tabs, obvs, but they should already have been detected.
  static const std::regex regexws("\\s+");
  lineParts = splitLine(line, regexws);
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
  static const std::regex regext("\\t+");
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
  } else if (format > 3 && format < 6) {
    nextSynthon = splitLine(nextLine, regexc);
  } else if (format > 5) {
    nextSynthon = splitLine(nextLine, regext);
  }
  if (nextSynthon.size() < 4) {
    throw std::runtime_error("Bad format for SynthonSpace file " + fileName +
                             " on line " + std::to_string(lineNum));
  }
  return nextSynthon;
}

int getSynthonNum(int format, const std::string &synthon) {
  int synthonNum{std::numeric_limits<int>::max()};
  if (format == 0 || format == 1 || format == 3 || format == 4 || format == 6 ||
      format == 7) {
    synthonNum = std::stoi(synthon);
  } else if (format == 2 || format == 5 || format == 8) {
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
  readStream(ifs, cancelled);
}

void SynthonSpace::readStream(std::istream &is, bool &cancelled) {
  int format = -1;
  std::string nextLine;
  int lineNum = 1;
  ControlCHandler::reset();

  while (!is.eof()) {
    if (ControlCHandler::getGotSignal()) {
      cancelled = true;
      return;
    }
    auto nextSynthon = readSynthonLine(is, lineNum, format, d_fileName);
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
    reaction->assessRingFormers();
    if (reaction->getNumRingFormers()) {
      d_hasRingFormer = true;
    }
    d_numProducts += reaction->getNumProducts();
    if (reaction->getSynthons().size() > d_maxNumSynthons) {
      d_maxNumSynthons = reaction->getSynthons().size();
    }
    reaction->initializeSearchOrders();
  }
  fillSynthonReactions();
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
  streamWrite(os, d_numConformers);
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
    const std::vector<std::uint64_t> &synthonPos, std::uint32_t version,
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
    tmp.second->readFromDBStream(is, version);
    tmp.first = tmp.second->getSmiles();
    synthons[i] = std::move(tmp);
  }
}

void threadedReadSynthons(
    const char *fileMap, const std::vector<std::uint64_t> &synthonPos,
    unsigned int numThreads, std::uint32_t version,
    std::vector<std::pair<std::string, std::unique_ptr<Synthon>>> &synthons) {
  size_t eachThread = 1 + (synthonPos.size() / numThreads);
  size_t start = 0;
  std::vector<std::thread> threads;
  for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
    threads.push_back(std::thread(readSynthons, start, start + eachThread,
                                  fileMap, std::ref(synthonPos), version,
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
  if (d_fileMajorVersion > 3010) {
    streamRead(is, d_numConformers);
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
  unsigned int numThreadsToUse = getNumThreadsToUse(numThreads);
  if (numThreadsToUse > 1) {
    threadedReadSynthons(fileMap.d_mappedMemory, synthonPos, numThreadsToUse,
                         d_fileMajorVersion, d_synthonPool);
  } else {
    readSynthons(0, numSynthons, fileMap.d_mappedMemory, synthonPos,
                 d_fileMajorVersion, d_synthonPool);
  }
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
  if (numThreadsToUse > 1) {
    threadedReadReactions(fileMap.d_mappedMemory, reactionPos, numThreadsToUse,
                          *this, d_fileMajorVersion, d_reactions);
  } else {
    readReactions(0, numReactions, fileMap.d_mappedMemory, reactionPos, *this,
                  d_fileMajorVersion, d_reactions);
  }
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
    reaction->assessRingFormers();
    if (reaction->getNumRingFormers()) {
      d_hasRingFormer = true;
    }
    if (reaction->getSynthons().size() > d_maxNumSynthons) {
      d_maxNumSynthons = reaction->getSynthons().size();
    }
  }
  fillSynthonReactions();
  BOOST_LOG(rdInfoLog) << "Number of synthons : " << d_synthonPool.size()
                       << " of which " << getNumSynthonsWithShapes()
                       << " have shapes.\n";
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
  enumerateToStream(os);
  os.close();
}

void SynthonSpace::enumerateToStream(std::ostream &os) const {
  for (const auto &[fst, snd] : d_reactions) {
    snd->enumerateToStream(os);
  }
}

unsigned int SynthonSpace::getMaxNumConnectors() const {
  unsigned int maxNumConns = 0;
  for (const auto &[id, synSet] : d_reactions) {
    maxNumConns =
        std::max(maxNumConns,
                 static_cast<unsigned int>(synSet->getConnectors().count()));
  }
  return maxNumConns;
}

bool SynthonSpace::hasFingerprints() const { return !d_fpType.empty(); }
unsigned int SynthonSpace::getNumConformers() const { return d_numConformers; }

void SynthonSpace::buildSynthonFingerprints(
    const FingerprintGenerator<std::uint64_t> &fpGen,
    unsigned int progressBarWidth) {
  BOOST_LOG(rdWarningLog) << "Building the fingerprints of "
                          << d_synthonPool.size()
                          << " synthons may take some time." << std::endl;

  if (const auto fpType = fpGen.infoString();
      fpType != d_fpType || !hasFingerprints()) {
    BOOST_LOG(rdWarningLog)
        << "Building the fingerprints may take some time." << std::endl;
    std::unique_ptr<ProgressBar> pbar;
    if (progressBarWidth) {
      pbar.reset(new ProgressBar(progressBarWidth, d_reactions.size()));
    }
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
      if (pbar) {
        pbar->increment();
      }
    }
    d_fpType = fpGen.infoString();
    for (auto &rxn : d_reactions) {
      rxn.second->initializeSearchOrders();
    }
  }
}

namespace {
void writeInterimFile(const SynthonSpace &space, const std::string &filename) {
  // Write to a temporary file and then move it onto filename.  That way
  // if there's a crash during the write, we reduce the chance of losing
  // everything.
  const std::string tempFile = filename + ".tmp";
  space.writeDBFile(tempFile);
  std::filesystem::rename(tempFile, filename);
}
}  // namespace

void SynthonSpace::buildSynthonShapes(bool &cancelled,
                                      const ShapeBuildParams &shapeParams) {
  if (d_numConformers == shapeParams.numConfs) {
    bool missingShapes = false;
    for (const auto &[id, synthon] : d_synthonPool) {
      if (!synthon->getShapes()) {
        missingShapes = true;
        break;
      }
    }
    if (missingShapes) {
      BOOST_LOG(rdWarningLog)
          << "Resuming building SynthonSpace shapes." << std::endl;
    } else {
      BOOST_LOG(rdWarningLog)
          << "SynthonSpace has already been built with " << shapeParams.numConfs
          << " conformers." << std::endl;
      return;
    }
  }
  BOOST_LOG(rdWarningLog) << "Building the conformers of "
                          << d_synthonPool.size()
                          << " synthons may take some time." << std::endl;
  if (d_synthonReactions.empty()) {
    fillSynthonReactions();
  }

  d_numConformers = shapeParams.numConfs;
  std::vector<unsigned int> doneRxns(d_synthonPool.size(), 0u);
  cancelled = false;
  std::vector<std::vector<std::unique_ptr<SampleMolRec>>> allSampleMols;
  buildSynthonSampleMolecules(shapeParams.maxSynthonAtoms, allSampleMols);
  // Do the big molecules first so there is less chance of sitting at the
  // end waiting for 1 big embedding to finish.
  std::ranges::sort(allSampleMols,
                    [](const auto &sm1, const auto &sm2) -> bool {
                      return sm1.front()->d_numAtoms > sm2.front()->d_numAtoms;
                    });

  bool interimWrite = true;
  if (shapeParams.interimWrites == 0 || shapeParams.interimFile.empty()) {
    interimWrite = false;
  }
  std::unique_ptr<ProgressBar> pbar;
  if (shapeParams.useProgressBar) {
    pbar.reset(
        new ProgressBar(shapeParams.useProgressBar, allSampleMols.size()));
  }

  while (!cancelled) {
    // Loop around until all synthons have some shapes or we've run out
    // of synthon/reaction combinations.
    std::vector<std::unique_ptr<SampleMolRec>> sampleMols;
    sampleMols.reserve(d_synthonPool.size());
    for (auto &allSampleMol : allSampleMols) {
      if (!allSampleMol.empty()) {
        // If we have a shape object in the synthon with some shapes, don't do
        // anything.  If there's a shape object but no shapes then probably
        // the embedding failed last time, so try again.  It might be the
        // fault of this synthon, but it might conceivably have been a problem
        // with whatever it was attached to, and another reaction might have
        // bolted on something more amenable.
        if (!allSampleMol.back()->d_synthon->getShapes() ||
            allSampleMol.back()->d_synthon->getShapes()->hasNoShapes()) {
          sampleMols.push_back(std::move(allSampleMol.back()));
          allSampleMol.pop_back();
        }
      }
      if (interimWrite && shapeParams.interimWrites == sampleMols.size()) {
        break;
      }
    }
    if (sampleMols.empty()) {
      break;
    }
    std::ranges::sort(sampleMols, [](const auto &a, const auto &b) -> bool {
      return a->d_numAtoms > b->d_numAtoms;
    });
    auto dgParams = DGeomHelpers::KDG;
    dgParams.numThreads = 1;
    dgParams.pruneRmsThresh = shapeParams.rmsThreshold;
    dgParams.randomSeed = shapeParams.randomSeed;
    dgParams.maxIterations = shapeParams.maxEmbedAttempts;
    dgParams.timeout = shapeParams.timeOut;
    auto numWithShapesB4 = getNumSynthonsWithShapes();
    details::makeShapesFromMols(sampleMols, dgParams, shapeParams, pbar);
    auto numWithShapes = getNumSynthonsWithShapes();
    if (interimWrite && numWithShapes > numWithShapesB4) {
      BOOST_LOG(rdInfoLog) << "Writing interim file " << shapeParams.interimFile
                           << " with " << numWithShapes << " shapes."
                           << std::endl;
      writeInterimFile(*this, shapeParams.interimFile);
    }
    if (ControlCHandler::getGotSignal()) {
      cancelled = true;
    }
  }
}

std::uint64_t SynthonSpace::getNumSynthonsWithShapes() const {
  size_t numShapes = 0;
  for (const auto &[smiles, synthon] : d_synthonPool) {
    if (synthon->getShapes() && !synthon->getShapes()->hasNoShapes()) {
      numShapes++;
    }
  }
  return numShapes;
}

void SynthonSpace::reportSynthonUsage(std::ostream &os) const {
  size_t maxSize = 0;
  for (const auto &it : d_synthonReactions) {
    if (it.second.size() > maxSize) {
      maxSize = it.second.size();
    }
  }
  os << "Max number of reactions for a synthon: " << maxSize << std::endl;
  std::vector<size_t> counts(maxSize + 1, 0);
  for (const auto &it : d_synthonReactions) {
    counts[it.second.size()]++;
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
      it != d_synthonPool.end() && it->first == smiles) {
    return it->second.get();
  }
  return nullptr;
}

void SynthonSpace::orderSynthonsForSearch() {
  for (auto &[name, reaction] : d_reactions) {
    reaction->initializeSearchOrders();
  }
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

void SynthonSpace::fillSynthonReactions() {
  for (auto &[id, rxn] : d_reactions) {
    auto rxnSynthons = rxn->getSynthons();
    for (auto &synthonSet : rxnSynthons) {
      for (auto &[id, synthon] : synthonSet) {
        auto smiles = synthon->getSmiles();
        if (auto it = d_synthonReactions.find(smiles);
            it != d_synthonReactions.end()) {
          synthon->updateMaxSynthonSetSize(rxn->getSynthons().size());
          it->second.push_back(rxn.get());
        } else {
          d_synthonReactions.insert(
              std::make_pair(smiles, std::vector<SynthonSet *>(1, rxn.get())));
          synthon->updateMaxSynthonSetSize(rxn->getSynthons().size());
        }
      }
    }
  }
}

void SynthonSpace::buildSynthonSampleMolecules(
    unsigned int maxSynthonAtoms,
    std::vector<std::vector<std::unique_ptr<SampleMolRec>>> &sampleMols) const {
  sampleMols.reserve(d_synthonReactions.size());
  for (const auto &[synthonSmiles, reactions] : d_synthonReactions) {
    auto synthon = getSynthonFromPool(synthonSmiles);
    if ((maxSynthonAtoms && synthon->getNumHeavyAtoms() > maxSynthonAtoms) ||
        (synthon->getShapes() && !synthon->getShapes()->hasNoShapes())) {
      continue;
    }
    std::vector<std::unique_ptr<SampleMolRec>> theseSamples;
    theseSamples.reserve(reactions.size());
    for (auto &reaction : reactions) {
      theseSamples.push_back(reaction->makeSampleMolecule(synthon));
      // In the unlikely event that we didn't get anything, drop it.
      if (!theseSamples.back()->d_numAtoms) {
        theseSamples.pop_back();
      }
    }
    std::ranges::sort(theseSamples, [](const auto &a, const auto &b) -> bool {
      return a->d_numAtoms > b->d_numAtoms;
    });
    sampleMols.push_back(std::move(theseSamples));
  }
}

void convertTextToDBFile(const std::string &inFilename,
                         const std::string &outFilename, bool &cancelled,
                         const FingerprintGenerator<std::uint64_t> *fpGen,
                         const ShapeBuildParams *options) {
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
  if (options) {
    synthSpace.buildSynthonShapes(cancelled, *options);
    if (cancelled) {
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
  std::string do_grouping() const override { return "\3"; }
  Char do_thousands_sep() const override { return ' '; }
};
}  // namespace

std::string formattedIntegerString(std::int64_t value) {
  std::ostringstream oss;
  oss.imbue(std::locale(oss.getloc(), new MyFacet<char>));
  oss << value;
  return oss.str();
}
}  // namespace SynthonSpaceSearch
}  // namespace RDKit
