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
#include <GraphMol/QueryOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

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

SubstructureResults SynthonSpace::substructureSearch(
    const ROMol &query, SynthonSpaceSearchParams params) {
  PRECONDITION(query.getNumAtoms() != 0, "Search query must contain atoms.");

  if (params.randomSample) {
    if (!d_randGen) {
      std::random_device rd;
      d_randGen = std::make_unique<std::mt19937>(rd());
    }
    if (params.randomSeed != -1) {
      d_randGen->seed(params.randomSeed);
    }
  }
  std::vector<std::unique_ptr<ROMol>> results;
  RDKit::MatchVectType dontCare;

  auto fragments = details::splitMolecule(query, params.maxBondSplits);
  std::vector<SynthonSpaceHitSet> allHits;
  size_t totHits = 0;
  for (auto &fragSet : fragments) {
    auto theseHits = searchFragSet(fragSet);
    if (!theseHits.empty()) {
      totHits += std::accumulate(
          theseHits.begin(), theseHits.end(), 0,
          [](const size_t prevVal, const SynthonSpaceHitSet &hs) -> size_t {
            return prevVal + hs.numHits;
          });
      allHits.insert(allHits.end(), theseHits.begin(), theseHits.end());
    }
  }

  if (params.buildHits) {
    std::sort(allHits.begin(), allHits.end(),
              [](const SynthonSpaceHitSet &hs1,
                 const SynthonSpaceHitSet &hs2) -> bool {
                if (hs1.reactionId == hs2.reactionId) {
                  return hs1.numHits < hs2.numHits;
                } else {
                  return hs1.reactionId < hs2.reactionId;
                }
              });
    // Keep track of the result names so we can weed out duplicates by
    // reaction and synthons.  Different splits may give rise to the same
    // synthon combination.  This will keep the same molecule produced via
    // different reactions which I think makes sense.  The resultsNames will
    // be accumulated even if the molecule itself doesn't make it into the
    // results set, for example if it isn't a random selection or it's
    // outside maxHits or hitStart.
    std::set<std::string> resultsNames;
    // The random sampling, if requested, doesn't always produce the required
    // number of hits in 1 sweep.  Also, totHits, being an upper bound, may
    // not be achievable.
    size_t tmpMaxHits(params.maxHits);
    int numSweeps = 0;
    while (results.size() < std::min(tmpMaxHits, totHits) &&
           numSweeps < params.numRandomSweeps) {
      buildHits(allHits, query, params, totHits, resultsNames, results);
      // totHits is an upper bound, so may not be reached.
      if (params.maxHits == -1 || !params.randomSample) {
        break;
      }
      numSweeps++;
    }
  }
  return SubstructureResults{std::move(results), totHits};
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

// Return a bitset giving the different connector types in this
// molecule.
boost::dynamic_bitset<> getConnectorPattern(
    const std::vector<std::unique_ptr<ROMol>> &fragSet) {
  boost::dynamic_bitset<> conns(MAX_CONNECTOR_NUM);
  for (const auto &frag : fragSet) {
    for (const auto &a : frag->atoms()) {
      if (!a->getAtomicNum() && a->getIsotope()) {
        conns.set(a->getIsotope() - 1);
      }
    }
  }
  return conns;
}

// Make pattern fps for the fragments and re-order fps and fragments into
// descending order of number of bits set in the fp (the "largest fragment
// heuristic").  The FP with the largest number of bits
// is the most likely to screen out a matching synthon set since smaller,
// less complex fragments are more likely to match something, so screen
// with that first.
std::vector<std::unique_ptr<ExplicitBitVect>> makePatternFPs(
    std::vector<std::unique_ptr<ROMol>> &molFrags) {
  std::vector<std::unique_ptr<ExplicitBitVect>> pattFPs;
  pattFPs.reserve(molFrags.size());
  for (const auto &frag : molFrags) {
    pattFPs.emplace_back(PatternFingerprintMol(*frag, 2048));
  }
  // Sort by descending number of bits set.
  std::vector<std::pair<size_t, ExplicitBitVect *>> fps(pattFPs.size());
  for (size_t i = 0; i < pattFPs.size(); ++i) {
    fps[i] = std::make_pair(i, pattFPs[i].get());
  }
  std::sort(fps.begin(), fps.end(),
            [](const std::pair<size_t, ExplicitBitVect *> &fp1,
               const std::pair<size_t, ExplicitBitVect *> &fp2) -> bool {
              return fp1.second->getNumOnBits() > fp2.second->getNumOnBits();
            });

  // Now put molFrags in the same order.
  std::vector<std::unique_ptr<ROMol>> newFrags(molFrags.size());
  std::vector<std::unique_ptr<ExplicitBitVect>> retFPs(molFrags.size());
  for (size_t i = 0; i < fps.size(); ++i) {
    newFrags[i] = std::move(molFrags[fps[i].first]);
    retFPs[i] = std::move(pattFPs[fps[i].first]);
  }
  molFrags = std::move(newFrags);
  return retFPs;
}

// Return copies of the mol fragments will all permutations of the connectors
// in the reaction onto the connectors in the fragments.
// E.g. if the reaction has 3 connectors, 1, 2 and 3 and the fragged mol has
// 2, return all permutations of 2 from 3.  It's ok if the fragged mol doesn't
// have all the connections in the reaction, although this may well result in
// a lot of hits.
std::vector<std::vector<std::unique_ptr<RWMol>>> getConnectorPermutations(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    const boost::dynamic_bitset<> &fragConns,
    const boost::dynamic_bitset<> &reactionConns) {
  std::vector<std::vector<std::unique_ptr<RWMol>>> connPerms;
  auto bitsToInts =
      [](const boost::dynamic_bitset<> &bits) -> std::vector<int> {
    std::vector<int> ints;
    for (size_t i = 0; i < bits.size(); ++i) {
      if (bits[i]) {
        ints.push_back(static_cast<int>(i));
      }
    }
    return ints;
  };
  auto numFragConns = fragConns.count();
  auto rConns = bitsToInts(reactionConns);
  auto perms = details::permMFromN(numFragConns, reactionConns.count());

  for (const auto &perm : perms) {
    connPerms.emplace_back();
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
  }

  return connPerms;
}

// Take the molFrags and flag those synthons that have pattern fingerprints
// where all the bits match with the fragment.  The pattern fingerprints are
// insensitive to isotope numbers, so this can be done on the initial
// fragmentation, without generating the combinations of connector numbers.
// Matches the pattFPs with the synthon sets in the order synthonOrder, but
// returns the bitsets in the original order.
std::vector<boost::dynamic_bitset<>> screenSynthonsWithFPs(
    const std::vector<std::unique_ptr<ExplicitBitVect>> &pattFPs,
    std::unique_ptr<SynthonSet> &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  std::vector<boost::dynamic_bitset<>> passedFPs;
  for (const auto &synthonSet : reaction->getSynthons()) {
    passedFPs.emplace_back(synthonSet.size());
  }

  boost::dynamic_bitset<> fragsMatched(synthonOrder.size());
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    const auto &synthonSet = reaction->getSynthons()[synthonOrder[i]];
    for (size_t j = 0; j < synthonSet.size(); ++j) {
      auto &synthon = synthonSet[j];
      if (AllProbeBitsMatch(*pattFPs[i], *synthon->getPattFP())) {
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
std::vector<boost::dynamic_bitset<>> getHitSynthons(
    std::vector<std::unique_ptr<RWMol>> &molFrags,
    const std::vector<boost::dynamic_bitset<>> &passedScreens,
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<unsigned int> &synthonOrder) {
  RDKit::MatchVectType dontCare;
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  for (const auto &synthonSet : reaction->getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }

  // The tests must be applied for all permutations of synthon list against
  // fragment.
  auto synthonOrders =
      details::permMFromN(molFrags.size(), reaction->getSynthons().size());

  boost::dynamic_bitset<> fragsMatched(synthonOrder.size());
  // Match the fragment to the synthon set in this order.
  for (size_t i = 0; i < synthonOrder.size(); ++i) {
    const auto &synthonsSet = reaction->getSynthons()[synthonOrder[i]];
    const auto &passedScreensSet = passedScreens[synthonOrder[i]];
    for (size_t j = 0; j < synthonsSet.size(); ++j) {
      if (passedScreensSet[j]) {
        auto &synthon = synthonsSet[j];
        if (SubstructMatch(*synthon->getMol(), *molFrags[i], dontCare)) {
          synthonsToUse[synthonOrder[i]][j] = true;
          fragsMatched[i] = true;
        }
      }
    }
    // if the fragment didn't match anything, the whole thing's a bust.
    if (!fragsMatched[i]) {
      synthonsToUse.clear();
      return synthonsToUse;
    }
  }
  // If all bits in one of the bitsets is unset, it means that nothing matched
  // that synthon.  If at least one of the bitsets has a set bit, all products
  // incorporating the synthon with no bits set must match the query so
  // should be used because the query matches products that don't incorporate
  // anything from 1 of the synthon lists.  For example, if the synthons are
  // [1*]Nc1c([2*])cccc1 and [1*]=CC=C[2*] and the query is c1ccccc1.
  bool someSet = std::any_of(
      synthonsToUse.begin(), synthonsToUse.end(),
      [](const boost::dynamic_bitset<> &bs) -> bool { return bs.any(); });
  if (someSet) {
    for (auto &rtu : synthonsToUse) {
      if (!rtu.count()) {
        rtu.set();
      }
    }
  }
  return synthonsToUse;
}

void buildConnectorRegions(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &connRegs,
    std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> &connRegFPs) {
  for (const auto &frag : molFrags) {
    auto fragConnRegs = getConnRegion(*frag);
    if (!fragConnRegs) {
      // There were no connector atoms.
      continue;
    }
    std::vector<std::unique_ptr<ROMol>> splitConnRegs;
    MolOps::getMolFrags(*fragConnRegs, splitConnRegs, false);
    connRegFPs.emplace_back();
    for (auto &cr : splitConnRegs) {
      connRegFPs.back().emplace_back(PatternFingerprintMol(*cr));
    }
    connRegs.push_back(std::move(splitConnRegs));
  }
}

// Return true if all the fragments have a connector region that matches
// something in the reaction, false otherwise.
bool checkConnectorRegions(
    const std::unique_ptr<SynthonSet> &reaction,
    const std::vector<std::vector<std::unique_ptr<ROMol>>> &connRegs,
    const std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>>
        &connRegFPs) {
  const auto &rxnConnRegs = reaction->getConnectorRegions();
  const auto &rxnConnRegsFP = reaction->getConnRegFP();
  MatchVectType dontCare;
  for (size_t i = 0; i < connRegFPs.size(); ++i) {
    bool connRegFound = false;
    for (size_t j = 0; j < connRegFPs[i].size(); ++j) {
      if (AllProbeBitsMatch(*connRegFPs[i][j], *rxnConnRegsFP)) {
        for (const auto &rxncr : rxnConnRegs) {
          if (SubstructMatch(*rxncr, *connRegs[i][j], dontCare)) {
            connRegFound = true;
            break;
          }
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
}  // namespace

std::vector<SynthonSpaceHitSet> SynthonSpace::searchFragSet(
    std::vector<std::unique_ptr<ROMol>> &fragSet) {
  std::vector<SynthonSpaceHitSet> results;

  auto pattFPs = makePatternFPs(fragSet);
  std::vector<std::vector<std::unique_ptr<ROMol>>> connRegs;
  std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> connRegFPs;
  std::vector<int> numFragConns;
  numFragConns.reserve(fragSet.size());
  for (const auto &frag : fragSet) {
    numFragConns.push_back(details::countConnections(MolToSmiles(*frag)));
  }

  auto conns = getConnectorPattern(fragSet);
  for (auto &it : d_reactions) {
    auto &reaction = it.second;
    // It can't be a hit if the number of fragments is more than the number
    // of synthon sets because some of the molecule won't be matched in any
    // of the potential products.  It can be less, in which case the unused
    // synthon set will be used completely, possibly resulting in a large
    // number of hits.
    if (fragSet.size() > reaction->getSynthons().size()) {
      continue;
    }

    // Check that all the frags have a connector region that matches something
    // in this reaction set.  Skip if not.
    if (connRegs.empty()) {
      buildConnectorRegions(fragSet, connRegs, connRegFPs);
    }
    if (!checkConnectorRegions(reaction, connRegs, connRegFPs)) {
      continue;
    }

    // Select only the synthons that have fingerprints that are a superset
    // of the fragment fingerprints.
    // Need to try all combinations of synthon orders.
    auto synthonOrders =
        details::permMFromN(pattFPs.size(), reaction->getSynthons().size());
    for (const auto &so : synthonOrders) {
      auto passedScreens = screenSynthonsWithFPs(pattFPs, reaction, so);
      // If none of the synthons passed the screens, move right along, nothing
      // to see.
      bool skip = std::all_of(
          passedScreens.begin(), passedScreens.end(),
          [](const boost::dynamic_bitset<> &bs) -> bool { return bs.none(); });
      if (skip) {
        continue;
      }

      // Get all the possible permutations of connector numbers compatible with
      // the number of synthon sets in this reaction.  So if the
      // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
      // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
      // that might be how they're labelled in the reaction database.
      auto connCombs =
          getConnectorPermutations(fragSet, conns, reaction->getConnectors());

      // Find all synthons that match the fragments with each connector
      // combination.
      for (auto &connComb : connCombs) {
        auto theseSynthons =
            getHitSynthons(connComb, passedScreens, reaction, so);
        if (!theseSynthons.empty()) {
          size_t numHits = std::accumulate(
              theseSynthons.begin(), theseSynthons.end(), 1,
              [](int prevRes, const boost::dynamic_bitset<> &s2) {
                return prevRes * s2.count();
              });
          if (numHits) {
            results.push_back(
                SynthonSpaceHitSet{reaction->getId(), theseSynthons, numHits});
          }
        }
      }
    }
  }
  return results;
}

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

namespace {
// class to step through all combinations of list of different sizes.
// returns (0,0,0), (0,0,1), (0,1,0) etc.
struct Stepper {
  explicit Stepper(std::vector<size_t> &sizes) : d_sizes(sizes) {
    d_currState = std::vector<size_t>(sizes.size(), 0);
  }
  void step() {
    // Don't do anything if we're at the end, but expect an infinite
    // loop if the user isn't wise to this.
    if (d_currState[0] == d_sizes[0]) {
      return;
    }
    std::int64_t i = static_cast<std::int64_t>(d_currState.size()) - 1;
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
  std::vector<size_t> d_currState;
  std::vector<size_t> d_sizes;
};
}  // namespace

void SynthonSpace::buildHits(
    const std::vector<SynthonSpaceHitSet> &hitsets, const ROMol &query,
    const SynthonSpaceSearchParams &params, const size_t totHits,
    std::set<std::string> &resultsNames,
    std::vector<std::unique_ptr<ROMol>> &results) const {
  if (hitsets.empty()) {
    return;
  }

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  const double randDiscrim =
      static_cast<double>(params.maxHits) / static_cast<double>(totHits);
  MatchVectType dontCare;

  for (const auto &hitset : hitsets) {
    const auto &synthonsToUse = hitset.synthonsToUse;
    auto synthons = getSynthonsToUse(synthonsToUse, hitset.reactionId);
    if (synthons.empty()) {
      return;
    }

    std::vector<size_t> numSynthons;
    numSynthons.reserve(synthons.size());
    for (auto &s : synthons) {
      numSynthons.push_back(s.size());
    }
    Stepper stepper(numSynthons);
    const size_t numReactions = synthonsToUse.size();
    MolzipParams mzparams;
    mzparams.label = MolzipLabel::Isotope;
    while (stepper.d_currState[0] != numSynthons[0]) {
      std::string combName =
          hitset.reactionId + "_" +
          synthons[0][stepper.d_currState[0]]->getProp<std::string>(
              common_properties::_Name);
      for (size_t i = 1; i < numReactions; ++i) {
        combName +=
            "_" + synthons[i][stepper.d_currState[i]]->getProp<std::string>(
                      common_properties::_Name);
      }
      if (!params.randomSample || params.maxHits == -1 ||
          (params.randomSample && dist(*d_randGen) < randDiscrim)) {
        if (resultsNames.insert(combName).second) {
          if (resultsNames.size() < static_cast<size_t>(params.hitStart)) {
            continue;
          }
          auto combMol = std::make_unique<ROMol>(
              ROMol(*synthons[0][stepper.d_currState[0]]));
          for (size_t i = 1; i < numReactions; ++i) {
            combMol.reset(
                combineMols(*combMol, *synthons[i][stepper.d_currState[i]]));
          }
          auto prod = molzip(*combMol, mzparams);
          MolOps::sanitizeMol(*dynamic_cast<RWMol *>(prod.get()));
          // Do a final check of the whole thing.  It can happen that the
          // fragments match synthons but the final product doesn't match, and
          // a key example is when the 2 synthons come together to form an
          // aromatic ring.  An aliphatic query can match the aliphatic synthon
          // so they are selected as a hit, but the final aromatic ring isn't
          // a match.  E.g. Cc1cccc(C(=O)N[1*])c1N=[2*] and c1ccoc1C(=[2*])[1*]
          // making Cc1cccc2c(=O)[nH]c(-c3ccco3)nc12.  The query c1ccc(CN)o1
          // when split is a match to the synthons (c1ccc(C[1*])o1 and [1*]N)
          // but the product the hydroxyquinazoline is aromatic, at least in
          // the RDKit model so the N in the query doesn't match.
          if (!SubstructMatch(*prod, query, dontCare)) {
            continue;
          }
          prod->setProp<std::string>(common_properties::_Name, combName);
          results.push_back(std::move(prod));
        }
      }
      if (results.size() == static_cast<size_t>(params.maxHits)) {
        return;
      }
      stepper.step();
    }
  }
}

std::vector<std::vector<ROMol *>> SynthonSpace::getSynthonsToUse(
    const std::vector<boost::dynamic_bitset<>> &synthonsToUse,
    const std::string &reaction_id) const {
  if (const auto &it = d_reactions.find(reaction_id); it == d_reactions.end()) {
    throw std::runtime_error("Reaction " + reaction_id +
                             "not in the reaction set.");
  }
  const auto &reaction = d_reactions.find(reaction_id)->second;

  std::vector<std::vector<ROMol *>> synthons(reaction->getSynthons().size(),
                                             std::vector<ROMol *>());
  for (size_t i = 0; i < synthonsToUse.size(); ++i) {
    for (size_t j = 0; j < synthonsToUse[i].size(); ++j) {
      if (synthonsToUse[i][j]) {
        synthons[i].push_back(reaction->getSynthons()[i][j]->getMol().get());
      }
    }
  }
  return synthons;
}
}  // namespace RDKit::SynthonSpaceSearch
