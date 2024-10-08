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
#include <GraphMol/HyperspaceSearch/Hyperspace.h>
#include <GraphMol/HyperspaceSearch/HyperspaceSubstructureSearch.h>
#include <GraphMol/HyperspaceSearch/ReactionSet.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace RDKit::HyperspaceSSSearch {
Hyperspace::Hyperspace(const std::string &fileName) : d_fileName(fileName) {
  readFile();
}

std::vector<std::unique_ptr<ROMol>> Hyperspace::search(
    const RDKit::ROMol &query, unsigned int maxBondSplits, int maxHits) {
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
  std::vector<HyperspaceHitSet> allHits;
  size_t totHits = 0;
  for (auto &fragSet : fragments) {
    auto theseHits = searchFragSet(fragSet);
    if (!theseHits.empty()) {
      totHits += std::reduce(
          theseHits.begin(), theseHits.end(), 0,
          [&](const size_t prevVal, const HyperspaceHitSet &hs) -> size_t {
            return prevVal + hs.numHits;
          });
      std::cout << "Total hits : " << totHits << std::endl;
      allHits.insert(allHits.end(), theseHits.begin(), theseHits.end());
    }
    // Just for safety, do a final check
    //    for (auto &hitMol : theseResults) {
    //      auto molName =
    //      hitMol->getProp<std::string>(common_properties::_Name); if
    //      (resultsNames.find(molName) == resultsNames.end() &&
    //          SubstructMatch(*hitMol, query, dontCare)) {
    //        results.emplace_back(std::unique_ptr<ROMol>(hitMol.release()));
    //        resultsNames.insert(molName);
    //      }
    //    }
    //    if (runningMaxHits != -1) {
    //      runningMaxHits -= theseResults.size();
    //      if (runningMaxHits < 1) {
    //        break;
    //      }
    //    }
  }
  std::sort(
      allHits.begin(), allHits.end(),
      [&](const HyperspaceHitSet &hs1, const HyperspaceHitSet &hs2) -> bool {
        if (hs1.reactionId == hs2.reactionId) {
          return hs1.numHits < hs2.numHits;
        } else {
          return hs1.reactionId < hs2.reactionId;
        }
      });

  buildHits(allHits, results, maxHits);
  return results;
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
  const static std::vector<std::regex> connRegexs{
      std::regex(R"(\[U\])"), std::regex(R"(\[Np\])"), std::regex(R"(\[Pu\])"),
      std::regex(R"(\[Am\])")};
  const static std::vector<std::string> repls{"[1*]", "[2*]", "[3*]", "[4*]"};
  for (size_t i = 0; i < connRegexs.size(); ++i) {
    smiles = std::regex_replace(smiles, connRegexs[i], repls[i]);
  }
}
}  // namespace

void Hyperspace::readFile() {
  std::ifstream ifs(d_fileName);
  // the first line is headers, check they're in the right format,
  // which is tab-separated:
  // SMILES	synton_id	synton#	reaction_id
  // or tab-separated:
  // SMILES	synton_id	synton#	reaction_id release
  // or
  // SMILES,synton_id,synton_role,reaction_id
  // Note that we keep the spelling "synton" from the original paper
  // and example input file, and allow any whitespace rather than just tab.
  std::string firstLine;
  getline(ifs, firstLine);
  //  std::cout << "parsing : " << firstLine << std::endl;
  std::regex regexz("[\\s,]+");

  auto lineParts = splitLine(firstLine, regexz);
  std::vector<std::vector<std::string>> firstLineOpts{
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id"},
      std::vector<std::string>{"SMILES", "synton_id", "synton#", "reaction_id",
                               "release"},
      std::vector<std::string>{"SMILES", "synton_id", "synton_role",
                               "reaction_id"}};
  int format = -1;
  for (size_t i = 0; i < firstLineOpts.size(); ++i) {
    if (lineParts == firstLineOpts[i]) {
      format = i;
    }
  }
  if (format == -1) {
    throw std::runtime_error("Bad format for hyperspace file " + d_fileName);
  }
  std::cout << "Format = " << format << std::endl;
  std::string nextLine;
  int lineNum = 1;

  while (getline(ifs, nextLine)) {
    ++lineNum;
    if (nextLine.empty() || nextLine[0] == '#') {
      continue;
    }
    auto nextReag = splitLine(nextLine, regexz);
    if (nextReag.size() < 4) {
      throw std::runtime_error("Bad format for hyperspace file " + d_fileName +
                               " on line " + std::to_string(lineNum));
    }
    //    std::cout << nextReag[0] << " : " << nextReag[1] << std::endl;
    if (auto it = d_reactions.find(nextReag[3]); it == d_reactions.end()) {
      d_reactions.insert(std::make_pair(
          nextReag[3], std::make_unique<ReactionSet>(nextReag[3])));
    }
    fixConnectors(nextReag[0]);
    auto &currReaction = d_reactions[nextReag[3]];
    size_t synthonNum{std::numeric_limits<size_t>::max()};
    if (format == 0 || format == 1) {
      synthonNum = std::stoi(nextReag[2]);
    } else if (format == 2) {
      // in this case it's a string "synton_2" etc.
      synthonNum = std::stoi(nextReag[2].substr(7));
    }
    currReaction->addReagent(synthonNum, nextReag[0], nextReag[1]);
  }
  for (auto &it : d_reactions) {
    auto &reaction = it.second;
    reaction->assignConnectorsUsed();
  }
}

namespace {
// Return a bitset giving the different connector types in this
// molecule.
boost::dynamic_bitset<> getConnectorPattern(
    const std::vector<std::unique_ptr<ROMol>> &fragSet) {
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

// Make pattern fps for the fragments and re-order fps and fragments into
// descending order of number of bits set in the fp (the "largest fragment
// heuristic").
std::vector<std::unique_ptr<ExplicitBitVect>> makePatternFPs(
    std::vector<std::unique_ptr<ROMol>> &molFrags) {
  std::vector<std::unique_ptr<ExplicitBitVect>> pattFPs;
  for (const auto &frag : molFrags) {
    pattFPs.emplace_back(
        std::unique_ptr<ExplicitBitVect>(PatternFingerprintMol(*frag, 2048)));
  }
  // Sort by descending number of bits set.  This is the "largest fragment
  // heuristic" from the paper.  The FP with the largest number of bits
  // is the most likely to screen out a matching reagent set since smaller,
  // less complex fragments are more likely to match something, so screen
  // with that first.
  std::vector<std::pair<size_t, ExplicitBitVect *>> fps(pattFPs.size());
  for (size_t i = 0; i < pattFPs.size(); ++i) {
    fps[i] = std::make_pair(i, pattFPs[i].get());
  }
  std::sort(fps.begin(), fps.end(),
            [&](const std::pair<size_t, ExplicitBitVect *> &fp1,
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
    connPerms.push_back(std::vector<std::unique_ptr<RWMol>>());
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

// Take the molFrags and flag those reagents that have pattern fingerprints
// where all the bits match with the fragment.  The pattern fingerprints are
// insensitive to isotope numbers, so this can be done on the initial
// fragmentation, without generating the combinations of connector numbers.
// Matches the pattFPs with the reagent sets in the order reagentOrder, but
// returns the bitsets in the original order.
std::vector<boost::dynamic_bitset<>> screenReagentsWithFPs(
    const std::vector<std::unique_ptr<ExplicitBitVect>> &pattFPs,
    std::unique_ptr<ReactionSet> &reaction,
    std::vector<unsigned int> reagentOrder) {
  std::vector<boost::dynamic_bitset<>> reagsToUse;
  //  std::cout << "screenReagentsWithFPs" << std::endl;
  std::vector<boost::dynamic_bitset<>> passedFPs;
  for (const auto &reagSet : reaction->reagents()) {
    passedFPs.push_back(boost::dynamic_bitset<>(reagSet.size()));
  }

  std::vector<boost::dynamic_bitset<>> thisPass;
  for (const auto &reagSet : reaction->reagents()) {
    thisPass.push_back(boost::dynamic_bitset<>(reagSet.size()));
  }
  boost::dynamic_bitset<> fragsMatched(reagentOrder.size());
  for (size_t i = 0; i < reagentOrder.size(); ++i) {
    const auto &reagSet = reaction->reagents()[reagentOrder[i]];
    for (size_t j = 0; j < reagSet.size(); ++j) {
      auto &reag = reagSet[j];
      if (AllProbeBitsMatch(*pattFPs[i], *reag->pattFP())) {
        thisPass[reagentOrder[i]][j] = true;
        fragsMatched[i] = true;
      }
    }
    // If nothing matched this fragment, the whole thing's a bust.
    if (!fragsMatched[i]) {
      //      std::cout << "Nothing matches frag " << i << " : "
      //                << pattFPs[i]->getNumOnBits() << std::endl;
      break;
    }
  }
  // If all the fragments had a match, these results are valid.
  if (fragsMatched.count() == fragsMatched.size()) {
    for (size_t i = 0; i < passedFPs.size(); ++i) {
      passedFPs[i] |= thisPass[i];
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
    std::vector<std::unique_ptr<RWMol>> &molFrags,
    const std::vector<boost::dynamic_bitset<>> &passedScreens,
    std::unique_ptr<ReactionSet> &reaction,
    const std::vector<unsigned int> &reagentOrder) {
  RDKit::MatchVectType dontCare;
  std::vector<boost::dynamic_bitset<>> reagsToUse;
  for (const auto &reagSet : reaction->reagents()) {
    reagsToUse.push_back(boost::dynamic_bitset<>(reagSet.size()));
  }

  // The tests must be applied for all permutations of reagent list against
  // fragment.
  auto reagentOrders =
      details::permMFromN(molFrags.size(), reaction->reagents().size());

  boost::dynamic_bitset<> fragsMatched(reagentOrder.size());
  // Match the fragment to the reagent set in this order.
  for (size_t i = 0; i < reagentOrder.size(); ++i) {
    const auto &reagSet = reaction->reagents()[reagentOrder[i]];
    const auto &passedScreensSet = passedScreens[reagentOrder[i]];
    for (size_t j = 0; j < reagSet.size(); ++j) {
      if (passedScreensSet[j]) {
        auto &reag = reagSet[j];
        if (SubstructMatch(*reag->mol(), *molFrags[i], dontCare)) {
          reagsToUse[reagentOrder[i]][j] = true;
          fragsMatched[i] = true;
        }
      }
    }
    // if the fragment didn't match anything, the whole thing's a bust.
    if (!fragsMatched[i]) {
      reagsToUse.clear();
      return reagsToUse;
    }
  }
  // If all bits in one of the bitsets is unset, it means that nothing matched
  // that reagent.  If at least one of the bitsets has a set bit, all products
  // incorporating the reagent with no bits set must match the query so
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
    connRegFPs.push_back(std::vector<std::unique_ptr<ExplicitBitVect>>());
    for (auto &cr : splitConnRegs) {
      connRegFPs.back().emplace_back(PatternFingerprintMol(*cr));
    }
    connRegs.push_back(std::move(splitConnRegs));
  }
}

// Return true if all the fragments have a connector region that matches
// something in the reaction, false otherwise.
bool checkConnectorRegions(
    std::unique_ptr<ReactionSet> &reaction,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &connRegs,
    std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> &connRegFPs) {
  const auto &rxnConnRegs = reaction->connectorRegions();
  const auto &rxnConnRegsFP = reaction->connRegFP();
  RDKit::MatchVectType dontCare;
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

std::vector<HyperspaceHitSet> Hyperspace::searchFragSet(
    std::vector<std::unique_ptr<ROMol>> &fragSet, int maxHits) {
  std::vector<HyperspaceHitSet> results;

  auto pattFPs = makePatternFPs(fragSet);
  std::vector<std::vector<std::unique_ptr<ROMol>>> connRegs;
  std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> connRegFPs;

  //  std::cout << "searchFragSet" << std::endl;
  auto conns = getConnectorPattern(fragSet);
  for (auto &it : d_reactions) {
    auto &reaction = it.second;
    //    std::cout << "Searching for " << fragSet.size() << " ::: ";
    //    for (size_t i = 0; i < fragSet.size(); ++i) {
    //      const auto &f = fragSet[i];
    //      std::cout << f->getNumAtoms() << " : " << MolToSmiles(*f) << " : "
    //                << MolToSmarts(*f) << " :: ";
    //      std::cout << f->getNumAtoms() << " : " << MolToSmiles(*f) << " : "
    //                << pattFPs[i]->getNumOnBits() << " :: ";
    //    }
    //    std::cout << " in " << reaction->id() << " : "
    //              << reaction->reagents().size() << std::endl;
    // It can't be a hit if the number of fragments is more than the number
    // of reagent sets because some of the molecule won't be matched in any
    // of the potential products.  It can be less, in which case the unused
    // reagent set will be used completely, possibly resulting in a large
    // number of hits.
    if (fragSet.size() > reaction->reagents().size()) {
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

    // Select only the reagents that have fingerprints that are a superset
    // of the fragment fingerprints.
    // Need to try all combinations of reagent orders.
    auto reagentOrders =
        details::permMFromN(pattFPs.size(), reaction->reagents().size());
    for (const auto &ro : reagentOrders) {
      auto passedScreens = screenReagentsWithFPs(pattFPs, reaction, ro);
      // If none of the reagents passed the screens, move right along, nothing
      // to see.
      bool skip = true;
      for (const auto &ps : passedScreens) {
        if (ps.count()) {
          skip = false;
          break;
        }
      }
      if (skip) {
        continue;
      }

      // Get all the possible permutations of connector numbers compatible with
      // the number of reagent sets in this reaction.  So if the
      // fragmented molecule is C[1*].N[2*] and there are 3 reagent sets
      // we also try C[2*].N[1*]  and C[3*].N[2*] because
      // that might be how they're labelled in the reaction database.
      auto connCombs =
          getConnectorPermutations(fragSet, conns, reaction->connectors());

      // Find all reagents that match the fragments with each connector
      // combination.
      for (auto &connComb : connCombs) {
        auto theseReagents =
            getHitReagents(connComb, passedScreens, reaction, ro);
        if (!theseReagents.empty()) {
          size_t numHits =
              std::reduce(theseReagents.begin(), theseReagents.end(), 1,
                          [&](int prevRes, const boost::dynamic_bitset<> &s2) {
                            return prevRes * s2.count();
                          });
          if (numHits) {
            for (size_t i = 0; i < theseReagents.size(); ++i) {
              std::cout << reaction->id() << " reagents " << i << " : "
                        << theseReagents[i].count() << " : " << theseReagents[i]
                        << std::endl;
            }
            results.push_back(
                HyperspaceHitSet{reaction->id(), theseReagents, numHits});
          }
        }
      }
    }
  }
  //  std::cout << "Leaving searchFragSet" << std::endl;
  return results;
}

void Hyperspace::writeToDBStream(const std::string &outFile) const {
  std::ofstream os(outFile, std::fstream::binary | std::fstream::trunc);
  streamWrite(os, d_reactions.size());
  for (const auto &rs : d_reactions) {
    rs.second->writeToDBStream(os);
  }
  os.close();
}

void Hyperspace::readFromDBStream(const std::string &inFile) {
  try {
    std::ifstream is(inFile, std::fstream::binary);
    size_t numRS;
    streamRead(is, numRS);
    for (size_t i = 0; i < numRS; ++i) {
      ReactionSet *rs = new ReactionSet;
      rs->readFromDBStream(is);
      d_reactions.insert(make_pair(rs->id(), rs));
    }
  } catch (std::exception &e) {
    std::cerr << "Error : " << e.what() << " for file " << inFile << "\n";
    exit(1);
  }
}

void Hyperspace::summarise(std::ostream &os) const {
  os << "Read from file " << d_fileName << "\n"
     << "Number of reactions : " << d_reactions.size() << "\n";
  for (const auto &reaction : d_reactions) {
    const auto &rxn = reaction.second;
    os << "Reaction name " << rxn->id() << "\n";
    for (size_t i = 0; i < rxn->reagents().size(); ++i) {
      os << "  Synthon set " << i << " has " << rxn->reagents()[i].size()
         << " reagents" << "\n";
    }
  }
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

void Hyperspace::buildHits(const std::vector<HyperspaceHitSet> &hitsets,
                           std::vector<std::unique_ptr<ROMol>> &results,
                           int maxHits) {
  if (hitsets.empty()) {
    return;
  }

  // Keep track of the result names so we can weed out duplicates by
  // reaction and reagents.  Different splits may give rise to the same
  // reagent combination.  This will keep the same molecule produced via
  // different reactions which I think makes sense.
  std::set<std::string> resultsNames;

  for (const auto &hitset : hitsets) {
    std::cout << "Build hits with " << hitset.reactionId << std::endl;
    const auto &reagentsToUse = hitset.reagsToUse;
    auto reags = getReagentsToUse(reagentsToUse, hitset.reactionId);
    if (reags.empty()) {
      //    std::cout << "Nothing to do" << std::endl;
      return;
    }

    std::vector<int> numReags;
    for (auto &r : reags) {
      std::cout << "Number of reagents : " << r.size() << std::endl;
      numReags.push_back(r.size());
    }
    Stepper stepper(numReags);
    const int numReactions = reagentsToUse.size();
    MolzipParams params;
    params.label = MolzipLabel::Isotope;
    while (stepper.d_currState[0] != numReags[0]) {
      std::string combName =
          hitset.reactionId + "_" +
          reags[0][stepper.d_currState[0]]->getProp<std::string>(
              common_properties::_Name);
      for (int i = 1; i < numReactions; ++i) {
        combName +=
            "_" + reags[i][stepper.d_currState[i]]->getProp<std::string>(
                      common_properties::_Name);
      }
      if (resultsNames.insert(combName).second) {
        std::unique_ptr<ROMol> combMol(
            new ROMol(*reags[0][stepper.d_currState[0]]));
        for (int i = 1; i < numReactions; ++i) {
          combMol.reset(
              combineMols(*combMol, *reags[i][stepper.d_currState[i]]));
        }
        auto prod = molzip(*combMol, params);
        prod->setProp<std::string>(common_properties::_Name, combName);
        MolOps::sanitizeMol(*static_cast<RWMol *>(prod.get()));
        results.push_back(std::move(prod));
      }
      // -1 means no limit, and size_t is unsigned, so the cast will make for
      // a very large number.
      if (results.size() == static_cast<size_t>(maxHits)) {
        return;
      }
      stepper.step();
    }
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

  std::vector<std::vector<ROMol *>> reags(reaction->reagents().size(),
                                          std::vector<ROMol *>());
  for (size_t i = 0; i < reagentsToUse.size(); ++i) {
    for (size_t j = 0; j < reagentsToUse[i].size(); ++j) {
      if (reagentsToUse[i][j]) {
        reags[i].push_back(reaction->reagents()[i][j]->mol().get());
      }
    }
  }
  return reags;
}
}  // namespace RDKit::HyperspaceSSSearch
