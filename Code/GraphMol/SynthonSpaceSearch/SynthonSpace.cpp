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
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace RDKit::SynthonSpaceSearch {

long SynthonSpace::numProducts() const {
  long totSize = 0;
  for (const auto &reaction : d_reactions) {
    const auto &rxn = reaction.second;
    size_t thisSize = 1;
    for (size_t i = 0; i < rxn->reagents().size(); ++i) {
      thisSize *= rxn->reagents()[i].size();
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
      d_randGen.reset(new std::mt19937(rd()));
    }
    if (params.randomSeed != -1) {
      d_randGen->seed(params.randomSeed);
    }
  }
  std::vector<std::unique_ptr<ROMol>> results;
  RDKit::MatchVectType dontCare;

  auto fragments = details::splitMolecule(query, params.maxBondSplits);
  //  for (const auto &frags : fragments) {
  //    std::cout << "****************************" << std::endl;
  //    for (const auto &f : frags) {
  //      std::cout << MolToSmiles(*f) << std::endl;
  //    }
  //  }
  std::vector<SynthonSpaceHitSet> allHits;
  size_t totHits = 0;
  for (auto &fragSet : fragments) {
    auto theseHits = searchFragSet(fragSet);
    if (!theseHits.empty()) {
      totHits += std::accumulate(
          theseHits.begin(), theseHits.end(), 0,
          [&](const size_t prevVal, const SynthonSpaceHitSet &hs) -> size_t {
            return prevVal + hs.numHits;
          });
      allHits.insert(allHits.end(), theseHits.begin(), theseHits.end());
    }
  }
  std::sort(allHits.begin(), allHits.end(),
            [&](const SynthonSpaceHitSet &hs1,
                const SynthonSpaceHitSet &hs2) -> bool {
              if (hs1.reactionId == hs2.reactionId) {
                return hs1.numHits < hs2.numHits;
              } else {
                return hs1.reactionId < hs2.reactionId;
              }
            });

  if (params.buildHits) {
    // Keep track of the result names so we can weed out duplicates by
    // reaction and reagents.  Different splits may give rise to the same
    // reagent combination.  This will keep the same molecule produced via
    // different reactions which I think makes sense.  The resultsNames will
    // be accumulated even if the molecule itself doesn't make it into the
    // results set, for example if it isn't a random selection or it's
    // outside maxHits or hitStart.
    std::set<std::string> resultsNames;
    // The random sampling, if requested, doesn't always produce the required
    // number of hits in 1 sweep.
    size_t tmpMaxHits(params.maxHits);
    int numSweeps = 0;
    while (results.size() < std::min(tmpMaxHits, totHits) && numSweeps < 10) {
      buildHits(allHits, query, params, totHits, resultsNames, results);
      // totHits is an upper bound, so may not be reached.
      if (params.maxHits == -1 || !params.randomSample) {
        break;
      }
      //      std::cout << "Sweep " << numSweeps << " : " << params.maxHits << "
      //      : "
      //                << totHits << " : " << resultsNames.size() << " : "
      //                << results.size() << std::endl;
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
  const static std::vector<std::regex> connRegexs{
      std::regex(R"(\[U\])"), std::regex(R"(\[Np\])"), std::regex(R"(\[Pu\])"),
      std::regex(R"(\[Am\])")};
  const static std::vector<std::string> repls{"[1*]", "[2*]", "[3*]", "[4*]"};
  for (size_t i = 0; i < connRegexs.size(); ++i) {
    smiles = std::regex_replace(smiles, connRegexs[i], repls[i]);
  }
}
}  // namespace

namespace {
// Return a bitset giving the different connector types in this
// molecule.
boost::dynamic_bitset<> getConnectorPattern(
    const std::vector<std::unique_ptr<ROMol>> &fragSet) {
  boost::dynamic_bitset<> conns(4);
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
          //          std::cout << reag->smiles() << " matched to "
          //                    << MolToSmiles(*molFrags[i], true, false, -1,
          //                    false)
          //                    << std::endl;
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

std::vector<SynthonSpaceHitSet> SynthonSpace::searchFragSet(
    std::vector<std::unique_ptr<ROMol>> &fragSet) {
  std::vector<SynthonSpaceHitSet> results;

  auto pattFPs = makePatternFPs(fragSet);
  std::vector<std::vector<std::unique_ptr<ROMol>>> connRegs;
  std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> connRegFPs;
  std::vector<int> numFragConns;
  for (const auto &frag : fragSet) {
    numFragConns.push_back(details::countConnections(MolToSmiles(*frag)));
  }

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
      //      if (!checkNumberOfConnections(numFragConns, reaction, ro)) {
      //        continue;
      //      }
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
          size_t numHits = std::accumulate(
              theseReagents.begin(), theseReagents.end(), 1,
              [&](int prevRes, const boost::dynamic_bitset<> &s2) {
                return prevRes * s2.count();
              });
          if (numHits) {
            //            for (size_t i = 0; i < theseReagents.size(); ++i) {
            //              std::cout << reaction->id() << " reagents " << i <<
            //              " : "
            //                        << theseReagents[i].count() << " : " <<
            //                        theseReagents[i]
            //                        << std::endl;
            //            }
            results.push_back(
                SynthonSpaceHitSet{reaction->id(), theseReagents, numHits});
          }
        }
      }
    }
  }
  //  std::cout << "Leaving searchFragSet" << std::endl;
  return results;
}

void SynthonSpace::readTextFile(const std::string &inFile) {
  d_fileName = inFile;
  std::ifstream ifs(d_fileName);
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
    throw std::runtime_error("Bad format for SynthonSpace file " + d_fileName);
  }
  std::string nextLine;
  int lineNum = 1;

  while (getline(ifs, nextLine)) {
    ++lineNum;
    if (nextLine.empty() || nextLine[0] == '#') {
      continue;
    }
    auto nextReag = splitLine(nextLine, regexz);
    if (nextReag.size() < 4) {
      throw std::runtime_error("Bad format for SynthonSpace file " +
                               d_fileName + " on line " +
                               std::to_string(lineNum));
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

  // Do some final processing.
  for (auto &it : d_reactions) {
    auto &reaction = it.second;
    reaction->assignConnectorsUsed();
  }
}

void SynthonSpace::writeDBFile(const std::string &outFile) const {
  std::ofstream os(outFile, std::fstream::binary | std::fstream::trunc);
  streamWrite(os, d_reactions.size());
  for (const auto &rs : d_reactions) {
    rs.second->writeToDBStream(os);
  }
  os.close();
}

void SynthonSpace::readDBFile(const std::string &inFile) {
  d_fileName = inFile;
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

void SynthonSpace::summarise(std::ostream &os) const {
  os << "Read from file " << d_fileName << "\n"
     << "Number of reactions : " << d_reactions.size() << "\n";
  size_t totSize = 0;
  for (const auto &reaction : d_reactions) {
    const auto &rxn = reaction.second;
    os << "Reaction name " << rxn->id() << "\n";
    size_t thisSize = 1;
    for (size_t i = 0; i < rxn->reagents().size(); ++i) {
      os << "  Synthon set " << i << " has " << rxn->reagents()[i].size()
         << " reagents" << "\n";
      thisSize *= rxn->reagents()[i].size();
    }
    totSize += thisSize;
  }
  os << "Approximate number of molecules in SynthonSpace : " << totSize
     << std::endl;
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

void SynthonSpace::buildHits(const std::vector<SynthonSpaceHitSet> &hitsets,
                             const ROMol &query,
                             const SynthonSpaceSearchParams &params,
                             size_t totHits,
                             std::set<std::string> &resultsNames,
                             std::vector<std::unique_ptr<ROMol>> &results) {
  if (hitsets.empty()) {
    return;
  }

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  double randDiscrim = double(params.maxHits) / double(totHits);
  RDKit::MatchVectType dontCare;

  for (const auto &hitset : hitsets) {
    //    std::cout << "Build hits with " << hitset.reactionId << std::endl;
    const auto &reagentsToUse = hitset.reagsToUse;
    auto reags = getReagentsToUse(reagentsToUse, hitset.reactionId);
    if (reags.empty()) {
      //    std::cout << "Nothing to do" << std::endl;
      return;
    }

    std::vector<int> numReags;
    for (auto &r : reags) {
      //      std::cout << "Number of reagents : " << r.size() << std::endl;
      numReags.push_back(r.size());
    }
    Stepper stepper(numReags);
    const int numReactions = reagentsToUse.size();
    MolzipParams mzparams;
    mzparams.label = MolzipLabel::Isotope;
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
      if (!params.randomSample || params.maxHits == -1 ||
          (params.randomSample && dist(*d_randGen) < randDiscrim)) {
        if (resultsNames.insert(combName).second) {
          if (resultsNames.size() < static_cast<size_t>(params.hitStart)) {
            continue;
          }
          std::unique_ptr<ROMol> combMol(
              new ROMol(*reags[0][stepper.d_currState[0]]));
          for (int i = 1; i < numReactions; ++i) {
            combMol.reset(
                combineMols(*combMol, *reags[i][stepper.d_currState[i]]));
          }
          //        std::cout << "Zipping " << MolToSmiles(*combMol) <<
          //        std::endl;
          auto prod = molzip(*combMol, mzparams);
          MolOps::sanitizeMol(*static_cast<RWMol *>(prod.get()));
          // Do a final check of the whole thing.  It can happen that the
          // fragments match synthons but the final product doesn't match, and
          // a key example is when the 2 synthons come together to form an
          // aromatic ring.  An aliphatic query can match the aliphatic synthon
          // so they are selected as a hit, but the final aromatic ring isn't
          // a match.  E.g. Cc1cccc(C(=O)N[1*])c1N=[2*] and c1ccoc1C(=[2*])[1*]
          // making Cc1cccc2c(=O)[nH]c(-c3ccco3)nc12.  The query c1ccc(CN)o1
          // when split is a match to the reagents (c1ccc(C[1*])o1 and [1*]N)
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

std::vector<std::vector<ROMol *>> SynthonSpace::getReagentsToUse(
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
}  // namespace RDKit::SynthonSpaceSearch
