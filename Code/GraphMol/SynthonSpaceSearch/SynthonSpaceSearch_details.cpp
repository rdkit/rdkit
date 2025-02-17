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
#include <list>
#include <memory>
#include <regex>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <RDGeneral/ControlCHandler.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <boost/fusion/container/vector/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>

namespace RDKit::SynthonSpaceSearch::details {

bool checkTimeOut(const TimePoint *endTime) {
  if (endTime != nullptr && Clock::now() > *endTime) {
    BOOST_LOG(rdWarningLog) << "Timed out.\n";
    return true;
  }
  return false;
}

// get a vector of vectors of unsigned ints that are all combinations of
// M items chosen from N e.g. all combinations of 3 bonds from a
// molecule. A modified form of the code in the first answer from
// https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
std::vector<std::vector<unsigned int>> combMFromN(const unsigned int m,
                                                  const unsigned int n) {
  std::string allN(m, 1);
  allN.resize(n, 0);
  std::vector<std::vector<unsigned int>> combs;
  do {
    combs.emplace_back();
    for (unsigned int i = 0; i < n; ++i) {
      if (allN[i]) {
        combs.back().push_back(i);
      }
    }
  } while (std::prev_permutation(allN.begin(), allN.end()));
  return combs;
}

std::vector<std::vector<unsigned int>> permMFromN(const unsigned int m,
                                                  const unsigned int n) {
  std::vector<std::vector<unsigned int>> perms;
  auto combs = combMFromN(m, n);
  for (auto &c : combs) {
    do {
      perms.push_back(c);
    } while (std::next_permutation(c.begin(), c.end()));
  }

  return perms;
}

// The fragmentation is valid if the 2 ends of each bond are in different
// fragments.  This assumes there are no ring-closing reactions in the
// library, which is probably ok.
bool checkConnectorsInDifferentFrags(
    const std::vector<std::unique_ptr<ROMol>> &molFrags, const int numSplits) {
  // Loop over the isotope numbers of the ends of the splits
  for (const auto &frag : molFrags) {
    for (int j = 1; j <= numSplits; ++j) {
      int dummyCount = 0;
      for (const auto &atom : frag->atoms()) {
        if (!atom->getAtomicNum() &&
            atom->getIsotope() == static_cast<unsigned int>(j)) {
          if (dummyCount) {
            return false;
          }
          ++dummyCount;
        }
      }
    }
  }
  return true;
}

bool checkConnectorsInDifferentFrags(const ROMol &mol,
                                     const VECT_INT_VECT &fragIdxs,
                                     const int numSplits) {
  int dummyAtoms[2 * MAX_CONNECTOR_NUM + 2] = {0};
  for (const auto &atom : mol.atoms()) {
    if (!atom->getAtomicNum()) {
      if (auto dummy = atom->getIsotope(); dummy <= MAX_CONNECTOR_NUM) {
        int pos = 2 * dummy;
        if (dummyAtoms[pos]) {
          dummyAtoms[pos + 1] = atom->getIdx() + 1;
        } else {
          dummyAtoms[pos] = atom->getIdx() + 1;
        }
      }
    }
  }
  for (const auto &fragIdx : fragIdxs) {
    for (int j = 0; j < numSplits; ++j) {
      if (dummyAtoms[2 * j]) {
        const int d1 = dummyAtoms[2 * j] - 1;
        const int d2 = dummyAtoms[2 * j + 1] - 1;
        int dummyCount = 0;
        for (const auto fi : fragIdx) {
          if (fi == d1 || fi == d2) {
            if (dummyCount) {
              return false;
            }
            ++dummyCount;
          }
        }
      }
    }
  }
  return true;
}

// Traverse the bonds from aromBond and return all the ones that are aromatic.
std::vector<const Bond *> getContiguousAromaticBonds(const ROMol &mol,
                                                     const Bond *aromBond) {
  std::vector<const Bond *> aromBonds(1, aromBond);
  std::list<const Bond *> toDo(1, aromBond);
  boost::dynamic_bitset<> done(mol.getNumBonds());
  done[aromBond->getIdx()] = true;
  while (!toDo.empty()) {
    const auto nextBond = toDo.front();
    toDo.pop_front();
    for (const auto nbr :
         make_iterator_range(mol.getAtomNeighbors(nextBond->getBeginAtom()))) {
      if (auto bond = mol.getBondBetweenAtoms(nextBond->getBeginAtomIdx(), nbr);
          !done[bond->getIdx()] && bond->getIsAromatic()) {
        aromBonds.push_back(bond);
        done[bond->getIdx()] = true;
        toDo.push_back(bond);
      }
    }
    for (const auto nbr :
         make_iterator_range(mol.getAtomNeighbors(nextBond->getEndAtom()))) {
      if (auto bond = mol.getBondBetweenAtoms(nextBond->getEndAtomIdx(), nbr);
          !done[bond->getIdx()] && bond->getIsAromatic()) {
        aromBonds.push_back(bond);
        done[bond->getIdx()] = true;
        toDo.push_back(bond);
      }
    }
  }
  return aromBonds;
}

namespace {
boost::dynamic_bitset<> flagRingBonds(const ROMol &mol) {
  const auto ringInfo = mol.getRingInfo();
  if (!ringInfo->isInitialized()) {
    // Query molecules don't seem to have the ring info generated on creation.
    MolOps::findSSSR(mol);
  }
  boost::dynamic_bitset<> ringBonds(mol.getNumBonds());
  for (const auto &r : ringInfo->bondRings()) {
    for (const auto b : r) {
      ringBonds.set(b);
    }
  }
  return ringBonds;
}

void addBondsToList(const ROMol &mol, const Atom *atom,
                    boost::dynamic_bitset<> &ringBonds,
                    boost::dynamic_bitset<> &doneAtoms,
                    std::list<const Atom *> &atoms,
                    std::vector<const Bond *> &ringBlock) {
  for (const auto nbond : mol.atomBonds(atom)) {
    if (ringBonds[nbond->getIdx()]) {
      ringBonds.set(nbond->getIdx(), false);
      ringBlock.push_back(nbond);
      const Atom *otherAtom = nbond->getOtherAtom(atom);
      if (!doneAtoms[otherAtom->getIdx()]) {
        atoms.push_back(otherAtom);
        doneAtoms.set(otherAtom->getIdx(), true);
      }
    }
  }
}

// Get all the contiguous ring blocks, e.g. in c1ccccc1Oc1cccc2[nH]ccc12
// get the benzene and indole as separate pieces.  Pass ringBonds by
// value as it gets changed and we'll be needing it again later.
std::vector<std::vector<const Bond *>> getRingBlocks(
    const ROMol &mol, boost::dynamic_bitset<> ringBonds) {
  std::vector<std::vector<const Bond *>> ringBlocks;
  while (ringBonds.count()) {
    for (const auto bond : mol.bonds()) {
      if (ringBonds[bond->getIdx()]) {
        ringBlocks.emplace_back(std::vector<const Bond *>{bond});
        ringBonds.set(bond->getIdx(), false);
        boost::dynamic_bitset<> doneAtoms(mol.getNumAtoms());
        std::list<const Atom *> toDo;
        addBondsToList(mol, bond->getBeginAtom(), ringBonds, doneAtoms, toDo,
                       ringBlocks.back());
        addBondsToList(mol, bond->getEndAtom(), ringBonds, doneAtoms, toDo,
                       ringBlocks.back());
        while (!toDo.empty()) {
          auto nextAtom = toDo.front();
          toDo.pop_front();
          addBondsToList(mol, nextAtom, ringBonds, doneAtoms, toDo,
                         ringBlocks.back());
        }
        std::sort(ringBlocks.back().begin(), ringBlocks.back().end(),
                  [](const Bond *b1, const Bond *b2) -> bool {
                    return b1->getIdx() < b2->getIdx();
                  });
      }
    }
  }
  std::cout << "ringBlocks size: " << ringBlocks.size() << std::endl;
  for (const auto &rb : ringBlocks) {
    for (auto bond : rb) {
      std::cout << bond->getIdx() << " ";
    }
    std::cout << std::endl;
  }
  return ringBlocks;
}

bool bondPairFragmentsBlock(
    size_t bondi, size_t bondj, unsigned int numAtoms,
    const std::vector<const Bond *> &ringBlock,
    std::vector<boost::dynamic_bitset<>> &ringAdjTable) {
  const Bond *bi = ringBlock[bondi];
  const Bond *bj = ringBlock[bondj];

  // Temporarily break the 2 bonds
  ringAdjTable[bi->getBeginAtomIdx()][bi->getEndAtomIdx()] = false;
  ringAdjTable[bi->getEndAtomIdx()][bi->getBeginAtomIdx()] = false;
  ringAdjTable[bj->getBeginAtomIdx()][bj->getEndAtomIdx()] = false;
  ringAdjTable[bj->getEndAtomIdx()][bj->getBeginAtomIdx()] = false;

  std::list<size_t> atoms(1, ringBlock[bondi]->getBeginAtomIdx());
  std::list<size_t> toDo(1, ringBlock[bondi]->getBeginAtomIdx());
  boost::dynamic_bitset<> doneAtom(ringAdjTable.size());
  doneAtom[ringBlock[bondi]->getBeginAtomIdx()] = true;
  while (!toDo.empty()) {
    const auto nextAtom = toDo.front();
    toDo.pop_front();
    const auto &theseConns = ringAdjTable[nextAtom];
    for (size_t i = 0; i < ringAdjTable.size(); i++) {
      if (theseConns[i] && !doneAtom[i]) {
        doneAtom[i] = true;
        toDo.push_back(i);
        atoms.push_back(i);
      }
    }
  }

  ringAdjTable[bi->getBeginAtomIdx()][bi->getEndAtomIdx()] = true;
  ringAdjTable[bi->getEndAtomIdx()][bi->getBeginAtomIdx()] = true;
  ringAdjTable[bj->getBeginAtomIdx()][bj->getEndAtomIdx()] = true;
  ringAdjTable[bj->getEndAtomIdx()][bj->getBeginAtomIdx()] = true;

  return atoms.size() < numAtoms;
}

void makeRingAtomAdjTable(const ROMol &mol,
                          const boost::dynamic_bitset<> &ringBonds,
                          std::vector<boost::dynamic_bitset<>> &ringAdjTable) {
  ringAdjTable = std::vector<boost::dynamic_bitset<>>(
      mol.getNumAtoms(), boost::dynamic_bitset<>(mol.getNumAtoms()));
  for (const auto bond : mol.bonds()) {
    if (ringBonds[bond->getIdx()]) {
      ringAdjTable[bond->getBeginAtomIdx()][bond->getEndAtomIdx()] = true;
      ringAdjTable[bond->getEndAtomIdx()][bond->getBeginAtomIdx()] = true;
    }
  }
}

// Take out any pairs of ring bonds that don't fragment the molecule.  These
// will be in fused ring systems where one of the bonds is in 1 sub-ring,
// one is in the other.
void findBondPairsThatFragment(
    const ROMol &mol, const boost::dynamic_bitset<> &ringBonds,
    std::vector<std::vector<const Bond *>> &ringBlocks,
    std::vector<std::pair<unsigned int, unsigned int>> &ringBondPairs) {
  std::vector<boost::dynamic_bitset<>> ringAdjTable;
  makeRingAtomAdjTable(mol, ringBonds, ringAdjTable);
  // If all the atoms in the bond are 2 connected, it's a simple ring so
  // nothing to do.
  for (const auto &ringBlock : ringBlocks) {
    bool ok = true;
    boost::dynamic_bitset<> blockAtoms(mol.getNumAtoms());
    for (const auto bond : ringBlock) {
      blockAtoms[bond->getBeginAtomIdx()] = true;
      blockAtoms[bond->getEndAtomIdx()] = true;
      if (ringAdjTable[bond->getBeginAtomIdx()].count() > 2 ||
          ringAdjTable[bond->getEndAtomIdx()].count() > 2) {
        ok = false;
      }
    }
    if (ok) {
      for (size_t i = 0; i < ringBlock.size() - 1; ++i) {
        for (size_t j = i + 1; j < ringBlock.size(); ++j) {
          ringBondPairs.emplace_back(
              std::make_pair(ringBlock[i]->getIdx(), ringBlock[j]->getIdx()));
        }
      }
    } else {
      // Need to check if each pair makes 2 fragments before adding.
      unsigned int numAtoms = blockAtoms.count();
      for (size_t i = 0; i < ringBlock.size() - 1; ++i) {
        for (size_t j = i + 1; j < ringBlock.size(); ++j) {
          if (bondPairFragmentsBlock(i, j, numAtoms, ringBlock, ringAdjTable)) {
            ringBondPairs.emplace_back(
                std::make_pair(ringBlock[i]->getIdx(), ringBlock[j]->getIdx()));
          }
        }
      }
    }
  }
}

void makeFragmentPairs(
    const ROMol &mol, const std::vector<int> &nonRingBonds,
    const std::vector<std::pair<int, int>> &ringBondPairs,
    std::map<std::string, std::vector<std::shared_ptr<ROMol>>> &fragSmis,
    unsigned int dummyStart,
    std::vector<std::vector<std::shared_ptr<ROMol>>> &fragments) {
  // Make a map of the original bond indices to the actual indices
  // in this fragment
  // std::cout << "makeFragmentPairs for " << MolToSmiles(mol) << std::endl;
  std::map<unsigned int, unsigned int> fragBonds;
  for (const auto bond : mol.bonds()) {
    // Bonds to dummy atoms won't have an origIdx.
    if (bond->hasProp("origIdx")) {
      fragBonds.insert(
          std::make_pair(bond->getProp<int>("origIdx"), bond->getIdx()));
    }
  }

  auto addFragmentation = [&](const std::unique_ptr<ROMol> &fragMol,
                              unsigned int ds) {
    std::string fragSmi(MolToSmiles(*fragMol));
    if (auto it = fragSmis.find(fragSmi); it != fragSmis.end()) {
      fragments.emplace_back(it->second);
    } else {
      std::vector<std::unique_ptr<ROMol>> frags;
      MolOps::getMolFrags(*fragMol, frags, false);
      if (frags.size() > 2) {
        std::cout << "not splitting " << MolToSmiles(mol) << " gives "
                  << frags.size() << " frags" << std::endl;
        return;
      }
      fragments.emplace_back();
      std::transform(
          frags.begin(), frags.end(), std::back_inserter(fragments.back()),
          [&](std::unique_ptr<ROMol> &frag) -> std::shared_ptr<ROMol> {
            frag->setProp<unsigned int>("dummyStart", ds);
            return std::shared_ptr<ROMol>(frag.release());
          });
      fragSmis.insert({fragSmi, fragments.back()});
    }
  };

  ++dummyStart;
  std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{
      1, {dummyStart, dummyStart}};
  for (auto bondIdx : nonRingBonds) {
    if (auto it = fragBonds.find(bondIdx); it != fragBonds.end()) {
      std::unique_ptr<ROMol> fragMol(MolFragmenter::fragmentOnBonds(
          mol, std::vector<unsigned int>{it->second}, true, &dummyLabels));
      addFragmentation(fragMol, dummyStart);
    }
  }

  dummyLabels.push_back(std::make_pair(dummyStart + 1, dummyStart + 1));
  for (const auto &bondPairIdxs : ringBondPairs) {
    auto it1 = fragBonds.find(bondPairIdxs.first);
    if (it1 == fragBonds.end()) {
      continue;
    }
    auto it2 = fragBonds.find(bondPairIdxs.second);
    if (it2 == fragBonds.end()) {
      continue;
    }
    std::unique_ptr<ROMol> fragMol(MolFragmenter::fragmentOnBonds(
        mol, std::vector<unsigned int>{it1->second, it2->second}, true,
        &dummyLabels));
    addFragmentation(fragMol, dummyStart + 1);
  }
}

unsigned int maxDummyStart(
    const std::vector<std::shared_ptr<ROMol>> &fragments) {
  unsigned int maxDS = 0;
  for (const auto &frag : fragments) {
    unsigned int dummyStart = frag->getProp<unsigned int>("dummyStart");
    if (dummyStart > maxDS) {
      maxDS = dummyStart;
    }
  }
  return maxDS;
}

void makeFragments(
    const ROMol &mol, const std::vector<unsigned int> &splitBonds,
    const std::vector<std::pair<unsigned int, unsigned int>> &dummyLabels,
    unsigned int maxBondSplits, const boost::dynamic_bitset<> &ringBonds,
    std::set<std::string> &fragSmis,
    std::vector<std::vector<std::shared_ptr<ROMol>>> &fragments) {
  std::cout << "splitting on ";
  for (unsigned int i = 0; i < splitBonds.size(); ++i) {
    std::cout << splitBonds[i] << " ";
  }
  std::cout << std::endl;
  // first, see how many fragments we're going to get. The ring bonds
  // are paired so they will split the same ring.
  int numRingBonds = 0;
  int numNonRingBonds = 0;
  for (auto i : splitBonds) {
    std::cout << i << " " << ringBonds[i] << std::endl;
    if (ringBonds[i]) {
      numRingBonds++;
    } else {
      numNonRingBonds++;
    }
  }
  unsigned int numFragsPoss = 1 + numNonRingBonds + (numRingBonds / 2);
  std::cout << "num frags poss : " << numFragsPoss << " : " << numNonRingBonds
            << " and " << numRingBonds << " vs " << maxBondSplits << " : "
            << MolToSmiles(mol) << std::endl;
  if (numFragsPoss > maxBondSplits) {
    return;
  }
  std::unique_ptr<ROMol> fragMol(
      MolFragmenter::fragmentOnBonds(mol, splitBonds, true, &dummyLabels));
  std::string fragSmi(MolToSmiles(*fragMol));
  if (fragSmis.insert(fragSmi).second) {
    std::vector<std::unique_ptr<ROMol>> molFrags;
    if (MolOps::getMolFrags(*fragMol, molFrags, false) <= maxBondSplits) {
      fragments.emplace_back();
      std::transform(
          molFrags.begin(), molFrags.end(),
          std::back_inserter(fragments.back()),
          [&](std::unique_ptr<ROMol> &frag) -> std::shared_ptr<ROMol> {
            return std::shared_ptr<ROMol>(frag.release());
          });
    }
  }
}

void addSplitBond(
    const std::pair<unsigned int, unsigned int> &bondPair,
    std::vector<unsigned int> &splitBonds,
    std::vector<std::pair<unsigned int, unsigned int>> &dummyLabels) {
  splitBonds.emplace_back(bondPair.first);
  splitBonds.emplace_back(bondPair.second);
  std::sort(splitBonds.begin(), splitBonds.end());
  splitBonds.erase(std::unique(splitBonds.begin(), splitBonds.end()),
                   splitBonds.end());
  dummyLabels.clear();
  for (unsigned int i = 0; i < static_cast<unsigned int>(splitBonds.size());
       ++i) {
    dummyLabels.emplace_back(std::make_pair(i + 1, i + 1));
  }
  // std::cout << "splitBonds: " << splitBonds.size() << " :: ";
  // for (unsigned int i = 0; i < static_cast<unsigned int>(splitBonds.size());
  //      ++i) {
  //   std::cout << splitBonds[i] << " ";
  // }
  // std::cout << std::endl;
  // std::cout << "dummyLabels: " << dummyLabels.size() << " :: ";
  // for (unsigned int i = 0; i < static_cast<unsigned int>(dummyLabels.size());
  //      ++i) {
  //   std::cout << dummyLabels[i].first << "," << dummyLabels[i].second << " ";
  // }
  // std::cout << std::endl;
}

// Build all combinations of maxBondSplits sets of bondPairs into splitBonds,
// removing any duplicate bonds.
void buildSplitBonds(
    const std::vector<std::pair<unsigned int, unsigned int>> &bondPairs,
    unsigned int maxBondSplits,
    std::vector<std::vector<unsigned int>> &splitBonds) {
  std::vector<unsigned int> nextSplits;
  splitBonds.reserve(maxBondSplits * maxBondSplits * bondPairs.size());
  for (unsigned int i = 1; i < maxBondSplits; ++i) {
    auto combs = combMFromN(i, static_cast<int>(bondPairs.size()));
    for (const auto &comb : combs) {
      nextSplits.clear();
      for (auto c : comb) {
        nextSplits.push_back(bondPairs[c].first);
        nextSplits.push_back(bondPairs[c].second);
      }
      std::sort(nextSplits.begin(), nextSplits.end());
      nextSplits.erase(std::unique(nextSplits.begin(), nextSplits.end()),
                       nextSplits.end());
      // Each split will need a connector num, so any split set that will
      // produce one higher than the SynthonSpace has been set up for is
      // a bust.  Splitting 3 rings each once will produce 4 fragments
      // and 6 broken bonds, for example.
      if (nextSplits.size() > MAX_CONNECTOR_NUM) {
        continue;
      }
      nextSplits.shrink_to_fit();
      splitBonds.push_back(nextSplits);
    }
  }
  std::sort(splitBonds.begin(), splitBonds.end());
  splitBonds.erase(std::unique(splitBonds.begin(), splitBonds.end()),
                   splitBonds.end());
  std::cout << "bonds to split with : " << splitBonds.size() << std::endl;
  for (const auto &sb : splitBonds) {
    std::cout << sb.size() << " :: ";
    for (auto s : sb) {
      std::cout << s << " ";
    }
    std::cout << std::endl;
  }
}
}  // namespace

std::vector<std::vector<std::shared_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxNumFrags, std::uint64_t maxNumFragSets,
    TimePoint *endTime, bool &timedOut) {
  if (maxNumFrags < 1) {
    maxNumFrags = 1;
  }
  maxNumFrags = std::min({maxNumFrags, MAX_CONNECTOR_NUM, query.getNumBonds()});

  auto ringBonds = flagRingBonds(query);

  // Now get all contiguous ring blocks
  auto ringBlocks = getRingBlocks(query, ringBonds);

  // Collect all the bond pairs that can fragment the molecule.
  std::vector<std::pair<unsigned int, unsigned int>> bondPairs;
  findBondPairsThatFragment(query, ringBonds, ringBlocks, bondPairs);
  // And all the non-ring bonds, which clearly can all make 2 fragments
  // when broken.  Put them in as pairs of the same value, for ease of
  // processing below.
  for (const auto b : query.bonds()) {
    if (!ringBonds[b->getIdx()]) {
      bondPairs.push_back({b->getIdx(), b->getIdx()});
    }
  }
  std::cout << "Number of splitters : " << bondPairs.size() << std::endl;
  for (const auto &pair : bondPairs) {
    std::cout << pair.first << " -> " << pair.second << std::endl;
  }
  std::vector<std::vector<std::shared_ptr<ROMol>>> fragments;
  // Keep the molecule itself (i.e. 0 splits).  It will probably produce
  // lots of hits but it is necessary if, for example, the query is a match
  // for a single synthon set.
  fragments.emplace_back();
  fragments.back().emplace_back(new ROMol(query));

  // Now do the splits.  Symmetrical molecules can give rise to the same
  // fragment set in different ways so keep track of what we've had to
  // avoid duplicates.
  std::set<std::string> fragSmis;
  std::vector<std::vector<unsigned int>> splitBonds;
  buildSplitBonds(bondPairs, maxNumFrags, splitBonds);
  std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
  for (unsigned int i = 1; i <= MAX_CONNECTOR_NUM; ++i) {
    dummyLabels.emplace_back(i, i);
  }

  for (const auto &sb : splitBonds) {
    makeFragments(query, sb, dummyLabels, maxNumFrags, ringBonds, fragSmis,
                  fragments);
  }

#if 0
  // To make maxBondSplits fragments, we need to do maxBondSplits-1 splits.
  for (size_t i = 0; i < bondPairs.size(); i++) {
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
    std::vector<unsigned int> splitBonds;
    addSplitBond(bondPairs[i], splitBonds, dummyLabels);
    makeFragments(query, splitBonds, dummyLabels, maxBondSplits, ringBonds,
                  fragSmis, fragments);
    if (fragments.size() > maxNumFragSets) {
      BOOST_LOG(rdWarningLog)
          << "Maximum number of fragments reached." << std::endl;
      break;
    }
    for (size_t j = 0; j < bondPairs.size(); j++) {
      splitBonds.clear();
      addSplitBond(bondPairs[i], splitBonds, dummyLabels);
      addSplitBond(bondPairs[j], splitBonds, dummyLabels);
      makeFragments(query, splitBonds, dummyLabels, maxBondSplits, ringBonds,
                    fragSmis, fragments);
      if (fragments.size() > maxNumFragSets) {
        BOOST_LOG(rdWarningLog)
            << "Maximum number of fragments reached." << std::endl;
        break;
      }
      for (size_t k = 0; k < bondPairs.size(); k++) {
        splitBonds.clear();
        addSplitBond(bondPairs[i], splitBonds, dummyLabels);
        addSplitBond(bondPairs[j], splitBonds, dummyLabels);
        addSplitBond(bondPairs[k], splitBonds, dummyLabels);
        makeFragments(query, splitBonds, dummyLabels, maxBondSplits, ringBonds,
                      fragSmis, fragments);
        if (fragments.size() > maxNumFragSets) {
          BOOST_LOG(rdWarningLog)
              << "Maximum number of fragments reached." << std::endl;
          break;
        }
      }
    }
  }
#endif
#if 0
  // To make all fragment pairs, use all pairs of nonRingBonds and all
  // ringBondPairs
  for (size_t splitNum = 1; splitNum < static_cast<size_t>(maxBondSplits);
       ++splitNum) {
    std::vector<std::vector<std::shared_ptr<ROMol>>> newFrags;
    for (auto &frags : fragments) {
      unsigned int dummyStart = maxDummyStart(frags);
      if (frags.size() != splitNum ||
          frags.size() == static_cast<size_t>(maxBondSplits)) {
        continue;
      }
      for (size_t i = 0; i < frags.size(); ++i) {
        // split frag i and make a new load of fragments with the products
        // of the split and the existing ones
        std::vector<std::vector<std::shared_ptr<ROMol>>> iNewFrags;
        makeFragmentPairs(*frags[i], nonRingBonds, ringBondPairs, fragSmis,
                          dummyStart, iNewFrags);
        for (const auto &iNewFrag : iNewFrags) {
          newFrags.emplace_back();
          for (size_t j = 0; j < frags.size(); ++j) {
            if (j != i) {
              newFrags.back().emplace_back(frags[j]);
            } else {
              newFrags.back().insert(newFrags.back().end(), iNewFrag.begin(),
                                     iNewFrag.end());
            }
          }
          if (newFrags.back().size() == 4) {
            std::cout << "BAWOOGA" << std::endl;
          }
        }
      }
    }
    fragments.insert(fragments.end(), newFrags.begin(), newFrags.end());
  }
  std::set<std::string> allSmis;
  std::sort(fragments.begin(), fragments.end());
  fragments.erase(std::unique(fragments.begin(), fragments.end()),
                  fragments.end());
  std::sort(fragments.begin(), fragments.end(),
            [](const std::vector<std::shared_ptr<ROMol>> &a,
               const std::vector<std::shared_ptr<ROMol>> &b) -> bool {
              return a.size() < b.size();
            });
#endif
  std::cout << "total number of frag sets : " << fragments.size() << std::endl;
  for (const auto &frags : fragments) {
    for (const auto &f : frags) {
      std::cout << MolToSmiles(*f) << ".";
    }
    std::cout << std::endl;
  }
  return fragments;
}

int countConnections(const ROMol &mol) {
  int res = 0;
  for (const auto atom : mol.atoms()) {
    if (!atom->getAtomicNum() && atom->getIsotope() >= 1 &&
        atom->getIsotope() <= MAX_CONNECTOR_NUM) {
      ++res;
    }
  }
  return res;
}

std::vector<boost::dynamic_bitset<>> getConnectorPatterns(
    const std::vector<std::shared_ptr<ROMol>> &mols) {
  std::vector<boost::dynamic_bitset<>> connPatterns(
      mols.size(), boost::dynamic_bitset<>(MAX_CONNECTOR_NUM + 1));
  for (size_t i = 0; i < mols.size(); i++) {
    for (const auto &a : mols[i]->atoms()) {
      if (!a->getAtomicNum() && a->getIsotope() <= MAX_CONNECTOR_NUM) {
        connPatterns[i].set(a->getIsotope());
      }
    }
  }
  return connPatterns;
}

boost::dynamic_bitset<> getConnectorPattern(
    const std::vector<std::shared_ptr<ROMol>> &fragSet) {
  boost::dynamic_bitset<> conns(MAX_CONNECTOR_NUM + 1);
  const auto connPatterns = getConnectorPatterns(fragSet);
  for (const auto &cp : connPatterns) {
    conns |= cp;
  }
  return conns;
}

namespace {
std::vector<int> bitsToInts(const boost::dynamic_bitset<> &bits) {
  std::vector<int> ints;
  for (size_t i = 0; i < bits.size(); ++i) {
    if (bits[i]) {
      ints.push_back(static_cast<int>(i));
    }
  }
  return ints;
}
}  // namespace

std::vector<std::vector<std::unique_ptr<ROMol>>> getConnectorPermutations(
    const std::vector<std::shared_ptr<ROMol>> &molFrags,
    const boost::dynamic_bitset<> &fragConns,
    const boost::dynamic_bitset<> &reactionConns) {
  std::vector<std::vector<std::unique_ptr<ROMol>>> connPerms;
  const auto numFragConns = fragConns.count();
  auto rConns = bitsToInts(reactionConns);
  const auto perms = permMFromN(numFragConns, reactionConns.count());

  for (const auto &perm : perms) {
    connPerms.emplace_back();
    // Copy the fragments and set the isotope numbers according to this
    // permutation.
    for (const auto &f : molFrags) {
      connPerms.back().emplace_back(new RWMol(*f));
      boost::dynamic_bitset<> atomDone(f->getNumAtoms());
      for (const auto atom : connPerms.back().back()->atoms()) {
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

std::vector<std::vector<boost::dynamic_bitset<>>> getConnectorPermutations(
    const std::vector<boost::dynamic_bitset<>> &fragConnPatts,
    const boost::dynamic_bitset<> &reactionConns) {
  boost::dynamic_bitset<> conns(MAX_CONNECTOR_NUM + 1);
  for (auto &fragConnPatt : fragConnPatts) {
    conns |= fragConnPatt;
  }

  const auto numFragConns = conns.count();
  auto rConns = bitsToInts(reactionConns);
  const auto perms = permMFromN(numFragConns, reactionConns.count());
  std::vector<std::vector<boost::dynamic_bitset<>>> retBitsets;
  for (const auto &perm : perms) {
    retBitsets.emplace_back();
    for (const auto &fragConnPatt : fragConnPatts) {
      boost::dynamic_bitset<> bs(MAX_CONNECTOR_NUM + 1);
      for (size_t i = 0; i < perm.size(); ++i) {
        if (fragConnPatt[i + 1]) {
          bs.set(perm[i] + 1);
        }
      }
      retBitsets.back().push_back(bs);
    }
  }
  return retBitsets;
}

void expandBitSet(std::vector<boost::dynamic_bitset<>> &bitSets) {
  const bool someSet = std::any_of(
      bitSets.begin(), bitSets.end(),
      [](const boost::dynamic_bitset<> &bs) -> bool { return bs.any(); });
  if (someSet) {
    for (auto &bs : bitSets) {
      if (!bs.count()) {
        bs.set();
      }
    }
  }
}

void bitSetsToVectors(const std::vector<boost::dynamic_bitset<>> &bitSets,
                      std::vector<std::vector<size_t>> &outVecs) {
  outVecs.resize(bitSets.size());
  for (size_t i = 0; i < bitSets.size(); ++i) {
    outVecs[i].reserve(bitSets[i].count());
    for (size_t j = 0; j < bitSets[i].size(); j++) {
      if (bitSets[i][j]) {
        outVecs[i].push_back(j);
      }
    }
  }
}

std::string buildProductName(const std::string &reactionId,
                             const std::vector<std::string> &fragIds) {
  std::string prodName = "";
  for (const auto &fragId : fragIds) {
    if (prodName != "") {
      prodName += ";";
    }
    prodName += fragId;
  }
  prodName += ";" + reactionId;
  return prodName;
}

std::unique_ptr<ROMol> buildProduct(const std::vector<const ROMol *> &synths) {
  MolzipParams mzparams;
  mzparams.label = MolzipLabel::Isotope;

  auto prodMol = std::make_unique<ROMol>(*synths.front());
  for (size_t i = 1; i < synths.size(); ++i) {
    prodMol.reset(combineMols(*prodMol, *synths[i]));
  }
  prodMol = molzip(*prodMol, mzparams);
  MolOps::sanitizeMol(*dynamic_cast<RWMol *>(prodMol.get()));

  return prodMol;
}
}  // namespace RDKit::SynthonSpaceSearch::details
