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
#include <list>
#include <memory>
#include <regex>
#include <thread>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <RDGeneral/ControlCHandler.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceHitSet.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <RDGeneral/RDThreads.h>

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
  int dummyAtoms[2 * MAX_CONNECTOR_NUM + 2] = {};
  for (const auto &atom : mol.atoms()) {
    if (!atom->getAtomicNum()) {
      if (const auto dummy = atom->getIsotope(); dummy <= MAX_CONNECTOR_NUM) {
        if (const int pos = 2 * dummy; dummyAtoms[pos]) {
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
      if (const Atom *otherAtom = nbond->getOtherAtom(atom);
          !doneAtoms[otherAtom->getIdx()]) {
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
          const auto nextAtom = toDo.front();
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
  return ringBlocks;
}

bool bondPairFragmentsBlock(
    const size_t bondi, const size_t bondj, const unsigned int numAtoms,
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
    const std::vector<std::vector<const Bond *>> &ringBlocks,
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
      const unsigned int numAtoms = blockAtoms.count();
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

void makeFragmentsForMol(
    const ROMol &mol, const std::vector<std::vector<unsigned int>> &splitBonds,
    size_t splitBondNum,
    const std::vector<std::pair<unsigned int, unsigned int>> &dummyLabels,
    const unsigned int maxNumFrags, const boost::dynamic_bitset<> &ringBonds,
    std::vector<std::pair<std::string, std::unique_ptr<ROMol>>> &fragments) {
  // first, see how many fragments we're going to get. The ring bonds
  // are paired so they will split the same ring.
  int numRingBonds = 0;
  int numNonRingBonds = 0;
  for (const auto sb : splitBonds[splitBondNum]) {
    if (ringBonds[sb]) {
      numRingBonds++;
    } else {
      numNonRingBonds++;
    }
  }
  if (const unsigned int numFragsPoss = 1 + numNonRingBonds + numRingBonds / 2;
      numFragsPoss > maxNumFrags) {
    return;
  }
  auto fragMol = MolFragmenter::fragmentOnBonds(mol, splitBonds[splitBondNum],
                                                true, &dummyLabels);
  const std::string fragSmi(MolToSmiles(*fragMol));
  fragments[splitBondNum] =
      std::pair<std::string, std::unique_ptr<ROMol>>(fragSmi, fragMol);
}

void doPartInitialFragmentation(
    const ROMol &mol, const std::vector<std::vector<unsigned int>> &splitBonds,
    const unsigned int maxNumFrags, const boost::dynamic_bitset<> &ringBonds,
    const TimePoint *endTime, std::atomic<std::int64_t> &mostRecentRingBond,
    std::int64_t lastRingBond,
    const std::vector<std::pair<unsigned int, unsigned int>> &dummyLabels,
    std::vector<std::pair<std::string, std::unique_ptr<ROMol>>> &tmpFrags) {
  int numTries = 100;
  bool timedOut = false;
  while (true) {
    std::int64_t thisRB = ++mostRecentRingBond;
    // std::cout << "ring bond " << thisRB << " of " << lastRingBond << " and "
    // << tmpFrags.size() << std::endl;
    if (thisRB > lastRingBond) {
      break;
    }
    makeFragmentsForMol(mol, splitBonds, thisRB, dummyLabels, maxNumFrags,
                        ringBonds, tmpFrags);
    --numTries;
    if (!numTries) {
      numTries = 100;
      timedOut = checkTimeOut(endTime);
      if (timedOut) {
        break;
      }
    }
  }
}

void doInitialFragmentation(
    const ROMol &mol, const std::vector<std::vector<unsigned int>> &splitBonds,
    const unsigned int maxNumFrags, const boost::dynamic_bitset<> &ringBonds,
    [[maybe_unused]] const int numThreads, const TimePoint *endTime,
    bool &timedOut,
    std::vector<std::pair<std::string, std::unique_ptr<ROMol>>> &tmpFrags) {
  std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
  for (unsigned int i = 1; i <= MAX_CONNECTOR_NUM; ++i) {
    dummyLabels.emplace_back(i, i);
  }

  // Now do the splits.  Symmetrical molecules can give rise to the same
  // fragment set in different ways so keep track of what we've had to
  // avoid duplicates.
  std::int64_t lastRingBond = splitBonds.size() - 1;
  std::atomic<std::int64_t> mostRecentRingBond = -1;
#if RDK_BUILD_THREADSAFE_SSS
  if (const auto numThreadsToUse = getNumThreadsToUse(numThreads);
      numThreads > 1) {
    std::vector<std::thread> threads;
    for (unsigned int i = 0U;
         i <
         std::min(static_cast<std::int64_t>(numThreadsToUse), lastRingBond + 1);
         ++i) {
      threads.push_back(std::thread(doPartInitialFragmentation, std::ref(mol),
                                    std::ref(splitBonds), maxNumFrags,
                                    std::ref(ringBonds), endTime,
                                    std::ref(mostRecentRingBond), lastRingBond,
                                    std::ref(dummyLabels), std::ref(tmpFrags)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    doPartInitialFragmentation(mol, splitBonds, maxNumFrags, ringBonds, endTime,
                               std::ref(mostRecentRingBond), lastRingBond,
                               dummyLabels, tmpFrags);
  }
#else
  doPartInitialFragmentation(mol, splitBonds, maxNumFrags, ringBonds, endTime,
                             std::ref(mostRecentRingBond), lastRingBond,
                             dummyLabels, tmpFrags);
#endif
  timedOut = details::checkTimeOut(endTime);
}

void doPartFinalFragmentation(
    const std::vector<std::pair<std::string, std::unique_ptr<ROMol>>> &tmpFrags,
    unsigned int maxNumFrags, const TimePoint *endTime,
    std::atomic<std::int64_t> &mostRecentFrag, std::int64_t lastFrag,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragments) {
  int numTries = 100;
  bool timedOut = false;
  while (true) {
    std::int64_t thisFrag = ++mostRecentFrag;
    // std::cout << "frag " << thisFrag << " of " << lastFrag << " and "
    // << tmpFrags.size() << std::endl;
    if (thisFrag > lastFrag) {
      break;
    }
    if (std::vector<std::unique_ptr<ROMol>> molFrags;
        MolOps::getMolFrags(*tmpFrags[thisFrag].second, molFrags, false) <=
        maxNumFrags) {
      // The first fragment was made from the whole query and is already in
      // fragments.
      fragments[thisFrag + 1] = std::move(molFrags);
    }
    --numTries;
    if (!numTries) {
      numTries = 100;
      timedOut = checkTimeOut(endTime);
      if (timedOut) {
        break;
      }
    }
  }
}

void doFinalFragmentation(
    const std::vector<std::pair<std::string, std::unique_ptr<ROMol>>> &tmpFrags,
    unsigned int maxNumFrags, [[maybe_unused]] int numThreads,
    const TimePoint *endTime, bool &timedOut,
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragments) {
  std::int64_t lastFrag = tmpFrags.size() - 1;
  std::atomic<std::int64_t> mostRecentFrag = -1;
#if RDK_BUILD_THREADSAFE_SSS
  if (const auto numThreadsToUse = getNumThreadsToUse(numThreads);
      numThreads > 1) {
    std::vector<std::thread> threads;
    for (unsigned int i = 0U;
         i < std::min(static_cast<std::int64_t>(numThreadsToUse), lastFrag + 1);
         ++i) {
      threads.push_back(std::thread(
          doPartFinalFragmentation, std::ref(tmpFrags), maxNumFrags, endTime,
          std::ref(mostRecentFrag), lastFrag, std::ref(fragments)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    doPartFinalFragmentation(tmpFrags, maxNumFrags, endTime, mostRecentFrag,
                             lastFrag, fragments);
  }
#else
  doPartFinalFragmentation(tmpFrags, maxNumFrags, endTime, mostRecentFrag,
                           lastFrag, fragments);
#endif
  timedOut = details::checkTimeOut(endTime);
}

// Build all combinations of maxBondSplits sets of bondPairs into splitBonds,
// removing any duplicate bonds.
void buildSplitBonds(
    const std::vector<std::pair<unsigned int, unsigned int>> &bondPairs,
    const unsigned int maxBondSplits,
    std::vector<std::vector<unsigned int>> &splitBonds) {
  std::vector<unsigned int> nextSplits;
  splitBonds.reserve(maxBondSplits * maxBondSplits * bondPairs.size());
  for (unsigned int i = 1; i < maxBondSplits; ++i) {
    auto combs = combMFromN(i, static_cast<int>(bondPairs.size()));
    for (const auto &comb : combs) {
      nextSplits.clear();
      for (const auto c : comb) {
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
}
}  // namespace

std::vector<std::vector<std::unique_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxNumFrags,
    const std::uint64_t maxNumFragSets, const TimePoint *endTime,
    const int numThreads, bool &timedOut) {
  if (maxNumFrags < 1) {
    maxNumFrags = 1;
  }
  maxNumFrags = std::min({maxNumFrags, MAX_CONNECTOR_NUM, query.getNumBonds()});

  auto ringBonds = flagRingBonds(query);

  // Now get all contiguous ring blocks
  const auto ringBlocks = getRingBlocks(query, ringBonds);

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

  std::vector<std::vector<unsigned int>> splitBonds;
  buildSplitBonds(bondPairs, maxNumFrags, splitBonds);
  std::vector<std::pair<std::string, std::unique_ptr<ROMol>>> tmpFrags(
      splitBonds.size());

  // First split leaves the fragments in the same molecule, and returns
  // the SMILES for it.
  doInitialFragmentation(query, splitBonds, maxNumFrags, ringBonds, numThreads,
                         endTime, timedOut, tmpFrags);
  std::vector<std::vector<std::unique_ptr<ROMol>>> fragments;
  if (timedOut || ControlCHandler::getGotSignal()) {
    return fragments;
  }

  // Keep unique SMILES onlyu
  std::sort(tmpFrags.begin(), tmpFrags.end(),
            [](const auto &lhs, const auto &rhs) -> bool {
              return lhs.first < rhs.first;
            });
  tmpFrags.erase(std::unique(tmpFrags.begin(), tmpFrags.end(),
                             [](const auto &lhs, const auto &rhs) -> bool {
                               return lhs.first == rhs.first;
                             }),
                 tmpFrags.end());
  if (tmpFrags.size() > maxNumFragSets) {
    tmpFrags.erase(tmpFrags.begin() + maxNumFragSets, tmpFrags.end());
  }

  // Keep the molecule itself (i.e. 0 splits).  It will probably produce
  // lots of hits but it is necessary if, for example, the query is a match
  // for a single synthon set.
  fragments.resize(tmpFrags.size() + 1);
  fragments.emplace_back();
  fragments.back().emplace_back(new ROMol(query));
  // And now split the molecules into the final fragments.
  doFinalFragmentation(tmpFrags, maxNumFrags, numThreads, endTime, timedOut,
                       fragments);

  fragments.erase(
      std::remove_if(fragments.begin(), fragments.end(),
                     [](const auto &fs) -> bool { return fs.empty(); }),
      fragments.end());
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
    const std::vector<std::unique_ptr<ROMol>> &mols) {
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
    const std::vector<std::unique_ptr<ROMol>> &fragSet) {
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

std::vector<std::vector<std::vector<std::pair<Atom *, unsigned int>>>>
getConnectorPermutations(const std::vector<std::unique_ptr<ROMol>> &molFrags,
                         const boost::dynamic_bitset<> &fragConns,
                         const boost::dynamic_bitset<> &reactionConns) {
  const auto numFragConns = fragConns.count();
  auto rConns = bitsToInts(reactionConns);
  const auto perms = permMFromN(numFragConns, reactionConns.count());

  std::vector<std::vector<std::vector<std::pair<Atom *, unsigned int>>>>
      fragConnPerms;
  fragConnPerms.reserve(perms.size());

  for (const auto &perm : perms) {
    fragConnPerms.emplace_back();
    // Copy the fragments and set the isotope numbers according to this
    // permutation.
    for (const auto &f : molFrags) {
      fragConnPerms.back().emplace_back();
      boost::dynamic_bitset<> atomDone(f->getNumAtoms());
      for (const auto atom : f->atoms()) {
        if (!atom->getAtomicNum()) {
          for (size_t i = 0; i < perm.size(); ++i) {
            if (!atomDone[atom->getIdx()] && atom->getIsotope() == i + 1) {
              fragConnPerms.back().back().emplace_back(atom, perm[i] + 1);
              atomDone[atom->getIdx()] = true;
            }
          }
        }
      }
    }
  }

  return fragConnPerms;
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

bool removeQueryAtoms(RWMol &mol) {
  bool didSomething = false;
  for (const Atom *atom : mol.atoms()) {
    if ((atom->getAtomicNum() || !atom->getIsotope()) && atom->hasQuery() &&
        atom->getQuery()->getDescription() != "AtomType") {
      std::unique_ptr<QueryAtom> qat;
      if (atom->getAtomicNum()) {
        qat.reset(new QueryAtom(atom->getAtomicNum()));
      } else {
        qat.reset(new QueryAtom());
        qat->setQuery(makeAAtomQuery());
      }
      mol.replaceAtom(atom->getIdx(), qat.get());
      didSomething = true;
    }
  }
  return didSomething;
}

std::unique_ptr<ROMol> buildConnRegion(const ROMol &mol) {
  boost::dynamic_bitset<> inFrag(mol.getNumAtoms());
  for (const auto a : mol.atoms()) {
    if (!a->getAtomicNum() && a->getIsotope()) {
      inFrag[a->getIdx()] = true;
      for (const auto &n1 : mol.atomNeighbors(a)) {
        inFrag[n1->getIdx()] = true;
        for (const auto &n2 : mol.atomNeighbors(n1)) {
          inFrag[n2->getIdx()] = true;
          for (const auto &n3 : mol.atomNeighbors(n2)) {
            inFrag[n3->getIdx()] = true;
          }
        }
      }
    }
  }
  if (!inFrag.count()) {
    return std::unique_ptr<RWMol>();
  }

  std::unique_ptr<RWMol> molCp(new RWMol(mol));
  molCp->beginBatchEdit();
  for (const auto aCp : molCp->atoms()) {
    if (!inFrag[aCp->getIdx()]) {
      molCp->removeAtom(aCp);
    } else {
      if (!aCp->getAtomicNum()) {
        if (aCp->getIsotope()) {
          aCp->setIsotope(1);
          if (aCp->hasQuery()) {
            aCp->expandQuery(makeAtomIsotopeQuery(1), Queries::COMPOSITE_OR);
          }
        }
      }
    }
  }
  molCp->commitBatchEdit();
  return molCp;
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

std::string buildProductName(
    const RDKit::SynthonSpaceSearch::SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &fragNums) {
  std::string prodName = "";
  for (size_t i = 0; i < fragNums.size(); ++i) {
    if (prodName != "") {
      prodName += ";";
    }
    prodName += hitset->synthonsToUse[i][fragNums[i]].first;
  }
  prodName += ";" + hitset->d_reaction->getId();
  return prodName;
}

std::unique_ptr<ROMol> buildProduct(
    const std::vector<const ROMol *> &synthons) {
  MolzipParams mzparams;
  mzparams.label = MolzipLabel::Isotope;

  auto prodMol = std::make_unique<ROMol>(*synthons.front());
  for (size_t i = 1; i < synthons.size(); ++i) {
    prodMol.reset(combineMols(*prodMol, *synthons[i]));
  }
  prodMol = molzip(*prodMol, mzparams);
  MolOps::sanitizeMol(*dynamic_cast<RWMol *>(prodMol.get()));

  return prodMol;
}

std::map<std::string, std::vector<ROMol *>> mapFragsBySmiles(
    std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets,
    bool &cancelled) {
  std::map<std::string, std::vector<ROMol *>> fragSmiToFrag;
  for (auto &fragSet : fragSets) {
    for (auto &frag : fragSet) {
      if (ControlCHandler::getGotSignal()) {
        cancelled = true;
        return fragSmiToFrag;
      }
      // For the fingerprints, ring info is required.
      unsigned int otf;
      sanitizeMol(*static_cast<RWMol *>(frag.get()), otf,
                  MolOps::SANITIZE_SYMMRINGS);
      std::string fragSmi = MolToSmiles(*frag);
      if (auto it = fragSmiToFrag.find(fragSmi); it == fragSmiToFrag.end()) {
        fragSmiToFrag.emplace(fragSmi, std::vector<ROMol *>(1, frag.get()));
      } else {
        it->second.emplace_back(frag.get());
      }
    }
  }
  return fragSmiToFrag;
}

}  // namespace RDKit::SynthonSpaceSearch::details
