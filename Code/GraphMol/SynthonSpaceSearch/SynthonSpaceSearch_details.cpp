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
#include <regex>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <RDGeneral/ControlCHandler.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>

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

std::vector<std::vector<std::unique_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxBondSplits, std::uint64_t maxNumFrags,
    TimePoint *endTime, bool &timedOut) {
  if (maxBondSplits < 1) {
    maxBondSplits = 1;
  }
  maxBondSplits =
      std::min({maxBondSplits, MAX_CONNECTOR_NUM, query.getNumBonds()});
  const auto ringInfo = query.getRingInfo();
  if (!ringInfo->isInitialized()) {
    // Query molecules don't seem to have the ring info generated on creation.
    MolOps::findSSSR(query);
  }
  boost::dynamic_bitset<> ringBonds(query.getNumBonds());
  for (const auto &r : ringInfo->bondRings()) {
    for (const auto b : r) {
      ringBonds.set(b);
    }
  }
  std::vector<std::vector<std::unique_ptr<ROMol>>> fragments;
  // Keep the molecule itself (i.e. 0 splits).  It will probably produce
  // lots of hits but it is necessary if, for example, the query is a match
  // for a single synthon set.
  fragments.emplace_back();
  fragments.back().emplace_back(new ROMol(query));

  // Now do the splits.  Symmetrical molecules can give rise to the same
  // fragment set in different ways so keep track of what we've had to
  // avoid duplicates.
  std::set<std::string> fragSmis;
  bool cancelled = false;
  timedOut = false;
  std::uint64_t numTries = 100;

  // Now do the splits.
  for (unsigned int i = 1; i <= maxBondSplits; ++i) {
    if (timedOut || cancelled) {
      break;
    }
    auto combs = combMFromN(i, static_cast<int>(query.getNumBonds()));
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
    for (unsigned int j = 1; j <= i; ++j) {
      dummyLabels.emplace_back(j, j);
    }
    for (auto &c : combs) {
      if (ControlCHandler::getGotSignal()) {
        cancelled = true;
        break;
      }
      --numTries;
      if (!numTries) {
        numTries = 100;
        timedOut = checkTimeOut(endTime);
        if (timedOut) {
          break;
        }
      }

      // don't break just 1 ring bond, as it can't create 2 fragments.  It
      // could be better than this, by checking that any number of ring
      // bonds are all in the same ring system.  Maybe look at that
      // if necessary for performance.  Triazoles can be created from 3
      // synthons, so breaking 3 ring bonds is ok.  This will still pass
      // through cases that break 2 ring bonds in separate ring systems,
      // such as a bond in each of the 2 phenyl rings in c1ccccc1c2ccccc2,
      // but they will be caught below.
      const auto numRingBonds = std::reduce(
          c.begin(), c.end(), 0, [&](const int prevRes, const int bondNum) {
            if (ringBonds[bondNum]) {
              return prevRes + 1;
            }
            return prevRes;
          });
      if (numRingBonds == 1) {
        continue;
      }
      std::unique_ptr<ROMol> fragMol(
          MolFragmenter::fragmentOnBonds(query, c, true, &dummyLabels));
      std::vector<std::unique_ptr<ROMol>> molFrags;
      // Must have been a ring-opening.
      if (const auto numFrags = MolOps::getMolFrags(*fragMol, molFrags, false);
          numFrags == 1) {
        continue;
      }
      if (checkConnectorsInDifferentFrags(molFrags, i)) {
        std::string fragSmi(MolToSmiles(*fragMol));
        if (!fragSmis.insert(fragSmi).second) {
          continue;
        }
        fragments.emplace_back(std::move(molFrags));
        if (fragments.size() > maxNumFrags) {
          BOOST_LOG(rdWarningLog)
              << "Maximum number of fragments reached." << std::endl;
          break;
        }
      }
    }
    if (fragments.size() > maxNumFrags) {
      break;
    }
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

std::vector<std::vector<std::unique_ptr<ROMol>>> getConnectorPermutations(
    const std::vector<std::unique_ptr<ROMol>> &molFrags,
    const boost::dynamic_bitset<> &fragConns,
    const boost::dynamic_bitset<> &reactionConns) {
  std::vector<std::vector<std::unique_ptr<ROMol>>> connPerms;
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

}  // namespace RDKit::SynthonSpaceSearch::details
