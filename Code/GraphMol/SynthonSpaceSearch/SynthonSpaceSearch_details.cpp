//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "SynthonSpace.h"

#include <algorithm>
#include <list>
#include <regex>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>

namespace RDKit::SynthonSpaceSearch::details {

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

std::vector<std::vector<unsigned int>> permMFromN(unsigned int m,
                                                  unsigned int n) {
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
    const std::vector<std::unique_ptr<ROMol>> &molFrags, int numSplits) {
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
      auto bond = mol.getBondBetweenAtoms(nextBond->getBeginAtomIdx(), nbr);
      if (!done[bond->getIdx()] && bond->getIsAromatic()) {
        aromBonds.push_back(bond);
        done[bond->getIdx()] = true;
        toDo.push_back(bond);
      }
    }
    for (const auto nbr :
         make_iterator_range(mol.getAtomNeighbors(nextBond->getEndAtom()))) {
      auto bond = mol.getBondBetweenAtoms(nextBond->getEndAtomIdx(), nbr);
      if (!done[bond->getIdx()] && bond->getIsAromatic()) {
        aromBonds.push_back(bond);
        done[bond->getIdx()] = true;
        toDo.push_back(bond);
      }
    }
  }
  return aromBonds;
}

// If a bond has been split in an aromatic ring, the bonds to dummy atoms
// will be aromatic.  However, the synthons that produce them must have
// non-aromatic bonds to the dummies, by definition.  An example is in
// triazole formation, where the synthons can be [2*]=N-N=C-[1*],
// [1*]N[3*] and [2*]=C[3*] which come together to give C1=NN=CN1 so long
// as there are appropriate substituents on the synthons.  Splitting triazole,
// however, gives [1*]cnn[3*], [2*]n[3*] and [1*]c[2*] (no explicit H on the
// nitrogen if the molecule is specified as a SMARTS originally).  These won't
// match the synthons, so no hit is found.  This function detects such
// situations and changes the bond types accordingly.  All aromatic bonds are
// set to single|double|aromatic because we can't know which kekule form would
// be appropriate for the synthons, and aromatic atoms are set to atomic number
// queries with the aromatic flag cleared.
void fixAromaticRingSplits(std::vector<std::unique_ptr<ROMol>> &molFrags) {
  for (auto &frag : molFrags) {
    auto buildQueryAtom = [](const Atom *atom) -> std::unique_ptr<QueryAtom> {
      std::unique_ptr<QueryAtom> nqa(new QueryAtom(atom->getAtomicNum()));
      if (!nqa->getAtomicNum()) {
        nqa->setIsotope(atom->getIsotope());
        nqa->expandQuery(makeAtomIsotopeQuery(atom->getIsotope()));
      }
      return nqa;
    };
    std::unique_ptr<RWMol> qmol;
    for (const auto atom : frag->atoms()) {
      // Allow for general dummy atoms in the query by only looking at atoms
      // where the isotope numbers have been set to the ones we're using.  For
      // these atoms, there should only be 1 bond.
      if (!atom->getAtomicNum() && atom->getIsotope() > 0 &&
          atom->getIsotope() < MAX_CONNECTOR_NUM + 1) {
        if (auto const fbond = (*frag)[*frag->getAtomBonds(atom).first];
            fbond->getIsAromatic()) {
          if (!qmol) {
            qmol = std::make_unique<RWMol>(*frag);
          }
          auto aromBonds = getContiguousAromaticBonds(*frag, fbond);
          for (const auto &ab : aromBonds) {
            const auto qab = qmol->getBondWithIdx(ab->getIdx());
            std::unique_ptr<QueryBond> qbond(new QueryBond(*qab));
            qbond->setQuery(makeSingleOrDoubleOrAromaticBondQuery());
            qmol->replaceBond(qab->getIdx(), qbond.get());

            const auto qba = qmol->getAtomWithIdx(ab->getBeginAtomIdx());
            auto nqba = buildQueryAtom(qba);
            qmol->replaceAtom(ab->getBeginAtomIdx(), nqba.get());

            const auto qea = qmol->getAtomWithIdx(ab->getEndAtomIdx());
            const auto nqea = buildQueryAtom(qea);
            qmol->replaceAtom(ab->getEndAtomIdx(), nqea.get());
          }
          break;
        }
      }
    }

    if (qmol) {
      frag = std::move(qmol);
    }
  }
}

std::vector<std::vector<std::unique_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxBondSplits) {
  if (maxBondSplits < 1) {
    maxBondSplits = 1;
  }
  maxBondSplits =
      std::min({maxBondSplits, MAX_CONNECTOR_NUM, query.getNumBonds()});
  const auto ringInfo = query.getRingInfo();
  boost::dynamic_bitset<> ringBonds(query.getNumBonds());
  for (const auto &r : ringInfo->bondRings()) {
    for (const auto b : r) {
      ringBonds.set(b);
    }
  }

  std::vector<std::vector<std::unique_ptr<ROMol>>> fragments;
  // Keep the molecule itself (i.e. 0 splits).  It will probably produce
  // lots of hits but one can imagine a use for it.
  fragments.emplace_back();
  fragments.back().emplace_back(new ROMol(query));

  // Now do the splits.
  for (unsigned int i = 1; i <= maxBondSplits; ++i) {
    auto combs = combMFromN(i, static_cast<int>(query.getNumBonds()));
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
    for (unsigned int j = 1; j <= i; ++j) {
      dummyLabels.emplace_back(j, j);
    }
    for (auto &c : combs) {
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
      auto numFrags = MolOps::getMolFrags(*fragMol, molFrags, false);
      // Must have been a ring-opening.
      if (numFrags == 1) {
        continue;
      }
      if (checkConnectorsInDifferentFrags(molFrags, i)) {
        fragments.emplace_back(std::move(molFrags));
      }
    }
  }
  return fragments;
}

int countConnections(const std::string &smiles) {
  static const std::regex conns(R"(\[[1-4]\*\])");
  return static_cast<int>(std::distance(
      std::sregex_token_iterator(smiles.begin(), smiles.end(), conns),
      std::sregex_token_iterator()));
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
  auto connPatterns = getConnectorPatterns(fragSet);
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

void expandBitSet(std::vector<boost::dynamic_bitset<>> &bitSets) {
  bool someSet = std::any_of(
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

}  // namespace RDKit::SynthonSpaceSearch::details