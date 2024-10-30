//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file contains an implementation of synthonspace substructure search
// similar to that described in
// 'Fast Substructure Search in Combinatorial Library Spaces',
// Thomas Liphardt and Thomas Sander,
// J. Chem. Inf. Model. 2023, 63, 16, 5133–5141
// https://doi.org/10.1021/acs.jcim.3c00290
//
// The algorithm allows the substructure searching of a very large library
// of structures that is described in synthon format (such as Enamine REAL)
// without enumerating the individual structures during the search process.
//
// It is not a direct implementation, as, for example, it uses a different
// fingerprint for the initial synthon screening.

#include <algorithm>
#include <list>
#include <regex>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>

namespace RDKit {
namespace SynthonSpaceSearch {
namespace details {

// get a vector of vectors of unsigned ints that are all combinations of
// M items chosen from N e.g. all combinations of 3 bonds from a
// molecule. A modified form of the code in the first answer from
// https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
std::vector<std::vector<unsigned int>> combMFromN(unsigned int m,
                                                  unsigned int n) {
  std::string allN(m, 1);
  allN.resize(n, 0);
  std::vector<std::vector<unsigned int>> combs;
  do {
    combs.push_back(std::vector<unsigned int>());
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
  //  std::cout << MolToSmiles(mol) << " : " << aromBond->getIdx() << " : "
  //            << aromBond->getBeginAtomIdx() << " ("
  //            << aromBond->getBeginAtom()->getAtomicNum()
  //            << ") : " << aromBond->getEndAtomIdx() << " ("
  //            << aromBond->getEndAtom()->getAtomicNum() << ")" << std::endl;
  std::vector<const Bond *> aromBonds(1, aromBond);
  std::list<const Bond *> toDo(1, aromBond);
  boost::dynamic_bitset<> done(mol.getNumBonds());
  done[aromBond->getIdx()] = true;
  while (!toDo.empty()) {
    auto nextBond = toDo.front();
    //    std::cout << "  next bond : " << " : " << nextBond->getIdx() << " : "
    //              << nextBond->getBeginAtomIdx() << " ("
    //              << nextBond->getBeginAtom()->getAtomicNum()
    //              << ") : " << nextBond->getEndAtomIdx() << " ("
    //              << nextBond->getEndAtom()->getAtomicNum() << ")" <<
    //              std::endl;
    toDo.pop_front();
    for (auto nbr :
         make_iterator_range(mol.getAtomNeighbors(nextBond->getBeginAtom()))) {
      auto bond = mol.getBondBetweenAtoms(nextBond->getBeginAtomIdx(), nbr);
      if (!done[bond->getIdx()] && bond->getIsAromatic()) {
        aromBonds.push_back(bond);
        done[bond->getIdx()] = true;
        toDo.push_back(bond);
      }
    }
    for (auto nbr :
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
  ROMol::OEDGE_ITER beg;
  for (auto &frag : molFrags) {
    auto buildQueryAtom = [](const Atom *atom) -> QueryAtom * {
      QueryAtom *nqa = new QueryAtom(atom->getAtomicNum());
      if (!nqa->getAtomicNum()) {
        nqa->setIsotope(atom->getIsotope());
        nqa->expandQuery(makeAtomIsotopeQuery(atom->getIsotope()));
      }
      return nqa;
    };
    std::unique_ptr<RWMol> qmol;
    for (auto atom : frag->atoms()) {
      // Allow for general dummy atoms in the query by only looking at atoms
      // where the isotope numbers have been set to the ones we're using.  For
      // these atoms, there should only be 1 bond.
      if (!atom->getAtomicNum() && atom->getIsotope() > 0 &&
          atom->getIsotope() < 5) {
        beg = frag->getAtomBonds(atom).first;
        auto fbond = (*frag)[*beg];
        if (fbond->getIsAromatic()) {
          if (!qmol) {
            qmol.reset(new RWMol(*frag));
          }
          auto aromBonds = getContiguousAromaticBonds(*frag, fbond);
          for (const auto &ab : aromBonds) {
            auto qab = qmol->getBondWithIdx(ab->getIdx());
            QueryBond *qbond = new QueryBond(*qab);
            qbond->setQuery(makeSingleOrDoubleOrAromaticBondQuery());
            qmol->replaceBond(qab->getIdx(), qbond);

            auto qba = qmol->getAtomWithIdx(ab->getBeginAtomIdx());
            QueryAtom *nqba = buildQueryAtom(qba);
            qmol->replaceAtom(ab->getBeginAtomIdx(), nqba);

            auto qea = qmol->getAtomWithIdx(ab->getEndAtomIdx());
            QueryAtom *nqea = buildQueryAtom(qea);
            qmol->replaceAtom(ab->getEndAtomIdx(), nqea);
          }
          break;
        }
      }
    }

    if (qmol) {
      frag.reset(qmol.release());
    }
  }
}

// Split the molecule into fragments.  maxBondSplits gives the maximum number
// of bonds to be used in each split.  There will a vector of vectors of
// molecules, 1 inner vector for each split i.e. maxBondSplits in total, the
// first with 1 split, the 2nd with 2 etc.  Each inner vector contains the
// fragments from a split molecule.
std::vector<std::vector<std::unique_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxBondSplits) {
  //  std::cout << "Splitting " << MolToSmiles(query) << " with " <<
  //  maxBondSplits
  //            << " bonds." << std::endl;

  if (maxBondSplits < 1) {
    maxBondSplits = 1;
  }
  if (maxBondSplits > 3) {
    maxBondSplits = 3;
  }
  if (maxBondSplits > query.getNumBonds()) {
    maxBondSplits = query.getNumBonds();
  }
  auto ringInfo = query.getRingInfo();
  boost::dynamic_bitset<> ringBonds(query.getNumBonds());
  for (const auto &r : ringInfo->bondRings()) {
    for (const auto b : r) {
      ringBonds.set(b);
    }
  }

  std::vector<std::vector<std::unique_ptr<ROMol>>> fragments;
  // Keep the molecule itself (i.e. 0 splits).  It will probably produce
  // lots of hits but one can imagine a use for it.
  fragments.push_back(std::vector<std::unique_ptr<ROMol>>());
  fragments.back().emplace_back(new ROMol(query));

  // Now do the splits.
  for (unsigned int i = 1; i <= maxBondSplits; ++i) {
    //    std::cout << "Splitting with up to " << i << " bonds" << std::endl;
    auto combs = combMFromN(i, static_cast<int>(query.getNumBonds()));
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
    for (unsigned int j = 1; j <= i; ++j) {
      dummyLabels.push_back(std::make_pair(j, j));
    }
    //    std::cout << "Number of possible splits : " << combs.size() <<
    //    std::endl;
    for (auto &c : combs) {
      //      for (auto &i : c) {
      //        std::cout << i << " ";
      //      }
      //      std::cout << std::endl;
      // don't break just 1 ring bond, as it can't create 2 fragments.  It
      // could be better than this, by checking that any number of ring
      // bonds are all in the same ring system.  Maybe look at that
      // if necessary for performance.  Triazoles can be created from 3
      // synthons, so breaking 3 ring bonds is ok.  This will still pass
      // through cases that break 2 ring bonds in separate ring systems,
      // such as a bond in each of the 2 phenyl rings in c1ccccc1c2ccccc2,
      // but they will be caught below.
      auto numRingBonds =
          std::reduce(c.begin(), c.end(), 0, [&](int prevRes, int bondNum) {
            if (ringBonds[bondNum]) {
              return prevRes + 1;
            } else {
              return prevRes;
            }
          });
      if (numRingBonds == 1) {
        continue;
      }
      std::unique_ptr<ROMol> fragMol(
          MolFragmenter::fragmentOnBonds(query, c, true, &dummyLabels));
      //      std::cout << MolToSmiles(*fragMol) << std::endl;
      std::vector<std::unique_ptr<ROMol>> molFrags;
      auto numFrags = MolOps::getMolFrags(*fragMol, molFrags, false);
      // Must have been a ring-opening.
      if (numFrags == 1) {
        continue;
      }
      if (checkConnectorsInDifferentFrags(molFrags, i)) {
        fixAromaticRingSplits(molFrags);
        fragments.emplace_back(std::move(molFrags));
      }
    }
    //    std::cout << "Number of valid splits : " << fragments.size() <<
    //    std::endl;
  }
  //  std::cout << "Fragments size : " << fragments.size() << std::endl;
  return fragments;
}

int countConnections(const std::string &smiles) {
  static const std::regex conns(R"(\[[1-4]\*\])");
  return static_cast<int>(std::distance(
      std::sregex_token_iterator(smiles.begin(), smiles.end(), conns),
      std::sregex_token_iterator()));
}
}  // namespace details

}  // namespace SynthonSpaceSearch
}  // namespace RDKit
