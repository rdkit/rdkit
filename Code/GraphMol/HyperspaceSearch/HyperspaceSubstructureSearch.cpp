//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file contains an implementation of hyperspace substructure search
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
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "Hyperspace.h"

namespace RDKit {
namespace HyperspaceSSSearch {
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

// Split the molecule into fragments.  maxBondSplits gives the maximum number
// of bonds to be used in each split.  There will 1 vector of molecules
// for each split i.e. maxBondSplits in total, the first with 1 split, the 2nd
// with 2 etc.  Each ROMol contains a split molecule containing the fragments
// after the split (the SMILES will be dot-connected with isotope-labelled
// dummies showing the split bonds).
std::vector<std::vector<std::unique_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxBondSplits) {
  std::cout << "Splitting " << MolToSmiles(query) << " with " << maxBondSplits
            << " bonds." << std::endl;

  if (maxBondSplits < 1) {
    maxBondSplits = 1;
  }
  if (maxBondSplits > 3) {
    maxBondSplits = 3;
  }
  std::vector<std::vector<std::unique_ptr<ROMol>>> fragments;
  auto ringInfo = query.getRingInfo();
  boost::dynamic_bitset<> ringBonds(query.getNumBonds());
  for (const auto &r : ringInfo->bondRings()) {
    for (const auto b : r) {
      ringBonds.set(b);
    }
  }
  std::vector<int> fragAtoms;

  // Keep the molecule itself (i.e. 0 splits).  It will probably produce
  // lots of hits but one can imagine a use for it.
  fragments.push_back(std::vector<std::unique_ptr<ROMol>>());
  fragments.back().emplace_back(new ROMol(query));

  // Now do the splits.
  for (unsigned int i = 1; i <= maxBondSplits; ++i) {
    std::cout << "Splitting with up to " << i << " bonds" << std::endl;
    auto combs = combMFromN(i, static_cast<int>(query.getNumBonds()));
    fragments.push_back(std::vector<std::unique_ptr<ROMol>>());
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
    for (unsigned int j = 1; j <= i; ++j) {
      dummyLabels.push_back(std::make_pair(j, j));
    }
    std::cout << "Number of possible splits : " << combs.size() << std::endl;
    for (auto &c : combs) {
      std::unique_ptr<ROMol> fragMol(
          MolFragmenter::fragmentOnBonds(query, c, true, &dummyLabels));
      auto numFrags = MolOps::getMolFrags(*fragMol, fragAtoms);
      if (numFrags == 1) {
        continue;
      }
      // the fragmentation is valid if the 2 ends of each bond are in different
      // fragments.
      bool fragsOk = true;
      // Loop over the isotope numbers of the ends of the splits
      for (unsigned int j = 1; j <= i; ++j) {
        // Loop over the different fragments
        for (unsigned int k = 0; k < j; ++k) {
          // Count the number of dummy atoms of the correct isotope
          // number in this fragment.
          int dummyCount = 0;
          for (unsigned int l = 0; l < fragMol->getNumAtoms(); ++l) {
            if (static_cast<unsigned int>(fragAtoms[l]) == k) {
              auto atom = fragMol->getAtomWithIdx(l);
              if (!atom->getAtomicNum() && atom->getIsotope() == j) {
                ++dummyCount;
              }
            }
          }
          // If it's more than 1 dummy of this isotope number in this fragment
          // it's a bad fragmentation, probably because a single bond in a
          // ring has been cut.
          if (dummyCount > 1) {
            fragsOk = false;
          }
        }
        if (!fragsOk) {
          break;
        }
      }
      if (fragsOk) {
        fragments.back().push_back(std::move(fragMol));
      }
    }
    std::cout << "Number of valid splits : " << fragments.back().size()
              << std::endl;
  }
  return fragments;
}

}  // namespace details

// Do a substructure search for query in the hyperspace.
std::vector<std::unique_ptr<ROMol>> SSSearch(const ROMol &query,
                                             unsigned int maxBondSplits,
                                             Hyperspace &hyperspace) {
  auto &reactions = hyperspace.reactions();
  for (const auto &r : reactions) {
    std::cout << "reaction " << r.first
              << "  number of synthons : " << r.second->d_reagents.size()
              << std::endl;
  }
  auto results = hyperspace.search(query, maxBondSplits);
  return results;
}

// overload taking the name of the hyperspace file.
std::vector<std::unique_ptr<ROMol>> SSSearch(const ROMol &query,
                                             unsigned int maxBondSplits,
                                             const std::string &libName) {
  std::cout << "Searching library " << libName << " for structures containing "
            << MolToSmiles(query) << std::endl;
  Hyperspace hyperspace(libName);
  std::cout << "Number of reactions : " << hyperspace.numReactions()
            << std::endl;
  auto &reactions = hyperspace.reactions();
  for (const auto &r : reactions) {
    std::cout << "reaction " << r.first
              << "  number of synthons : " << r.second->d_reagents.size()
              << std::endl;
  }
  auto results = hyperspace.search(query, maxBondSplits);
  return results;
}

}  // namespace HyperspaceSSSearch
}  // namespace RDKit
