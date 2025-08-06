//
//  Copyright (C) 2025 David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <numeric>
#include <unordered_set>

#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

namespace RDKit {
namespace {
bool atomsMatch(const Atom *atom1, const Atom *atom2) {
  if (atom1->getAtomicNum() == atom2->getAtomicNum() &&
      atom1->getIsAromatic() == atom2->getIsAromatic()) {
    return true;
  }
  return false;
}

// Return true if the 2 pairs of atoms are both not bonded, or both
// bonded by a bond of the same type.
bool bondsMatch(const ROMol &mol1, unsigned int atom11, unsigned int atom12,
                const ROMol &mol2, unsigned int atom21, unsigned int atom22) {
  if (atom11 == atom12 || atom21 == atom22) {
    return false;
  }
  auto bond1 = mol1.getBondBetweenAtoms(atom11, atom12);
  auto bond2 = mol2.getBondBetweenAtoms(atom21, atom22);
  if (!bond1 && !bond2) {
    return true;
  }
  if (bond1 && bond2) {
    if (bond1->getBondType() == bond2->getBondType()) {
      return true;
    }
  }
  return false;
}

void buildAtomPairs(
    const ROMol &mol1, const ROMol &mol2,
    std::vector<std::pair<unsigned int, unsigned int>> &atomPairs) {
  atomPairs.reserve(mol1.getNumAtoms() * mol2.getNumAtoms());
  for (const auto at1 : mol1.atoms()) {
    for (const auto at2 : mol2.atoms()) {
      if (atomsMatch(at1, at2)) {
        atomPairs.push_back(std::make_pair(at1->getIdx(), at2->getIdx()));
      }
    }
  }
}

void buildCorrespondenceGraph(
    const std::vector<std::pair<unsigned int, unsigned int>> &atomPairs,
    const ROMol &mol1, const ROMol &mol2,
    std::vector<std::unordered_set<unsigned int>> &corrGraph) {
  const auto distMat1 = MolOps::getDistanceMat(mol1);
  const auto distMat2 = MolOps::getDistanceMat(mol2);
  for (size_t i = 0U; i < atomPairs.size() - 1; ++i) {
    for (size_t j = i + 1; j < atomPairs.size(); ++j) {
      if (!bondsMatch(mol1, atomPairs[i].first, atomPairs[j].first, mol2,
                      atomPairs[i].second, atomPairs[j].second)) {
        continue;
      }
      auto d1 = distMat1[atomPairs[i].first * mol1.getNumAtoms() +
                         atomPairs[j].first];
      auto d2 = distMat2[atomPairs[i].second * mol2.getNumAtoms() +
                         atomPairs[j].second];
      if (fabs(d1 - d2) < 1.0e-6) {
        corrGraph[i].insert(j);
        corrGraph[j].insert(i);
      }
    }
  }
}

void bron_kerbosch(
    std::unordered_set<unsigned int> &clique,
    std::unordered_set<unsigned int> &remaining,
    std::unordered_set<unsigned int> &done,
    const std::vector<std::unordered_set<unsigned int>> &corrGraph,
    std::vector<std::vector<unsigned int>> &maxCliques) {
#if 0
  std::cout << "\nNext Round" << std::endl;
  std::cout << "clique : ";
  for (auto c : clique) {
    std::cout << c << " ";
  }
  std::cout << "  remaining : ";
  for (auto c : remaining) {
    std::cout << c << " ";
  }
  std::cout << "  done : ";
  for (auto c : done) {
    std::cout << c << " ";
  }
  std::cout << std::endl;
#endif
  if (remaining.empty() && done.empty()) {
    // std::cout << "CLIQUE : " << clique.size() << " :: ";
    // for (auto c : clique) {
    //   std::cout << c << " ";
    // }
    // std::cout << std::endl;
    if (maxCliques.empty()) {
      maxCliques.push_back({clique.begin(), clique.end()});
    } else {
      if (clique.size() > maxCliques.front().size()) {
        maxCliques.clear();
      } else if (clique.size() < maxCliques.front().size()) {
        return;
      }
      maxCliques.push_back({clique.begin(), clique.end()});
    }
    return;
  }
  // select a pivot node from the union of remaining and done that
  // has maximum degree
  std::unordered_set<unsigned int> rem_union_done(remaining.begin(),
                                                  remaining.end());
  rem_union_done.insert(done.begin(), done.end());
  auto pivot = *max_element(rem_union_done.begin(), rem_union_done.end(),
                            [&](const auto &a, const auto &b) -> bool {
                              return corrGraph[a].size() < corrGraph[b].size();
                            });
  // std::cout << "pivot : " << pivot << std::endl;
  // Make a set of possible vertices to add to the clique.  These are
  // vertices in 'remaining' that are not neighbours of 'pivot'.
  std::vector<unsigned int> possibles;
  for (auto r : remaining) {
    if (corrGraph[pivot].find(r) == corrGraph[pivot].end()) {
      possibles.push_back(r);
    }
  }
#if 0
  std::sort(possibles.begin(), possibles.end());
  std::cout << "possibles : ";
  for (auto p : possibles) {
    std::cout << p << " ";
  }
  std::cout << std::endl;
#endif
  for (auto v : possibles) {
    // std::cout << "vertex : " << v << std::endl;
    // Make a new clique from the current clique + v
    std::unordered_set<unsigned int> newClique(clique.begin(), clique.end());
    newClique.insert(v);
    // The new remaining are members of remaining that are neighbours of vertex.
    std::unordered_set<unsigned int> newRemaining;
    for (auto r : remaining) {
      if (corrGraph[v].find(r) != corrGraph[v].end()) {
        newRemaining.insert(r);
      }
    }
    // Likewise, the new done are members of done that are neighbours of vertex.
    std::unordered_set<unsigned int> newDone;
    for (auto d : done) {
      if (corrGraph[v].find(d) != corrGraph[v].end()) {
        newDone.insert(d);
      }
    }

    // Step
    bron_kerbosch(newClique, newRemaining, newDone, corrGraph, maxCliques);
    // Move vertex from remaining to done
    remaining.erase(v);
    done.insert(v);
  }
}
}  // namespace

void TwoMolMCSS(const ROMol &mol1, const ROMol &mol2,
                std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
                    &maxCliques) {
  std::vector<std::pair<unsigned int, unsigned int>> atomPairs;
  buildAtomPairs(mol1, mol2, atomPairs);
  if (atomPairs.empty()) {
    return;
  }
  std::vector<std::unordered_set<unsigned int>> corrGraph(
      atomPairs.size(), std::unordered_set<unsigned int>());
  buildCorrespondenceGraph(atomPairs, mol1, mol2, corrGraph);
  // std::vector<std::unordered_set<unsigned int>> corrGraph(
  //     6, std::unordered_set<unsigned int>());
  // corrGraph[0].insert(1);
  // corrGraph[0].insert(4);
  // corrGraph[1].insert(0);
  // corrGraph[1].insert(2);
  // corrGraph[1].insert(4);
  // corrGraph[2].insert(1);
  // corrGraph[2].insert(3);
  // corrGraph[3].insert(2);
  // corrGraph[3].insert(4);
  // corrGraph[3].insert(5);
  //
  // corrGraph[4].insert(0);
  // corrGraph[4].insert(1);
  // corrGraph[4].insert(3);
  //
  // corrGraph[5].insert(3);

  std::unordered_set<unsigned int> clique;
  std::unordered_set<unsigned int> remaining(corrGraph.size());
  for (unsigned int i = 0U; i < corrGraph.size(); ++i) {
    remaining.insert(i);
  }
  std::unordered_set<unsigned int> done;
  std::vector<std::vector<unsigned int>> rawMaxCliques;
  bron_kerbosch(clique, remaining, done, corrGraph, rawMaxCliques);
  maxCliques.reserve(rawMaxCliques.size());
  for (const auto &rmc : rawMaxCliques) {
    std::vector<std::pair<unsigned int, unsigned int>> newClique;
    newClique.reserve(rmc.size());
    for (const auto cm : rmc) {
      newClique.push_back(atomPairs[cm]);
    }
    std::ranges::sort(newClique);
    maxCliques.push_back(std::move(newClique));
  }
}

std::string makeSMARTSFromMCSS(
    const ROMol &mol1,
    const std::vector<std::pair<unsigned int, unsigned int>> &clique) {
  RWMol qmol(mol1);
  std::unordered_set<unsigned int> atomsToKeep;
  std::transform(clique.begin(), clique.end(),
                 std::inserter(atomsToKeep, atomsToKeep.begin()),
                 [](const auto &a) -> unsigned int { return a.first; });
  qmol.beginBatchEdit();
  for (auto atom : qmol.atoms()) {
    if (atomsToKeep.find(atom->getIdx()) == atomsToKeep.end()) {
      qmol.removeAtom(atom);
    }
  }
  qmol.commitBatchEdit();
  return MolToSmarts(qmol);
}
}  // namespace RDKit