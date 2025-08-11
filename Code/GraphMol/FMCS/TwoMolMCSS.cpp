//
//  Copyright (C) 2025 David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "GraphMol/CIPLabeler/Descriptor.h"

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

// Returns 0 if one bond is nullptr.
// Returns 'c' if they are the same type (for the cEdges) or
// 'd' if they aren't or are both nullptr (for the dEdges).
char bondsMatch(const Bond *bond1, const Bond *bond2) {
  if (!bond1 && !bond2) {
    return 'd';
  }
  if (bond1 && bond2) {
    if (bond1->getBondType() == bond2->getBondType()) {
      return 'c';
    }
    return 'd';
  }
  return 0;
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
    std::vector<std::unordered_set<unsigned int>> &corrGraph,
    std::vector<std::unordered_set<unsigned int>> &cEdges,
    std::vector<std::unordered_set<unsigned int>> &dEdges) {
  std::vector<std::vector<const Bond *>> bond1s(
      mol1.getNumAtoms(),
      std::vector<const Bond *>(mol1.getNumAtoms(), nullptr));
  for (const auto b : mol1.bonds()) {
    bond1s[b->getBeginAtomIdx()][b->getEndAtomIdx()] = b;
    bond1s[b->getEndAtomIdx()][b->getBeginAtomIdx()] = b;
  }
  std::vector<std::vector<const Bond *>> bond2s(
      mol2.getNumAtoms(),
      std::vector<const Bond *>(mol2.getNumAtoms(), nullptr));
  for (const auto b : mol2.bonds()) {
    bond2s[b->getBeginAtomIdx()][b->getEndAtomIdx()] = b;
    bond2s[b->getEndAtomIdx()][b->getBeginAtomIdx()] = b;
  }
  for (size_t i = 0U; i < atomPairs.size() - 1; ++i) {
    for (size_t j = i + 1; j < atomPairs.size(); ++j) {
      if (atomPairs[i].first == atomPairs[j].first ||
          atomPairs[i].second == atomPairs[j].second) {
        continue;
      }
      auto bond1 = bond1s[atomPairs[i].first][atomPairs[j].first];
      auto bond2 = bond2s[atomPairs[i].second][atomPairs[j].second];
      auto bm = bondsMatch(bond1, bond2);
      if (!bm) {
        continue;
      }
      corrGraph[i].insert(j);
      corrGraph[j].insert(i);
      if (bm == 'c') {
        cEdges[i].insert(j);
        cEdges[j].insert(i);
      } else if (bm == 'd') {
        dEdges[i].insert(j);
        dEdges[j].insert(i);
      }
#if 0
      auto d1 = distMat1[atomPairs[i].first * mol1.getNumAtoms() +
                         atomPairs[j].first];
      auto d2 = distMat2[atomPairs[i].second * mol2.getNumAtoms() +
                         atomPairs[j].second];
      if (fabs(d1 - d2) < 1.0e-6) {
        corrGraph[i].insert(j);
        corrGraph[j].insert(i);
      }
#endif
    }
  }
}
void buildCorrespondenceGraph(
    const std::vector<std::pair<unsigned int, unsigned int>> &atomPairs,
    const ROMol &mol1, const ROMol &mol2,
    std::vector<std::vector<char>> &corrGraph) {
  std::vector<std::vector<const Bond *>> bond1s(
      mol1.getNumAtoms(),
      std::vector<const Bond *>(mol1.getNumAtoms(), nullptr));
  for (const auto b : mol1.bonds()) {
    bond1s[b->getBeginAtomIdx()][b->getEndAtomIdx()] = b;
    bond1s[b->getEndAtomIdx()][b->getBeginAtomIdx()] = b;
  }
  std::vector<std::vector<const Bond *>> bond2s(
      mol2.getNumAtoms(),
      std::vector<const Bond *>(mol2.getNumAtoms(), nullptr));
  for (const auto b : mol2.bonds()) {
    bond2s[b->getBeginAtomIdx()][b->getEndAtomIdx()] = b;
    bond2s[b->getEndAtomIdx()][b->getBeginAtomIdx()] = b;
  }
  for (size_t i = 0U; i < atomPairs.size() - 1; ++i) {
    for (size_t j = i + 1; j < atomPairs.size(); ++j) {
      if (atomPairs[i].first == atomPairs[j].first ||
          atomPairs[i].second == atomPairs[j].second) {
        continue;
      }
      auto bond1 = bond1s[atomPairs[i].first][atomPairs[j].first];
      auto bond2 = bond2s[atomPairs[i].second][atomPairs[j].second];
      auto bm = bondsMatch(bond1, bond2);
      corrGraph[i][j] = bm;
      corrGraph[j][i] = bm;
#if 0
      auto d1 = distMat1[atomPairs[i].first * mol1.getNumAtoms() +
                         atomPairs[j].first];
      auto d2 = distMat2[atomPairs[i].second * mol2.getNumAtoms() +
                         atomPairs[j].second];
      if (fabs(d1 - d2) < 1.0e-6) {
        corrGraph[i].insert(j);
        corrGraph[j].insert(i);
      }
#endif
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
    // Step if we could get something larger than we have already.
    if (maxCliques.empty() ||
        newRemaining.size() + newClique.size() > maxCliques.front().size()) {
      bron_kerbosch(newClique, newRemaining, newDone, corrGraph, maxCliques);
    }
    // Move vertex from remaining to done
    remaining.erase(v);
    done.insert(v);
  }
}

// Notation from Koch, Theoretical Computer Science 250 (2001) 1–30
// Algorithm 3
// c is the current clique
// p is the set of nodes that can be added to c. They are neighbours of
// the last node in c (u) via c-edges.
// d is the set of nodes that cannot be added to c because they are
// neighbours of u via d-edges.
// s is the set of nodes already tried as u in this step in the
// recursion.
// t is the set of nodes that have already been used to initialise a
// clique or have appeared in a clique already
void enumerate_c_cliques(
    std::vector<unsigned int> &c, const std::unordered_set<unsigned int> &p,
    const std::unordered_set<unsigned int> &d,
    std::unordered_set<unsigned int> &s, std::unordered_set<unsigned int> &t,
    const std::vector<std::unordered_set<unsigned int>> &corrGraph,
    const std::vector<std::unordered_set<unsigned int>> &cEdges,
    const std::vector<std::unordered_set<unsigned int>> &dEdges,
    std::vector<std::vector<unsigned int>> &maxCliques) {
#if 0
  std::cout << "c ";
  for (auto cm : c) {
    std::cout << cm << " ";
  }
  std::cout << ":: p ";
  for (auto pm : p) {
    std::cout << pm << " ";
  }
  std::cout << ":: d ";
  for (auto pm : d) {
    std::cout << pm << " ";
  }
  std::cout << ":: s ";
  for (auto pm : s) {
    std::cout << pm << " ";
  }
  std::cout << ":: t ";
  for (auto pm : t) {
    std::cout << pm << " ";
  }
  std::cout << std::endl;
#endif
  if (!maxCliques.empty() && c.size() + corrGraph.size() - p.size() - s.size() <
                                 maxCliques.front().size()) {
    return;
  }
  if (p.empty() && s.empty() && c.size() > 1) {
    // Doing this speeds up the search enormously, because it means that
    // we won't start a new clique with a node that has already been in a
    // clique.  It does mean that symmetrically equivalent cliques won't
    // be returned.
    t.insert(c.begin(), c.end());
    // std::cout << "CLIQUE : " << c.size() << " :: ";
    // for (auto cm : c) {
    //   std::cout << cm << " ";
    // }
    // std::cout << std::endl;
    if (maxCliques.empty()) {
      // std::cout << "New MAX CLIQUE : " << c.size() << " :: ";
      maxCliques.push_back({c.begin(), c.end()});
    } else {
      if (c.size() > maxCliques.front().size()) {
        maxCliques.clear();
      } else if (c.size() < maxCliques.front().size()) {
        return;
      }
      // std::cout << "New MAX CLIQUE : " << c.size() << " :: ";
      // std::cout << std::endl;
      maxCliques.push_back(c);
    }
    return;
  }
  for (auto u : p) {
    if (s.contains(u)) {
      continue;
    }
    // Form a new P which is the current p minus this u.
    std::unordered_set<unsigned int> newP(p);
    newP.erase(u);
    // newD and newS start out as copies of the incoming.
    std::unordered_set<unsigned int> newD(d);
    std::unordered_set<unsigned int> newS(s);
    // n is the set of all neighbours of u.
    const auto &n = corrGraph[u];
    // Now loop over all nodes that aren't connected by a d-edge to the current
    // clique and include in newP all of those that are neighbours of u
    // via a c-edge and take then out of newD.
    for (auto v : d) {
      if (cEdges[u].contains(v)) {
        // If v is in t, we've already tried it as a root of the clique,
        // so any clique that contains it has already been found.  Therefore,
        // put it in s, not to be considered further.
        if (t.contains(v)) {
          newS.insert(v);
        } else {
          newP.insert(v);
        }
        newD.erase(v);
      }
    }
    // the new clique is the current clique plus u.
    std::vector<unsigned int> newC(c);
    newC.push_back(u);

    // We now recurse by using only the members of newP, newD and newS
    // that are neighbours of u
    std::unordered_set<unsigned int> onwardP;
    for (auto pe : newP) {
      if (n.contains(pe)) {
        onwardP.insert(pe);
      }
    }
    std::unordered_set<unsigned int> onwardD;
    for (auto de : newD) {
      if (n.contains(de)) {
        onwardD.insert(de);
      }
    }
    std::unordered_set<unsigned int> onwardS;
    for (auto se : newS) {
      if (n.contains(se)) {
        onwardS.insert(se);
      }
    }
    enumerate_c_cliques(newC, onwardP, onwardD, onwardS, t, corrGraph, cEdges,
                        dEdges, maxCliques);
    s.insert(u);
  }
}

void enumerate_z_cliques(
    std::vector<unsigned int> &c, const std::unordered_set<unsigned int> &p,
    const std::unordered_set<unsigned int> &d,
    std::unordered_set<unsigned int> &s, std::unordered_set<unsigned int> &t,
    const std::vector<std::unordered_set<unsigned int>> &corrGraph,
    const std::vector<std::unordered_set<unsigned int>> &cEdges,
    const std::vector<std::unordered_set<unsigned int>> &dEdges,
    std::vector<std::vector<unsigned int>> &maxCliques) {
#if 0
  std::cout << "c ";
  for (auto cm : c) {
    std::cout << cm << " ";
  }
  std::cout << ":: p ";
  for (auto pm : p) {
    std::cout << pm << " ";
  }
  std::cout << ":: d ";
  for (auto pm : d) {
    std::cout << pm << " ";
  }
  std::cout << ":: s ";
  for (auto pm : s) {
    std::cout << pm << " ";
  }
  std::cout << ":: t ";
  for (auto pm : t) {
    std::cout << pm << " ";
  }
  std::cout << std::endl;
#endif
  if (p.empty()) {
    if (s.empty() && c.size() > 1) {
      // Doing this speeds up the search enormously, because it means that
      // we won't start a new clique with a node that has already been in a
      // clique.  It does mean that symmetrically equivalent cliques won't
      // be returned.
      t.insert(c.begin(), c.end());
      // std::cout << "CLIQUE : " << c.size() << " :: ";
      // for (auto cm : c) {
      //   std::cout << cm << " ";
      // }
      // std::cout << std::endl;
      if (maxCliques.empty()) {
        // std::cout << "New MAX CLIQUE : " << c.size() << " :: ";
        maxCliques.push_back({c.begin(), c.end()});
      } else {
        if (c.size() > maxCliques.front().size()) {
          maxCliques.clear();
        } else if (c.size() < maxCliques.front().size()) {
          return;
        }
        // std::cout << "New MAX CLIQUE : " << c.size() << " :: ";
        // std::cout << std::endl;
        maxCliques.push_back(c);
      }
    }
    return;
  }
  auto ut = *max_element(p.begin(), p.end(),
                         [&](const auto &a, const auto &b) -> bool {
                           return corrGraph[a].size() < corrGraph[b].size();
                         });
  for (auto ui : p) {
    // if ui is adjacent to ut
    // OR if ui is connected via a c-path to a vertex in D that is
    // not adjacent to ut
    bool ok1 = false;
    bool ok2 = !corrGraph[ui].contains(ut);
    if (!ok2) {
      for (auto de : d) {
        if (cEdges[ui].contains(de) && !corrGraph[de].contains(ut)) {
          ok1 = true;
          break;
        }
      }
    }
    if (ok2 || ok1) {
      // Form a new P which is the current p minus this ui.
      std::unordered_set<unsigned int> newP(p);
      newP.erase(ui);
      // newD and newS start out as copies of the incoming.
      std::unordered_set<unsigned int> newD(d);
      std::unordered_set<unsigned int> newS(s);
      // n is the set of all neighbours of u.
      const auto &n = corrGraph[ui];
      for (auto v : d) {
        if (p.contains(v)) {
          newP.insert(v);
        } else if (d.contains(v)) {
          // can v be added to P?
          if (cEdges[ui].contains(v)) {
            if (t.contains(v)) {
              newS.insert(v);
            } else {
              newP.insert(v);
            }
            newD.erase(v);
          } else if (s.contains(v)) {
            newS.erase(v);
          }
        }
      }
      // the new clique is the current clique plus u.
      std::vector<unsigned int> newC(c);
      newC.push_back(ui);

      // We now recurse by using only the members of newP, newD and newS
      // that are neighbours of u
      std::unordered_set<unsigned int> onwardP;
      for (auto pe : newP) {
        if (n.contains(pe)) {
          onwardP.insert(pe);
        }
      }
      std::unordered_set<unsigned int> onwardD;
      for (auto de : newD) {
        if (n.contains(de)) {
          onwardD.insert(de);
        }
      }
      std::unordered_set<unsigned int> onwardS;
      for (auto se : newS) {
        if (n.contains(se)) {
          onwardS.insert(se);
        }
      }
      enumerate_z_cliques(newC, onwardP, onwardD, onwardS, t, corrGraph, cEdges,
                          dEdges, maxCliques);
      s.insert(ui);
    }
  }
}

void enumerate_z_cliques(std::vector<unsigned int> &c,
                         const std::unordered_set<unsigned int> &p,
                         const std::unordered_set<unsigned int> &d,
                         std::unordered_set<unsigned int> &s,
                         std::unordered_set<unsigned int> &t,
                         const std::vector<std::vector<char>> &corrGraph,
                         std::vector<std::vector<unsigned int>> &maxCliques) {
#if 0
  std::cout << "c ";
  for (auto cm : c) {
    std::cout << cm << " ";
  }
  std::cout << ":: p ";
  for (auto pm : p) {
    std::cout << pm << " ";
  }
  std::cout << ":: d ";
  for (auto pm : d) {
    std::cout << pm << " ";
  }
  std::cout << ":: s ";
  for (auto pm : s) {
    std::cout << pm << " ";
  }
  std::cout << ":: t ";
  for (auto pm : t) {
    std::cout << pm << " ";
  }
  std::cout << std::endl;
#endif
  if (p.empty()) {
    if (s.empty() && c.size() > 1) {
      // Doing this speeds up the search enormously, because it means that
      // we won't start a new clique with a node that has already been in a
      // clique.  It does mean that symmetrically equivalent cliques won't
      // be returned.
      t.insert(c.begin(), c.end());
      // std::cout << "CLIQUE : " << c.size() << " :: ";
      // for (auto cm : c) {
      //   std::cout << cm << " ";
      // }
      // std::cout << std::endl;
      if (maxCliques.empty()) {
        // std::cout << "New MAX CLIQUE : " << c.size() << " :: ";
        maxCliques.push_back({c.begin(), c.end()});
      } else {
        if (c.size() > maxCliques.front().size()) {
          maxCliques.clear();
        } else if (c.size() < maxCliques.front().size()) {
          return;
        }
        // std::cout << "New MAX CLIQUE : " << c.size() << " :: ";
        // std::cout << std::endl;
        maxCliques.push_back(c);
      }
    }
    return;
  }
  auto ut = *max_element(p.begin(), p.end(),
                         [&](const auto &a, const auto &b) -> bool {
                           return corrGraph[a].size() < corrGraph[b].size();
                         });
  for (auto ui : p) {
    // if ui is adjacent to ut
    // OR if ui is connected via a c-path to a vertex in D that is
    // not adjacent to ut
    bool ok1 = false;
    bool ok2 = !corrGraph[ui][ut];
    if (!ok2) {
      for (auto de : d) {
        if (!corrGraph[de][ut] && corrGraph[ui][de] == 'c') {
          ok1 = true;
          break;
        }
      }
    }
    if (ok2 || ok1) {
      // Form a new P which is the current p minus this ui.
      std::unordered_set<unsigned int> newP(p);
      newP.erase(ui);
      // newD and newS start out as copies of the incoming.
      std::unordered_set<unsigned int> newD(d);
      std::unordered_set<unsigned int> newS(s);
      // n is the set of all neighbours of u.
      const auto &n = corrGraph[ui];
      for (auto v : d) {
        if (p.contains(v)) {
          newP.insert(v);
        } else if (d.contains(v)) {
          // can v be added to P?
          if (corrGraph[ui][v] == 'c') {
            if (t.contains(v)) {
              newS.insert(v);
            } else {
              newP.insert(v);
            }
            newD.erase(v);
          } else if (s.contains(v)) {
            newS.erase(v);
          }
        }
      }
      // the new clique is the current clique plus u.
      std::vector<unsigned int> newC(c);
      newC.push_back(ui);

      // We now recurse by using only the members of newP, newD and newS
      // that are neighbours of u
      std::unordered_set<unsigned int> onwardP;
      for (auto pe : newP) {
        if (n[pe]) {
          onwardP.insert(pe);
        }
      }
      std::unordered_set<unsigned int> onwardD;
      for (auto de : newD) {
        if (n[de]) {
          onwardD.insert(de);
        }
      }
      std::unordered_set<unsigned int> onwardS;
      for (auto se : newS) {
        if (n[se]) {
          onwardS.insert(se);
        }
      }
      enumerate_z_cliques(newC, onwardP, onwardD, onwardS, t, corrGraph,
                          maxCliques);
      s.insert(ui);
    }
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
#if 0
  int i = 0;
  for (const auto &ap : atomPairs) {
    std::cout << i++ << " : " << ap.first << " -> " << ap.second << std::endl;
  }
#endif
  std::vector<std::vector<char>> corrGraph(
      atomPairs.size(), std::vector<char>(atomPairs.size(), 0));
  buildCorrespondenceGraph(atomPairs, mol1, mol2, corrGraph);

#if 0
  std::vector<std::unordered_set<unsigned int>> oldCorrGraph(
      atomPairs.size(), std::unordered_set<unsigned int>());
  // cEdges are edges in the correspondence graph where the 2 atoms
  // in each pair are bonded
  std::vector<std::unordered_set<unsigned int>> cEdges(
      atomPairs.size(), std::unordered_set<unsigned int>());
  // cEdges are the opposite - atoms in both pairs are not bonded.
  std::vector<std::unordered_set<unsigned int>> dEdges(
      atomPairs.size(), std::unordered_set<unsigned int>());
  buildCorrespondenceGraph(atomPairs, mol1, mol2, oldCorrGraph, cEdges, dEdges);
  auto printGraph =
      [](const std::vector<std::unordered_set<unsigned int>> &g) -> void {
    unsigned int i = 0;
    for (const auto &l : g) {
      std::cout << i++ << " :: " << l.size() << " :: ";
      for (const auto &u : l) {
        std::cout << u << " ";
      }
      std::cout << std::endl;
    }
  };
  std::cout << "Corr Graph" << std::endl;
  printGraph(corrGraph);
  std::cout << "cEdges" << std::endl;
  printGraph(cEdges);
  std::cout << "dEdges" << std::endl;
  printGraph(dEdges);
#endif
  std::vector<unsigned int> clique;
  std::vector<std::vector<unsigned int>> rawMaxCliques;
#if 0
  std::unordered_set<unsigned int> remaining(corrGraph.size());
  for (unsigned int i = 0U; i < corrGraph.size(); ++i) {
    remaining.insert(i);
  }
  std::unordered_set<unsigned int> done;
  bron_kerbosch(clique, remaining, done, corrGraph, rawMaxCliques);
#else
  // Running enumerate_z_cliques over each node in the corrGraph in turn.
  std::unordered_set<unsigned int> t;
  std::unordered_set<unsigned int> p;
  std::unordered_set<unsigned int> d;
  std::unordered_set<unsigned int> s;
  for (size_t u = 0U; u < corrGraph.size(); ++u) {
    // std::cout << "starting with " << u << std::endl;
    p.clear();
    d.clear();
    s.clear();
    for (size_t v = 0; v < corrGraph[u].size(); ++v) {
      if (corrGraph[v][u] == 'c') {
        if (t.contains(v)) {
          s.insert(v);
        } else {
          p.insert(v);
        }
      } else if (corrGraph[v][u] == 'd') {
        d.insert(v);
      }
    }
    clique.clear();
    clique.push_back(u);
    enumerate_z_cliques(clique, p, d, s, t, corrGraph, rawMaxCliques);
    t.insert(u);
  }
#endif
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