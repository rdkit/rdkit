//
//  Copyright (C) 2025 David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "TwoMolMCSS.h"

#include "GraphMol/CIPLabeler/Descriptor.h"

#include <algorithm>
#include <numeric>
#include <unordered_set>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/FMCS/MatchTable.h>
#include <GraphMol/FMCS/chatgpt.h>
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

// Returns 0 if one bond is nullptr and the other isn't
// Returns 'c' if they are the same type (for the c-Edges) or
// 'd' if they aren't or are both nullptr (for the d-Edges).
char bondsMatch(const Bond *bond1, const Bond *bond2,
                const FMCS::MatchTable &bondMatchTable) {
  // If they're both not bonded or 1 is bonded and the other isn't,
  // they're d-Edges.
  if ((!bond1 && !bond2) || (bond1 && !bond2) || (!bond1 && bond2)) {
    return 'd';
  }
  // If they're both bonded and the same type, they're c-Edges.
  // If they're both bonded and not the same type, they're not a match
  // so return 0.
  if (bond1 && bond2) {
    if (bondMatchTable.at(bond1->getIdx(), bond2->getIdx())) {
      return 'c';
    }
    return 0;
  }
  // Default is not a match.
  return 0;
}

void buildAtomPairs(
    const ROMol &mol1, const ROMol &mol2,
    const FMCS::MatchTable &atomMatchTable,
    std::vector<std::pair<unsigned int, unsigned int>> &atomPairs) {
  atomPairs.reserve(mol1.getNumAtoms() * mol2.getNumAtoms());
  for (const auto at1 : mol1.atoms()) {
    for (const auto at2 : mol2.atoms()) {
      if (atomMatchTable.at(at1->getIdx(), at2->getIdx())) {
        atomPairs.push_back(std::make_pair(at1->getIdx(), at2->getIdx()));
      }
    }
  }
}

void buildCorrespondenceGraph(
    const std::vector<std::pair<unsigned int, unsigned int>> &atomPairs,
    const ROMol &mol1, const ROMol &mol2,
    const FMCS::MatchTable &bondMatchTable,
    std::vector<std::vector<char>> &corrGraph,
    std::vector<unsigned int> &corrGraphNumConns,
    std::vector<std::pair<unsigned int, unsigned int>> &atomPairStarts,
    std::vector<std::vector<unsigned int>> &cEdges) {
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

  for (size_t i = 1; i < atomPairs.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      if (atomPairs[i].first == atomPairs[j].first ||
          atomPairs[i].second == atomPairs[j].second) {
        continue;
      }
      auto bond1 = bond1s[atomPairs[i].first][atomPairs[j].first];
      auto bond2 = bond2s[atomPairs[i].second][atomPairs[j].second];
      auto bm = bondsMatch(bond1, bond2, bondMatchTable);
      corrGraph[i][j] = bm;
      corrGraph[j][i] = bm;
      if (bm) {
        corrGraphNumConns[i]++;
        corrGraphNumConns[j]++;
      }
      if (bm == 'c') {
        cEdges[i].push_back(j);
        cEdges[j].push_back(i);
      }
    }
  }
  atomPairStarts.resize(mol1.getNumAtoms());
  for (size_t i = 0U; i < atomPairs.size(); ++i) {
    atomPairStarts[atomPairs[i].first].first++;
  }
  for (size_t i = 1U; i < atomPairStarts.size(); ++i) {
    atomPairStarts[i].second =
        atomPairStarts[i - 1].first + atomPairStarts[i - 1].second;
  }
  atomPairStarts.erase(
      std::remove_if(atomPairStarts.begin(), atomPairStarts.end(),
                     [](const auto &a) -> bool { return a.first == 0; }),
      atomPairStarts.end());
  std::ranges::sort(atomPairStarts, [](const auto &a, const auto &b) -> bool {
    return a.first < b.first;
  });
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
  auto pivot =
      *std::max_element(rem_union_done.begin(), rem_union_done.end(),
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

bool inline inVector(const std::vector<unsigned int> &v, unsigned int val) {
  return std::find(v.begin(), v.end(), val) != v.end();
}

// Returns true if there's a path along c-Edges from ui to a vertex in d that
// isn't connected to ut.
bool cPathToNonUtNbor(unsigned int ui, unsigned int ut,
                      const std::vector<unsigned int> &d,
                      const std::vector<std::vector<char>> &corrGraph,
                      const std::vector<std::vector<unsigned int>> &cEdges,
                      std::unique_ptr<boost::dynamic_bitset<>> &inD) {
  const auto &nUt = corrGraph[ut];
  boost::dynamic_bitset<> tried(corrGraph.size());
  if (!inD) {
    inD.reset(new boost::dynamic_bitset<>(corrGraph.size()));
    for (auto dm : d) {
      (*inD)[dm] = true;
    }
  }
  for (auto dm : d) {
    if (nUt[dm]) {
      continue;
    }
    // Is there a c-path in the corrGraph from ui to dm?
    // If they're directly connected it's easy
    if (corrGraph[ui][dm] == 'c') {
      // std::cout << "return true" << std::endl;
      return true;
    }

    // Otherwise it's a breadth-first search.
    std::queue<unsigned int> q;
    q.push(ui);

    tried.reset();
    tried[ui] = true;
    while (!q.empty()) {
      auto x = q.front();
      q.pop();
      for (auto xn : cEdges[x]) {
        if (corrGraph[x][xn] != 'c' || tried[xn] || !(*inD)[xn]) {
          continue;
        }
        // It's in the path and in d, so if it's not a n'bour of ut
        // return true.
        if (!nUt[xn]) {
          // std::cout << "return true" << std::endl;
          return true;
        }
        q.push(xn);
        tried[xn] = true;
      }
    }
  }
  // std::cout << "return false" << std::endl;
  return false;
}

void enumerate_z_cliques(std::vector<unsigned int> &c,
                         const std::vector<unsigned int> p,
                         const std::vector<unsigned int> &d,
                         std::vector<unsigned int> &s,
                         const std::vector<char> &t,
                         const std::vector<std::vector<char>> &corrGraph,
                         const std::vector<unsigned int> &corrGraphNumConns,
                         const std::vector<std::vector<unsigned int>> &cEdges,
                         std::vector<std::vector<unsigned int>> &maxCliques) {
#if 0
  if (c.size() == 1) {
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
      std::cout << int(pm) << " ";
    }
    std::cout << std::endl;
  }
#endif
#if 0
  if (c.size() == 1) {
    std::cout << "NEW : " << c[0] << std::endl;
  }
  std::cout << "c : " << c.size() << " p : " << p.size() << " d : " << d.size()
            << " s : " << s.size() << " t : " << t.size() << " : "
            << "cg : " << corrGraph.size() << " : ";
  if (maxCliques.empty()) {
    std::cout << " NONE" << std::endl;
  } else {
    std::cout << "max clique : " << maxCliques.front().size() << std::endl;
  }
#endif
  if (p.empty()) {
    if (s.empty() && c.size() > 1) {
      // std::cout << "CLIQUE : " << c.size() << " :: ";
      // for (auto cm : c) {
      //   std::cout << cm << " ";
      // }
      // std::cout << std::endl;
      if (maxCliques.empty()) {
        // std::cout << "New MAX CLIQUE : " << c.size() << std::endl;
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
  auto ut = *std::max_element(
      p.begin(), p.end(), [&](const auto &a, const auto &b) -> bool {
        return corrGraphNumConns[a] < corrGraphNumConns[b];
      });
  // std::cout << "pivot node " << ut << std::endl;
  std::unique_ptr<boost::dynamic_bitset<>> inD;
  for (auto ui : p) {
    if (inVector(s, ui)) {
      continue;
    }
    // if ui is adjacent to ut
    // OR if ui is connected via a c-path to a vertex in D that is
    // not adjacent to ut
    bool ok1 = false;
    bool ok2 = !corrGraph[ui][ut];
    if (!ok2) {
      ok1 = cPathToNonUtNbor(ui, ut, d, corrGraph, cEdges, inD);
    }
    if (ok2 || ok1) {
      const auto &n = corrGraph[ui];
      // Form a new P which is the members of the current p that are
      // neighbours of ui
      std::vector<unsigned int> newP;
      for (auto pe : p) {
        if (n[pe]) {
          newP.push_back(pe);
        }
      }
      // std::cout << "init newP : " << newP.size() << " : ";
      // for (auto pe : newP) {
      //   std::cout << pe << " ";
      // }
      // std::cout << std::endl;
      // and the same for newS
      std::vector<unsigned int> newS;
      for (auto se : s) {
        if (n[se]) {
          newS.push_back(se);
        }
      }
      // The newP, newS and newD should only contain neighbours of ui.
      std::vector<unsigned int> newD;
      newD.reserve(newD.size());
      for (auto v : d) {
        bool addD = true;
        if (inVector(p, v) && n[v]) {
          // std::cout << "add " << v << " to newP 1" << std::endl;
          newP.push_back(v);
        } else {
          // can v be added to P?
          if (n[v] == 'c') {
            if (t[v] and n[v]) {
              newS.push_back(v);
            } else if (n[v]) {
              // std::cout << "add " << v << " to newP 2" << std::endl;
              newP.push_back(v);
            }
            addD = false;
          } else if (n[v] && inVector(s, v)) {
            newS.push_back(v);
          }
        }
        if (n[v] && addD) {
          newD.push_back(v);
        }
      }
      // std::cout << "inter newP : " << newP.size() << " : ";
      // for (auto pe : newP) {
      //   std::cout << pe << " ";
      // }
      // std::cout << std::endl;
      // the new clique is the current clique plus u.
      // std::vector<unsigned int> newC(c);
      // newC.push_back(ui);
      c.push_back(ui);
      enumerate_z_cliques(c, newP, newD, newS, t, corrGraph, corrGraphNumConns,
                          cEdges, maxCliques);
      c.pop_back();
      s.push_back(ui);
    }
  }
}
}  // namespace

void TwoMolMCSS(const ROMol &mol1, const ROMol &mol2, unsigned int minMCSSSize,
                const FMCS::MatchTable &atomMatchTable,
                const FMCS::MatchTable &bondMatchTable, bool uniquify,
                std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
                    &maxCliques) {
  std::vector<std::pair<unsigned int, unsigned int>> atomPairs;
  // std::cout << "minMCSSSize = " << minMCSSSize << std::endl;
  if (minMCSSSize > std::min(mol1.getNumAtoms(), mol2.getNumAtoms())) {
    return;
  }
  buildAtomPairs(mol1, mol2, atomMatchTable, atomPairs);
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
  std::vector<unsigned int> corrGraphNumConns(atomPairs.size(), 0);
  std::vector<std::pair<unsigned int, unsigned int>> atomPairStarts;
  std::vector<std::vector<unsigned int>> cEdges(atomPairs.size());
  buildCorrespondenceGraph(atomPairs, mol1, mol2, bondMatchTable, corrGraph,
                           corrGraphNumConns, atomPairStarts, cEdges);

#if 0
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
  // Running enumerate_z_cliques over each node in the corrGraph in turn.
#if 1
  std::vector<char> t(corrGraph.size(), 0);
  std::vector<unsigned int> p;
  std::vector<unsigned int> d;
  std::vector<unsigned int> s;
  // Start each round of enumerate_z_cliques at an atom pair given by
  // the atomPairStarts.
  for (size_t aps = 0U; aps < atomPairStarts.size(); ++aps) {
    for (unsigned int i = 0u; i < atomPairStarts[aps].first; ++i) {
      size_t u = atomPairStarts[aps].second + i;
      p.clear();
      d.clear();
      s.clear();
      for (size_t v = 0; v < corrGraph[u].size(); ++v) {
        if (corrGraph[v][u] == 'c') {
          if (t[v]) {
            s.push_back(v);
          } else {
            p.push_back(v);
          }
        } else if (corrGraph[v][u] == 'd') {
          d.push_back(v);
        }
      }
      clique.clear();
      clique.push_back(u);
      enumerate_z_cliques(clique, p, d, s, t, corrGraph, corrGraphNumConns,
                          cEdges, rawMaxCliques);
      t[u] = 1;
      if (!rawMaxCliques.empty() &&
          rawMaxCliques.front().size() > minMCSSSize) {
        minMCSSSize = rawMaxCliques.front().size();
      }
      // If we have tried minMCSSSize atoms from mol1 as start points
      // and haven't found a clique of the requisite size, we won't be
      // able to from now on in.
      if (aps + minMCSSSize > atomPairStarts.size()) {
        break;
      }
    }
    if (aps + minMCSSSize > atomPairStarts.size()) {
      break;
    }
  }
  // std::cout << "Number of cliques: " << rawMaxCliques.size() << std::endl;
  // for (auto &c : rawMaxCliques) {
  //   std::ranges::sort(c);
  // }
  // std::ranges::sort(rawMaxCliques);
  // for (const auto &clique : rawMaxCliques) {
  //   std::cout << clique.size() << " :";
  //   for (const auto &c : clique) {
  //     std::cout << " " << c;
  //   }
  //   std::cout << std::endl;
  // }
#endif
  if (rawMaxCliques.empty()) {
    return;
  } else {
    if (rawMaxCliques.front().size() < minMCSSSize) {
      return;
    }
  }
  maxCliques.clear();
  maxCliques.reserve(rawMaxCliques.size());
  std::vector<std::vector<unsigned int>> cliqueAtoms;
  cliqueAtoms.reserve(rawMaxCliques.size());
  for (const auto &rmc : rawMaxCliques) {
    if (uniquify) {
      std::vector<unsigned int> theseCliqueAtoms;
      theseCliqueAtoms.resize(2 * rmc.size());
      for (size_t i = 0; i < rmc.size(); ++i) {
        theseCliqueAtoms[i] = atomPairs[rmc[i]].first;
        theseCliqueAtoms[i + rmc.size()] = atomPairs[rmc[i]].second;
      }
      std::sort(theseCliqueAtoms.begin(),
                theseCliqueAtoms.begin() + rmc.size());
      std::sort(theseCliqueAtoms.begin() + rmc.size(), theseCliqueAtoms.end());
      if (std::find(cliqueAtoms.begin(), cliqueAtoms.end(), theseCliqueAtoms) ==
          cliqueAtoms.end()) {
        cliqueAtoms.emplace_back(std::move(theseCliqueAtoms));
      } else {
        continue;
      }
    }
    std::vector<std::pair<unsigned int, unsigned int>> newClique;
    newClique.reserve(rmc.size());
    for (const auto cm : rmc) {
      newClique.push_back(atomPairs[cm]);
    }
    std::ranges::sort(newClique);
    // std::cout << makeSMARTSFromMCSS(mol1, newClique) << std::endl;
    maxCliques.push_back(std::move(newClique));
  }
#if 0
  KochBK<unsigned int> g;
  for (unsigned int i = 0; i < corrGraph.size(); ++i) {
    for (unsigned int j = 0; j < corrGraph[i].size(); ++j) {
      if (corrGraph[i][j]) {
        g.add_edge(i, j, corrGraph[i][j]);
      }
    }
  }
  auto cliques = g.enumerate_c_cliques();
  std::ranges::sort(cliques, [](const auto &a, const auto &b) -> bool {
    return a.size() < b.size();
  });
  std::cout << "c-cliques:\n";
  for (auto &C : cliques) {
    std::cout << C.size() << " { ";
    for (auto &v : C) std::cout << v << " ";
    std::cout << "}\n";
  }
#endif
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