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
#include <chrono>
#include <numeric>
#include <queue>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/FMCS/MatchTable.h>
#if 0
// This is ChatGPT's reference implementation.  It's useful to include
// for debugging purposes.
#include <GraphMol/FMCS/chatgpt.h>
#endif
#include <GraphMol/SmilesParse/SmartsWrite.h>

using Clock = std::chrono::steady_clock;
using TimePoint = std::chrono::time_point<Clock>;

namespace RDKit {
namespace {
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
  // If they're both bonded and not the same type they count as not
  // connected in the correspondence graph but not d-Edges as they
  // are clearly connected in the molecules.
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
  // Make a list of where each atom pair starts for each atom in molecule 1,
  // along with the number of pairs that atom has with atoms in molecule 2.
  // Sort this into ascending number of pairs.  The start nodes in the
  // correspondence graph will be taken in this order which can reduce the
  // search time on average.
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
  // Sort the cEdges into ascending number of connections.  This makes
  // cPathToNonUtNbor a bit faster.
  for (auto &ce : cEdges) {
    std::ranges::sort(ce, [&](const auto &a, const auto &b) -> bool {
      return cEdges[a].size() < cEdges[b].size();
    });
  }
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
  std::queue<unsigned int> q;
  q.push(ui);
  tried.reset();
  tried[ui] = true;
  while (!q.empty()) {
    auto x = q.front();
    q.pop();
    for (auto xn : cEdges[x]) {
      if (tried[xn] || !(*inD)[xn]) {
        continue;
      }
      // It's in the path and in d, so if it's not a n'bour of ut
      // return true.
      if (!nUt[xn]) {
        return true;
      }
      q.push(xn);
      tried[xn] = true;
    }
  }
  return false;
}

bool inline inVector(const std::vector<unsigned int> &v, unsigned int val) {
  return std::find(v.begin(), v.end(), val) != v.end();
}

// Notation from Koch, Theoretical Computer Science 250 (2001) 1–30
// Algorithm 5.
// c is the current clique
// p is the set of nodes that can be added to c. They are neighbours of
// the last node in c (u) via c-edges.
// d is the set of nodes that cannot be added to c because they are
// neighbours of u via d-edges.
// s is the set of nodes already tried as u in this step in the
// recursion.
// t is the set of nodes that have already been used to initialise a
// clique
void enumerate_z_cliques(std::vector<unsigned int> &c,
                         const std::vector<unsigned int> p,
                         const std::vector<unsigned int> &d,
                         std::vector<unsigned int> &s,
                         const std::vector<char> &t,
                         const std::vector<std::vector<char>> &corrGraph,
                         const std::vector<unsigned int> &corrGraphNumConns,
                         const std::vector<std::vector<unsigned int>> &cEdges,
                         std::vector<std::vector<unsigned int>> &maxCliques) {
  if (p.empty()) {
    if (s.empty() && c.size() > 1) {
      if (maxCliques.empty()) {
        maxCliques.push_back({c.begin(), c.end()});
      } else {
        if (c.size() > maxCliques.front().size()) {
          maxCliques.clear();
        } else if (c.size() < maxCliques.front().size()) {
          return;
        }
        maxCliques.push_back(c);
      }
    }
    return;
  }
  // Select a pivot node ut.  Koch doesn't say how to do this, but the BK
  // algorithm uses a maxmimally connected one so do that here too.
  auto ut = *std::max_element(
      p.begin(), p.end(), [&](const auto &a, const auto &b) -> bool {
        return corrGraphNumConns[a] < corrGraphNumConns[b];
      });
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
      // Form a new p which is the members of the current p that are
      // neighbours of ui
      std::vector<unsigned int> newP;
      for (auto pe : p) {
        if (n[pe]) {
          newP.push_back(pe);
        }
      }
      // and the same for newS
      std::vector<unsigned int> newS;
      for (auto se : s) {
        if (n[se]) {
          newS.push_back(se);
        }
      }
      // Now form newD and fiddle with newS and newP at the same time.  This
      // is the nub of Koch's algorithm.  It's somewhat re-worked from her
      // description which was a bit inefficient.
      std::vector<unsigned int> newD;
      newD.reserve(newD.size());
      for (auto v : d) {
        bool addD = true;
        if (inVector(p, v) && n[v]) {
          newP.push_back(v);
        } else {
          // can v be added to P?
          if (n[v] == 'c') {
            if (t[v] and n[v]) {
              newS.push_back(v);
            } else if (n[v]) {
              newP.push_back(v);
            }
            addD = false;
          } else if (n[v] && inVector(s, v)) {
            newS.push_back(v);
          }
        }
        if (addD && n[v]) {
          newD.push_back(v);
        }
      }
      // the new clique is the current clique plus ui.
      c.push_back(ui);
      enumerate_z_cliques(c, newP, newD, newS, t, corrGraph, corrGraphNumConns,
                          cEdges, maxCliques);
      // Set it up for the next round, marking ui as used.
      c.pop_back();
      s.push_back(ui);
    }
  }
}
}  // namespace

void TwoMolMCSS(const ROMol &mol1, const ROMol &mol2, unsigned int minMCSSSize,
                const FMCS::MatchTable &atomMatchTable,
                const FMCS::MatchTable &bondMatchTable, bool uniquify,
                unsigned int timeOut,
                std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
                    &maxCliques) {
  const TimePoint *endTime = nullptr;
  TimePoint endTimePt;
  if (timeOut > 0) {
    endTimePt = Clock::now() + std::chrono::seconds(timeOut);
    endTime = &endTimePt;
  }

  std::vector<std::pair<unsigned int, unsigned int>> atomPairs;
  if (minMCSSSize > std::min(mol1.getNumAtoms(), mol2.getNumAtoms())) {
    return;
  }
  buildAtomPairs(mol1, mol2, atomMatchTable, atomPairs);
  if (atomPairs.empty()) {
    return;
  }
  std::vector<std::vector<char>> corrGraph(
      atomPairs.size(), std::vector<char>(atomPairs.size(), 0));
  std::vector<unsigned int> corrGraphNumConns(atomPairs.size(), 0);
  std::vector<std::pair<unsigned int, unsigned int>> atomPairStarts;
  std::vector<std::vector<unsigned int>> cEdges(atomPairs.size());
  buildCorrespondenceGraph(atomPairs, mol1, mol2, bondMatchTable, corrGraph,
                           corrGraphNumConns, atomPairStarts, cEdges);

  std::vector<unsigned int> clique;
  std::vector<std::vector<unsigned int>> rawMaxCliques;
  bool timedOut = false;
  // Running enumerate_z_cliques over each node in the corrGraph in turn.
  std::vector<char> t(corrGraph.size(), 0);
  std::vector<unsigned int> p;
  std::vector<unsigned int> d;
  std::vector<unsigned int> s;
  // Start each round of enumerate_z_cliques at an atom pair given by
  // the atomPairStarts.
  for (size_t aps = 0U; aps < atomPairStarts.size(); ++aps) {
    for (unsigned int i = 0u; i < atomPairStarts[aps].first; ++i) {
      if (endTime != nullptr && Clock::now() > *endTime) {
        BOOST_LOG(rdWarningLog) << "Timed out.\n";
        timedOut = true;
        break;
      }

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
    }
    // If the number of start atoms used so far + current minMCSSSize
    // exceeds the number of atomPairStarts then the clique can't
    // be improved on, so it's ok to stop.
    if (timedOut || aps + minMCSSSize > atomPairStarts.size()) {
      break;
    }
  }

  if (rawMaxCliques.empty() || rawMaxCliques.front().size() < minMCSSSize) {
    return;
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
    maxCliques.push_back(std::move(newClique));
  }
#if 0
  // This is using ChatGPT's reference implementation.  It's useful to
  // have it for debugging purposes.
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
}  // namespace RDKit
