#pragma once
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>
#include <stdexcept>
// This is an implementation of the algorithm used in TwoMolMCSS.h
// the ChatGPT 5.0 produced for me.  It is inefficient but correct
// and therefore a useful reference.  It was also pretty awesome to
// watch it being created.  The prompt was:
// "I would like an implementation of the koch modification of the
// Bron Kerbosch algorithm, specifically algorithm 5."
// It had to think about it for just under 2 minutes.

/**
 * KochBK - Koch's modification of Bron–Kerbosch algorithm (Algorithm 5).
 *
 * Reference: I. Koch, "Enumerating all connected maximal common subgraphs in
 * two graphs," Theoretical Computer Science 250(1–2):1–30, 2001.
 *
 * Graph model:
 *  - Undirected, simple graph.
 *  - Each edge labeled either 'c' or 'd'.
 *  - Enumerates all maximal cliques spanned by c-edges ("c-cliques").
 *
 * Template parameter: Vertex type (must be hashable, e.g., int, string).
 */
template <typename Vertex>
class KochBK {
 public:
  using VSet = std::unordered_set<Vertex>;
  using Clique = std::unordered_set<Vertex>;

  // ---- Graph construction ----
  void add_edge(const Vertex &u, const Vertex &v, char t = 'c') {
    if (u == v) return;
    if (t != 'c' && t != 'd')
      throw std::invalid_argument("edge type must be 'c' or 'd'");
    V.insert(u);
    V.insert(v);
    neigh[u].insert(v);
    neigh[v].insert(u);
    if (t == 'c') {
      c_neigh[u].insert(v);
      c_neigh[v].insert(u);
    } else {
      d_neigh[u].insert(v);
      d_neigh[v].insert(u);
    }
  }

  // ---- Public API ----
  std::vector<Clique> enumerate_c_cliques() {
    std::vector<Clique> out;
    VSet T;
    for (const Vertex &u : V) {
      VSet P, D, S;
      for (const Vertex &v : neigh[u]) {
        if (c_neigh[u].count(v)) {
          if (T.count(v))
            S.insert(v);
          else
            P.insert(v);
        } else if (d_neigh[u].count(v)) {
          D.insert(v);
        }
      }
      alg5({u}, P, D, S, T, out);
      T.insert(u);
    }
    return out;
  }

 private:
  // ---- Graph storage ----
  VSet V;
  std::unordered_map<Vertex, VSet> neigh;    // all neighbors
  std::unordered_map<Vertex, VSet> c_neigh;  // c-neighbors only
  std::unordered_map<Vertex, VSet> d_neigh;  // d-neighbors only

  // ---- Helpers ----
  bool adjacent(const Vertex &u, const Vertex &v) const {
    auto it = neigh.find(u);
    return it != neigh.end() && it->second.count(v);
  }
  bool adjacent_c(const Vertex &u, const Vertex &v) const {
    auto it = c_neigh.find(u);
    return it != c_neigh.end() && it->second.count(v);
  }

  bool cpath_exists_to_non_ut_neighbor(const Vertex &start_u, const Vertex &ut,
                                       const VSet &D) const {
    if (D.empty()) return false;
    std::unordered_set<Vertex> visited;
    std::queue<Vertex> q;
    visited.insert(start_u);
    q.push(start_u);

    while (!q.empty()) {
      Vertex x = q.front();
      q.pop();
      auto it = c_neigh.find(x);
      if (it == c_neigh.end()) continue;
      for (const Vertex &w : it->second) {
        if (w == start_u) continue;
        if (!D.count(w) || visited.count(w)) continue;
        if (!adjacent(ut, w)) return true;
        visited.insert(w);
        q.push(w);
      }
    }
    return false;
  }

  // ---- Algorithm 5 ----
  void alg5(Clique C, VSet P, VSet D, VSet S, const VSet &T,
            std::vector<Clique> &out) {
    if (P.empty()) {
      if (S.empty()) {
        out.push_back(C);
      }
      return;
    }

    Vertex ut = *P.begin();  // pivot
    std::vector<Vertex> Pvec(P.begin(), P.end());

    for (const Vertex &ui : Pvec) {
      bool cond =
          (!adjacent(ui, ut)) || cpath_exists_to_non_ut_neighbor(ui, ut, D);
      if (!cond) continue;

      P.erase(ui);

      VSet P0 = P;
      VSet D0 = D;
      VSet S0 = S;
      VSet N = neigh[ui];

      for (const Vertex &v : std::vector<Vertex>(D0.begin(), D0.end())) {
        if (P.count(v)) {
          P0.insert(v);
          D0.erase(v);
        } else if (D.count(v)) {
          if (adjacent_c(v, ui)) {
            if (T.count(v))
              S0.insert(v);
            else
              P0.insert(v);
            D0.erase(v);
          }
        } else if (S.count(v)) {
          S0.insert(v);
        }
      }

      // Restrict sets to neighbors of ui
      VSet newP, newD, newS;
      for (const Vertex &x : P0)
        if (N.count(x)) newP.insert(x);
      for (const Vertex &x : D0)
        if (N.count(x)) newD.insert(x);
      for (const Vertex &x : S0)
        if (N.count(x)) newS.insert(x);

      Clique Cnext = C;
      Cnext.insert(ui);
      alg5(Cnext, newP, newD, newS, T, out);

      S.insert(ui);
    }
  }
};
