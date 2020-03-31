//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <queue>
#include <sstream>
#include <vector>

#include "BaseMol.hpp"
#include "Edge.hpp"
#include "TooManyNodesException.hpp"

namespace RDKit {
namespace NewCIPLabelling {

template <typename A, typename B> class Node;

namespace {

/**
 * Upper limit on the size of the digraph, stops out of memory error with a
 * more graceful failure. 0=Infinite
 */
const int MAX_NODE_COUNT = 100000;

/**
 * Used for debugging only, 0=Infinite
 */
const int MAX_NODE_DIST = 0;
} // namespace

/**
 * A class to hold directed acyclic graphs representing the molecule.
 *
 * The root of the DAG is one of the foci of the configuration for
 * which the label is being calculated. The tmproot may be set to
 * other nodes that may become relevant in the calculation.
 *
 */
template <typename A, typename B> class Digraph {
private:
  BaseMol<A, B> *mol = nullptr;
  Node<A, B> *root = nullptr;
  Node<A, B> *tmproot = nullptr;
  int numNodes = 0;
  A rule6Ref = nullptr;

  void addEdge(Node<A, B> *beg, B bond, Node<A, B> *end) {
    auto e = new Edge<A, B>(beg, end, bond);
    beg->add(e);
    end->add(e);
  }

public:
  Digraph() = default;

  Digraph(BaseMol<A, B> *mol) : mol{mol} {}

  Digraph(BaseMol<A, B> *mol, A atom) : mol{mol} { init(atom); }

  ~Digraph() { delete getCurrRoot(); }

  BaseMol<A, B> *getMol() const { return mol; };

  Node<A, B> *getRoot() const { return root; };

  Node<A, B> *getCurrRoot() const {
    return tmproot == nullptr ? root : tmproot;
  }

  int getNumNodes() const { return numNodes; }

  std::vector<Node<A, B> *> getNodes(A atom) const {
    auto result = std::vector<Node<A, B> *>{};
    auto queue = std::queue<Node<A, B> *>{};
    queue.push(getCurrRoot());

    while (!queue.empty()) {
      auto node = queue.front();
      queue.pop();

      if (atom == node->getAtom()) {
        result.push_back(node);
      }
      for (const auto &e : node->getEdges()) {
        if (!e->isBeg(node)) {
          continue;
        }
        queue.push(e->getEnd());
      }
    }
    return result;
  }

  /**
   * Access the reference atom for Rule 6 (if one is set).
   */
  A getRule6Ref() const { return rule6Ref; }

  /**
   * Used exclusively for Rule 6, we set one atom as the reference.
   * @param ref reference atom
   */
  void setRule6Ref(A ref) { rule6Ref = ref; }

  Node<A, B> *init(A atom) {
    auto visit = std::vector<char>(mol->getNumAtoms());
    visit[mol->getAtomIdx(atom)] = 1;

    auto atomic_num_num = mol->getAtomicNum(atom);
    auto atomic_num_den = 1;
    auto dist = 1;
    auto flags = 0x0;

    root = new Node<A, B>(this, std::move(visit), atom, atomic_num_num,
                          atomic_num_den, dist, flags);

    ++numNodes;

    return root;
  }

  /**
   * Sets the root node of this digraph by flipping the directions
   * of edges as required.
   *
   * @param newroot the new root
   */
  void changeRoot(Node<A, B> *newroot) {
    auto queue = std::queue<Node<A, B> *>();
    queue.push(newroot);

    auto toflip = std::vector<Edge<A, B> *>();
    while (!queue.empty()) {
      auto node = queue.front();
      queue.pop();

      for (const auto &e : node->getEdges()) {
        if (e->isEnd(node)) {
          toflip.push_back(e);
          queue.push(e->getBeg());
        }
      }
    }
    for (auto &e : toflip) {
      e->flip();
    }
    tmproot = newroot;
  }

  void expand(Node<A, B> *beg) {
    const A atom = beg->getAtom();
    const auto &edges = beg->getEdges();
    const B prev = edges.size() > 0 && !edges[0]->isBeg(beg)
                       ? edges[0]->getBond()
                       : nullptr;

    if (MAX_NODE_DIST > 0 && beg->getDistance() > MAX_NODE_DIST) {
      return;
    }
    if (MAX_NODE_COUNT > 0 && numNodes >= MAX_NODE_COUNT) {
      std::stringstream errmsg;
      errmsg << "Digraph generation failed: more than " << MAX_NODE_COUNT
             << "nodes found.";
      throw TooManyNodesException(errmsg.str());
    }

    // create 'explicit' nodes
    for (const auto &bond : mol->getBonds(atom)) {
      const A nbr = mol->getOther(bond, atom);
      const int nbrIdx = mol->getAtomIdx(nbr);
      const int bord = mol->getBondOrder(bond);
      const int virtual_nodes = bord - 1;

      if (beg->visit[nbrIdx] == 0) {
        auto end = beg->newChild(nbrIdx, nbr);
        ++numNodes;
        addEdge(beg, bond, end);

        // duplicate nodes for bond orders (except for root atoms...)
        // for example >S=O
        if (root != beg) {
          if (mol->getCharge(atom) < 0 &&
              mol->getFractionalAtomicNum(atom).denominator() > 1) {
            end = beg->newBondDuplicateChild(nbrIdx, nbr);
            ++numNodes;
            addEdge(beg, bond, end);
          } else {
            for (int i = 0; i < virtual_nodes; ++i) {
              end = beg->newBondDuplicateChild(nbrIdx, nbr);
              ++numNodes;
              addEdge(beg, bond, end);
            }
          }
        }
      } else if (bond == prev) { // bond order expansion (backwards)
        if (root->getAtom() != nbr) {
          for (int i = 0; i < virtual_nodes; ++i) {
            auto end = beg->newBondDuplicateChild(nbrIdx, nbr);
            ++numNodes;
            addEdge(beg, bond, end);
          }
        }
      } else { // ring closures
        auto end = beg->newRingDuplicateChild(nbrIdx, nbr);
        ++numNodes;
        addEdge(beg, bond, end);

        if (mol->getCharge(atom) < 0 &&
            mol->getFractionalAtomicNum(atom).denominator() > 1) {
          end = beg->newBondDuplicateChild(nbrIdx, nbr);
          ++numNodes;
          addEdge(beg, bond, end);
        } else {
          for (int i = 0; i < virtual_nodes; ++i) {
            end = beg->newBondDuplicateChild(nbrIdx, nbr);
            ++numNodes;
            addEdge(beg, bond, end);
          }
        }
      }
    }

    // Create implicit hydrogen nodes
    const int hcnt = mol->getNumHydrogens(atom);
    for (int i = 0; i < hcnt; ++i) {
      auto end = beg->newImplicitHydrogenChild();
      ++numNodes;
      addEdge(beg, nullptr, end);
    }
  }

  void expandAll() {
    auto queue = std::queue<Node<A, B> *>{};
    queue.push(root);
    while (!queue.empty()) {
      auto node = queue.front();
      queue.pop();
      for (const auto &e : node->getEdges()) {
        if (!e->isBeg(node)) {
          continue;
        }
        if (!e->getEnd()->isTerminal()) {
          queue.push(e->getEnd());
        }
      }
    }
  }
};

} // namespace NewCIPLabelling
} // namespace RDKit
