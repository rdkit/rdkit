//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <sstream>

#include "Digraph.h"
#include "CIPMol.h"
#include "Node.h"
#include "Edge.h"

namespace RDKit {
namespace CIPLabeler {

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

void Digraph::addEdge(Node *beg, Bond *bond, Node *end) {
  auto e = new Edge(beg, end, bond);
  beg->add(e);
  end->add(e);
}

Digraph::Digraph() = default;

Digraph::Digraph(CIPMol *mol) : mol{mol} {}

Digraph::Digraph(CIPMol *mol, Atom *atom) : mol{mol} { init(atom); }

Digraph::~Digraph() { delete getCurrRoot(); }

CIPMol *Digraph::getMol() const { return mol; };

Node *Digraph::getRoot() const { return root; };

Node *Digraph::getCurrRoot() const {
  return tmproot == nullptr ? root : tmproot;
}

int Digraph::getNumNodes() const { return numNodes; }

std::vector<Node *> Digraph::getNodes(Atom *atom) const {
  auto result = std::vector<Node *>{};
  auto queue = std::vector<Node *>({getCurrRoot()});

  for (auto pos = 0u; pos < queue.size(); ++pos) {
    const auto node = queue[pos];

    if (atom == node->getAtom()) {
      result.push_back(node);
    }
    for (const auto &e : node->getEdges()) {
      if (!e->isBeg(node)) {
        continue;
      }
      queue.push_back(e->getEnd());
    }
  }
  return result;
}

/**
 * Access the reference atom for Rule 6 (if one is set).
 */
Atom *Digraph::getRule6Ref() const { return rule6Ref; }

/**
 * Used exclusively for Rule 6, we set one atom as the reference.
 * @param ref reference atom
 */
void Digraph::setRule6Ref(Atom *ref) { rule6Ref = ref; }

Node *Digraph::init(Atom *atom) {
  PRECONDITION(atom, "cannot init digraph on a nullptr")
  auto visit = std::vector<char>(mol->getNumAtoms());
  visit[atom->getIdx()] = 1;

  auto dist = 1;
  auto flags = 0x0;
  auto frac = Fraction(atom->getAtomicNum());

  root = new Node(this, std::move(visit), atom, std::move(frac), dist, flags);

  ++numNodes;

  return root;
}

/**
 * Sets the root node of this digraph by flipping the directions
 * of edges as required.
 *
 * @param newroot the new root
 */
void Digraph::changeRoot(Node *newroot) {
  auto queue = std::vector<Node *>({newroot});

  auto toflip = std::vector<Edge *>();
  for (auto pos = 0u; pos < queue.size(); ++pos) {
    const auto node = queue[pos];

    for (const auto &e : node->getEdges()) {
      if (e->isEnd(node)) {
        toflip.push_back(e);
        queue.push_back(e->getBeg());
      }
    }
  }
  for (auto &e : toflip) {
    e->flip();
  }
  tmproot = newroot;
}

void Digraph::expand(Node *beg) {
  const auto &atom = beg->getAtom();
  const auto &edges = beg->getEdges();
  const auto &prev =
      edges.size() > 0 && !edges[0]->isBeg(beg) ? edges[0]->getBond() : nullptr;

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
    const auto &nbr = bond->getOtherAtom(atom);
    const int nbrIdx = nbr->getIdx();
    const int bord = mol->getBondOrder(bond);
    const int virtual_nodes = bord - 1;

    if (beg->visit[nbrIdx] == 0) {
      auto end = beg->newChild(nbrIdx, nbr);
      ++numNodes;
      addEdge(beg, bond, end);

      // duplicate nodes for bond orders (except for root atoms...)
      // for example >S=O
      if (root != beg) {
        if (atom->getFormalCharge() < 0 &&
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

      if (atom->getFormalCharge() < 0 &&
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
  const int hcnt = atom->getTotalNumHs();
  for (int i = 0; i < hcnt; ++i) {
    auto end = beg->newImplicitHydrogenChild();
    ++numNodes;
    addEdge(beg, nullptr, end);
  }
}

void Digraph::expandAll() {
  auto queue = std::vector<Node *>({root});

  for (auto pos = 0u; pos < queue.size(); ++pos) {
    const auto node = queue[pos];

    for (const auto &e : node->getEdges()) {
      if (!e->isBeg(node)) {
        continue;
      }
      if (!e->getEnd()->isTerminal()) {
        queue.push_back(e->getEnd());
      }
    }
  }
}

} // namespace CIPLabeler
} // namespace RDKit
