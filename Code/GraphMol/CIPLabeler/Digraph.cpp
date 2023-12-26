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

#include <list>
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
}  // namespace

Node &Digraph::addNode(std::vector<char> &&visit, Atom *atom,
                       boost::rational<int> &&frac, int dist, int flags) {
  d_nodes.emplace_back(this, std::move(visit), atom, std::move(frac), dist,
                       flags);
  return d_nodes.back();
}

void Digraph::addEdge(Node *beg, Bond *bond, Node *end) {
  d_edges.emplace_back(beg, end, bond);
  auto &e = d_edges.back();
  beg->add(&e);
  end->add(&e);
}

Digraph::Digraph(const CIPMol &mol, Atom *atom, bool atropisomerMode)
    : d_mol{mol} {
  PRECONDITION(atom, "cannot init digraph on a nullptr")

  auto visit = std::vector<char>(d_mol.getNumAtoms());
  visit[atom->getIdx()] = 1;

  auto dist = 1;
  auto flags = 0x0;
  auto atomic_num = atom->getAtomicNum();

  dp_root = &addNode(std::move(visit), atom, atomic_num, dist, flags);
  dp_origin = dp_root;
  d_atropisomerMode = atropisomerMode;
}

const CIPMol &Digraph::getMol() const { return d_mol; };

Node *Digraph::getOriginalRoot() const { return dp_origin; };

Node *Digraph::getCurrentRoot() const { return dp_root; }

int Digraph::getNumNodes() const { return d_nodes.size(); }

std::vector<Node *> Digraph::getNodes(Atom *atom) const {
  std::vector<Node *> result;
  auto queue = std::list<Node *>({getCurrentRoot()});

  for (const auto &node : queue) {
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
Atom *Digraph::getRule6Ref() const { return dp_rule6Ref; }

/**
 * Used exclusively for Rule 6, we set one atom as the reference.
 * @param ref reference atom
 */
void Digraph::setRule6Ref(Atom *ref) { dp_rule6Ref = ref; }

/**
 * Sets the root node of this digraph by flipping the directions
 * of edges as required.
 *
 * @param newroot the new root
 */
void Digraph::changeRoot(Node *newroot) {
  std::vector<Edge *> toflip;
  auto queue = std::list<Node *>({newroot});
  for (const auto &node : queue) {
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
  dp_root = newroot;
}

void Digraph::expand(Node *beg) {
  const auto &atom = beg->getAtom();
  const auto &edges = beg->getEdges();
  const auto &prev =
      edges.size() > 0 && !edges[0]->isBeg(beg) ? edges[0]->getBond() : nullptr;

  if (MAX_NODE_DIST > 0 && beg->getDistance() > MAX_NODE_DIST) {
    return;
  }
  if (MAX_NODE_COUNT > 0 && d_nodes.size() >= MAX_NODE_COUNT) {
    std::stringstream errmsg;
    errmsg << "Digraph generation failed: more than " << MAX_NODE_COUNT
           << "nodes found.";
    throw TooManyNodesException(errmsg.str());
  }

  // create 'explicit' nodes
  for (const auto &bond : d_mol.getBonds(atom)) {
    const auto &nbr = bond->getOtherAtom(atom);
    const int nbrIdx = nbr->getIdx();
    const int bord = d_mol.getBondOrder(bond);
    const int virtual_nodes = bord - 1;

    if (!beg->isVisited(nbrIdx)) {
      auto end = beg->newChild(nbrIdx, nbr);
      addEdge(beg, bond, end);

      // duplicate nodes for bond orders (except for root atoms...)
      // for example >S=O
      if (dp_origin != beg || d_atropisomerMode) {
        if (atom->getFormalCharge() < 0 &&
            d_mol.getFractionalAtomicNum(atom).denominator() > 1) {
          end = beg->newBondDuplicateChild(nbrIdx, nbr);
          addEdge(beg, bond, end);
        } else {
          for (int i = 0; i < virtual_nodes; ++i) {
            end = beg->newBondDuplicateChild(nbrIdx, nbr);
            addEdge(beg, bond, end);
          }
        }
      }
    } else if (bond == prev) {  // bond order expansion (backwards)
      if (dp_origin->getAtom() != nbr || d_atropisomerMode) {
        for (int i = 0; i < virtual_nodes; ++i) {
          auto end = beg->newBondDuplicateChild(nbrIdx, nbr);
          addEdge(beg, bond, end);
        }
      }
    } else {  // ring closures
      auto end = beg->newRingDuplicateChild(nbrIdx, nbr);
      addEdge(beg, bond, end);

      if (atom->getFormalCharge() < 0 &&
          d_mol.getFractionalAtomicNum(atom).denominator() > 1) {
        end = beg->newBondDuplicateChild(nbrIdx, nbr);
        addEdge(beg, bond, end);
      } else {
        for (int i = 0; i < virtual_nodes; ++i) {
          end = beg->newBondDuplicateChild(nbrIdx, nbr);
          addEdge(beg, bond, end);
        }
      }
    }
  }

  // Create implicit hydrogen nodes
  const int hcnt = atom->getTotalNumHs();
  for (int i = 0; i < hcnt; ++i) {
    auto end = beg->newImplicitHydrogenChild();
    addEdge(beg, nullptr, end);
  }
}

}  // namespace CIPLabeler
}  // namespace RDKit
