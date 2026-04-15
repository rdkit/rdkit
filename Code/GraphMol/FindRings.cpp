//
//  Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Rings.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDLog.h>

#include <boost/dynamic_bitset.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace {
using namespace RDKit;

// normalizes a ring by rotating/reversing it so that the first atom
// is the one with the smallest index, and the second atom is the neighbor
// to the first one that again has the smallest index.
// This change should have a small performance footprint while it helps
// keeping test results consistent when making changes to ring detection.
void normalize_ring(std::vector<int> &ring) {
  auto newStart = std::ranges::min_element(ring);
  std::ranges::rotate(ring, newStart);

  if (ring.back() < ring[1]) {
    // we don't need to move the central element!
    auto numPairsToMove = (ring.size() - 1) / 2;
    auto front = ring.begin() + 1;
    std::swap_ranges(front, front + numPairsToMove, ring.rbegin());
  }
}

std::vector<int> rdl_cycle_to_atom_ring(RDL_cycle *cycle) {
  std::vector<int> ring;
  ring.reserve(cycle->weight);

  // Edges in a cycle are not returned in iteration order.
  // so we need to take care of that while we convert them
  // into an atom ring.
  boost::dynamic_bitset<> unseen_edges(cycle->weight);
  unseen_edges.set();

  ring.push_back(cycle->edges[0][0]);
  ring.push_back(cycle->edges[0][1]);
  unseen_edges.set(0, false);

  while (ring.size() < cycle->weight) {
    // Note we don't want to close the cycle: that would
    // add the initial atom at the end too.
    for (auto edgeIdx = unseen_edges.find_first();
         edgeIdx != boost::dynamic_bitset<>::npos;
         edgeIdx = unseen_edges.find_next(edgeIdx)) {
      auto edge = cycle->edges[edgeIdx];
      for (auto j = 0; j < 2; ++j) {
        if (static_cast<unsigned int>(ring.back()) == edge[j]) {
          ring.push_back(edge[1 - j]);
          if (ring.size() == cycle->weight) {
            unseen_edges.reset();
          } else {
            unseen_edges.set(edgeIdx, false);
          }
          break;
        }
      }
    }
  }

  // For consistency, normalize the ring
  normalize_ring(ring);

  return ring;
}
}  // namespace

namespace RingUtils {

void convertToBonds(const INT_VECT &ring, INT_VECT &bondRing,
                    const ROMol &mol) {
  const auto rsiz = rdcast<unsigned int>(ring.size());
  bondRing.resize(rsiz);
  for (unsigned int i = 0; i < (rsiz - 1); i++) {
    const Bond *bnd = mol.getBondBetweenAtoms(ring[i], ring[i + 1]);
    if (!bnd) {
      throw ValueErrorException("expected bond not found");
    }
    bondRing[i] = bnd->getIdx();
  }
  // bond from last to first atom
  const Bond *bnd = mol.getBondBetweenAtoms(ring[rsiz - 1], ring[0]);
  if (!bnd) {
    throw ValueErrorException("expected bond not found");
  }

  bondRing[rsiz - 1] = bnd->getIdx();
}

void convertToBonds(const VECT_INT_VECT &res, VECT_INT_VECT &brings,
                    const ROMol &mol) {
  brings.reserve(res.size());
  for (const auto &ring : res) {
    INT_VECT bring;
    convertToBonds(ring, bring, mol);
    brings.push_back(bring);
  }
}

auto ringComparer = [](const auto &v1, const auto &v2) {
  if (v1.size() == v2.size()) {
    return v1 < v2;
  }
  return v1.size() < v2.size();
};
}  // end of namespace RingUtils

namespace FindRings {
using namespace RDKit;

void storeRingInfo(const ROMol &mol, const INT_VECT &ring) {
  INT_VECT bondIndices;
  RingUtils::convertToBonds(ring, bondIndices, mol);
  mol.getRingInfo()->addRing(ring, bondIndices);
}

void storeRingsInfo(const ROMol &mol, const VECT_INT_VECT &rings) {
  mol.getRingInfo()->preallocate(rings.size(), rings.size());
  for (const auto &ring : rings) {
    storeRingInfo(mol, ring);
  }
}
}  // namespace FindRings

namespace RDKit {
namespace MolOps {
int findSSSR(const ROMol &mol, VECT_INT_VECT *res, bool includeDativeBonds,
             bool includeHydrogenBonds) {
  if (!res) {
    VECT_INT_VECT rings;
    return findSSSR(mol, rings, includeDativeBonds, includeHydrogenBonds);
  } else {
    return findSSSR(mol, (*res), includeDativeBonds, includeHydrogenBonds);
  }
}

int findSSSR(const ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds,
             bool includeHydrogenBonds) {
  res.clear();
  auto ringInfo = mol.getRingInfo();
  if (ringInfo->isInitialized()) {
    ringInfo->reset();
  }
  ringInfo->initialize(FIND_RING_TYPE_SSSR);
  ringInfo->preallocate(mol.getNumAtoms(), mol.getNumBonds());

  findRingFamilies(mol, includeDativeBonds, includeHydrogenBonds);
  auto urfdata = mol.getRingInfo()->dp_urfData.get();

  RDL_cycle **sssr = nullptr;
  auto sssrSize = RDL_getSSSR(urfdata, &sssr);
  if (sssrSize == RDL_INVALID_RESULT) {
    throw ValueErrorException("Failed finding a SSSR for the mol.");
  }

  res.reserve(sssrSize);
  for (unsigned int i = 0; i < sssrSize; ++i) {
    auto ring = rdl_cycle_to_atom_ring(sssr[i]);
    res.push_back(std::move(ring));
  }
  RDL_deleteCycles(sssr, sssrSize);

  std::ranges::sort(res, RingUtils::ringComparer);

  FindRings::storeRingsInfo(mol, res);

  // update the ring memberships of atoms and bonds in the molecule:
  // store the SSSR rings on the molecule as a property
  // we will ignore any existing SSSRs on the molecule - simply overwrite
  return rdcast<int>(res.size());
}

int symmetrizeSSSR(ROMol &mol, bool includeDativeBonds,
                   bool includeHydrogenBonds) {
  VECT_INT_VECT tmp;
  return symmetrizeSSSR(mol, tmp, includeDativeBonds, includeHydrogenBonds);
};

int symmetrizeSSSR(ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds,
                   bool includeHydrogenBonds) {
  // FIX: need to set flag here the symmetrization has been done in order to
  // avoid repeating this work
  findSSSR(mol, res, includeDativeBonds, includeHydrogenBonds);

  // reinit as SYMM_SSSR
  auto ringInfo = mol.getRingInfo();
  ringInfo->initialize(FIND_RING_TYPE_SYMM_SSSR);

  // get the bond rings (we just calculated them when getting the SSSR.
  // copy is intentional, as we'll probably be adding more rings
  auto bondsssrs = ringInfo->bondRings();

  //
  // For each ring in the URFs, figure out if it could replace a single
  // ring in the SSSR. A ring could be swapped out if:
  //
  // * They are the same size
  // * The replacement doesn't remove any bonds from the union of the bonds
  //   in the SSSR.
  //
  // The latter can be checked by determining if the SSSR ring is the unique
  // provider of any ring bond. If it is, the replacement ring must also
  // provide that bond.
  //
  // May miss extra rings that would need to swap two (or three...) rings
  // to be included.

  // counts of each bond
  std::vector<int> bondCounts(mol.getNumBonds(), 0);
  for (const auto &r : bondsssrs) {
    for (const auto &b : r) {
      bondCounts[b] += 1;
    }
  }

  // This was also initialized when we got the SSSR.
  auto urfdata = mol.getRingInfo()->dp_urfData.get();

  INT_VECT extraRing;
  for (unsigned int i = 0; i < RDL_getNofURF(urfdata); ++i) {
    auto it = RDL_getRCyclesForURFIterator(urfdata, i);
    while (!RDL_cycleIteratorAtEnd(it)) {
      auto *cycle = RDL_cycleIteratorGetCycle(it);
      auto extraAtomRing = rdl_cycle_to_atom_ring(cycle);
      RingUtils::convertToBonds(extraAtomRing, extraRing, mol);

      if (std::ranges::find(bondsssrs, extraRing) !=
          bondsssrs.end()) {
        // already in the SSSR, skip
        RDL_deleteCycle(cycle);
        RDL_cycleIteratorNext(it);
        continue;
      }

      for (const auto &ring : bondsssrs) {
        if (ring.size() != cycle->weight) {
          continue;
        }

        // If `ring` is the only provider of some bond, extraRing must also
        // provide that bond.
        bool shareBond = false;
        bool replacesAllUniqueBonds = true;
        for (auto &bondID : ring) {
          const int bondCount = bondCounts[bondID];
          if (bondCount == 1 || !shareBond) {
            auto position = std::ranges::find(extraRing, bondID);
            if (position != extraRing.end()) {
              shareBond = true;
            } else if (bondCount == 1) {
              // 1 means `ring` is the only ring in the SSSR to provide this
              // bond, and extraRing did not provide it (so extraRing is not an
              // acceptable substitution in the SSSR for ring)
              replacesAllUniqueBonds = false;
            }
          }
        }

        if (shareBond && replacesAllUniqueBonds) {
          // small optimization: only convert the RDL cycle into
          // an atom ring if we are going to use it.
          res.push_back(extraAtomRing);
          FindRings::storeRingInfo(mol, extraAtomRing);
          break;
        }
      }

      RDL_deleteCycle(cycle);
      RDL_cycleIteratorNext(it);
    }
    RDL_deleteCycleIterator(it);
  }

  return rdcast<int>(res.size());
}

namespace {
void _DFS(const ROMol &mol, const Atom *atom, INT_VECT &atomColors,
          std::vector<const Atom *> &traversalOrder, VECT_INT_VECT &res,
          const Atom *fromAtom = nullptr) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(atomColors[atom->getIdx()] == 0, "bad color");
  atomColors[atom->getIdx()] = 1;
  traversalOrder.push_back(atom);

  for (const auto nbr : mol.atomNeighbors(atom)) {
    unsigned int nbrIdx = nbr->getIdx();
    if (atomColors[nbrIdx] == 0) {
      if (nbr->getDegree() < 2) {
        atomColors[nbr->getIdx()] = 2;
      } else {
        _DFS(mol, nbr, atomColors, traversalOrder, res, atom);
      }
    } else if (atomColors[nbrIdx] == 1) {
      if (fromAtom && nbrIdx != fromAtom->getIdx()) {
        INT_VECT cycle;
        auto lastElem =
            std::find(traversalOrder.rbegin(), traversalOrder.rend(), atom);
        for (auto rIt = lastElem;  // traversalOrder.rbegin();
             rIt != traversalOrder.rend() && (*rIt)->getIdx() != nbrIdx;
             ++rIt) {
          cycle.push_back((*rIt)->getIdx());
        }
        cycle.push_back(nbrIdx);
        res.push_back(cycle);
      }
    }
  }
  atomColors[atom->getIdx()] = 2;
  traversalOrder.pop_back();
}
}  // end of anonymous namespace
void fastFindRings(const ROMol &mol) {
  if (mol.getRingInfo()->isInitialized()) {
    mol.getRingInfo()->reset();
  }

  mol.getRingInfo()->initialize(FIND_RING_TYPE_FAST);

  VECT_INT_VECT res;
  res.resize(0);

  unsigned int nats = mol.getNumAtoms();

  INT_VECT atomColors(nats, 0);

  for (unsigned int i = 0; i < nats; ++i) {
    if (atomColors[i]) {
      continue;
    }
    if (mol.getAtomWithIdx(i)->getDegree() < 2) {
      atomColors[i] = 2;
      continue;
    }
    std::vector<const Atom *> traversalOrder;
    _DFS(mol, mol.getAtomWithIdx(i), atomColors, traversalOrder, res);
  }

  FindRings::storeRingsInfo(mol, res);
}

void findRingFamilies(const ROMol &mol, bool includeDativeBonds,
                      bool includeHydrogenBonds) {
  if (mol.getRingInfo()->isInitialized()) {
    // return if we've done this before
    if (mol.getRingInfo()->areRingFamiliesInitialized()) {
      return;
    }
  } else {
    mol.getRingInfo()->initialize();
  }

  // RDL_calculate fails and returns null if the graph is empty,
  // just trick it into not freaking out by giving it a fake atom
  auto numAtoms = mol.getNumAtoms();
  if (numAtoms == 0) {
    numAtoms = 1;
  }

  RDL_graph *graph = RDL_initNewGraph(numAtoms);
  for (auto cbi : mol.bonds()) {
    if (auto bt = cbi->getBondType();
        bt == Bond::ZERO || (!includeDativeBonds && isDative(bt)) ||
        (!includeHydrogenBonds && bt == Bond::HYDROGEN)) {
      continue;
    }

    RDL_addUEdge(graph, cbi->getBeginAtomIdx(), cbi->getEndAtomIdx());
  }
  RDL_data *urfdata = RDL_calculate(graph);
  if (urfdata == nullptr) {
    RDL_deleteGraph(graph);
    mol.getRingInfo()->dp_urfData.reset();
    throw ValueErrorException("Cannot get URFs");
  }
  mol.getRingInfo()->dp_urfData.reset(urfdata, &RDL_deleteData);

  for (unsigned int i = 0; i < RDL_getNofURF(urfdata); ++i) {
    RDL_node *nodes = nullptr;
    unsigned nNodes = RDL_getNodesForURF(urfdata, i, &nodes);
    if (nNodes == RDL_INVALID_RESULT) {
      free(nodes);
      throw ValueErrorException("Cannot get URF nodes");
    }
    RDL_edge *edges = nullptr;
    unsigned nEdges = RDL_getEdgesForURF(urfdata, i, &edges);
    if (nEdges == RDL_INVALID_RESULT) {
      free(nodes);
      free(edges);
      throw ValueErrorException("Cannot get URF edges");
    }
    INT_VECT nvect(nNodes), evect(nEdges);
    for (unsigned int ridx = 0; ridx < nNodes; ++ridx) {
      nvect[ridx] = nodes[ridx];
    }
    for (unsigned int ridx = 0; ridx < nEdges; ++ridx) {
      unsigned int bidx = edges[ridx][0];
      unsigned int eidx = edges[ridx][1];
      evect[ridx] = mol.getBondBetweenAtoms(bidx, eidx)->getIdx();
    }
    mol.getRingInfo()->addRingFamily(nvect, evect);
    free(nodes);
    free(edges);
  }
}
}  // namespace MolOps

}  // namespace RDKit
