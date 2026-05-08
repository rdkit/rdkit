
//
//  Copyright (C) 2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>

namespace RDKit {

using AtomsIterator = CXXAtomIterator<MolGraph, Atom *>;
using AtomNeighborsIterator =
    CXXAtomIterator<MolGraph, Atom *, MolGraph::adjacency_iterator>;
using BondsIterator = CXXBondIterator<MolGraph, Bond *>;
using AtomBondsIterator =
    CXXBondIterator<MolGraph, Bond *, MolGraph::out_edge_iterator>;

template <typename IterT=AtomsIterator>
struct AtomSeqHolder {
  IterT iter;
  ROMol *mol;
  AtomSeqHolder(ROMol &mol) : iter(mol.atoms()), mol(&mol) {};
  AtomSeqHolder(ROMol &mol, const IterT &iter) : iter(iter), mol(&mol) {};
  size_t size() const { return iter.size(); }
  typename IterT::CXXAtomIter begin() { return iter.begin(); }
  typename IterT::CXXAtomIter end() { return iter.end(); }
  Atom *operator[](int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= mol->getNumAtoms()) {
      throw IndexErrorException(idx);
    }
    return mol->getAtomWithIdx(idx);
  }
};

template <typename IterT=BondsIterator>
struct BondSeqHolder {
  IterT iter;
  ROMol *mol;
  BondSeqHolder(ROMol &mol) : iter(mol.bonds()), mol(&mol) {};
  BondSeqHolder(ROMol &mol, const IterT &iter) : iter(iter), mol(&mol) {};
  size_t size() const { return iter.size(); }
  typename IterT::CXXBondIter begin() { return iter.begin(); }
  typename IterT::CXXBondIter end() { return iter.end(); }
  Bond *operator[](int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= mol->getNumBonds()) {
      throw IndexErrorException(idx);
    }
    return mol->getBondWithIdx(idx);
  }
};
}  // namespace RDKit