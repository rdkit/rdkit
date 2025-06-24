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
#pragma once

#include <memory>

#include <GraphMol/RDKitBase.h>

#include "Descriptor.h"
#include "Mancude.h"

namespace RDKit {

namespace CIPLabeler {

template <typename T, typename U>
class CIPMolSpan {
 public:
  class CIPMolIter {
   public:
    CIPMolIter() = delete;
    CIPMolIter(ROMol &mol, U pos) : d_mol{mol}, d_pos{std::move(pos)} {}

    T &operator*() {
      d_current = d_mol[*d_pos];
      return d_current;
    }

    CIPMolIter &operator++() {
      ++d_pos;
      return *this;
    }

    bool operator!=(const CIPMolIter &it) const { return d_pos != it.d_pos; }

   private:
    ROMol &d_mol;
    U d_pos;
    T d_current = nullptr;
  };

 public:
  CIPMolSpan() = delete;
  CIPMolSpan(ROMol &mol, std::pair<U, U> &&itr)
      : d_mol{mol},
        d_istart{std::move(itr.first)},
        d_iend{std::move(itr.second)} {}

  CIPMolIter begin() { return {d_mol, d_istart}; }
  CIPMolIter end() { return {d_mol, d_iend}; }

 private:
  ROMol &d_mol;
  const U d_istart;
  const U d_iend;
};

class CIPMol {
 public:
  CIPMol() = delete;

  explicit CIPMol(ROMol &mol);

  // Average atomic number with other atoms that are in an
  // aromatic ring with this one.
  boost::rational<int> getFractionalAtomicNum(Atom *atom) const;

  unsigned getNumAtoms() const;

  unsigned getNumBonds() const;

  Atom *getAtom(int idx) const;

  CXXAtomIterator<MolGraph, Atom *> atoms() const;

  Bond *getBond(int idx) const;

  CIPMolSpan<Bond *, ROMol::OEDGE_ITER> getBonds(Atom *atom) const;

  CIPMolSpan<Atom *, ROMol::ADJ_ITER> getNeighbors(Atom *atom) const;

  bool isInRing(Bond *bond) const;

  // Integer bond order of a kekulized molecule
  // Dative bonds get bond order 0.
  int getBondOrder(Bond *bond) const;

 private:
  ROMol &d_mol;
  std::vector<RDKit::Bond::BondType> d_kekulized_bonds;
  std::vector<boost::rational<int>> d_atomnums;
};

}  // namespace CIPLabeler
}  // namespace RDKit
