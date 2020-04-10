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

#include "Descriptor.hpp"
#include "Mancude.h"

namespace RDKit {

namespace CIPLabeler {

template <typename T, typename U> class CIPMolIterator {
public:
  class CIPMolIter {
  private:
    ROMol *mol;
    U pos;
    T current = nullptr;

  public:
    CIPMolIter() = delete;
    CIPMolIter(ROMol *mol, U pos) : mol{mol}, pos{pos} {}

    T &operator*() {
      current = (*mol)[*pos];
      return current;
    }

    CIPMolIter &operator++() {
      ++pos;
      return *this;
    }

    bool operator!=(const CIPMolIter &it) const { return pos != it.pos; }
  };

private:
  ROMol *mol;
  U istart;
  U iend;

public:
  CIPMolIterator() = delete;
  CIPMolIterator(ROMol *mol, std::pair<U, U> itr)
      : mol{mol}, istart{itr.first}, iend{itr.second} {}

  CIPMolIter begin() { return {mol, istart}; }
  CIPMolIter end() { return {mol, iend}; }
};

class CIPMol {
private:
  ROMol *mol;
  std::unique_ptr<RWMol> kekulized_mol = nullptr;

  std::vector<Fraction> atomnums;

public:
  explicit CIPMol(ROMol *mol);

  Fraction getFractionalAtomicNum(Atom *atom) const;

  unsigned getNumAtoms() const;

  unsigned getNumBonds() const;

  Atom *getAtom(int idx) const;

  CXXAtomIterator<MolGraph, Atom *> atoms() const;

  Bond *getBond(int idx) const;

  CIPMolIterator<Bond *, ROMol::OEDGE_ITER> getBonds(Atom *atom) const;

  CIPMolIterator<Atom *, ROMol::ADJ_ITER> getNeighbors(Atom *atom) const;

  bool isInRing(Bond *bond) const;

  int getBondOrder(Bond *bond) const;

  void setAtomDescriptor(Atom *atom, const std::string &key, Descriptor desc);

  void setBondDescriptor(Bond *bond, const std::string &key, Descriptor desc);
};

} // namespace CIPLabeler
} // namespace RDKit