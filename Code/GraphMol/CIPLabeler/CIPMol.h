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
  public:
    CIPMolIter() = delete;
    CIPMolIter(ROMol *mol, U pos) : dp_mol{mol}, d_pos{pos} {}

    T &operator*() {
      d_current = (*dp_mol)[*d_pos];
      return d_current;
    }

    CIPMolIter &operator++() {
      ++d_pos;
      return *this;
    }

    bool operator!=(const CIPMolIter &it) const { return d_pos != it.d_pos; }

  private:
    ROMol *dp_mol;
    U d_pos;
    T d_current = nullptr;
  };

public:
  CIPMolIterator() = delete;
  CIPMolIterator(ROMol *mol, std::pair<U, U> itr)
      : dp_mol{mol}, d_istart{itr.first}, d_iend{itr.second} {}

  CIPMolIter begin() { return {dp_mol, d_istart}; }
  CIPMolIter end() { return {dp_mol, d_iend}; }

private:
  ROMol *dp_mol;
  U d_istart;
  U d_iend;
};

class CIPMol {
public:
  explicit CIPMol(ROMol *mol);

  boost::rational<int> getFractionalAtomicNum(Atom *atom) const;

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

private:
  ROMol *dp_mol;
  std::unique_ptr<RWMol> dp_kekulized_mol = nullptr;

  std::vector<boost::rational<int>> d_atomnums;
};

} // namespace CIPLabeler
} // namespace RDKit