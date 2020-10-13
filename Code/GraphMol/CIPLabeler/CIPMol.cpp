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

#include <GraphMol/MolOps.h>

#include "CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

CIPMol::CIPMol(ROMol &mol) : d_mol{mol} {}

boost::rational<int> CIPMol::getFractionalAtomicNum(Atom *atom) const {
  PRECONDITION(atom, "bad atom")
  if (d_atomnums.empty())
    const_cast<CIPMol *>(this)->d_atomnums = calcFracAtomNums(*this);
  return d_atomnums[atom->getIdx()];
}

unsigned CIPMol::getNumAtoms() const { return d_mol.getNumAtoms(); }

unsigned CIPMol::getNumBonds() const { return d_mol.getNumBonds(); };

Atom *CIPMol::getAtom(int idx) const { return d_mol.getAtomWithIdx(idx); };

Bond *CIPMol::getBond(int idx) const { return d_mol.getBondWithIdx(idx); };

CIPMolSpan<Bond *, ROMol::OEDGE_ITER> CIPMol::getBonds(Atom *atom) const {
  PRECONDITION(atom, "bad atom")
  return {d_mol, d_mol.getAtomBonds(atom)};
}

CIPMolSpan<Atom *, ROMol::ADJ_ITER> CIPMol::getNeighbors(Atom *atom) const {
  PRECONDITION(atom, "bad atom")
  return {d_mol, d_mol.getAtomNeighbors(atom)};
}

bool CIPMol::isInRing(Bond *bond) const {
  PRECONDITION(bond, "bad bond")
  const auto rings = d_mol.getRingInfo();

  if (!rings->isInitialized()) {
    MolOps::fastFindRings(d_mol);
  }

  return rings->numBondRings(bond->getIdx()) != 0u;
};

int CIPMol::getBondOrder(Bond *bond) const {
  PRECONDITION(bond, "bad bond")
  if (dp_kekulized_mol == nullptr) {
    auto tmp = new RWMol(d_mol);
    MolOps::Kekulize(*tmp);
    const_cast<CIPMol *>(this)->dp_kekulized_mol.reset(tmp);
  }

  const auto kekulized_bond = dp_kekulized_mol->getBondWithIdx(bond->getIdx());

  // Dative bonds might need to be considered with a different bond order
  // for the end atom at the end of the bond.
  switch (kekulized_bond->getBondType()) {
  case Bond::ZERO:
  case Bond::HYDROGEN:
  case Bond::DATIVE:
  case Bond::DATIVEL:
  case Bond::DATIVER:
    return 0;
  case Bond::SINGLE:
    return 1;
  case Bond::DOUBLE:
    return 2;
  case Bond::TRIPLE:
    return 3;
  case Bond::QUADRUPLE:
    return 4;
  case Bond::QUINTUPLE:
    return 5;
  case Bond::HEXTUPLE:
    return 6;
  default:
    throw std::runtime_error("Non integer-order bonds are not allowed.");
  }
};

} // namespace CIPLabeler
} // namespace RDKit
