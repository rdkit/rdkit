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

CIPMol::CIPMol(ROMol *mol) : mol{mol} {}

Fraction CIPMol::getFractionalAtomicNum(Atom *atom) const {
  PRECONDITION(atom, "bad atom")
  if (atomnums.empty())
    const_cast<CIPMol *>(this)->atomnums = calcFracAtomNums(this);
  return atomnums[atom->getIdx()];
}

unsigned CIPMol::getNumAtoms() const { return mol->getNumAtoms(); }

unsigned CIPMol::getNumBonds() const { return mol->getNumBonds(); };

Atom *CIPMol::getAtom(int idx) const { return mol->getAtomWithIdx(idx); };

CXXAtomIterator<MolGraph, Atom *> CIPMol::atoms() const { return mol->atoms(); }

Bond *CIPMol::getBond(int idx) const { return mol->getBondWithIdx(idx); };

CIPMolIterator<Bond *, ROMol::OEDGE_ITER> CIPMol::getBonds(Atom *atom) const {
  return {mol, mol->getAtomBonds(atom)};
}

CIPMolIterator<Atom *, ROMol::ADJ_ITER> CIPMol::getNeighbors(Atom *atom) const {
  return {mol, mol->getAtomNeighbors(atom)};
}

bool CIPMol::isInRing(Bond *bond) const {
  const auto rings = mol->getRingInfo();

  if (!rings->isInitialized()) {
    MolOps::findSSSR(*mol);
  }

  return rings->numBondRings(bond->getIdx()) != 0u;
};

int CIPMol::getBondOrder(Bond *bond) const {
  if (kekulized_mol == nullptr) {
    auto tmp = new RWMol(*mol);
    MolOps::Kekulize(*tmp);
    const_cast<CIPMol *>(this)->kekulized_mol.reset(tmp);
  }

  const auto kekulized_bond = kekulized_mol->getBondWithIdx(bond->getIdx());
  switch (kekulized_bond->getBondType()) {
  case Bond::ZERO:
  case Bond::HYDROGEN:
  case Bond::DATIVE:
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

void CIPMol::setAtomDescriptor(Atom *atom, const std::string &key,
                               Descriptor desc) {
  if (key == CIP_LABEL_KEY) {
    switch (desc) {
    case Descriptor::NONE:
      throw std::runtime_error("Received an invalid as Atom Descriptor");
    case Descriptor::seqTrans:
    case Descriptor::seqCis:
    case Descriptor::E:
    case Descriptor::Z:
    case Descriptor::M:
    case Descriptor::P:
    case Descriptor::m:
    case Descriptor::p:
    case Descriptor::SP_4:
    case Descriptor::TBPY_5:
    case Descriptor::OC_6:
      throw std::runtime_error(
          "Received a Descriptor that is not supported for atoms");
    default:
      atom->setProp(common_properties::_CIPCode, to_string(desc));
    }
  } else {
    std::stringstream ss;
    ss << "Setting property key '" << key << "' is not implemented for atoms.";
    throw std::runtime_error(ss.str());
  }
}

void CIPMol::setBondDescriptor(Bond *bond, const std::string &key,
                               Descriptor desc) {
  // This is incorrect, but right now, we need don't know the anchors to set
  // seqCis/TRANSlabels

  if (key == CIP_LABEL_KEY) {
    switch (desc) {
    case Descriptor::NONE:
      throw std::runtime_error("Received an invalid as Bond Descriptor");
    case Descriptor::R:
    case Descriptor::S:
    case Descriptor::r:
    case Descriptor::s:
    case Descriptor::SP_4:
    case Descriptor::TBPY_5:
    case Descriptor::OC_6:
      throw std::runtime_error(
          "Received a Descriptor that is not supported for bonds");
    default:
      bond->setProp(common_properties::_CIPCode, to_string(desc));
    }
  } else {
    std::stringstream ss;
    ss << "Setting property key '" << key << "' is not implemented for bonds.";
    throw std::runtime_error(ss.str());
  }
}

} // namespace CIPLabeler
} // namespace RDKit