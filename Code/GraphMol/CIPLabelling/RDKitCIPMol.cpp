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

#include <boost/range/iterator_range.hpp>

#include "RDKitCIPMol.hpp"
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace CIPLabelling {

RDKitCIPMol::RDKitCIPMol(ROMol *mol) : mol{mol} {}

int RDKitCIPMol::getNumAtoms() const { return mol->getNumAtoms(); }

int RDKitCIPMol::getNumBonds() const { return mol->getNumBonds(); };

RdkA RDKitCIPMol::getAtom(int idx) const { return mol->getAtomWithIdx(idx); };

int RDKitCIPMol::getAtomIdx(RdkA atom) const { return atom->getIdx(); };

std::vector<RdkA> RDKitCIPMol::atoms() const {
  std::vector<RdkA> atoms;
  atoms.reserve(mol->getNumAtoms());
  for (auto &atom : mol->atoms()) {
    atoms.push_back(atom);
  }
  return atoms;
}

RdkB RDKitCIPMol::getBond(int idx) const { return mol->getBondWithIdx(idx); };

int RDKitCIPMol::getBondIdx(RdkB bond) const { return bond->getIdx(); };

std::vector<RdkB> RDKitCIPMol::getBonds(RdkA atom) const {
  std::vector<RdkB> bonds;
  bonds.reserve(atom->getDegree());
  for (const auto &ibnd : boost::make_iterator_range(mol->getAtomBonds(atom))) {
    bonds.push_back((*mol)[ibnd]);
  }
  return bonds;
}

RdkA RDKitCIPMol::getOther(RdkB bond, RdkA atom) const {
  return bond->getOtherAtom(atom);
};

RdkA RDKitCIPMol::getBeg(RdkB bond) const { return bond->getBeginAtom(); };

RdkA RDKitCIPMol::getEnd(RdkB bond) const { return bond->getEndAtom(); };

bool RDKitCIPMol::isInRing(RdkB bond) const {
  const auto rings = mol->getRingInfo();
  return rings->numBondRings(bond->getIdx()) != 0u;
};

int RDKitCIPMol::getAtomicNum(RdkA atom) const {
  if (atom == nullptr) {
    return 1;
  }
  return atom->getAtomicNum();
};

int RDKitCIPMol::getNumHydrogens(RdkA atom) const {
  return atom->getTotalNumHs();
};

#if 0
***** NOT IMPLEMENTED YET *****
Descriptor RDKitCIPMol::getAtomDescriptor(RdkA atom,
                                          const std::string& key) const {}

Descriptor RDKitCIPMol::getBondDescriptor(RdkB bond,
                                          const std::string& key) const {}
#endif

int RDKitCIPMol::getMassNum(RdkA atom) const {
  if (atom == nullptr) {
    return 0;
  }
  return atom->getIsotope();
};

int RDKitCIPMol::getCharge(RdkA atom) const { return atom->getFormalCharge(); };

int RDKitCIPMol::getBondOrder(RdkB bond) const {
  if (kekulized_mol == nullptr) {
    auto tmp = new RWMol(*mol);
    MolOps::Kekulize(*tmp);
    const_cast<RDKitCIPMol *>(this)->kekulized_mol.reset(tmp);
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

void RDKitCIPMol::setAtomDescriptor(RdkA atom, const std::string &key,
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

void RDKitCIPMol::setBondDescriptor(RdkB bond, const std::string &key,
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

} // namespace CIPLabelling
} // namespace RDKit