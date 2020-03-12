//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <sstream>

#include <boost/range/iterator_range.hpp>

#include "RDKitCipMol.hpp"
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>

namespace RDKit
{
namespace NewCIPLabelling
{

RDKitCipMol::RDKitCipMol(ROMol* mol) : mol{mol} {}

int RDKitCipMol::getNumAtoms() const
{
    return mol->getNumAtoms();
}

int RDKitCipMol::getNumBonds() const
{
    return mol->getNumBonds();
};

RdkA RDKitCipMol::getAtom(int idx) const
{
    return mol->getAtomWithIdx(idx);
};

int RDKitCipMol::getAtomIdx(RdkA atom) const
{
    return atom->getIdx();
};

std::vector<RdkA> RDKitCipMol::atoms() const
{
    std::vector<RdkA> atoms;
    atoms.reserve(mol->getNumAtoms());
    for (auto& atom : mol->atoms()) {
        atoms.push_back(atom);
    }
    return atoms;
}

RdkB RDKitCipMol::getBond(int idx) const
{
    return mol->getBondWithIdx(idx);
};

int RDKitCipMol::getBondIdx(RdkB bond) const
{
    return bond->getIdx();
};

std::vector<RdkB> RDKitCipMol::getBonds(RdkA atom) const
{
    std::vector<RdkB> bonds;
    bonds.reserve(atom->getDegree());
    for (const auto& ibnd :
         boost::make_iterator_range(mol->getAtomBonds(atom))) {
        bonds.push_back((*mol)[ibnd]);
    }
    return bonds;
}

RdkA RDKitCipMol::getOther(RdkB bond, RdkA atom) const
{
    return bond->getOtherAtom(atom);
};

RdkA RDKitCipMol::getBeg(RdkB bond) const
{
    return bond->getBeginAtom();
};

RdkA RDKitCipMol::getEnd(RdkB bond) const
{
    return bond->getEndAtom();
};

bool RDKitCipMol::isInRing(RdkB bond) const
{
    const auto rings = mol->getRingInfo();
    return rings->numBondRings(bond->getIdx()) != 0u;
};

int RDKitCipMol::getAtomicNum(RdkA atom) const
{
    if (atom == nullptr) {
        return 1;
    }
    return atom->getAtomicNum();
};

int RDKitCipMol::getNumHydrogens(RdkA atom) const
{
    return atom->getTotalNumHs();
};

#if 0
***** NOT IMPLEMENTED YET *****
Descriptor RDKitCipMol::getAtomDescriptor(RdkA atom,
                                          const std::string& key) const {}

Descriptor RDKitCipMol::getBondDescriptor(RdkB bond,
                                          const std::string& key) const {}
#endif

int RDKitCipMol::getMassNum(RdkA atom) const
{
    if (atom == nullptr) {
        return 0;
    }
    return atom->getIsotope();
};

int RDKitCipMol::getCharge(RdkA atom) const
{
    return atom->getFormalCharge();
};

int RDKitCipMol::getBondOrder(RdkB bond) const
{
    auto order = bond->getBondTypeAsDouble();
    return static_cast<int>(std::ceil(order));
};

void RDKitCipMol::setAtomDescriptor(RdkA atom, const std::string& key,
                                    Descriptor desc)
{
    if (key == CIP_LABEL_KEY) {
        std::stringstream ss;
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
            ss << desc;
        }
        atom->setProp(common_properties::_CIPCode, ss.str());
    } else {
        std::stringstream ss;
        ss << "Setting property key '" << key
           << "' is not implemented for atoms.";
        throw std::runtime_error(ss.str());
    }
}

void RDKitCipMol::setBondDescriptor(RdkB bond, const std::string& key,
                                    Descriptor desc)
{
    // This is incorrect, but right now, we need don't know the anchors to set
    // seqCis/TRANSlabels

    if (key == CIP_LABEL_KEY) {
        std::stringstream ss;
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
            ss << desc;
        }
        bond->setProp(common_properties::_CIPCode, ss.str());
    } else {
        std::stringstream ss;
        ss << "Setting property key '" << key
           << "' is not implemented for bonds.";
        throw std::runtime_error(ss.str());
    }
}

} // namespace NewCIPLabelling
} // namespace RDKit