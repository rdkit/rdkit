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
#pragma once

#include <ostream>

namespace RDKit
{
namespace NewCIPLabelling
{

static const std::string CIP_LABEL_KEY = "cip.label";
static const std::string CONF_INDEX = "conf.index";

/**
 * Defines a descriptor which can be assigned to an atom to indicate the type of
 * chirality (if there is any). Each descriptor defines it's general @{link
 * Type} which can be useful when comparing centres of different geometry.
 *
 */
enum class Descriptor {
    NONE, // Unspecified
    UNKNOWN,
    ns, // Other
    /**
     * Tetrahedral
     */
    R,
    S,
    r,
    s,
    /**
     * Cis/Trans
     */
    seqTrans,
    seqCis,
    E,
    Z,
    /* Axial */
    M,
    P,
    m,
    p,

    SP_4,
    TBPY_5,
    OC_6
};

static std::ostream& operator<<(std::ostream& os, const Descriptor& desc)
{
    switch (desc) {
    case Descriptor::NONE:
        os << "NONE";
        break;
    case Descriptor::UNKNOWN:
        os << "UNKNOWN";
        break;
    case Descriptor::ns:
        os << "ns";
        break;
    case Descriptor::R:
        os << 'R';
        break;
    case Descriptor::S:
        os << 'S';
        break;
    case Descriptor::r:
        os << 'r';
        break;
    case Descriptor::s:
        os << 's';
        break;
    case Descriptor::seqTrans:
        os << "seqTrans";
        break;
    case Descriptor::seqCis:
        os << "seqCis";
        break;
    case Descriptor::E:
        os << 'E';
        break;
    case Descriptor::Z:
        os << 'Z';
        break;
    case Descriptor::M:
        os << 'M';
        break;
    case Descriptor::P:
        os << 'P';
        break;
    case Descriptor::m:
        os << 'm';
        break;
    case Descriptor::p:
        os << 'p';
        break;
    case Descriptor::SP_4:
        os << "SP_4";
        break;
    case Descriptor::TBPY_5:
        os << "TBPY_5";
        break;
    case Descriptor::OC_6:
        os << "OC_6";
        break;
    }
    return os;
}

} // namespace NewCIPLabelling
} // namespace RDKit