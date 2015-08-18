// $Id$
//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include <map>
#include <vector>
#include <string>
#include <stdexcept>
#include <boost/cstdint.hpp>
#include "../RDKitBase.h"

namespace RDKit
{
 namespace MolHash
 {
    typedef boost::uint32_t HashCodeType;

    HashCodeType generateMoleculeHashCode(const ROMol &mol,
        const std::vector<unsigned> *atomsToUse=0,
        const std::vector<unsigned> *bondsToUse=0,  // ?? listed bonds between/to/from excluded atom(s) ??
        const std::vector<boost::uint32_t> *atomCodes =0,
        const std::vector<boost::uint32_t> *bondCodes =0);

    enum CodeFlags  // bitwise flags to combine and compute atom/bond codes
    {
        CF_NO_LABELS= 0x0000,
        CF_ELEMENT  = 0x0001,
        CF_CHARGE   = 0x0002,
        CF_VALENCE  = 0x0004,
        CF_ISOTOPE  = 0x0008,
        CF_ATOM_CHIRALITY= 0x0010,
        CF_ATOM_AROMATIC = 0x0020,
        CF_ATOM_ALL = 0x00FF,
        CF_BOND_ORDER         = 0x0100, // ignore AROMATIZATION if corresponding flag is not specified
        CF_BOND_AROMATIZATION = 0x0200,
        CF_BOND_TYPE_EXACT    = CF_BOND_ORDER | CF_BOND_AROMATIZATION, // exact type value with aromatic
        CF_BOND_CHIRALITY     = 0x0400, // include bond chirality information into bond code
        CF_BOND_IN_RING       = 0x0800,
        CF_BOND_ALL           = 0xFF00,
        CF_ALL = 0xFFFF,
    };

    void fillAtomBondCodes(const ROMol &mol, boost::uint64_t flags     // CodeFlags constants combination
                          , std::vector<boost::uint32_t> *atomCodes    // NULL is allowed
                          , std::vector<boost::uint32_t> *bondCodes);  // NULL is allowed

#pragma pack(push,1)
    struct HashSet
    {
        boost::uint16_t  Version;
        boost::uint16_t  Reserved;
        boost::uint16_t  NumAtoms;
        boost::uint16_t  NumBonds;
        boost::uint32_t  FormulaCRC32;
        HashCodeType     NonChiralAtomsHash;
        HashCodeType     NonChiralBondsHash;
        HashCodeType     ChiralAtomsHash;
        HashCodeType     ChiralBondsHash;
        HashCodeType     ChiralityHash;
    public:
        HashSet() { memset(this, 0, sizeof(*this));}
    };
#pragma pack(pop)

    void generateMoleculeHashSet(const ROMol &mol, HashSet& res,
        const std::vector<unsigned> *atomsToUse=0,
        const std::vector<unsigned> *bondsToUse=0);

    std::string generateMoleculeHashSet(const ROMol &mol,
                        const std::vector<unsigned> *atomsToUse=0,
                        const std::vector<unsigned> *bondsToUse=0);

    std::string encode(const void* bin, size_t size);   // binary data to Base64 encoded string

}}
