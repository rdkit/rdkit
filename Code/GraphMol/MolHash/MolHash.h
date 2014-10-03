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
#include "../RDKitBase.h"

namespace RDKit
{
 namespace MolHash
 {
    typedef unsigned long HashCodeT;    // 32 bits !!!
//    typedef unsigned long long HashCodeT;

    HashCodeT generateMoleculeHashCode(const ROMol &mol,
        const std::vector<unsigned> *atomsToUse=0,
        const std::vector<unsigned> *bondsToUse=0,  // ?? listed bonds between/to/from excluded atom(s) ??
        const std::vector<unsigned> *atomCodes =0,
        const std::vector<unsigned> *bondCodes =0);

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

    void fillAtomBondCodes(const ROMol &mol, unsigned flags     // CodeFlags constants combination
                          , std::vector<unsigned> *atomCodes    // NULL is allowed
                          , std::vector<unsigned> *bondCodes);  // NULL is allowed

#pragma pack(push,1)
    struct HashSet
    {
        unsigned short  Version;
        unsigned short  Reserved;
        unsigned short  NumAtoms;
        unsigned short  NumBonds;
        unsigned long   FormulaCRC32;
        HashCodeT       NonChiralAtomsHash;
        HashCodeT       NonChiralBondsHash;
        HashCodeT       ChiralAtomsHash;
        HashCodeT       ChiralBondsHash;
        HashCodeT       ChiralityHash;
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
