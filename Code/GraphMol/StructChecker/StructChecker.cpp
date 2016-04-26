//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "StructChecker.h"

namespace RDKit {
 namespace StructureCheck {

    unsigned StructChecker::checkMolStructure(RWMol &mol)const {
        unsigned flags = StructureFlags::NO_CHANGE; // == 0. return value

        if ( 0 != Options.MaxMolSize && mol.getNumAtoms() > 0
            &&(mol.getNumAtoms() > Options.MaxMolSize || mol.getNumBonds() > Options.MaxMolSize))
            return StructureFlags::BAD_MOLECULE;
/* uses SDL text
        if (Options.ConvertAtomTexts)
        {
            if(!convertAtomAliases(mol))
                flags |= ALIAS_CONVERSION_FAILED;
            else
                flags |= TRANSFORMED;
        }

        if (Options.ConvertSText)
            ;//new_data_list = ConvertSTEXTToData(mol, new_data_list);
*/
        return flags;
    }

 }// namespace StructureCheck
} // namespace RDKit
