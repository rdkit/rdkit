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
#include "Pattern.h"

namespace RDKit {
 namespace StructureCheck {

    unsigned StructChecker::checkMolStructure(RWMol &mol)const {
        unsigned flags = NO_CHANGE; // == 0. return value

        if (0 != Options.MaxMolSize
            && (mol.getNumAtoms() > Options.MaxMolSize || mol.getNumBonds() > Options.MaxMolSize)) {
            return SIZE_CHECK_FAILED;
        }
        mol.getRingInfo()->initialize();

/* it uses SDL text
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
        if (!Options.AugmentedAtomPairs.empty()) {
            if (TransformAugmentedAtoms(mol, Options.AugmentedAtomPairs))
                flags |= TRANSFORMED;
        }

/* TODO:
        stereo_result = DubiousStereochemistry(mp);
        if (FixDubious3DMolecule(mp) & CONVERTED_TO_2D)
        {
            stereo_result = TRUE;
            flags |= DUBIOUS_STEREO_REMOVED;
        }
*/
        if (0 != (flags & TRANSFORMED)) {   // sanitaze molecule
//???? .............. ????
            mol.getRingInfo()->initialize();
        }
        return flags;
    }

 }// namespace StructureCheck
} // namespace RDKit
