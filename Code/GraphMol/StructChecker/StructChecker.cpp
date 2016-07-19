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
#include "Stereo.h"

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

        unsigned stereo_result = DubiousStereochemistry(mol);
        if (0!=(FixDubious3DMolecule(mol) & CONVERTED_TO_2D))
        {
            stereo_result = 1;
            flags |= DUBIOUS_STEREO_REMOVED;
        }
/*
        if (remove_minor_fragments)
        {
            // Add MF_PRE and MW_PRE data fields
            new_data_list = AddMWMF(new_data_list, mp, "MW_PRE", "MF_PRE");
            new_mp = StripSmallFragments(CopyMolecule(mp), &fragments_found);
            if (new_mp)      // new molecule data structure has been allocated
            {                // => need to overwrite components pointed to by input pointer
                             //    and deallocate skeleton container object
                             // deallocate objects pointed to by components of input (*mp)
                FreeMoleculeChildObjects(mp);
                // shallow copy (*new_mp) into 'emptied' container
                (*mp) = (*new_mp);
                // deallocate new container struct
                MyFree((char *)new_mp);
            }
            if (fragments_found) result |= FRAGMENTS_FOUND;
            // Add MF_POST and MW_POST data fields
            new_data_list = AddMWMF(new_data_list, mp, "MW_POST", "MF_POST");
        }

        for (i = 0; i<ntautomers; i++)         // do tautomer standardization 
        {
            fprintf(stderr, "tautomerizing with rule %d\n", i);
            for (j = 0;j<3;j++)      // limit to 3 run per rule
            {
                tmp = ApplyTautomer(mp, from_tautomer[i], to_tautomer[i]);
                if (!tmp) break;
                result |= TAUTOMER_TRANSFORMED;
                sprintf(msg_buffer,
                    "%10s: has been tautomerized with rule '%s'",
                    mp->name, from_tautomer[i]->name);
                AddMsgToList(msg_buffer);
            }
        }
*/
/*        if (!IsNULL(data_list) && !IsNULL(new_data_list))
        {        // append new data list if any
            for (dph = data_list; !IsNULL(dph->next); dph = dph->next)
                ;
            dph->next = new_data_list;
        }
*/
        if (stereo_result == EITHER_BOND_FOUND) {  // looks for EITHER bonds
            flags |= EITHER_WARNING;
            RemoveDubiousStereochemistry(mol);
            flags |= DUBIOUS_STEREO_REMOVED;
        }
        else if (stereo_result > EITHER_BOND_FOUND) { // more severe errors
            flags |= STEREO_ERROR;
            if (Options.CheckStereo)
                flags |= BAD_MOLECULE;
            else
            {
                RemoveDubiousStereochemistry(mol);
                flags |= DUBIOUS_STEREO_REMOVED;
            }
        }


// ...............


// the end:
        if (0 != (flags & TRANSFORMED)) {   // sanitaze molecule
//???? .............. ????
            mol.getRingInfo()->initialize();
        }
        return flags;
    }

 }// namespace StructureCheck
} // namespace RDKit
