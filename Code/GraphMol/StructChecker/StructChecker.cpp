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
#include "StripSmallFragments.h"

namespace RDKit {
 namespace StructureCheck {

    unsigned StructChecker::checkMolStructure(RWMol &mol)const {
        unsigned flags = NO_CHANGE; // == 0. return value

        if (0 != Options.MaxMolSize
            && (mol.getNumAtoms() > Options.MaxMolSize || mol.getNumBonds() > Options.MaxMolSize)) {
            return SIZE_CHECK_FAILED;
        }
        if( ! mol.getRingInfo()->isInitialized())
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

        if (Options.RemoveMinorFragments)
        {
            AddMWMF(mol, "MW_PRE");            // Add  mol mass data field
            if (StripSmallFragments(mol)) {
                flags |= FRAGMENTS_FOUND;
            }
            AddMWMF(mol, "MW_POST");            // Add  mol mass data field

        }
/*
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
            if (Options.CheckStereo) {
                flags |= BAD_MOLECULE;
            }
            else
            {
                RemoveDubiousStereochemistry(mol);
                flags |= DUBIOUS_STEREO_REMOVED;
            }
        }
//line 1612
/*
        if (charges_read && TotalCharge(mol) != Options.DesiredCharge)
        {
            bool tmp = RechargeMolecule(mp, desired_charge, &ndeprot, &nrefine);
            if (mp->symbol_lists || mp->prop_lines)
                tmp = false;
            charge_bad = !tmp;
            flags |= (charge_bad ? BAD_MOLECULE : RECHARGED);
        }
*/
        double clashLimit = 0.15; // see void SetCollisionLimit(int percent)
        if (Options.CheckCollisions && AtomClash(mol, clashLimit))
            flags |= ATOM_CLASH;
/*
        if (! Options.GoodAtoms.empty() && ! CheckAtoms(mol, Options.GoodAtoms))
            flags |= ATOM_CHECK_FAILED;
*/
        if (Options.CheckStereo && ! CheckStereo(mol))
            flags |= STEREO_ERROR;

//        if (Options.GroupsToSGroups)
//            ConvertGroupsToSGroups(mol);
/*
//line 1630
        stereo_bad = FALSE;
        for (i = 0; i<nstereopat; i++)
        {
            ssp = stereo_patterns[i];
            tmp = ForceStereoTemplate(mp, ssp);
            if (tmp == (-1))
            {
                flags |= STEREO_FORCED_BAD; // problem enforcing stereochemistry of 'ssp->name'
            }
            else if (tmp == 15) // "STEREO_FORCED"
            {
                flags |= STEREO_TRANSFORMED;    // stereochemistry of 'ssp->name' enforced",
            }
        }
//line 1655
        for (i = 0; i<npat; i++)         // do template cleaning 
        {
            ssp = patterns[i];
            if (TemplateClean(mol, ssp))
            {
                result |= TEMPLATE_TRANSFORMED; // has been cleaned with template 'ssp->name'
            }
        }
//line 1669
        for (i = 0; i<nrpat; i++)         // do template rotation
        {
            ssp = rotate_patterns[i];
            if (TemplateRotate(mol, ssp))
            {
                result |= TEMPLATE_TRANSFORMED; // has been rotated by template 'ssp->name'
            }
        }
   }
*/

// the end:
        if (0 != (flags & TRANSFORMED)) {   // sanitaze molecule
//???? .............. ????
            if (mol.getRingInfo()->isInitialized())
                mol.getRingInfo()->reset();
            mol.getRingInfo()->initialize();
        }
        return flags;
    }

 }// namespace StructureCheck
} // namespace RDKit
