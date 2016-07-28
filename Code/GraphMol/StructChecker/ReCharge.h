//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

namespace RDKit {
 namespace StructureCheck {

     /*
     * Removes hydrogens from *mp until desired_charge is reached. The
     * positions for hydrogen removal are selected by "acidity" combined
     * with a refinement algorithm. It returns TRUE if molecule could be
     * neutralized and FALSE if any problem were encountered.
     * *ndeprot and *nrefine are set to the number of deprotonations
     * and refinement cycles performed.
     */
     bool RechargeMolecule(RWMol &mol, int desired_charge);

     /*
     * Returns the total charge of all atoms in molecule.
     */
     int  TotalCharge     (const ROMol &mol);
 }
}

