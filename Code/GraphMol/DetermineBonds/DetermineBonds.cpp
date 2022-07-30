//
//  Copyright (C) 2022 Sreya Gogineni and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "DetermineBonds.h"
#include <GraphMol/RDKitBase.h>
//#include "EHTTools.h"
#include <iostream>

namespace RDKit {


void determineConnectivity(RWMol &mol, bool useHuckel, int charge, double covFactor) { // accept reference
    unsigned int numAtoms = mol.getNumAtoms();
    for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
            mol.removeBond(i, j);
            mol.getAtomWithIdx(i)->setNoImplicit(false);
            mol.getAtomWithIdx(j)->setNoImplicit(false);
        }
    }
    if (useHuckel) {
        connectivityHuckel(mol, charge);
    } else {
        connectivityVdW(mol, covFactor);
    }
}

void connectivityHuckel(RWMol &mol, int charge) {
    unsigned int numAtoms = mol.getNumAtoms();
    mol.getAtomWithIdx(0)->setFormalCharge(charge);
    //
}

void connectivityVdW(RWMol &mol, double covFactor) {
    unsigned int numAtoms = mol.getNumAtoms();
    double *distMat = MolOps::get3DDistanceMat(mol);
    
    double rcov[numAtoms];
    for (unsigned int i = 0; i < numAtoms; i++) {
        rcov[i] = PeriodicTable::getTable()->getRcovalent(mol.getAtomWithIdx(i)->getAtomicNum());
        rcov[i] *= covFactor;
    }

    for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
            if (distMat[i * numAtoms + j] <= (rcov[i] + rcov[j])) {
                mol.addBond(i, j, Bond::BondType::SINGLE);
            }
        }
    }
}

} // namespace RDKit


