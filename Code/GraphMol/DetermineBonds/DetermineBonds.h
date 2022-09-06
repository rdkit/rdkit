//
//  Copyright (C) 2022 Sreya Gogineni and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_DETERMINEBONDS_H
#define RD_DETERMINEBONDS_H
#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
// ! assigns atomic connectivity to a molecule using atomic coordinates, disregarding pre-existing bonds
/*!
    \param mol is the molecule of interest; it must have a 3D conformer
    \param useHueckel (optional) if this is  \c true, the extended Hueckel theorem will be used to determine connectivity rather than the van der Waals method
    \param charge (optional) the charge of the molecule; it must be provided if the Hueckel method is used and charge is non-zero
    \param covFactor (optional) the factor with which to multiply each covalent radius if the van der Waals method is used
 */
void determineConnectivity(RWMol &mol, bool useHueckel=false, int charge=0, double covFactor=1.3);

void determineBondOrder(RWMol &mol, int charge=0, bool allowChargedFragments=true, bool embedChiral=true, bool useAtomMap=false);

}

#endif
